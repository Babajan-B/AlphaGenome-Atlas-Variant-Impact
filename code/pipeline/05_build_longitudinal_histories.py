#!/usr/bin/env python3
"""Aggregate release observations into variant histories and frozen Objective 1 cohorts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import duckdb


OUT_SCHEMA = pa.schema([
    ("canonical_id", pa.string()), ("variant_id", pa.string()), ("geneinfo", pa.string()),
    ("chrom", pa.string()), ("norm_pos", pa.int64()), ("norm_ref", pa.string()), ("norm_alt", pa.string()),
    ("first_observed", pa.string()), ("last_observed", pa.string()),
    ("first_class", pa.string()), ("last_class", pa.string()),
    ("first_vus_date", pa.string()), ("last_vus_date", pa.string()),
    ("first_reclassification_date", pa.string()), ("trajectory", pa.string()),
    ("observed_releases", pa.int32()), ("class_changes", pa.int32()),
    ("reversions", pa.int32()), ("representation_count", pa.int32()),
    ("high_confidence", pa.bool_()),
])


def review_high(value: str) -> bool:
    text = (value or "").lower()
    return "expert_panel" in text or "practice_guideline" in text


def summarize(rows: list[dict]) -> dict:
    rows.sort(key=lambda r: r["release_date"])
    classes = [r["classification_group"] for r in rows]
    compact = [classes[0]]
    for value in classes[1:]:
        if value != compact[-1]:
            compact.append(value)
    first_vus = next((r["release_date"] for r in rows if r["classification_group"] == "VUS"), None)
    vus_rows = [r for r in rows if r["classification_group"] == "VUS"]
    last_vus = vus_rows[-1]["release_date"] if vus_rows else None
    after_vus = [r for r in rows if first_vus and r["release_date"] > first_vus and r["classification_group"] in {"P_LP", "B_LB", "CONFLICTING"}]
    reclass = after_vus[0] if after_vus else None
    if first_vus:
        if reclass:
            trajectory = f"VUS_TO_{reclass['classification_group']}"
        elif rows[-1]["classification_group"] == "VUS":
            trajectory = "PERSISTENT_VUS"
        else:
            trajectory = "VUS_OTHER_UNRESOLVED"
    elif all(c == "P_LP" for c in classes):
        trajectory = "STABLE_P_LP"
    elif all(c == "B_LB" for c in classes):
        trajectory = "STABLE_B_LB"
    else:
        trajectory = "OTHER"
    reversions = max(0, len(compact) - 2)
    return {
        "canonical_id": rows[0]["canonical_id"], "variant_id": rows[0]["variant_id"],
        "geneinfo": next((r["geneinfo"] for r in reversed(rows) if r["geneinfo"]), ""),
        "chrom": rows[0]["chrom"], "norm_pos": rows[0]["norm_pos"],
        "norm_ref": rows[0]["norm_ref"], "norm_alt": rows[0]["norm_alt"],
        "first_observed": rows[0]["release_date"], "last_observed": rows[-1]["release_date"],
        "first_class": classes[0], "last_class": classes[-1], "first_vus_date": first_vus,
        "last_vus_date": last_vus, "first_reclassification_date": reclass["release_date"] if reclass else None,
        "trajectory": trajectory, "observed_releases": len(rows), "class_changes": len(compact) - 1,
        "reversions": reversions, "representation_count": len({(r['source_pos'],r['source_ref'],r['source_alt']) for r in rows}),
        "high_confidence": bool(reclass and review_high(reclass.get("review_status", ""))),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("/path/to/project_root/AG-VUS"))
    args = parser.parse_args()
    source = args.root / "data/processed/release_observations"
    fragments = list(source.glob("release_date=*/observations.parquet"))
    if len(fragments) != 420:
        raise RuntimeError(f"Expected 420 release partitions, found {len(fragments)}")
    out_dir = args.root / "data/processed/longitudinal"
    out_dir.mkdir(parents=True, exist_ok=True)
    output = out_dir / "variant_histories.parquet"
    writer = pq.ParquetWriter(output, OUT_SCHEMA, compression="zstd")
    database = out_dir / "objective1_aggregation.duckdb"
    connection = duckdb.connect(str(database))
    connection.execute(f"PRAGMA temp_directory='{args.root / 'tmp/duckdb'}'")
    connection.execute("PRAGMA threads=16")
    parquet_glob = str(source / "release_date=*" / "observations.parquet")
    reader = connection.execute(
        "SELECT * FROM read_parquet(?, hive_partitioning=true) "
        "ORDER BY canonical_id, release_date", [parquet_glob]
    ).fetch_record_batch(rows_per_batch=100_000)
    batch, current_key, current = [], None, []
    try:
        for record_batch in reader:
            for row in record_batch.to_pylist():
                key = row["canonical_id"]
                if current_key is not None and key != current_key:
                    batch.append(summarize(current))
                    current = []
                    if len(batch) >= 50_000:
                        writer.write_table(pa.Table.from_pylist(batch, schema=OUT_SCHEMA)); batch.clear()
                current_key = key
                current.append(row)
        if current:
            batch.append(summarize(current))
        if batch:
            writer.write_table(pa.Table.from_pylist(batch, schema=OUT_SCHEMA))
    finally:
        writer.close()
        connection.close()

    counts = {}
    keep = {"VUS_TO_P_LP", "VUS_TO_B_LB", "PERSISTENT_VUS", "VUS_TO_CONFLICTING", "STABLE_P_LP", "STABLE_B_LB"}
    frozen_path = out_dir / "objective1_frozen_cohorts.parquet"
    frozen_writer = pq.ParquetWriter(frozen_path, OUT_SCHEMA, compression="zstd")
    variants = frozen_rows = 0
    try:
        for record_batch in pq.ParquetFile(output).iter_batches(batch_size=100_000):
            variants += record_batch.num_rows
            trajectory_index = OUT_SCHEMA.get_field_index("trajectory")
            for value in record_batch.column(trajectory_index).to_pylist():
                counts[value] = counts.get(value, 0) + 1
            table_batch = pa.Table.from_batches([record_batch], schema=OUT_SCHEMA)
            mask = pc.is_in(table_batch["trajectory"], value_set=pa.array(sorted(keep)))
            selected = table_batch.filter(mask)
            if selected.num_rows:
                frozen_writer.write_table(selected)
                frozen_rows += selected.num_rows
    finally:
        frozen_writer.close()
    summary = {"status": "complete", "release_partitions": len(fragments), "variants": variants,
               "frozen_rows": frozen_rows, "trajectory_counts": counts,
               "normalization_boundary": "split ALT + minimal representation; non-ACGTN excluded; reference-context left-alignment audit remains separate"}
    (args.root / "reports/objective1_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
