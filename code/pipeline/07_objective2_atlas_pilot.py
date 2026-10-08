#!/usr/bin/env python3
"""Inventory frozen-cohort eligibility and checkpoint a bounded Atlas pilot."""

import argparse
import hashlib
import importlib.metadata
import json
import os
import pickle
import signal
import time
from contextlib import contextmanager
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import duckdb
import numpy as np
from alphagenome.atlas import atlas
from alphagenome.data import genome


SCORERS = ["AVI_SCORE", "AVI_SCORE_MODEL_FEATURES", "AVI_SCORE_FEATURE_IMPORTANCE"]


def now():
    return datetime.now(timezone.utc).isoformat()


@contextmanager
def deadline(seconds):
    def expired(signum, frame):
        raise TimeoutError("Atlas request deadline exceeded")
    previous = signal.signal(signal.SIGALRM, expired)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def atomic_json(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("/path/to/project_root/AG-VUS"))
    parser.add_argument("--per-cohort", type=int, default=20)
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--use-existing-inventory", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.per_cohort <= 20:
        raise ValueError("Pilot limit is 1-20 variants per cohort")
    root = args.root
    qc = json.loads((root / "logs/OBJECTIVE1.COMPLETE").read_text())
    if qc.get("status") != "complete":
        raise RuntimeError("Objective 1 completion gate failed")
    source = root / "data/processed/longitudinal/objective1_frozen_cohorts_ref_validated.parquet"
    output = root / "annotations/alphagenome/objective2_pilot"
    output.mkdir(parents=True, exist_ok=True)
    if args.use_existing_inventory:
        inventory = json.loads((output / "eligibility_inventory.json").read_text())
        rows = json.loads((output / "pilot_manifest.json").read_text())
        if sum(r["n"] for r in inventory) != qc["frozen_cohort_total"]:
            raise RuntimeError("Saved inventory does not reconcile with final QC")
        if not rows or len(rows) > 120 or any(r["eligibility_reason"] != "eligible_snv" for r in rows):
            raise RuntimeError("Invalid saved pilot manifest")
        run_queries(rows, source, output, qc)
        return
    con = duckdb.connect()
    con.execute("SET memory_limit='2GB'")
    con.execute("SET threads=1")
    spill = output / "spill"
    spill.mkdir(exist_ok=True)
    con.execute("SET temp_directory = ?", [str(spill)])
    source_sql = str(source).replace("'", "''")
    con.execute(f"CREATE VIEW cohorts AS SELECT * FROM read_parquet('{source_sql}')")
    con.execute("""
        CREATE VIEW eligibility AS
        SELECT canonical_id, trajectory, chrom, norm_pos, norm_ref, norm_alt,
               molecular_consequence, geneinfo,
               CASE
                 WHEN reference_status <> 'valid' OR reference_status IS NULL THEN 'reference_invalid'
                 WHEN chrom NOT IN ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y') THEN 'unsupported_contig'
                 WHEN length(norm_ref) <> 1 OR length(norm_alt) <> 1 THEN 'non_snv'
                 WHEN norm_ref NOT IN ('A','C','G','T') OR norm_alt NOT IN ('A','C','G','T') THEN 'ambiguous_allele'
                 WHEN norm_ref = norm_alt THEN 'no_change'
                 ELSE 'eligible_snv'
               END AS eligibility_reason
        FROM cohorts
    """)
    inventory = [dict(zip(["trajectory", "reason", "n"], row)) for row in con.execute(
        "SELECT trajectory, eligibility_reason, count(*) FROM eligibility GROUP BY ALL ORDER BY 1,2"
    ).fetchall()]
    if sum(r["n"] for r in inventory) != qc["frozen_cohort_total"]:
        raise RuntimeError("Eligibility inventory does not reconcile with Objective 1 QC")
    eligibility_tmp = output / "eligibility.tmp.parquet"
    safe_path = str(eligibility_tmp).replace("'", "''")
    con.execute(f"COPY eligibility TO '{safe_path}' (FORMAT PARQUET, COMPRESSION ZSTD)")
    eligibility_tmp.replace(output / "eligibility.parquet")
    selected = con.execute("""
        SELECT * EXCLUDE (rank) FROM (
          SELECT *, row_number() OVER (
              PARTITION BY trajectory ORDER BY md5(canonical_id), canonical_id
          ) AS rank
          FROM eligibility WHERE eligibility_reason = 'eligible_snv'
        ) WHERE rank <= ? ORDER BY trajectory, canonical_id
    """, [args.per_cohort])
    columns = [r[0] for r in selected.description]
    rows = [dict(zip(columns, r)) for r in selected.fetchall()]
    con.close()
    atomic_json(output / "eligibility_inventory.json", inventory)
    atomic_json(output / "pilot_manifest.json", rows)
    print(f"inventory_rows={sum(r['n'] for r in inventory)} pilot_variants={len(rows)}", flush=True)
    if args.inventory_only:
        return
    run_queries(rows, source, output, qc)


def run_queries(rows, source, output, qc):

    api_key = os.environ.get("ALPHA_GENOME_API_KEY")
    if not api_key:
        raise RuntimeError("API key is not configured")
    client = atlas.create(api_key)
    print("api_stage=scorer_metadata", flush=True)
    with deadline(60):
        metadata = client.scorer_metadata()
    missing = set(SCORERS) - set(metadata)
    atomic_json(output / "scorer_inventory.json", sorted(metadata))
    atomic_json(output / "scorer_metadata.json", {name: repr(metadata[name]) for name in SCORERS if name in metadata})
    if missing:
        raise RuntimeError(f"Missing required Atlas scorers: {sorted(missing)}")
    records = output / "variant_results"
    records.mkdir(exist_ok=True)
    package = importlib.metadata.version("alphagenome")
    results = []
    for index, row in enumerate(rows, 1):
        ident = row["canonical_id"]
        token = hashlib.sha256(ident.encode()).hexdigest()[:24]
        record_path = records / f"{token}.json"
        if record_path.exists():
            old = json.loads(record_path.read_text())
            if old.get("canonical_id") == ident and old.get("status") in {"annotated", "no_finite_avi"}:
                results.append(old)
                print(f"variant={index}/{len(rows)} checkpoint=reused", flush=True)
                continue
        result_record = {
            "canonical_id": ident, "trajectory": row["trajectory"],
            "queried_at_utc": now(), "alphagenome_version": package,
            "requested_scorers": SCORERS, "status": "api_error", "attempts": 0,
            "response_files": {}, "scorer_shapes": {}, "finite_avi_values": 0,
        }
        variant = genome.Variant(
            chromosome="chr" + row["chrom"], position=int(row["norm_pos"]),
            reference_bases=row["norm_ref"], alternate_bases=row["norm_alt"],
        )
        response = None
        for attempt in range(1, 4):
            result_record["attempts"] = attempt
            try:
                with deadline(45):
                    response = client.query_variant(variant, requested_scorers=SCORERS)
                break
            except Exception as exc:
                # Do not log exception messages that might contain credentials.
                result_record["error_type"] = type(exc).__name__
                if attempt < 3:
                    time.sleep(2 ** attempt)
        if response is not None:
            try:
                for scorer, data in response.items():
                    # Atlas AnnData metadata can contain Python objects that
                    # H5AD cannot serialize. Preserve the exact trusted response.
                    destination = records / f"{token}.{scorer}.pickle"
                    temporary = destination.with_suffix(".pickle.tmp")
                    with temporary.open("wb") as handle:
                        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
                    temporary.replace(destination)
                    result_record["response_files"][scorer] = str(destination)
                    result_record["scorer_shapes"][scorer] = list(data.shape)
                    result_record.setdefault("scorer_schema", {})[scorer] = {
                        "obs_columns": list(map(str, data.obs.columns)),
                        "var_columns": list(map(str, data.var.columns)),
                        "layers": list(map(str, data.layers.keys())),
                        "uns_keys": list(map(str, data.uns.keys())),
                    }
                    if scorer == "AVI_SCORE" and data.X is not None:
                        values = data.X.toarray() if hasattr(data.X, "toarray") else np.asarray(data.X)
                        result_record["finite_avi_values"] = int(np.isfinite(values).sum())
                result_record["status"] = "annotated" if result_record["finite_avi_values"] else "no_finite_avi"
                result_record.pop("error_type", None)
            except Exception as exc:
                result_record["status"] = "serialization_error"
                result_record["error_type"] = type(exc).__name__
        atomic_json(record_path, result_record)
        results.append(result_record)
        print(f"variant={index}/{len(rows)} status={result_record['status']}", flush=True)
        time.sleep(0.5)
    counts = Counter(r["status"] for r in results)
    summary = {
        "status": "pilot_complete" if not counts["api_error"] and not counts["serialization_error"] else "pilot_complete_with_errors",
        "completed_at_utc": now(), "alphagenome_version": package,
        "input_cohort": str(source), "reference_qc": qc["frozen_cohort_total"],
        "pilot_size": len(rows), "status_counts": dict(counts),
        "by_cohort": {cohort: dict(Counter(r["status"] for r in results if r["trajectory"] == cohort))
                      for cohort in sorted({r["trajectory"] for r in rows})},
        "scientific_boundary": "Coverage pilot only; no performance claims or calibrated probability interpretation; AVI and AlphaMissense are not independent evidence",
    }
    atomic_json(output / "pilot_summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)
    if counts["api_error"] or counts["serialization_error"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
