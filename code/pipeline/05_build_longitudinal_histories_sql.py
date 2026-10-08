#!/usr/bin/env python3
"""Build Objective 1 longitudinal cohorts with DuckDB SQL.

This replaces the original Python row loop, which did not complete at full
ClinVar scale. Intermediate tables live in a restartable DuckDB database.
"""

from __future__ import annotations

import argparse
import json
import hashlib
from pathlib import Path
from urllib.parse import quote

import duckdb
import pyarrow.parquet as pq


KEEP = (
    "VUS_TO_P_LP",
    "VUS_TO_B_LB",
    "PERSISTENT_VUS",
    "VUS_TO_CONFLICTING",
    "STABLE_P_LP",
    "STABLE_B_LB",
)


def sql_string(path: Path) -> str:
    return str(path).replace("'", "''")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("/path/to/project_root/AG-VUS"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-limit", default="48GB")
    parser.add_argument("--chrom")
    parser.add_argument("--expected-releases", type=int, default=420)
    args = parser.parse_args()

    root = args.root
    source = root / "data/processed/release_observations"
    fragments = list(source.glob("release_date=*/observations.parquet"))
    if len(fragments) != args.expected_releases:
        raise RuntimeError(f"Expected {args.expected_releases} release partitions, found {len(fragments)}")

    out_dir = root / "data/processed/longitudinal"
    if args.chrom is not None:
        out_dir = out_dir / "chromosome_parts_v2" / quote(args.chrom, safe="")
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = out_dir / "spill"
    temp_dir.mkdir(parents=True, exist_ok=True)
    database = out_dir / "objective1_aggregation_sql.duckdb"
    histories_tmp = out_dir / "variant_histories.sql.tmp.parquet"
    histories = out_dir / "variant_histories.parquet"
    frozen_tmp = out_dir / "objective1_frozen_cohorts.sql.tmp.parquet"
    frozen = out_dir / "objective1_frozen_cohorts.parquet"
    summary_path = root / "reports/objective1_summary.json"
    if args.chrom is not None:
        summary_path = out_dir / "summary.json"
    fingerprint = hashlib.sha256("\n".join(
        f"{p.parent.name}:{p.stat().st_size}:{p.stat().st_mtime_ns}" for p in sorted(fragments)
    ).encode()).hexdigest()
    if args.chrom is not None and summary_path.exists() and histories.exists() and frozen.exists():
        previous = json.loads(summary_path.read_text())
        if previous.get("input_fingerprint") == fingerprint and previous.get("pipeline_version") == 2:
            if (pq.ParquetFile(histories).metadata.num_rows == previous["variants"] and
                    pq.ParquetFile(frozen).metadata.num_rows == previous["frozen_rows"]):
                print(f"chrom={args.chrom} checkpoint=reused", flush=True)
                return

    for path in (histories_tmp, frozen_tmp):
        path.unlink(missing_ok=True)

    con = duckdb.connect(str(database))
    con.execute(f"PRAGMA threads={args.threads}")
    con.execute("SET memory_limit = ?", [args.memory_limit])
    con.execute(f"PRAGMA temp_directory='{sql_string(temp_dir)}'")
    con.execute("PRAGMA preserve_insertion_order=false")
    parquet_glob = sql_string(source / "release_date=*" / "observations.parquet")
    chrom_filter = "" if args.chrom is None else " AND chrom = '" + args.chrom.replace("'", "''") + "'"

    # Collapse any duplicate representations within the same dated release so
    # every canonical variant has one deterministic state per release.
    print("stage=release_state", flush=True)
    con.execute(
        f"""
        CREATE OR REPLACE TABLE release_state AS
        SELECT
            canonical_id,
            CAST(release_date AS VARCHAR) AS release_date,
            arg_max(variant_id, variant_id) AS variant_id,
            arg_max(geneinfo, length(geneinfo)) AS geneinfo,
            arg_max(clnhgvs, length(clnhgvs)) AS clnhgvs,
            arg_max(molecular_consequence, length(molecular_consequence)) AS molecular_consequence,
            chrom,
            norm_pos,
            norm_ref,
            norm_alt,
            arg_max(source_pos, variant_id) AS source_pos,
            arg_max(source_ref, variant_id) AS source_ref,
            arg_max(source_alt, variant_id) AS source_alt,
            CASE
                WHEN count(DISTINCT classification_group) > 1 THEN 'CONFLICTING'
                ELSE min(classification_group)
            END AS classification_group,
            string_agg(DISTINCT review_status, '|') AS review_status
        FROM read_parquet('{parquet_glob}', hive_partitioning=true)
        WHERE normalization_status <> 'non_acgtn_allele' {chrom_filter}
        GROUP BY canonical_id, release_date, chrom, norm_pos, norm_ref, norm_alt
        """
    )
    print("stage=ordered_state", flush=True)
    con.execute(
        """
        CREATE OR REPLACE TABLE ordered_state AS
        SELECT *,
            lag(classification_group) OVER (
                PARTITION BY canonical_id ORDER BY release_date
            ) AS previous_class
        FROM release_state
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE variant_dates AS
        SELECT canonical_id,
            min(release_date) FILTER (WHERE classification_group = 'VUS') AS first_vus_date,
            max(release_date) FILTER (WHERE classification_group = 'VUS') AS last_vus_date
        FROM release_state
        GROUP BY canonical_id
        """
    )
    con.execute(
        """
        CREATE OR REPLACE TABLE first_reclassification AS
        SELECT s.canonical_id,
            arg_min(s.release_date, s.release_date) AS first_reclassification_date,
            arg_min(s.classification_group, s.release_date) AS reclassification_group,
            arg_min(s.review_status, s.release_date) AS reclassification_review_status
        FROM release_state s
        JOIN variant_dates d USING (canonical_id)
        WHERE d.first_vus_date IS NOT NULL
          AND s.release_date > d.first_vus_date
          AND s.classification_group IN ('P_LP', 'B_LB', 'CONFLICTING')
        GROUP BY s.canonical_id
        """
    )

    print("stage=history_export", flush=True)
    histories_query = """
        WITH aggregate_state AS (
            SELECT
                canonical_id,
                arg_min(variant_id, release_date) AS variant_id,
                arg_max(geneinfo, release_date) AS geneinfo,
                arg_max(clnhgvs, release_date) AS clnhgvs,
                arg_max(molecular_consequence, release_date) AS molecular_consequence,
                arg_min(chrom, release_date) AS chrom,
                arg_min(norm_pos, release_date) AS norm_pos,
                arg_min(norm_ref, release_date) AS norm_ref,
                arg_min(norm_alt, release_date) AS norm_alt,
                min(release_date) AS first_observed,
                max(release_date) AS last_observed,
                arg_min(classification_group, release_date) AS first_class,
                arg_max(classification_group, release_date) AS last_class,
                count(*)::INTEGER AS observed_releases,
                sum(CASE WHEN previous_class IS NOT NULL AND previous_class <> classification_group THEN 1 ELSE 0 END)::INTEGER AS class_changes,
                greatest(sum(CASE WHEN previous_class IS NOT NULL AND previous_class <> classification_group THEN 1 ELSE 0 END) - 1, 0)::INTEGER AS reversions,
                count(DISTINCT struct_pack(pos := source_pos, ref := source_ref, alt := source_alt))::INTEGER AS representation_count,
                count(*) FILTER (WHERE classification_group <> 'P_LP') = 0 AS always_p_lp,
                count(*) FILTER (WHERE classification_group <> 'B_LB') = 0 AS always_b_lb
            FROM ordered_state
            GROUP BY canonical_id
        )
        SELECT
            a.canonical_id, a.variant_id, a.geneinfo, a.clnhgvs, a.molecular_consequence,
            a.chrom, a.norm_pos, a.norm_ref, a.norm_alt,
            a.first_observed, a.last_observed, a.first_class, a.last_class,
            d.first_vus_date, d.last_vus_date, r.first_reclassification_date,
            CASE
                WHEN d.first_vus_date IS NOT NULL AND r.reclassification_group IS NOT NULL
                    THEN 'VUS_TO_' || r.reclassification_group
                WHEN d.first_vus_date IS NOT NULL AND a.last_class = 'VUS' THEN 'PERSISTENT_VUS'
                WHEN d.first_vus_date IS NOT NULL THEN 'VUS_OTHER_UNRESOLVED'
                WHEN a.always_p_lp THEN 'STABLE_P_LP'
                WHEN a.always_b_lb THEN 'STABLE_B_LB'
                ELSE 'OTHER'
            END AS trajectory,
            a.observed_releases, a.class_changes, a.reversions, a.representation_count,
            coalesce(
                contains(lower(r.reclassification_review_status), 'expert_panel') OR
                contains(lower(r.reclassification_review_status), 'practice_guideline'),
                false
            ) AS high_confidence
        FROM aggregate_state a
        JOIN variant_dates d USING (canonical_id)
        LEFT JOIN first_reclassification r USING (canonical_id)
    """
    con.execute(
        f"COPY ({histories_query}) TO '{sql_string(histories_tmp)}' "
        "(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)"
    )
    histories_tmp.replace(histories)

    print("stage=cohort_export", flush=True)
    keep_sql = ", ".join(f"'{value}'" for value in KEEP)
    con.execute(
        f"COPY (SELECT * FROM read_parquet('{sql_string(histories)}') "
        f"WHERE trajectory IN ({keep_sql})) TO '{sql_string(frozen_tmp)}' "
        "(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 100000)"
    )
    frozen_tmp.replace(frozen)

    trajectory_rows = con.execute(
        f"SELECT trajectory, count(*) FROM read_parquet('{sql_string(histories)}') "
        "GROUP BY trajectory ORDER BY trajectory"
    ).fetchall()
    summary = {
        "status": "complete",
        "engine": "DuckDB SQL",
        "pipeline_version": 2,
        "input_fingerprint": fingerprint,
        "chromosome": args.chrom,
        "release_partitions": len(fragments),
        "variants": pq.ParquetFile(histories).metadata.num_rows,
        "frozen_rows": pq.ParquetFile(frozen).metadata.num_rows,
        "trajectory_counts": {name: count for name, count in trajectory_rows},
        "normalization_boundary": (
            "split ALT + minimal representation; non-ACGTN excluded; "
            "GRCh38 reference validation follows as a separate gated job"
        ),
    }
    temporary_summary = summary_path.with_suffix(".json.tmp")
    temporary_summary.write_text(json.dumps(summary, indent=2) + "\n")
    temporary_summary.replace(summary_path)
    con.close()
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
