#!/usr/bin/env python3
"""Validate normalized REF alleles against GRCh38 and freeze final Objective 1 cohorts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pysam


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("/path/to/project_root/AG-VUS"))
    args = parser.parse_args()
    source = args.root / "data/processed/longitudinal/variant_histories.parquet"
    reference = args.root / "data/reference/GRCh38_no_alt_analysis_set.fna"
    output_dir = args.root / "data/processed/longitudinal"
    audit_path = output_dir / "reference_validation_audit.parquet"
    final_path = output_dir / "objective1_frozen_cohorts_ref_validated.parquet"
    fasta = pysam.FastaFile(str(reference))
    audit_schema = pa.schema([
        ("canonical_id", pa.string()), ("reference_status", pa.string()),
        ("observed_ref", pa.string()), ("expected_ref", pa.string()),
    ])
    history_file = pq.ParquetFile(source)
    history_schema = history_file.schema_arrow.append(pa.field("reference_status", pa.string()))
    audit_writer = pq.ParquetWriter(audit_path, audit_schema, compression="zstd")
    final_writer = pq.ParquetWriter(final_path, history_schema, compression="zstd")
    counts = {"valid": 0, "mismatch": 0, "contig_missing": 0, "fetch_error": 0}
    cohort_trajectories = {
        "VUS_TO_P_LP",
        "VUS_TO_B_LB",
        "PERSISTENT_VUS",
        "VUS_TO_CONFLICTING",
        "STABLE_P_LP",
        "STABLE_B_LB",
    }
    cohort_counts = {trajectory: 0 for trajectory in sorted(cohort_trajectories)}
    try:
        for batch in history_file.iter_batches(batch_size=50_000):
            table = pa.Table.from_batches([batch])
            audits, statuses = [], []
            for row in table.select(["canonical_id", "chrom", "norm_pos", "norm_ref"]).to_pylist():
                contig = row["chrom"] if row["chrom"].startswith("chr") else f"chr{row['chrom']}"
                try:
                    expected = fasta.fetch(contig, row["norm_pos"] - 1, row["norm_pos"] - 1 + len(row["norm_ref"])).upper()
                    status = "valid" if expected == row["norm_ref"].upper() else "mismatch"
                except KeyError:
                    expected, status = "", "contig_missing"
                except Exception:
                    expected, status = "", "fetch_error"
                counts[status] += 1
                statuses.append(status)
                audits.append({"canonical_id": row["canonical_id"], "reference_status": status,
                               "observed_ref": row["norm_ref"], "expected_ref": expected})
            audit_writer.write_table(pa.Table.from_pylist(audits, schema=audit_schema))
            with_status = table.append_column("reference_status", pa.array(statuses))
            trajectories = table.column("trajectory").to_pylist()
            keep = [
                status == "valid" and trajectory in cohort_trajectories
                for status, trajectory in zip(statuses, trajectories)
            ]
            selected = with_status.filter(pa.array(keep))
            if selected.num_rows:
                final_writer.write_table(selected)
                for trajectory in selected.column("trajectory").to_pylist():
                    cohort_counts[trajectory] += 1
    finally:
        audit_writer.close(); final_writer.close(); fasta.close()
    summary = {
        "status": "complete", "reference": str(reference), "reference_validation": counts,
        "retention_rate": round(counts["valid"] / max(1, sum(counts.values())), 6),
        "frozen_cohort_counts": cohort_counts,
        "frozen_cohort_total": sum(cohort_counts.values()),
        "final_cohorts": str(final_path), "audit": str(audit_path),
        "scientific_boundary": "ALT splitting, parsimonious normalization, and GRCh38 REF validation completed; repeat-context left alignment is not claimed beyond ClinVar source normalization.",
    }
    (args.root / "reports/objective1_final_qc.json").write_text(json.dumps(summary, indent=2))
    (args.root / "logs/OBJECTIVE1.COMPLETE").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
