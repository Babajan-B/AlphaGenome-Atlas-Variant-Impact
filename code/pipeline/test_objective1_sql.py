#!/usr/bin/env python3
"""Small regression test of SQL trajectories against the original definitions."""

import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def module(filename):
    spec = importlib.util.spec_from_file_location(filename, Path(__file__).parent / filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    parser = module("04_parse_clinvar_releases.py")
    original = module("05_build_longitudinal_histories.py")
    sequences = [
        ["VUS", "P_LP", "P_LP"], ["VUS", "B_LB", "VUS"],
        ["VUS", "VUS", "VUS"], ["VUS", "CONFLICTING", "CONFLICTING"],
        ["P_LP", "P_LP", "P_LP"], ["B_LB", "B_LB", "B_LB"],
        ["VUS", "OTHER", "OTHER"], ["OTHER", "P_LP", "B_LB"],
    ]
    with tempfile.TemporaryDirectory(prefix="agvus_sql_test_") as directory:
        root = Path(directory)
        (root / "reports").mkdir()
        expected = {}
        for index, classes in enumerate(sequences):
            rows = []
            for offset, value in enumerate(classes):
                row = {field.name: "" for field in parser.SCHEMA}
                row.update(
                    release_date=f"2026010{offset + 1}", variant_id=str(index),
                    chrom=str(index % 2 + 1), source_pos=index + 1, source_ref="A", source_alt="C",
                    norm_pos=index + 1, norm_ref="A", norm_alt="C",
                    canonical_id=f"GRCh38_{index % 2 + 1}_{index + 1}_A_C",
                    classification_group=value, review_status="reviewed_by_expert_panel",
                    geneinfo="TEST:1", clnhgvs="test", molecular_consequence="missense_variant",
                    normalization_status="minimal_normalized",
                )
                rows.append(row)
            expected[rows[0]["canonical_id"]] = original.summarize(rows)
            for row in rows:
                output = root / "data/processed/release_observations" / f"release_date={row['release_date']}" / "observations.parquet"
                output.parent.mkdir(parents=True, exist_ok=True)
                existing = pq.ParquetFile(output).read().to_pylist() if output.exists() else []
                pq.write_table(pa.Table.from_pylist(existing + [row], schema=parser.SCHEMA), output)
        subprocess.run([
            sys.executable, str(Path(__file__).parent / "05_build_longitudinal_batched.py"),
            "--root", str(root), "--expected-releases", "3", "--threads", "1",
        ], check=True)
        actual = pq.ParquetFile(root / "data/processed/longitudinal/variant_histories.parquet").read().to_pylist()
        assert len(actual) == len(expected)
        for row in actual:
            for key, value in expected[row["canonical_id"]].items():
                assert row[key] == value, (row["canonical_id"], key, row[key], value)
        assert pq.ParquetFile(root / "data/processed/longitudinal/objective1_frozen_cohorts.parquet").metadata.num_rows == 6
        subprocess.run([
            sys.executable, str(Path(__file__).parent / "05_build_longitudinal_batched.py"),
            "--root", str(root), "--expected-releases", "3", "--threads", "1",
        ], check=True)
        print("SQL_REGRESSION_PASS: 8 trajectories match original definitions")


if __name__ == "__main__":
    main()
