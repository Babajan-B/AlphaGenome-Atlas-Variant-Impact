#!/usr/bin/env python3
"""Test eligibility and deterministic sampling without API calls or credentials."""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def main():
    with tempfile.TemporaryDirectory(prefix="agvus_obj2_test_") as directory:
        root = Path(directory)
        source = root / "data/processed/longitudinal"
        source.mkdir(parents=True)
        (root / "logs").mkdir()
        rows = []
        for index, (chrom, ref, alt) in enumerate([
            ("1", "A", "C"), ("2", "G", "T"), ("X", "A", "AT"),
            ("Y", "N", "C"), ("MT", "A", "C"), ("1", "A", "A"),
        ]):
            rows.append(dict(canonical_id=str(index), trajectory="PERSISTENT_VUS",
                             chrom=chrom, norm_pos=1, norm_ref=ref, norm_alt=alt,
                             molecular_consequence="", geneinfo="",
                             reference_status="valid"))
        pq.write_table(pa.Table.from_pylist(rows), source / "objective1_frozen_cohorts_ref_validated.parquet")
        (root / "logs/OBJECTIVE1.COMPLETE").write_text(json.dumps({"status": "complete", "frozen_cohort_total": len(rows)}))
        command = [sys.executable, str(Path(__file__).parent / "07_objective2_atlas_pilot.py"),
                   "--root", str(root), "--inventory-only", "--per-cohort", "1"]
        subprocess.run(command, check=True)
        out = root / "annotations/alphagenome/objective2_pilot"
        counts = {r["reason"]: r["n"] for r in json.loads((out / "eligibility_inventory.json").read_text())}
        assert counts == {"eligible_snv": 2, "non_snv": 1, "ambiguous_allele": 1,
                          "unsupported_contig": 1, "no_change": 1}, counts
        first = json.loads((out / "pilot_manifest.json").read_text())
        assert len(first) == 1 and first[0]["eligibility_reason"] == "eligible_snv"
        subprocess.run(command, check=True)
        assert first == json.loads((out / "pilot_manifest.json").read_text())
        print("OBJECTIVE2_INVENTORY_TEST_PASS")


if __name__ == "__main__":
    main()
