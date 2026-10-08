#!/usr/bin/env python3
"""Parse one shard of dated ClinVar VCFs into minimal, audited Parquet partitions."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


INFO_RE = re.compile(r"(?:^|;)([A-Za-z0-9_]+)=([^;]*)")


def info_dict(text: str) -> dict[str, str]:
    return {key: value for key, value in INFO_RE.findall(text)}


def harmonize(value: str) -> str:
    text = value.lower().replace("%2c", ",")
    if "conflict" in text:
        return "CONFLICTING"
    has_path = "pathogenic" in text
    has_benign = "benign" in text
    if has_path and has_benign:
        return "CONFLICTING"
    if has_path:
        return "P_LP"
    if has_benign:
        return "B_LB"
    if "uncertain" in text or "vus" in text:
        return "VUS"
    return "OTHER"


def minimal_normalize(pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    ref, alt = ref.upper(), alt.upper()
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


SCHEMA = pa.schema([
    ("release_date", pa.string()), ("variant_id", pa.string()),
    ("chrom", pa.string()), ("source_pos", pa.int64()),
    ("source_ref", pa.string()), ("source_alt", pa.string()),
    ("norm_pos", pa.int64()), ("norm_ref", pa.string()), ("norm_alt", pa.string()),
    ("canonical_id", pa.string()), ("classification_raw", pa.string()),
    ("classification_group", pa.string()), ("review_status", pa.string()),
    ("geneinfo", pa.string()), ("clnhgvs", pa.string()), ("molecular_consequence", pa.string()),
    ("normalization_status", pa.string()),
])


def flush(rows: list[dict], writer: pq.ParquetWriter) -> None:
    if rows:
        writer.write_table(pa.Table.from_pylist(rows, schema=SCHEMA))
        rows.clear()


def parse_release(vcf: Path, date: str, output: Path) -> dict:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".parquet.tmp")
    rows: list[dict] = []
    records = alleles = symbolic = 0
    writer = pq.ParquetWriter(temporary, SCHEMA, compression="zstd")
    try:
        with gzip.open(vcf, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line or line[0] == "#":
                    continue
                fields = line.rstrip("\n").split("\t", 8)
                if len(fields) < 8:
                    continue
                chrom, pos_text, record_id, ref, alt_text, _, _, info_text = fields[:8]
                records += 1
                info = info_dict(info_text)
                for alt in alt_text.split(","):
                    alleles += 1
                    if alt.startswith("<") or alt == "*":
                        symbolic += 1
                        continue
                    pos, norm_ref, norm_alt = minimal_normalize(int(pos_text), ref, alt)
                    status = "minimal_normalized"
                    if not set(norm_ref + norm_alt) <= set("ACGTN"):
                        status = "non_acgtn_allele"
                    variant_id = info.get("CLNVID") or info.get("ALLELEID") or record_id
                    rows.append({
                        "release_date": date, "variant_id": str(variant_id), "chrom": chrom.removeprefix("chr"),
                        "source_pos": int(pos_text), "source_ref": ref, "source_alt": alt,
                        "norm_pos": pos, "norm_ref": norm_ref, "norm_alt": norm_alt,
                        "canonical_id": f"GRCh38_{chrom.removeprefix('chr')}_{pos}_{norm_ref}_{norm_alt}",
                        "classification_raw": info.get("CLNSIG", ""),
                        "classification_group": harmonize(info.get("CLNSIG", "")),
                        "review_status": info.get("CLNREVSTAT", ""), "geneinfo": info.get("GENEINFO", ""),
                        "clnhgvs": info.get("CLNHGVS", ""), "molecular_consequence": info.get("MC", ""),
                        "normalization_status": status,
                    })
                    if len(rows) >= 100_000:
                        flush(rows, writer)
        flush(rows, writer)
    finally:
        writer.close()
    temporary.replace(output)
    return {"release_date": date, "vcf": str(vcf), "output": str(output), "records": records,
            "alleles": alleles, "symbolic_excluded": symbolic, "status": "ok"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("/path/to/project_root/AG-VUS"))
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    args = parser.parse_args()
    manifest = args.root / "reports/clinvar_all_download_status.tsv"
    with manifest.open() as handle:
        releases = [row for row in csv.DictReader(handle, delimiter="\t") if row["status"] == "verified"]
    if len(releases) != 420:
        raise RuntimeError(f"Expected 420 verified releases, found {len(releases)}")
    selected = releases[args.shard_index::args.shard_count]
    report_path = args.root / "reports/objective1_parse" / f"shard_{args.shard_index:02d}.jsonl"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("a", encoding="utf-8") as report:
        for row in selected:
            date = row["release_date"]
            output = args.root / "data/processed/release_observations" / f"release_date={date}" / "observations.parquet"
            if output.exists():
                result = {"release_date": date, "output": str(output), "status": "already_complete"}
            else:
                result = parse_release(Path(row["vcf"]), date, output)
            report.write(json.dumps(result) + "\n")
            report.flush()
            print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
