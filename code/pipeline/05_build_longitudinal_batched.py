#!/usr/bin/env python3
"""Run memory-capped chromosome rebuilds serially and merge their checkpoints."""

import argparse
import json
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from urllib.parse import quote

import duckdb
import pyarrow.parquet as pq


def merge(parts, destination):
    temporary = destination.with_suffix(".merge.tmp.parquet")
    writer = None
    count = 0
    try:
        for part in parts:
            source = pq.ParquetFile(part)
            if writer is None:
                writer = pq.ParquetWriter(temporary, source.schema_arrow, compression="zstd")
            for batch in source.iter_batches(batch_size=50000):
                writer.write_batch(batch)
                count += batch.num_rows
    finally:
        if writer is not None:
            writer.close()
    if writer is None:
        raise RuntimeError("No chromosome parts were produced")
    if pq.ParquetFile(temporary).metadata.num_rows != count:
        raise RuntimeError("Merged Parquet row count mismatch")
    temporary.replace(destination)
    return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("/path/to/project_root/AG-VUS"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-limit", default="48GB")
    parser.add_argument("--expected-releases", type=int, default=420)
    args = parser.parse_args()
    start = time.monotonic()
    source = args.root / "data/processed/release_observations"
    fragments = list(source.glob("release_date=*/observations.parquet"))
    if len(fragments) != args.expected_releases:
        raise RuntimeError("Release partition count mismatch")
    output = args.root / "data/processed/longitudinal"
    output.mkdir(parents=True, exist_ok=True)
    (args.root / "reports").mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute("SET memory_limit='1GB'")
    con.execute("SET threads=2")
    chromosomes = [r[0] for r in con.execute(
        "SELECT DISTINCT chrom FROM read_parquet(?, hive_partitioning=true) ORDER BY chrom",
        [str(source / "release_date=*" / "observations.parquet")],
    ).fetchall()]
    con.close()
    if not chromosomes or any(chrom is None or not str(chrom) for chrom in chromosomes):
        raise RuntimeError("Empty or invalid chromosome inventory")
    parts = []
    counts = Counter()
    for index, chrom in enumerate(chromosomes, 1):
        chrom = str(chrom)
        before = time.monotonic()
        print(f"chromosome={chrom} progress={index}/{len(chromosomes)} started", flush=True)
        subprocess.run([
            sys.executable, str(Path(__file__).parent / "05_build_longitudinal_histories_sql.py"),
            "--root", str(args.root), "--chrom", chrom,
            "--threads", str(args.threads), "--memory-limit", args.memory_limit,
            "--expected-releases", str(args.expected_releases),
        ], check=True)
        part = output / "chromosome_parts_v2" / quote(chrom, safe="")
        summary = json.loads((part / "summary.json").read_text())
        if summary.get("status") != "complete":
            raise RuntimeError(f"Chromosome {chrom} lacks a complete checkpoint")
        counts.update(summary["trajectory_counts"])
        parts.append(part)
        # These are this pipeline's disposable intermediate database files.
        for filename in ("objective1_aggregation_sql.duckdb", "objective1_aggregation_sql.duckdb.wal"):
            (part / filename).unlink(missing_ok=True)
        print(f"chromosome={chrom} completed seconds={time.monotonic() - before:.1f}", flush=True)
    variants = merge([p / "variant_histories.parquet" for p in parts], output / "variant_histories.parquet")
    frozen = merge([p / "objective1_frozen_cohorts.parquet" for p in parts], output / "objective1_frozen_cohorts.parquet")
    if variants != sum(counts.values()):
        raise RuntimeError("Chromosome trajectory counts do not reconcile with merged rows")
    summary = {
        "status": "complete", "engine": "serial chromosome-batched DuckDB SQL",
        "pipeline_version": 2, "release_partitions": len(fragments),
        "chromosomes": chromosomes, "variants": variants, "frozen_rows": frozen,
        "trajectory_counts": dict(counts), "elapsed_seconds": time.monotonic() - start,
        "memory_limit_per_worker": args.memory_limit, "worker_threads": args.threads,
        "normalization_boundary": "ALT splitting and parsimonious normalization; GRCh38 REF validation follows; full left alignment is not claimed",
    }
    destination = args.root / "reports/objective1_summary.json"
    temporary = destination.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(summary, indent=2) + "\n")
    temporary.replace(destination)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
