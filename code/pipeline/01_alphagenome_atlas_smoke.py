#!/usr/bin/env python3
"""Validate Atlas authentication and make one small, reproducible query."""

from __future__ import annotations

import json
import os
import platform
from datetime import datetime, timezone
from pathlib import Path

import alphagenome
from alphagenome.atlas import atlas
from alphagenome.data import genome


def main() -> None:
    api_key = os.environ.get("ALPHA_GENOME_API_KEY")
    if not api_key:
        raise SystemExit("ALPHA_GENOME_API_KEY is not set")

    root = Path(os.environ.get("AGVUS_ROOT", "/path/to/project_root/AG-VUS"))
    out_dir = root / "annotations" / "alphagenome" / "smoke_test"
    out_dir.mkdir(parents=True, exist_ok=True)

    client = atlas.create(api_key)
    metadata = client.scorer_metadata()
    scorer_names = sorted(metadata)
    if not scorer_names:
        raise RuntimeError("Atlas returned no scorer metadata")

    preferred = next(
        (name for name in scorer_names if "avi" in name.lower()), scorer_names[0]
    )
    variant = genome.Variant(
        chromosome="chr22",
        position=36_201_698,
        reference_bases="A",
        alternate_bases="C",
    )
    result = client.query_variant(variant, requested_scorers=[preferred])

    summary = {
        "queried_at_utc": datetime.now(timezone.utc).isoformat(),
        "python": platform.python_version(),
        "alphagenome": getattr(alphagenome, "__version__", "unknown"),
        "scorer_count": len(scorer_names),
        "requested_scorer": preferred,
        "variant": "chr22:36201698:A:C",
        "returned_keys": sorted(result),
        "status": "passed",
    }
    (out_dir / "atlas_smoke_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (out_dir / "scorer_names.txt").write_text(
        "\n".join(scorer_names) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

