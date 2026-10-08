#!/usr/bin/env python3
"""Audit trusted Atlas pilot artifacts and export raw AVI matrix entries."""
import hashlib
import json
import pickle
from collections import Counter
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq


def main():
    root = Path("/path/to/project_root/AG-VUS")
    output = root / "annotations/alphagenome/objective2_pilot"
    summary = json.loads((output / "pilot_summary.json").read_text())
    manifest = json.loads((output / "pilot_manifest.json").read_text())
    if summary["status"] != "pilot_complete" or summary["pilot_size"] != len(manifest):
        raise RuntimeError("Pilot has not completed successfully")
    counts = Counter()
    scores = []
    for row in manifest:
        token = hashlib.sha256(row["canonical_id"].encode()).hexdigest()[:24]
        record = json.loads((output / "variant_results" / f"{token}.json").read_text())
        if record["canonical_id"] != row["canonical_id"] or record["status"] != "annotated":
            raise RuntimeError("Missing annotated pilot row")
        for scorer in ("AVI_SCORE", "AVI_SCORE_MODEL_FEATURES", "AVI_SCORE_FEATURE_IMPORTANCE"):
            file = Path(record["response_files"][scorer])
            if file.parent != output / "variant_results" or file.name != f"{token}.{scorer}.pickle":
                raise RuntimeError("Unexpected artifact path")
            # Load only our own API-generated trusted files, never user uploads.
            with file.open("rb") as handle:
                data = pickle.load(handle)
            if list(data.shape) != record["scorer_shapes"][scorer]:
                raise RuntimeError("Response shape mismatch")
            counts[scorer] += 1
            if scorer == "AVI_SCORE":
                values = data.X.toarray() if hasattr(data.X, "toarray") else np.asarray(data.X)
                if int(np.isfinite(values).sum()) != record["finite_avi_values"]:
                    raise RuntimeError("AVI finite-value count mismatch")
                for i, j in np.ndindex(values.shape):
                    value = float(values[i, j])
                    scores.append({
                        "canonical_id": row["canonical_id"], "trajectory": row["trajectory"],
                        "scorer": scorer, "obs_index": i, "feature_index": j,
                        "obs_id": str(data.obs_names[i]), "feature_id": str(data.var_names[j]),
                        "value": value if np.isfinite(value) else None,
                        "feature_metadata": json.dumps({str(k): repr(v) for k, v in data.var.iloc[j].items()}),
                        "queried_at_utc": record["queried_at_utc"],
                        "alphagenome_version": record["alphagenome_version"],
                    })
    destination = output / "pilot_avi_scores.parquet"
    temporary = destination.with_suffix(".tmp.parquet")
    pq.write_table(pa.Table.from_pylist(scores), temporary, compression="zstd")
    temporary.replace(destination)
    audit = {"status": "verified", "variants": len(manifest), "scorer_artifact_counts": dict(counts),
             "avi_matrix_entries": len(scores), "avi_table": str(destination),
             "boundary": "Raw AVI matrix entries; no averaging, calibration claim or clinical classification"}
    (output / "pilot_artifact_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    (root / "logs/OBJECTIVE2_ATLAS_PILOT.COMPLETE").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
