#!/bin/bash
set -euo pipefail
ROOT=/path/to/project_root/AG-VUS
OUT="$ROOT/annotations/alphagenome/objective2_pilot"
exec >"$ROOT/logs/objective2_login_queries.log" 2>&1
echo "started=$(date --iso-8601=seconds) host=$(hostname)"
trap 'echo "worker_exit=$? ended=$(date --iso-8601=seconds)"' EXIT
set -a
source /path/to/AG-VUS/config/alphagenome.env
set +a
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
flock -n "$OUT/query_worker.lock" timeout 18000 \
  "$ROOT/environment/alphagenome-py310/bin/python" -u \
  /path/to/AG-VUS/scripts/07_objective2_atlas_pilot.py --root "$ROOT" --use-existing-inventory
