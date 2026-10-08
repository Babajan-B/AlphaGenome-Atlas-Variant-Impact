#!/usr/bin/env bash
set -euo pipefail

ROOT=/path/to/project_root/AG-VUS
ENV_FILE=/path/to/AG-VUS/config/alphagenome.env
PY="$ROOT/environment/alphagenome-py310/bin/python"

export AGVUS_ROOT="$ROOT"
source "$ENV_FILE"

/path/to/AG-VUS/scripts/00_inventory_existing_data.sh
/path/to/AG-VUS/scripts/02_link_existing_clinvar.sh
"$PY" /path/to/AG-VUS/scripts/01_alphagenome_atlas_smoke.py \
  > "$ROOT/logs/alphagenome_atlas_smoke.log" 2>&1

echo "Phase 0 complete"

