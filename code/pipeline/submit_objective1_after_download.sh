#!/usr/bin/env bash
set -euo pipefail
ROOT=/path/to/project_root/AG-VUS
DOWNLOAD_PID_FILE="$ROOT/logs/objective1_full_timeline.pid"
LOG="$ROOT/logs/objective1_submit_after_download.log"
exec >>"$LOG" 2>&1

if [[ -f "$DOWNLOAD_PID_FILE" ]]; then
  pid=$(cat "$DOWNLOAD_PID_FILE")
  while kill -0 "$pid" 2>/dev/null; do sleep 300; done
fi

read -r verified failed < <("$ROOT/environment/alphagenome-py310/bin/python" - "$ROOT/reports/clinvar_all_download_status.tsv" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
verified = sum(row.get("status", "").strip() == "verified" for row in rows)
print(verified, len(rows) - verified)
PY
)
if [[ "$verified" -ne 420 || "$failed" -ne 0 ]]; then
  echo "Download gate failed: verified=$verified failed=$failed"
  exit 1
fi

until "$ROOT/environment/alphagenome-py310/bin/python" -c 'import duckdb' 2>/dev/null; do
  echo "Waiting for DuckDB dependency installation"
  sleep 60
done
until [[ -f "$ROOT/logs/GRCH38_REFERENCE.READY" ]]; do
  echo "Waiting for indexed GRCh38 reference"
  sleep 60
done

job_ids=()
for shard in $(seq 0 7); do
  job_ids+=("$(/opt/pbs/bin/qsub -v SHARD_INDEX=$shard,SHARD_COUNT=8 /path/to/AG-VUS/scripts/objective1_parse_shard.pbs)")
done
dependency=$(IFS=:; echo "${job_ids[*]}")
aggregate=$(/opt/pbs/bin/qsub -W depend=afterok:"$dependency" /path/to/AG-VUS/scripts/objective1_aggregate.pbs)
validate=$(/opt/pbs/bin/qsub -W depend=afterok:"$aggregate" /path/to/AG-VUS/scripts/objective1_validate.pbs)
printf 'parse_jobs=%s\naggregate_job=%s\nvalidate_job=%s\n' "$dependency" "$aggregate" "$validate" | tee "$ROOT/reports/objective1_job_chain.txt"
