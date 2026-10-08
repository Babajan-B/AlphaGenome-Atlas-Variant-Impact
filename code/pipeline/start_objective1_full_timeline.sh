#!/usr/bin/env bash
set -euo pipefail

ROOT=/path/to/project_root/AG-VUS
PY="$ROOT/environment/alphagenome-py310/bin/python"
MONTHLY_PID_FILE="$ROOT/logs/clinvar_monthly_download.pid"
LOG="$ROOT/logs/objective1_full_timeline.log"
PID_FILE="$ROOT/logs/objective1_full_timeline.pid"

if test -f "$PID_FILE" && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
  echo "Objective 1 timeline process already running as PID $(cat "$PID_FILE")"
  exit 0
fi

nohup bash -c '
  set -euo pipefail
  root=/path/to/project_root/AG-VUS
  monthly_pid_file="$root/logs/clinvar_monthly_download.pid"
  if test -f "$monthly_pid_file"; then
    monthly_pid=$(cat "$monthly_pid_file")
    while kill -0 "$monthly_pid" 2>/dev/null; do sleep 300; done
  fi
  exec "$root/environment/alphagenome-py310/bin/python" \
    /path/to/AG-VUS/scripts/03_download_monthly_clinvar.py \
    --selection all
' > "$LOG" 2>&1 < /dev/null &

echo $! > "$PID_FILE"
echo "Objective 1 timeline process started as PID $(cat "$PID_FILE")"
echo "It will wait for the monthly download and then fetch only missing releases."

