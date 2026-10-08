#!/usr/bin/env bash
set -euo pipefail

ROOT=/path/to/project_root/AG-VUS
PY="$ROOT/environment/alphagenome-py310/bin/python"
LOG="$ROOT/logs/clinvar_monthly_download.log"
PID="$ROOT/logs/clinvar_monthly_download.pid"

if test -f "$PID" && kill -0 "$(cat "$PID")" 2>/dev/null; then
  echo "Download already running as PID $(cat "$PID")"
  exit 0
fi

nohup "$PY" /path/to/AG-VUS/scripts/03_download_monthly_clinvar.py \
  > "$LOG" 2>&1 < /dev/null &
echo $! > "$PID"
echo "Started ClinVar monthly download as PID $(cat "$PID")"
echo "Log: $LOG"

