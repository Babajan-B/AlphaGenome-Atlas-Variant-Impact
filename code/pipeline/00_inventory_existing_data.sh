#!/usr/bin/env bash
set -euo pipefail

ROOT="${AGVUS_ROOT:-/path/to/project_root/AG-VUS}"
SOURCE=/path/to/project_root/KAUBioMED_LLM/data/raw/clinvar
OUT="$ROOT/reports/data_inventory.tsv"

mkdir -p "$(dirname "$OUT")"
printf 'bytes\tmodified\tsha256\tpath\n' > "$OUT"

find "$SOURCE" -type f -print0 | sort -z | while IFS= read -r -d '' path; do
  bytes=$(stat -c '%s' "$path")
  modified=$(stat -c '%y' "$path" | cut -d. -f1)
  checksum=$(sha256sum "$path" | cut -d' ' -f1)
  printf '%s\t%s\t%s\t%s\n' "$bytes" "$modified" "$checksum" "$path" >> "$OUT"
done

echo "Wrote $OUT"

