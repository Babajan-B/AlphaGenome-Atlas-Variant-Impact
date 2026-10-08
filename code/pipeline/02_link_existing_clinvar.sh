#!/usr/bin/env bash
set -euo pipefail

ROOT="${AGVUS_ROOT:-/path/to/project_root/AG-VUS}"
SOURCE=/path/to/project_root/KAUBioMED_LLM/data/raw/clinvar
TARGET="$ROOT/data/raw/clinvar/existing"

mkdir -p "$TARGET"
for path in \
  "$SOURCE/clinvar_GRCh38.vcf.gz" \
  "$SOURCE/clinvar_GRCh38.vcf.gz.tbi" \
  "$SOURCE/variant_summary.txt.gz" \
  "$SOURCE/submission_summary.txt.gz" \
  "$SOURCE/archive_2.0/2023/clinvar_20231226.vcf.gz" \
  "$SOURCE/archive_2.0/2025/clinvar_20251221.vcf.gz"; do
  test -f "$path" || { echo "Missing: $path" >&2; exit 1; }
  ln -sfn "$path" "$TARGET/$(basename "$path")"
done

echo "Linked existing ClinVar assets into $TARGET"

