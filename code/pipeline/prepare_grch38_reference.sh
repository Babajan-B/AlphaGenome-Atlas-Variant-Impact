#!/usr/bin/env bash
set -euo pipefail
ROOT=/path/to/project_root/AG-VUS
DIR="$ROOT/data/reference"
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mkdir -p "$DIR"
cd "$DIR"
if [[ ! -s GRCh38_no_alt_analysis_set.fna.gz ]]; then
  curl -fL --retry 5 --retry-delay 10 "$URL" -o GRCh38_no_alt_analysis_set.fna.gz.part
  mv GRCh38_no_alt_analysis_set.fna.gz.part GRCh38_no_alt_analysis_set.fna.gz
fi
if [[ ! -s GRCh38_no_alt_analysis_set.fna ]]; then
  gzip -dc GRCh38_no_alt_analysis_set.fna.gz > GRCh38_no_alt_analysis_set.fna.part
  mv GRCh38_no_alt_analysis_set.fna.part GRCh38_no_alt_analysis_set.fna
fi
"$ROOT/environment/alphagenome-py310/bin/python" - <<'PY'
import pysam
pysam.faidx('/path/to/project_root/AG-VUS/data/reference/GRCh38_no_alt_analysis_set.fna')
PY
touch "$ROOT/logs/GRCH38_REFERENCE.READY"
