# Artificial intelligence-based AlphaGenome scores and ClinVar VUS reclassification: data and code

This repository benchmarks the AlphaGenome Variant Impact (AVI) score, released with the [AlphaGenome Atlas](https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/), against the later reclassification of ClinVar variants of uncertain significance (VUS).

Data and code for: Banaganapalli B. *Artificial Intelligence-Based AlphaGenome Scores and the Reclassification of ClinVar Variants of Uncertain Significance.* Manuscript submitted for publication.

Contact: Babajan Banaganapalli, Department of Genetic Medicine, Faculty of Medicine; Princess Al-Jawhara Center of Excellence in Research of Hereditary Disorders; Center of Excellence in High Performance Computing (Aziz); King Abdulaziz University, Jeddah, Saudi Arabia (bbabajan@kau.edu.sa).

## Contents

### `data/analysed_cohort_42383.csv.gz`

One row per analysed single-nucleotide variant (GRCh38). Each was a ClinVar VUS later reclassified to benign/likely benign or pathogenic/likely pathogenic across 420 dated ClinVar releases (2017–2026).

| Column | Meaning |
|---|---|
| `canonical_id` | `GRCh38_<chrom>_<pos>_<ref>_<alt>` |
| `chrom`, `pos`, `ref`, `alt` | GRCh38 coordinate (1-based) and alleles |
| `gene_symbol` | ClinVar GENEINFO symbol(s), `;`-separated |
| `trajectory` | `VUS_TO_B_LB` or `VUS_TO_P_LP` |
| `later_pathogenic` | Outcome: 1 = later pathogenic/likely pathogenic, 0 = later benign/likely benign |
| `first_vus_date`, `first_reclassification_date` | ClinVar release dates (YYYYMMDD) |
| `class_changes`, `reversions` | Number of classification changes and reversions in the history |
| `consequence_class`, `molecular_consequence` | Collapsed class used in the paper, and ClinVar MC field |
| `avi_raw` | AlphaGenome Variant Impact score, raw (before PHRED scaling). AlphaGenome output: see Licence |
| `avi_partition` | AVI data split by chromosome: `training`, `validation`, or `held_out` (primary analysis) |
| `alphagenome_sdk_version`, `queried_at_utc` | Query provenance |

Comparator scores (AlphaMissense, CADD, REVEL, SpliceAI) are not included because of their licences. They can be re-attached by exact match on `chrom`, `pos`, `ref`, `alt` from the sources cited in the paper.

### `results/`

- `results_v2.json`: every statistic reported in the paper.
- `circularity_by_date.json`: the four-predictor comparison split at the AlphaMissense release date (Supplementary Table S6).
- `S*.csv`: the stratified result tables behind the Supplementary Tables.

### `code/`

- `pipeline/`: release download and parsing, longitudinal reconstruction, reference validation, AlphaGenome querying, comparator annotation, and the PBS job files. Paths are written as `/path/to/project_root` and `/path/to/AG-VUS`. The AlphaGenome API key is read from the `ALPHA_GENOME_API_KEY` environment variable and is not included.
- `analysis/`: statistics (`stats_v2.py`, `noncoding.py`) and the label-circularity check (`circularity_by_date.py`).

## Reproducing the analysis

```bash
pip install pandas numpy scipy scikit-learn matplotlib pyarrow
```

- `results/` holds every reported statistic, so the paper's numbers can be checked without rerunning anything.
- `code/analysis/` (statistics and circularity check) expects the analysed cohort with comparator scores re-attached. AlphaMissense, CADD, REVEL and SpliceAI scores are not redistributed here; join them to `data/analysed_cohort_42383.csv.gz` on `chrom`, `pos`, `ref`, `alt`, then edit the input paths at the top of each script.
- `code/pipeline/` was run on a PBS cluster. Paths are written as `/path/to/project_root`. AlphaGenome queries need an API key in the `ALPHA_GENOME_API_KEY` environment variable.

## Licence

- Code: MIT (see `LICENSE`).
- Data and result tables: CC BY 4.0, **except the `avi_raw` column**.
- **AlphaGenome output notice.** The `avi_raw` column of `data/analysed_cohort_42383.csv.gz`, and any value derived from it, is AlphaGenome output obtained through the AlphaGenome API. It is provided under and subject to the AlphaGenome Output Terms of Use (https://deepmind.google.com/science/alphagenome/output-terms). By using this information, you agree to the AlphaGenome Output Terms of Use. It is not covered by CC BY 4.0. See `data/NOTICE.txt`.
- When using the AVI scores, please cite AlphaGenome (Avsec Ž *et al.* Advancing regulatory variant effect prediction with AlphaGenome. *Nature* 2026; 649: 1206–1218. https://doi.org/10.1038/s41586-025-10014-0) and the AlphaGenome Atlas (Cheng J *et al.* AlphaGenome Atlas: in silico mutagenesis of the entire human genome improves prioritization and interpretation of non-coding variants. *medRxiv* 2026. https://doi.org/10.64898/2026.09.16.26363192).

## Citation

If you use these data or code, please cite the article (citation to be added on publication).
