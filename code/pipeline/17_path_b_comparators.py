"""Extract full primary-source AlphaMissense/REVEL coverage and compare identical cohorts."""
import argparse
import gzip
import hashlib
import json
import zipfile
from pathlib import Path

import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, precision_recall_curve

from path_b_statistics import paired_delong
from importlib import import_module
metrics = import_module('16_path_b_af_analysis').metrics
collapse = import_module('14_path_b_stratification').consequence
evaluate = import_module('14_path_b_stratification').evaluate


def extract(frame, source, predictor):
    matches = []
    if predictor == 'am':
        digest = hashlib.md5()
        with source.open('rb') as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(chunk)
        if digest.hexdigest() != '9fd167735f16a1b87da6eb3e4c25fcb5':
            raise RuntimeError('AlphaMissense official MD5 mismatch')
        names = ['chrom', 'pos', 'ref', 'alt', 'genome', 'uniprot_id', 'transcript_id', 'protein_variant', 'am_pathogenicity', 'am_class']
        handles = [gzip.open(source, 'rt')]
    else:
        archive = zipfile.ZipFile(source)
        handles = [archive.open(name) for name in archive.namelist()
                   if not name.endswith('/') and not name.startswith('__MACOSX/')]
        if not handles:
            raise RuntimeError('No REVEL CSV members')
    keys = ['chrom', 'pos', 'ref', 'alt']
    for handle in handles:
        options = dict(sep='\t', names=names, comment='#') if predictor == 'am' else {}
        for chunk in pd.read_csv(handle, chunksize=100000, dtype=str, **options):
            if predictor == 'revel':
                chunk = chunk.rename(columns={'chr': 'chrom', 'grch38_pos': 'pos', 'REVEL': 'revel_score'})
            chunk['chrom'] = chunk.chrom.str.removeprefix('chr')
            chunk['pos'] = pd.to_numeric(chunk.pos, errors='coerce')
            chunk = chunk.dropna(subset=['pos'])
            hit = frame[keys].merge(chunk, on=keys, how='inner')
            if len(hit):
                matches.append(hit)
        handle.close()
    if predictor == 'revel':
        archive.close()
    if not matches:
        raise RuntimeError('No exact matches from primary source')
    return pd.concat(matches, ignore_index=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    output = root / 'manuscript/path_b/comparators'
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.read_parquet(root / 'reports/objective2_benchmark/variant_scores.parquet')
    parts = frame.canonical_id.str.extract(r'^GRCh38_([0-9XY]+)_(\d+)_([ACGT])_([ACGT])$')
    for name, value in zip(['chrom', 'pos', 'ref', 'alt'], [parts[0], parts[1].astype(int), parts[2], parts[3]]):
        frame[name] = value
    keys = ['chrom', 'pos', 'ref', 'alt']
    audit = {}
    for predictor, filename, column in [('am', 'AlphaMissense_hg38.tsv.gz', 'am_pathogenicity'),
                                         ('revel', 'revel-v1.3_all_chromosomes.zip', 'revel_score')]:
        matches = extract(frame, root / 'data/path_b_sources' / filename, predictor)
        matches[column] = pd.to_numeric(matches[column], errors='raise')
        matches.to_parquet(output / f'{predictor}_all_exact_transcript_matches.parquet', index=False)
        counts = matches.groupby(keys)[column].nunique()
        ambiguous = counts[counts > 1].reset_index()[keys]
        filtered = matches.merge(ambiguous.assign(ambiguous=True), on=keys, how='left')
        filtered = filtered[filtered.ambiguous.isna()].drop_duplicates(keys)
        frame = frame.merge(filtered[keys + [column]], on=keys, how='left', validate='one_to_one')
        audit[predictor] = {'exact_alleles': len(counts), 'ambiguous_score_alleles_abstained': len(ambiguous),
                            'scored_variants': int(frame[column].notna().sum())}
    for dataset, column in [('cadd', 'cadd_phred'), ('spliceai', 'spliceai_max_ds')]:
        path = root / 'annotations/path_b' / dataset
        if not (path / 'COMPLETE.json').exists():
            raise RuntimeError(f'Comparator extraction incomplete: {dataset}')
        scores = pd.read_parquet(path / 'all_variants.parquet')
        frame = frame.merge(scores[['canonical_id', column]], on='canonical_id', how='left', validate='one_to_one')
        audit[dataset] = {'scored_variants': int(frame[column].notna().sum())}
    frame['consequence_class'] = frame.molecular_consequence.map(collapse)
    frame.to_parquet(output / 'full_primary_source_comparator_scores.parquet', index=False)
    report = {'coverage': audit, 'comparisons': {}, 'boundary': 'Post hoc exact GRCh38 allele joins; ambiguous transcript scores abstain. Ensembl SpliceAI raw MANE v1.4 D50 has transcript-limited coverage. Missing scores remain missing.'}
    for name, columns, selected in [('four_way', ['avi', 'am_pathogenicity', 'cadd_phred', 'revel_score'], frame),
                                     ('splice_donor_acceptor', ['avi', 'cadd_phred', 'spliceai_max_ds'], frame[frame.consequence_class.isin(['Splice donor', 'Splice acceptor'])])]:
        selected = selected.dropna(subset=columns)
        report['comparisons'][name] = {column: evaluate(selected.assign(avi=selected[column]), 2000) for column in columns}
        report['comparisons'][name]['paired_delong_vs_avi'] = {column: paired_delong(selected.label, selected.avi, selected[column]) for column in columns[1:]}
        pairs = report['comparisons'][name]['paired_delong_vs_avi']
        valid = sorted([(key, value['p_value']) for key, value in pairs.items() if value.get('p_value') is not None], key=lambda item: item[1])
        running = 0
        for rank, (key, p_value) in enumerate(valid):
            running = max(running, min(1, (len(valid) - rank) * p_value))
            pairs[key]['holm_p_value'] = running
        if selected.label.nunique() == 2:
            for column in columns:
                fpr, tpr, thresholds = roc_curve(selected.label, selected[column])
                pd.DataFrame({'fpr': fpr, 'tpr': tpr, 'threshold': thresholds}).to_csv(output / f'{name}_{column}_roc.csv', index=False)
                precision, recall, thresholds = precision_recall_curve(selected.label, selected[column])
                pd.DataFrame({'precision': precision, 'recall': recall, 'threshold': np.r_[thresholds, np.nan]}).to_csv(output / f'{name}_{column}_pr.csv', index=False)
        selected.to_parquet(output / f'{name}_matched_variants.parquet', index=False)
    (output / 'comparator_results.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps(audit), flush=True)


if __name__ == '__main__':
    main()
