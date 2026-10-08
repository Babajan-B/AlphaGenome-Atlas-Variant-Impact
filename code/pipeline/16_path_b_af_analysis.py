"""Paired AF baselines, AF strata and deterministic frequency-matched analyses."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve, precision_recall_curve

from path_b_statistics import paired_delong
from importlib import import_module
collapse = import_module('14_path_b_stratification').consequence
PARTITIONS = import_module('14_path_b_stratification').PARTITIONS


def metrics(frame, column):
    result = {'n': len(frame), 'events': int(frame.label.sum())}
    if frame.label.nunique() < 2:
        return dict(result, status='single_class_or_empty')
    return dict(result, prevalence=float(frame.label.mean()),
                auroc=float(roc_auc_score(frame.label, frame[column])),
                auprc=float(average_precision_score(frame.label, frame[column])))


def comparison(frame):
    if frame.label.nunique() < 2:
        return {'n': len(frame), 'status': 'single_class_or_empty'}
    return {'avi': metrics(frame, 'avi'), 'af_baseline': metrics(frame, 'af_baseline'),
            'paired_delong': paired_delong(frame.label, frame.avi, frame.af_baseline)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    scores = pd.read_parquet(args.root / 'reports/objective2_benchmark/variant_scores.parquet')
    output = args.root / 'manuscript/path_b/frequency_analysis'
    output.mkdir(parents=True, exist_ok=True)
    scores['consequence_class'] = scores.molecular_consequence.map(collapse)
    report = {'epsilon': 1e-12, 'matching_seed': 20260918, 'datasets': {},
              'boundary': 'Post hoc retrospective robustness analyses. AUROC differences are not fractions of explained biological signal. Unobserved variants are not proof of adequate gnomAD coverage.'}
    for dataset, fields in [('genomes', ['fafmax_faf95_max', 'AF']),
                            ('joint', ['fafmax_faf95_max_joint', 'AF_joint'])]:
        source = args.root / 'annotations/path_b' / dataset
        if not (source / 'COMPLETE.json').exists():
            raise RuntimeError(f'Frequency extraction incomplete: {dataset}')
        frequencies = pd.read_parquet(source / 'all_variants.parquet')
        frame = scores.merge(frequencies, on='canonical_id', how='left', validate='one_to_one')
        if len(frame) != 42383 or frame.extraction_status.isna().any():
            raise RuntimeError('Frequency joins do not reconcile')
        report['datasets'][dataset] = {}
        for field in fields:
            selected = frame.copy()
            absent = selected.extraction_status.eq('absent_after_successful_query')
            observed_missing = (~absent) & selected[field].isna()
            selected['frequency'] = selected[field].mask(absent, 0)
            selected = selected[~observed_missing].copy()
            if not selected.frequency.between(0, 1).all():
                raise RuntimeError('Invalid allele frequency')
            selected['af_baseline'] = -np.log10(selected.frequency + 1e-12)
            selected['af_band'] = pd.cut(selected.frequency, [0, 1e-5, 1e-4, 1e-3, np.inf],
                                        right=False, labels=['observed_<1e-5', '1e-5_to_<1e-4', '1e-4_to_<1e-3', '>=1e-3']).astype(str)
            selected.loc[selected.extraction_status.eq('absent_after_successful_query'), 'af_band'] = 'absent'
            item = {'observed_missing_excluded': int(observed_missing.sum()),
                    'absent_imputed_zero': int(absent.sum()), 'overall': comparison(selected),
                    'missense': comparison(selected[selected.consequence_class == 'Missense']),
                    'outside_training_validation_chromosomes': comparison(selected[selected.chrom.isin(PARTITIONS['outside_AVI_training_validation_chromosomes'])]),
                    'af_bands': {key: comparison(subset) for key, subset in selected.groupby('af_band')}}
            selected['match_bin'] = 'observed_zero'
            positive_frequency = selected.frequency > 0
            if positive_frequency.any():
                bins = pd.qcut(selected.loc[positive_frequency, 'frequency'], q=10, duplicates='drop')
                selected.loc[positive_frequency, 'match_bin'] = bins.astype(str)
            selected.loc[selected.extraction_status.eq('absent_after_successful_query'), 'match_bin'] = 'absent'
            matches = []
            for _, subset in selected.groupby('match_bin'):
                size = min(int((subset.label == 0).sum()), int((subset.label == 1).sum()))
                if size:
                    matches.extend([subset[subset.label == label].sample(size, random_state=20260918) for label in [0, 1]])
            if matches:
                matched = pd.concat(matches)
                item['af_matched'] = comparison(matched)
                matched.to_parquet(output / f'{dataset}_{field}_af_matched.parquet', index=False)
            selected.to_parquet(output / f'{dataset}_{field}_variants.parquet', index=False)
            for name in ['avi', 'af_baseline']:
                fpr, tpr, thresholds = roc_curve(selected.label, selected[name])
                pd.DataFrame({'fpr': fpr, 'tpr': tpr, 'threshold': thresholds}).to_csv(output / f'{dataset}_{field}_{name}_roc.csv', index=False)
            report['datasets'][dataset][field] = item
            print(json.dumps({'dataset': dataset, 'field': field, 'overall': item['overall']}), flush=True)
    (output / 'frequency_results.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
