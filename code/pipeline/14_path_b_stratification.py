"""Reproducible post hoc Path B chromosome and consequence analyses."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

PARTITIONS = {
    'AVI_training_chromosomes': {'1', '4', '7', '8', '10', '13', '15'},
    'AVI_validation_chromosomes': {'2', '5', '11', '14', '17', '20', '22', 'X'},
    'outside_AVI_training_validation_chromosomes': {'3', '6', '9', '12', '16', '18', '19', '21', 'Y'},
}
PRIORITY = [
    ('SO:0001587', 'Nonsense'), ('SO:0001575', 'Splice donor'),
    ('SO:0001574', 'Splice acceptor'), ('SO:0001578', 'Stop lost'),
    ('SO:0001582', 'Start lost'), ('SO:0001583', 'Missense'),
    ('SO:0001819', 'Synonymous'), ('SO:0001627', 'Intronic'),
    ('SO:0001623', "5' UTR"), ('SO:0001624', "3' UTR"),
    ('SO:0001619', 'Non-coding transcript'),
]


def consequence(value):
    for token, name in PRIORITY:
        if token in str(value):
            return name
    return 'Other' if pd.notna(value) and value else 'Unannotated'


def evaluate(frame, nboot):
    y = frame.label.to_numpy(dtype=int)
    scores = frame.avi.to_numpy(dtype=float)
    result = {'n': len(y), 'events': int(y.sum()), 'non_events': int(len(y) - y.sum())}
    if len(y) == 0 or len(np.unique(y)) != 2:
        return dict(result, status='not_evaluable')
    result.update(prevalence=float(y.mean()), auroc=float(roc_auc_score(y, scores)),
                  auprc=float(average_precision_score(y, scores)))
    if nboot:
        positive, negative = np.flatnonzero(y), np.flatnonzero(y == 0)
        rng = np.random.default_rng(20260918)
        boot = np.empty((nboot, 2))
        for i in range(nboot):
            idx = np.r_[rng.choice(positive, len(positive), replace=True),
                        rng.choice(negative, len(negative), replace=True)]
            boot[i] = roc_auc_score(y[idx], scores[idx]), average_precision_score(y[idx], scores[idx])
        result['auroc_ci95'] = np.percentile(boot[:, 0], [2.5, 97.5]).tolist()
        result['auprc_ci95'] = np.percentile(boot[:, 1], [2.5, 97.5]).tolist()
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--bootstrap', type=int, default=2000)
    args = parser.parse_args()
    source = args.root / 'results/objective2_benchmark/variant_scores_with_external_predictors.parquet'
    if not source.exists():
        source = args.root / 'reports/objective2_benchmark/variant_scores.parquet'
    output = args.root / 'manuscript/path_b/stratification'
    output.mkdir(parents=True, exist_ok=True)
    frame = pd.read_parquet(source)
    if len(frame) != 42383 or frame.canonical_id.nunique() != len(frame):
        raise RuntimeError('Unexpected cohort')
    coords = frame.canonical_id.str.extract(r'^GRCh38_([0-9XY]+)_(\d+)_([ACGT])_([ACGT])$')
    if coords.isna().any().any():
        raise RuntimeError('Unexpected canonical variant representation')
    frame['chrom'] = coords[0]
    frame['consequence_class'] = frame.molecular_consequence.map(consequence)
    frame['label'] = frame.trajectory.eq('VUS_TO_P_LP').astype(int)
    rows = []
    for partition, chromosomes in PARTITIONS.items():
        subset = frame[frame.chrom.isin(chromosomes)]
        for scope, selected in [('all', subset), ('missense', subset[subset.consequence_class == 'Missense'])]:
            item = dict(partition=partition, scope=scope, **evaluate(selected, args.bootstrap))
            rows.append(item)
            print(json.dumps(item), flush=True)
    pd.DataFrame(rows).to_csv(output / 'chromosome_partition_metrics.csv', index=False)
    consequences = [dict(consequence=key, **evaluate(subset, args.bootstrap))
                    for key, subset in frame.groupby('consequence_class')]
    table = pd.DataFrame(consequences)
    table[(table.events >= 20) & (table.non_events >= 20)].to_csv(output / 'consequence_main.csv', index=False)
    table[(table.events < 20) | (table.non_events < 20)].to_csv(output / 'consequence_supplementary.csv', index=False)
    audit = {'source': str(source), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
             'bootstrap': args.bootstrap, 'seed': 20260918, 'variants': len(frame),
             'missense_variants': int((frame.consequence_class == 'Missense').sum()),
             'partition_counts_reconcile': sum(row['n'] for row in rows if row['scope'] == 'all') == len(frame),
             'chrY_variants': int((frame.chrom == 'Y').sum()),
             'boundary': 'Post hoc robustness analysis. Outside train/validation chromosomes is not the exact published coordinate-filtered test set.'}
    (output / 'audit.json').write_text(json.dumps(audit, indent=2) + '\n')


if __name__ == '__main__':
    main()
