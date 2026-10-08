"""Audit trusted AVI artifacts and perform a retrospective, unfitted benchmark."""
import hashlib
import json
import pickle
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve, precision_recall_curve

ROOT = Path('/path/to/project_root/AG-VUS')


def metrics(frame):
    result = {'n': len(frame), 'positive': int(frame.label.sum())}
    if frame.label.nunique() < 2:
        return dict(result, status='single_class_not_evaluable')
    return dict(result, auroc=float(roc_auc_score(frame.label, frame.avi)),
                auprc=float(average_precision_score(frame.label, frame.avi)),
                prevalence=float(frame.label.mean()))


def main():
    batch = ROOT / 'annotations/alphagenome/objective2_reclassified_full'
    if not (ROOT / 'logs/OBJECTIVE2_RECLASSIFIED_ATLAS.COMPLETE').exists():
        raise RuntimeError('Full annotation completion gate missing')
    output = ROOT / 'reports/objective2_benchmark'
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((batch / 'manifest.json').read_text())
    if len(manifest) != 42383 or len({r['canonical_id'] for r in manifest}) != 42383:
        raise RuntimeError('Manifest count or uniqueness failure')
    rows = []
    for index, row in enumerate(manifest, 1):
        token = hashlib.sha256(row['canonical_id'].encode()).hexdigest()[:24]
        record = json.loads((batch / 'variant_results' / (token + '.json')).read_text())
        if record['canonical_id'] != row['canonical_id'] or record['status'] != 'annotated':
            raise RuntimeError('Identity/status audit failed')
        files = record['response_files']
        for scorer in ('AVI_SCORE', 'AVI_SCORE_MODEL_FEATURES', 'AVI_SCORE_FEATURE_IMPORTANCE'):
            path = Path(files[scorer])
            if path.parent.parent.parent != ROOT / 'annotations/alphagenome':
                raise RuntimeError('Untrusted artifact directory')
            if path.name != f'{token}.{scorer}.pickle' or not path.is_file():
                raise RuntimeError('Missing or unexpected artifact')
        # Only unpickle our own recorded API responses, never external uploads.
        with Path(files['AVI_SCORE']).open('rb') as handle:
            data = pickle.load(handle)
        values = data.X.toarray() if hasattr(data.X, 'toarray') else np.asarray(data.X)
        if list(data.shape) != record['scorer_shapes']['AVI_SCORE'] or values.size != 1:
            raise RuntimeError('AVI must have exactly one score; refusing arbitrary averaging')
        score = float(values.ravel()[0])
        if not np.isfinite(score):
            raise RuntimeError('Nonfinite AVI score')
        rows.append({'canonical_id': row['canonical_id'], 'trajectory': row['trajectory'],
                     'avi': score, 'avi_obs_id': str(data.obs_names[0]),
                     'avi_feature_id': str(data.var_names[0]),
                     'queried_at_utc': record['queried_at_utc'],
                     'alphagenome_version': record['alphagenome_version']})
        if index % 1000 == 0:
            print(f'audited={index}/{len(manifest)}', flush=True)
    scores = pd.DataFrame(rows)
    con = duckdb.connect()
    con.execute("SET memory_limit='4GB'")
    con.execute('SET threads=2')
    history = con.execute('''SELECT canonical_id, molecular_consequence, geneinfo,
        first_vus_date, first_reclassification_date, class_changes, reversions
        FROM read_parquet(?) WHERE trajectory IN ('VUS_TO_P_LP','VUS_TO_B_LB')''',
        [str(ROOT / 'data/processed/longitudinal/objective1_frozen_cohorts_ref_validated.parquet')]).df()
    con.close()
    frame = scores.merge(history, on='canonical_id', how='left', validate='one_to_one', indicator=True)
    if not frame['_merge'].eq('both').all():
        raise RuntimeError('History join incomplete')
    frame = frame.drop(columns='_merge')
    frame['label'] = frame.trajectory.eq('VUS_TO_P_LP').astype(int)
    frame['reclassification_year'] = pd.to_datetime(frame.first_reclassification_date, errors='raise').dt.year
    temporary = output / 'variant_scores.tmp.parquet'
    frame.to_parquet(temporary, index=False)
    temporary.replace(output / 'variant_scores.parquet')
    groups = {'overall': metrics(frame), 'by_consequence': {}, 'by_reclassification_year': {}}
    for key, subset in frame.groupby('molecular_consequence', dropna=False):
        groups['by_consequence'][str(key)] = metrics(subset)
    for key, subset in frame.groupby('reclassification_year', dropna=False):
        groups['by_reclassification_year'][str(key)] = metrics(subset)
    stable = frame[(frame.reversions == 0) & (frame.class_changes == 1)]
    groups['no_reversion_single_change_sensitivity'] = metrics(stable)
    fpr, tpr, thresholds = roc_curve(frame.label, frame.avi)
    pd.DataFrame({'fpr': fpr, 'tpr': tpr, 'threshold': thresholds}).to_csv(output / 'roc_points.csv', index=False)
    precision, recall, thresholds = precision_recall_curve(frame.label, frame.avi)
    pd.DataFrame({'precision': precision, 'recall': recall,
                  'threshold': np.r_[thresholds, np.nan]}).to_csv(output / 'pr_points.csv', index=False)
    groups['boundary'] = ('Retrospective fixed-predictor evaluation, not prospective temporal validation. '
        'No fitting or threshold selection on these outcomes. Predictor training-overlap audit pending. '
        'AVI is not assumed calibrated: Brier/ECE omitted. AVI includes AlphaMissense-related evidence; '
        'do not count them independently. Consequence groups retain source annotations exactly.')
    (output / 'metrics.json').write_text(json.dumps(groups, indent=2, allow_nan=False) + '\n')
    (ROOT / 'logs/OBJECTIVE2_AVI_BENCHMARK.COMPLETE').write_text('Retrospective AVI benchmark exported\n')
    print(json.dumps(groups['overall']), flush=True)


if __name__ == '__main__':
    main()
