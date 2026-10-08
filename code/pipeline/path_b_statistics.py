"""Paired DeLong statistics using placement values, retaining tied scores."""
import numpy as np
from scipy.stats import norm, rankdata


def paired_delong(y, first, second):
    y = np.asarray(y, dtype=int)
    predictions = np.vstack([first, second]).astype(float)
    if not np.isfinite(predictions).all() or not set(y) <= {0, 1}:
        raise ValueError('Nonfinite scores or nonbinary labels')
    m, n = int(y.sum()), int((y == 0).sum())
    if min(m, n) < 2:
        return {'status': 'fewer_than_two_per_class'}
    order = np.argsort(-y, kind='stable')
    sorted_scores = predictions[:, order]
    tx = np.array([rankdata(s[:m]) for s in sorted_scores])
    ty = np.array([rankdata(s[m:]) for s in sorted_scores])
    tz = np.array([rankdata(s) for s in sorted_scores])
    aucs = tz[:, :m].sum(axis=1) / (m * n) - (m + 1) / (2 * n)
    covariance = np.cov((tz[:, :m] - tx) / n) / m + np.cov(1 - (tz[:, m:] - ty) / m) / n
    difference = float(aucs[0] - aucs[1])
    variance = float(covariance[0, 0] + covariance[1, 1] - 2 * covariance[0, 1])
    if variance <= 1e-15:
        return {'auroc_difference': difference, 'status': 'zero_variance',
                'p_value': 1.0 if abs(difference) < 1e-12 else None}
    se = float(np.sqrt(variance))
    return {'auroc_difference': difference, 'se': se,
            'p_value': float(2 * norm.sf(abs(difference / se))),
            'difference_ci95': [difference - 1.96 * se, difference + 1.96 * se]}
