#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Non-coding / synonymous comparison: AVI against the genome-wide comparators
in the variant classes where AlphaMissense and REVEL return nothing.
Appends to results_v2.json.
"""
import json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve

warnings.filterwarnings("ignore")
OUT = "/path/to/output"
SEED = 20260918
NONCODE = ["Synonymous", "Intronic", "3' UTR", "5' UTR", "Non-coding transcript"]


def _midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x); T = np.zeros(N); i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N); T2[J] = T
    return T2


def delong(y, s1, s2):
    y = np.asarray(y, float); o = (-y).argsort(kind="mergesort"); m = int(y.sum())
    P = np.vstack([np.asarray(s1, float)[o], np.asarray(s2, float)[o]]); n = P.shape[1] - m
    tx = np.array([_midrank(P[r, :m]) for r in range(2)])
    ty = np.array([_midrank(P[r, m:]) for r in range(2)])
    tz = np.array([_midrank(P[r, :]) for r in range(2)])
    auc = tz[:, :m].sum(1) / m / n - (m + 1) / 2 / n
    v01 = (tz[:, :m] - tx) / n; v10 = 1 - (tz[:, m:] - ty) / m
    cov = np.cov(v01) / m + np.cov(v10) / n
    diff = float(auc[0] - auc[1]); var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    se = float(np.sqrt(var))
    return diff, se, float(2 * stats.norm.sf(abs(diff / se))), (diff - 1.96 * se, diff + 1.96 * se)


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); a = np.empty_like(p); run = 0.0
    for r, i in enumerate(o):
        run = max(run, (len(p) - r) * p[i]); a[i] = min(1.0, run)
    return a


def boot(y, s, nboot=2000, seed=SEED):
    y = np.asarray(y); s = np.asarray(s)
    pos = np.flatnonzero(y == 1); neg = np.flatnonzero(y == 0)
    rng = np.random.default_rng(seed)
    r = np.empty(nboot); p = np.empty(nboot)
    for i in range(nboot):
        idx = np.concatenate([rng.choice(pos, pos.size, True), rng.choice(neg, neg.size, True)])
        r[i] = roc_auc_score(y[idx], s[idx]); p[i] = average_precision_score(y[idx], s[idx])
    return np.percentile(r, [2.5, 97.5]), np.percentile(p, [2.5, 97.5])


def met(y, s, nboot=2000):
    y = np.asarray(y); s = np.asarray(s); prev = float(y.mean())
    d = dict(n=int(len(y)), events=int(y.sum()), non_events=int((y == 0).sum()),
             prevalence=round(prev, 4))
    if y.sum() < 2 or (y == 0).sum() < 2:
        d["status"] = "too_few"; return d
    rc, pc = boot(y, s, nboot)
    d.update(auroc=round(float(roc_auc_score(y, s)), 4), auroc_lo=round(rc[0], 4),
             auroc_hi=round(rc[1], 4), auprc=round(float(average_precision_score(y, s)), 4),
             auprc_lo=round(pc[0], 4), auprc_hi=round(pc[1], 4),
             auprc_baseline=round(prev, 4),
             auprc_lift=round(float(average_precision_score(y, s)) - prev, 4))
    return d


df = pd.read_parquet(f"{OUT}/analysis_frame.parquet")
nc = df[df.cc.isin(NONCODE)].copy()
R = json.load(open(f"{OUT}/results_v2.json"))

block = dict(
    classes_included=NONCODE,
    n=int(len(nc)), events=int(nc.label.sum()), prevalence=round(float(nc.label.mean()), 4),
    alphamissense_coverage=round(float(nc.am_pathogenicity.notna().mean()), 4),
    revel_coverage=round(float(nc.revel_score.notna().mean()), 4),
    cadd_coverage=round(float(nc.cadd_phred.notna().mean()), 4),
    spliceai_coverage=round(float(nc.spliceai_max_ds.notna().mean()), 4),
)

y = nc.label.to_numpy()
block["pooled"] = {
    "avi": met(y, nc.avi), "cadd_phred": met(y, nc.cadd_phred)}
sm = nc.spliceai_max_ds.notna().to_numpy()
block["pooled"]["spliceai_max_ds"] = met(y[sm], nc.spliceai_max_ds[sm])
block["spliceai_subset_n"] = int(sm.sum())

pairs = []
for b in ["cadd_phred", "spliceai_max_ds"]:
    if b == "spliceai_max_ds":
        d, se, p, ci = delong(y[sm], nc.avi[sm], nc.spliceai_max_ds[sm])
    else:
        d, se, p, ci = delong(y, nc.avi, nc.cadd_phred)
    pairs.append(dict(a="AlphaGenome AVI", b=b, delta=round(d, 4), se=round(se, 5),
                      ci=[round(ci[0], 4), round(ci[1], 4)], p=p))
for k, a in zip(pairs, holm([x["p"] for x in pairs])):
    k["p_holm"] = float(a)
block["pooled_pairwise_vs_avi"] = pairs

# held-out chromosomes only
h = nc[nc.partition.eq("held_out")]
yh = h.label.to_numpy()
block["held_out"] = {"avi": met(yh, h.avi), "cadd_phred": met(yh, h.cadd_phred)}
dh, seh, ph, cih = delong(yh, h.avi, h.cadd_phred)
block["held_out_avi_vs_cadd"] = dict(delta=round(dh, 4), se=round(seh, 5),
                                     ci=[round(cih[0], 4), round(cih[1], 4)], p=ph)

# per class
byc = []
for cc in NONCODE:
    g = nc[nc.cc.eq(cc)]
    y2 = g.label.to_numpy()
    row = dict(consequence=cc, n=int(len(g)), events=int(y2.sum()),
               prevalence=round(float(y2.mean()), 4) if len(g) else None)
    if y2.sum() >= 2 and (y2 == 0).sum() >= 2:
        row["avi"] = met(y2, g.avi, 1000)
        row["cadd_phred"] = met(y2, g.cadd_phred, 1000)
        d, se, p, ci = delong(y2, g.avi, g.cadd_phred)
        row["avi_vs_cadd"] = dict(delta=round(d, 4), ci=[round(ci[0], 4), round(ci[1], 4)], p=p)
    else:
        row["status"] = "too_few_events"
    byc.append(row)
block["by_class"] = byc

# ROC and precision-recall curve coordinates
cur = {}
for k, arr, mask in [("avi", nc.avi, np.ones(len(nc), bool)),
                     ("cadd_phred", nc.cadd_phred, np.ones(len(nc), bool)),
                     ("spliceai_max_ds", nc.spliceai_max_ds, sm)]:
    f2, t2, _ = roc_curve(y[mask], arr[mask])
    p2, r2, _ = precision_recall_curve(y[mask], arr[mask])
    cur[k] = dict(fpr=f2.tolist(), tpr=t2.tolist(), rec=r2.tolist(), pre=p2.tolist())
json.dump(cur, open(f"{OUT}/noncoding_curves.json", "w"))
nc.to_parquet(f"{OUT}/noncoding_frame.parquet")

R["noncoding_comparison"] = block
json.dump(R, open(f"{OUT}/results_v2.json", "w"), indent=2, default=float)
pd.DataFrame([dict(scope="pooled", predictor=k, **v) for k, v in block["pooled"].items()]
             + [dict(scope="held_out", predictor=k, **v) for k, v in block["held_out"].items()]
             ).to_csv(f"{OUT}/T_noncoding.csv", index=False)

print(json.dumps({k: v for k, v in block.items() if k != "by_class"}, indent=2, default=float))
print("\nBY CLASS")
for r in byc:
    if "avi" in r:
        print(f"  {r['consequence']:22s} n={r['n']:6,d} ev={r['events']:4d} prev={r['prevalence']:.4f} "
              f"AVI={r['avi']['auroc']:.4f} ({r['avi']['auroc_lo']:.4f}-{r['avi']['auroc_hi']:.4f})  "
              f"CADD={r['cadd_phred']['auroc']:.4f}  d={r['avi_vs_cadd']['delta']:+.4f} p={r['avi_vs_cadd']['p']:.3g}")
    else:
        print(f"  {r['consequence']:22s} n={r['n']:6,d} ev={r['events']:4d}  {r['status']}")
