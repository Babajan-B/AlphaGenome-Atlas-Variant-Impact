#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Project 9 v2 statistics: independent recomputation of all Path B results
including the frequency, comparator and sensitivity analyses.

Every figure quoted in the manuscript is produced here from the Path B
parquet outputs, so the manuscript numbers and the archived tables cannot
drift apart.
"""
import json, os, warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve, precision_recall_curve

warnings.filterwarnings("ignore")
P = "/path/to/project_root/manuscript/path_b/"   # Path B outputs from code/pipeline
OUT = "/path/to/output"
os.makedirs(OUT, exist_ok=True)
NBOOT = 2000
SEED = 20260918

TRAIN_CHR = {"1", "4", "7", "8", "10", "13", "15"}
VAL_CHR = {"2", "5", "11", "14", "17", "20", "22", "X"}

SO_PRIORITY = [
    ("SO:0001587|nonsense", "Nonsense"),
    ("SO:0001575|splice_donor_variant", "Splice donor"),
    ("SO:0001574|splice_acceptor_variant", "Splice acceptor"),
    ("SO:0001578|stop_lost", "Stop lost"),
    ("SO:0001582|initiator_codon_variant", "Start lost"),
    ("SO:0001583|missense_variant", "Missense"),
    ("SO:0001819|synonymous_variant", "Synonymous"),
    ("SO:0001627|intron_variant", "Intronic"),
    ("SO:0001623|5_prime_UTR_variant", "5' UTR"),
    ("SO:0001624|3_prime_UTR_variant", "3' UTR"),
    ("SO:0001619|non-coding_transcript_variant", "Non-coding transcript"),
]


def collapse(mc):
    if not isinstance(mc, str) or mc == "":
        return "Unannotated"
    for token, label in SO_PRIORITY:
        if token in mc:
            return label
    return "Other"


def norm_chr(c):
    return str(c).replace("chr", "")


def partition(c):
    c = norm_chr(c)
    return "training" if c in TRAIN_CHR else ("validation" if c in VAL_CHR else "held_out")


# ---------------------------------------------------------------- DeLong
def _midrank(x):
    J = np.argsort(x); Z = x[J]; N = len(x); T = np.zeros(N, float); i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N, float); T2[J] = T
    return T2


def _fast_delong(preds_sorted, m):
    n = preds_sorted.shape[1] - m
    pos = preds_sorted[:, :m]; neg = preds_sorted[:, m:]
    k = preds_sorted.shape[0]
    tx = np.empty([k, m]); ty = np.empty([k, n]); tz = np.empty([k, m + n])
    for r in range(k):
        tx[r, :] = _midrank(pos[r, :]); ty[r, :] = _midrank(neg[r, :])
        tz[r, :] = _midrank(preds_sorted[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - float(m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx[:, :]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:, :]) / m
    return aucs, np.atleast_2d(np.cov(v01) / m + np.cov(v10) / n)


def delong_pair(y, s1, s2):
    y = np.asarray(y, float)
    order = (-y).argsort(kind="mergesort"); m = int(y.sum())
    preds = np.vstack([np.asarray(s1, float)[order], np.asarray(s2, float)[order]])
    aucs, cov = _fast_delong(preds, m)
    diff = float(aucs[0] - aucs[1])
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    if var <= 0:
        return diff, np.nan, np.nan, (np.nan, np.nan)
    se = float(np.sqrt(var)); z = diff / se
    return diff, se, float(2 * stats.norm.sf(abs(z))), (diff - 1.96 * se, diff + 1.96 * se)


def holm(pvals):
    pv = np.asarray(pvals, float); order = np.argsort(pv)
    adj = np.empty_like(pv); run = 0.0
    for rank, idx in enumerate(order):
        run = max(run, (len(pv) - rank) * pv[idx])
        adj[idx] = min(1.0, run)
    return adj


def boot_ci(y, s, nboot=NBOOT, seed=SEED):
    y = np.asarray(y); s = np.asarray(s)
    pos = np.flatnonzero(y == 1); neg = np.flatnonzero(y == 0)
    if pos.size == 0 or neg.size == 0:
        return (np.nan, np.nan), (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    r = np.empty(nboot); p = np.empty(nboot)
    for i in range(nboot):
        idx = np.concatenate([rng.choice(pos, pos.size, True), rng.choice(neg, neg.size, True)])
        r[i] = roc_auc_score(y[idx], s[idx]); p[i] = average_precision_score(y[idx], s[idx])
    return tuple(np.percentile(r, [2.5, 97.5])), tuple(np.percentile(p, [2.5, 97.5]))


def metrics(y, s, nboot=NBOOT, label=""):
    """AUROC/AUPRC with CIs, prevalence baseline and AUPRC lift over baseline."""
    y = np.asarray(y); s = np.asarray(s)
    prev = float(y.mean())
    d = dict(n=int(len(y)), events=int(y.sum()),
             non_events=int((y == 0).sum()), prevalence=round(prev, 4))
    if y.sum() == 0 or (y == 0).sum() == 0:
        d["status"] = "single_class_not_evaluable"
        return d
    roc = float(roc_auc_score(y, s)); prc = float(average_precision_score(y, s))
    rc, pc = boot_ci(y, s, nboot)
    d.update(auroc=round(roc, 4), auroc_lo=round(rc[0], 4), auroc_hi=round(rc[1], 4),
             auprc=round(prc, 4), auprc_lo=round(pc[0], 4), auprc_hi=round(pc[1], 4),
             auprc_baseline=round(prev, 4), auprc_lift=round(prc - prev, 4))
    return d


R = {"seed": SEED, "bootstrap_replicates": NBOOT}

# ============================================================ load
comp = pd.read_parquet(P + "comparators/full_primary_source_comparator_scores.parquet")
comp["cc"] = comp["molecular_consequence"].map(collapse)
comp["partition"] = comp["chrom"].map(partition)
comp["gene"] = comp["geneinfo"].fillna("").str.split("|").str[0].str.split(":").str[0].replace("", "Unassigned")
y_all = comp["label"].to_numpy(); s_all = comp["avi"].to_numpy()

R["cohort"] = dict(n=int(len(comp)), events=int(y_all.sum()),
                   non_events=int((y_all == 0).sum()), prevalence=round(float(y_all.mean()), 4))

# ============================================================ 1. chromosome partitions (PRIMARY)
part_rows = []
for pt in ["held_out", "validation", "training"]:
    for scope, mask in [("all", comp.partition.eq(pt)),
                        ("missense", comp.partition.eq(pt) & comp.cc.eq("Missense"))]:
        g = comp[mask]
        part_rows.append(dict(partition=pt, scope=scope,
                              **metrics(g.label, g.avi, label=f"{pt}/{scope}")))
part = pd.DataFrame(part_rows)
part.to_csv(f"{OUT}/T_chromosome_partitions.csv", index=False)
R["chromosome_partitions"] = part_rows
R["partition_chromosomes"] = {
    "training": sorted(TRAIN_CHR, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x)),
    "validation": sorted(VAL_CHR, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x)),
    "held_out": sorted({norm_chr(c) for c in comp.chrom} - TRAIN_CHR - VAL_CHR,
                       key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))}
R["chrY_variants"] = int((comp.chrom.map(norm_chr) == "Y").sum())

# held-out is the primary analysis population
held = comp[comp.partition.eq("held_out")]
yh = held.label.to_numpy(); sh = held.avi.to_numpy()
R["primary_held_out"] = metrics(yh, sh, label="held_out_primary")

# distributions on held-out
pos_s, neg_s = sh[yh == 1], sh[yh == 0]
u, p_mw = stats.mannwhitneyu(pos_s, neg_s, alternative="two-sided")
R["distributions_held_out"] = dict(
    pos_median=round(float(np.median(pos_s)), 3), pos_q1=round(float(np.percentile(pos_s, 25)), 3),
    pos_q3=round(float(np.percentile(pos_s, 75)), 3),
    neg_median=round(float(np.median(neg_s)), 3), neg_q1=round(float(np.percentile(neg_s, 25)), 3),
    neg_q3=round(float(np.percentile(neg_s, 75)), 3),
    mannwhitney_p="<0.001" if p_mw < 1e-3 else f"{p_mw:.3g}",
    cliffs_delta=round(2 * roc_auc_score(yh, sh) - 1, 4))

# ============================================================ 2. consequence strata + prevalence baseline
cons_rows = []
for cc, g in comp.groupby("cc"):
    cons_rows.append(dict(consequence=cc, **metrics(g.label, g.avi,
                                                    nboot=1000 if len(g) >= 50 else 200, label=cc)))
cons = pd.DataFrame(cons_rows).sort_values("n", ascending=False)
cons["main_table"] = (cons["events"] >= 20) & (cons["non_events"] >= 20)
cons.to_csv(f"{OUT}/T_consequence_with_baseline.csv", index=False)
R["by_consequence"] = cons.to_dict("records")

# consequence-class prior benchmark (leave-one-out) and stratified concordance
grp = comp.groupby("cc")["label"]
cnt = grp.transform("count").to_numpy(); tot = grp.transform("sum").to_numpy()
prior = (tot - y_all) / np.maximum(cnt - 1, 1)
diff, se, pv, ci = delong_pair(y_all, s_all, prior)
R["consequence_prior"] = dict(**metrics(y_all, prior, label="consequence_prior"),
                              avi_minus_prior=round(diff, 4), delong_p=pv,
                              delong_ci=[round(ci[0], 4), round(ci[1], 4)])
num = den = 0.0
for cc, g in comp.groupby("cc"):
    yy = g.label.to_numpy()
    if yy.sum() and (yy == 0).sum():
        w = yy.sum() * (yy == 0).sum()
        num += roc_auc_score(yy, g.avi) * w; den += w
R["consequence_stratified_auroc"] = round(num / den, 4)

# ============================================================ 3. frequency decomposition
freq = {}
for tag, fn, fld in [("genomes_faf95", "frequency_analysis/genomes_fafmax_faf95_max_variants.parquet",
                      "fafmax_faf95_max"),
                     ("joint_faf95", "frequency_analysis/joint_fafmax_faf95_max_joint_variants.parquet",
                      "fafmax_faf95_max_joint")]:
    d = pd.read_parquet(P + fn)
    d["cc"] = d["molecular_consequence"].map(collapse)
    d["partition"] = d["chrom"].map(partition)
    blk = dict(field=fld, n=int(len(d)),
               absent_imputed=int((d["af_band"] == "absent").sum()) if "af_band" in d else None)
    blk["avi_overall"] = metrics(d.label, d.avi, label="avi_overall")
    blk["af_overall"] = metrics(d.label, d.af_baseline, label="af_overall")
    df_, se_, p_, ci_ = delong_pair(d.label, d.avi, d.af_baseline)
    blk["avi_vs_af"] = dict(difference=round(df_, 4), se=round(se_, 5), p=p_,
                            ci=[round(ci_[0], 4), round(ci_[1], 4)])
    h = d[d.partition.eq("held_out")]
    blk["avi_held_out"] = metrics(h.label, h.avi, label="avi_held_out")
    blk["af_held_out"] = metrics(h.label, h.af_baseline, label="af_held_out")
    # AF bands
    bands = []
    for b, g in d.groupby("af_band"):
        bands.append(dict(band=b, avi=metrics(g.label, g.avi, nboot=500, label=f"avi_{b}"),
                          af=metrics(g.label, g.af_baseline, nboot=500, label=f"af_{b}")))
    blk["af_bands"] = bands
    freq[tag] = blk

# AF-matched, and the doubly-controlled (AF-matched AND consequence-stratified) analysis
for tag, fn in [("genomes_faf95", "frequency_analysis/genomes_fafmax_faf95_max_af_matched.parquet"),
                ("joint_faf95", "frequency_analysis/joint_fafmax_faf95_max_joint_af_matched.parquet")]:
    m = pd.read_parquet(P + fn)
    m["cc"] = m["molecular_consequence"].map(collapse)
    freq[tag]["af_matched"] = metrics(m.label, m.avi, label="af_matched")
    freq[tag]["af_matched_af_baseline"] = metrics(m.label, m.af_baseline, label="af_matched_af")
    num = den = 0.0; rows = []
    for cc, g in m.groupby("cc"):
        yy = g.label.to_numpy()
        if yy.sum() and (yy == 0).sum() and len(g) >= 50:
            a = roc_auc_score(yy, g.avi); w = yy.sum() * (yy == 0).sum()
            num += a * w; den += w
            rows.append(dict(consequence=cc, n=int(len(g)), events=int(yy.sum()),
                             prevalence=round(float(yy.mean()), 4), auroc=round(float(a), 4)))
    freq[tag]["af_matched_consequence_stratified_auroc"] = round(num / den, 4)
    freq[tag]["af_matched_within_consequence"] = rows
    freq[tag]["af_matched_composition"] = {k: int(v) for k, v in m["cc"].value_counts().items()}
R["frequency"] = freq

# ============================================================ 4. four-way comparator comparison
PRED = {"avi": "AlphaGenome AVI", "am_pathogenicity": "AlphaMissense",
        "cadd_phred": "CADD", "revel_score": "REVEL"}
cols = list(PRED)
fw_src = pd.read_parquet(P + "comparators/four_way_matched_variants.parquet")
fw_src["cc"] = fw_src["molecular_consequence"].map(collapse)
fw_src["partition"] = fw_src["chrom"].map(partition)
yf = fw_src["label"].to_numpy()
fw = {c: dict(label=PRED[c], **metrics(yf, fw_src[c], label=PRED[c])) for c in cols}

pairs = []
for i in range(len(cols)):
    for j in range(i + 1, len(cols)):
        a, b = cols[i], cols[j]
        dd, se, pv, ci = delong_pair(yf, fw_src[a], fw_src[b])
        pairs.append(dict(a=PRED[a], b=PRED[b], delta=round(dd, 4), se=round(se, 5),
                          ci=[round(ci[0], 4), round(ci[1], 4)], p=pv))
for k, adj in zip(pairs, holm([x["p"] for x in pairs])):
    k["p_holm"] = float(adj)
corr = fw_src[cols].corr(method="spearman")
corr.to_csv(f"{OUT}/T_predictor_spearman.csv")

# coverage from the full table
cov = {c: dict(nonmissing=int(comp[c].notna().sum()),
               pct=round(float(comp[c].notna().mean()), 4))
       for c in ["avi", "am_pathogenicity", "revel_score", "cadd_phred", "spliceai_max_ds"]}
n_missense = int(comp.cc.eq("Missense").sum())
for c in ["am_pathogenicity", "revel_score"]:
    cov[c]["pct_of_missense"] = round(float(comp.loc[comp.cc.eq("Missense"), c].notna().mean()), 4)

R["four_way"] = dict(n=int(len(fw_src)), events=int(yf.sum()),
                     prevalence=round(float(yf.mean()), 4),
                     composition={k: int(v) for k, v in fw_src["cc"].value_counts().items()},
                     predictors=fw, pairwise=pairs,
                     spearman={f"{PRED[a]} vs {PRED[b]}": round(float(corr.loc[a, b]), 4)
                               for i, a in enumerate(cols) for b in cols[i + 1:]},
                     coverage=cov, n_missense_in_cohort=n_missense)

# four-way restricted to held-out chromosomes
fh = fw_src[fw_src.partition.eq("held_out")]
R["four_way_held_out"] = dict(n=int(len(fh)), events=int(fh.label.sum()),
                              prevalence=round(float(fh.label.mean()), 4),
                              predictors={c: dict(label=PRED[c],
                                                  **metrics(fh.label, fh[c], nboot=1000, label=PRED[c]))
                                          for c in cols})
fhp = []
for b in ["am_pathogenicity", "revel_score", "cadd_phred"]:
    dd, se, pv, ci = delong_pair(fh.label, fh["avi"], fh[b])
    fhp.append(dict(a="AlphaGenome AVI", b=PRED[b], delta=round(dd, 4),
                    ci=[round(ci[0], 4), round(ci[1], 4)], p=pv))
for k, adj in zip(fhp, holm([x["p"] for x in fhp])):
    k["p_holm"] = float(adj)
R["four_way_held_out"]["pairwise_vs_avi"] = fhp

# ============================================================ 5. splice (exploratory, underpowered)
sp = pd.read_parquet(P + "comparators/splice_donor_acceptor_matched_variants.parquet")
ys = sp["label"].to_numpy()
spm = {c: dict(label=lab, **metrics(ys, sp[c], label=lab))
       for c, lab in [("avi", "AlphaGenome AVI"), ("cadd_phred", "CADD"),
                      ("spliceai_max_ds", "SpliceAI max delta score")]}
sppairs = []
for b in ["cadd_phred", "spliceai_max_ds"]:
    dd, se, pv, ci = delong_pair(ys, sp["avi"], sp[b])
    sppairs.append(dict(a="AlphaGenome AVI", b=spm[b]["label"], delta=round(dd, 4),
                        ci=[round(ci[0], 4), round(ci[1], 4)], p=pv))
for k, adj in zip(sppairs, holm([x["p"] for x in sppairs])):
    k["p_holm"] = float(adj)
R["splice_exploratory"] = dict(n=int(len(sp)), events=int(ys.sum()),
                               non_events=int((ys == 0).sum()),
                               prevalence=round(float(ys.mean()), 4),
                               predictors=spm, pairwise_vs_avi=sppairs,
                               interpretation=("Inconclusive: no significant difference detected; "
                                               "equivalence not established."))

# ============================================================ 6. operating points on the held-out set
fpr, tpr, thr = roc_curve(yh, sh)
rows = []
for target in [0.80, 0.90, 0.95]:
    i = int(np.argmin(np.abs(tpr - target))); t = thr[i]
    rows.append(("sensitivity ~%.0f%%" % (target * 100), t))
j = int(np.argmax(tpr - fpr)); rows.append(("Youden index maximum", thr[j]))
for q in [0.90, 0.95]:
    rows.append((f"top {100*(1-q):.0f}% of scores", float(np.quantile(sh, q))))
op = []
for crit, t in rows:
    pred = sh >= t
    tp = int(((pred) & (yh == 1)).sum()); fp = int(((pred) & (yh == 0)).sum())
    fn = int(((~pred) & (yh == 1)).sum()); tn = int(((~pred) & (yh == 0)).sum())
    op.append(dict(criterion=crit, threshold=round(float(t), 3),
                   flagged=int(pred.sum()), flagged_pct=round(float(pred.mean()), 4),
                   sensitivity=round(tp / (tp + fn), 4), specificity=round(tn / (tn + fp), 4),
                   ppv=round(tp / (tp + fp), 4) if tp + fp else None,
                   npv=round(tn / (tn + fn), 4) if tn + fn else None,
                   tp=tp, fp=fp, fn=fn, tn=tn))
pd.DataFrame(op).to_csv(f"{OUT}/T_operating_points_held_out.csv", index=False)
R["operating_points_held_out"] = op

t90 = float(np.quantile(sh, 0.90))
R["enrichment_held_out"] = dict(
    top_decile_event_rate=round(float(yh[sh >= t90].mean()), 4),
    cohort_prevalence=round(float(yh.mean()), 4),
    fold=round(float(yh[sh >= t90].mean() / yh.mean()), 2),
    lower_half_event_rate=round(float(yh[sh < float(np.quantile(sh, 0.5))].mean()), 5))

# ============================================================ 7. sensitivity analyses (held-out)
sens = {}
m1 = (held.reversions == 0) & (held.class_changes == 1)
sens["single_change_no_reversion"] = metrics(held[m1].label, held[m1].avi, nboot=1000)
pa = held.cc.isin(["Missense", "Nonsense", "Splice donor", "Splice acceptor", "Start lost", "Stop lost"])
sens["protein_altering"] = metrics(held[pa].label, held[pa].avi, nboot=1000)
gs = held.sample(frac=1.0, random_state=7).drop_duplicates("gene")
sens["one_variant_per_gene"] = metrics(gs.label, gs.avi, nboot=1000)
top10 = set(held.gene.value_counts().head(10).index)
ex = held[~held.gene.isin(top10)]
sens["excluding_top10_genes"] = metrics(ex.label, ex.avi, nboot=1000)
R["sensitivity_held_out"] = sens

# gene-level on held-out
gl = []
for gname, g in held.groupby("gene"):
    d = dict(gene=gname, n=int(len(g)), events=int(g.label.sum()),
             prevalence=round(float(g.label.mean()), 4), median_avi=round(float(g.avi.median()), 4))
    if g.label.sum() and (g.label == 0).sum():
        d["auroc"] = round(float(roc_auc_score(g.label, g.avi)), 4)
    gl.append(d)
genes = pd.DataFrame(gl).sort_values("n", ascending=False)
genes.to_csv(f"{OUT}/T_gene_level_held_out.csv", index=False)
ev = genes.dropna(subset=["auroc"])
R["genes_held_out"] = dict(n_genes=int(len(genes)),
                           n_ge20=int((genes.n >= 20).sum()),
                           n_ge20_evaluable=int(((genes.n >= 20) & genes.auroc.notna()).sum()),
                           median_auroc_ge20=round(float(ev.loc[ev.n >= 20, "auroc"].median()), 4),
                           q1=round(float(ev.loc[ev.n >= 20, "auroc"].quantile(.25)), 4),
                           q3=round(float(ev.loc[ev.n >= 20, "auroc"].quantile(.75)), 4))

# ============================================================ 8. ROC and precision-recall curve coordinates
np.save(f"{OUT}/roc_held_out.npy", np.vstack(roc_curve(yh, sh)[:2]))
pre, rec, _ = precision_recall_curve(yh, sh)
np.save(f"{OUT}/pr_held_out.npy", np.vstack([rec, pre]))
curves = {}
for c in cols:
    f2, t2, _ = roc_curve(yf, fw_src[c]); p2, r2, _ = precision_recall_curve(yf, fw_src[c])
    curves[c] = dict(fpr=f2.tolist(), tpr=t2.tolist(), rec=r2.tolist(), pre=p2.tolist())
json.dump(curves, open(f"{OUT}/fourway_curves.json", "w"))
comp[["cc", "gene", "label", "avi", "partition", "chrom", "reclassification_year",
      "am_pathogenicity", "revel_score", "cadd_phred", "spliceai_max_ds"]].to_parquet(
    f"{OUT}/analysis_frame.parquet")
fw_src.to_parquet(f"{OUT}/fourway_frame.parquet")
for tag, fn in [("genomes", "frequency_analysis/genomes_fafmax_faf95_max_variants.parquet")]:
    pd.read_parquet(P + fn).to_parquet(f"{OUT}/freq_{tag}.parquet")
pd.read_parquet(P + "frequency_analysis/genomes_fafmax_faf95_max_af_matched.parquet").to_parquet(
    f"{OUT}/freq_genomes_matched.parquet")

json.dump(R, open(f"{OUT}/results_v2.json", "w"), indent=2, default=float)

# ---------------------------------------------------------------- report
print("COHORT:", R["cohort"])
print("\nCHROMOSOME PARTITIONS")
print(part[["partition", "scope", "n", "events", "prevalence", "auroc", "auroc_lo", "auroc_hi",
            "auprc", "auprc_baseline", "auprc_lift"]].to_string(index=False))
print("\nPRIMARY (held out):", json.dumps(R["primary_held_out"], default=float))
print("\nCONSEQUENCE (with prevalence baseline)")
print(cons[["consequence", "n", "events", "prevalence", "auroc", "auroc_lo", "auroc_hi",
            "auprc", "auprc_baseline", "auprc_lift", "main_table"]].to_string(index=False))
print("\nCONSEQUENCE PRIOR:", json.dumps(R["consequence_prior"], default=float))
print("CONSEQUENCE-STRATIFIED AUROC:", R["consequence_stratified_auroc"])
print("\nFREQUENCY")
for tag, b in freq.items():
    print(f"  [{tag}] n={b['n']} avi={b['avi_overall']['auroc']} af={b['af_overall']['auroc']} "
          f"diff={b['avi_vs_af']['difference']} | held-out avi={b['avi_held_out']['auroc']} "
          f"af={b['af_held_out']['auroc']} | AF-matched avi={b['af_matched']['auroc']} "
          f"af={b['af_matched_af_baseline']['auroc']} | matched+cons-stratified="
          f"{b['af_matched_consequence_stratified_auroc']}")
print("\nFOUR-WAY n=%d events=%d prev=%.4f comp=%s" % (
    R["four_way"]["n"], R["four_way"]["events"], R["four_way"]["prevalence"],
    R["four_way"]["composition"]))
for c in cols:
    v = fw[c]
    print(f"   {v['label']:18s} AUROC {v['auroc']:.4f} ({v['auroc_lo']:.4f}-{v['auroc_hi']:.4f})"
          f"  AUPRC {v['auprc']:.4f} (baseline {v['auprc_baseline']:.4f}, lift {v['auprc_lift']:+.4f})")
for k in pairs:
    print(f"   {k['a']} - {k['b']}: d={k['delta']:+.4f} CI{k['ci']} p={k['p']:.3g} holm={k['p_holm']:.3g}")
print("\nFOUR-WAY held-out n=%d" % R["four_way_held_out"]["n"])
for c in cols:
    v = R["four_way_held_out"]["predictors"][c]
    print(f"   {v['label']:18s} AUROC {v['auroc']:.4f} ({v['auroc_lo']:.4f}-{v['auroc_hi']:.4f})")
for k in fhp:
    print(f"   AVI - {k['b']}: d={k['delta']:+.4f} CI{k['ci']} holm={k['p_holm']:.3g}")
print("\nSPLICE:", json.dumps(R["splice_exploratory"], indent=1, default=float)[:1200])
print("\nCOVERAGE:", json.dumps(cov, default=float))
print("\nOPERATING POINTS (held out)"); print(pd.DataFrame(op).to_string(index=False))
print("\nENRICHMENT:", R["enrichment_held_out"])
print("\nSENSITIVITY:", json.dumps(sens, default=float))
print("\nGENES:", R["genes_held_out"])
