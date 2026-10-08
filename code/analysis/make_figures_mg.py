#!/usr/bin/env python3
"""Redesigned figures for the Mammalian Genome submission.

Springer artwork rules applied: 174 mm single-column width, height <= 234 mm,
Arial 8-11 pt, lowercase bold panel letters, no titles inside the artwork,
RGB, 600 dpi TIFF (combination art) plus a 300 dpi PNG preview, and every
series distinguished by line style or marker shape as well as colour.

Curves are recomputed from the Path B source parquet files; every reported
statistic (AUROC, CI, p) is read from results_v2.json or the supplementary
CSVs so the figures cannot drift from the manuscript.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from sklearn.metrics import roc_curve, precision_recall_curve

ROOT = Path(__file__).resolve().parents[2]          # project9_ag_vus
PB = ROOT / "manuscript" / "path_b"
V2 = ROOT / "manuscript_v2"
OUT = ROOT / "mammalian_submission" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

R = json.load(open(V2 / "analysis" / "results_v2.json"))
df = pd.read_parquet(PB / "comparators" / "full_primary_source_comparator_scores.parquet")
fw = pd.read_parquet(PB / "comparators" / "four_way_matched_variants.parquet")
fq = pd.read_parquet(PB / "frequency_analysis" / "genomes_fafmax_faf95_max_variants.parquet")
genes = pd.read_csv(V2 / "supplementary" / "S7_gene_level_held_out.csv")
nc_cls = pd.read_csv(V2 / "supplementary" / "S9b_noncoding_by_class.csv")

HELD = {str(c) for c in R["partition_chromosomes"]["held_out"]}
df["held"] = df.chrom.astype(str).isin(HELD)
fw["held"] = fw.chrom.astype(str).isin(HELD)
fq["held"] = fq.chrom.astype(str).isin(HELD)

# Okabe-Ito palette; each predictor also gets its own line style and marker
BLUE, VERM, GREEN, ORANGE, SKY, PURPLE, GREY = (
    "#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#7F7F7F")
PRED = {  # key: (label, colour, linestyle, marker)
    "avi": ("AlphaGenome AVI", BLUE, "-", "o"),
    "revel_score": ("REVEL", GREEN, (0, (5, 1.5)), "s"),
    "am_pathogenicity": ("AlphaMissense", ORANGE, (0, (1, 1.2)), "D"),
    "cadd_phred": ("CADD", VERM, (0, (4, 1.5, 1, 1.5)), "^"),
    "spliceai_max_ds": ("SpliceAI", PURPLE, (0, (5, 1.5)), "v"),
}

MM = 1 / 25.4
WIDTH = 174 * MM
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "legend.frameon": False, "legend.handlelength": 2.6,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 1.4, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.facecolor": "white", "figure.facecolor": "white",
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.default": "regular",
})


def save(fig, n):
    from PIL import Image
    tmp = OUT / f"Fig{n}_rgba.png"
    fig.savefig(tmp, dpi=600)
    Image.open(tmp).convert("RGB").save(OUT / f"Fig{n}.tif", dpi=(600, 600), compression="tiff_lzw")
    tmp.unlink()
    fig.savefig(OUT / f"Fig{n}.png", dpi=300)
    plt.close(fig)
    print("wrote", f"Fig{n}")


def letter(ax, s, x=-0.17, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=11, fontweight="bold",
            va="bottom", ha="left")


def roc_axes(ax):
    ax.plot([0, 1], [0, 1], color="0.6", lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.005); ax.set_aspect("equal")
    ax.set_xticks([0, .25, .5, .75, 1]); ax.set_yticks([0, .25, .5, .75, 1])
    ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1"])
    ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1"])
    ax.set_xlabel("False-positive rate (1 − specificity)")
    ax.set_ylabel("True-positive rate (sensitivity)")


def pr_axes(ax, prev):
    ax.axhline(prev, color="0.6", lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1.005); ax.set_aspect("equal")
    ax.set_xticks([0, .25, .5, .75, 1]); ax.set_yticks([0, .25, .5, .75, 1])
    ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1"])
    ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1"])
    ax.set_xlabel("Recall (sensitivity)"); ax.set_ylabel("Precision (PPV)")


def curve(ax, y, s, key, label, kind="roc", lw=1.4):
    _, col, ls, _ = PRED[key]
    m = ~np.isnan(s)
    if kind == "roc":
        x1, y1, _ = roc_curve(y[m], s[m])
    else:
        p, r, _ = precision_recall_curve(y[m], s[m]); x1, y1 = r, p
    ax.plot(x1, y1, color=col, ls=ls, lw=lw, label=label, zorder=3)


def ci_row(ax, y, est, lo, hi, color, marker, mfc=None, ms=4.5):
    ax.errorbar(est, y, xerr=[[est - lo], [hi - est]], fmt=marker, color=color,
                mfc=mfc if mfc else color, mec=color, ms=ms, elinewidth=1.0,
                capsize=2.2, capthick=1.0, zorder=3)


# =================================================================== Fig 1
held = df[df.held]
y, s = held.label.to_numpy(), held.avi.to_numpy()
PH = R["primary_held_out"]
fig = plt.figure(figsize=(WIDTH, 150 * MM), layout="constrained")
gs = fig.add_gridspec(2, 2, )

ax = fig.add_subplot(gs[0, 0])
roc_axes(ax)
fpr, tpr, thr = roc_curve(y, s)
ax.plot(fpr, tpr, color=BLUE, lw=1.6, zorder=3)
op = R["operating_points_held_out"]
you = [o for o in op if "youden" in o["criterion"].lower()][0]
ax.plot(1 - you["specificity"], you["sensitivity"], "o", ms=5.5, mfc="white",
        mec=VERM, mew=1.4, zorder=4)
ax.annotate(f"Youden cut-off {you['threshold']:.3f}\nsens. {you['sensitivity']:.2f}, "
            f"spec. {you['specificity']:.2f}",
            xy=(1 - you["specificity"], you["sensitivity"]), xytext=(0.40, 0.22),
            fontsize=8, color=VERM, arrowprops=dict(arrowstyle="-", color=VERM, lw=0.8))
ax.text(0.97, 0.04, f"AUROC {PH['auroc']:.3f}\n95% CI {PH['auroc_lo']:.3f}–{PH['auroc_hi']:.3f}",
        transform=ax.transAxes, ha="right", va="bottom")
letter(ax, "a", x=-0.30)

ax = fig.add_subplot(gs[0, 1])
pr_axes(ax, PH["prevalence"])
p, r, _ = precision_recall_curve(y, s)
ax.plot(r, p, color=BLUE, lw=1.6, zorder=3)
ax.text(0.03, PH["prevalence"] + 0.025, f"Prevalence {PH['prevalence']:.3f}", color="0.35")
ax.text(0.03, 0.30, f"AUPRC {PH['auprc']:.3f}\n95% CI {PH['auprc_lo']:.3f}–{PH['auprc_hi']:.3f}\n"
        f"Lift +{PH['auprc_lift']:.3f}", va="bottom")
letter(ax, "b", x=-0.30)

ax = fig.add_subplot(gs[1, 0])
bins = np.linspace(np.floor(s.min() * 4) / 4, np.ceil(s.max() * 4) / 4, 41)
ax.hist(s[y == 0], bins=bins, density=True, color=SKY, alpha=0.75, label="Later benign/LB",
        edgecolor="white", lw=0.3)
ax.hist(s[y == 1], bins=bins, density=True, histtype="step", color=VERM, lw=1.4,
        hatch="////", label="Later pathogenic/LP")
ax.axvline(you["threshold"], color="0.25", lw=0.9, ls=(0, (3, 2)))
ax.text(you["threshold"] - 0.05, ax.get_ylim()[1] * 0.97, "Youden cut-off", fontsize=8, va="top", ha="right")
ax.set_xlabel("AVI score"); ax.set_ylabel("Density")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=8, handlelength=1.4,
          columnspacing=1.0)
letter(ax, "c", x=-0.30)

ax = fig.add_subplot(gs[1, 1])
pr_ = pd.DataFrame(R["chromosome_partitions"])
order = [("held_out", "Held out"), ("validation", "Validation"), ("training", "Training")]
for i, (pt, nm) in enumerate(order):
    for sc, dy, col, mk, mfc in [("all", -0.15, BLUE, "o", None), ("missense", 0.15, VERM, "s", "white")]:
        rr = pr_[(pr_.partition == pt) & (pr_.scope == sc)].iloc[0]
        ci_row(ax, i + dy, rr.auroc, rr.auroc_lo, rr.auroc_hi, col, mk, mfc)
ax.set_yticks(range(3)); ax.set_yticklabels([nm for _, nm in order])
ax.invert_yaxis(); ax.set_ylim(2.6, -0.6)
ax.set_xlim(0.88, 0.97); ax.set_xticks([0.88, 0.90, 0.92, 0.94, 0.96])
ax.set_xlabel("AUROC (95% CI)"); ax.set_ylabel("Chromosome partition")
ax.legend(handles=[Line2D([], [], color=BLUE, marker="o", ls="none", label="All variants"),
                   Line2D([], [], color=VERM, marker="s", mfc="white", ls="none", label="Missense only")],
          loc="lower right", fontsize=8, handlelength=1, borderaxespad=0.2)
letter(ax, "d", x=-0.30)
save(fig, 1)

# =================================================================== Fig 2
bc = pd.DataFrame(R["by_consequence"])
bc = bc[bc.main_table == True].sort_values("auroc")  # noqa: E712
CP = R["consequence_prior"]
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 78 * MM), layout="constrained")
ax = axs[0]
for i, rr in enumerate(bc.itertuples()):
    ci_row(ax, i, rr.auroc, rr.auroc_lo, rr.auroc_hi, BLUE, "o")
ax.axvline(0.5, color="0.6", lw=0.8)
ax.axvline(CP["auroc"], color="0.3", lw=0.9, ls=(0, (1, 1.5)))
ax.text(CP["auroc"] + 0.008, len(bc) - 0.45, f"Consequence\nclass alone\n{CP['auroc']:.3f}",
        fontsize=8, va="top", color="0.25")
ax.set_yticks(range(len(bc)))
ax.set_yticklabels([f"{c}\n{n:,} variants; {p*100:.1f}% P/LP"
                    for c, n, p in zip(bc.consequence, bc.n, bc.prevalence)])
ax.set_xlim(0.45, 1.0); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_ylim(-0.6, len(bc) - 0.4)
ax.set_xlabel("AUROC (95% CI)")
letter(ax, "a", x=-0.62, y=1.04)

ax = axs[1]
for i, rr in enumerate(bc.itertuples()):
    ax.plot([rr.auprc_baseline, rr.auprc], [i, i], color="0.55", lw=1.2, zorder=1)
    ax.plot(rr.auprc_baseline, i, "o", mfc="white", mec="0.25", ms=5, zorder=2)
    ax.plot(rr.auprc, i, "D", color=ORANGE, mec="0.15", mew=0.5, ms=5.5, zorder=3)
    ax.text(1.03, i, f"+{rr.auprc_lift:.3f}", va="center", ha="left", fontsize=8)
ax.text(1.03, len(bc) - 0.4, "Lift", fontsize=8, fontweight="bold", ha="left", va="bottom")
ax.set_yticks(range(len(bc))); ax.set_yticklabels([])
ax.set_xlim(0, 1.0); ax.set_ylim(-0.6, len(bc) - 0.4)
ax.set_xlabel("Precision–recall value")
ax.legend(handles=[Line2D([], [], marker="o", mfc="white", mec="0.25", ls="none", label="Prevalence baseline"),
                   Line2D([], [], marker="D", color=ORANGE, mec="0.15", ls="none", label="AUPRC")],
          loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=2, fontsize=8, handlelength=1,
          columnspacing=1.2)
letter(ax, "b", x=-0.06, y=1.04)
save(fig, 2)

# =================================================================== Fig 3
G = R["frequency"]["genomes_faf95"]
fig = plt.figure(figsize=(WIDTH, 160 * MM), layout="constrained")
top, bot = fig.subfigures(2, 1, height_ratios=[1.15, 1])
axa, axb = top.subplots(1, 2, gridspec_kw=dict(width_ratios=[1, 1.1]))
ax = axa
roc_axes(ax)
fh = fq[fq.held]
curve(ax, fh.label.to_numpy(), fh.avi.to_numpy(), "avi", f"AVI  {G['avi_held_out']['auroc']:.3f}")
fpr, tpr, _ = roc_curve(fh.label, fh.af_baseline)
ax.plot(fpr, tpr, color=ORANGE, ls=(0, (5, 1.5)), lw=1.4,
        label=f"gnomAD FAF  {G['af_held_out']['auroc']:.3f}", zorder=2)
ax.legend(title="AUROC", loc="lower right", fontsize=8, title_fontsize=8)
letter(ax, "a", x=-0.30)

ax = axb
band_lbl = {"absent": "Absent", "observed_<1e-5": r"$<10^{-5}$", "1e-5_to_<1e-4": r"$10^{-5}$–$10^{-4}$",
            "1e-4_to_<1e-3": r"$10^{-4}$–$10^{-3}$", ">=1e-3": r"$\geq10^{-3}$"}
bands = {b["band"]: b for b in G["af_bands"]}
for i, k in enumerate(band_lbl):
    b = bands[k]
    ci_row(ax, i - 0.15, b["avi"]["auroc"], b["avi"]["auroc_lo"], b["avi"]["auroc_hi"], BLUE, "o")
    if k == "absent":
        ax.plot(0.5, i + 0.15, marker="s", mfc="white", mec=ORANGE, ms=4.5, ls="none")
    else:
        ci_row(ax, i + 0.15, b["af"]["auroc"], b["af"]["auroc_lo"], b["af"]["auroc_hi"],
               ORANGE, "s", "white")
ax.set_yticks(range(5))
ax.set_yticklabels([f"{band_lbl[k]}\n{bands[k]['avi']['n']:,} ({bands[k]['avi']['events']:,} P/LP)"
                    for k in band_lbl])
ax.invert_yaxis(); ax.axvline(0.5, color="0.6", lw=0.8)
ax.set_xlim(0.4, 1.0); ax.set_xlabel("AUROC (95% CI)")
ax.legend(handles=[Line2D([], [], color=BLUE, marker="o", ls="none", label="AVI"),
                   Line2D([], [], color=ORANGE, marker="s", mfc="white", ls="none", label="FAF baseline")],
          loc="upper center", bbox_to_anchor=(0.58, 1.0), fontsize=8, handlelength=1, borderaxespad=0.2)
letter(ax, "b", x=-0.62)

ax = bot.subplots()
am = G["af_matched"]
from sklearn.metrics import roc_auc_score
ao = G["avi_overall"]
ladder = [
    ("No control, whole cohort", "42,383", roc_auc_score(df.label, df.avi), None, None),
    ("Consequence-stratified, whole cohort", "42,383", R["consequence_stratified_auroc"], None, None),
    ("No control, FAF-evaluable set", f"{ao['n']:,}", ao["auroc"], ao["auroc_lo"], ao["auroc_hi"]),
    ("Allele-frequency matched", f"{am['n']:,}", am["auroc"], am["auroc_lo"], am["auroc_hi"]),
    ("Frequency matched and consequence-stratified", f"{am['n']:,}",
     G["af_matched_consequence_stratified_auroc"], None, None),
]
for i, (nm, n, est, lo, hi) in enumerate(ladder):
    if lo is not None:
        ci_row(ax, i, est, lo, hi, GREEN, "o", ms=6)
    else:
        ax.plot(est, i, "o", color=GREEN, ms=6, zorder=3)
    ax.plot([0.5, est], [i, i], color=GREEN, lw=0.8, alpha=0.5, zorder=1)
    ax.text(est + 0.006, i - 0.18, f"{est:.3f}", va="bottom", fontsize=8)
ax.set_yticks(range(len(ladder)))
ax.set_yticklabels([f"{nm} (n = {n})" for nm, n, *_ in ladder])
ax.invert_yaxis(); ax.set_xlim(0.5, 1.0); ax.set_ylim(len(ladder) - 0.4, -0.6)
ax.set_xlabel("AVI AUROC (axis starts at chance, 0.5)")
letter(ax, "c", x=-0.62, y=1.0)
save(fig, 3)

# =================================================================== Fig 4
FW = R["four_way"]; FH = R["four_way_held_out"]
keys = ["revel_score", "am_pathogenicity", "avi", "cadd_phred"]
fig = plt.figure(figsize=(WIDTH, 160 * MM), layout="constrained")
gs = fig.add_gridspec(2, 2, )
yv = fw.label.to_numpy()
ax = fig.add_subplot(gs[0, 0]); roc_axes(ax)
for k in keys:
    curve(ax, yv, fw[k].to_numpy(), k, f"{PRED[k][0]}  {FW['predictors'][k]['auroc']:.3f}",
          lw=1.8 if k == "avi" else 1.3)
ax.legend(title="AUROC", loc="lower right", fontsize=8, title_fontsize=8)
letter(ax, "a", x=-0.30)
ax = fig.add_subplot(gs[0, 1]); pr_axes(ax, FW["prevalence"])
for k in keys:
    curve(ax, yv, fw[k].to_numpy(), k, f"{PRED[k][0]}  {FW['predictors'][k]['auprc']:.3f}", "pr",
          lw=1.8 if k == "avi" else 1.3)
ax.legend(title="AUPRC", loc="lower left", fontsize=8, title_fontsize=8)
ax.text(0.97, FW["prevalence"] - 0.03, f"Prevalence {FW['prevalence']:.3f}", ha="right", va="top", color="0.35")
letter(ax, "b", x=-0.30)

ax = fig.add_subplot(gs[1, 0])
for i, k in enumerate(keys):
    _, col, _, mk = PRED[k]
    a, h = FW["predictors"][k], FH["predictors"][k]
    ci_row(ax, i - 0.16, a["auroc"], a["auroc_lo"], a["auroc_hi"], col, mk)
    ci_row(ax, i + 0.16, h["auroc"], h["auroc_lo"], h["auroc_hi"], col, mk, "white")
ax.set_yticks(range(4)); ax.set_yticklabels([PRED[k][0] for k in keys])
ax.invert_yaxis(); ax.set_xlim(0.80, 0.94); ax.set_xlabel("AUROC (95% CI)")
ax.legend(handles=[Line2D([], [], color="0.2", marker="o", ls="none",
                          label=f"All matched (n = {FW['n']:,})"),
                   Line2D([], [], color="0.2", marker="o", mfc="white", ls="none",
                          label=f"Held-out chr. (n = {FH['n']:,})")],
          loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, fontsize=8, handlelength=1)
letter(ax, "c", x=-0.42, y=1.12)

ax = fig.add_subplot(gs[1, 1])
pw_all = {p["b"]: p for p in FW["pairwise"] if p["a"] == "AlphaGenome AVI"}
pw_h = {p["b"]: p for p in FH["pairwise_vs_avi"]}
comps = [("REVEL", "revel_score"), ("AlphaMissense", "am_pathogenicity"), ("CADD", "cadd_phred")]


def stars(p):
    return "ns" if p >= 0.05 else "*" if p >= 0.01 else "**" if p >= 1e-3 else "***"


for i, (nm, k) in enumerate(comps):
    col, mk = PRED[k][1], PRED[k][3]
    for dy, pw, mfc in [(-0.16, pw_all[nm], None), (0.16, pw_h[nm], "white")]:
        # plot comparator minus AVI so that bars right of zero favour the comparator
        d, lo, hi = -pw["delta"], -pw["ci"][1], -pw["ci"][0]
        ci_row(ax, i + dy, d, lo, hi, col, mk, mfc)
        ax.text(0.10, i + dy, stars(pw["p_holm"]), va="center", fontsize=8)
ax.axvline(0, color="0.4", lw=0.8)
ax.set_yticks(range(3)); ax.set_yticklabels([f"{nm} − AVI" for nm, _ in comps])
ax.invert_yaxis(); ax.set_xlim(-0.10, 0.115)
ax.set_xlabel("Difference in AUROC (95% CI)")
ax.text(0.02, -0.75, "Favours comparator →", fontsize=8, color="0.3")
ax.text(-0.02, -0.75, "← Favours AVI", fontsize=8, color="0.3", ha="right")
ax.set_ylim(2.5, -0.95)
letter(ax, "d", x=-0.42, y=1.12)
save(fig, 4)

# =================================================================== Fig 5
NC = R["noncoding_comparison"]
ncd = df[df.consequence_class.isin(NC["classes_included"])]
yv = ncd.label.to_numpy()
ks = ["avi", "spliceai_max_ds", "cadd_phred"]
fig = plt.figure(figsize=(WIDTH, 160 * MM), layout="constrained")
gs = fig.add_gridspec(2, 2, )
ax = fig.add_subplot(gs[0, 0]); roc_axes(ax)
for k in ks:
    curve(ax, yv, ncd[k].to_numpy(), k, f"{PRED[k][0]}  {NC['pooled'][k]['auroc']:.3f}",
          lw=1.8 if k == "avi" else 1.3)
ax.legend(title="AUROC", loc="lower right", fontsize=8, title_fontsize=8)
letter(ax, "a", x=-0.30)
ax = fig.add_subplot(gs[0, 1]); pr_axes(ax, NC["prevalence"])
for k in ks:
    curve(ax, yv, ncd[k].to_numpy(), k, f"{PRED[k][0]}  {NC['pooled'][k]['auprc']:.3f}", "pr",
          lw=1.8 if k == "avi" else 1.3)
ax.legend(title="AUPRC", loc="lower left", bbox_to_anchor=(0.0, 0.07), fontsize=8, title_fontsize=8)
ax.text(0.03, NC["prevalence"] + 0.012, f"Prevalence {NC['prevalence']:.3f}", ha="left", color="0.35")
letter(ax, "b", x=-0.30)

ax = fig.add_subplot(gs[1, 0])
cls = nc_cls[nc_cls.status.isna()].reset_index(drop=True)
for i, rr in cls.iterrows():
    ci_row(ax, i - 0.16, rr.avi_auroc, rr.avi_auroc_lo, rr.avi_auroc_hi, BLUE, "o")
    ci_row(ax, i + 0.16, rr.cadd_auroc, rr.cadd_auroc_lo, rr.cadd_auroc_hi, VERM, "^", "white")
    ax.text(1.03, i, stars(rr.p), va="center", fontsize=8)
ax.set_yticks(range(len(cls)))
ax.set_yticklabels([f"{c.replace(chr(39), '′')}\n{n:,} ({e} P/LP)" for c, n, e in
                    zip(cls.consequence, cls.n, cls.events)])
ax.invert_yaxis(); ax.set_xlim(0.5, 1.08); ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
ax.set_xlabel("AUROC (95% CI)")
ax.legend(handles=[Line2D([], [], color=BLUE, marker="o", ls="none", label="AVI"),
                   Line2D([], [], color=VERM, marker="^", mfc="white", ls="none", label="CADD")],
          loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=8, handlelength=1)
letter(ax, "c", x=-0.42)

ax = fig.add_subplot(gs[1, 1])
cov = [("avi", 1.0), ("cadd_phred", NC["cadd_coverage"]), ("spliceai_max_ds", NC["spliceai_coverage"]),
       ("revel_score", NC["revel_coverage"]), ("am_pathogenicity", NC["alphamissense_coverage"])]
hatches = ["", "", "", "////", "////"]
for i, ((k, v), h) in enumerate(zip(cov, hatches)):
    col = PRED[k][1] if not h else "white"
    ax.barh(i, v * 100, color=col, edgecolor=PRED[k][1] if h else "0.15", hatch=h, lw=0.8, height=0.62)
    ax.text(v * 100 + 2, i, f"{v*100:.1f}%" if v >= 0.001 else f"{v*100:.2f}%", va="center", fontsize=8)
ax.set_yticks(range(len(cov))); ax.set_yticklabels([PRED[k][0] for k, _ in cov])
ax.invert_yaxis(); ax.set_xlim(0, 118); ax.set_xticks([0, 25, 50, 75, 100])
ax.set_xlabel(f"Variants scored (% of {NC['n']:,})")
letter(ax, "d", x=-0.42)
save(fig, 5)

# =================================================================== Fig 6
GH = R["genes_held_out"]
g = genes[(genes.n >= 20) & genes.auroc.notna()]
fig, axs = plt.subplots(1, 2, figsize=(WIDTH, 80 * MM), layout="constrained",
                        gridspec_kw=dict(width_ratios=[1, 1.15]))
ax = axs[0]
sz = 8 + 40 * (g.events / g.events.max())
ax.scatter(g.n, g.auroc, s=sz, facecolor=PURPLE, edgecolor="0.2", lw=0.4, alpha=0.8, zorder=3)
ax.set_xscale("log")
ax.axhline(GH["median_auroc_ge20"], color=VERM, lw=1.0, ls=(0, (4, 2)))
ax.axhline(0.5, color="0.6", lw=0.8, ls=(0, (1, 1.5)))
ax.legend(handles=[Line2D([], [], color=VERM, ls=(0, (4, 2)), label=f"Median {GH['median_auroc_ge20']:.3f}"),
                   Line2D([], [], color="0.6", ls=(0, (1, 1.5)), label="Chance"),
                   Line2D([], [], marker="o", ls="none", mfc=PURPLE, mec="0.2", label="Gene")],
          loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=8, handlelength=2,
          columnspacing=1.0)
from matplotlib.ticker import FixedLocator, NullFormatter
ax.xaxis.set_major_locator(FixedLocator([20, 50, 100, 200, 500]))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xticklabels(["20", "50", "100", "200", "500"])
ax.set_xlabel("Variants in gene (log scale)"); ax.set_ylabel("Within-gene AUROC")
ax.set_ylim(0.3, 1.04)
letter(ax, "a", x=-0.22)

ax = axs[1]
SH = R["sensitivity_held_out"]
rows = [("Primary (held-out)", PH, "o", BLUE, None),
        ("Single change, no reversion", SH["single_change_no_reversion"], "o", BLUE, None),
        ("Excluding 10 largest genes", SH["excluding_top10_genes"], "o", BLUE, None),
        ("One variant per gene", SH["one_variant_per_gene"], "o", BLUE, None),
        ("Protein-altering only", SH["protein_altering"], "o", BLUE, None),
        ("Consequence class alone", CP, "s", GREY, "white")]
for i, (nm, d, mk, col, mfc) in enumerate(rows):
    ci_row(ax, i, d["auroc"], d["auroc_lo"], d["auroc_hi"], col, mk, mfc)
ax.axvline(PH["auroc"], color=VERM, lw=0.9, ls=(0, (4, 2)))
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([f"{nm}\n(n = {d['n']:,})" for nm, d, *_ in rows])
ax.invert_yaxis(); ax.set_xlim(0.6, 1.0); ax.set_xlabel("AUROC (95% CI)")
letter(ax, "b", x=-0.62)
save(fig, 6)
