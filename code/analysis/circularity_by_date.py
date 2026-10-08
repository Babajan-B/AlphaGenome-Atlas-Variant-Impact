#!/usr/bin/env python3
"""Label-circularity check for the four-predictor comparison.

ClinVar submitters could use AlphaMissense as PP3/BP4 evidence only after its
public release (Cheng et al., Science, 19 September 2023). If AVI's missense
deficit came from that circularity, it should vanish among variants whose first
reclassification predates the release. Same DeLong, Holm and outcome-stratified
bootstrap code and seed as manuscript_v2/analysis/stats_v2.py.
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
src = (ROOT / "manuscript_v2" / "analysis" / "stats_v2.py").read_text()
helpers = src[:src.index('R = {"seed"')]
helpers = helpers.replace("os.makedirs(OUT, exist_ok=True)", "")
ns = {}
exec(compile(helpers, "stats_v2_helpers", "exec"), ns)
metrics, delong_pair, holm = ns["metrics"], ns["delong_pair"], ns["holm"]

fw = pd.read_parquet(ROOT / "manuscript" / "path_b" / "comparators" / "four_way_matched_variants.parquet")
fw["date"] = pd.to_datetime(fw.first_reclassification_date.astype(str), format="%Y%m%d")
CUT = pd.Timestamp("2023-09-19")
PRED = {"avi": "AlphaGenome AVI", "am_pathogenicity": "AlphaMissense",
        "revel_score": "REVEL", "cadd_phred": "CADD"}

out = {"cutoff": str(CUT.date()), "periods": {}}
for name, sub in [("before_alphamissense", fw[fw.date < CUT]),
                  ("after_alphamissense", fw[fw.date >= CUT])]:
    y = sub.label.to_numpy()
    block = {"n": int(len(sub)), "events": int(y.sum()),
             "first_date": str(sub.date.min().date()), "last_date": str(sub.date.max().date()),
             "predictors": {k: metrics(y, sub[k].to_numpy()) for k in PRED}}
    pw = []
    for k in ["am_pathogenicity", "revel_score", "cadd_phred"]:
        d, se, p, ci = delong_pair(y, sub["avi"], sub[k])
        pw.append(dict(comparator=PRED[k], delta=round(d, 4), ci=[round(ci[0], 4), round(ci[1], 4)], p=p))
    for row, ph in zip(pw, holm([r["p"] for r in pw])):
        row["p_holm"] = float(ph)
    block["avi_minus_comparator"] = pw
    out["periods"][name] = block

(HERE / "circularity_by_date.json").write_text(json.dumps(out, indent=2))
for name, b in out["periods"].items():
    print(name, b["n"], b["events"], b["first_date"], b["last_date"])
    for k, m in b["predictors"].items():
        print(f"   {PRED[k]:16s} {m['auroc']:.3f} ({m['auroc_lo']:.3f}-{m['auroc_hi']:.3f})  AUPRC {m['auprc']:.3f}")
    for r in b["avi_minus_comparator"]:
        print(f"   AVI-{r['comparator']:14s} {r['delta']:+.3f} ({r['ci'][0]:+.3f} to {r['ci'][1]:+.3f}) p={r['p']:.2g} holm={r['p_holm']:.2g}")
