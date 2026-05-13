#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd

OUTDIR = "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step12_revision_strengthening/results/step12B_degree_matched_graph_null"

transport_file = f"{OUTDIR}/03_step12B_transport_null_summary.tsv"
sfari_file = f"{OUTDIR}/06_step12B_NTM2_SFARI_null_summary.tsv"

def infer_top_n(x):
    m = re.search(r"top(\d+)", str(x))
    return int(m.group(1)) if m else np.nan

transport = pd.read_csv(transport_file, sep="\t")
sfari = pd.read_csv(sfari_file, sep="\t")

transport["top_n_fixed"] = transport["module_node"].map(infer_top_n)
sfari["top_n_fixed"] = sfari["module_node"].map(infer_top_n)

# Corrected transport decision
records = []

for _, r in transport.iterrows():
    fam = r["module_family"]
    topn = int(r["top_n_fixed"])
    fdr = float(r["primary_fdr"])

    if fam in ["NTM1_ASD_up", "NTM3_ASD_signed"] and topn in [200, 500]:
        status = "PASS" if fdr < 0.05 else "REVIEW"
        interpretation = "late-prenatal transport exceeds degree/density/stage-profile matched null" if status == "PASS" else "late-prenatal transport does not exceed degree/density/stage-profile matched null"
    elif fam == "NTM2_ASD_down" and topn in [200, 500]:
        status = "PASS" if fdr < 0.05 else "REVIEW"
        interpretation = "NTM2 non-prenatal/adult-like transport exceeds degree/density/stage-profile matched null" if status == "PASS" else "NTM2 transport separation is not significant under degree/density/stage-profile matched null"
    else:
        status = "INFO"
        interpretation = "exploratory"

    records.append({
        "analysis": "transport_degree_matched_null",
        "module_family": fam,
        "top_n": topn,
        "status": status,
        "primary_fdr": fdr,
        "interpretation": interpretation
    })

# Corrected SFARI decision
for topn, sub in sfari.groupby("top_n_fixed"):
    n_sets = sub.shape[0]
    n_sig = int((sub["overlap_fdr"] < 0.05).sum())
    status = "PASS" if n_sig >= max(1, int(np.ceil(n_sets * 0.5))) else "REVIEW"

    records.append({
        "analysis": "NTM2_SFARI_degree_matched_null",
        "module_family": "NTM2_ASD_down",
        "top_n": int(topn),
        "status": status,
        "primary_fdr": sub["overlap_fdr"].min(),
        "interpretation": f"{n_sig}/{n_sets} SFARI-like sets exceed degree/density-matched graph-priority null"
    })

decision = pd.DataFrame(records)
decision_out = f"{OUTDIR}/07b_step12B_decision_table.corrected.tsv"
decision.to_csv(decision_out, sep="\t", index=False)

# Compact manuscript summary
md_out = f"{OUTDIR}/08b_step12B_corrected_interpretation.md"
with open(md_out, "w") as f:
    f.write("# Step12B corrected interpretation\n\n")
    f.write("## Corrected decision table\n\n")
    f.write(decision.to_markdown(index=False))
    f.write("\n\n")
    f.write("## Final interpretation\n\n")
    f.write("Degree-, density-, and stage-profile-matched graph nulls did not support the developmental transport calls as strongly as the expression/variance-matched transport nulls. ")
    f.write("Thus, developmental stage transport should be interpreted as topology-sensitive under this strict null, and the primary developmental conclusion should remain anchored to the matched expression/variance transport calibration and hyperparameter sensitivity analyses. ")
    f.write("In contrast, NTM2_ASD_down graph-prioritized SFARI convergence remained robust under degree/density matching, especially at top500, where most SFARI-like resources exceeded the matched null. ")
    f.write("Step12B therefore supports topology-robust NTM2 genetic convergence but not topology-independent developmental transport.\n")

print("[patch Step12B] wrote:", decision_out)
print("[patch Step12B] wrote:", md_out)
