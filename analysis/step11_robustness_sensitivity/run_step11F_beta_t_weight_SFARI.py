#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

OUTDIR = "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step11_robustness_sensitivity/results"

GRAPH_NODES = "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step04_graph_construction/results/09_step04C_native_graph_nodes.tsv"
GRAPH_EDGES = "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step04_graph_construction/results/10_step04C_native_graph_edges.tsv.gz"
MODULE_WEIGHT_FILE = "/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv"

def clean_gene(x):
    x = str(x).strip()
    x = x.replace("gene:", "")
    x = x.replace("GENE:", "")
    x = x.split("|")[-1] if x.startswith("gene|") else x
    return x.strip()

def clean_risk(x):
    return str(x).strip().replace("risk_set:", "")

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    q = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if valid.sum() == 0:
        return q
    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(pv)
    qv = np.empty(m)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, ranked[i] * m / (i + 1))
        qv[order[i]] = prev
    q[valid] = np.minimum(qv, 1.0)
    return q

def hypergeom_enrichment(top_genes, risk_genes, universe):
    top_genes = set(map(clean_gene, top_genes)) & set(universe)
    risk_genes = set(map(clean_gene, risk_genes)) & set(universe)
    universe = set(universe)

    M = len(universe)
    n = len(risk_genes)
    N = len(top_genes)
    x = len(top_genes & risk_genes)

    if M == 0 or n == 0 or N == 0:
        return x, N, n, M, np.nan, np.nan

    p = hypergeom.sf(x - 1, M, n, N)
    a = x
    b = N - x
    c = n - x
    d = M - a - b - c
    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return x, N, n, M, odds, p

def find_col(df, candidates):
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None

print("[Step11F v3] Reading graph nodes")
nodes = pd.read_csv(GRAPH_NODES, sep="\t", low_memory=False)
node_col = find_col(nodes, ["node_id", "node", "id", "name", "node_name"])
type_col = find_col(nodes, ["node_type", "type", "class", "category"])

nodes[node_col] = nodes[node_col].astype(str)

if type_col is not None:
    gene_nodes_raw = nodes.loc[
        nodes[type_col].astype(str).str.lower().eq("gene"),
        node_col
    ].astype(str).tolist()
else:
    gene_nodes_raw = [
        x for x in nodes[node_col].astype(str).tolist()
        if not x.startswith("module:") and not x.startswith("stage:") and not x.startswith("risk_set:")
    ]

universe = sorted(set(clean_gene(x) for x in gene_nodes_raw if clean_gene(x) not in ["", "nan", "None"]))
print("[Step11F v3] Universe genes:", len(universe))

print("[Step11F v3] Reading graph edges")
edges = pd.read_csv(GRAPH_EDGES, sep="\t", compression="gzip", low_memory=False)

source_col = find_col(edges, ["source", "src", "from", "node1", "from_node"])
target_col = find_col(edges, ["target", "dst", "to", "node2", "to_node"])
etype_col = find_col(edges, ["edge_type", "type", "relation", "edge_class"])

edges[source_col] = edges[source_col].astype(str)
edges[target_col] = edges[target_col].astype(str)
edges[etype_col] = edges[etype_col].astype(str).str.replace("-", "_").str.replace(" ", "_")

risk_sets = {}
for _, r in edges.iterrows():
    et = str(r[etype_col]).lower()
    if "risk" not in et:
        continue

    s = str(r[source_col])
    t = str(r[target_col])

    if s.startswith("risk_set:"):
        risk_sets.setdefault(clean_risk(s), set()).add(clean_gene(t))
    elif t.startswith("risk_set:"):
        risk_sets.setdefault(clean_risk(t), set()).add(clean_gene(s))

sfari_sets = {k: sorted(v & set(universe)) for k, v in risk_sets.items() if "sfari" in k.lower()}

print("[Step11F v3] SFARI sets:")
for k, v in sfari_sets.items():
    print("  ", k, len(v))

print("[Step11F v3] Reading module weights")
mw = pd.read_csv(MODULE_WEIGHT_FILE, sep="\t", low_memory=False)

print("[Step11F v3] Module weight columns:")
print(list(mw.columns))

gene_col = find_col(mw, ["gene_symbol_fixed", "gene_symbol", "gene", "symbol"])
program_col = find_col(mw, ["program", "module_family", "module", "module_name"])
topn_col = find_col(mw, ["top_n", "topn", "n_top"])
weight_col = find_col(mw, ["weight", "module_weight", "signed_weight", "score"])
beta_col = find_col(mw, ["beta_ASD_vs_Control", "beta", "estimate", "effect", "logFC", "logfc"])
t_col = find_col(mw, ["t", "t_stat", "stat", "statistic"])

required = {
    "gene_col": gene_col,
    "program_col": program_col,
    "topn_col": topn_col,
    "weight_col": weight_col,
    "beta_col": beta_col,
    "t_col": t_col
}
print("[Step11F v3] Detected columns:")
for k, v in required.items():
    print("  ", k, "=", v)

for k in ["gene_col", "program_col", "topn_col", "weight_col"]:
    if required[k] is None:
        raise RuntimeError(f"Required column not detected: {k}")

mw2 = pd.DataFrame()
mw2["gene"] = mw[gene_col].map(clean_gene)
mw2["program"] = mw[program_col].astype(str)
mw2["top_n"] = pd.to_numeric(mw[topn_col], errors="coerce")
mw2["current_weight"] = pd.to_numeric(mw[weight_col], errors="coerce")

if beta_col is not None:
    mw2["beta_signed"] = pd.to_numeric(mw[beta_col], errors="coerce")
else:
    mw2["beta_signed"] = np.nan

if t_col is not None:
    mw2["t_signed"] = pd.to_numeric(mw[t_col], errors="coerce")
else:
    mw2["t_signed"] = np.nan

mw2 = mw2[
    mw2["program"].eq("NTM2_ASD_down") &
    mw2["top_n"].isin([200, 500]) &
    mw2["gene"].isin(universe)
].copy()

print("[Step11F v3] NTM2 module rows in graph universe:")
print(mw2.groupby(["program", "top_n"]).size())

# Add rank-normalized versions.
for base in ["current_weight", "beta_signed", "t_signed"]:
    rank_col = base + "_rank_normalized"
    mw2[rank_col] = np.nan
    for key, idx in mw2.groupby(["program", "top_n"]).groups.items():
        vals = pd.to_numeric(mw2.loc[idx, base], errors="coerce").fillna(0).values
        sign = np.sign(vals)
        sign[sign == 0] = 1
        ranks = pd.Series(np.abs(vals)).rank(method="average").values
        denom = np.nanmax(ranks)
        mw2.loc[idx, rank_col] = sign * ranks / denom if denom > 0 else sign

schemes = [
    "current_weight",
    "current_weight_rank_normalized",
    "beta_signed",
    "beta_signed_rank_normalized",
    "t_signed",
    "t_signed_rank_normalized",
]

records = []
diag_records = []

for (program, top_n), sub in mw2.groupby(["program", "top_n"]):
    for scheme in schemes:
        tmp = sub[["gene", scheme]].copy()
        tmp[scheme] = pd.to_numeric(tmp[scheme], errors="coerce")
        tmp = tmp[np.isfinite(tmp[scheme])].copy()

        n_nonzero = int((np.abs(tmp[scheme]) > 0).sum())
        diag_records.append({
            "program": program,
            "top_n": int(top_n),
            "weight_scheme": scheme,
            "n_genes_in_universe": tmp.shape[0],
            "n_nonzero_weights": n_nonzero,
            "min_weight": tmp[scheme].min() if tmp.shape[0] else np.nan,
            "max_weight": tmp[scheme].max() if tmp.shape[0] else np.nan,
        })

        if tmp.shape[0] == 0 or n_nonzero == 0:
            continue

        tmp = tmp.sort_values(scheme, key=lambda x: np.abs(x), ascending=False)

        for cutoff in [200, 500]:
            top_genes = tmp["gene"].head(min(cutoff, tmp.shape[0])).tolist()

            for risk_name, risk_genes in sfari_sets.items():
                x, N, n, M, odds, p = hypergeom_enrichment(top_genes, risk_genes, universe)
                records.append({
                    "module_family": program,
                    "top_n": int(top_n),
                    "weight_scheme": scheme,
                    "gene_cutoff": cutoff,
                    "risk_set": risk_name,
                    "overlap": x,
                    "tested_gene_n": N,
                    "risk_gene_n": n,
                    "universe_n": M,
                    "odds_ratio_approx": odds,
                    "p_value": p,
                })

diag = pd.DataFrame(diag_records)
diag_out = f"{OUTDIR}/12_step11F_v3_beta_t_weight_diagnostics.tsv"
diag.to_csv(diag_out, sep="\t", index=False)

res = pd.DataFrame(records)
if res.shape[0] > 0:
    res["fdr_within_scheme"] = np.nan
    for key, idx in res.groupby(["module_family", "top_n", "weight_scheme", "gene_cutoff"]).groups.items():
        res.loc[idx, "fdr_within_scheme"] = bh_fdr(res.loc[idx, "p_value"].values)

out1 = f"{OUTDIR}/13_step11F_v3_NTM2_beta_t_SFARI_enrichment.tsv"
res.to_csv(out1, sep="\t", index=False)

if res.shape[0] > 0:
    summary = (
        res.groupby(["module_family", "top_n", "weight_scheme", "gene_cutoff"], dropna=False)
        .agg(
            n_risk_sets=("risk_set", "count"),
            n_fdr_lt_0_05=("fdr_within_scheme", lambda x: int(np.sum(np.asarray(x, dtype=float) < 0.05))),
            median_overlap=("overlap", "median"),
            median_or=("odds_ratio_approx", "median"),
            min_p=("p_value", "min"),
            median_p=("p_value", "median"),
        )
        .reset_index()
    )
else:
    summary = pd.DataFrame()

out2 = f"{OUTDIR}/14_step11F_v3_NTM2_beta_t_SFARI_summary.tsv"
summary.to_csv(out2, sep="\t", index=False)

md = f"{OUTDIR}/15_step11F_v3_overall_summary.md"
with open(md, "w") as f:
    f.write("# Step11F v3 namespace-fixed beta/t-weight SFARI sensitivity summary\n\n")

    f.write("## Purpose\n\n")
    f.write("This analysis fixes gene-ID namespace mismatch by normalizing graph nodes, module genes, and risk-set genes to plain gene symbols before testing whether NTM2 SFARI convergence depends on the original signed weight definition.\n\n")

    f.write("## Diagnostics\n\n")
    f.write(f"- Graph universe genes after cleaning: {len(universe)}\n")
    f.write(f"- SFARI-like risk sets detected: {len(sfari_sets)}\n")
    for k, v in sfari_sets.items():
        f.write(f"  - {k}: {len(v)} genes in graph universe\n")
    f.write("\n")

    f.write("## Weight diagnostics\n\n")
    f.write(diag.to_markdown(index=False))
    f.write("\n\n")

    f.write("## Enrichment summary\n\n")
    if summary.shape[0] > 0:
        f.write(summary.to_markdown(index=False))
    else:
        f.write("No enrichment records were generated. Check gene overlap diagnostics.")
    f.write("\n\n")

    f.write("## Interpretation\n\n")
    f.write("If current, beta-based, t-statistic-based, and rank-normalized weights all retain SFARI enrichment for NTM2, this supports robustness of NTM2 genetic convergence to module-weight definition. If beta/t weights weaken while current/rank-compressed weights remain positive, Step11F should be reported as a partial sensitivity analysis focused on NTM2 genetic convergence rather than developmental transport.\n")

print("[Step11F v3] Wrote:", diag_out)
print("[Step11F v3] Wrote:", out1)
print("[Step11F v3] Wrote:", out2)
print("[Step11F v3] Wrote:", md)
