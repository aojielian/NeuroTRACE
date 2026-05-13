#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step12D1: gnomAD v2.1.1 held-out constraint validation for NeuroTRACE.

Purpose:
  Test whether NeuroTRACE native module genes and graph-prioritized genes
  converge on held-out gnomAD v2.1.1 constrained genes.

Inputs:
  - gnomAD constraint gene sets from Step12D0
  - native module weights
  - graph-priority gene tables

Outputs:
  03_step12D1_constraint_validation_enrichment.tsv
  04_step12D1_constraint_validation_summary.tsv
  05_step12D1_decision_table.tsv
  06_step12D1_overall_summary.md
"""

import argparse
import os
import re
import glob
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def clean_gene(x):
    x = str(x).strip()
    x = x.replace("gene:", "")
    x = x.replace("GENE:", "")
    return x


def find_col(df, candidates, required=False):
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    if required:
        raise ValueError(f"Missing column from {candidates}; available={list(df.columns)}")
    return None


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


def infer_module_family(x):
    x = str(x)
    if "NTM1" in x:
        return "NTM1_ASD_up"
    if "NTM2" in x:
        return "NTM2_ASD_down"
    if "NTM3" in x:
        return "NTM3_ASD_signed"
    return x


def infer_top_n(x):
    m = re.search(r"top[_|\\-]?(\\d+)", str(x), flags=re.I)
    if m:
        return int(m.group(1))
    return np.nan


def hypergeom_enrichment(top_genes, set_genes, universe):
    top_genes = set(clean_gene(x) for x in top_genes) & set(universe)
    set_genes = set(clean_gene(x) for x in set_genes) & set(universe)
    universe = set(universe)

    M = len(universe)
    n = len(set_genes)
    N = len(top_genes)
    x = len(top_genes & set_genes)

    if M == 0 or n == 0 or N == 0:
        return x, N, n, M, np.nan, np.nan

    p = hypergeom.sf(x - 1, M, n, N)

    a = x
    b = N - x
    c = n - x
    d = M - a - b - c

    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return x, N, n, M, odds, p


def load_constraint_sets(gene_set_file):
    gs = pd.read_csv(gene_set_file, sep="\t", low_memory=False)
    sets = {
        k: set(sub["gene"].astype(str).map(clean_gene))
        for k, sub in gs.groupby("gene_set")
    }
    return sets


def load_native_modules(module_weight_file):
    mw = pd.read_csv(module_weight_file, sep="\t", low_memory=False)

    gene_col = find_col(mw, ["gene_symbol_fixed", "gene_symbol", "gene", "symbol"], required=True)
    program_col = find_col(mw, ["program", "module_family", "module"], required=True)
    topn_col = find_col(mw, ["top_n", "topn"], required=True)
    weight_col = find_col(mw, ["weight", "module_weight", "signed_weight", "score"])

    mw["gene_clean"] = mw[gene_col].astype(str).map(clean_gene)
    mw["module_family"] = mw[program_col].astype(str).map(infer_module_family)
    mw["top_n_clean"] = pd.to_numeric(mw[topn_col], errors="coerce")

    if weight_col is not None:
        mw["weight_abs"] = pd.to_numeric(mw[weight_col], errors="coerce").abs()
    else:
        mw["weight_abs"] = 1.0

    records = []
    for (module, top_n), sub in mw.groupby(["module_family", "top_n_clean"]):
        if module not in ["NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed"]:
            continue
        if int(top_n) not in [50, 100, 200, 500]:
            continue
        genes = sub.sort_values("weight_abs", ascending=False)["gene_clean"].dropna().astype(str).tolist()
        records.append({
            "source_type": "native_module_genes",
            "source_file": module_weight_file,
            "module_family": module,
            "top_n": int(top_n),
            "priority_cutoff": int(top_n),
            "genes": genes
        })

    return records, set(mw["gene_clean"].dropna().astype(str))


def discover_priority_files(project_root):
    patterns = [
        "**/*gene_graph_priority*.tsv",
        "**/*gene_graph_priority*.tsv.gz",
        "**/*graph_priority_genes*.tsv",
        "**/*graph_priority_genes*.tsv.gz",
        "**/*priority_genes*.tsv",
        "**/*priority_genes*.tsv.gz",
    ]

    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(project_root, pat), recursive=True))

    files = sorted(set(files))

    usable = []
    for p in files:
        try:
            if p.endswith(".gz"):
                df = pd.read_csv(p, sep="\t", compression="gzip", nrows=20, low_memory=False)
            else:
                df = pd.read_csv(p, sep="\t", nrows=20, low_memory=False)
            cols = [c.lower() for c in df.columns]
            if any("gene" == c or c.endswith("gene") or "gene_node" == c for c in cols) and any("program" in c or "module" in c for c in cols):
                usable.append(p)
        except Exception:
            pass

    return usable


def load_priority_records(priority_file):
    if priority_file.endswith(".gz"):
        df = pd.read_csv(priority_file, sep="\t", compression="gzip", low_memory=False)
    else:
        df = pd.read_csv(priority_file, sep="\t", low_memory=False)

    gene_col = find_col(df, ["gene", "gene_symbol", "gene_node"], required=True)
    program_col = find_col(df, ["program", "module_family", "module", "module_node"], required=True)
    topn_col = find_col(df, ["top_n", "topn"])
    rank_col = find_col(df, ["priority_rank", "rank"])
    score_col = find_col(df, ["graph_priority_score", "priority_score", "score"])

    df["gene_clean"] = df[gene_col].astype(str).map(clean_gene)
    df["module_family"] = df[program_col].astype(str).map(infer_module_family)

    if topn_col is not None:
        df["top_n_clean"] = pd.to_numeric(df[topn_col], errors="coerce")
    else:
        df["top_n_clean"] = df[program_col].map(infer_top_n)

    if rank_col is not None:
        df["rank_clean"] = pd.to_numeric(df[rank_col], errors="coerce")
        df = df.sort_values(["module_family", "top_n_clean", "rank_clean"])
    elif score_col is not None:
        df["score_clean"] = pd.to_numeric(df[score_col], errors="coerce")
        df = df.sort_values(["module_family", "top_n_clean", "score_clean"], ascending=[True, True, False])
    else:
        df["rank_clean"] = df.groupby(["module_family", "top_n_clean"]).cumcount() + 1

    records = []
    for (module, top_n), sub in df.groupby(["module_family", "top_n_clean"], dropna=False):
        if module not in ["NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed"]:
            continue
        if not np.isfinite(top_n):
            continue
        top_n = int(top_n)

        all_genes = sub["gene_clean"].dropna().astype(str).tolist()

        for cutoff in [200, 500, 1000]:
            genes = all_genes[:min(cutoff, len(all_genes))]
            records.append({
                "source_type": "graph_priority",
                "source_file": priority_file,
                "module_family": module,
                "top_n": top_n,
                "priority_cutoff": cutoff,
                "genes": genes
            })

    return records, set(df["gene_clean"].dropna().astype(str))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene_sets", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step12_revision_strengthening/results/step12D_gnomad_constraint/01_step12D0_gnomad_constraint_gene_sets.tsv")
    parser.add_argument("--module_weights", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv")
    parser.add_argument("--project_root", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project")
    parser.add_argument("--outdir", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step12_revision_strengthening/results/step12D_gnomad_constraint")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    constraint_sets = load_constraint_sets(args.gene_sets)

    native_records, native_universe = load_native_modules(args.module_weights)

    priority_files = discover_priority_files(args.project_root)

    priority_records = []
    priority_universe = set()

    for p in priority_files:
        recs, uni = load_priority_records(p)
        priority_records.extend(recs)
        priority_universe |= uni

    priority_file_index = pd.DataFrame({"priority_file": priority_files})
    priority_file_index.to_csv(
        os.path.join(args.outdir, "02_step12D1_discovered_priority_files.tsv"),
        sep="\t",
        index=False
    )

    all_records = native_records + priority_records

    # Use graph-priority universe when available; otherwise native universe.
    base_universe = priority_universe if len(priority_universe) > 0 else native_universe

    # Intersect with gnomAD known genes.
    gnomad_genes = set()
    for s in constraint_sets.values():
        gnomad_genes |= s

    universe = set(base_universe) & set(gnomad_genes)

    results = []

    for rec in all_records:
        genes = rec["genes"]

        for set_name, set_genes in constraint_sets.items():
            x, N, n, M, odds, p = hypergeom_enrichment(genes, set_genes, universe)

            results.append({
                "source_type": rec["source_type"],
                "source_file": rec["source_file"],
                "module_family": rec["module_family"],
                "top_n": rec["top_n"],
                "priority_cutoff": rec["priority_cutoff"],
                "constraint_gene_set": set_name,
                "overlap": x,
                "tested_gene_n": N,
                "constraint_set_n": n,
                "universe_n": M,
                "odds_ratio": odds,
                "p_value": p
            })

    res = pd.DataFrame(results)

    # FDR within source_type/source_file/priority_cutoff across all constraint tests.
    res["fdr_global"] = bh_fdr(res["p_value"].values)

    res["fdr_within_source_cutoff"] = np.nan
    for key, idx in res.groupby(["source_type", "source_file", "priority_cutoff"]).groups.items():
        res.loc[idx, "fdr_within_source_cutoff"] = bh_fdr(res.loc[idx, "p_value"].values)

    out1 = os.path.join(args.outdir, "03_step12D1_constraint_validation_enrichment.tsv")
    res.to_csv(out1, sep="\t", index=False)

    summary = (
        res.groupby(["source_type", "module_family", "top_n", "priority_cutoff"], dropna=False)
        .agg(
            n_constraint_sets=("constraint_gene_set", "count"),
            n_fdr_lt_0_05=("fdr_within_source_cutoff", lambda x: int(np.sum(np.asarray(x, dtype=float) < 0.05))),
            median_overlap=("overlap", "median"),
            median_or=("odds_ratio", "median"),
            min_p=("p_value", "min"),
            min_fdr=("fdr_within_source_cutoff", "min")
        )
        .reset_index()
        .sort_values(["source_type", "module_family", "top_n", "priority_cutoff"])
    )

    out2 = os.path.join(args.outdir, "04_step12D1_constraint_validation_summary.tsv")
    summary.to_csv(out2, sep="\t", index=False)

    decision_records = []
    focus = res[
        res["constraint_gene_set"].isin([
            "gnomAD_v2.1.1_LOEUF_lt_0.35",
            "gnomAD_v2.1.1_pLI_gt_0.90",
            "gnomAD_v2.1.1_misZ_top_10pct",
            "gnomAD_v2.1.1_LOEUF_lowest_10pct"
        ])
    ].copy()

    for (source_type, module, top_n, cutoff), sub in focus.groupby(["source_type", "module_family", "top_n", "priority_cutoff"], dropna=False):
        n_sets = sub.shape[0]
        n_sig = int((sub["fdr_within_source_cutoff"] < 0.05).sum())
        min_fdr = sub["fdr_within_source_cutoff"].min()
        median_or = sub["odds_ratio"].median()

        if module == "NTM2_ASD_down" and int(top_n) in [200, 500]:
            status = "PASS" if n_sig >= 1 and median_or > 1 else "REVIEW"
        else:
            status = "INFO" if n_sig == 0 else "SUPPORTIVE"

        decision_records.append({
            "source_type": source_type,
            "module_family": module,
            "top_n": top_n,
            "priority_cutoff": cutoff,
            "status": status,
            "n_primary_constraint_sets": n_sets,
            "n_primary_sets_fdr_lt_0_05": n_sig,
            "min_fdr": min_fdr,
            "median_or": median_or
        })

    decision = pd.DataFrame(decision_records)
    out3 = os.path.join(args.outdir, "05_step12D1_decision_table.tsv")
    decision.to_csv(out3, sep="\t", index=False)

    md = os.path.join(args.outdir, "06_step12D1_overall_summary.md")
    with open(md, "w") as f:
        f.write("# Step12D1 gnomAD v2.1.1 held-out constraint validation\n\n")

        f.write("## Purpose\n\n")
        f.write("This analysis evaluates whether NeuroTRACE native module genes and graph-prioritized genes converge on held-out gnomAD v2.1.1 gene constraint resources. These resources were not added as graph nodes and were used only for post hoc validation.\n\n")

        f.write("## Inputs\n\n")
        f.write(f"- Constraint gene sets: `{args.gene_sets}`\n")
        f.write(f"- Module weights: `{args.module_weights}`\n")
        f.write(f"- Discovered priority files: {len(priority_files)}\n")
        f.write(f"- Enrichment universe size: {len(universe)}\n\n")

        f.write("## Discovered priority files\n\n")
        f.write(priority_file_index.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Constraint validation summary\n\n")
        f.write(summary.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Decision table\n\n")
        f.write(decision.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Interpretation guide\n\n")
        f.write("- PASS for NTM2 indicates that at least one primary gnomAD constraint set is significantly enriched after FDR correction and median OR is greater than 1.\n")
        f.write("- This analysis should be described as held-out constraint validation because gnomAD constraint resources are not used as graph nodes or priors.\n")
        f.write("- If NTM2 top500 is stronger than top200, report top500 as the robust constraint-convergence threshold.\n")

    print("[Step12D1] Wrote:", out1)
    print("[Step12D1] Wrote:", out2)
    print("[Step12D1] Wrote:", out3)
    print("[Step12D1] Wrote:", md)


if __name__ == "__main__":
    main()
