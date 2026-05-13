#!/usr/bin/env python3

import re
import math
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
STEP = BASE / "neurotrace_algorithm_project/step10_independent_genetic_risk"
OUT = STEP / "results"
OUT.mkdir(parents=True, exist_ok=True)

SATTERSTROM_FILE = BASE / "Satterstrom_102_ASD_risk_genes.txt"

GRAPH_NODES = BASE / "neurotrace_algorithm_project/step04_graph_construction/results/09_step04C_native_graph_nodes.tsv"
NO_RISK_PRIORITY = BASE / "neurotrace_algorithm_project/step10_risk_free_ablation/results/93_step10A_no_risk_gene_graph_priority.tsv.gz"
MODULE_FILE = BASE / "neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv"

CUTOFFS = [50, 100, 200, 500]
MAIN_TOP_N = [200, 500]


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def clean_gene(x):
    if x is None or pd.isna(x):
        return ""
    return str(x).strip().upper()


def parse_gene_file(path):
    genes = []
    with open(path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            tokens = re.split(r"[\s,;|]+", line)
            for t in tokens:
                g = clean_gene(t)
                if not g:
                    continue
                if g in {"GENE", "GENES", "SYMBOL", "GENE_SYMBOL", "HGNC", "NA", "NULL"}:
                    continue
                if re.match(r"^[A-Z0-9][A-Z0-9._-]{1,30}$", g):
                    genes.append(g)

    return sorted(set(genes))


def hypergeom_sf(k, K, n, N):
    if k <= 0:
        return 1.0
    upper = min(K, n)
    if k > upper:
        return 0.0

    def log_choose(a, b):
        if b < 0 or b > a:
            return -np.inf
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)

    logs = []
    denom = log_choose(N, n)
    for i in range(k, upper + 1):
        logs.append(log_choose(K, i) + log_choose(N - K, n - i) - denom)

    m = max(logs)
    return float(np.exp(m) * np.sum(np.exp(np.array(logs) - m)))


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    out = np.empty(n, dtype=float)
    prev = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        val = p[idx] * n / rank
        prev = min(prev, val)
        out[idx] = min(prev, 1.0)
    return out


def enrichment(test_genes, risk_genes, universe):
    universe = set(universe)
    test = set(test_genes).intersection(universe)
    risk = set(risk_genes).intersection(universe)

    N = len(universe)
    n = len(test)
    K = len(risk)
    k = len(test.intersection(risk))

    p = hypergeom_sf(k, K, n, N) if N > 0 and n > 0 and K > 0 else 1.0

    a = k
    b = n - k
    c = K - k
    d = N - K - b
    OR = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

    return {
        "overlap_n": k,
        "test_n": n,
        "risk_n_in_universe": K,
        "universe_n": N,
        "odds_ratio": OR,
        "p_value": p,
        "overlap_genes": ";".join(sorted(test.intersection(risk)))
    }


def main():
    print(f"[{now()}] Step10G2 Satterstrom102 independent genetic-risk enrichment started")

    for f in [SATTERSTROM_FILE, GRAPH_NODES, NO_RISK_PRIORITY, MODULE_FILE]:
        if not Path(f).exists():
            raise FileNotFoundError(f"Missing input: {f}")

    nodes = pd.read_csv(GRAPH_NODES, sep="\t")
    universe = set(nodes.loc[nodes["node_type"] == "gene", "label"].map(clean_gene))
    universe = {g for g in universe if g}

    satterstrom = parse_gene_file(SATTERSTROM_FILE)

    inventory = pd.DataFrame([{
        "geneset": "Satterstrom_102_ASD_FDR0.1",
        "category": "independent_genetic_risk",
        "path": str(SATTERSTROM_FILE),
        "n_raw_genes": len(satterstrom),
        "n_in_graph_universe": len(set(satterstrom).intersection(universe)),
        "n_missing_from_graph_universe": len(set(satterstrom) - universe),
        "missing_genes": ";".join(sorted(set(satterstrom) - universe))
    }])
    inventory.to_csv(OUT / "01_step10G2_Satterstrom102_inventory.tsv", sep="\t", index=False)

    # 1. No-risk graph-priority enrichment
    gp = pd.read_csv(NO_RISK_PRIORITY, sep="\t", compression="gzip")
    gp["gene"] = gp["gene"].map(clean_gene)
    gp["top_n"] = gp["top_n"].astype(str)

    rows = []
    for (program, top_n), sub in gp.groupby(["program", "top_n"]):
        sub = sub.sort_values("priority_rank")
        for cutoff in CUTOFFS:
            test_genes = sub.head(cutoff)["gene"].tolist()
            e = enrichment(test_genes, satterstrom, universe)
            rows.append({
                "analysis": "no_risk_graph_priority",
                "program": program,
                "top_n": top_n,
                "priority_cutoff": cutoff,
                "geneset": "Satterstrom_102_ASD_FDR0.1",
                "geneset_category": "independent_genetic_risk",
                **e
            })

    priority_enrich = pd.DataFrame(rows)
    priority_enrich["fdr"] = bh_fdr(priority_enrich["p_value"].values)
    priority_enrich = priority_enrich.sort_values(["fdr", "p_value", "program", "top_n", "priority_cutoff"])
    priority_enrich.to_csv(OUT / "02_step10G2_no_risk_graph_priority_Satterstrom102_enrichment.tsv", sep="\t", index=False)

    # 2. Native module gene enrichment
    mod = pd.read_csv(MODULE_FILE, sep="\t")
    mod["gene"] = mod["gene_symbol_fixed"].map(clean_gene)
    mod = mod[mod["top_n"].isin(MAIN_TOP_N)].copy()

    rows = []
    for (program, top_n), sub in mod.groupby(["program", "top_n"]):
        test_genes = sub["gene"].tolist()
        e = enrichment(test_genes, satterstrom, universe)
        rows.append({
            "analysis": "native_module_genes",
            "program": program,
            "top_n": top_n,
            "priority_cutoff": "",
            "geneset": "Satterstrom_102_ASD_FDR0.1",
            "geneset_category": "independent_genetic_risk",
            **e
        })

    module_enrich = pd.DataFrame(rows)
    module_enrich["fdr"] = bh_fdr(module_enrich["p_value"].values)
    module_enrich = module_enrich.sort_values(["fdr", "p_value", "program", "top_n"])
    module_enrich.to_csv(OUT / "03_step10G2_native_module_Satterstrom102_enrichment.tsv", sep="\t", index=False)

    # 3. Combine and summarize
    combined = pd.concat([priority_enrich, module_enrich], ignore_index=True)
    combined = combined.sort_values(["fdr", "p_value", "analysis", "program", "top_n", "priority_cutoff"])
    combined.to_csv(OUT / "04_step10G2_combined_Satterstrom102_enrichment.tsv", sep="\t", index=False)

    summary = combined.groupby(["analysis"]).agg(
        n_tests=("p_value", "size"),
        n_fdr_lt_005=("fdr", lambda x: int((x < 0.05).sum())),
        min_fdr=("fdr", "min"),
        max_or=("odds_ratio", "max"),
        max_overlap=("overlap_n", "max")
    ).reset_index()
    summary.to_csv(OUT / "05_step10G2_Satterstrom102_summary.tsv", sep="\t", index=False)

    top_hits = combined.head(30).copy()
    top_hits.to_csv(OUT / "06_step10G2_top_Satterstrom102_hits.tsv", sep="\t", index=False)

    with open(OUT / "07_step10G2_Satterstrom102_enrichment_summary.md", "w") as f:
        f.write("# NeuroTRACE Step10G2 Satterstrom 102 ASD risk gene enrichment summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Purpose\n")
        f.write("Test independent genetic-risk convergence using Satterstrom et al. 102 ASD risk genes at FDR <= 0.1. This is independent of SFARI and ArkingLab canonical transcriptomic modules.\n\n")
        f.write("## Gene-set inventory\n\n")
        f.write(inventory.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Summary\n\n")
        f.write(summary.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Top hits\n\n")
        f.write(top_hits.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Interpretation\n")
        if (combined["fdr"] < 0.05).any():
            f.write("Satterstrom independent ASD risk gene enrichment is detected in at least one NeuroTRACE output. This supports independent genetic-risk convergence beyond SFARI.\n")
        else:
            f.write("No FDR-significant Satterstrom enrichment was detected. Interpret SFARI enrichment as stronger than independent Satterstrom convergence, or consider larger/alternative genetic-risk resources.\n")

    print(f"[{now()}] Step10G2 done")
    print(f"[{now()}] Satterstrom raw genes: {len(satterstrom)}")
    print(f"[{now()}] Satterstrom genes in graph universe: {len(set(satterstrom).intersection(universe))}")
    print(f"[{now()}] Results: {OUT}")


if __name__ == "__main__":
    main()
