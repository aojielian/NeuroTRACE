#!/usr/bin/env python3

import math
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
STEP = BASE / "neurotrace_algorithm_project/step10_zhou2022_genesets_fixed"
OUT = STEP / "results"
GS = STEP / "gene_sets"

GRAPH_NODES = BASE / "neurotrace_algorithm_project/step04_graph_construction/results/09_step04C_native_graph_nodes.tsv"
NO_RISK_PRIORITY = BASE / "neurotrace_algorithm_project/step10_risk_free_ablation/results/93_step10A_no_risk_gene_graph_priority.tsv.gz"
MODULE_FILE = BASE / "neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv"

OUT.mkdir(parents=True, exist_ok=True)

CUTOFFS = [50, 100, 200, 500]
MAIN_TOP_N = [200, 500]

GENESETS_TO_USE = [
    "Zhou_Known_ASD_NDD_618",
    "Zhou_Known_ASD_NDD_high_pLI",
    "Zhou_Known_ASD_NDD_LOEUF_top",
    "Zhou_Known_ASD_NDD_GeneScore_1_or_1S",
    "Zhou_DeNovoWEST_stage1_p_lt_0.001",
    "Zhou_DeNovoWEST_stage1_FDR_lt_0.1",
    "Zhou_DeNovoWEST_stage1_exomewide",
    "Zhou_DeNovoWEST_stage1_LoF_or_Dmis_driven",
    "Zhou_TDT_prioritized_260",
    "Zhou_TDT_prioritized_autosomal",
    "Zhou_TDT_overtransmitted_HC_LoF",
    "Zhou_Selected_404_combined",
    "Zhou_Selected_159_deNovo",
    "Zhou_Selected_245_inherited_TDT",
    "Zhou_Combined_DeNovoWEST_FDR_lt_0.1",
    "Zhou_Combined_DeNovoWEST_exomewide",
    "Zhou_Meta_391_all_selected",
    "Zhou_Meta_inherited_plus_deNovo",
    "Zhou_Meta_novel",
    "Zhou_Meta_known",
    "Zhou_dnLoF_enriched_96_constrained",
    "Zhou_dnLoF_enriched_ASDdnSignif",
    "Zhou_dnLoF_enriched_SFARI_category",
    "Zhou_constrained_background_5754",
    "Zhou_ExACpLI",
    "Zhou_Archetype_A1_neurotransmission",
    "Zhou_Archetype_A2_chromatin_modification",
    "Zhou_Archetype_A3_RNA_processing",
    "Zhou_Archetype_A4_vesicle_mediated_transport",
    "Zhou_Archetype_A5_MAPK_signaling_migration",
    "Zhou_Archetype_A6_cytoskeleton_mitosis",
]

def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def clean_gene(x):
    if x is None or pd.isna(x):
        return ""
    y = str(x).strip().upper()
    if y in {"", "NA", "NAN", "NULL", "---", "."}:
        return ""
    return y

def read_genes(path):
    out = []
    with open(path, "r", errors="replace") as f:
        for line in f:
            g = clean_gene(line)
            if g:
                out.append(g)
    return sorted(set(out))

def geneset_category(name):
    if "Archetype" in name:
        return "mechanistic_archetype"
    if "TDT" in name or "inherited" in name:
        return "inherited_risk"
    if "DeNovo" in name or "dnLoF" in name:
        return "de_novo_risk"
    if "Meta" in name or "Selected" in name:
        return "combined_risk"
    if "Known_ASD_NDD" in name:
        return "known_ASD_NDD"
    if "constrained" in name or "ExAC" in name or "LOEUF" in name:
        return "constraint_background"
    return "other"

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

def enrichment(test_genes, set_genes, universe):
    universe = set(universe)
    test = set(test_genes).intersection(universe)
    gs = set(set_genes).intersection(universe)

    N = len(universe)
    n = len(test)
    K = len(gs)
    k = len(test.intersection(gs))

    pval = hypergeom_sf(k, K, n, N) if N and n and K else 1.0

    a = k
    b = n - k
    c = K - k
    d = N - K - b
    OR = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

    return {
        "overlap_n": k,
        "test_n": n,
        "geneset_n_in_universe": K,
        "universe_n": N,
        "odds_ratio": OR,
        "p_value": pval,
        "overlap_genes": ";".join(sorted(test.intersection(gs))),
    }

def main():
    print(f"[{now()}] Step10G4 Zhou2022 enrichment started")

    for f in [GRAPH_NODES, NO_RISK_PRIORITY, MODULE_FILE]:
        if not Path(f).exists():
            raise FileNotFoundError(f"Missing input: {f}")
    if not GS.exists():
        raise FileNotFoundError(f"Missing gene set directory: {GS}")

    nodes = pd.read_csv(GRAPH_NODES, sep="\t")
    universe = set(nodes.loc[nodes["node_type"] == "gene", "label"].map(clean_gene))
    universe = {g for g in universe if g}

    genesets = {}
    inv_rows = []
    for name in GENESETS_TO_USE:
        path = GS / f"{name}.txt"
        if not path.exists():
            continue
        genes = read_genes(path)
        if not genes:
            continue
        genesets[name] = genes
        inv_rows.append({
            "geneset": name,
            "category": geneset_category(name),
            "path": str(path),
            "n_raw_genes": len(genes),
            "n_in_graph_universe": len(set(genes).intersection(universe)),
        })

    inv = pd.DataFrame(inv_rows).sort_values(["category", "geneset"])
    inv.to_csv(OUT / "09_step10G4_Zhou_geneset_inventory.tsv", sep="\t", index=False)

    # no-risk graph-priority enrichment
    gp = pd.read_csv(NO_RISK_PRIORITY, sep="\t", compression="gzip")
    gp["gene"] = gp["gene"].map(clean_gene)
    gp["top_n"] = gp["top_n"].astype(str)

    rows = []
    for (program, top_n), sub in gp.groupby(["program", "top_n"]):
        sub = sub.sort_values("priority_rank")
        for cutoff in CUTOFFS:
            test_genes = sub.head(cutoff)["gene"].tolist()
            for gs_name, gs_genes in genesets.items():
                e = enrichment(test_genes, gs_genes, universe)
                rows.append({
                    "analysis": "no_risk_graph_priority",
                    "program": program,
                    "top_n": top_n,
                    "priority_cutoff": cutoff,
                    "geneset": gs_name,
                    "geneset_category": geneset_category(gs_name),
                    **e,
                })

    priority = pd.DataFrame(rows)
    priority["fdr"] = bh_fdr(priority["p_value"].values)
    priority = priority.sort_values(["fdr", "p_value", "program", "top_n", "priority_cutoff", "geneset"])
    priority.to_csv(OUT / "10_step10G4_no_risk_graph_priority_Zhou_enrichment.tsv", sep="\t", index=False)

    # native module enrichment
    mod = pd.read_csv(MODULE_FILE, sep="\t")
    mod["gene"] = mod["gene_symbol_fixed"].map(clean_gene)
    mod = mod[mod["top_n"].isin(MAIN_TOP_N)].copy()

    rows = []
    for (program, top_n), sub in mod.groupby(["program", "top_n"]):
        test_genes = sub["gene"].tolist()
        for gs_name, gs_genes in genesets.items():
            e = enrichment(test_genes, gs_genes, universe)
            rows.append({
                "analysis": "native_module_genes",
                "program": program,
                "top_n": top_n,
                "priority_cutoff": "",
                "geneset": gs_name,
                "geneset_category": geneset_category(gs_name),
                **e,
            })

    module = pd.DataFrame(rows)
    module["fdr"] = bh_fdr(module["p_value"].values)
    module = module.sort_values(["fdr", "p_value", "program", "top_n", "geneset"])
    module.to_csv(OUT / "11_step10G4_native_module_Zhou_enrichment.tsv", sep="\t", index=False)

    combined = pd.concat([priority, module], ignore_index=True)
    combined = combined.sort_values(["fdr", "p_value", "analysis", "program", "top_n", "priority_cutoff", "geneset"])
    combined.to_csv(OUT / "12_step10G4_combined_Zhou_enrichment.tsv", sep="\t", index=False)

    top_hits = combined.head(100)
    top_hits.to_csv(OUT / "13_step10G4_top_Zhou_enrichment_hits.tsv", sep="\t", index=False)

    cat_summary = combined.groupby(["analysis", "geneset_category"]).agg(
        n_tests=("p_value", "size"),
        n_fdr_lt_005=("fdr", lambda x: int((x < 0.05).sum())),
        min_fdr=("fdr", "min"),
        max_or=("odds_ratio", "max"),
        max_overlap=("overlap_n", "max")
    ).reset_index()
    cat_summary.to_csv(OUT / "14_step10G4_Zhou_category_summary.tsv", sep="\t", index=False)

    program_summary = combined.groupby(["analysis", "program", "geneset_category"]).agg(
        n_tests=("p_value", "size"),
        n_fdr_lt_005=("fdr", lambda x: int((x < 0.05).sum())),
        min_fdr=("fdr", "min"),
        max_or=("odds_ratio", "max"),
        max_overlap=("overlap_n", "max")
    ).reset_index()
    program_summary.to_csv(OUT / "15_step10G4_Zhou_program_category_summary.tsv", sep="\t", index=False)

    with open(OUT / "16_step10G4_Zhou_enrichment_summary.md", "w") as f:
        f.write("# NeuroTRACE Step10G4 Zhou 2022 enrichment summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Purpose\n")
        f.write("Test NeuroTRACE convergence with fixed Zhou 2022 gene sets, including known ASD/NDD genes, de novo risk sets, inherited/TDT genes, combined meta-analysis genes, constraint sets, and mechanistic archetypes.\n\n")
        f.write("## Gene-set inventory\n\n")
        f.write(inv.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Category summary\n\n")
        f.write(cat_summary.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Program-category summary\n\n")
        f.write(program_summary.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Top hits\n\n")
        f.write(top_hits.head(30).to_markdown(index=False))
        f.write("\n\n")
        f.write("## Interpretation\n")
        f.write("Use these results to assess convergence with independent ASD/NDD genetic-risk, inherited/de novo risk, constraint, and mechanistic archetype resources. Interpret results with graph-universe coverage.\n")

    print(f"[{now()}] Step10G4 done")
    print(f"[{now()}] Gene sets used: {len(genesets)}")
    print(f"[{now()}] Results: {OUT}")

if __name__ == "__main__":
    main()
