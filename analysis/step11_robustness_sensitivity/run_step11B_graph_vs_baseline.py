#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step11B: Compare NeuroTRACE full graph diffusion with simpler real-data baselines.

Purpose
-------
This step tests whether the graph-diffusion layer adds value beyond simpler
developmental mapping baselines.

Compared methods
----------------
1. direct_module_stage
   Uses only direct module-stage edges already encoded in the native graph.

2. gene_stage_projection
   Uses a simple module-gene by stage-gene projection, without PPR diffusion.

3. graph_diffusion_only
   Uses PPR probability mass on developmental stage nodes.

4. full_graph
   Combines direct module-stage alignment and graph diffusion.

5. full_graph_no_risk
   Same as full_graph but removes risk-gene prior edges before diffusion.

Additional gene-priority comparison
-----------------------------------
For NTM2_ASD_down, compare SFARI enrichment among:
  - native module genes
  - full graph-prioritized genes
  - no-risk graph-prioritized genes

Outputs
-------
01_step11B_stage_method_comparison.tsv
02_step11B_main_call_summary.tsv
03_step11B_method_rank_summary.tsv
04_step11B_gene_priority_SFARI_comparison.tsv
05_step11B_overall_summary.md
06_step11B_manuscript_ready_text.md
"""

import argparse
import os
import re
import math
import numpy as np
import pandas as pd
from collections import defaultdict

try:
    from scipy import sparse
    from scipy.stats import hypergeom
except Exception as e:
    raise RuntimeError("This script requires scipy, numpy, and pandas.") from e


def log(msg):
    print(f"[Step11B] {msg}", flush=True)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def read_table_auto(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing input file: {path}")
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def find_col(df, candidates, required=True, label="column"):
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    if required:
        raise ValueError(
            f"Could not find {label}. Tried {candidates}. Available columns: {list(df.columns)}"
        )
    return None


def normalize_edge_type(x):
    if pd.isna(x):
        return "unknown"
    x = str(x).strip().replace("-", "_").replace(" ", "_")
    return x


def clean_stage_name(x):
    return str(x).replace("stage:", "")


def zscore(x):
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return x
    sd = np.nanstd(x)
    if (not np.isfinite(sd)) or sd == 0:
        return np.zeros_like(x)
    return (x - np.nanmean(x)) / sd


def softmax(x, temperature=0.75):
    x = np.asarray(x, dtype=float)
    if len(x) == 0:
        return x
    x = x / temperature
    x = x - np.nanmax(x)
    ex = np.exp(x)
    denom = np.nansum(ex)
    if denom <= 0 or not np.isfinite(denom):
        return np.ones_like(x) / len(x)
    return ex / denom


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
    qv = np.empty(m, dtype=float)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        val = ranked[i] * m / (i + 1)
        prev = min(prev, val)
        qv[order[i]] = prev
    q[valid] = np.minimum(qv, 1.0)
    return q


def infer_module_family(module_name):
    x = str(module_name).lower()
    if "ntm1" in x:
        return "NTM1_ASD_up"
    if "ntm2" in x:
        return "NTM2_ASD_down"
    if "ntm3" in x:
        return "NTM3_ASD_signed"
    return str(module_name)


def infer_top_n(module_name):
    x = str(module_name)
    m = re.search(r"top[_\|\-]?(\d+)", x, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))
    m = re.search(r"(\d+)$", x)
    if m:
        return int(m.group(1))
    return np.nan


def guess_node_classes(nodes, node_col, type_col=None):
    classes = {}
    stage_tokens = [
        "early_prenatal", "mid_prenatal", "late_prenatal",
        "infancy", "infant", "child", "adolescence", "adolescent",
        "adulthood", "adult", "prenatal", "postnatal"
    ]

    for _, row in nodes.iterrows():
        node = str(row[node_col])
        low = node.lower()
        t = ""
        if type_col is not None and type_col in nodes.columns:
            t = str(row[type_col]).lower()

        if "module" in t or low.startswith("module:") or low.startswith("ntm"):
            cls = "module"
        elif "stage" in t or low.startswith("stage:") or any(s in low for s in stage_tokens):
            cls = "stage"
        elif "risk" in t or low.startswith("risk_set:") or "sfari" in low:
            cls = "risk"
        elif "gene" in t:
            cls = "gene"
        else:
            cls = "gene"

        classes[node] = cls

    return classes


def detect_signed_weight_col(edges, weight_col):
    """
    Try to identify a signed module-gene weight column.
    If none exists, use edge weight and record that signed information is unavailable.
    """
    candidates = [
        "signed_weight", "module_weight", "raw_weight", "weight_signed",
        "signed_score", "gene_weight", "ntm_weight", "score_signed", "beta_weight"
    ]
    lower = {c.lower(): c for c in edges.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()], True
    return weight_col, False


def edge_weight_value(row, col):
    try:
        x = float(row[col])
        if np.isfinite(x):
            return x
        return 0.0
    except Exception:
        return 0.0


def extract_direct_module_stage(edges, source_col, target_col, weight_col, edge_type_col, module_nodes, stage_nodes):
    module_set = set(module_nodes)
    stage_set = set(stage_nodes)
    direct = defaultdict(float)

    for _, r in edges.iterrows():
        et = normalize_edge_type(r[edge_type_col])
        if "module_stage" not in et and not ("module" in et and "stage" in et):
            continue

        s = str(r[source_col])
        t = str(r[target_col])
        w = edge_weight_value(r, weight_col)

        if s in module_set and t in stage_set:
            direct[(s, t)] += w
        elif t in module_set and s in stage_set:
            direct[(t, s)] += w

    return direct


def extract_module_gene(edges, source_col, target_col, weight_col, signed_weight_col, edge_type_col, module_nodes, gene_nodes):
    module_set = set(module_nodes)
    gene_set = set(gene_nodes)
    mg = defaultdict(dict)

    for _, r in edges.iterrows():
        et = normalize_edge_type(r[edge_type_col])
        if "module_gene" not in et and not ("module" in et and "gene" in et):
            continue

        s = str(r[source_col])
        t = str(r[target_col])

        if s in module_set and t in gene_set:
            m, g = s, t
        elif t in module_set and s in gene_set:
            m, g = t, s
        else:
            continue

        w = edge_weight_value(r, signed_weight_col)
        mg[m][g] = mg[m].get(g, 0.0) + w

    return mg


def extract_stage_gene(edges, source_col, target_col, weight_col, edge_type_col, stage_nodes, gene_nodes):
    stage_set = set(stage_nodes)
    gene_set = set(gene_nodes)
    sg = defaultdict(dict)

    for _, r in edges.iterrows():
        et = normalize_edge_type(r[edge_type_col])
        if "stage_gene" not in et and not ("stage" in et and "gene" in et):
            continue

        s = str(r[source_col])
        t = str(r[target_col])

        if s in stage_set and t in gene_set:
            st, g = s, t
        elif t in stage_set and s in gene_set:
            st, g = t, s
        else:
            continue

        w = edge_weight_value(r, weight_col)
        sg[st][g] = sg[st].get(g, 0.0) + w

    return sg


def extract_risk_sets_from_edges(edges, source_col, target_col, edge_type_col, node_classes):
    risk_sets = defaultdict(set)
    for _, r in edges.iterrows():
        et = normalize_edge_type(r[edge_type_col])
        if "risk" not in et:
            continue
        s = str(r[source_col])
        t = str(r[target_col])
        cs = node_classes.get(s, "unknown")
        ct = node_classes.get(t, "unknown")
        if cs == "risk" and ct == "gene":
            risk_sets[s].add(t)
        elif ct == "risk" and cs == "gene":
            risk_sets[t].add(s)
    return dict(risk_sets)


def build_transition(edges, all_nodes, source_col, target_col, weight_col, edge_type_col, edge_scales):
    node_to_idx = {n: i for i, n in enumerate(all_nodes)}
    rows, cols, vals = [], [], []

    for _, r in edges.iterrows():
        s = str(r[source_col])
        t = str(r[target_col])
        if s not in node_to_idx or t not in node_to_idx:
            continue

        et = normalize_edge_type(r[edge_type_col])
        scale = edge_scales.get(et, edge_scales.get("default", 1.0))
        if scale <= 0:
            continue

        w = abs(edge_weight_value(r, weight_col)) * scale
        if w <= 0:
            continue

        i = node_to_idx[s]
        j = node_to_idx[t]
        rows.extend([i, j])
        cols.extend([j, i])
        vals.extend([w, w])

    n = len(all_nodes)
    A = sparse.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    row_sum = np.asarray(A.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    P = sparse.diags(1.0 / row_sum) @ A
    return P


def run_ppr(P, seed_idx, alpha=0.35, max_iter=120, tol=1e-10):
    n = P.shape[0]
    e = np.zeros(n)
    e[seed_idx] = 1.0
    p = e.copy()
    for _ in range(max_iter):
        new_p = alpha * e + (1 - alpha) * (P.T @ p)
        if np.linalg.norm(new_p - p, ord=1) < tol:
            p = new_p
            break
        p = new_p
    s = p.sum()
    if s > 0:
        p = p / s
    return p


def normalized_dot(a_dict, b_dict):
    genes = sorted(set(a_dict.keys()) & set(b_dict.keys()))
    if len(genes) == 0:
        return 0.0, 0

    a = np.array([a_dict[g] for g in genes], dtype=float)
    b = np.array([b_dict[g] for g in genes], dtype=float)

    denom = math.sqrt(float(np.sum(a * a))) * math.sqrt(float(np.sum(b * b)))
    if denom == 0 or not np.isfinite(denom):
        return 0.0, len(genes)

    return float(np.sum(a * b) / denom), len(genes)


def entropy(p):
    p = np.asarray(p, dtype=float)
    p = p[p > 0]
    if len(p) == 0:
        return np.nan
    return float(-np.sum(p * np.log(p)))


def summarize_stage_scores(raw_scores, stage_nodes, method, module, expected_stage, temperature):
    raw = np.asarray(raw_scores, dtype=float)
    z = zscore(raw)
    prob = softmax(z, temperature=temperature)

    order = np.argsort(prob)[::-1]
    top_stage = stage_nodes[order[0]]
    top_prob = float(prob[order[0]])
    second_prob = float(prob[order[1]]) if len(order) > 1 else np.nan
    margin = top_prob - second_prob if np.isfinite(second_prob) else np.nan

    expected_prob = np.nan
    expected_rank = np.nan
    for rank, idx in enumerate(order, start=1):
        if clean_stage_name(stage_nodes[idx]).lower() == expected_stage.lower():
            expected_prob = float(prob[idx])
            expected_rank = rank
            break

    records = []
    for st, raw_val, z_val, p_val in zip(stage_nodes, raw, z, prob):
        records.append({
            "method": method,
            "module_node": module,
            "module_family": infer_module_family(module),
            "top_n": infer_top_n(module),
            "stage": st,
            "stage_clean": clean_stage_name(st),
            "raw_score": raw_val,
            "z_score": z_val,
            "stage_probability": p_val,
            "top_stage": top_stage,
            "top_stage_clean": clean_stage_name(top_stage),
            "top_probability": top_prob,
            "second_probability": second_prob,
            "top_margin": margin,
            "entropy": entropy(prob),
            "expected_stage": expected_stage,
            "expected_stage_probability": expected_prob,
            "expected_stage_rank": expected_rank,
            "expected_stage_is_top": int(clean_stage_name(top_stage).lower() == expected_stage.lower())
        })
    return records


def hypergeom_enrichment(top_genes, risk_genes, universe_genes):
    top_genes = set(top_genes) & set(universe_genes)
    risk_genes = set(risk_genes) & set(universe_genes)
    universe_genes = set(universe_genes)

    M = len(universe_genes)
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--nodes",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step04_graph_construction/results/09_step04C_native_graph_nodes.tsv"
    )
    parser.add_argument(
        "--edges",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step04_graph_construction/results/10_step04C_native_graph_edges.tsv.gz"
    )
    parser.add_argument(
        "--outdir",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step11_robustness_sensitivity/results"
    )
    parser.add_argument("--expected_stage", default="late_prenatal")
    parser.add_argument("--temperature", type=float, default=0.75)
    parser.add_argument("--alpha", type=float, default=0.35)
    parser.add_argument("--priority_top_gene_k", type=int, default=500)
    parser.add_argument("--main_module_regex", default=r"NTM[123].*(top[_\|\-]?(200|500)|200|500)")
    args = parser.parse_args()

    ensure_dir(args.outdir)

    log(f"Reading nodes: {args.nodes}")
    nodes = read_table_auto(args.nodes)
    log(f"Reading edges: {args.edges}")
    edges = read_table_auto(args.edges)

    node_col = find_col(nodes, ["node_id", "node", "id", "name", "node_name"], label="node id")
    type_col = find_col(nodes, ["node_type", "type", "class", "category"], required=False, label="node type")

    source_col = find_col(edges, ["source", "src", "from", "node1", "from_node"], label="edge source")
    target_col = find_col(edges, ["target", "dst", "to", "node2", "to_node"], label="edge target")
    weight_col = find_col(edges, ["weight", "edge_weight", "w", "score"], required=False, label="edge weight")
    edge_type_col = find_col(edges, ["edge_type", "type", "relation", "edge_class"], required=False, label="edge type")

    if weight_col is None:
        weight_col = "__unit_weight__"
        edges[weight_col] = 1.0
    if edge_type_col is None:
        edge_type_col = "__edge_type__"
        edges[edge_type_col] = "unknown"

    signed_weight_col, signed_available = detect_signed_weight_col(edges, weight_col)

    nodes[node_col] = nodes[node_col].astype(str)
    edges[source_col] = edges[source_col].astype(str)
    edges[target_col] = edges[target_col].astype(str)
    edges[edge_type_col] = edges[edge_type_col].map(normalize_edge_type)

    node_classes = guess_node_classes(nodes, node_col, type_col)
    all_nodes = list(nodes[node_col].astype(str))
    node_to_idx = {n: i for i, n in enumerate(all_nodes)}

    module_nodes = [n for n in all_nodes if node_classes.get(n) == "module"]
    stage_nodes = [n for n in all_nodes if node_classes.get(n) == "stage"]
    gene_nodes = [n for n in all_nodes if node_classes.get(n) == "gene"]

    main_re = re.compile(args.main_module_regex, flags=re.IGNORECASE)
    main_modules = [m for m in module_nodes if main_re.search(m)]
    if len(main_modules) == 0:
        main_modules = module_nodes

    log(f"Detected {len(gene_nodes)} gene nodes, {len(module_nodes)} module nodes, {len(stage_nodes)} stage nodes.")
    log(f"Main modules: {len(main_modules)}")
    log(f"Stage nodes: {stage_nodes}")
    log(f"Module-gene signed weight column used: {signed_weight_col}; signed_available={signed_available}")

    direct_ms = extract_direct_module_stage(
        edges, source_col, target_col, weight_col, edge_type_col, main_modules, stage_nodes
    )
    module_gene = extract_module_gene(
        edges, source_col, target_col, weight_col, signed_weight_col, edge_type_col, main_modules, gene_nodes
    )
    stage_gene = extract_stage_gene(
        edges, source_col, target_col, weight_col, edge_type_col, stage_nodes, gene_nodes
    )
    risk_sets = extract_risk_sets_from_edges(edges, source_col, target_col, edge_type_col, node_classes)
    sfari_sets = {k: v for k, v in risk_sets.items() if "sfari" in k.lower()}

    log(f"Direct module-stage pairs: {len(direct_ms)}")
    log(f"Module-gene modules with edges: {len(module_gene)}")
    log(f"Stage-gene stages with edges: {len(stage_gene)}")
    log(f"SFARI-like risk sets: {len(sfari_sets)}")

    default_scales = {
        "gene_gene": 0.35,
        "gene_knn": 0.35,
        "knn": 0.35,
        "module_gene": 1.00,
        "stage_gene": 0.80,
        "module_stage": 1.25,
        "risk_gene_prior": 0.30,
        "risk_gene": 0.30,
        "risk": 0.30,
        "default": 1.0
    }

    no_risk_scales = dict(default_scales)
    no_risk_scales["risk_gene_prior"] = 0.0
    no_risk_scales["risk_gene"] = 0.0
    no_risk_scales["risk"] = 0.0

    log("Building full graph transition matrix.")
    P_full = build_transition(edges, all_nodes, source_col, target_col, weight_col, edge_type_col, default_scales)
    log("Building no-risk graph transition matrix.")
    P_no_risk = build_transition(edges, all_nodes, source_col, target_col, weight_col, edge_type_col, no_risk_scales)

    stage_indices = [node_to_idx[s] for s in stage_nodes]
    gene_indices = [node_to_idx[g] for g in gene_nodes]

    stage_records = []
    gene_priority_records = []

    for module in main_modules:
        if module not in node_to_idx:
            continue

        # 1. Direct module-stage baseline.
        direct_scores = np.array([direct_ms.get((module, s), 0.0) for s in stage_nodes], dtype=float)
        stage_records.extend(
            summarize_stage_scores(
                direct_scores, stage_nodes, "direct_module_stage",
                module, args.expected_stage, args.temperature
            )
        )

        # 2. Simple gene-stage projection baseline.
        projection_scores = []
        projection_overlap = []
        for st in stage_nodes:
            score, n_overlap = normalized_dot(module_gene.get(module, {}), stage_gene.get(st, {}))
            projection_scores.append(score)
            projection_overlap.append(n_overlap)

        projection_records = summarize_stage_scores(
            projection_scores, stage_nodes, "gene_stage_projection",
            module, args.expected_stage, args.temperature
        )
        for rec, n_ov in zip(projection_records, projection_overlap):
            rec["gene_stage_overlap_n"] = n_ov
            rec["module_gene_signed_available"] = signed_available
            rec["module_gene_weight_column"] = signed_weight_col
        stage_records.extend(projection_records)

        # 3. Full graph diffusion only.
        ppr_full = run_ppr(P_full, node_to_idx[module], alpha=args.alpha)
        graph_scores = np.array([ppr_full[i] for i in stage_indices], dtype=float)
        stage_records.extend(
            summarize_stage_scores(
                graph_scores, stage_nodes, "graph_diffusion_only",
                module, args.expected_stage, args.temperature
            )
        )

        # 4. Full graph = direct z + graph z.
        full_scores = 0.50 * zscore(direct_scores) + 0.50 * zscore(graph_scores)
        stage_records.extend(
            summarize_stage_scores(
                full_scores, stage_nodes, "full_graph",
                module, args.expected_stage, args.temperature
            )
        )

        # 5. Full graph no-risk.
        ppr_no_risk = run_ppr(P_no_risk, node_to_idx[module], alpha=args.alpha)
        no_risk_graph_scores = np.array([ppr_no_risk[i] for i in stage_indices], dtype=float)
        full_no_risk_scores = 0.50 * zscore(direct_scores) + 0.50 * zscore(no_risk_graph_scores)
        stage_records.extend(
            summarize_stage_scores(
                full_no_risk_scores, stage_nodes, "full_graph_no_risk",
                module, args.expected_stage, args.temperature
            )
        )

        # Gene-priority enrichment comparison for SFARI.
        if len(sfari_sets) > 0:
            native_genes_sorted = sorted(
                module_gene.get(module, {}).items(),
                key=lambda kv: abs(kv[1]),
                reverse=True
            )
            native_genes = [g for g, w in native_genes_sorted]

            full_gene_scores = pd.DataFrame({
                "gene": gene_nodes,
                "score": [ppr_full[i] for i in gene_indices]
            }).sort_values("score", ascending=False)
            full_genes = full_gene_scores["gene"].head(args.priority_top_gene_k).tolist()

            no_risk_gene_scores = pd.DataFrame({
                "gene": gene_nodes,
                "score": [ppr_no_risk[i] for i in gene_indices]
            }).sort_values("score", ascending=False)
            no_risk_genes = no_risk_gene_scores["gene"].head(args.priority_top_gene_k).tolist()

            gene_sets_to_test = {
                "native_module_genes": native_genes,
                "full_graph_priority": full_genes,
                "no_risk_graph_priority": no_risk_genes
            }

            for source_name, genes in gene_sets_to_test.items():
                if source_name == "native_module_genes":
                    genes = genes[:min(len(genes), args.priority_top_gene_k)]

                for risk_name, risk_genes in sfari_sets.items():
                    x, N, n, M, odds, pval = hypergeom_enrichment(genes, risk_genes, gene_nodes)
                    gene_priority_records.append({
                        "module_node": module,
                        "module_family": infer_module_family(module),
                        "module_top_n": infer_top_n(module),
                        "gene_set_source": source_name,
                        "priority_top_gene_k": args.priority_top_gene_k,
                        "risk_set": risk_name,
                        "overlap": x,
                        "tested_gene_n": N,
                        "risk_gene_n": n,
                        "universe_n": M,
                        "odds_ratio_approx": odds,
                        "p_value": pval
                    })

    stage_df = pd.DataFrame(stage_records)

    # Add main-call interpretation.
    stage_df["is_main_stage_record"] = stage_df["stage_clean"].str.lower() == args.expected_stage.lower()
    stage_df["main_call_expected"] = "none"
    stage_df.loc[
        stage_df["module_family"].isin(["NTM1_ASD_up", "NTM3_ASD_signed"]),
        "main_call_expected"
    ] = "late_prenatal"
    stage_df.loc[
        stage_df["module_family"].eq("NTM2_ASD_down"),
        "main_call_expected"
    ] = "not_late_prenatal"

    out_stage = os.path.join(args.outdir, "01_step11B_stage_method_comparison.tsv")
    stage_df.to_csv(out_stage, sep="\t", index=False)

    # Main call summary by method/module.
    top_only = stage_df[stage_df["stage_clean"] == stage_df["top_stage_clean"]].copy()

    main_summary = (
        top_only
        .groupby(["method", "module_node", "module_family", "top_n"], dropna=False)
        .agg(
            top_stage=("top_stage_clean", "first"),
            top_probability=("top_probability", "first"),
            second_probability=("second_probability", "first"),
            top_margin=("top_margin", "first"),
            entropy=("entropy", "first"),
            expected_stage_probability=("expected_stage_probability", "first"),
            expected_stage_rank=("expected_stage_rank", "first"),
            expected_stage_is_top=("expected_stage_is_top", "first")
        )
        .reset_index()
    )

    def pass_status(row):
        fam = row["module_family"]
        topn = row["top_n"]
        if fam in ["NTM1_ASD_up", "NTM3_ASD_signed"] and topn in [200, 500]:
            return "PASS" if row["top_stage"] == args.expected_stage else "REVIEW"
        if fam == "NTM2_ASD_down" and topn in [200, 500]:
            return "PASS" if row["top_stage"] != args.expected_stage else "REVIEW"
        return "INFO"

    main_summary["main_call_status"] = main_summary.apply(pass_status, axis=1)

    out_main = os.path.join(args.outdir, "02_step11B_main_call_summary.tsv")
    main_summary.to_csv(out_main, sep="\t", index=False)

    # Method rank/aggregate summary.
    method_summary = (
        main_summary
        .groupby("method", dropna=False)
        .agg(
            n_main_tests=("main_call_status", "count"),
            pass_fraction=("main_call_status", lambda x: np.mean(np.array(x) == "PASS")),
            median_top_margin=("top_margin", "median"),
            mean_top_margin=("top_margin", "mean"),
            median_entropy=("entropy", "median"),
            n_late_top=("top_stage", lambda x: np.sum(np.array(x) == args.expected_stage)),
            n_mid_top=("top_stage", lambda x: np.sum(np.array(x) == "mid_prenatal")),
            n_adolescence_top=("top_stage", lambda x: np.sum(np.array(x) == "adolescence")),
            n_adulthood_top=("top_stage", lambda x: np.sum(np.array(x) == "adulthood"))
        )
        .reset_index()
    )

    out_method = os.path.join(args.outdir, "03_step11B_method_rank_summary.tsv")
    method_summary.to_csv(out_method, sep="\t", index=False)

    # Gene-priority SFARI comparison.
    if len(gene_priority_records) > 0:
        gp_df = pd.DataFrame(gene_priority_records)
        gp_df["fdr_within_source"] = np.nan
        for source_name, idx in gp_df.groupby("gene_set_source").groups.items():
            gp_df.loc[idx, "fdr_within_source"] = bh_fdr(gp_df.loc[idx, "p_value"].values)

        gp_df["fdr_global"] = bh_fdr(gp_df["p_value"].values)

        out_gp = os.path.join(args.outdir, "04_step11B_gene_priority_SFARI_comparison.tsv")
        gp_df.to_csv(out_gp, sep="\t", index=False)
    else:
        gp_df = pd.DataFrame()
        out_gp = os.path.join(args.outdir, "04_step11B_gene_priority_SFARI_comparison.tsv")
        gp_df.to_csv(out_gp, sep="\t", index=False)

    # Summaries for manuscript.
    summary_path = os.path.join(args.outdir, "05_step11B_overall_summary.md")
    text_path = os.path.join(args.outdir, "06_step11B_manuscript_ready_text.md")

    # Compare full graph vs direct method margins.
    margin_compare = main_summary.pivot_table(
        index=["module_node", "module_family", "top_n"],
        columns="method",
        values="top_margin",
        aggfunc="first"
    ).reset_index()
    if "full_graph" in margin_compare.columns and "direct_module_stage" in margin_compare.columns:
        margin_compare["full_minus_direct_margin"] = (
            margin_compare["full_graph"] - margin_compare["direct_module_stage"]
        )
    margin_compare_path = os.path.join(args.outdir, "07_step11B_full_vs_direct_margin.tsv")
    margin_compare.to_csv(margin_compare_path, sep="\t", index=False)

    # Gene-priority compact summary.
    if gp_df.shape[0] > 0:
        gp_summary = (
            gp_df
            .groupby(["module_node", "module_family", "module_top_n", "gene_set_source"], dropna=False)
            .agg(
                n_risk_sets=("risk_set", "count"),
                n_fdr_global_lt_0_05=("fdr_global", lambda x: np.sum(np.asarray(x, dtype=float) < 0.05)),
                median_overlap=("overlap", "median"),
                median_or=("odds_ratio_approx", "median"),
                min_p=("p_value", "min"),
                median_p=("p_value", "median")
            )
            .reset_index()
        )
    else:
        gp_summary = pd.DataFrame()

    gp_summary_path = os.path.join(args.outdir, "08_step11B_gene_priority_SFARI_summary.tsv")
    gp_summary.to_csv(gp_summary_path, sep="\t", index=False)

    with open(summary_path, "w") as f:
        f.write("# Step11B full graph diffusion versus simpler baselines\n\n")
        f.write("## Purpose\n\n")
        f.write("This analysis compares full NeuroTRACE graph diffusion with direct module-stage alignment, simple module-gene by stage-gene projection, graph diffusion alone, and risk-free full graph diffusion.\n\n")

        f.write("## Input and configuration\n\n")
        f.write(f"- Nodes: `{args.nodes}`\n")
        f.write(f"- Edges: `{args.edges}`\n")
        f.write(f"- PPR alpha: {args.alpha}\n")
        f.write(f"- Stage softmax temperature: {args.temperature}\n")
        f.write(f"- Module-gene signed weight column used: `{signed_weight_col}`\n")
        f.write(f"- Signed module-gene weights detected: {signed_available}\n")
        f.write(f"- Priority gene cutoff: {args.priority_top_gene_k}\n\n")

        f.write("## Method-level summary\n\n")
        f.write(method_summary.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Main call summary\n\n")
        f.write(main_summary.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Full graph versus direct margin\n\n")
        f.write(margin_compare.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Gene-priority SFARI summary\n\n")
        if gp_summary.shape[0] > 0:
            f.write(gp_summary.to_markdown(index=False))
            f.write("\n\n")
        else:
            f.write("No SFARI-like risk sets were detected.\n\n")

        f.write("## Interpretation guide\n\n")
        f.write("- If full_graph and full_graph_no_risk preserve the expected NTM1/NTM3 late-prenatal calls while NTM2 remains non-prenatal, the graph layer supports the same biological separation as the direct alignment.\n")
        f.write("- If graph-prioritized genes show stronger or comparable SFARI enrichment than native module genes, the gene-priority layer provides added convergence beyond a raw module list.\n")
        f.write("- If direct and full graph results are highly concordant, the conservative interpretation is that graph diffusion preserves and contextualizes the direct alignment rather than replacing it.\n")

    with open(text_path, "w") as f:
        f.write("# Step11B manuscript-ready text\n\n")
        f.write("To quantify the incremental value of the graph layer, we compared full NeuroTRACE diffusion with simpler developmental baselines, including direct module-stage alignment, a module-gene by stage-gene projection, graph diffusion alone, and a risk-free full-graph model. ")
        f.write("The main interpretive criterion was whether NTM1_ASD_up and NTM3_ASD_signed retained late-prenatal top-stage alignment while NTM2_ASD_down remained separated from the prenatal-transport modules. ")
        f.write("We also compared SFARI enrichment among native module genes, full graph-prioritized genes, and no-risk graph-prioritized genes to evaluate whether graph-based prioritization added convergence beyond the original module gene lists.\n\n")
        f.write("After running this analysis, insert the exact method-level and gene-priority results from `05_step11B_overall_summary.md`.\n")

    log("Finished Step11B.")
    for p in [
        out_stage, out_main, out_method, out_gp, summary_path,
        text_path, margin_compare_path, gp_summary_path
    ]:
        log(f"Wrote: {p}")


if __name__ == "__main__":
    main()
