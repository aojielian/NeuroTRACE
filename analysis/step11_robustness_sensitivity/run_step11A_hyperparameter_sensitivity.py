#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step11A: NeuroTRACE hyperparameter sensitivity analysis

Purpose
-------
Test whether the main NeuroTRACE graph-transport conclusions are stable across:
  1. Personalized PageRank restart alpha
  2. Stage softmax temperature
  3. Edge-type scaling factors
  4. Risk-edge inclusion/removal

Primary questions
-----------------
  Q1. Do NTM1_ASD_up and NTM3_ASD_signed top200/top500 continue to favor late_prenatal?
  Q2. Does NTM2_ASD_down remain distinct from the late-prenatal transport pattern?
  Q3. Does NTM2 no-risk SFARI enrichment remain stable when risk-gene prior edges are removed?
  Q4. Are conclusions robust rather than driven by a single fixed parameter setting?

Inputs
------
Default inputs are based on the current NeuroTRACE project structure:
  - native graph nodes
  - native graph edges
  - optional output directory

Outputs
-------
  01_step11A_param_grid_stage_probabilities.tsv
  02_step11A_top_stage_stability.tsv
  03_step11A_late_prenatal_stability.tsv
  04_step11A_no_risk_SFARI_stability.tsv
  05_step11A_overall_summary.md
  06_step11A_top_stage_matrix.tsv
  07_step11A_main_call_pass_fail.tsv

Notes
-----
This script is intentionally defensive about column names because upstream graph
tables may have slightly different naming conventions.
"""

import argparse
import gzip
import math
import os
import re
import sys
import itertools
from collections import defaultdict

import numpy as np
import pandas as pd

try:
    from scipy import sparse
    from scipy.stats import hypergeom
except Exception as e:
    raise RuntimeError(
        "This script requires scipy. Please run it with the Python environment "
        "that contains scipy, numpy, and pandas."
    ) from e


# -----------------------------
# Utility functions
# -----------------------------

def log(msg):
    print(f"[Step11A] {msg}", flush=True)


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
            f"Could not find {label}. Tried: {candidates}. "
            f"Available columns: {list(df.columns)}"
        )
    return None


def zscore(x):
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return x
    sd = np.nanstd(x)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(x)
    return (x - np.nanmean(x)) / sd


def softmax(x, temperature=1.0):
    x = np.asarray(x, dtype=float)
    if temperature <= 0:
        raise ValueError("temperature must be positive")
    x = x / temperature
    x = x - np.nanmax(x)
    ex = np.exp(x)
    denom = np.nansum(ex)
    if denom == 0 or not np.isfinite(denom):
        return np.ones_like(x) / len(x)
    return ex / denom


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        q[order[i]] = prev
    return np.minimum(q, 1.0)


def normalize_edge_type(x):
    if pd.isna(x):
        return "unknown"
    x = str(x)
    x = x.strip()
    x = x.replace("-", "_")
    x = x.replace(" ", "_")
    return x


def guess_node_classes(nodes, node_col, type_col=None):
    """
    Return dict node_id -> broad class:
      gene, module, stage, risk, unknown
    """
    classes = {}
    stage_patterns = [
        "early_prenatal", "mid_prenatal", "late_prenatal",
        "infant", "child", "adolescent", "adult", "prenatal",
        "postnatal"
    ]

    for _, row in nodes.iterrows():
        node = str(row[node_col])
        low = node.lower()

        t = ""
        if type_col is not None and type_col in nodes.columns:
            t = str(row[type_col]).lower()

        cls = "unknown"

        if "module" in t or low.startswith("ntm") or "asd_up" in low or "asd_down" in low or "asd_signed" in low:
            cls = "module"
        elif "stage" in t or any(p in low for p in stage_patterns):
            cls = "stage"
        elif "risk" in t or "sfari" in low or "satterstrom" in low or "zhou" in low:
            cls = "risk"
        elif "gene" in t:
            cls = "gene"

        classes[node] = cls

    # Fallback: nodes not classified as module/stage/risk become genes.
    for node, cls in list(classes.items()):
        if cls == "unknown":
            classes[node] = "gene"

    return classes


def build_sparse_transition(edges, all_nodes, source_col, target_col, weight_col, edge_type_col, edge_scales):
    """
    Build row-normalized sparse transition matrix from undirected weighted edges.
    Returns:
      P: row-normalized sparse matrix
    """
    node_to_idx = {n: i for i, n in enumerate(all_nodes)}
    rows = []
    cols = []
    vals = []

    for _, r in edges.iterrows():
        s = str(r[source_col])
        t = str(r[target_col])
        if s not in node_to_idx or t not in node_to_idx:
            continue

        et = normalize_edge_type(r[edge_type_col]) if edge_type_col else "unknown"
        scale = edge_scales.get(et, edge_scales.get("default", 1.0))
        if scale <= 0:
            continue

        try:
            w = float(r[weight_col])
        except Exception:
            w = 1.0

        if not np.isfinite(w):
            continue

        w = abs(w) * scale
        if w <= 0:
            continue

        i = node_to_idx[s]
        j = node_to_idx[t]

        # Treat as undirected for PPR.
        rows.extend([i, j])
        cols.extend([j, i])
        vals.extend([w, w])

    n = len(all_nodes)
    A = sparse.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    row_sum = np.asarray(A.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    inv = sparse.diags(1.0 / row_sum)
    P = inv @ A
    return P


def run_ppr(P, seed_idx, alpha=0.35, max_iter=120, tol=1e-10):
    n = P.shape[0]
    e = np.zeros(n, dtype=float)
    e[seed_idx] = 1.0
    p = e.copy()
    for _ in range(max_iter):
        new_p = alpha * e + (1.0 - alpha) * (P.T @ p)
        if np.linalg.norm(new_p - p, ord=1) < tol:
            p = new_p
            break
        p = new_p
    s = p.sum()
    if s > 0:
        p = p / s
    return p


def extract_direct_module_stage(edges, source_col, target_col, weight_col, edge_type_col, module_nodes, stage_nodes):
    """
    Extract direct module-stage edge weight matrix if available.
    Missing pairs get 0.
    """
    module_set = set(module_nodes)
    stage_set = set(stage_nodes)
    direct = defaultdict(float)

    for _, r in edges.iterrows():
        et = normalize_edge_type(r[edge_type_col]) if edge_type_col else ""
        if "module_stage" not in et and not ("module" in et and "stage" in et):
            continue

        s = str(r[source_col])
        t = str(r[target_col])
        try:
            w = float(r[weight_col])
        except Exception:
            w = 0.0

        if s in module_set and t in stage_set:
            direct[(s, t)] += w
        elif t in module_set and s in stage_set:
            direct[(t, s)] += w

    return direct


def extract_risk_sets_from_edges(edges, source_col, target_col, edge_type_col, node_classes):
    """
    Extract risk-set genes from risk_gene_prior edges in the original graph.
    Returns dict risk_set_name -> set(gene_ids).
    """
    risk_sets = defaultdict(set)

    for _, r in edges.iterrows():
        et = normalize_edge_type(r[edge_type_col]) if edge_type_col else ""
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


def infer_module_family(module_name):
    x = module_name.lower()
    if "ntm1" in x:
        return "NTM1_ASD_up"
    if "ntm2" in x:
        return "NTM2_ASD_down"
    if "ntm3" in x:
        return "NTM3_ASD_signed"
    return module_name


def infer_top_n(module_name):
    x = str(module_name)
    m = re.search(r"top[_-]?(\d+)", x, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))
    m = re.search(r"_(\d+)$", x)
    if m:
        return int(m.group(1))
    return np.nan


def hypergeom_enrichment(top_genes, risk_genes, universe_genes):
    top_genes = set(top_genes) & set(universe_genes)
    risk_genes = set(risk_genes) & set(universe_genes)
    universe_genes = set(universe_genes)

    M = len(universe_genes)
    n = len(risk_genes)
    N = len(top_genes)
    x = len(top_genes & risk_genes)

    if M == 0 or n == 0 or N == 0:
        return {
            "overlap": x,
            "top_n": N,
            "risk_n": n,
            "universe_n": M,
            "odds_ratio_approx": np.nan,
            "p_value": np.nan
        }

    p = hypergeom.sf(x - 1, M, n, N)

    a = x
    b = N - x
    c = n - x
    d = M - a - b - c
    odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

    return {
        "overlap": x,
        "top_n": N,
        "risk_n": n,
        "universe_n": M,
        "odds_ratio_approx": odds,
        "p_value": p
    }


# -----------------------------
# Main
# -----------------------------

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
    parser.add_argument("--top_gene_k", type=int, default=500)
    parser.add_argument("--max_iter", type=int, default=120)
    parser.add_argument("--tol", type=float, default=1e-10)
    parser.add_argument("--main_module_regex", default=r"NTM[123].*(top[_-]?(200|500)|200|500)")
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
        log("No edge weight column found. Creating unit weights.")
        weight_col = "__unit_weight__"
        edges[weight_col] = 1.0

    if edge_type_col is None:
        log("No edge type column found. Creating unknown edge type.")
        edge_type_col = "__edge_type__"
        edges[edge_type_col] = "unknown"

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
    risk_nodes = [n for n in all_nodes if node_classes.get(n) == "risk"]

    main_re = re.compile(args.main_module_regex, flags=re.IGNORECASE)
    main_module_nodes = [m for m in module_nodes if main_re.search(m)]
    if len(main_module_nodes) == 0:
        main_module_nodes = module_nodes

    log(f"Detected nodes: total={len(all_nodes)}, genes={len(gene_nodes)}, modules={len(module_nodes)}, stages={len(stage_nodes)}, risks={len(risk_nodes)}")
    log(f"Main modules selected: {len(main_module_nodes)}")
    log(f"Stage nodes: {stage_nodes}")

    direct_ms = extract_direct_module_stage(edges, source_col, target_col, weight_col, edge_type_col, module_nodes, stage_nodes)
    risk_sets = extract_risk_sets_from_edges(edges, source_col, target_col, edge_type_col, node_classes)

    sfari_keys = [k for k in risk_sets.keys() if "sfari" in k.lower()]
    log(f"Detected risk sets from graph: {len(risk_sets)}")
    if sfari_keys:
        log("SFARI-like risk sets: " + ", ".join(sfari_keys))
    else:
        log("No SFARI-like risk set detected from risk edges. SFARI stability table will be empty.")

    # Hyperparameter grid.
    alphas = [0.15, 0.25, 0.35, 0.50, 0.70]
    temperatures = [0.50, 0.75, 1.00, 1.50]
    gene_gene_scales = [0.20, 0.35, 0.50]
    module_stage_scales = [0.75, 1.00, 1.25, 1.50]
    risk_scales = [0.00, 0.10, 0.30, 0.50]

    # Keep module-gene and stage-gene stable in this first grid to avoid exploding combinations.
    module_gene_scale = 1.00
    stage_gene_scale = 0.80

    param_records = []
    stage_prob_records = []
    sfari_records = []

    total_grid = len(alphas) * len(temperatures) * len(gene_gene_scales) * len(module_stage_scales) * len(risk_scales)
    log(f"Total parameter settings: {total_grid}")

    grid_id = 0

    for alpha, temp, gg_scale, ms_scale, risk_scale in itertools.product(
        alphas, temperatures, gene_gene_scales, module_stage_scales, risk_scales
    ):
        grid_id += 1

        edge_scales = {
            "gene_gene": gg_scale,
            "gene_knn": gg_scale,
            "knn": gg_scale,
            "module_gene": module_gene_scale,
            "stage_gene": stage_gene_scale,
            "module_stage": ms_scale,
            "risk_gene_prior": risk_scale,
            "risk_gene": risk_scale,
            "risk": risk_scale,
            "default": 1.0
        }

        log(f"Grid {grid_id}/{total_grid}: alpha={alpha}, temp={temp}, gene_gene={gg_scale}, module_stage={ms_scale}, risk={risk_scale}")

        P = build_sparse_transition(
            edges=edges,
            all_nodes=all_nodes,
            source_col=source_col,
            target_col=target_col,
            weight_col=weight_col,
            edge_type_col=edge_type_col,
            edge_scales=edge_scales
        )

        stage_indices = [node_to_idx[s] for s in stage_nodes]
        gene_indices = [node_to_idx[g] for g in gene_nodes]

        for module in main_module_nodes:
            if module not in node_to_idx:
                continue

            ppr = run_ppr(P, node_to_idx[module], alpha=alpha, max_iter=args.max_iter, tol=args.tol)

            graph_stage = np.array([ppr[i] for i in stage_indices], dtype=float)
            graph_z = zscore(graph_stage)

            # Direct module-stage edge component if available.
            direct_stage = np.array([direct_ms.get((module, s), 0.0) for s in stage_nodes], dtype=float)
            direct_z = zscore(direct_stage)

            # Combine graph diffusion and direct stage information.
            # The module-stage edge scale already changes the graph component; this 50:50 combination
            # avoids dropping direct alignment when graph scores are very small.
            combined = 0.50 * graph_z + 0.50 * direct_z

            probs = softmax(combined, temperature=temp)
            top_idx = int(np.argmax(probs))
            top_stage = stage_nodes[top_idx]
            expected_idx = None
            expected_prob = np.nan
            for i, s in enumerate(stage_nodes):
                if s.lower() == args.expected_stage.lower():
                    expected_idx = i
                    expected_prob = probs[i]

            module_family = infer_module_family(module)
            top_n = infer_top_n(module)

            for s, raw_g, dz, gz, comb, prob in zip(stage_nodes, direct_stage, direct_z, graph_z, combined, probs):
                stage_prob_records.append({
                    "grid_id": grid_id,
                    "alpha": alpha,
                    "temperature": temp,
                    "gene_gene_scale": gg_scale,
                    "module_gene_scale": module_gene_scale,
                    "stage_gene_scale": stage_gene_scale,
                    "module_stage_scale": ms_scale,
                    "risk_gene_prior_scale": risk_scale,
                    "module_node": module,
                    "module_family": module_family,
                    "top_n": top_n,
                    "stage": s,
                    "direct_stage_weight": raw_g,
                    "direct_stage_z": dz,
                    "graph_stage_z": gz,
                    "combined_stage_score": comb,
                    "stage_probability": prob,
                    "top_stage": top_stage,
                    "is_top_stage": int(s == top_stage),
                    "expected_stage": args.expected_stage,
                    "expected_stage_probability": expected_prob,
                    "expected_stage_is_top": int(top_stage.lower() == args.expected_stage.lower())
                })

            # No-risk SFARI stability for NTM2 only, evaluated when risk edges are removed.
            if risk_scale == 0.0 and "ntm2" in module.lower() and sfari_keys:
                gene_scores = pd.DataFrame({
                    "gene": gene_nodes,
                    "ppr_score": [ppr[i] for i in gene_indices]
                }).sort_values("ppr_score", ascending=False)

                top_genes = gene_scores["gene"].head(args.top_gene_k).tolist()

                for risk_name in sfari_keys:
                    enr = hypergeom_enrichment(top_genes, risk_sets[risk_name], gene_nodes)
                    sfari_records.append({
                        "grid_id": grid_id,
                        "alpha": alpha,
                        "temperature": temp,
                        "gene_gene_scale": gg_scale,
                        "module_gene_scale": module_gene_scale,
                        "stage_gene_scale": stage_gene_scale,
                        "module_stage_scale": ms_scale,
                        "risk_gene_prior_scale": risk_scale,
                        "module_node": module,
                        "module_family": module_family,
                        "top_n": top_n,
                        "risk_set": risk_name,
                        "top_gene_k": args.top_gene_k,
                        **enr
                    })

    stage_df = pd.DataFrame(stage_prob_records)
    stage_out = os.path.join(args.outdir, "01_step11A_param_grid_stage_probabilities.tsv")
    stage_df.to_csv(stage_out, sep="\t", index=False)

    # Top-stage stability summary.
    top_df = (
        stage_df[stage_df["is_top_stage"] == 1]
        .groupby(["module_node", "module_family", "top_n", "top_stage"], dropna=False)
        .size()
        .reset_index(name="n_parameter_settings")
    )
    total_by_module = (
        stage_df[stage_df["is_top_stage"] == 1]
        .groupby(["module_node"], dropna=False)
        .size()
        .reset_index(name="n_total_parameter_settings")
    )
    top_df = top_df.merge(total_by_module, on="module_node", how="left")
    top_df["top_stage_fraction"] = top_df["n_parameter_settings"] / top_df["n_total_parameter_settings"]
    top_out = os.path.join(args.outdir, "02_step11A_top_stage_stability.tsv")
    top_df.to_csv(top_out, sep="\t", index=False)

    # Expected-stage stability summary.
    expected_df = (
        stage_df[stage_df["stage"].str.lower() == args.expected_stage.lower()]
        .groupby(["module_node", "module_family", "top_n"], dropna=False)
        .agg(
            n_parameter_settings=("grid_id", "nunique"),
            expected_stage_top_fraction=("expected_stage_is_top", "mean"),
            expected_stage_probability_median=("stage_probability", "median"),
            expected_stage_probability_mean=("stage_probability", "mean"),
            expected_stage_probability_min=("stage_probability", "min"),
            expected_stage_probability_max=("stage_probability", "max")
        )
        .reset_index()
    )
    exp_out = os.path.join(args.outdir, "03_step11A_late_prenatal_stability.tsv")
    expected_df.to_csv(exp_out, sep="\t", index=False)

    # SFARI no-risk stability table.
    if len(sfari_records) > 0:
        sfari_df = pd.DataFrame(sfari_records)
        sfari_df["fdr"] = np.nan
        for risk_name, idx in sfari_df.groupby("risk_set").groups.items():
            p = sfari_df.loc[idx, "p_value"].values
            valid = np.isfinite(p)
            q = np.full(len(p), np.nan)
            if valid.sum() > 0:
                q[valid] = bh_fdr(p[valid])
            sfari_df.loc[idx, "fdr"] = q

        sfari_out = os.path.join(args.outdir, "04_step11A_no_risk_SFARI_stability.tsv")
        sfari_df.to_csv(sfari_out, sep="\t", index=False)
    else:
        sfari_df = pd.DataFrame()
        sfari_out = os.path.join(args.outdir, "04_step11A_no_risk_SFARI_stability.tsv")
        sfari_df.to_csv(sfari_out, sep="\t", index=False)

    # Matrix for quick heatmap-like inspection.
    matrix = (
        stage_df[stage_df["is_top_stage"] == 1]
        .pivot_table(
            index=["module_node", "module_family", "top_n"],
            columns="top_stage",
            values="grid_id",
            aggfunc="nunique",
            fill_value=0
        )
        .reset_index()
    )
    matrix_out = os.path.join(args.outdir, "06_step11A_top_stage_matrix.tsv")
    matrix.to_csv(matrix_out, sep="\t", index=False)

    # PASS/FAIL criteria for main conclusions.
    pass_records = []
    for _, r in expected_df.iterrows():
        fam = str(r["module_family"])
        topn = r["top_n"]
        frac = float(r["expected_stage_top_fraction"])
        medp = float(r["expected_stage_probability_median"])

        criterion = "not_applicable"
        status = "INFO"

        if fam in ["NTM1_ASD_up", "NTM3_ASD_signed"] and topn in [200, 500]:
            criterion = f"{args.expected_stage} top-stage fraction >= 0.80"
            status = "PASS" if frac >= 0.80 else "REVIEW"
        elif fam == "NTM2_ASD_down" and topn in [200, 500]:
            criterion = f"{args.expected_stage} should not dominate; top-stage fraction <= 0.50"
            status = "PASS" if frac <= 0.50 else "REVIEW"

        pass_records.append({
            "module_node": r["module_node"],
            "module_family": fam,
            "top_n": topn,
            "criterion": criterion,
            "status": status,
            "expected_stage_top_fraction": frac,
            "expected_stage_probability_median": medp,
            "n_parameter_settings": int(r["n_parameter_settings"])
        })

    pass_df = pd.DataFrame(pass_records)
    pass_out = os.path.join(args.outdir, "07_step11A_main_call_pass_fail.tsv")
    pass_df.to_csv(pass_out, sep="\t", index=False)

    # Overall summary markdown.
    summary_out = os.path.join(args.outdir, "05_step11A_overall_summary.md")
    with open(summary_out, "w") as f:
        f.write("# Step11A NeuroTRACE hyperparameter sensitivity summary\n\n")
        f.write("## Purpose\n\n")
        f.write("This analysis evaluates whether the main NeuroTRACE graph-transport conclusions are stable across PPR restart alpha, stage softmax temperature, gene-gene edge scaling, module-stage edge scaling, and risk-gene prior edge scaling.\n\n")

        f.write("## Input graph\n\n")
        f.write(f"- Nodes: `{args.nodes}`\n")
        f.write(f"- Edges: `{args.edges}`\n")
        f.write(f"- Total nodes: {len(all_nodes)}\n")
        f.write(f"- Gene nodes: {len(gene_nodes)}\n")
        f.write(f"- Module nodes: {len(module_nodes)}\n")
        f.write(f"- Stage nodes: {len(stage_nodes)}\n")
        f.write(f"- Risk nodes: {len(risk_nodes)}\n")
        f.write(f"- Parameter settings: {total_grid}\n\n")

        f.write("## Main late-prenatal stability table\n\n")
        if expected_df.shape[0] > 0:
            f.write(expected_df.to_markdown(index=False))
            f.write("\n\n")
        else:
            f.write("No expected-stage records were generated. Check stage node naming.\n\n")

        f.write("## PASS/REVIEW table\n\n")
        if pass_df.shape[0] > 0:
            f.write(pass_df.to_markdown(index=False))
            f.write("\n\n")

        f.write("## No-risk SFARI stability\n\n")
        if sfari_df.shape[0] > 0:
            sfari_summary = (
                sfari_df
                .groupby(["module_node", "module_family", "top_n", "risk_set"], dropna=False)
                .agg(
                    n_tests=("grid_id", "count"),
                    median_overlap=("overlap", "median"),
                    median_or=("odds_ratio_approx", "median"),
                    min_p=("p_value", "min"),
                    median_p=("p_value", "median"),
                    fdr_lt_0_05_fraction=("fdr", lambda x: np.mean(np.asarray(x, dtype=float) < 0.05))
                )
                .reset_index()
            )
            f.write(sfari_summary.to_markdown(index=False))
            f.write("\n\n")
        else:
            f.write("No SFARI-like risk set was detected from graph risk edges, or no NTM2 no-risk tests were generated.\n\n")

        f.write("## Interpretation guide\n\n")
        f.write("- For NTM1_ASD_up and NTM3_ASD_signed top200/top500, robust support is defined as late_prenatal remaining the top stage in >=80% of parameter settings.\n")
        f.write("- For NTM2_ASD_down top200/top500, robust separation is defined as late_prenatal not dominating the top-stage calls.\n")
        f.write("- For no-risk SFARI stability, robust support is defined as median OR > 1 and FDR < 0.05 in a large fraction of risk-free parameter settings.\n")
        f.write("- REVIEW does not mean failure; it means the module should be inspected before manuscript-level interpretation.\n")

    log("Finished Step11A.")
    log(f"Wrote: {stage_out}")
    log(f"Wrote: {top_out}")
    log(f"Wrote: {exp_out}")
    log(f"Wrote: {sfari_out}")
    log(f"Wrote: {summary_out}")
    log(f"Wrote: {matrix_out}")
    log(f"Wrote: {pass_out}")


if __name__ == "__main__":
    main()
