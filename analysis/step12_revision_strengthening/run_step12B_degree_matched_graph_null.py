#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step12B: Degree-matched / neighborhood-density-matched graph nulls.

Purpose
-------
Reviewer-facing robustness analysis to test whether NeuroTRACE stage transport
and NTM2 SFARI convergence are driven by graph degree or neighborhood density.

Main tests
----------
1. Degree-matched developmental transport null for NTM1_ASD_up and
   NTM3_ASD_signed top200/top500.
2. Degree-matched non-prenatal / adult-like null for NTM2_ASD_down top200/top500.
3. Degree-matched graph-priority SFARI null for NTM2_ASD_down top200/top500.

Matching variables
------------------
For each observed module gene, sample a replacement gene matched by:
  - graph weighted degree bin
  - local gene-gene neighborhood density bin
  - developmental stage-profile mean bin
  - developmental stage-profile variance bin

Outputs
-------
01_step12B_observed_transport.tsv
02_step12B_degree_matched_random_transport.tsv.gz
03_step12B_transport_null_summary.tsv
04_step12B_NTM2_observed_SFARI.tsv
05_step12B_NTM2_degree_matched_random_SFARI.tsv.gz
06_step12B_NTM2_SFARI_null_summary.tsv
07_step12B_decision_table.tsv
08_step12B_overall_summary.md
"""

import argparse
import os
import re
import math
from collections import defaultdict

import numpy as np
import pandas as pd

from scipy import sparse
from scipy.stats import hypergeom


def log(x):
    print(f"[Step12B] {x}", flush=True)


def ensure_dir(x):
    os.makedirs(x, exist_ok=True)


def read_table(path):
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def find_col(df, candidates, required=False):
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    if required:
        raise ValueError(f"Missing required column from {candidates}; available={list(df.columns)}")
    return None


def clean_gene(x):
    x = str(x).strip()
    x = x.replace("gene:", "")
    x = x.replace("GENE:", "")
    return x


def clean_stage(x):
    return str(x).replace("stage:", "")


def norm_type(x):
    if pd.isna(x):
        return "unknown"
    return str(x).strip().replace("-", "_").replace(" ", "_")


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    if (not np.isfinite(sd)) or sd == 0:
        return np.zeros_like(x)
    return (x - np.nanmean(x)) / sd


def softmax(x, temperature=0.75):
    x = np.asarray(x, dtype=float) / temperature
    x = x - np.nanmax(x)
    ex = np.exp(x)
    denom = np.nansum(ex)
    if denom <= 0 or not np.isfinite(denom):
        return np.ones(len(x)) / len(x)
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
    return int(m.group(1)) if m else np.nan


def edge_weight(row, col):
    try:
        v = float(row[col])
        if np.isfinite(v):
            return v
    except Exception:
        pass
    return 0.0


def guess_classes(nodes, node_col, type_col=None):
    out = {}
    stage_tokens = [
        "early_prenatal", "mid_prenatal", "late_prenatal",
        "adolescence", "adulthood", "adult", "prenatal"
    ]
    for _, r in nodes.iterrows():
        node = str(r[node_col])
        low = node.lower()
        t = str(r[type_col]).lower() if type_col and type_col in nodes.columns else ""

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

        out[node] = cls

    return out


def extract_module_gene(edges, source_col, target_col, weight_col, edge_type_col, module_nodes, gene_nodes):
    module_set = set(module_nodes)
    gene_set = set(gene_nodes)
    out = defaultdict(list)

    for idx, r in edges.iterrows():
        et = norm_type(r[edge_type_col])
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

        out[m].append({
            "edge_index": idx,
            "module": m,
            "gene": g,
            "weight_abs": abs(edge_weight(r, weight_col))
        })

    return out


def extract_stage_gene(edges, source_col, target_col, weight_col, edge_type_col, stage_nodes, gene_nodes):
    stage_set = set(stage_nodes)
    gene_set = set(gene_nodes)
    out = defaultdict(dict)

    for _, r in edges.iterrows():
        et = norm_type(r[edge_type_col])
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

        out[st][g] = out[st].get(g, 0.0) + edge_weight(r, weight_col)

    return out


def extract_direct_module_stage(edges, source_col, target_col, weight_col, edge_type_col, module_nodes, stage_nodes):
    module_set = set(module_nodes)
    stage_set = set(stage_nodes)
    out = defaultdict(float)

    for idx, r in edges.iterrows():
        et = norm_type(r[edge_type_col])
        if "module_stage" not in et and not ("module" in et and "stage" in et):
            continue

        s = str(r[source_col])
        t = str(r[target_col])
        w = edge_weight(r, weight_col)

        if s in module_set and t in stage_set:
            out[(s, t)] += w
        elif t in module_set and s in stage_set:
            out[(t, s)] += w

    return out


def extract_risk_sets(edges, source_col, target_col, edge_type_col, classes):
    risk_sets = defaultdict(set)

    for _, r in edges.iterrows():
        et = norm_type(r[edge_type_col])
        if "risk" not in et:
            continue

        s = str(r[source_col])
        t = str(r[target_col])
        cs = classes.get(s, "unknown")
        ct = classes.get(t, "unknown")

        if cs == "risk" and ct == "gene":
            risk_sets[s].add(t)
        elif ct == "risk" and cs == "gene":
            risk_sets[t].add(s)

    return dict(risk_sets)


def build_gene_graph_metrics(edges, source_col, target_col, weight_col, edge_type_col, gene_nodes):
    gene_set = set(gene_nodes)
    degree = defaultdict(float)
    unweighted_degree = defaultdict(int)
    neighbor_weight = defaultdict(list)

    for _, r in edges.iterrows():
        et = norm_type(r[edge_type_col])
        if not ("gene_gene" in et or "gene_knn" in et or et == "knn" or ("gene" in et and "gene" in et)):
            # Still allow generic gene-gene by source/target class below.
            pass

        s = str(r[source_col])
        t = str(r[target_col])
        if s not in gene_set or t not in gene_set:
            continue

        w = abs(edge_weight(r, weight_col))
        if w <= 0:
            continue

        degree[s] += w
        degree[t] += w
        unweighted_degree[s] += 1
        unweighted_degree[t] += 1
        neighbor_weight[s].append(w)
        neighbor_weight[t].append(w)

    records = []
    for g in gene_nodes:
        wlist = neighbor_weight.get(g, [])
        records.append({
            "gene": g,
            "weighted_degree": degree.get(g, 0.0),
            "unweighted_degree": unweighted_degree.get(g, 0),
            "local_density": float(np.mean(wlist)) if len(wlist) else 0.0
        })
    return pd.DataFrame(records)


def stage_profile_metrics(gene_nodes, stage_gene, stage_nodes):
    records = []
    for g in gene_nodes:
        vals = []
        for st in stage_nodes:
            vals.append(stage_gene.get(st, {}).get(g, 0.0))
        vals = np.asarray(vals, dtype=float)
        records.append({
            "gene": g,
            "stage_profile_mean": float(np.mean(vals)),
            "stage_profile_var": float(np.var(vals)),
            "stage_profile_absmean": float(np.mean(np.abs(vals)))
        })
    return pd.DataFrame(records)


def add_bins(df, n_bins=5):
    df = df.copy()
    for col in ["weighted_degree", "local_density", "stage_profile_mean", "stage_profile_var"]:
        vals = pd.to_numeric(df[col], errors="coerce").fillna(0)
        if len(np.unique(vals)) <= 1:
            df[col + "_bin"] = 0
        else:
            df[col + "_bin"] = pd.qcut(vals.rank(method="first"), q=n_bins, labels=False, duplicates="drop")
    df["match_bin"] = (
        df["weighted_degree_bin"].astype(str) + "_" +
        df["local_density_bin"].astype(str) + "_" +
        df["stage_profile_mean_bin"].astype(str) + "_" +
        df["stage_profile_var_bin"].astype(str)
    )
    return df


def sample_matched_genes(obs_genes, gene_features, rng):
    gene_to_bin = dict(zip(gene_features["gene"], gene_features["match_bin"]))
    bin_to_genes = {b: sub["gene"].tolist() for b, sub in gene_features.groupby("match_bin")}
    all_genes = gene_features["gene"].tolist()

    sampled = []
    for g in obs_genes:
        b = gene_to_bin.get(g, None)
        pool = bin_to_genes.get(b, all_genes)
        if len(pool) == 0:
            pool = all_genes
        sampled.append(rng.choice(pool))

    # Prefer unique genes while preserving size where possible.
    unique_sampled = list(dict.fromkeys(sampled))
    if len(unique_sampled) < len(obs_genes):
        pool = list(set(all_genes) - set(unique_sampled))
        need = len(obs_genes) - len(unique_sampled)
        if len(pool) >= need:
            unique_sampled += rng.choice(pool, size=need, replace=False).tolist()
        else:
            unique_sampled += rng.choice(all_genes, size=need, replace=True).tolist()

    return unique_sampled[:len(obs_genes)]


def weighted_stage_projection(gene_weights, stage_gene, stage_nodes):
    scores = []
    overlaps = []
    for st in stage_nodes:
        sg = stage_gene.get(st, {})
        common = sorted(set(gene_weights.keys()) & set(sg.keys()))
        overlaps.append(len(common))

        if len(common) == 0:
            scores.append(0.0)
            continue

        a = np.asarray([gene_weights[g] for g in common], dtype=float)
        b = np.asarray([sg[g] for g in common], dtype=float)
        denom = math.sqrt(np.sum(a * a)) * math.sqrt(np.sum(b * b))
        scores.append(float(np.sum(a * b) / denom) if denom > 0 else 0.0)

    return np.asarray(scores), overlaps


def summarize_transport(module, stage_nodes, scores, temperature=0.75):
    stage_clean = [clean_stage(x) for x in stage_nodes]
    z = zscore(scores)
    p = softmax(z, temperature=temperature)
    order = np.argsort(p)[::-1]

    def prob_of(s):
        return float(sum(prob for st, prob in zip(stage_clean, p) if st == s))

    return {
        "module_node": module,
        "module_family": infer_module_family(module),
        "top_n": infer_top_n(module),
        "top_stage": stage_clean[order[0]],
        "top_probability": float(p[order[0]]),
        "second_stage": stage_clean[order[1]] if len(order) > 1 else "NA",
        "second_probability": float(p[order[1]]) if len(order) > 1 else np.nan,
        "top_margin": float(p[order[0]] - p[order[1]]) if len(order) > 1 else np.nan,
        "late_prenatal_probability": prob_of("late_prenatal"),
        "mid_prenatal_probability": prob_of("mid_prenatal"),
        "adolescence_probability": prob_of("adolescence"),
        "adulthood_probability": prob_of("adulthood"),
        "adolescence_plus_adulthood_probability": prob_of("adolescence") + prob_of("adulthood")
    }


def build_transition(edges, all_nodes, source_col, target_col, weight_col, edge_type_col,
                     scales, drop_edge_indices=None, add_edges=None):
    if drop_edge_indices is None:
        drop_edge_indices = set()
    else:
        drop_edge_indices = set(drop_edge_indices)

    if add_edges is None:
        add_edges = []

    node_to_idx = {n: i for i, n in enumerate(all_nodes)}

    rows, cols, vals = [], [], []

    for idx, r in edges.iterrows():
        if idx in drop_edge_indices:
            continue

        s = str(r[source_col])
        t = str(r[target_col])
        if s not in node_to_idx or t not in node_to_idx:
            continue

        et = norm_type(r[edge_type_col])
        scale = scales.get(et, scales.get("default", 1.0))
        if scale <= 0:
            continue

        w = abs(edge_weight(r, weight_col)) * scale
        if w <= 0:
            continue

        i = node_to_idx[s]
        j = node_to_idx[t]

        rows.extend([i, j])
        cols.extend([j, i])
        vals.extend([w, w])

    for s, t, w, et in add_edges:
        if s not in node_to_idx or t not in node_to_idx:
            continue
        scale = scales.get(et, scales.get("default", 1.0))
        w = abs(float(w)) * scale
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
        new_p = alpha * e + (1.0 - alpha) * (P.T @ p)
        if np.linalg.norm(new_p - p, ord=1) < tol:
            p = new_p
            break
        p = new_p
    s = p.sum()
    if s > 0:
        p = p / s
    return p


def hypergeom_enrichment(top_genes, risk_genes, universe):
    top_genes = set(top_genes) & set(universe)
    risk_genes = set(risk_genes) & set(universe)
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


def empirical_p_high(obs, null):
    null = np.asarray(null, dtype=float)
    return (1 + np.sum(null >= obs)) / (len(null) + 1)


def empirical_p_low(obs, null):
    null = np.asarray(null, dtype=float)
    return (1 + np.sum(null <= obs)) / (len(null) + 1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step04_graph_construction/results/09_step04C_native_graph_nodes.tsv")
    parser.add_argument("--edges", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step04_graph_construction/results/10_step04C_native_graph_edges.tsv.gz")
    parser.add_argument("--outdir", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step12_revision_strengthening/results/step12B_degree_matched_graph_null")
    parser.add_argument("--n_random", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260509)
    parser.add_argument("--temperature", type=float, default=0.75)
    parser.add_argument("--alpha", type=float, default=0.35)
    parser.add_argument("--priority_top_gene_k", type=int, default=500)
    args = parser.parse_args()

    ensure_dir(args.outdir)
    rng = np.random.default_rng(args.seed)

    log(f"Reading nodes: {args.nodes}")
    nodes = read_table(args.nodes)
    log(f"Reading edges: {args.edges}")
    edges = read_table(args.edges)

    node_col = find_col(nodes, ["node_id", "node", "id", "name", "node_name"], required=True)
    type_col = find_col(nodes, ["node_type", "type", "class", "category"])
    source_col = find_col(edges, ["source", "src", "from", "node1", "from_node"], required=True)
    target_col = find_col(edges, ["target", "dst", "to", "node2", "to_node"], required=True)
    weight_col = find_col(edges, ["weight", "edge_weight", "w", "score"])
    edge_type_col = find_col(edges, ["edge_type", "type", "relation", "edge_class"])

    if weight_col is None:
        edges["__unit_weight__"] = 1.0
        weight_col = "__unit_weight__"
    if edge_type_col is None:
        edges["__edge_type__"] = "unknown"
        edge_type_col = "__edge_type__"

    nodes[node_col] = nodes[node_col].astype(str)
    edges[source_col] = edges[source_col].astype(str)
    edges[target_col] = edges[target_col].astype(str)
    edges[edge_type_col] = edges[edge_type_col].map(norm_type)

    classes = guess_classes(nodes, node_col, type_col)

    all_nodes = list(nodes[node_col].astype(str))
    node_to_idx = {n: i for i, n in enumerate(all_nodes)}

    module_nodes = [n for n in all_nodes if classes.get(n) == "module"]
    stage_nodes = [n for n in all_nodes if classes.get(n) == "stage"]
    gene_nodes = [n for n in all_nodes if classes.get(n) == "gene"]

    main_modules = [
        m for m in module_nodes
        if re.search(r"NTM[123].*top(200|500)", m)
    ]
    if len(main_modules) == 0:
        main_modules = [
            m for m in module_nodes
            if "NTM" in m and ("200" in m or "500" in m)
        ]

    main_modules = sorted(main_modules)

    log(f"Genes: {len(gene_nodes)}")
    log(f"Stages: {[clean_stage(x) for x in stage_nodes]}")
    log(f"Main modules: {main_modules}")

    module_gene = extract_module_gene(edges, source_col, target_col, weight_col, edge_type_col, main_modules, gene_nodes)
    stage_gene = extract_stage_gene(edges, source_col, target_col, weight_col, edge_type_col, stage_nodes, gene_nodes)
    direct_ms = extract_direct_module_stage(edges, source_col, target_col, weight_col, edge_type_col, main_modules, stage_nodes)
    risk_sets = extract_risk_sets(edges, source_col, target_col, edge_type_col, classes)
    sfari_sets = {k: v for k, v in risk_sets.items() if "sfari" in k.lower()}

    log(f"SFARI-like risk sets: {len(sfari_sets)}")

    gene_graph = build_gene_graph_metrics(edges, source_col, target_col, weight_col, edge_type_col, gene_nodes)
    gene_stage = stage_profile_metrics(gene_nodes, stage_gene, stage_nodes)
    gene_features = gene_graph.merge(gene_stage, on="gene", how="left")
    gene_features = add_bins(gene_features, n_bins=5)

    feature_out = os.path.join(args.outdir, "00_step12B_gene_degree_density_stage_features.tsv")
    gene_features.to_csv(feature_out, sep="\t", index=False)

    # Risk-free graph scales for gene-priority test.
    no_risk_scales = {
        "gene_gene": 0.35,
        "gene_knn": 0.35,
        "knn": 0.35,
        "module_gene": 1.00,
        "stage_gene": 0.80,
        "module_stage": 1.25,
        "risk_gene_prior": 0.0,
        "risk_gene": 0.0,
        "risk": 0.0,
        "default": 1.0
    }

    # Drop original module-gene and module-stage edges for randomized module insertion.
    module_set = set(main_modules)
    drop_module_edges = set()
    for idx, r in edges.iterrows():
        s = str(r[source_col])
        t = str(r[target_col])
        et = norm_type(r[edge_type_col])
        if (s in module_set or t in module_set) and (
            "module_gene" in et or "module_stage" in et or ("module" in et and ("gene" in et or "stage" in et))
        ):
            drop_module_edges.add(idx)

    observed_transport = []
    random_transport = []

    observed_sfari = []
    random_sfari = []

    stage_indices = [node_to_idx[s] for s in stage_nodes]
    gene_indices = [node_to_idx[g] for g in gene_nodes]

    for module in main_modules:
        log(f"Processing module: {module}")
        mg = module_gene.get(module, [])
        if len(mg) == 0:
            log(f"WARNING: no module-gene edges for {module}")
            continue

        obs_genes = [x["gene"] for x in mg]
        obs_weights = np.asarray([x["weight_abs"] for x in mg], dtype=float)
        obs_gene_weights = {g: w for g, w in zip(obs_genes, obs_weights)}

        # Observed direct projection based on module genes and stage-gene profiles.
        obs_direct, obs_overlap = weighted_stage_projection(obs_gene_weights, stage_gene, stage_nodes)
        obs_t = summarize_transport(module, stage_nodes, obs_direct, temperature=args.temperature)
        obs_t["record_type"] = "observed_direct_projection"
        obs_t["n_module_genes"] = len(obs_genes)
        obs_t["mean_stage_overlap"] = float(np.mean(obs_overlap))
        observed_transport.append(obs_t)

        # Observed no-risk graph-priority SFARI.
        add_edges_obs = [(module, g, float(w), "module_gene") for g, w in zip(obs_genes, obs_weights)]
        # Add module-stage edges from direct projection, shifted nonnegative.
        shift = obs_direct - np.min(obs_direct)
        if np.max(shift) > 0:
            shift = shift / np.max(shift)
        for st, w in zip(stage_nodes, shift):
            add_edges_obs.append((module, st, float(w), "module_stage"))

        P_obs = build_transition(
            edges, all_nodes, source_col, target_col, weight_col, edge_type_col,
            no_risk_scales, drop_edge_indices=drop_module_edges, add_edges=add_edges_obs
        )
        ppr_obs = run_ppr(P_obs, node_to_idx[module], alpha=args.alpha)
        gene_scores_obs = pd.DataFrame({
            "gene": gene_nodes,
            "score": [ppr_obs[i] for i in gene_indices]
        }).sort_values("score", ascending=False)

        top_obs = gene_scores_obs["gene"].head(args.priority_top_gene_k).tolist()

        for risk_name, risk_genes in sfari_sets.items():
            x, N, n, M, odds, p = hypergeom_enrichment(top_obs, risk_genes, gene_nodes)
            observed_sfari.append({
                "module_node": module,
                "module_family": infer_module_family(module),
                "top_n": infer_top_n(module),
                "risk_set": risk_name,
                "overlap": x,
                "tested_gene_n": N,
                "risk_gene_n": n,
                "universe_n": M,
                "odds_ratio_approx": odds,
                "p_value": p
            })

        # Random modules.
        for rep in range(1, args.n_random + 1):
            if rep % 200 == 0:
                log(f"{module}: random {rep}/{args.n_random}")

            rand_genes = sample_matched_genes(obs_genes, gene_features, rng)
            rand_weights = rng.permutation(obs_weights)

            rand_gene_weights = {g: w for g, w in zip(rand_genes, rand_weights)}

            rand_direct, rand_overlap = weighted_stage_projection(rand_gene_weights, stage_gene, stage_nodes)
            rand_t = summarize_transport(module, stage_nodes, rand_direct, temperature=args.temperature)
            rand_t["record_type"] = "degree_density_stage_matched_random"
            rand_t["replicate"] = rep
            rand_t["n_module_genes"] = len(rand_genes)
            rand_t["mean_stage_overlap"] = float(np.mean(rand_overlap))
            random_transport.append(rand_t)

            # Only do graph-priority SFARI null for NTM2 to limit runtime and focus claim.
            if "NTM2" in module and len(sfari_sets) > 0:
                add_edges = [(module, g, float(w), "module_gene") for g, w in zip(rand_genes, rand_weights)]

                shift_r = rand_direct - np.min(rand_direct)
                if np.max(shift_r) > 0:
                    shift_r = shift_r / np.max(shift_r)
                for st, w in zip(stage_nodes, shift_r):
                    add_edges.append((module, st, float(w), "module_stage"))

                P_rand = build_transition(
                    edges, all_nodes, source_col, target_col, weight_col, edge_type_col,
                    no_risk_scales, drop_edge_indices=drop_module_edges, add_edges=add_edges
                )
                ppr_rand = run_ppr(P_rand, node_to_idx[module], alpha=args.alpha)

                gene_scores_rand = pd.DataFrame({
                    "gene": gene_nodes,
                    "score": [ppr_rand[i] for i in gene_indices]
                }).sort_values("score", ascending=False)

                top_rand = gene_scores_rand["gene"].head(args.priority_top_gene_k).tolist()

                for risk_name, risk_genes in sfari_sets.items():
                    x, N, n, M, odds, p = hypergeom_enrichment(top_rand, risk_genes, gene_nodes)
                    random_sfari.append({
                        "module_node": module,
                        "module_family": infer_module_family(module),
                        "top_n": infer_top_n(module),
                        "replicate": rep,
                        "risk_set": risk_name,
                        "overlap": x,
                        "tested_gene_n": N,
                        "risk_gene_n": n,
                        "universe_n": M,
                        "odds_ratio_approx": odds,
                        "p_value": p
                    })

    obs_df = pd.DataFrame(observed_transport)
    rand_df = pd.DataFrame(random_transport)
    obs_sfari_df = pd.DataFrame(observed_sfari)
    rand_sfari_df = pd.DataFrame(random_sfari)

    obs_out = os.path.join(args.outdir, "01_step12B_observed_transport.tsv")
    rand_out = os.path.join(args.outdir, "02_step12B_degree_matched_random_transport.tsv.gz")
    obs_df.to_csv(obs_out, sep="\t", index=False)
    rand_df.to_csv(rand_out, sep="\t", index=False, compression="gzip")

    obs_sfari_out = os.path.join(args.outdir, "04_step12B_NTM2_observed_SFARI.tsv")
    rand_sfari_out = os.path.join(args.outdir, "05_step12B_NTM2_degree_matched_random_SFARI.tsv.gz")
    obs_sfari_df.to_csv(obs_sfari_out, sep="\t", index=False)
    rand_sfari_df.to_csv(rand_sfari_out, sep="\t", index=False, compression="gzip")

    # Transport null summary.
    trans_records = []
    for _, obs in obs_df.iterrows():
        module = obs["module_node"]
        null = rand_df[rand_df["module_node"].eq(module)]

        if null.shape[0] == 0:
            continue

        fam = obs["module_family"]
        top_n = obs["top_n"]

        rec = {
            "module_node": module,
            "module_family": fam,
            "top_n": top_n,
            "observed_top_stage": obs["top_stage"],
            "observed_late_prenatal_probability": obs["late_prenatal_probability"],
            "null_late_prenatal_mean": null["late_prenatal_probability"].mean(),
            "null_late_prenatal_median": null["late_prenatal_probability"].median(),
            "empirical_p_late_prenatal_high": empirical_p_high(obs["late_prenatal_probability"], null["late_prenatal_probability"]),
            "empirical_p_late_prenatal_low": empirical_p_low(obs["late_prenatal_probability"], null["late_prenatal_probability"]),
            "observed_adolescence_plus_adulthood_probability": obs["adolescence_plus_adulthood_probability"],
            "null_adolescence_plus_adulthood_mean": null["adolescence_plus_adulthood_probability"].mean(),
            "null_adolescence_plus_adulthood_median": null["adolescence_plus_adulthood_probability"].median(),
            "empirical_p_adolescence_plus_adulthood_high": empirical_p_high(
                obs["adolescence_plus_adulthood_probability"],
                null["adolescence_plus_adulthood_probability"]
            )
        }

        if fam in ["NTM1_ASD_up", "NTM3_ASD_signed"]:
            rec["primary_test"] = "late_prenatal_high"
            rec["primary_empirical_p"] = rec["empirical_p_late_prenatal_high"]
        elif fam == "NTM2_ASD_down":
            rec["primary_test"] = "late_prenatal_low_or_adultlike_high"
            rec["primary_empirical_p"] = min(
                rec["empirical_p_late_prenatal_low"],
                rec["empirical_p_adolescence_plus_adulthood_high"]
            )
        else:
            rec["primary_test"] = "NA"
            rec["primary_empirical_p"] = np.nan

        trans_records.append(rec)

    trans_summary = pd.DataFrame(trans_records)
    if trans_summary.shape[0] > 0:
        trans_summary["primary_fdr"] = bh_fdr(trans_summary["primary_empirical_p"].values)

    trans_summary_out = os.path.join(args.outdir, "03_step12B_transport_null_summary.tsv")
    trans_summary.to_csv(trans_summary_out, sep="\t", index=False)

    # SFARI null summary for NTM2.
    sfari_records = []
    if obs_sfari_df.shape[0] > 0 and rand_sfari_df.shape[0] > 0:
        for _, obs in obs_sfari_df.iterrows():
            if obs["module_family"] != "NTM2_ASD_down":
                continue

            null = rand_sfari_df[
                rand_sfari_df["module_node"].eq(obs["module_node"]) &
                rand_sfari_df["risk_set"].eq(obs["risk_set"])
            ]

            if null.shape[0] == 0:
                continue

            rec = {
                "module_node": obs["module_node"],
                "module_family": obs["module_family"],
                "top_n": obs["top_n"],
                "risk_set": obs["risk_set"],
                "observed_overlap": obs["overlap"],
                "null_overlap_mean": null["overlap"].mean(),
                "null_overlap_median": null["overlap"].median(),
                "empirical_p_overlap_high": empirical_p_high(obs["overlap"], null["overlap"]),
                "observed_or": obs["odds_ratio_approx"],
                "null_or_mean": null["odds_ratio_approx"].mean(),
                "null_or_median": null["odds_ratio_approx"].median(),
                "empirical_p_or_high": empirical_p_high(obs["odds_ratio_approx"], null["odds_ratio_approx"]),
                "observed_hypergeom_p": obs["p_value"]
            }
            sfari_records.append(rec)

    sfari_summary = pd.DataFrame(sfari_records)
    if sfari_summary.shape[0] > 0:
        sfari_summary["overlap_fdr"] = bh_fdr(sfari_summary["empirical_p_overlap_high"].values)
        sfari_summary["or_fdr"] = bh_fdr(sfari_summary["empirical_p_or_high"].values)

    sfari_summary_out = os.path.join(args.outdir, "06_step12B_NTM2_SFARI_null_summary.tsv")
    sfari_summary.to_csv(sfari_summary_out, sep="\t", index=False)

    # Decision table.
    decisions = []

    for _, row in trans_summary.iterrows():
        fam = row["module_family"]
        top_n = row["top_n"]
        fdr = row.get("primary_fdr", np.nan)

        if fam in ["NTM1_ASD_up", "NTM3_ASD_signed"] and top_n in [200, 500]:
            status = "PASS" if np.isfinite(fdr) and fdr < 0.05 else "REVIEW"
            interpretation = "late-prenatal transport exceeds degree/density/stage-profile matched random modules"
        elif fam == "NTM2_ASD_down" and top_n in [200, 500]:
            status = "PASS" if np.isfinite(fdr) and fdr < 0.05 else "REVIEW"
            interpretation = "NTM2 is depleted for late-prenatal transport or enriched for adolescence/adulthood versus matched random modules"
        else:
            status = "INFO"
            interpretation = "exploratory"

        decisions.append({
            "analysis": "transport_degree_matched_null",
            "module_family": fam,
            "top_n": top_n,
            "status": status,
            "primary_fdr": fdr,
            "interpretation": interpretation
        })

    if sfari_summary.shape[0] > 0:
        # Collapse NTM2 SFARI summary by top_n.
        for top_n, sub in sfari_summary.groupby("top_n"):
            n_sets = sub.shape[0]
            n_pass = int(np.sum(sub["overlap_fdr"] < 0.05))
            decisions.append({
                "analysis": "NTM2_SFARI_degree_matched_null",
                "module_family": "NTM2_ASD_down",
                "top_n": top_n,
                "status": "PASS" if n_pass >= max(1, math.ceil(n_sets * 0.5)) else "REVIEW",
                "primary_fdr": sub["overlap_fdr"].min(),
                "interpretation": f"{n_pass}/{n_sets} SFARI-like sets exceed degree/density matched graph-priority null"
            })

    decision_df = pd.DataFrame(decisions)
    decision_out = os.path.join(args.outdir, "07_step12B_decision_table.tsv")
    decision_df.to_csv(decision_out, sep="\t", index=False)

    # Markdown summary.
    md_out = os.path.join(args.outdir, "08_step12B_overall_summary.md")
    with open(md_out, "w") as f:
        f.write("# Step12B degree-matched graph null summary\n\n")

        f.write("## Purpose\n\n")
        f.write("This analysis tests whether NeuroTRACE transport and NTM2 gene-priority SFARI convergence are driven by graph degree or local neighborhood density. Random modules were matched to observed module genes by weighted degree, local density, developmental stage-profile mean, and developmental stage-profile variance.\n\n")

        f.write("## Configuration\n\n")
        f.write(f"- Random modules per observed module: {args.n_random}\n")
        f.write(f"- PPR alpha for gene-priority null: {args.alpha}\n")
        f.write(f"- Stage softmax temperature: {args.temperature}\n")
        f.write("- Matching variables: weighted degree, local density, stage-profile mean, stage-profile variance\n")
        f.write("- Graph-priority SFARI test: NTM2 only, risk-free graph setting\n\n")

        f.write("## Transport null summary\n\n")
        if trans_summary.shape[0] > 0:
            f.write(trans_summary.to_markdown(index=False))
        else:
            f.write("No transport summary generated.")
        f.write("\n\n")

        f.write("## NTM2 SFARI degree-matched null summary\n\n")
        if sfari_summary.shape[0] > 0:
            f.write(sfari_summary.to_markdown(index=False))
        else:
            f.write("No SFARI null summary generated.")
        f.write("\n\n")

        f.write("## Decision table\n\n")
        if decision_df.shape[0] > 0:
            f.write(decision_df.to_markdown(index=False))
        else:
            f.write("No decision table generated.")
        f.write("\n\n")

        f.write("## Interpretation guide\n\n")
        f.write("- PASS for transport means observed stage probability remains stronger than degree/density/stage-profile matched random modules.\n")
        f.write("- PASS for NTM2 SFARI means observed graph-prioritized SFARI overlap exceeds the degree/density-matched graph-priority null for at least half of SFARI-like sets.\n")
        f.write("- REVIEW should be interpreted as a topology-sensitive boundary, not necessarily failure of the primary expression/variance-matched analyses.\n")

    log("Finished Step12B")
    for p in [
        feature_out,
        obs_out,
        rand_out,
        trans_summary_out,
        obs_sfari_out,
        rand_sfari_out,
        sfari_summary_out,
        decision_out,
        md_out
    ]:
        log(f"Wrote: {p}")


if __name__ == "__main__":
    main()
