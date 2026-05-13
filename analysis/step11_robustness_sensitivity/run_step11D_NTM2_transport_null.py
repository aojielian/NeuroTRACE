#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step11D: Matched-random transport calibration for NTM2_ASD_down.

Purpose
-------
Test whether NTM2_ASD_down is specifically separated from late-prenatal transport
and whether its adolescence/adulthood transport is stronger than matched random
modules.

This step complements the existing Step10B transport null for NTM1/NTM3.

Outputs
-------
01_step11D_NTM2_observed_transport.tsv
02_step11D_NTM2_random_transport.tsv.gz
03_step11D_NTM2_transport_null_summary.tsv
04_step11D_NTM2_top_stage_null_summary.tsv
05_step11D_overall_summary.md
"""

import argparse
import os
import re
import gzip
import math
import numpy as np
import pandas as pd
from collections import defaultdict

try:
    from scipy import sparse
except Exception as e:
    raise RuntimeError("This script requires scipy, numpy, and pandas.") from e


def log(x):
    print(f"[Step11D] {x}", flush=True)


def ensure_dir(x):
    os.makedirs(x, exist_ok=True)


def read_table(path):
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    return pd.read_csv(path, sep="\t", low_memory=False)


def find_col(df, candidates, required=True):
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    if required:
        raise ValueError(f"Missing column from candidates {candidates}. Available: {list(df.columns)}")
    return None


def norm_type(x):
    if pd.isna(x):
        return "unknown"
    return str(x).strip().replace("-", "_").replace(" ", "_")


def clean_stage(x):
    return str(x).replace("stage:", "")


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


def edge_weight(row, col):
    try:
        v = float(row[col])
        if np.isfinite(v):
            return v
    except Exception:
        pass
    return 0.0


def infer_top_n(x):
    m = re.search(r"top[_|\\-]?(\\d+)", str(x), flags=re.IGNORECASE)
    return int(m.group(1)) if m else np.nan


def infer_module_family(x):
    x = str(x)
    if "NTM1" in x:
        return "NTM1_ASD_up"
    if "NTM2" in x:
        return "NTM2_ASD_down"
    if "NTM3" in x:
        return "NTM3_ASD_signed"
    return x


def guess_classes(nodes, node_col, type_col=None):
    out = {}
    stage_tokens = ["early_prenatal", "mid_prenatal", "late_prenatal", "adolescence", "adulthood", "adult"]
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


def detect_signed_col(edges, weight_col):
    candidates = [
        "signed_weight", "module_weight", "raw_weight", "weight_signed",
        "signed_score", "gene_weight", "ntm_weight", "score_signed", "beta_weight"
    ]
    lower = {c.lower(): c for c in edges.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return weight_col


def extract_module_gene(edges, source_col, target_col, weight_col, signed_col, edge_type_col, module_nodes, gene_nodes):
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
            "weight_abs": abs(edge_weight(r, weight_col)),
            "weight_signed": edge_weight(r, signed_col)
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

    for _, r in edges.iterrows():
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
        rows += [i, j]
        cols += [j, i]
        vals += [w, w]

    for s, t, w, et in add_edges:
        if s not in node_to_idx or t not in node_to_idx:
            continue
        scale = scales.get(et, scales.get("default", 1.0))
        w = abs(float(w)) * scale
        if w <= 0:
            continue
        i = node_to_idx[s]
        j = node_to_idx[t]
        rows += [i, j]
        cols += [j, i]
        vals += [w, w]

    n = len(all_nodes)
    A = sparse.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    row_sum = np.asarray(A.sum(axis=1)).ravel()
    row_sum[row_sum == 0] = 1.0
    return sparse.diags(1.0 / row_sum) @ A


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
        a = np.array([gene_weights[g] for g in common], dtype=float)
        b = np.array([sg[g] for g in common], dtype=float)
        denom = math.sqrt(np.sum(a * a)) * math.sqrt(np.sum(b * b))
        scores.append(float(np.sum(a * b) / denom) if denom > 0 else 0.0)
    return np.array(scores), overlaps


def summarize_probs(module, label, stage_nodes, direct_scores, graph_scores, temperature):
    combined = 0.5 * zscore(direct_scores) + 0.5 * zscore(graph_scores)
    probs = softmax(combined, temperature=temperature)

    stage_clean = [clean_stage(s) for s in stage_nodes]
    order = np.argsort(probs)[::-1]
    top_stage = stage_clean[order[0]]

    def prob_of(stage_name):
        return float(sum(p for s, p in zip(stage_clean, probs) if s == stage_name))

    late = prob_of("late_prenatal")
    mid = prob_of("mid_prenatal")
    adolescence = prob_of("adolescence")
    adulthood = prob_of("adulthood")
    adult_like = adolescence + adulthood

    return {
        "module_node": module,
        "module_family": infer_module_family(module),
        "top_n": infer_top_n(module),
        "record_type": label,
        "top_stage": top_stage,
        "late_prenatal_probability": late,
        "mid_prenatal_probability": mid,
        "adolescence_probability": adolescence,
        "adulthood_probability": adulthood,
        "adolescence_plus_adulthood_probability": adult_like,
        "top_probability": float(probs[order[0]]),
        "second_probability": float(probs[order[1]]) if len(order) > 1 else np.nan,
        "top_margin": float(probs[order[0]] - probs[order[1]]) if len(order) > 1 else np.nan
    }


def make_gene_bins(gene_nodes, stage_gene, stage_nodes, n_bins=5):
    rows = []
    for g in gene_nodes:
        vals = np.array([stage_gene.get(st, {}).get(g, 0.0) for st in stage_nodes], dtype=float)
        rows.append({
            "gene": g,
            "stage_mean": float(np.mean(vals)),
            "stage_var": float(np.var(vals))
        })
    df = pd.DataFrame(rows)
    df["mean_bin"] = pd.qcut(df["stage_mean"].rank(method="first"), q=n_bins, labels=False, duplicates="drop")
    df["var_bin"] = pd.qcut(df["stage_var"].rank(method="first"), q=n_bins, labels=False, duplicates="drop")
    df["bin"] = df["mean_bin"].astype(str) + "_" + df["var_bin"].astype(str)
    return df


def sample_matched_genes(obs_genes, bin_df, rng):
    gene_to_bin = dict(zip(bin_df["gene"], bin_df["bin"]))
    bin_to_genes = {b: sub["gene"].tolist() for b, sub in bin_df.groupby("bin")}

    sampled = []
    for g in obs_genes:
        b = gene_to_bin.get(g, None)
        pool = bin_to_genes.get(b, bin_df["gene"].tolist())
        sampled.append(rng.choice(pool))

    # If duplicated, refill from global pool to keep approximate size.
    sampled_unique = list(dict.fromkeys(sampled))
    if len(sampled_unique) < len(obs_genes):
        pool = list(set(bin_df["gene"].tolist()) - set(sampled_unique))
        need = len(obs_genes) - len(sampled_unique)
        if len(pool) >= need:
            sampled_unique += rng.choice(pool, size=need, replace=False).tolist()
        else:
            sampled_unique += rng.choice(bin_df["gene"].tolist(), size=need, replace=True).tolist()
    return sampled_unique[:len(obs_genes)]


def empirical_p_high(obs, null):
    null = np.asarray(null, dtype=float)
    return (1.0 + np.sum(null >= obs)) / (len(null) + 1.0)


def empirical_p_low(obs, null):
    null = np.asarray(null, dtype=float)
    return (1.0 + np.sum(null <= obs)) / (len(null) + 1.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nodes", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step04_graph_construction/results/09_step04C_native_graph_nodes.tsv")
    parser.add_argument("--edges", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step04_graph_construction/results/10_step04C_native_graph_edges.tsv.gz")
    parser.add_argument("--outdir", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step11_robustness_sensitivity/results")
    parser.add_argument("--n_random", type=int, default=1000)
    parser.add_argument("--alpha", type=float, default=0.35)
    parser.add_argument("--temperature", type=float, default=0.75)
    parser.add_argument("--seed", type=int, default=20260509)
    args = parser.parse_args()

    ensure_dir(args.outdir)
    rng = np.random.default_rng(args.seed)

    nodes = read_table(args.nodes)
    edges = read_table(args.edges)

    node_col = find_col(nodes, ["node_id", "node", "id", "name", "node_name"])
    type_col = find_col(nodes, ["node_type", "type", "class", "category"], required=False)
    source_col = find_col(edges, ["source", "src", "from", "node1", "from_node"])
    target_col = find_col(edges, ["target", "dst", "to", "node2", "to_node"])
    weight_col = find_col(edges, ["weight", "edge_weight", "w", "score"], required=False)
    edge_type_col = find_col(edges, ["edge_type", "type", "relation", "edge_class"], required=False)

    if weight_col is None:
        edges["__unit_weight__"] = 1.0
        weight_col = "__unit_weight__"
    if edge_type_col is None:
        edges["__edge_type__"] = "unknown"
        edge_type_col = "__edge_type__"

    signed_col = detect_signed_col(edges, weight_col)

    nodes[node_col] = nodes[node_col].astype(str)
    edges[source_col] = edges[source_col].astype(str)
    edges[target_col] = edges[target_col].astype(str)
    edges[edge_type_col] = edges[edge_type_col].map(norm_type)

    classes = guess_classes(nodes, node_col, type_col)
    all_nodes = list(nodes[node_col])
    node_to_idx = {n: i for i, n in enumerate(all_nodes)}

    module_nodes = [n for n in all_nodes if classes.get(n) == "module"]
    stage_nodes = [n for n in all_nodes if classes.get(n) == "stage"]
    gene_nodes = [n for n in all_nodes if classes.get(n) == "gene"]
    stage_indices = [node_to_idx[s] for s in stage_nodes]

    ntm2_modules = [m for m in module_nodes if "NTM2" in m and ("top200" in m or "top500" in m)]
    ntm2_modules = sorted(ntm2_modules, key=infer_top_n)

    log(f"NTM2 modules: {ntm2_modules}")
    log(f"Stages: {[clean_stage(s) for s in stage_nodes]}")
    log(f"Genes: {len(gene_nodes)}")
    log(f"Random modules per observed module: {args.n_random}")

    module_gene = extract_module_gene(edges, source_col, target_col, weight_col, signed_col, edge_type_col, module_nodes, gene_nodes)
    stage_gene = extract_stage_gene(edges, source_col, target_col, weight_col, edge_type_col, stage_nodes, gene_nodes)
    direct_ms = extract_direct_module_stage(edges, source_col, target_col, weight_col, edge_type_col, module_nodes, stage_nodes)

    bin_df = make_gene_bins(gene_nodes, stage_gene, stage_nodes, n_bins=5)

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

    observed_records = []
    random_records = []

    for module in ntm2_modules:
        obs_edges = module_gene.get(module, [])
        if len(obs_edges) == 0:
            log(f"No module-gene edges for {module}; skipping")
            continue

        obs_genes = [x["gene"] for x in obs_edges]
        obs_weights_signed = np.array([x["weight_signed"] for x in obs_edges], dtype=float)
        obs_weights_abs = np.abs(np.array([x["weight_abs"] for x in obs_edges], dtype=float))
        obs_gene_weights = {g: w for g, w in zip(obs_genes, obs_weights_signed)}

        drop_idx = [x["edge_index"] for x in obs_edges]
        # Also drop module-stage edges so random modules do not inherit observed direct stage edges.
        for idx, r in edges.iterrows():
            et = norm_type(r[edge_type_col])
            if "module_stage" in et or ("module" in et and "stage" in et):
                s = str(r[source_col])
                t = str(r[target_col])
                if s == module or t == module:
                    drop_idx.append(idx)

        # Observed run: keep observed module-gene edges, but use no-risk graph.
        P_obs = build_transition(edges, all_nodes, source_col, target_col, weight_col, edge_type_col,
                                 no_risk_scales, drop_edge_indices=[])
        ppr_obs = run_ppr(P_obs, node_to_idx[module], alpha=args.alpha)
        graph_obs = np.array([ppr_obs[i] for i in stage_indices], dtype=float)
        direct_obs = np.array([direct_ms.get((module, st), 0.0) for st in stage_nodes], dtype=float)
        obs_summary = summarize_probs(module, "observed", stage_nodes, direct_obs, graph_obs, args.temperature)
        observed_records.append(obs_summary)

        # Random modules.
        for rep in range(1, args.n_random + 1):
            rand_genes = sample_matched_genes(obs_genes, bin_df, rng)

            # Preserve observed weight distribution but permute assignment.
            perm_signed = rng.permutation(obs_weights_signed)
            perm_abs = np.abs(rng.permutation(obs_weights_abs))

            rand_gene_weights = {g: w for g, w in zip(rand_genes, perm_signed)}
            direct_rand, overlap_counts = weighted_stage_projection(rand_gene_weights, stage_gene, stage_nodes)

            add_edges = []
            for g, w in zip(rand_genes, perm_abs):
                add_edges.append((module, g, float(w), "module_gene"))

            # Add random direct module-stage edges derived from projection scores.
            # Shift scores to nonnegative for adjacency; direct signed signal is still used in combined scoring.
            direct_shift = direct_rand - np.min(direct_rand)
            if np.max(direct_shift) > 0:
                direct_shift = direct_shift / np.max(direct_shift)
            else:
                direct_shift = np.zeros_like(direct_shift)

            for st, w in zip(stage_nodes, direct_shift):
                add_edges.append((module, st, float(w), "module_stage"))

            P_rand = build_transition(edges, all_nodes, source_col, target_col, weight_col, edge_type_col,
                                      no_risk_scales, drop_edge_indices=drop_idx, add_edges=add_edges)
            ppr_rand = run_ppr(P_rand, node_to_idx[module], alpha=args.alpha)
            graph_rand = np.array([ppr_rand[i] for i in stage_indices], dtype=float)

            rand_summary = summarize_probs(module, "random", stage_nodes, direct_rand, graph_rand, args.temperature)
            rand_summary["replicate"] = rep
            rand_summary["n_module_genes"] = len(obs_genes)
            rand_summary["mean_stage_overlap"] = float(np.mean(overlap_counts))
            random_records.append(rand_summary)

            if rep % 200 == 0:
                log(f"{module}: random replicate {rep}/{args.n_random}")

    obs_df = pd.DataFrame(observed_records)
    rand_df = pd.DataFrame(random_records)

    obs_out = os.path.join(args.outdir, "01_step11D_NTM2_observed_transport.tsv")
    rand_out = os.path.join(args.outdir, "02_step11D_NTM2_random_transport.tsv.gz")
    obs_df.to_csv(obs_out, sep="\t", index=False)
    rand_df.to_csv(rand_out, sep="\t", index=False, compression="gzip")

    summary_records = []
    top_records = []

    for module in ntm2_modules:
        obs = obs_df[obs_df["module_node"] == module].iloc[0]
        null = rand_df[rand_df["module_node"] == module].copy()

        rec = {
            "module_node": module,
            "module_family": infer_module_family(module),
            "top_n": infer_top_n(module),
            "observed_top_stage": obs["top_stage"],
            "observed_late_prenatal_probability": obs["late_prenatal_probability"],
            "null_late_prenatal_mean": null["late_prenatal_probability"].mean(),
            "null_late_prenatal_median": null["late_prenatal_probability"].median(),
            "empirical_p_late_prenatal_low": empirical_p_low(obs["late_prenatal_probability"], null["late_prenatal_probability"]),
            "empirical_p_late_prenatal_high": empirical_p_high(obs["late_prenatal_probability"], null["late_prenatal_probability"]),
            "observed_adolescence_plus_adulthood_probability": obs["adolescence_plus_adulthood_probability"],
            "null_adolescence_plus_adulthood_mean": null["adolescence_plus_adulthood_probability"].mean(),
            "null_adolescence_plus_adulthood_median": null["adolescence_plus_adulthood_probability"].median(),
            "empirical_p_adolescence_plus_adulthood_high": empirical_p_high(obs["adolescence_plus_adulthood_probability"], null["adolescence_plus_adulthood_probability"]),
            "observed_top_margin": obs["top_margin"],
            "null_top_margin_median": null["top_margin"].median()
        }

        # Decision: strict adult-like if adult/adolescence high is significant;
        # non-prenatal if late-prenatal is not top and empirical low is small.
        if (obs["top_stage"] in ["adolescence", "adulthood"]) and rec["empirical_p_adolescence_plus_adulthood_high"] <= 0.05:
            rec["decision"] = "adult_like_enriched"
        elif (obs["top_stage"] in ["adolescence", "adulthood"]) and rec["empirical_p_late_prenatal_low"] <= 0.05:
            rec["decision"] = "late_prenatal_depleted"
        elif obs["top_stage"] in ["adolescence", "adulthood"]:
            rec["decision"] = "non_prenatal_descriptive"
        else:
            rec["decision"] = "review"

        summary_records.append(rec)

        top_counts = null["top_stage"].value_counts().reset_index()
        top_counts.columns = ["null_top_stage", "n_random"]
        top_counts["module_node"] = module
        top_counts["top_n"] = infer_top_n(module)
        top_counts["fraction"] = top_counts["n_random"] / len(null)
        top_records.append(top_counts)

    summary_df = pd.DataFrame(summary_records)
    top_df = pd.concat(top_records, ignore_index=True) if top_records else pd.DataFrame()

    summary_out = os.path.join(args.outdir, "03_step11D_NTM2_transport_null_summary.tsv")
    top_out = os.path.join(args.outdir, "04_step11D_NTM2_top_stage_null_summary.tsv")
    summary_df.to_csv(summary_out, sep="\t", index=False)
    top_df.to_csv(top_out, sep="\t", index=False)

    md_out = os.path.join(args.outdir, "05_step11D_overall_summary.md")
    with open(md_out, "w") as f:
        f.write("# Step11D NTM2 transport null calibration\n\n")
        f.write("## Purpose\n\n")
        f.write("This analysis tests whether NTM2_ASD_down is specifically separated from late-prenatal transport and whether its adolescence/adulthood transport exceeds matched random modules.\n\n")
        f.write("## Configuration\n\n")
        f.write(f"- Random modules per NTM2 threshold: {args.n_random}\n")
        f.write("- Matching variables: BrainSpan graph genes binned by stage-profile mean and variance\n")
        f.write("- Weight structure: observed NTM2 signed/absolute weight distributions permuted across matched random genes\n")
        f.write("- Graph mode: no-risk graph\n")
        f.write(f"- PPR alpha: {args.alpha}\n")
        f.write(f"- Stage softmax temperature: {args.temperature}\n\n")

        f.write("## Summary table\n\n")
        f.write(summary_df.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Null top-stage distribution\n\n")
        f.write(top_df.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Interpretation guide\n\n")
        f.write("- `adult_like_enriched`: NTM2 top stage is adolescence/adulthood and adolescence+adulthood probability is higher than matched random modules.\n")
        f.write("- `late_prenatal_depleted`: NTM2 top stage is adolescence/adulthood and late-prenatal probability is lower than matched random modules.\n")
        f.write("- `non_prenatal_descriptive`: NTM2 is non-prenatal, but the empirical null does not support a strong enrichment/depletion claim.\n")
        f.write("- `review`: inspect before manuscript interpretation.\n")

    log("Finished Step11D.")
    for p in [obs_out, rand_out, summary_out, top_out, md_out]:
        log(f"Wrote: {p}")


if __name__ == "__main__":
    main()
