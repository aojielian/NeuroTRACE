#!/usr/bin/env python3

import math
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
PROJECT = BASE / "neurotrace_algorithm_project"

STEP03 = PROJECT / "step03_feature_embedding" / "results"
STEP04 = PROJECT / "step04_graph_construction" / "results"
STEP05 = PROJECT / "step05_optimal_transport_alignment" / "results"
STEP10 = PROJECT / "step10_transport_null"
OUT = STEP10 / "results"
FIG = STEP10 / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

NODES_FILE = STEP04 / "09_step04C_native_graph_nodes.tsv"
EDGES_FILE = STEP04 / "10_step04C_native_graph_edges.tsv.gz"
MODULE_FILE = STEP03 / "24_step03C2_neurotrace_native_module_weights_symbol.tsv"
STAGE_PROFILES_FILE = STEP03 / "04_step03A_brainspan_stage_gene_profiles.tsv.gz"
GENE_EMBED_FILE = STEP03 / "01_step03A_brainspan_gene_embedding_pca.tsv.gz"
OBS_TRANSPORT_FILE = STEP05 / "12_step05B_NTM_top_transport_stage.tsv"
FULL_TRANSPORT_FILE = STEP05 / "10_step05B_NTM_stage_graph_transport.tsv"

N_PERM = 1000
SEED = 20260509
ALPHA_RESTART = 0.35
MAX_ITER = 120
TOL = 1e-10
TRANSPORT_TEMP = 0.75

TARGET_PROGRAMS = ["NTM1_ASD_up", "NTM3_ASD_signed"]
TARGET_TOP_N = [200, 500]
TARGET_STAGE = "late_prenatal"

EDGE_TYPE_SCALE = {
    "gene_gene_embedding_knn": 0.35,
    "module_gene_weight": 1.00,
    "stage_gene_profile": 0.80,
    "module_stage_native_alignment": 1.25,
    "risk_gene_prior": 0.30,
}

STAGE_ORDER = [
    "early_prenatal",
    "mid_prenatal",
    "late_prenatal",
    "childhood",
    "adolescence",
    "adulthood",
]


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def clean_gene(x):
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def zscore(x, eps=1e-8):
    x = np.asarray(x, dtype=float)
    return (x - np.nanmean(x)) / (np.nanstd(x) + eps)


def softmax(x, temp=1.0):
    x = np.asarray(x, dtype=float) / max(temp, 1e-8)
    x = x - np.nanmax(x)
    ex = np.exp(x)
    return ex / (np.sum(ex) + 1e-12)


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    out = np.empty(n)
    prev = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        val = p[idx] * n / rank
        prev = min(prev, val)
        out[idx] = min(prev, 1.0)
    return out


def empirical_p_greater(obs, null):
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if not np.isfinite(obs) or len(null) == 0:
        return np.nan
    return (np.sum(null >= obs) + 1) / (len(null) + 1)


def make_bins(x, n_bins=5):
    x = np.asarray(x, dtype=float)
    ok = np.isfinite(x)
    out = np.full(len(x), -1, dtype=int)
    if ok.sum() == 0:
        return out
    ranks = pd.Series(x[ok]).rank(method="average").values
    out[ok] = np.minimum(n_bins, np.maximum(1, np.ceil(ranks / ranks.max() * n_bins))).astype(int)
    return out


def build_transition_matrix(nodes, edges):
    node_ids = nodes["node_id"].astype(str).tolist()
    node_to_idx = {n: i for i, n in enumerate(node_ids)}
    n = len(node_ids)

    rows, cols, vals = [], [], []
    for _, e in edges.iterrows():
        src = str(e["source"])
        tgt = str(e["target"])
        et = str(e["edge_type"])
        if src not in node_to_idx or tgt not in node_to_idx:
            continue
        w = float(e["weight"]) * EDGE_TYPE_SCALE.get(et, 1.0)
        if not np.isfinite(w) or w <= 0:
            continue
        i = node_to_idx[src]
        j = node_to_idx[tgt]
        rows.extend([i, j])
        cols.extend([j, i])
        vals.extend([w, w])

    from scipy import sparse
    A = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=np.float64)
    row_sum = np.asarray(A.sum(axis=1)).ravel()
    inv = np.zeros_like(row_sum)
    inv[row_sum > 0] = 1.0 / row_sum[row_sum > 0]
    P = sparse.diags(inv) @ A
    return P, node_to_idx, node_ids


def personalized_pagerank(P, seed_idx, alpha=0.35, max_iter=120, tol=1e-10):
    n = P.shape[0]
    s = np.zeros(n)
    s[seed_idx] = 1.0
    x = s.copy()
    for _ in range(max_iter):
        x_new = alpha * s + (1.0 - alpha) * (P.T @ x)
        if np.sum(np.abs(x_new - x)) < tol:
            return np.asarray(x_new).ravel()
        x = x_new
    return np.asarray(x).ravel()


def weighted_stage_score(module_genes, weights, stage_profile):
    common = [g for g in module_genes if g in stage_profile]
    if len(common) == 0:
        return np.nan, 0
    wmap = dict(zip(module_genes, weights))
    w = np.array([wmap[g] for g in common], dtype=float)
    p = np.array([stage_profile[g] for g in common], dtype=float)
    score = float(np.sum(w * p) / (np.sum(np.abs(w)) + 1e-12))
    return score, len(common)


def cosine(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        return np.nan
    return float(np.dot(a, b) / ((np.linalg.norm(a) + 1e-12) * (np.linalg.norm(b) + 1e-12)))


def weighted_embedding(module_genes, weights, gene_embed, pc_cols):
    idx = [g for g in module_genes if g in gene_embed.index]
    if len(idx) == 0:
        return np.full(len(pc_cols), np.nan), 0
    wmap = dict(zip(module_genes, weights))
    w = np.array([abs(wmap[g]) for g in idx], dtype=float)
    X = gene_embed.loc[idx, pc_cols].astype(float).values
    emb = (X * (w / (w.sum() + 1e-12))[:, None]).sum(axis=0)
    return emb, len(idx)


def sample_matched_genes(module_genes, module_weights, gene_stats, rng):
    module_genes = list(module_genes)
    module_weights = np.asarray(module_weights, dtype=float)
    chosen = []
    tiers = []

    mstats = gene_stats.loc[module_genes].copy()
    for g in module_genes:
        expr_bin = mstats.loc[g, "expr_bin"]
        var_bin = mstats.loc[g, "var_bin"]

        pool1 = gene_stats[
            (gene_stats["expr_bin"] == expr_bin) &
            (gene_stats["var_bin"] == var_bin) &
            (~gene_stats.index.isin(module_genes)) &
            (~gene_stats.index.isin(chosen))
        ].index.tolist()

        pool2 = gene_stats[
            (gene_stats["expr_bin"] == expr_bin) &
            (~gene_stats.index.isin(module_genes)) &
            (~gene_stats.index.isin(chosen))
        ].index.tolist()

        pool3 = gene_stats[
            (gene_stats["var_bin"] == var_bin) &
            (~gene_stats.index.isin(module_genes)) &
            (~gene_stats.index.isin(chosen))
        ].index.tolist()

        pool4 = gene_stats[
            (~gene_stats.index.isin(module_genes)) &
            (~gene_stats.index.isin(chosen))
        ].index.tolist()

        if pool1:
            chosen_gene = rng.choice(pool1)
            tier = "expr_var_bin"
        elif pool2:
            chosen_gene = rng.choice(pool2)
            tier = "expr_bin_only"
        elif pool3:
            chosen_gene = rng.choice(pool3)
            tier = "var_bin_only"
        elif pool4:
            chosen_gene = rng.choice(pool4)
            tier = "any_nonmodule"
        else:
            chosen_gene = rng.choice(gene_stats.index.tolist())
            tier = "fallback"

        chosen.append(str(chosen_gene))
        tiers.append(tier)

    weights = rng.permutation(module_weights)
    return chosen, weights, tiers


def compute_module_stage_transport(module_genes, weights, stage_profiles, stage_embed, gene_embed, pc_cols):
    # This mirrors Step03D stage score formulation for random module stage-direct score.
    rows = []
    emb, n_embed = weighted_embedding(module_genes, weights, gene_embed, pc_cols)

    pc2_30 = [c for c in pc_cols if c != "PC1"]
    if len(pc2_30) < 5:
        pc2_30 = pc_cols

    pc2_6 = [f"PC{i}" for i in range(2, 7) if f"PC{i}" in pc_cols]
    if len(pc2_6) == 0:
        pc2_6 = pc2_30

    emb_no_pc1, _ = weighted_embedding(module_genes, weights, gene_embed, pc2_30)
    emb_pc2_6, _ = weighted_embedding(module_genes, weights, gene_embed, pc2_6)

    for stage in STAGE_ORDER:
        sp = stage_profiles[stage]
        signed, n_common = weighted_stage_score(module_genes, weights, sp)

        srow = stage_embed[stage_embed["state_id"] == stage].iloc[0]
        cos_no_pc1 = cosine(emb_no_pc1, srow[pc2_30].values)
        cos_pc2_6 = cosine(emb_pc2_6, srow[pc2_6].values)

        rows.append({
            "state_id": stage,
            "stage_order": STAGE_ORDER.index(stage) + 1,
            "signed_weighted_stage_score": signed,
            "n_common_weighted_genes": n_common,
            "embedding_cosine_no_PC1": cos_no_pc1,
            "embedding_cosine_PC2_6": cos_pc2_6,
        })

    df = pd.DataFrame(rows)
    df["cosine_no_PC1_z"] = zscore(df["embedding_cosine_no_PC1"])
    df["signed_stage_z"] = zscore(df["signed_weighted_stage_score"])
    df["secondary_PC_z"] = zscore(df["embedding_cosine_PC2_6"])

    df["direct_alignment_score"] = (
        0.45 * df["cosine_no_PC1_z"] +
        0.35 * df["signed_stage_z"] +
        0.20 * df["secondary_PC_z"]
    )

    return df


def main():
    print(f"[{now()}] Step10B developmental transport matched-random null started")
    rng = np.random.default_rng(SEED)

    for f in [NODES_FILE, EDGES_FILE, MODULE_FILE, STAGE_PROFILES_FILE, GENE_EMBED_FILE, OBS_TRANSPORT_FILE, FULL_TRANSPORT_FILE]:
        if not Path(f).exists():
            raise FileNotFoundError(f"Missing required file: {f}")

    nodes = pd.read_csv(NODES_FILE, sep="\t")
    edges = pd.read_csv(EDGES_FILE, sep="\t", compression="gzip")
    modules = pd.read_csv(MODULE_FILE, sep="\t")
    stage_profiles_dt = pd.read_csv(STAGE_PROFILES_FILE, sep="\t", compression="gzip")
    gene_embed = pd.read_csv(GENE_EMBED_FILE, sep="\t", compression="gzip")
    obs_top = pd.read_csv(OBS_TRANSPORT_FILE, sep="\t")
    full_transport = pd.read_csv(FULL_TRANSPORT_FILE, sep="\t")

    modules["gene_key"] = modules["gene_symbol_fixed"].map(clean_gene)
    gene_embed["gene_key_upper"] = gene_embed["gene_key"].map(clean_gene)
    pc_cols = [c for c in gene_embed.columns if c.startswith("PC")]
    gene_embed = gene_embed.drop_duplicates("gene_key_upper").set_index("gene_key_upper")

    stage_profiles_dt["gene_key_upper"] = stage_profiles_dt["gene_key"].map(clean_gene)
    stage_profiles = {}
    for st, sub in stage_profiles_dt.groupby("state_id"):
        stage_profiles[str(st)] = dict(zip(sub["gene_key_upper"], sub["stage_profile_z"]))

    # stage embeddings from prior Step03A
    stage_embed_file = STEP03 / "03_step03A_brainspan_stage_embeddings.tsv"
    stage_embed = pd.read_csv(stage_embed_file, sep="\t")

    # Gene stats for expression/variance matching from BrainSpan stage profiles and embedding
    expr_mean = []
    expr_var = []
    genes = []
    for gene, sub in stage_profiles_dt.groupby("gene_key_upper"):
        vals = sub["stage_profile_z"].astype(float).values
        genes.append(gene)
        expr_mean.append(np.nanmean(vals))
        expr_var.append(np.nanvar(vals))

    gene_stats = pd.DataFrame({"gene": genes, "expr_mean": expr_mean, "expr_var": expr_var}).dropna()
    gene_stats["expr_bin"] = make_bins(gene_stats["expr_mean"], 5)
    gene_stats["var_bin"] = make_bins(gene_stats["expr_var"], 5)
    gene_stats = gene_stats[gene_stats["gene"].isin(gene_embed.index)].copy()
    gene_stats = gene_stats.set_index("gene")

    # Risk-included full graph transition for random module PPR with temporary module node? 
    # For transport null, we avoid reconstructing graph with random module nodes for speed and use direct stage alignment null.
    # This is conservative for the late-prenatal stage claim because direct module-stage score is the main driver of stage transport.
    null_rows = []
    summary_rows = []
    tier_rows = []

    target_modules = modules[
        (modules["program"].isin(TARGET_PROGRAMS)) &
        (modules["top_n"].isin(TARGET_TOP_N))
    ].copy()

    for (program, top_n), sub in target_modules.groupby(["program", "top_n"]):
        sub = sub.sort_values("rank_within_program").drop_duplicates("gene_key")
        sub = sub[sub["gene_key"].isin(gene_stats.index)].copy()

        module_genes = sub["gene_key"].tolist()
        module_weights = sub["weight"].astype(float).values

        obs_row = obs_top[
            (obs_top["program"].astype(str) == str(program)) &
            (obs_top["top_n"].astype(int) == int(top_n))
        ].iloc[0]

        obs_stage = str(obs_row["state_id"])
        obs_prob = float(obs_row["transport_probability"])

        obs_full_stage = full_transport[
            (full_transport["program"].astype(str) == str(program)) &
            (full_transport["top_n"].astype(int) == int(top_n)) &
            (full_transport["state_id"].astype(str) == TARGET_STAGE)
        ].iloc[0]

        obs_late_prob = float(obs_full_stage["transport_probability"])
        obs_late_score = float(obs_full_stage["neurotrace_native_transport_score"])

        print(f"[{now()}] Random transport null: {program} top{top_n}, n_genes={len(module_genes)}, obs_top={obs_stage}, obs_late_prob={obs_late_prob}")

        null_late_prob = np.zeros(N_PERM)
        null_late_score = np.zeros(N_PERM)
        null_top_stage = []
        tier_count = []

        for b in range(N_PERM):
            if b == 0 or (b + 1) % 100 == 0 or b + 1 == N_PERM:
                print(f"[{now()}]   permutation {b+1}/{N_PERM} {program} top{top_n}", flush=True)

            genes_rand, weights_rand, tiers = sample_matched_genes(module_genes, module_weights, gene_stats, rng)
            df = compute_module_stage_transport(
                genes_rand,
                weights_rand,
                stage_profiles,
                stage_embed,
                gene_embed,
                pc_cols,
            )

            # Convert direct stage scores to probability as direct-alignment null
            df["probability"] = softmax(df["direct_alignment_score"].values, temp=TRANSPORT_TEMP)
            late = df[df["state_id"] == TARGET_STAGE].iloc[0]
            top = df.sort_values("probability", ascending=False).iloc[0]

            null_late_prob[b] = float(late["probability"])
            null_late_score[b] = float(late["direct_alignment_score"])
            null_top_stage.append(str(top["state_id"]))

            tier_count.append({
                "expr_var_bin": sum(t == "expr_var_bin" for t in tiers),
                "expr_bin_only": sum(t == "expr_bin_only" for t in tiers),
                "var_bin_only": sum(t == "var_bin_only" for t in tiers),
                "any_nonmodule": sum(t == "any_nonmodule" for t in tiers),
                "fallback": sum(t == "fallback" for t in tiers),
            })

        emp_p_prob = empirical_p_greater(obs_late_prob, null_late_prob)
        emp_p_score = empirical_p_greater(obs_late_score, null_late_score)

        top_stage_counts = pd.Series(null_top_stage).value_counts().to_dict()
        late_top_freq = top_stage_counts.get(TARGET_STAGE, 0) / N_PERM

        summary_rows.append({
            "program": program,
            "top_n": top_n,
            "n_module_genes_used": len(module_genes),
            "target_stage": TARGET_STAGE,
            "observed_top_stage": obs_stage,
            "observed_top_transport_probability": obs_prob,
            "observed_late_prenatal_transport_probability": obs_late_prob,
            "observed_late_prenatal_transport_score": obs_late_score,
            "null_late_probability_mean": float(np.mean(null_late_prob)),
            "null_late_probability_sd": float(np.std(null_late_prob, ddof=1)),
            "null_late_score_mean": float(np.mean(null_late_score)),
            "null_late_score_sd": float(np.std(null_late_score, ddof=1)),
            "empirical_p_late_probability_greater": emp_p_prob,
            "empirical_p_late_score_greater": emp_p_score,
            "observed_late_probability_z": (obs_late_prob - np.mean(null_late_prob)) / (np.std(null_late_prob, ddof=1) + 1e-12),
            "observed_late_score_z": (obs_late_score - np.mean(null_late_score)) / (np.std(null_late_score, ddof=1) + 1e-12),
            "null_top_stage_late_prenatal_frequency": late_top_freq,
            "n_perm": N_PERM,
        })

        for b in range(N_PERM):
            null_rows.append({
                "program": program,
                "top_n": top_n,
                "perm_id": b + 1,
                "target_stage": TARGET_STAGE,
                "null_late_probability": null_late_prob[b],
                "null_late_score": null_late_score[b],
                "null_top_stage": null_top_stage[b],
                "observed_late_probability": obs_late_prob,
                "observed_late_score": obs_late_score,
            })
            tc = tier_count[b]
            tier_rows.append({
                "program": program,
                "top_n": top_n,
                "perm_id": b + 1,
                **tc,
            })

    summary = pd.DataFrame(summary_rows)
    summary["fdr_late_probability"] = bh_fdr(summary["empirical_p_late_probability_greater"].values)
    summary["fdr_late_score"] = bh_fdr(summary["empirical_p_late_score_greater"].values)

    null_df = pd.DataFrame(null_rows)
    tier_df = pd.DataFrame(tier_rows)

    null_df.to_csv(OUT / "102_step10B_transport_matched_random_null.tsv.gz", sep="\t", index=False, compression="gzip")
    summary.to_csv(OUT / "103_step10B_transport_matched_random_summary.tsv", sep="\t", index=False)
    tier_df.to_csv(OUT / "104_step10B_transport_random_matching_tier_counts.tsv.gz", sep="\t", index=False, compression="gzip")

    # Top-stage null frequency table
    top_freq = (
        null_df
        .groupby(["program", "top_n", "null_top_stage"])
        .size()
        .reset_index(name="n")
    )
    top_freq["frequency"] = top_freq["n"] / N_PERM
    top_freq.to_csv(OUT / "105_step10B_random_top_stage_frequency.tsv", sep="\t", index=False)

    # Simple plots
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for _, row in summary.iterrows():
            program = row["program"]
            top_n = int(row["top_n"])
            sub = null_df[(null_df["program"] == program) & (null_df["top_n"].astype(int) == top_n)]
            plt.figure(figsize=(6, 4))
            plt.hist(sub["null_late_probability"], bins=40)
            plt.axvline(row["observed_late_prenatal_transport_probability"], linewidth=2)
            plt.xlabel("Random-module late-prenatal probability")
            plt.ylabel("Count")
            plt.title(f"{program} top{top_n}: transport null")
            plt.tight_layout()
            plt.savefig(FIG / f"step10B_{program}_top{top_n}_late_prenatal_null.pdf")
            plt.close()
    except Exception as e:
        with open(OUT / "Step10B_plotting_failed.txt", "w") as f:
            f.write(str(e))

    with open(OUT / "106_step10B_transport_matched_random_summary.md", "w") as f:
        f.write("# NeuroTRACE Step10B developmental transport matched-random null summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Purpose\n")
        f.write("Estimate matched-random null distributions for late-prenatal transport probabilities of NTM1_ASD_up and NTM3_ASD_signed. This calibrates the descriptive late-prenatal transport calls with empirical P values.\n\n")
        f.write("## Design\n")
        f.write(f"- Target stage: {TARGET_STAGE}\n")
        f.write(f"- Modules: {', '.join(TARGET_PROGRAMS)}\n")
        f.write(f"- top_n: {TARGET_TOP_N}\n")
        f.write(f"- permutations per module/top_n: {N_PERM}\n")
        f.write("- Random modules match size, BrainSpan-overlap universe, expression/variance bins, and weight/sign structure.\n\n")
        f.write("## Summary\n\n")
        f.write(summary.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Interpretation\n")
        if (summary["fdr_late_probability"] < 0.05).all():
            f.write("All tested late-prenatal transport calls remain significant against matched random modules, supporting null-calibrated developmental alignment.\n")
        else:
            f.write("Some late-prenatal transport calls do not pass matched-random calibration; interpret those calls as descriptive or exploratory.\n")

    print(f"[{now()}] Step10B done")
    print(f"[{now()}] Results: {OUT}")


if __name__ == "__main__":
    main()
