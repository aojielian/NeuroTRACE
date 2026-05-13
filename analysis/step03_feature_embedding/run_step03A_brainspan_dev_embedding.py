#!/usr/bin/env python3

import os
import re
import gzip
import math
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
PROJECT = BASE / "neurotrace_algorithm_project"
STEP = PROJECT / "step03_feature_embedding"
OUT = STEP / "results"
FIG = STEP / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

PATH_FILE = BASE / "neurotrace_algorithm_project/step01_resource_manifest/results/09_step01B_preferred_paths.sh"

N_PCS = 30
MIN_STAGE_SAMPLES = 5
TOP_GENES_FOR_STAGE_EMBED = 500
SEED = 20260508

def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def read_path_file(path_file):
    d = {}
    with open(path_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            v = v.strip().strip('"').strip("'")
            d[k.strip()] = v
    return d

def zscore_rows(mat, eps=1e-8):
    mat = np.asarray(mat, dtype=float)
    mu = np.nanmean(mat, axis=1, keepdims=True)
    sd = np.nanstd(mat, axis=1, keepdims=True)
    return (mat - mu) / (sd + eps)

def zscore_vec(x, eps=1e-8):
    x = np.asarray(x, dtype=float)
    return (x - np.nanmean(x)) / (np.nanstd(x) + eps)

def safe_gene_key(x):
    if pd.isna(x):
        return ""
    return str(x).strip().upper()

def weighted_embedding(gene_scores, gene_embedding, top_n=500):
    scores = np.asarray(gene_scores, dtype=float)
    if np.all(~np.isfinite(scores)):
        return np.repeat(np.nan, gene_embedding.shape[1])
    scores = np.nan_to_num(scores, nan=0.0)
    idx = np.argsort(-np.abs(scores))[:min(top_n, len(scores))]
    w = np.abs(scores[idx])
    if w.sum() <= 0:
        return np.repeat(np.nan, gene_embedding.shape[1])
    w = w / w.sum()
    return (gene_embedding[idx, :] * w[:, None]).sum(axis=0)

def cosine(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        return np.nan
    return float(np.dot(a, b) / ((np.linalg.norm(a) + 1e-12) * (np.linalg.norm(b) + 1e-12)))

def signed_weighted_stage_score(weights, stage_profile):
    w = np.asarray(weights, dtype=float)
    p = np.asarray(stage_profile, dtype=float)
    denom = np.sum(np.abs(w)) + 1e-12
    return float(np.sum(w * p) / denom)

def main():
    print(f"[{now()}] Step03A BrainSpan developmental embedding started")

    paths = read_path_file(PATH_FILE)

    required = [
        "NEUROTRACE_BRAINSPAN_RRSD_EXPRESSION",
        "NEUROTRACE_BRAINSPAN_RRSD_SAMPLES",
        "NEUROTRACE_BRAINSPAN_RRSD_GENE_INFO",
        "NEUROTRACE_DEVBRIDGE_SIGNATURE_WEIGHTS_V3",
        "NEUROTRACE_SFARI_GENE_RELEASE_20260322_EXPORT_20260427"
    ]

    for k in required:
        if k not in paths:
            raise RuntimeError(f"Missing path variable in {PATH_FILE}: {k}")
        if not Path(paths[k]).exists():
            raise FileNotFoundError(f"{k} path does not exist: {paths[k]}")

    expr_path = paths["NEUROTRACE_BRAINSPAN_RRSD_EXPRESSION"]
    sample_path = paths["NEUROTRACE_BRAINSPAN_RRSD_SAMPLES"]
    gene_info_path = paths["NEUROTRACE_BRAINSPAN_RRSD_GENE_INFO"]
    weight_path = paths["NEUROTRACE_DEVBRIDGE_SIGNATURE_WEIGHTS_V3"]
    sfari_path = paths["NEUROTRACE_SFARI_GENE_RELEASE_20260322_EXPORT_20260427"]

    print(f"[{now()}] Reading BrainSpan expression: {expr_path}")
    expr = pd.read_csv(expr_path, sep="\t")
    if "feature_id" not in expr.columns:
        raise RuntimeError("BrainSpan expression must contain feature_id column")

    print(f"[{now()}] Reading BrainSpan samples: {sample_path}")
    samples = pd.read_csv(sample_path, sep="\t")

    print(f"[{now()}] Reading BrainSpan gene info: {gene_info_path}")
    gene_info = pd.read_csv(gene_info_path, sep="\t")

    print(f"[{now()}] Reading DevBridge signature weights: {weight_path}")
    weights = pd.read_csv(weight_path, sep="\t")

    print(f"[{now()}] Reading SFARI: {sfari_path}")
    sfari = pd.read_csv(sfari_path)

    # Expression cleanup
    expr["gene_key"] = expr["feature_id"].map(safe_gene_key)
    expr = expr[expr["gene_key"] != ""].copy()
    expr = expr.drop_duplicates("gene_key", keep="first")

    sample_cols = [c for c in expr.columns if c.startswith("BS_")]
    if len(sample_cols) == 0:
        sample_cols = [c for c in expr.columns if c != "feature_id" and c != "gene_key"]

    sample_ids = set(samples["sample_id"].astype(str)) if "sample_id" in samples.columns else set()
    sample_cols = [c for c in sample_cols if c in sample_ids]

    if len(sample_cols) < 20:
        raise RuntimeError(f"Too few matched BrainSpan sample columns: {len(sample_cols)}")

    samples = samples[samples["sample_id"].astype(str).isin(sample_cols)].copy()
    samples["sample_id"] = samples["sample_id"].astype(str)

    # Order sample metadata to expression columns
    samples = samples.set_index("sample_id").loc[sample_cols].reset_index()

    stage_col = "stage_group" if "stage_group" in samples.columns else None
    if stage_col is None:
        raise RuntimeError("BrainSpan sample metadata lacks stage_group")

    X = expr[sample_cols].apply(pd.to_numeric, errors="coerce").values
    finite_gene = np.isfinite(X).mean(axis=1) >= 0.95
    X = X[finite_gene, :]
    genes = expr.loc[finite_gene, "feature_id"].astype(str).values
    gene_keys = expr.loc[finite_gene, "gene_key"].astype(str).values

    row_var = np.nanvar(X, axis=1)
    keep_var = row_var > 1e-8
    X = X[keep_var, :]
    genes = genes[keep_var]
    gene_keys = gene_keys[keep_var]

    print(f"[{now()}] Expression matrix after filtering: genes={X.shape[0]}, samples={X.shape[1]}")

    Xz = zscore_rows(X)
    Xz = np.nan_to_num(Xz, nan=0.0)

    # SVD gene embedding
    print(f"[{now()}] Computing SVD gene embedding")
    U, S, Vt = np.linalg.svd(Xz, full_matrices=False)
    n_pc = min(N_PCS, U.shape[1])
    gene_embedding = U[:, :n_pc] * S[:n_pc]

    pc_cols = [f"PC{i+1}" for i in range(n_pc)]

    gene_embed_df = pd.DataFrame(gene_embedding, columns=pc_cols)
    gene_embed_df.insert(0, "gene_key", gene_keys)
    gene_embed_df.insert(0, "feature_id", genes)
    gene_embed_df.to_csv(OUT / "01_step03A_brainspan_gene_embedding_pca.tsv.gz", sep="\t", index=False, compression="gzip")

    # Variance explained
    var = S**2
    ve = var / var.sum()
    ve_df = pd.DataFrame({
        "PC": [f"PC{i+1}" for i in range(len(ve))],
        "variance_explained": ve,
        "cumulative_variance_explained": np.cumsum(ve)
    })
    ve_df.to_csv(OUT / "02_step03A_brainspan_pca_variance_explained.tsv", sep="\t", index=False)

    # Developmental stage profiles and embeddings
    print(f"[{now()}] Computing stage profiles and embeddings")
    stage_rows = []
    stage_profile_rows = []

    sample_stage = samples[["sample_id", stage_col]].copy()
    sample_stage[stage_col] = sample_stage[stage_col].astype(str)

    for stage, ss in sample_stage.groupby(stage_col):
        cols = ss["sample_id"].tolist()
        idx = [sample_cols.index(c) for c in cols if c in sample_cols]
        if len(idx) < MIN_STAGE_SAMPLES:
            continue

        profile = Xz[:, idx].mean(axis=1)
        emb = weighted_embedding(profile, gene_embedding, top_n=TOP_GENES_FOR_STAGE_EMBED)

        row = {
            "state_id": stage,
            "state_type": "BrainSpan_stage_group",
            "n_samples": len(idx)
        }
        for i, val in enumerate(emb):
            row[f"PC{i+1}"] = val
        stage_rows.append(row)

        prof = pd.DataFrame({
            "state_id": stage,
            "feature_id": genes,
            "gene_key": gene_keys,
            "stage_profile_z": profile
        })
        stage_profile_rows.append(prof)

    stage_embed = pd.DataFrame(stage_rows)
    stage_embed.to_csv(OUT / "03_step03A_brainspan_stage_embeddings.tsv", sep="\t", index=False)

    stage_profiles = pd.concat(stage_profile_rows, ignore_index=True)
    stage_profiles.to_csv(OUT / "04_step03A_brainspan_stage_gene_profiles.tsv.gz", sep="\t", index=False, compression="gzip")

    # DevBridge adult program embeddings
    print(f"[{now()}] Computing adult program embeddings")
    weights["gene_key"] = weights["gene_symbol"].map(safe_gene_key)

    gene_index = {g: i for i, g in enumerate(gene_keys)}
    program_rows = []
    program_vector_rows = []

    group_cols = ["program", "top_n"] if "top_n" in weights.columns else ["program"]

    for key, sub in weights.groupby(group_cols):
        if isinstance(key, tuple):
            program, top_n = key
        else:
            program, top_n = key, ""

        vec = np.zeros(len(gene_keys), dtype=float)
        matched = 0
        for _, r in sub.iterrows():
            g = r["gene_key"]
            if g in gene_index:
                vec[gene_index[g]] = float(r["weight"])
                matched += 1

        if matched == 0:
            continue

        emb = weighted_embedding(vec, gene_embedding, top_n=min(TOP_GENES_FOR_STAGE_EMBED, max(50, matched)))
        row = {
            "program": program,
            "top_n": top_n,
            "n_weight_genes": len(sub),
            "n_matched_universe": matched,
            "match_rate": matched / max(len(sub), 1)
        }
        for i, val in enumerate(emb):
            row[f"PC{i+1}"] = val
        program_rows.append(row)

        tmp = pd.DataFrame({
            "program": program,
            "top_n": top_n,
            "feature_id": genes,
            "gene_key": gene_keys,
            "program_weight": vec
        })
        tmp = tmp[tmp["program_weight"] != 0].copy()
        program_vector_rows.append(tmp)

    program_embed = pd.DataFrame(program_rows)
    program_embed.to_csv(OUT / "05_step03A_adult_program_embeddings.tsv", sep="\t", index=False)

    if program_vector_rows:
        program_vectors = pd.concat(program_vector_rows, ignore_index=True)
    else:
        program_vectors = pd.DataFrame(columns=["program", "top_n", "feature_id", "gene_key", "program_weight"])
    program_vectors.to_csv(OUT / "06_step03A_adult_program_weight_vectors_matched.tsv", sep="\t", index=False)

    # SFARI risk prior
    print(f"[{now()}] Mapping SFARI to embedding universe")
    sfari_gene_col = "gene-symbol" if "gene-symbol" in sfari.columns else None
    if sfari_gene_col is None:
        raise RuntimeError("SFARI file lacks gene-symbol column")

    sfari["gene_key"] = sfari[sfari_gene_col].map(safe_gene_key)

    sfari_score_col = "gene-score" if "gene-score" in sfari.columns else None
    syndromic_col = "syndromic" if "syndromic" in sfari.columns else None

    universe = pd.DataFrame({
        "feature_id": genes,
        "gene_key": gene_keys
    })

    sfari_small = sfari.copy()
    keep_cols = ["gene_key", sfari_gene_col]
    if sfari_score_col:
        keep_cols.append(sfari_score_col)
    if syndromic_col:
        keep_cols.append(syndromic_col)
    sfari_small = sfari_small[keep_cols].drop_duplicates("gene_key")

    risk = universe.merge(sfari_small, on="gene_key", how="left")
    risk["is_SFARI_all"] = risk[sfari_gene_col].notna().astype(int)

    if sfari_score_col:
        score_num = pd.to_numeric(risk[sfari_score_col], errors="coerce")
        risk["is_SFARI_score_le1"] = ((score_num <= 1) & score_num.notna()).astype(int)
        risk["is_SFARI_score_le2"] = ((score_num <= 2) & score_num.notna()).astype(int)
        risk["is_SFARI_score_le3"] = ((score_num <= 3) & score_num.notna()).astype(int)
    else:
        risk["is_SFARI_score_le1"] = 0
        risk["is_SFARI_score_le2"] = 0
        risk["is_SFARI_score_le3"] = 0

    if syndromic_col:
        syn = risk[syndromic_col].astype(str).str.lower()
        risk["is_SFARI_syndromic"] = syn.isin(["1", "true", "yes", "y"]).astype(int)
    else:
        risk["is_SFARI_syndromic"] = 0

    risk["is_SFARI_S_plus_1"] = ((risk["is_SFARI_syndromic"] == 1) | (risk["is_SFARI_score_le1"] == 1)).astype(int)
    risk["is_SFARI_S_plus_1_2"] = ((risk["is_SFARI_syndromic"] == 1) | (risk["is_SFARI_score_le2"] == 1)).astype(int)

    risk.to_csv(OUT / "07_step03A_sfari_prior_on_embedding_universe.tsv", sep="\t", index=False)

    # Program-to-stage preliminary alignment
    print(f"[{now()}] Computing preliminary program-to-stage alignment")
    stage_embed_mat = stage_embed[pc_cols].values
    prog_embed_mat = program_embed[pc_cols].values

    stage_profile_map = {}
    for stage, sub in stage_profiles.groupby("state_id"):
        stage_profile_map[stage] = sub.set_index("gene_key")["stage_profile_z"].to_dict()

    align_rows = []

    for i, prow in program_embed.iterrows():
        pvec_embed = prog_embed_mat[i, :]
        program = prow["program"]
        top_n = prow["top_n"]

        pweights = program_vectors[
            (program_vectors["program"] == program) &
            (program_vectors["top_n"].astype(str) == str(top_n))
        ][["gene_key", "program_weight"]].copy()

        wmap = dict(zip(pweights["gene_key"], pweights["program_weight"]))

        for j, srow in stage_embed.iterrows():
            state = srow["state_id"]
            svec_embed = stage_embed_mat[j, :]
            cos = cosine(pvec_embed, svec_embed)

            common = [g for g in wmap.keys() if g in stage_profile_map[state]]
            if common:
                w = np.array([wmap[g] for g in common], dtype=float)
                prof = np.array([stage_profile_map[state][g] for g in common], dtype=float)
                signed_score = signed_weighted_stage_score(w, prof)
                abs_score = float(np.sum(np.abs(w) * prof) / (np.sum(np.abs(w)) + 1e-12))
            else:
                signed_score = np.nan
                abs_score = np.nan

            align_rows.append({
                "program": program,
                "top_n": top_n,
                "state_id": state,
                "state_type": "BrainSpan_stage_group",
                "n_program_genes_matched": int(prow["n_matched_universe"]),
                "n_common_weighted_genes": len(common),
                "embedding_cosine": cos,
                "signed_weighted_stage_score": signed_score,
                "abs_weighted_stage_score": abs_score
            })

    align = pd.DataFrame(align_rows)
    align["embedding_cosine_rank_within_program"] = align.groupby(["program", "top_n"])["embedding_cosine"].rank(
        ascending=False,
        method="min"
    )
    align["signed_score_rank_within_program"] = align.groupby(["program", "top_n"])["signed_weighted_stage_score"].rank(
        ascending=False,
        method="min"
    )
    align.to_csv(OUT / "08_step03A_program_to_brainspan_stage_alignment.tsv", sep="\t", index=False)

    # Summary tables
    overlap_rows = [
        {"item": "BrainSpan_embedding_genes", "n": len(gene_keys)},
        {"item": "BrainSpan_samples", "n": len(sample_cols)},
        {"item": "BrainSpan_stage_groups_used", "n": stage_embed.shape[0]},
        {"item": "Adult_program_embeddings", "n": program_embed.shape[0]},
        {"item": "Adult_program_weighted_gene_rows_matched", "n": program_vectors.shape[0]},
        {"item": "SFARI_all_in_universe", "n": int(risk["is_SFARI_all"].sum())},
        {"item": "SFARI_score_le1_in_universe", "n": int(risk["is_SFARI_score_le1"].sum())},
        {"item": "SFARI_score_le2_in_universe", "n": int(risk["is_SFARI_score_le2"].sum())},
        {"item": "SFARI_S_plus_1_in_universe", "n": int(risk["is_SFARI_S_plus_1"].sum())},
        {"item": "SFARI_S_plus_1_2_in_universe", "n": int(risk["is_SFARI_S_plus_1_2"].sum())}
    ]
    pd.DataFrame(overlap_rows).to_csv(OUT / "09_step03A_universe_overlap_summary.tsv", sep="\t", index=False)

    top_align = align.sort_values(
        ["program", "top_n", "embedding_cosine_rank_within_program", "signed_score_rank_within_program"]
    ).groupby(["program", "top_n"]).head(5)
    top_align.to_csv(OUT / "10_step03A_top5_stage_alignment_per_program.tsv", sep="\t", index=False)

    # Optional heatmap figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        heat = align.pivot_table(
            index="program",
            columns="state_id",
            values="embedding_cosine",
            aggfunc="mean"
        )

        plt.figure(figsize=(max(8, 0.55 * heat.shape[1]), max(5, 0.35 * heat.shape[0])))
        plt.imshow(heat.values, aspect="auto")
        plt.xticks(range(heat.shape[1]), heat.columns, rotation=45, ha="right")
        plt.yticks(range(heat.shape[0]), heat.index)
        plt.colorbar(label="Embedding cosine")
        plt.title("Step03A preliminary adult program-to-BrainSpan stage alignment")
        plt.tight_layout()
        plt.savefig(FIG / "01_step03A_program_to_stage_embedding_cosine_heatmap.pdf")
        plt.close()
    except Exception as e:
        with open(OUT / "Step03A_plotting_failed.txt", "w") as f:
            f.write(str(e) + "\n")

    with open(OUT / "11_step03A_brainspan_embedding_summary.md", "w") as f:
        f.write("# NeuroTRACE Step03A BrainSpan developmental embedding summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Inputs\n")
        f.write(f"- BrainSpan RRSD expression: `{expr_path}`\n")
        f.write(f"- BrainSpan sample metadata: `{sample_path}`\n")
        f.write(f"- BrainSpan gene info: `{gene_info_path}`\n")
        f.write(f"- DevBridge signature weights: `{weight_path}`\n")
        f.write(f"- SFARI gene prior: `{sfari_path}`\n\n")
        f.write("## Output summary\n")
        for r in overlap_rows:
            f.write(f"- {r['item']}: {r['n']}\n")
        f.write("\n## Main outputs\n")
        f.write("- `01_step03A_brainspan_gene_embedding_pca.tsv.gz`\n")
        f.write("- `03_step03A_brainspan_stage_embeddings.tsv`\n")
        f.write("- `04_step03A_brainspan_stage_gene_profiles.tsv.gz`\n")
        f.write("- `05_step03A_adult_program_embeddings.tsv`\n")
        f.write("- `07_step03A_sfari_prior_on_embedding_universe.tsv`\n")
        f.write("- `08_step03A_program_to_brainspan_stage_alignment.tsv`\n")
        f.write("- `10_step03A_top5_stage_alignment_per_program.tsv`\n\n")
        f.write("## Interpretation\n")
        f.write("Step03A creates a real-data PCA fallback embedding from BrainSpan cortical developmental expression. This is the first real-data representation layer for NeuroTRACE and can be replaced or augmented later by foundation-model gene embeddings.\n")

    print(f"[{now()}] Step03A done")
    print(f"[{now()}] Results: {OUT}")
    print(f"[{now()}] Figures: {FIG}")

if __name__ == "__main__":
    main()
