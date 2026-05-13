#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
STEP = BASE / "neurotrace_algorithm_project" / "step03_feature_embedding"
OUT = STEP / "results"
FIG = STEP / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

PATH_FILE = BASE / "neurotrace_algorithm_project/step01_resource_manifest/results/09_step01B_preferred_paths.sh"

GENE_EMBED = OUT / "01_step03A_brainspan_gene_embedding_pca.tsv.gz"
STAGE_EMBED = OUT / "03_step03A_brainspan_stage_embeddings.tsv"
STAGE_PROFILES = OUT / "04_step03A_brainspan_stage_gene_profiles.tsv.gz"
PROGRAM_EMBED = OUT / "05_step03A_adult_program_embeddings.tsv"
PROGRAM_WEIGHTS = OUT / "06_step03A_adult_program_weight_vectors_matched.tsv"
SFARI_UNIVERSE = OUT / "07_step03A_sfari_prior_on_embedding_universe.tsv"
STEP03A_ALIGN = OUT / "08_step03A_program_to_brainspan_stage_alignment.tsv"

STAGE_ORDER = [
    "early_prenatal",
    "mid_prenatal",
    "late_prenatal",
    "childhood",
    "adolescence",
    "adulthood"
]

PC_SETS = {
    "all_PC1_30": list(range(1, 31)),
    "PC2_30_no_PC1": list(range(2, 31)),
    "PC3_30_no_PC1_PC2": list(range(3, 31)),
    "PC1_only_age_axis": [1],
    "PC2_6_core_secondary": list(range(2, 7))
}


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_path_file(path_file):
    d = {}
    with open(path_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip().strip('"').strip("'")
    return d


def safe_gene_key(x):
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def zscore(x, eps=1e-8):
    x = np.asarray(x, dtype=float)
    return (x - np.nanmean(x)) / (np.nanstd(x) + eps)


def cosine(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        return np.nan
    return float(np.dot(a, b) / ((np.linalg.norm(a) + 1e-12) * (np.linalg.norm(b) + 1e-12)))


def robust_syndromic_parse(x):
    """
    SFARI export 的 syndromic 字段可能不是 1/true/yes。
    这里宽松识别常见编码；同时保留 raw value 供审计。
    """
    if pd.isna(x):
        return 0
    s = str(x).strip().lower()
    if s in ["1", "true", "yes", "y", "syndromic", "syndromic gene", "s"]:
        return 1
    if s in ["0", "false", "no", "n", "", "nan", "none", "not syndromic"]:
        return 0
    # 如果字段是类似 "Syndromic / Functional" 这类长字符串
    if "syndromic" in s and "not" not in s:
        return 1
    return 0


def weighted_stage_score(program_sub, stage_profile_sub):
    merged = program_sub.merge(
        stage_profile_sub[["gene_key", "stage_profile_z"]],
        on="gene_key",
        how="inner"
    )
    if merged.empty:
        return {
            "n_common_weighted_genes": 0,
            "signed_weighted_stage_score": np.nan,
            "abs_weighted_stage_score": np.nan,
            "positive_weight_stage_score": np.nan,
            "negative_weight_stage_score": np.nan
        }

    w = merged["program_weight"].astype(float).values
    p = merged["stage_profile_z"].astype(float).values

    signed = float(np.sum(w * p) / (np.sum(np.abs(w)) + 1e-12))
    abs_score = float(np.sum(np.abs(w) * p) / (np.sum(np.abs(w)) + 1e-12))

    pos = w > 0
    neg = w < 0

    pos_score = float(np.sum(w[pos] * p[pos]) / (np.sum(np.abs(w[pos])) + 1e-12)) if pos.sum() > 0 else np.nan
    neg_score = float(np.sum(np.abs(w[neg]) * p[neg]) / (np.sum(np.abs(w[neg])) + 1e-12)) if neg.sum() > 0 else np.nan

    return {
        "n_common_weighted_genes": int(merged.shape[0]),
        "signed_weighted_stage_score": signed,
        "abs_weighted_stage_score": abs_score,
        "positive_weight_stage_score": pos_score,
        "negative_weight_stage_score": neg_score
    }


def pc_columns_for_set(df, pc_list):
    cols = [f"PC{i}" for i in pc_list if f"PC{i}" in df.columns]
    return cols


def main():
    print(f"[{now()}] Step03B refined alignment diagnostics started")

    for f in [PATH_FILE, GENE_EMBED, STAGE_EMBED, STAGE_PROFILES, PROGRAM_EMBED, PROGRAM_WEIGHTS, SFARI_UNIVERSE, STEP03A_ALIGN]:
        if not Path(f).exists():
            raise FileNotFoundError(f"Missing required file: {f}")

    paths = read_path_file(PATH_FILE)
    sfari_raw_path = paths.get("NEUROTRACE_SFARI_GENE_RELEASE_20260322_EXPORT_20260427", "")

    gene_embed = pd.read_csv(GENE_EMBED, sep="\t")
    stage_embed = pd.read_csv(STAGE_EMBED, sep="\t")
    stage_profiles = pd.read_csv(STAGE_PROFILES, sep="\t")
    program_embed = pd.read_csv(PROGRAM_EMBED, sep="\t")
    program_weights = pd.read_csv(PROGRAM_WEIGHTS, sep="\t")
    sfari_universe = pd.read_csv(SFARI_UNIVERSE, sep="\t")
    step03a_align = pd.read_csv(STEP03A_ALIGN, sep="\t")

    # ---------- 1. SFARI raw syndromic audit ----------
    sfari_audit_rows = []
    sfari_fix_summary_rows = []

    if sfari_raw_path and Path(sfari_raw_path).exists():
        sfari_raw = pd.read_csv(sfari_raw_path)
        syndromic_col = "syndromic" if "syndromic" in sfari_raw.columns else None
        gene_col = "gene-symbol" if "gene-symbol" in sfari_raw.columns else None
        score_col = "gene-score" if "gene-score" in sfari_raw.columns else None

        if syndromic_col:
            tmp = (
                sfari_raw[syndromic_col]
                .astype(str)
                .fillna("NA")
                .value_counts(dropna=False)
                .reset_index()
            )
            tmp.columns = ["syndromic_raw_value", "n"]
            tmp.to_csv(OUT / "12_step03B_sfari_syndromic_raw_value_counts.tsv", sep="\t", index=False)

        if gene_col:
            sfari_raw["gene_key"] = sfari_raw[gene_col].map(safe_gene_key)
            if syndromic_col:
                sfari_raw["syndromic_fixed"] = sfari_raw[syndromic_col].map(robust_syndromic_parse)
            else:
                sfari_raw["syndromic_fixed"] = 0

            if score_col:
                sfari_raw["gene_score_num"] = pd.to_numeric(sfari_raw[score_col], errors="coerce")
            else:
                sfari_raw["gene_score_num"] = np.nan

            universe = sfari_universe[["feature_id", "gene_key"]].drop_duplicates()
            fixed = universe.merge(
                sfari_raw[["gene_key", gene_col, score_col, syndromic_col, "syndromic_fixed", "gene_score_num"]].drop_duplicates("gene_key"),
                on="gene_key",
                how="left"
            )

            fixed["is_SFARI_all_fixed"] = fixed[gene_col].notna().astype(int)
            fixed["is_SFARI_score_le1_fixed"] = ((fixed["gene_score_num"] <= 1) & fixed["gene_score_num"].notna()).astype(int)
            fixed["is_SFARI_score_le2_fixed"] = ((fixed["gene_score_num"] <= 2) & fixed["gene_score_num"].notna()).astype(int)
            fixed["is_SFARI_score_le3_fixed"] = ((fixed["gene_score_num"] <= 3) & fixed["gene_score_num"].notna()).astype(int)
            fixed["is_SFARI_syndromic_fixed"] = fixed["syndromic_fixed"].fillna(0).astype(int)
            fixed["is_SFARI_S_plus_1_fixed"] = ((fixed["is_SFARI_syndromic_fixed"] == 1) | (fixed["is_SFARI_score_le1_fixed"] == 1)).astype(int)
            fixed["is_SFARI_S_plus_1_2_fixed"] = ((fixed["is_SFARI_syndromic_fixed"] == 1) | (fixed["is_SFARI_score_le2_fixed"] == 1)).astype(int)

            fixed.to_csv(OUT / "13_step03B_sfari_prior_on_embedding_universe_fixed.tsv", sep="\t", index=False)

            for col in [
                "is_SFARI_all_fixed",
                "is_SFARI_score_le1_fixed",
                "is_SFARI_score_le2_fixed",
                "is_SFARI_score_le3_fixed",
                "is_SFARI_syndromic_fixed",
                "is_SFARI_S_plus_1_fixed",
                "is_SFARI_S_plus_1_2_fixed"
            ]:
                sfari_fix_summary_rows.append({
                    "sfari_set": col,
                    "n_in_embedding_universe": int(fixed[col].sum())
                })

    if not sfari_fix_summary_rows:
        sfari_fix_summary_rows.append({
            "sfari_set": "sfari_raw_audit_failed_or_missing",
            "n_in_embedding_universe": -1
        })

    pd.DataFrame(sfari_fix_summary_rows).to_csv(
        OUT / "14_step03B_sfari_fixed_set_counts.tsv",
        sep="\t",
        index=False
    )

    # ---------- 2. Recompute refined embedding cosine across PC sets ----------
    stage_ids = stage_embed["state_id"].astype(str).tolist()

    refined_rows = []

    for _, prow in program_embed.iterrows():
        program = prow["program"]
        top_n = prow["top_n"]

        psub = program_weights[
            (program_weights["program"].astype(str) == str(program)) &
            (program_weights["top_n"].astype(str) == str(top_n))
        ].copy()

        for _, srow in stage_embed.iterrows():
            state_id = str(srow["state_id"])
            sprof = stage_profiles[stage_profiles["state_id"].astype(str) == state_id].copy()
            score_dict = weighted_stage_score(psub, sprof)

            base = {
                "program": program,
                "top_n": top_n,
                "state_id": state_id,
                "stage_order": STAGE_ORDER.index(state_id) + 1 if state_id in STAGE_ORDER else 999,
                "n_program_genes": int(psub.shape[0]),
                **score_dict
            }

            for pc_set_name, pc_list in PC_SETS.items():
                cols_p = pc_columns_for_set(program_embed, pc_list)
                cols_s = pc_columns_for_set(stage_embed, pc_list)
                cols = [c for c in cols_p if c in cols_s]
                if len(cols) == 0:
                    val = np.nan
                else:
                    val = cosine(prow[cols].astype(float).values, srow[cols].astype(float).values)
                base[f"embedding_cosine__{pc_set_name}"] = val

            refined_rows.append(base)

    refined = pd.DataFrame(refined_rows)

    # ranks
    for col in [c for c in refined.columns if c.startswith("embedding_cosine__")]:
        refined[f"{col}_rank"] = refined.groupby(["program", "top_n"])[col].rank(
            ascending=False,
            method="min"
        )

    for col in [
        "signed_weighted_stage_score",
        "abs_weighted_stage_score",
        "positive_weight_stage_score",
        "negative_weight_stage_score"
    ]:
        refined[f"{col}_rank"] = refined.groupby(["program", "top_n"])[col].rank(
            ascending=False,
            method="min"
        )

    # composite scores: avoid PC1-dominated allPC as primary
    refined["cosine_no_PC1_z"] = refined.groupby(["program", "top_n"])["embedding_cosine__PC2_30_no_PC1"].transform(lambda x: (x - x.mean()) / (x.std(ddof=0) + 1e-8))
    refined["signed_stage_z"] = refined.groupby(["program", "top_n"])["signed_weighted_stage_score"].transform(lambda x: (x - x.mean()) / (x.std(ddof=0) + 1e-8))
    refined["secondary_PC_z"] = refined.groupby(["program", "top_n"])["embedding_cosine__PC2_6_core_secondary"].transform(lambda x: (x - x.mean()) / (x.std(ddof=0) + 1e-8))

    refined["neurotrace_stage_alignment_score_v1"] = (
        0.45 * refined["cosine_no_PC1_z"] +
        0.35 * refined["signed_stage_z"] +
        0.20 * refined["secondary_PC_z"]
    )

    refined["neurotrace_stage_alignment_rank_v1"] = refined.groupby(["program", "top_n"])["neurotrace_stage_alignment_score_v1"].rank(
        ascending=False,
        method="min"
    )

    refined = refined.sort_values(["program", "top_n", "neurotrace_stage_alignment_rank_v1", "stage_order"])
    refined.to_csv(OUT / "15_step03B_refined_program_to_stage_alignment.tsv", sep="\t", index=False)

    # ---------- 3. Compare all-PC vs no-PC1 top state ----------
    compare_rows = []
    for (program, top_n), sub in refined.groupby(["program", "top_n"]):
        sub = sub.copy()

        def winner(metric):
            rr = sub.sort_values([metric, "stage_order"], ascending=[False, True]).iloc[0]
            return rr["state_id"], rr[metric]

        all_state, all_val = winner("embedding_cosine__all_PC1_30")
        nopc1_state, nopc1_val = winner("embedding_cosine__PC2_30_no_PC1")
        signed_state, signed_val = winner("signed_weighted_stage_score")
        final_state, final_val = winner("neurotrace_stage_alignment_score_v1")

        compare_rows.append({
            "program": program,
            "top_n": top_n,
            "all_PC_winner": all_state,
            "all_PC_cosine": all_val,
            "no_PC1_winner": nopc1_state,
            "no_PC1_cosine": nopc1_val,
            "signed_score_winner": signed_state,
            "signed_score": signed_val,
            "final_refined_winner": final_state,
            "final_refined_score": final_val,
            "allPC_equals_noPC1": int(all_state == nopc1_state),
            "allPC_equals_final": int(all_state == final_state),
            "noPC1_equals_final": int(nopc1_state == final_state)
        })

    compare = pd.DataFrame(compare_rows)
    compare.to_csv(OUT / "16_step03B_alignment_winner_sensitivity.tsv", sep="\t", index=False)

    # ---------- 4. Program-level freeze table ----------
    top1 = refined[refined["neurotrace_stage_alignment_rank_v1"] == 1].copy()
    top1 = top1.sort_values(["program", "top_n"])

    freeze_rows = []
    for _, r in top1.iterrows():
        if r["state_id"] in ["early_prenatal", "mid_prenatal", "late_prenatal"]:
            dev_window = "prenatal"
        elif r["state_id"] in ["childhood", "adolescence"]:
            dev_window = "postnatal_developmental"
        elif r["state_id"] == "adulthood":
            dev_window = "adult_like"
        else:
            dev_window = "unknown"

        freeze_rows.append({
            "program": r["program"],
            "top_n": r["top_n"],
            "refined_top_stage": r["state_id"],
            "developmental_window": dev_window,
            "neurotrace_stage_alignment_score_v1": r["neurotrace_stage_alignment_score_v1"],
            "embedding_cosine_no_PC1": r["embedding_cosine__PC2_30_no_PC1"],
            "embedding_cosine_all_PC": r["embedding_cosine__all_PC1_30"],
            "signed_weighted_stage_score": r["signed_weighted_stage_score"],
            "positive_weight_stage_score": r["positive_weight_stage_score"],
            "negative_weight_stage_score": r["negative_weight_stage_score"],
            "n_common_weighted_genes": r["n_common_weighted_genes"],
            "interpretation_flag": (
                "robust_same_as_allPC" if compare[
                    (compare["program"].astype(str) == str(r["program"])) &
                    (compare["top_n"].astype(str) == str(r["top_n"]))
                ]["allPC_equals_final"].iloc[0] == 1 else "PC1_sensitive_refined"
            )
        })

    freeze = pd.DataFrame(freeze_rows)
    freeze.to_csv(OUT / "17_step03B_program_stage_alignment_freeze.tsv", sep="\t", index=False)

    # ---------- 5. Top3 refined table ----------
    top3 = refined[refined["neurotrace_stage_alignment_rank_v1"] <= 3].copy()
    top3.to_csv(OUT / "18_step03B_top3_refined_stage_alignment_per_program.tsv", sep="\t", index=False)

    # ---------- 6. Plot ----------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        heat = refined.pivot_table(
            index=["program", "top_n"],
            columns="state_id",
            values="neurotrace_stage_alignment_score_v1",
            aggfunc="mean"
        )

        heat = heat[[s for s in STAGE_ORDER if s in heat.columns]]

        plt.figure(figsize=(max(8, 0.65 * heat.shape[1]), max(5, 0.35 * heat.shape[0])))
        plt.imshow(heat.values, aspect="auto")
        plt.xticks(range(heat.shape[1]), heat.columns, rotation=45, ha="right")
        plt.yticks(range(heat.shape[0]), [f"{a}|top{b}" for a, b in heat.index])
        plt.colorbar(label="Refined NeuroTRACE stage alignment score")
        plt.title("Step03B refined program-to-BrainSpan stage alignment")
        plt.tight_layout()
        plt.savefig(FIG / "02_step03B_refined_program_to_stage_heatmap.pdf")
        plt.close()

        # sensitivity counts
        comp_counts = compare[["allPC_equals_noPC1", "allPC_equals_final", "noPC1_equals_final"]].sum()
        plt.figure(figsize=(6, 4))
        plt.bar(comp_counts.index, comp_counts.values)
        plt.ylabel("Number of program/top_n combinations")
        plt.title("Step03B alignment winner sensitivity")
        plt.xticks(rotation=30, ha="right")
        plt.tight_layout()
        plt.savefig(FIG / "03_step03B_alignment_winner_sensitivity.pdf")
        plt.close()

    except Exception as e:
        with open(OUT / "Step03B_plotting_failed.txt", "w") as f:
            f.write(str(e) + "\n")

    # ---------- 7. Summary ----------
    n_pc_sensitive = int((compare["allPC_equals_final"] == 0).sum())
    n_total = int(compare.shape[0])

    with open(OUT / "19_step03B_refined_alignment_summary.md", "w") as f:
        f.write("# NeuroTRACE Step03B refined BrainSpan alignment diagnostics summary\n\n")
        f.write(f"Generated: {now()}\n\n")

        f.write("## Purpose\n")
        f.write("Step03B refines Step03A program-to-stage alignment by diagnosing PC1 sensitivity, recomputing cosine similarity after excluding PC1, combining no-PC1 embedding similarity with signed stage-profile scores, and auditing SFARI syndromic parsing.\n\n")

        f.write("## Key diagnostics\n")
        f.write(f"- Program/top_n combinations tested: {n_total}\n")
        f.write(f"- PC1-sensitive combinations where all-PC winner differs from refined winner: {n_pc_sensitive}\n")
        f.write("- Refined score: 0.45 * no-PC1 embedding z-score + 0.35 * signed stage score z-score + 0.20 * PC2-6 embedding z-score\n\n")

        f.write("## SFARI fixed set counts\n\n")
        f.write(pd.DataFrame(sfari_fix_summary_rows).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Program-stage freeze table\n\n")
        f.write(freeze.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Interpretation\n")
        f.write("If Step03A all-PC alignment is dominated by PC1, Step03B should be preferred for manuscript interpretation. The refined alignment is more appropriate for NeuroTRACE because it reduces the influence of the dominant age axis and integrates signed program expression against developmental stage profiles.\n")

    print(f"[{now()}] Step03B done")
    print(f"[{now()}] Results: {OUT}")
    print(f"[{now()}] Figures: {FIG}")


if __name__ == "__main__":
    main()
