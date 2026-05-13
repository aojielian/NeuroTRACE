#!/usr/bin/env python3

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

GENE_EMBED = OUT / "01_step03A_brainspan_gene_embedding_pca.tsv.gz"
STAGE_EMBED = OUT / "03_step03A_brainspan_stage_embeddings.tsv"
STAGE_PROFILES = OUT / "04_step03A_brainspan_stage_gene_profiles.tsv.gz"
NTM_MODULES = OUT / "24_step03C2_neurotrace_native_module_weights_symbol.tsv"
NTM_DE_SYMBOL = OUT / "23_step03C2_native_DE_symbol_level.tsv.gz"
SFARI_FIXED = OUT / "13_step03B_sfari_prior_on_embedding_universe_fixed.tsv"

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
    "PC2_6_core_secondary": list(range(2, 7)),
    "PC1_only_age_axis": [1],
}

MIN_OVERLAP_FOR_HIGH_CONF = 50
MIN_OVERLAP_RATE_FOR_HIGH_CONF = 0.25


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def zscore(x, eps=1e-8):
    x = np.asarray(x, dtype=float)
    return (x - np.nanmean(x)) / (np.nanstd(x) + eps)


def cosine(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)):
        return np.nan
    return float(np.dot(a, b) / ((np.linalg.norm(a) + 1e-12) * (np.linalg.norm(b) + 1e-12)))


def clean_gene(x):
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def pc_columns_for_set(df, pc_list):
    return [f"PC{i}" for i in pc_list if f"PC{i}" in df.columns]


def weighted_embedding(module_sub, gene_embed_sub, pc_cols):
    m = module_sub.merge(
        gene_embed_sub[["gene_key_upper"] + pc_cols],
        on="gene_key_upper",
        how="inner"
    )
    if m.empty:
        return np.repeat(np.nan, len(pc_cols)), 0

    w = m["weight"].astype(float).values
    W = np.abs(w)
    if W.sum() <= 0:
        return np.repeat(np.nan, len(pc_cols)), m.shape[0]

    X = m[pc_cols].astype(float).values
    emb = (X * (W / W.sum())[:, None]).sum(axis=0)
    return emb, m.shape[0]


def weighted_stage_score(module_sub, stage_profile_sub):
    m = module_sub.merge(
        stage_profile_sub[["gene_key_upper", "stage_profile_z"]],
        on="gene_key_upper",
        how="inner"
    )
    if m.empty:
        return {
            "n_common_weighted_genes": 0,
            "signed_weighted_stage_score": np.nan,
            "abs_weighted_stage_score": np.nan,
            "positive_weight_stage_score": np.nan,
            "negative_weight_stage_score": np.nan,
            "positive_gene_n": 0,
            "negative_gene_n": 0,
        }

    w = m["weight"].astype(float).values
    p = m["stage_profile_z"].astype(float).values

    signed = float(np.sum(w * p) / (np.sum(np.abs(w)) + 1e-12))
    abs_score = float(np.sum(np.abs(w) * p) / (np.sum(np.abs(w)) + 1e-12))

    pos = w > 0
    neg = w < 0

    pos_score = float(np.sum(w[pos] * p[pos]) / (np.sum(np.abs(w[pos])) + 1e-12)) if pos.sum() > 0 else np.nan
    neg_score = float(np.sum(np.abs(w[neg]) * p[neg]) / (np.sum(np.abs(w[neg])) + 1e-12)) if neg.sum() > 0 else np.nan

    return {
        "n_common_weighted_genes": int(m.shape[0]),
        "signed_weighted_stage_score": signed,
        "abs_weighted_stage_score": abs_score,
        "positive_weight_stage_score": pos_score,
        "negative_weight_stage_score": neg_score,
        "positive_gene_n": int(pos.sum()),
        "negative_gene_n": int(neg.sum()),
    }


def hypergeom_sf(k, K, n, N):
    import math
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


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    out = np.empty(n)
    prev = 1.0
    for r in range(n, 0, -1):
        idx = order[r - 1]
        val = p[idx] * n / r
        prev = min(prev, val)
        out[idx] = min(prev, 1.0)
    return out


def main():
    print(f"[{now()}] Step03D native NTM-to-BrainSpan alignment started")

    for f in [GENE_EMBED, STAGE_EMBED, STAGE_PROFILES, NTM_MODULES, NTM_DE_SYMBOL, SFARI_FIXED]:
        if not f.exists():
            raise FileNotFoundError(f"Missing required input: {f}")

    gene_embed = pd.read_csv(GENE_EMBED, sep="\t")
    stage_embed = pd.read_csv(STAGE_EMBED, sep="\t")
    stage_profiles = pd.read_csv(STAGE_PROFILES, sep="\t")
    modules = pd.read_csv(NTM_MODULES, sep="\t")
    de_symbol = pd.read_csv(NTM_DE_SYMBOL, sep="\t")
    sfari = pd.read_csv(SFARI_FIXED, sep="\t")

    gene_embed["gene_key_upper"] = gene_embed["gene_key"].map(clean_gene)
    stage_profiles["gene_key_upper"] = stage_profiles["gene_key"].map(clean_gene)

    if "gene_symbol_fixed" not in modules.columns:
        raise RuntimeError("NTM module file lacks gene_symbol_fixed column.")

    modules["gene_key_upper"] = modules["gene_symbol_fixed"].map(clean_gene)
    modules["module_id"] = modules["program"].astype(str) + "|top" + modules["top_n"].astype(str)

    sfari["gene_key_upper"] = sfari["gene_key"].map(clean_gene)

    pc_cols_all = [c for c in [f"PC{i}" for i in range(1, 31)] if c in gene_embed.columns]
    if len(pc_cols_all) < 5:
        raise RuntimeError("Too few PC columns in BrainSpan gene embedding.")

    # ---------- module overlap summary ----------
    overlap_rows = []
    brain_genes = set(gene_embed["gene_key_upper"])
    sfari_sets = [
        "is_SFARI_all_fixed",
        "is_SFARI_score_le1_fixed",
        "is_SFARI_score_le2_fixed",
        "is_SFARI_syndromic_fixed",
        "is_SFARI_S_plus_1_fixed",
        "is_SFARI_S_plus_1_2_fixed",
    ]
    sfari_sets = [c for c in sfari_sets if c in sfari.columns]

    for (program, top_n), sub in modules.groupby(["program", "top_n"]):
        genes = set(sub["gene_key_upper"])
        ov = genes.intersection(brain_genes)
        row = {
            "program": program,
            "top_n": top_n,
            "n_module_genes": len(genes),
            "n_overlap_brainspan": len(ov),
            "overlap_rate_brainspan": len(ov) / max(len(genes), 1),
        }
        for ss in sfari_sets:
            risk_genes = set(sfari.loc[sfari[ss] == 1, "gene_key_upper"])
            row[f"overlap_{ss}"] = len(genes.intersection(risk_genes))
        overlap_rows.append(row)

    overlap = pd.DataFrame(overlap_rows)
    overlap.to_csv(OUT / "31_step03D_NTM_brainspan_overlap_summary.tsv", sep="\t", index=False)

    # ---------- module embeddings ----------
    module_embed_rows = []
    for (program, top_n), sub in modules.groupby(["program", "top_n"]):
        row = {
            "program": program,
            "top_n": top_n,
            "module_id": f"{program}|top{top_n}",
            "n_module_genes": sub["gene_key_upper"].nunique(),
        }

        for pc_set_name, pc_list in PC_SETS.items():
            pc_cols = pc_columns_for_set(gene_embed, pc_list)
            emb, n_overlap = weighted_embedding(sub, gene_embed, pc_cols)
            row[f"n_overlap_embedding__{pc_set_name}"] = n_overlap
            for i, c in enumerate(pc_cols):
                row[f"{pc_set_name}__{c}"] = emb[i]

        module_embed_rows.append(row)

    module_embed = pd.DataFrame(module_embed_rows)
    module_embed.to_csv(OUT / "32_step03D_NTM_module_embeddings.tsv", sep="\t", index=False)

    # ---------- NTM-to-stage alignment ----------
    align_rows = []

    for _, mrow in module_embed.iterrows():
        program = mrow["program"]
        top_n = mrow["top_n"]
        module_id = mrow["module_id"]

        msub = modules[
            (modules["program"].astype(str) == str(program)) &
            (modules["top_n"].astype(str) == str(top_n))
        ].copy()

        for _, srow in stage_embed.iterrows():
            state_id = str(srow["state_id"])
            sprof = stage_profiles[stage_profiles["state_id"].astype(str) == state_id].copy()

            base = {
                "program": program,
                "top_n": top_n,
                "module_id": module_id,
                "state_id": state_id,
                "stage_order": STAGE_ORDER.index(state_id) + 1 if state_id in STAGE_ORDER else 999,
                "n_module_genes": msub["gene_key_upper"].nunique(),
            }

            score_dict = weighted_stage_score(msub, sprof)
            base.update(score_dict)

            for pc_set_name, pc_list in PC_SETS.items():
                pc_cols = pc_columns_for_set(gene_embed, pc_list)
                mvec = np.array([mrow.get(f"{pc_set_name}__{c}", np.nan) for c in pc_cols], dtype=float)
                svec = srow[pc_cols].astype(float).values
                base[f"embedding_cosine__{pc_set_name}"] = cosine(mvec, svec)
                base[f"n_overlap_embedding__{pc_set_name}"] = int(mrow.get(f"n_overlap_embedding__{pc_set_name}", 0))

            align_rows.append(base)

    align = pd.DataFrame(align_rows)

    for col in [c for c in align.columns if c.startswith("embedding_cosine__")]:
        align[f"{col}_rank"] = align.groupby(["program", "top_n"])[col].rank(ascending=False, method="min")

    for col in [
        "signed_weighted_stage_score",
        "abs_weighted_stage_score",
        "positive_weight_stage_score",
        "negative_weight_stage_score",
    ]:
        align[f"{col}_rank"] = align.groupby(["program", "top_n"])[col].rank(ascending=False, method="min")

    align["cosine_no_PC1_z"] = align.groupby(["program", "top_n"])["embedding_cosine__PC2_30_no_PC1"].transform(zscore)
    align["signed_stage_z"] = align.groupby(["program", "top_n"])["signed_weighted_stage_score"].transform(zscore)
    align["secondary_PC_z"] = align.groupby(["program", "top_n"])["embedding_cosine__PC2_6_core_secondary"].transform(zscore)

    align["neurotrace_native_stage_alignment_score_v1"] = (
        0.45 * align["cosine_no_PC1_z"] +
        0.35 * align["signed_stage_z"] +
        0.20 * align["secondary_PC_z"]
    )

    align["neurotrace_native_stage_alignment_rank_v1"] = align.groupby(["program", "top_n"])["neurotrace_native_stage_alignment_score_v1"].rank(
        ascending=False,
        method="min"
    )

    align["alignment_confidence"] = np.where(
        (align["n_common_weighted_genes"] >= MIN_OVERLAP_FOR_HIGH_CONF) &
        (align["n_common_weighted_genes"] / align["n_module_genes"] >= MIN_OVERLAP_RATE_FOR_HIGH_CONF),
        "usable_overlap",
        "low_overlap_caution"
    )

    align = align.sort_values(["program", "top_n", "neurotrace_native_stage_alignment_rank_v1", "stage_order"])
    align.to_csv(OUT / "33_step03D_NTM_to_brainspan_stage_alignment.tsv", sep="\t", index=False)

    # ---------- top tables ----------
    top1 = align[align["neurotrace_native_stage_alignment_rank_v1"] == 1].copy()
    top1 = top1.sort_values(["program", "top_n"])

    def window_label(stage):
        if stage in ["early_prenatal", "mid_prenatal", "late_prenatal"]:
            return "prenatal"
        if stage in ["childhood", "adolescence"]:
            return "postnatal_developmental"
        if stage == "adulthood":
            return "adult_like"
        return "unknown"

    top1["developmental_window"] = top1["state_id"].map(window_label)

    top1.to_csv(OUT / "34_step03D_NTM_top_stage_freeze.tsv", sep="\t", index=False)

    top3 = align[align["neurotrace_native_stage_alignment_rank_v1"] <= 3].copy()
    top3.to_csv(OUT / "35_step03D_NTM_top3_stage_alignment.tsv", sep="\t", index=False)

    # ---------- all-PC sensitivity ----------
    sens_rows = []
    for (program, top_n), sub in align.groupby(["program", "top_n"]):
        def winner(metric):
            r = sub.sort_values([metric, "stage_order"], ascending=[False, True]).iloc[0]
            return r["state_id"], r[metric]

        all_state, all_val = winner("embedding_cosine__all_PC1_30")
        nopc1_state, nopc1_val = winner("embedding_cosine__PC2_30_no_PC1")
        signed_state, signed_val = winner("signed_weighted_stage_score")
        final_state, final_val = winner("neurotrace_native_stage_alignment_score_v1")

        sens_rows.append({
            "program": program,
            "top_n": top_n,
            "all_PC_winner": all_state,
            "all_PC_cosine": all_val,
            "no_PC1_winner": nopc1_state,
            "no_PC1_cosine": nopc1_val,
            "signed_score_winner": signed_state,
            "signed_score": signed_val,
            "final_native_winner": final_state,
            "final_native_score": final_val,
            "allPC_equals_final": int(all_state == final_state),
            "noPC1_equals_final": int(nopc1_state == final_state),
        })

    sens = pd.DataFrame(sens_rows)
    sens.to_csv(OUT / "36_step03D_NTM_alignment_winner_sensitivity.tsv", sep="\t", index=False)

    # ---------- SFARI enrichment among NTM modules against BrainSpan universe ----------
    risk_rows = []
    universe_genes = set(gene_embed["gene_key_upper"])
    universe_n = len(universe_genes)

    for ss in sfari_sets:
        risk_genes = set(sfari.loc[sfari[ss] == 1, "gene_key_upper"]).intersection(universe_genes)
        K = len(risk_genes)

        for (program, top_n), sub in modules.groupby(["program", "top_n"]):
            genes = set(sub["gene_key_upper"]).intersection(universe_genes)
            n = len(genes)
            k = len(genes.intersection(risk_genes))
            pval = hypergeom_sf(k, K, n, universe_n) if n > 0 and K > 0 else 1.0
            a = k
            b = n - k
            c = K - k
            d = universe_n - K - b
            OR = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
            risk_rows.append({
                "program": program,
                "top_n": top_n,
                "risk_set": ss,
                "overlap_n": k,
                "module_genes_in_brainspan_universe": n,
                "risk_genes_in_brainspan_universe": K,
                "brainspan_universe_n": universe_n,
                "odds_ratio": OR,
                "p_value": pval,
            })

    risk = pd.DataFrame(risk_rows)
    risk["fdr"] = bh_fdr(risk["p_value"].values)
    risk = risk.sort_values(["fdr", "p_value", "program", "top_n"])
    risk.to_csv(OUT / "37_step03D_NTM_SFARI_enrichment_within_brainspan_universe.tsv", sep="\t", index=False)

    # ---------- plot ----------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        heat = align.pivot_table(
            index=["program", "top_n"],
            columns="state_id",
            values="neurotrace_native_stage_alignment_score_v1",
            aggfunc="mean"
        )
        heat = heat[[s for s in STAGE_ORDER if s in heat.columns]]

        plt.figure(figsize=(max(8, 0.75 * heat.shape[1]), max(5, 0.36 * heat.shape[0])))
        plt.imshow(heat.values, aspect="auto")
        plt.xticks(range(heat.shape[1]), heat.columns, rotation=45, ha="right")
        plt.yticks(range(heat.shape[0]), [f"{a}|top{b}" for a, b in heat.index])
        plt.colorbar(label="Native NTM stage alignment score")
        plt.title("Step03D NeuroTRACE-native modules to BrainSpan stages")
        plt.tight_layout()
        plt.savefig(FIG / "05_step03D_NTM_to_BrainSpan_stage_heatmap.pdf")
        plt.close()

        p = top1.copy()
        labels = p["program"].astype(str) + "|top" + p["top_n"].astype(str)
        plt.figure(figsize=(10, 4))
        plt.bar(labels, p["n_common_weighted_genes"])
        plt.xticks(rotation=45, ha="right")
        plt.ylabel("Overlap genes used for alignment")
        plt.title("Step03D NTM-BrainSpan overlap supporting top-stage calls")
        plt.tight_layout()
        plt.savefig(FIG / "06_step03D_NTM_overlap_for_top_stage.pdf")
        plt.close()

    except Exception as e:
        with open(OUT / "Step03D_plotting_failed.txt", "w") as f:
            f.write(str(e) + "\n")

    # ---------- summary ----------
    n_low = int((top1["alignment_confidence"] == "low_overlap_caution").sum())
    n_total = int(top1.shape[0])

    with open(OUT / "38_step03D_NTM_brainspan_alignment_summary.md", "w") as f:
        f.write("# NeuroTRACE Step03D native NTM-to-BrainSpan alignment summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Purpose\n")
        f.write("Step03D projects NeuroTRACE-native Gandal-derived NTM modules to BrainSpan developmental stages. This step does not use DevMap/DevBridge program weights and should replace Prog12/Prog13/Prog14-based prototype analyses in the real-data branch.\n\n")

        f.write("## Inputs\n")
        f.write(f"- NTM symbol modules: `{NTM_MODULES}`\n")
        f.write(f"- BrainSpan gene embedding: `{GENE_EMBED}`\n")
        f.write(f"- BrainSpan stage profiles: `{STAGE_PROFILES}`\n")
        f.write(f"- SFARI fixed prior: `{SFARI_FIXED}`\n\n")

        f.write("## Overlap warning\n")
        f.write(f"- Top-stage calls: {n_total}\n")
        f.write(f"- Low-overlap caution calls: {n_low}\n")
        f.write("Alignment should be interpreted together with `n_common_weighted_genes` and `alignment_confidence`.\n\n")

        f.write("## Top native NTM stage calls\n\n")
        f.write(top1[[
            "program", "top_n", "state_id", "developmental_window",
            "neurotrace_native_stage_alignment_score_v1",
            "n_common_weighted_genes",
            "n_module_genes",
            "alignment_confidence",
            "embedding_cosine__PC2_30_no_PC1",
            "signed_weighted_stage_score"
        ]].to_markdown(index=False))
        f.write("\n\n")

        f.write("## SFARI enrichment within BrainSpan universe: top rows\n\n")
        f.write(risk.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Interpretation\n")
        f.write("Step03D is the first real-data NeuroTRACE-native developmental alignment output. Its strongest claims should be limited to modules and top_n settings with adequate BrainSpan overlap. Low-overlap modules should be kept as exploratory until additional adult cohorts or a larger developmental reference universe are incorporated.\n")

    print(f"[{now()}] Step03D done")
    print(f"[{now()}] NTM modules tested: {modules[['program','top_n']].drop_duplicates().shape[0]}")
    print(f"[{now()}] Top-stage low-overlap caution calls: {n_low}/{n_total}")
    print(f"[{now()}] Results: {OUT}")
    print(f"[{now()}] Figures: {FIG}")


if __name__ == "__main__":
    main()
