#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step11G: BrainSpan universe coverage and nested-universe sensitivity.

Purpose
-------
Address the concern that the current NeuroTRACE developmental graph is restricted
to a 5,000-gene BrainSpan embedding universe. This script quantifies:

1. Coverage of NTM module genes across nested BrainSpan universes.
2. Coverage of SFARI, Satterstrom, and available Zhou gene-set resources.
3. Whether direct developmental stage alignment for NTM1/NTM3 and NTM2
   is stable across nested universes within the available 5k BrainSpan graph.
4. Whether larger BrainSpan raw files appear to be available for future 8k/10k
   reconstruction.

Important interpretation
------------------------
This does not fully replace a true 10k BrainSpan graph rebuild unless original
BrainSpan expression inputs are available. It provides a coverage-aware and
nested-universe sensitivity analysis within the current 5k embedding universe.
"""

import argparse
import os
import re
import glob
import math
import numpy as np
import pandas as pd


def log(msg):
    print(f"[Step11G] {msg}", flush=True)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def read_table(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    if path.endswith(".csv"):
        return pd.read_csv(path, low_memory=False)
    if path.endswith(".xlsx"):
        return pd.read_excel(path)
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
    if "|" in x and x.lower().startswith("gene|"):
        x = x.split("|")[-1]
    return x.strip()


def clean_stage(x):
    return str(x).strip().replace("stage:", "")


def zscore(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    if not np.isfinite(sd) or sd == 0:
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


def infer_module_family(x):
    x = str(x)
    if "NTM1" in x:
        return "NTM1_ASD_up"
    if "NTM2" in x:
        return "NTM2_ASD_down"
    if "NTM3" in x:
        return "NTM3_ASD_signed"
    return x


def weighted_stage_projection(module_weights, stage_profile_df):
    """
    module_weights: DataFrame with gene, weight
    stage_profile_df: long table with gene, stage, stage_profile_z
    """
    merged = stage_profile_df.merge(module_weights, on="gene", how="inner")
    if merged.shape[0] == 0:
        return pd.DataFrame()

    records = []
    for stage, sub in merged.groupby("stage_clean"):
        w = pd.to_numeric(sub["weight"], errors="coerce").values.astype(float)
        s = pd.to_numeric(sub["stage_profile_z"], errors="coerce").values.astype(float)

        ok = np.isfinite(w) & np.isfinite(s)
        w = w[ok]
        s = s[ok]

        if len(w) == 0:
            raw = np.nan
            n = 0
        else:
            denom = math.sqrt(np.sum(w * w)) * math.sqrt(np.sum(s * s))
            raw = float(np.sum(w * s) / denom) if denom > 0 else 0.0
            n = len(w)

        records.append({
            "stage_clean": stage,
            "raw_alignment_score": raw,
            "n_common_weighted_genes": n
        })

    out = pd.DataFrame(records)
    out["z_alignment_score"] = zscore(out["raw_alignment_score"].values)
    out["stage_probability"] = softmax(out["z_alignment_score"].values, temperature=0.75)
    return out


def load_gene_set_file(path):
    """
    Robustly read gene sets from txt/tsv/csv/gz. Returns dict name -> set(genes).
    """
    out = {}
    base = os.path.basename(path)

    try:
        if path.endswith(".gz"):
            df = pd.read_csv(path, sep=None, engine="python", compression="gzip", low_memory=False)
        elif path.endswith(".csv"):
            df = pd.read_csv(path, low_memory=False)
        else:
            # Try tabular first.
            df = pd.read_csv(path, sep=None, engine="python", low_memory=False)

        # If table has obvious gene column.
        gene_col = find_col(df, ["gene", "gene_symbol", "symbol", "Gene", "genes", "gene_name"])
        set_col = find_col(df, ["gene_set", "set", "pathway", "name", "category", "risk_set"])

        if gene_col is not None:
            if set_col is not None:
                for set_name, sub in df.groupby(set_col):
                    out[f"{base}:{set_name}"] = set(clean_gene(x) for x in sub[gene_col].dropna().astype(str))
            else:
                out[base] = set(clean_gene(x) for x in df[gene_col].dropna().astype(str))
            return out

        # Otherwise flatten all string cells as possible genes for small files.
        if df.shape[1] == 1:
            out[base] = set(clean_gene(x) for x in df.iloc[:, 0].dropna().astype(str))
            return out

    except Exception:
        pass

    # Plain text fallback.
    try:
        genes = []
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = re.split(r"[,\t; ]+", line)
                for p in parts:
                    p = clean_gene(p)
                    if p and re.match(r"^[A-Za-z0-9._-]+$", p):
                        genes.append(p)
        if genes:
            out[base] = set(genes)
    except Exception:
        pass

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--embedding",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step03_feature_embedding/results/01_step03A_brainspan_gene_embedding_pca.tsv.gz",
    )
    parser.add_argument(
        "--stage_profiles",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step03_feature_embedding/results/04_step03A_brainspan_stage_gene_profiles.tsv.gz",
    )
    parser.add_argument(
        "--module_weights",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step03_feature_embedding/results/24_step03C2_neurotrace_native_module_weights_symbol.tsv",
    )
    parser.add_argument(
        "--sfari_fixed",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step03_feature_embedding/results/13_step03B_sfari_prior_on_embedding_universe_fixed.tsv",
    )
    parser.add_argument(
        "--satterstrom",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/Satterstrom_102_ASD_risk_genes.txt",
    )
    parser.add_argument(
        "--zhou_gene_set_dir",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step10_zhou2022_genesets_fixed/gene_sets",
    )
    parser.add_argument(
        "--project_root",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project",
    )
    parser.add_argument(
        "--outdir",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step11_robustness_sensitivity/results",
    )
    args = parser.parse_args()

    ensure_dir(args.outdir)

    log("Reading BrainSpan embedding")
    emb = read_table(args.embedding)
    gene_col_emb = find_col(emb, ["gene_key", "gene", "feature_id", "symbol"], required=True)
    emb["gene"] = emb[gene_col_emb].map(clean_gene)

    log("Reading BrainSpan stage profiles")
    sp = read_table(args.stage_profiles)
    gene_col_sp = find_col(sp, ["gene_key", "gene", "feature_id", "symbol"], required=True)
    stage_col = find_col(sp, ["state_id", "stage", "stage_id"], required=True)
    val_col = find_col(sp, ["stage_profile_z", "profile_z", "z", "value"], required=True)

    sp = sp[[gene_col_sp, stage_col, val_col]].copy()
    sp.columns = ["gene", "stage", "stage_profile_z"]
    sp["gene"] = sp["gene"].map(clean_gene)
    sp["stage_clean"] = sp["stage"].map(clean_stage)
    sp["stage_profile_z"] = pd.to_numeric(sp["stage_profile_z"], errors="coerce")

    # Rank current 5k universe by variance of stage profiles.
    gene_var = (
        sp.groupby("gene", dropna=False)
        .agg(stage_profile_var=("stage_profile_z", "var"),
             stage_profile_absmean=("stage_profile_z", lambda x: np.nanmean(np.abs(pd.to_numeric(x, errors="coerce")))))
        .reset_index()
    )
    gene_var["stage_profile_var"] = pd.to_numeric(gene_var["stage_profile_var"], errors="coerce").fillna(0)
    gene_var["stage_profile_absmean"] = pd.to_numeric(gene_var["stage_profile_absmean"], errors="coerce").fillna(0)

    # Restrict to genes in embedding.
    emb_genes = set(emb["gene"].dropna().astype(str))
    gene_var = gene_var[gene_var["gene"].isin(emb_genes)].copy()
    gene_var = gene_var.sort_values(["stage_profile_var", "stage_profile_absmean"], ascending=False)
    gene_var["rank_within_5k"] = np.arange(1, gene_var.shape[0] + 1)

    available_n = gene_var.shape[0]
    requested_sizes = [1000, 2000, 3000, 4000, 5000]
    universe_sizes = [n for n in requested_sizes if n <= available_n]
    if available_n not in universe_sizes:
        universe_sizes.append(available_n)
    universe_sizes = sorted(set(universe_sizes))

    log(f"Available BrainSpan embedding genes: {available_n}")
    log(f"Nested universe sizes: {universe_sizes}")

    # Load module weights.
    log("Reading module weights")
    mw = read_table(args.module_weights)
    gene_col_mw = find_col(mw, ["gene_symbol_fixed", "gene_symbol", "gene", "symbol"], required=True)
    program_col = find_col(mw, ["program", "module_family", "module"], required=True)
    topn_col = find_col(mw, ["top_n", "topn"], required=True)
    weight_col = find_col(mw, ["weight", "program_weight", "module_weight", "signed_weight"], required=True)

    mw2 = pd.DataFrame()
    mw2["gene"] = mw[gene_col_mw].map(clean_gene)
    mw2["program"] = mw[program_col].astype(str)
    mw2["module_family"] = mw2["program"].map(infer_module_family)
    mw2["top_n"] = pd.to_numeric(mw[topn_col], errors="coerce")
    mw2["weight"] = pd.to_numeric(mw[weight_col], errors="coerce")

    mw2 = mw2[
        mw2["module_family"].isin(["NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed"]) &
        mw2["top_n"].isin([50, 100, 200, 500]) &
        np.isfinite(mw2["weight"])
    ].copy()

    # Load SFARI fixed as risk sets.
    risk_sets = {}
    if os.path.exists(args.sfari_fixed):
        sf = read_table(args.sfari_fixed)
        sf_gene_col = find_col(sf, ["gene_key", "gene", "gene-symbol", "symbol", "feature_id"])
        if sf_gene_col is not None:
            sf["gene"] = sf[sf_gene_col].map(clean_gene)
            for c in sf.columns:
                if c.startswith("is_SFARI"):
                    risk_sets[c] = set(sf.loc[pd.to_numeric(sf[c], errors="coerce").fillna(0).astype(int).eq(1), "gene"].dropna().astype(str))
                    if len(risk_sets[c]) == 0:
                        risk_sets[c] = set(sf.loc[sf[c].astype(str).str.upper().isin(["TRUE", "T", "1"]), "gene"].dropna().astype(str))

    # Satterstrom.
    if os.path.exists(args.satterstrom):
        gs = load_gene_set_file(args.satterstrom)
        for k, v in gs.items():
            risk_sets[f"Satterstrom102:{k}"] = set(v)

    # Zhou gene sets.
    if os.path.exists(args.zhou_gene_set_dir):
        for path in sorted(glob.glob(os.path.join(args.zhou_gene_set_dir, "*"))):
            if os.path.isfile(path) and re.search(r"\.(txt|tsv|csv|gz)$", path):
                gs = load_gene_set_file(path)
                for k, v in gs.items():
                    risk_sets[f"Zhou2022:{k}"] = set(v)

    # Keep only non-empty sets.
    risk_sets = {k: set(clean_gene(g) for g in v if clean_gene(g)) for k, v in risk_sets.items() if len(v) > 0}
    log(f"Risk/gene sets loaded: {len(risk_sets)}")

    # Larger BrainSpan file audit.
    candidate_records = []
    patterns = [
        "**/*BrainSpan*",
        "**/*brainspan*",
        "**/*stage*profile*",
        "**/*expression*",
    ]
    seen = set()
    for pat in patterns:
        for path in glob.glob(os.path.join(args.project_root, pat), recursive=True):
            if path in seen or not os.path.isfile(path):
                continue
            seen.add(path)
            try:
                size_mb = os.path.getsize(path) / (1024 * 1024)
                base = os.path.basename(path)
                if size_mb < 0.001:
                    continue
                candidate_records.append({
                    "path": path,
                    "file_name": base,
                    "size_mb": size_mb,
                    "mtime": os.path.getmtime(path),
                })
            except Exception:
                pass

    candidate_df = pd.DataFrame(candidate_records).sort_values("size_mb", ascending=False) if candidate_records else pd.DataFrame()
    candidate_out = os.path.join(args.outdir, "00_step11G_larger_brainspan_candidate_files.tsv")
    candidate_df.to_csv(candidate_out, sep="\t", index=False)

    # Nested universe analyses.
    coverage_records = []
    stage_records = []
    risk_records = []

    for n in universe_sizes:
        universe = set(gene_var.head(n)["gene"].astype(str))

        # Module coverage and stage projection.
        for (module_family, top_n), sub in mw2.groupby(["module_family", "top_n"], dropna=False):
            module_genes = set(sub["gene"].astype(str))
            in_universe = module_genes & universe

            coverage_records.append({
                "universe_size": n,
                "module_family": module_family,
                "top_n": int(top_n),
                "module_gene_n": len(module_genes),
                "module_genes_in_universe": len(in_universe),
                "module_coverage_fraction": len(in_universe) / len(module_genes) if module_genes else np.nan,
            })

            sub_u = sub[sub["gene"].isin(universe)].copy()
            if sub_u.shape[0] == 0:
                continue

            sp_u = sp[sp["gene"].isin(universe)].copy()
            align = weighted_stage_projection(sub_u[["gene", "weight"]], sp_u)

            if align.shape[0] == 0:
                continue

            align["universe_size"] = n
            align["module_family"] = module_family
            align["top_n"] = int(top_n)

            top_row = align.sort_values("stage_probability", ascending=False).iloc[0]
            align["top_stage"] = top_row["stage_clean"]
            align["top_probability"] = top_row["stage_probability"]

            for _, r in align.iterrows():
                stage_records.append(r.to_dict())

        # Risk/gene-set coverage.
        for set_name, genes in risk_sets.items():
            genes_clean = set(clean_gene(g) for g in genes if clean_gene(g))
            in_universe = genes_clean & universe
            risk_records.append({
                "universe_size": n,
                "gene_set": set_name,
                "gene_set_n": len(genes_clean),
                "genes_in_universe": len(in_universe),
                "coverage_fraction": len(in_universe) / len(genes_clean) if genes_clean else np.nan,
            })

    coverage_df = pd.DataFrame(coverage_records)
    stage_df = pd.DataFrame(stage_records)
    risk_df = pd.DataFrame(risk_records)

    coverage_out = os.path.join(args.outdir, "01_step11G_module_coverage_by_nested_universe.tsv")
    stage_out = os.path.join(args.outdir, "02_step11G_stage_alignment_by_nested_universe.tsv")
    risk_out = os.path.join(args.outdir, "03_step11G_risk_resource_coverage_by_nested_universe.tsv")

    coverage_df.to_csv(coverage_out, sep="\t", index=False)
    stage_df.to_csv(stage_out, sep="\t", index=False)
    risk_df.to_csv(risk_out, sep="\t", index=False)

    # Main call summary.
    if stage_df.shape[0] > 0:
        top_stage_df = (
            stage_df[stage_df["stage_clean"].eq(stage_df["top_stage"])]
            .groupby(["module_family", "top_n", "top_stage"], dropna=False)
            .agg(n_universes=("universe_size", "nunique"),
                 min_universe=("universe_size", "min"),
                 max_universe=("universe_size", "max"),
                 median_top_probability=("top_probability", "median"))
            .reset_index()
        )

        main_summary = (
            stage_df[stage_df["stage_clean"].eq(stage_df["top_stage"])]
            .pivot_table(index=["module_family", "top_n"],
                         columns="top_stage",
                         values="universe_size",
                         aggfunc="nunique",
                         fill_value=0)
            .reset_index()
        )
    else:
        top_stage_df = pd.DataFrame()
        main_summary = pd.DataFrame()

    top_out = os.path.join(args.outdir, "04_step11G_top_stage_stability_by_nested_universe.tsv")
    matrix_out = os.path.join(args.outdir, "05_step11G_top_stage_matrix_by_nested_universe.tsv")
    top_stage_df.to_csv(top_out, sep="\t", index=False)
    main_summary.to_csv(matrix_out, sep="\t", index=False)

    # Compact risk summary for key sets.
    key_patterns = ["SFARI_all", "SFARI_score_le1", "SFARI_S_plus_1_2", "Satterstrom", "Zhou"]
    risk_df["is_key_resource"] = risk_df["gene_set"].apply(lambda x: any(p.lower() in str(x).lower() for p in key_patterns))
    risk_key = risk_df[risk_df["is_key_resource"]].copy()
    risk_key_out = os.path.join(args.outdir, "06_step11G_key_risk_resource_coverage.tsv")
    risk_key.to_csv(risk_key_out, sep="\t", index=False)

    # Decision table for main module calls.
    decision_records = []
    n_universes = len(universe_sizes)
    if stage_df.shape[0] > 0:
        top_calls = stage_df[stage_df["stage_clean"].eq(stage_df["top_stage"])].copy()

        for (module_family, top_n), sub in top_calls.groupby(["module_family", "top_n"]):
            late_count = int((sub["top_stage"] == "late_prenatal").sum())
            mid_count = int((sub["top_stage"] == "mid_prenatal").sum())
            post_count = int(sub["top_stage"].isin(["adolescence", "adulthood"]).sum())

            if module_family in ["NTM1_ASD_up", "NTM3_ASD_signed"] and top_n in [200, 500]:
                status = "PASS" if late_count / n_universes >= 0.80 else "REVIEW"
                criterion = "late_prenatal fraction >= 0.80"
            elif module_family == "NTM2_ASD_down" and top_n in [200, 500]:
                status = "PASS" if late_count == 0 else "REVIEW"
                criterion = "late_prenatal count == 0"
            else:
                status = "INFO"
                criterion = "exploratory threshold"

            decision_records.append({
                "module_family": module_family,
                "top_n": int(top_n),
                "n_universes": n_universes,
                "late_prenatal_top_count": late_count,
                "mid_prenatal_top_count": mid_count,
                "postnatal_top_count": post_count,
                "criterion": criterion,
                "status": status,
            })

    decision_df = pd.DataFrame(decision_records)
    decision_out = os.path.join(args.outdir, "07_step11G_decision_table.tsv")
    decision_df.to_csv(decision_out, sep="\t", index=False)

    # Markdown summary.
    md_out = os.path.join(args.outdir, "08_step11G_overall_summary.md")
    with open(md_out, "w") as f:
        f.write("# Step11G BrainSpan universe coverage and nested-universe sensitivity\n\n")

        f.write("## Purpose\n\n")
        f.write("This analysis quantifies how the current 5,000-gene BrainSpan embedding universe affects module coverage, genetic-resource coverage, and developmental stage alignment. It also tests nested-universe sensitivity within the available 5k universe.\n\n")

        f.write("## Inputs\n\n")
        f.write(f"- BrainSpan embedding: `{args.embedding}`\n")
        f.write(f"- BrainSpan stage profiles: `{args.stage_profiles}`\n")
        f.write(f"- Module weights: `{args.module_weights}`\n")
        f.write(f"- SFARI fixed universe table: `{args.sfari_fixed}`\n")
        f.write(f"- Satterstrom file: `{args.satterstrom}`\n")
        f.write(f"- Zhou gene-set directory: `{args.zhou_gene_set_dir}`\n\n")

        f.write("## Universe design\n\n")
        f.write(f"- Available 5k embedding genes after cleaning: {available_n}\n")
        f.write(f"- Nested universe sizes tested: {', '.join(map(str, universe_sizes))}\n")
        f.write("- Nested universes were ranked by variance across BrainSpan stage-profile z-scores within the current 5k embedding.\n\n")

        f.write("## Main decision table\n\n")
        if decision_df.shape[0] > 0:
            f.write(decision_df.to_markdown(index=False))
        else:
            f.write("No decision table generated.")
        f.write("\n\n")

        f.write("## Top-stage matrix\n\n")
        if main_summary.shape[0] > 0:
            f.write(main_summary.to_markdown(index=False))
        else:
            f.write("No top-stage matrix generated.")
        f.write("\n\n")

        f.write("## Key risk-resource coverage\n\n")
        if risk_key.shape[0] > 0:
            # Show compact subset for 5k and smallest universe.
            show = risk_key[risk_key["universe_size"].isin([min(universe_sizes), max(universe_sizes)])].copy()
            f.write(show.head(80).to_markdown(index=False))
        else:
            f.write("No key risk resources detected.")
        f.write("\n\n")

        f.write("## Larger BrainSpan file audit\n\n")
        if candidate_df.shape[0] > 0:
            f.write("Potential BrainSpan/expression files were detected. Inspect `00_step11G_larger_brainspan_candidate_files.tsv` before attempting 8k/10k reconstruction.\n\n")
            f.write(candidate_df.head(20).to_markdown(index=False))
        else:
            f.write("No obvious larger BrainSpan expression files were detected under the project root.\n")
        f.write("\n\n")

        f.write("## Interpretation guide\n\n")
        f.write("- If NTM1/NTM3 top200/top500 retain late_prenatal across most nested universes, this supports stability within the current 5k embedding.\n")
        f.write("- If NTM2 top200/top500 remains non-late-prenatal, this supports separation from prenatal-transport modules.\n")
        f.write("- This analysis does not fully answer 10k-universe robustness unless raw BrainSpan expression files are available and a larger graph is reconstructed.\n")
        f.write("- If larger BrainSpan files are available, Step11G2 should rebuild 8k/10k embeddings; otherwise, report this as a quantified coverage limitation.\n")

    log("Finished Step11G.")
    for p in [
        candidate_out,
        coverage_out,
        stage_out,
        risk_out,
        top_out,
        matrix_out,
        risk_key_out,
        decision_out,
        md_out,
    ]:
        log(f"Wrote: {p}")


if __name__ == "__main__":
    main()
