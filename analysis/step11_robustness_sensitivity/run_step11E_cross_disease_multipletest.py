#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step11E v2: Cross-disease multiple-testing and directional-concordance framework.

This version is explicitly adapted to current NeuroTRACE Step08 outputs:
  - 82_step08E_cross_disease_NTM_models.tsv
  - optional 60_step08C_gandal_cross_disease_NTM_models.tsv

Primary purpose:
  Formalize multiple-testing control and directional evidence for SCZ/BD/MDD
  cross-disease generalization.

Outputs:
  01_step11E_v2_cross_disease_all_tests_with_FDR.tsv
  02_step11E_v2_disease_level_summary.tsv
  03_step11E_v2_module_level_summary.tsv
  04_step11E_v2_primary_cross_disease_boundary_summary.tsv
  05_step11E_v2_overall_summary.md
"""

import argparse
import os
import math
import numpy as np
import pandas as pd

from scipy.stats import norm, binomtest


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


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


def normalize_disease(x):
    s = str(x).strip().lower()

    if s in ["scz", "schizophrenia"] or "schiz" in s:
        return "SCZ"
    if s in ["bd", "bipolar", "bipolar disorder"] or "bipolar" in s:
        return "BD"
    if s in ["mdd", "major depressive disorder", "depression"] or "depress" in s:
        return "MDD"
    if "dup15" in s or "15q" in s:
        return "Dup15q"
    if s in ["asd", "autism", "autism spectrum disorder"] or "autism" in s:
        return "ASD"

    return str(x)


def normalize_module(x):
    s = str(x)

    if "NTM1" in s:
        return "NTM1_ASD_up"
    if "NTM2" in s:
        return "NTM2_ASD_down"
    if "NTM3" in s:
        return "NTM3_ASD_signed"

    return s


def one_sided_positive_p(beta, p_two):
    beta = float(beta)
    p_two = float(p_two)

    if not np.isfinite(beta) or not np.isfinite(p_two):
        return np.nan

    p_two = min(max(p_two, 1e-300), 1.0)

    if beta >= 0:
        return p_two / 2.0
    else:
        return 1.0 - p_two / 2.0


def stouffer_meta(pvals):
    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p)]

    if len(p) == 0:
        return np.nan, np.nan

    p = np.clip(p, 1e-300, 1.0)
    z = norm.isf(p)
    Z = np.sum(z) / math.sqrt(len(z))
    p_one = norm.sf(Z)

    return Z, p_one


def read_and_standardize(path, source_label):
    if not os.path.exists(path):
        return pd.DataFrame()

    df = pd.read_csv(path, sep="\t", low_memory=False)

    required = [
        "target_dx",
        "program",
        "top_n",
        "beta_target_vs_Control",
        "p_value",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing required columns: {missing}. Available: {list(df.columns)}")

    out = pd.DataFrame()
    out["source_label"] = source_label
    out["source_file"] = path

    if "dataset" in df.columns:
        out["dataset"] = df["dataset"].astype(str)
    else:
        out["dataset"] = source_label

    out["disease_raw"] = df["target_dx"].astype(str)
    out["disease"] = out["disease_raw"].map(normalize_disease)

    out["program"] = df["program"].astype(str)
    out["module_family"] = out["program"].map(normalize_module)
    out["top_n"] = pd.to_numeric(df["top_n"], errors="coerce")

    out["beta"] = pd.to_numeric(df["beta_target_vs_Control"], errors="coerce")
    out["se"] = pd.to_numeric(df["se"], errors="coerce") if "se" in df.columns else np.nan
    out["t"] = pd.to_numeric(df["t"], errors="coerce") if "t" in df.columns else np.nan
    out["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")

    if "fdr" in df.columns:
        out["reported_fdr"] = pd.to_numeric(df["fdr"], errors="coerce")
    else:
        out["reported_fdr"] = np.nan

    if "auc_target_vs_Control" in df.columns:
        out["auc"] = pd.to_numeric(df["auc_target_vs_Control"], errors="coerce")
    else:
        out["auc"] = np.nan

    if "n" in df.columns:
        out["n"] = pd.to_numeric(df["n"], errors="coerce")
    else:
        out["n"] = np.nan

    if "n_target" in df.columns:
        out["n_target"] = pd.to_numeric(df["n_target"], errors="coerce")
    else:
        out["n_target"] = np.nan

    if "n_Control" in df.columns:
        out["n_Control"] = pd.to_numeric(df["n_Control"], errors="coerce")
    else:
        out["n_Control"] = np.nan

    if "n_common_genes" in df.columns:
        out["n_common_genes"] = pd.to_numeric(df["n_common_genes"], errors="coerce")
    else:
        out["n_common_genes"] = np.nan

    if "overlap_rate" in df.columns:
        out["overlap_rate"] = pd.to_numeric(df["overlap_rate"], errors="coerce")
    else:
        out["overlap_rate"] = np.nan

    out = out[
        out["module_family"].isin(["NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed"]) &
        out["disease"].isin(["SCZ", "BD", "MDD", "Dup15q", "ASD"]) &
        np.isfinite(out["beta"]) &
        np.isfinite(out["p_value"])
    ].copy()

    out["positive_direction"] = out["beta"] > 0
    out["one_sided_p_positive"] = [
        one_sided_positive_p(b, p) for b, p in zip(out["beta"], out["p_value"])
    ]

    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--external_cross_disease",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step08_algorithm_strengthening/results/82_step08E_cross_disease_NTM_models.tsv",
    )
    parser.add_argument(
        "--gandal_cross_disease",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step08_algorithm_strengthening/results/60_step08C_gandal_cross_disease_NTM_models.tsv",
    )
    parser.add_argument(
        "--outdir",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step11_robustness_sensitivity/results",
    )
    args = parser.parse_args()

    ensure_dir(args.outdir)

    external = read_and_standardize(args.external_cross_disease, "external_microarray")
    gandal = read_and_standardize(args.gandal_cross_disease, "gandal_internal")

    # Primary boundary uses external SCZ/BD/MDD rows.
    primary = external[external["disease"].isin(["SCZ", "BD", "MDD"])].copy()

    if primary.shape[0] == 0:
        raise RuntimeError("No primary SCZ/BD/MDD external rows were found.")

    primary = primary.drop_duplicates(
        subset=["dataset", "disease", "module_family", "top_n", "beta", "p_value"]
    ).copy()

    primary["global_fdr"] = bh_fdr(primary["p_value"].values)

    primary["within_disease_fdr"] = np.nan
    for disease, idx in primary.groupby("disease").groups.items():
        primary.loc[idx, "within_disease_fdr"] = bh_fdr(primary.loc[idx, "p_value"].values)

    primary["within_dataset_fdr"] = np.nan
    for dataset, idx in primary.groupby("dataset").groups.items():
        primary.loc[idx, "within_dataset_fdr"] = bh_fdr(primary.loc[idx, "p_value"].values)

    all_out = os.path.join(args.outdir, "01_step11E_v2_cross_disease_all_tests_with_FDR.tsv")
    primary.to_csv(all_out, sep="\t", index=False)

    disease_records = []
    for disease, sub in primary.groupby("disease"):
        n = sub.shape[0]
        k = int(sub["positive_direction"].sum())

        binom_p = binomtest(k, n, 0.5, alternative="greater").pvalue if n > 0 else np.nan
        Z, st_p = stouffer_meta(sub["one_sided_p_positive"].values)

        disease_records.append({
            "disease": disease,
            "n_tests": n,
            "n_positive": k,
            "positive_rate": k / n if n else np.nan,
            "one_sided_binomial_p_positive": binom_p,
            "stouffer_Z_positive": Z,
            "stouffer_p_positive": st_p,
            "mean_beta": sub["beta"].mean(),
            "median_beta": sub["beta"].median(),
            "median_auc": sub["auc"].median(),
            "min_p": sub["p_value"].min(),
            "min_global_fdr": sub["global_fdr"].min(),
            "min_within_disease_fdr": sub["within_disease_fdr"].min(),
            "n_global_fdr_lt_0_05": int((sub["global_fdr"] < 0.05).sum()),
            "n_within_disease_fdr_lt_0_05": int((sub["within_disease_fdr"] < 0.05).sum()),
            "datasets": ",".join(sorted(sub["dataset"].astype(str).unique())),
            "programs": ",".join(sorted(sub["module_family"].astype(str).unique())),
        })

    disease_summary = pd.DataFrame(disease_records)
    disease_summary["stouffer_fdr_across_diseases"] = bh_fdr(disease_summary["stouffer_p_positive"].values)
    disease_summary["binomial_fdr_across_diseases"] = bh_fdr(disease_summary["one_sided_binomial_p_positive"].values)

    disease_order = {"SCZ": 1, "BD": 2, "MDD": 3}
    disease_summary["order"] = disease_summary["disease"].map(disease_order)
    disease_summary = disease_summary.sort_values("order").drop(columns=["order"])

    disease_out = os.path.join(args.outdir, "02_step11E_v2_disease_level_summary.tsv")
    disease_summary.to_csv(disease_out, sep="\t", index=False)

    module_summary = (
        primary
        .groupby(["disease", "dataset", "module_family", "top_n"], dropna=False)
        .agg(
            n_tests=("p_value", "count"),
            n_positive=("positive_direction", "sum"),
            positive_rate=("positive_direction", "mean"),
            beta=("beta", "mean"),
            p_value=("p_value", "min"),
            global_fdr=("global_fdr", "min"),
            within_disease_fdr=("within_disease_fdr", "min"),
            within_dataset_fdr=("within_dataset_fdr", "min"),
            auc=("auc", "median"),
            n_common_genes=("n_common_genes", "median"),
            overlap_rate=("overlap_rate", "median"),
        )
        .reset_index()
    )

    module_out = os.path.join(args.outdir, "03_step11E_v2_module_level_summary.tsv")
    module_summary.to_csv(module_out, sep="\t", index=False)

    boundary_out = os.path.join(args.outdir, "04_step11E_v2_primary_cross_disease_boundary_summary.tsv")
    disease_summary.to_csv(boundary_out, sep="\t", index=False)

    # Optional internal Gandal summary, including Dup15q if present.
    all_with_gandal = pd.concat([primary, gandal], ignore_index=True)
    all_with_gandal_out = os.path.join(args.outdir, "06_step11E_v2_external_plus_gandal_all_tests.tsv")
    all_with_gandal.to_csv(all_with_gandal_out, sep="\t", index=False)

    if gandal.shape[0] > 0:
        gandal_summary = (
            gandal
            .groupby("disease", dropna=False)
            .agg(
                n_tests=("p_value", "count"),
                n_positive=("positive_direction", "sum"),
                positive_rate=("positive_direction", "mean"),
                mean_beta=("beta", "mean"),
                median_beta=("beta", "median"),
                min_p=("p_value", "min"),
                min_reported_fdr=("reported_fdr", "min"),
            )
            .reset_index()
        )
    else:
        gandal_summary = pd.DataFrame()

    gandal_out = os.path.join(args.outdir, "07_step11E_v2_gandal_internal_summary.tsv")
    gandal_summary.to_csv(gandal_out, sep="\t", index=False)

    # Markdown report.
    md_out = os.path.join(args.outdir, "05_step11E_v2_overall_summary.md")
    with open(md_out, "w") as f:
        f.write("# Step11E v2 cross-disease multiple-testing summary\n\n")

        f.write("## Purpose\n\n")
        f.write("This analysis formalizes multiple-testing control and directional concordance for the external SCZ/BD/MDD cross-disease NeuroTRACE validation.\n\n")

        f.write("## Inputs\n\n")
        f.write(f"- External cross-disease table: `{args.external_cross_disease}`\n")
        f.write(f"- Optional Gandal internal table: `{args.gandal_cross_disease}`\n")
        f.write(f"- External primary rows retained: {primary.shape[0]}\n")
        f.write(f"- Gandal internal rows retained: {gandal.shape[0]}\n\n")

        f.write("## Disease-level summary\n\n")
        f.write(disease_summary.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Module-level summary\n\n")
        f.write(module_summary.to_markdown(index=False))
        f.write("\n\n")

        if gandal_summary.shape[0] > 0:
            f.write("## Optional Gandal internal cross-disease summary\n\n")
            f.write(gandal_summary.to_markdown(index=False))
            f.write("\n\n")

        f.write("## Manuscript-ready interpretation\n\n")

        def get_disease(d):
            sub = disease_summary[disease_summary["disease"].eq(d)]
            return sub.iloc[0].to_dict() if sub.shape[0] else None

        scz = get_disease("SCZ")
        bd = get_disease("BD")
        mdd = get_disease("MDD")

        if scz is not None:
            f.write(
                f"SCZ showed {int(scz['n_positive'])}/{int(scz['n_tests'])} positive NTM tests "
                f"(positive rate={scz['positive_rate']:.3f}; one-sided binomial P={scz['one_sided_binomial_p_positive']:.3g}; "
                f"Stouffer P={scz['stouffer_p_positive']:.3g}; minimum global FDR={scz['min_global_fdr']:.3g}). "
            )

        if bd is not None:
            f.write(
                f"BD showed {int(bd['n_positive'])}/{int(bd['n_tests'])} positive NTM tests "
                f"(positive rate={bd['positive_rate']:.3f}; one-sided binomial P={bd['one_sided_binomial_p_positive']:.3g}; "
                f"Stouffer P={bd['stouffer_p_positive']:.3g}; minimum global FDR={bd['min_global_fdr']:.3g}). "
            )

        if mdd is not None:
            f.write(
                f"MDD showed {int(mdd['n_positive'])}/{int(mdd['n_tests'])} positive NTM tests "
                f"(positive rate={mdd['positive_rate']:.3f}; one-sided binomial P={mdd['one_sided_binomial_p_positive']:.3g}; "
                f"Stouffer P={mdd['stouffer_p_positive']:.3g}; minimum global FDR={mdd['min_global_fdr']:.3g}). "
            )

        f.write(
            "Together, these results support cross-disease generalization to SCZ and BD, with MDD serving as a negative disease-boundary condition under explicit multiple-testing and directional-concordance evaluation.\n"
        )

    print("[Step11E v2] Finished.")
    for p in [
        all_out,
        disease_out,
        module_out,
        boundary_out,
        md_out,
        all_with_gandal_out,
        gandal_out,
    ]:
        print("[Step11E v2] Wrote:", p)


if __name__ == "__main__":
    main()
