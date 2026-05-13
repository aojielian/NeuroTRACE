#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step12D0: Inspect gnomAD v2.1.1 constraint metrics and generate held-out constraint gene sets.

Input:
  /gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz

Outputs:
  00_step12D0_gnomad_file_columns.tsv
  01_step12D0_gnomad_constraint_gene_sets.tsv
  02_step12D0_gnomad_constraint_summary.md
"""

import argparse
import os
import numpy as np
import pandas as pd


def find_col(df, candidates, required=False):
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    if required:
        raise ValueError(f"Missing required column from {candidates}. Available columns: {list(df.columns)}")
    return None


def clean_gene(x):
    x = str(x).strip()
    x = x.replace("gene:", "")
    x = x.replace("GENE:", "")
    return x


def safe_numeric(s):
    return pd.to_numeric(s, errors="coerce")


def add_set(records, df, set_name, mask, gene_col, description):
    sub = df.loc[mask & df[gene_col].notna(), gene_col].astype(str).map(clean_gene)
    sub = sorted(set([g for g in sub if g and g.lower() != "nan"]))
    for g in sub:
        records.append({
            "gene_set": set_name,
            "gene": g,
            "description": description,
            "n_genes_in_set": len(sub)
        })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--constraint_file",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
    )
    parser.add_argument(
        "--outdir",
        default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step12_revision_strengthening/results/step12D_gnomad_constraint"
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("[Step12D0] Reading:", args.constraint_file, flush=True)
    df = pd.read_csv(args.constraint_file, sep="\t", compression="gzip", low_memory=False)

    print("[Step12D0] Shape:", df.shape, flush=True)
    print("[Step12D0] Columns:")
    for i, c in enumerate(df.columns, 1):
        print(i, c)

    gene_col = find_col(df, ["gene", "gene_symbol", "symbol", "hgnc_symbol"], required=True)
    gene_id_col = find_col(df, ["gene_id", "ensg", "ensembl_gene_id"])
    pli_col = find_col(df, ["pLI", "pli"])
    loeuf_col = find_col(df, ["oe_lof_upper", "loeuf", "LOEUF", "lof.oe_ci.upper", "oe_lof_upper_bin"])
    oe_lof_col = find_col(df, ["oe_lof", "obs_exp_lof", "lof.oe"])
    mis_z_col = find_col(df, ["mis_z", "misZ", "missense_z", "mis_z_score"])
    syn_z_col = find_col(df, ["syn_z", "synZ", "syn_z_score"])
    constraint_flag_col = find_col(df, ["constraint_flag", "flag", "flags"])

    detected = {
        "gene_col": gene_col,
        "gene_id_col": gene_id_col,
        "pLI_col": pli_col,
        "LOEUF_col": loeuf_col,
        "oe_lof_col": oe_lof_col,
        "mis_z_col": mis_z_col,
        "syn_z_col": syn_z_col,
        "constraint_flag_col": constraint_flag_col,
    }

    col_summary = []
    for c in df.columns:
        s = df[c]
        col_summary.append({
            "column": c,
            "dtype": str(s.dtype),
            "n_non_missing": int(s.notna().sum()),
            "n_unique": int(s.nunique(dropna=True)),
            "example_values": ",".join(s.dropna().astype(str).head(5).tolist())
        })

    col_summary = pd.DataFrame(col_summary)
    col_summary_out = os.path.join(args.outdir, "00_step12D0_gnomad_file_columns.tsv")
    col_summary.to_csv(col_summary_out, sep="\t", index=False)

    # Clean and prepare numeric columns.
    df["_gene_clean"] = df[gene_col].astype(str).map(clean_gene)

    if pli_col is not None:
        df["_pLI"] = safe_numeric(df[pli_col])
    else:
        df["_pLI"] = np.nan

    if loeuf_col is not None:
        df["_LOEUF"] = safe_numeric(df[loeuf_col])
    else:
        df["_LOEUF"] = np.nan

    if oe_lof_col is not None:
        df["_oe_lof"] = safe_numeric(df[oe_lof_col])
    else:
        df["_oe_lof"] = np.nan

    if mis_z_col is not None:
        df["_mis_z"] = safe_numeric(df[mis_z_col])
    else:
        df["_mis_z"] = np.nan

    if syn_z_col is not None:
        df["_syn_z"] = safe_numeric(df[syn_z_col])
    else:
        df["_syn_z"] = np.nan

    records = []

    # Primary v2.1.1 recommended constrained sets.
    if df["_LOEUF"].notna().sum() > 0:
        add_set(
            records, df, "gnomAD_v2.1.1_LOEUF_lt_0.35",
            df["_LOEUF"] < 0.35,
            "_gene_clean",
            "Established gnomAD v2.1.1 LoF-constrained genes using LOEUF/oe_lof_upper < 0.35"
        )

        add_set(
            records, df, "gnomAD_v2.1.1_LOEUF_lt_0.20",
            df["_LOEUF"] < 0.20,
            "_gene_clean",
            "Stricter LoF-constrained genes using LOEUF/oe_lof_upper < 0.20"
        )

        q10 = df["_LOEUF"].quantile(0.10)
        q20 = df["_LOEUF"].quantile(0.20)

        add_set(
            records, df, "gnomAD_v2.1.1_LOEUF_lowest_10pct",
            df["_LOEUF"] <= q10,
            "_gene_clean",
            f"Lowest 10 percent of LOEUF/oe_lof_upper values, cutoff={q10:.4g}"
        )

        add_set(
            records, df, "gnomAD_v2.1.1_LOEUF_lowest_20pct",
            df["_LOEUF"] <= q20,
            "_gene_clean",
            f"Lowest 20 percent of LOEUF/oe_lof_upper values, cutoff={q20:.4g}"
        )

    if df["_pLI"].notna().sum() > 0:
        add_set(
            records, df, "gnomAD_v2.1.1_pLI_gt_0.90",
            df["_pLI"] > 0.90,
            "_gene_clean",
            "LoF-intolerant genes using pLI > 0.90"
        )

        add_set(
            records, df, "gnomAD_v2.1.1_pLI_gt_0.99",
            df["_pLI"] > 0.99,
            "_gene_clean",
            "Very high pLI genes using pLI > 0.99"
        )

    if df["_mis_z"].notna().sum() > 0:
        q90 = df["_mis_z"].quantile(0.90)
        q95 = df["_mis_z"].quantile(0.95)

        add_set(
            records, df, "gnomAD_v2.1.1_misZ_top_10pct",
            df["_mis_z"] >= q90,
            "_gene_clean",
            f"Top 10 percent missense Z genes, cutoff={q90:.4g}"
        )

        add_set(
            records, df, "gnomAD_v2.1.1_misZ_top_5pct",
            df["_mis_z"] >= q95,
            "_gene_clean",
            f"Top 5 percent missense Z genes, cutoff={q95:.4g}"
        )

    if df["_syn_z"].notna().sum() > 0:
        q90 = df["_syn_z"].quantile(0.90)
        add_set(
            records, df, "gnomAD_v2.1.1_synZ_top_10pct",
            df["_syn_z"] >= q90,
            "_gene_clean",
            f"Top 10 percent synonymous Z genes, cutoff={q90:.4g}; mainly used as negative/sensitivity control"
        )

    gene_sets = pd.DataFrame(records)
    gene_sets_out = os.path.join(args.outdir, "01_step12D0_gnomad_constraint_gene_sets.tsv")
    gene_sets.to_csv(gene_sets_out, sep="\t", index=False)

    set_summary = (
        gene_sets.groupby(["gene_set", "description"], dropna=False)
        .agg(n_genes=("gene", "nunique"))
        .reset_index()
        .sort_values("gene_set")
    )

    set_summary_out = os.path.join(args.outdir, "01b_step12D0_gnomad_constraint_gene_set_summary.tsv")
    set_summary.to_csv(set_summary_out, sep="\t", index=False)

    md_out = os.path.join(args.outdir, "02_step12D0_gnomad_constraint_summary.md")
    with open(md_out, "w") as f:
        f.write("# Step12D0 gnomAD v2.1.1 constraint file inspection\n\n")

        f.write("## Input\n\n")
        f.write(f"- Constraint file: `{args.constraint_file}`\n")
        f.write(f"- Rows: {df.shape[0]}\n")
        f.write(f"- Columns: {df.shape[1]}\n\n")

        f.write("## Detected key columns\n\n")
        for k, v in detected.items():
            f.write(f"- {k}: `{v}`\n")
        f.write("\n")

        f.write("## First five rows\n\n")
        f.write(df.head(5).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Generated gene sets\n\n")
        f.write(set_summary.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Recommended use in Step12D\n\n")
        f.write("- Primary held-out constraint set: `gnomAD_v2.1.1_LOEUF_lt_0.35`\n")
        f.write("- Sensitivity sets: `pLI_gt_0.90`, `LOEUF_lowest_10pct`, `misZ_top_10pct`\n")
        f.write("- Synonymous Z set, if generated, should be interpreted as a control/sensitivity set rather than a disease-relevant constraint set.\n")

    print("[Step12D0] Wrote:", col_summary_out)
    print("[Step12D0] Wrote:", gene_sets_out)
    print("[Step12D0] Wrote:", set_summary_out)
    print("[Step12D0] Wrote:", md_out)


if __name__ == "__main__":
    main()
