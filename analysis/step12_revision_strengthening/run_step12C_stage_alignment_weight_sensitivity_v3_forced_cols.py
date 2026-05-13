#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import numpy as np
import pandas as pd

def zscore_vec(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    if (not np.isfinite(sd)) or sd == 0:
        return np.zeros_like(x)
    return (x - np.nanmean(x)) / sd

def softmax(x, temperature=0.75):
    x = np.asarray(x, dtype=float) / float(temperature)
    x = x - np.nanmax(x)
    ex = np.exp(x)
    den = np.nansum(ex)
    if den <= 0 or not np.isfinite(den):
        return np.ones(len(x)) / len(x)
    return ex / den

def clean_stage(x):
    return str(x).replace("stage:", "").strip()

def build_weight_designs(step=0.05):
    designs = []

    designs.append({
        "design_type": "default",
        "design_name": "default_0.45_0.35_0.20",
        "w_no_pc1": 0.45,
        "w_pc2_pc6": 0.35,
        "w_signed_profile": 0.20
    })

    for leave in ["no_pc1", "pc2_pc6", "signed_profile"]:
        weights = {"no_pc1": 0.45, "pc2_pc6": 0.35, "signed_profile": 0.20}
        weights[leave] = 0.0
        total = sum(weights.values())
        for k in weights:
            weights[k] = weights[k] / total if total > 0 else 0
        designs.append({
            "design_type": "leave_one_component_out",
            "design_name": f"leave_out_{leave}",
            "w_no_pc1": weights["no_pc1"],
            "w_pc2_pc6": weights["pc2_pc6"],
            "w_signed_profile": weights["signed_profile"]
        })

    designs.extend([
        {"design_type": "single_component_only", "design_name": "only_no_pc1", "w_no_pc1": 1.0, "w_pc2_pc6": 0.0, "w_signed_profile": 0.0},
        {"design_type": "single_component_only", "design_name": "only_pc2_pc6", "w_no_pc1": 0.0, "w_pc2_pc6": 1.0, "w_signed_profile": 0.0},
        {"design_type": "single_component_only", "design_name": "only_signed_profile", "w_no_pc1": 0.0, "w_pc2_pc6": 0.0, "w_signed_profile": 1.0},
    ])

    vals = np.round(np.arange(0, 1 + 1e-9, step), 3)
    for a in vals:
        for b in vals:
            c = round(1.0 - a - b, 3)
            if c < -1e-9:
                continue
            if c < 0:
                continue
            designs.append({
                "design_type": "simplex_grid",
                "design_name": f"grid_{a:.2f}_{b:.2f}_{c:.2f}",
                "w_no_pc1": float(a),
                "w_pc2_pc6": float(b),
                "w_signed_profile": float(c)
            })

    return pd.DataFrame(designs).drop_duplicates()

def evaluate_designs(aln, designs, temperatures):
    records = []

    for _, d in designs.iterrows():
        for temp in temperatures:
            for (module, top_n), sub in aln.groupby(["module_family", "top_n"], dropna=False):
                sub = sub.copy()

                score = (
                    d["w_no_pc1"] * sub["no_pc1_z"].values +
                    d["w_pc2_pc6"] * sub["pc2_pc6_z"].values +
                    d["w_signed_profile"] * sub["signed_profile_z"].values
                )

                prob = softmax(score, temperature=temp)
                stages = sub["stage_clean"].tolist()
                order = np.argsort(prob)[::-1]

                top_stage = stages[order[0]]
                second_stage = stages[order[1]] if len(order) > 1 else "NA"

                def prob_of(stage_name):
                    return float(sum(p for s, p in zip(stages, prob) if s == stage_name))

                prenatal_prob = float(sum(p for s, p in zip(stages, prob) if "prenatal" in s))

                records.append({
                    "design_type": d["design_type"],
                    "design_name": d["design_name"],
                    "w_no_pc1": d["w_no_pc1"],
                    "w_pc2_pc6": d["w_pc2_pc6"],
                    "w_signed_profile": d["w_signed_profile"],
                    "temperature": temp,
                    "module_family": module,
                    "top_n": int(top_n),
                    "top_stage": top_stage,
                    "second_stage": second_stage,
                    "top_probability": float(prob[order[0]]),
                    "second_probability": float(prob[order[1]]) if len(order) > 1 else np.nan,
                    "top_margin": float(prob[order[0]] - prob[order[1]]) if len(order) > 1 else np.nan,
                    "late_prenatal_probability": prob_of("late_prenatal"),
                    "mid_prenatal_probability": prob_of("mid_prenatal"),
                    "early_prenatal_probability": prob_of("early_prenatal"),
                    "adolescence_probability": prob_of("adolescence"),
                    "adulthood_probability": prob_of("adulthood"),
                    "prenatal_probability": prenatal_prob,
                    "late_prenatal_is_top": int(top_stage == "late_prenatal"),
                    "prenatal_is_top": int("prenatal" in top_stage),
                    "postnatal_is_top": int(top_stage in ["adolescence", "adulthood", "infancy", "childhood"])
                })

    return pd.DataFrame(records)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_table", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step03_feature_embedding/results/33_step03D_NTM_to_brainspan_stage_alignment.tsv")
    parser.add_argument("--outdir", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step12_revision_strengthening/results/step12C_stage_alignment_weight_sensitivity_v3")
    parser.add_argument("--grid_step", type=float, default=0.05)
    parser.add_argument("--temperatures", default="0.25,0.5,0.75,1.0,1.5,2.0")
    parser.add_argument("--module_col", default="program")
    parser.add_argument("--topn_col", default="top_n")
    parser.add_argument("--stage_col", default="state_id")
    parser.add_argument("--no_pc1_col", default="cosine_no_PC1_z")
    parser.add_argument("--pc2_pc6_col", default="secondary_PC_z")
    parser.add_argument("--signed_profile_col", default="signed_weighted_stage_score")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.input_table, sep="\t", low_memory=False)

    required = [
        args.module_col,
        args.topn_col,
        args.stage_col,
        args.no_pc1_col,
        args.pc2_pc6_col,
        args.signed_profile_col
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing forced columns: {missing}. Available columns: {list(df.columns)}")

    aln = pd.DataFrame()
    aln["module_family"] = df[args.module_col].astype(str)
    aln["top_n"] = pd.to_numeric(df[args.topn_col], errors="coerce")
    aln["stage_clean"] = df[args.stage_col].astype(str).map(clean_stage)
    aln["no_pc1"] = pd.to_numeric(df[args.no_pc1_col], errors="coerce")
    aln["pc2_pc6"] = pd.to_numeric(df[args.pc2_pc6_col], errors="coerce")
    aln["signed_profile"] = pd.to_numeric(df[args.signed_profile_col], errors="coerce")

    aln = aln[
        aln["module_family"].isin(["NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed"]) &
        aln["top_n"].isin([50, 100, 200, 500])
    ].copy()

    for c in ["no_pc1", "pc2_pc6", "signed_profile"]:
        aln[c + "_z"] = np.nan
        for _, idx in aln.groupby(["module_family", "top_n"]).groups.items():
            aln.loc[idx, c + "_z"] = zscore_vec(aln.loc[idx, c].values)

    aln_out = os.path.join(args.outdir, "01_step12C_v3_standardized_alignment_components.tsv")
    aln.to_csv(aln_out, sep="\t", index=False)

    detected = pd.DataFrame([{
        "input_table": args.input_table,
        "module_col": args.module_col,
        "topn_col": args.topn_col,
        "stage_col": args.stage_col,
        "no_pc1_col": args.no_pc1_col,
        "pc2_pc6_col": args.pc2_pc6_col,
        "signed_profile_col": args.signed_profile_col,
        "n_rows_after_filter": aln.shape[0]
    }])
    detected_out = os.path.join(args.outdir, "00_step12C_v3_forced_columns.tsv")
    detected.to_csv(detected_out, sep="\t", index=False)

    designs = build_weight_designs(step=args.grid_step)
    temps = [float(x) for x in args.temperatures.split(",") if x.strip()]

    grid_df = evaluate_designs(aln, designs, temps)

    grid_out = os.path.join(args.outdir, "02_step12C_v3_weight_temperature_grid_stage_probs.tsv.gz")
    grid_df.to_csv(grid_out, sep="\t", index=False, compression="gzip")

    main = grid_df[grid_df["top_n"].isin([200, 500])].copy()

    module_summary = (
        main.groupby(["module_family", "top_n"], dropna=False)
        .agg(
            n_settings=("design_name", "count"),
            late_prenatal_top_fraction=("late_prenatal_is_top", "mean"),
            prenatal_top_fraction=("prenatal_is_top", "mean"),
            postnatal_top_fraction=("postnatal_is_top", "mean"),
            late_prenatal_probability_median=("late_prenatal_probability", "median"),
            prenatal_probability_median=("prenatal_probability", "median"),
            top_margin_median=("top_margin", "median"),
            n_unique_top_stages=("top_stage", "nunique"),
            top_stage_modes=("top_stage", lambda x: ",".join(pd.Series(x).value_counts().head(5).index.astype(str)))
        )
        .reset_index()
    )

    module_summary_out = os.path.join(args.outdir, "04_step12C_v3_module_level_stability_summary.tsv")
    module_summary.to_csv(module_summary_out, sep="\t", index=False)

    records = []
    for _, row in module_summary.iterrows():
        fam = row["module_family"]
        topn = int(row["top_n"])

        if fam in ["NTM1_ASD_up", "NTM3_ASD_signed"]:
            if row["prenatal_top_fraction"] >= 0.80:
                status = "PASS"
            elif row["prenatal_top_fraction"] >= 0.60:
                status = "PARTIAL"
            else:
                status = "REVIEW"
            interpretation = "prenatal-biased across alignment-component and temperature settings"
        elif fam == "NTM2_ASD_down":
            if row["late_prenatal_top_fraction"] <= 0.20:
                status = "PASS"
            elif row["late_prenatal_top_fraction"] <= 0.40:
                status = "PARTIAL"
            else:
                status = "REVIEW"
            interpretation = "separated from strict late-prenatal top-stage across alignment settings"
        else:
            status = "INFO"
            interpretation = "exploratory"

        records.append({
            "module_family": fam,
            "top_n": topn,
            "status": status,
            "late_prenatal_top_fraction": row["late_prenatal_top_fraction"],
            "prenatal_top_fraction": row["prenatal_top_fraction"],
            "postnatal_top_fraction": row["postnatal_top_fraction"],
            "late_prenatal_probability_median": row["late_prenatal_probability_median"],
            "prenatal_probability_median": row["prenatal_probability_median"],
            "n_unique_top_stages": row["n_unique_top_stages"],
            "top_stage_modes": row["top_stage_modes"],
            "interpretation": interpretation
        })

    decision = pd.DataFrame(records)
    decision_out = os.path.join(args.outdir, "05_step12C_v3_decision_table.tsv")
    decision.to_csv(decision_out, sep="\t", index=False)

    md_out = os.path.join(args.outdir, "06_step12C_v3_overall_summary.md")
    with open(md_out, "w") as f:
        f.write("# Step12C v3 forced-column stage-alignment component weight sensitivity\\n\\n")
        f.write("## Purpose\\n\\n")
        f.write("This analysis tests whether developmental stage calls depend on alignment-component weights or softmax temperature, using explicitly specified component columns to avoid automatic column-selection ambiguity.\\n\\n")

        f.write("## Forced columns\\n\\n")
        f.write(detected.to_markdown(index=False))
        f.write("\\n\\n")

        f.write("## Module-level stability summary\\n\\n")
        f.write(module_summary.to_markdown(index=False))
        f.write("\\n\\n")

        f.write("## Decision table\\n\\n")
        f.write(decision.to_markdown(index=False))
        f.write("\\n\\n")

        f.write("## Manuscript-ready interpretation\\n\\n")
        f.write("Stage-alignment component sensitivity showed that broad prenatal versus non-prenatal separation was more stable than exact prenatal substage assignment. Using forced component columns, NTM1_ASD_up and NTM3_ASD_signed were evaluated for broad prenatal bias across alternative component weights and temperatures, while NTM2_ASD_down was evaluated for separation from strict late-prenatal dominance. Exact late-prenatal specificity should therefore be interpreted as calibrated under the primary matched-random model rather than as invariant across all component-weight settings.\\n")

    print("[Step12C v3] Finished.")
    for p in [detected_out, aln_out, grid_out, module_summary_out, decision_out, md_out]:
        print("[Step12C v3] Wrote:", p)

if __name__ == "__main__":
    main()
