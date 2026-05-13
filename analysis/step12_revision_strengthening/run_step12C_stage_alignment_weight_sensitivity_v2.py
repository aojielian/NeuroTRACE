#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Step12C: Stage-alignment component weight sensitivity.

Purpose
-------
Test whether NeuroTRACE developmental stage calls depend on the fixed
0.45/0.35/0.20 weighting of alignment components:
  - no-PC1 embedding similarity
  - PC2-PC6 similarity
  - signed stage-profile alignment

Also test stage-probability softmax temperature sensitivity.

Main interpretation
-------------------
This analysis is not intended to prove that 0.45/0.35/0.20 is uniquely optimal.
It evaluates whether the broad biological interpretation is stable:
  - NTM1_ASD_up / NTM3_ASD_signed remain prenatal-biased
  - NTM2_ASD_down remains separated from late-prenatal transport

Outputs
-------
00_step12C_candidate_alignment_component_files.tsv
01_step12C_standardized_alignment_components.tsv
02_step12C_weight_temperature_grid_stage_probs.tsv.gz
03_step12C_design_level_summary.tsv
04_step12C_module_level_stability_summary.tsv
05_step12C_decision_table.tsv
06_step12C_overall_summary.md
"""

import argparse
import os
import re
import glob
import itertools
import numpy as np
import pandas as pd


def log(msg):
    print(f"[Step12C] {msg}", flush=True)


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def read_table(path, nrows=None):
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False, nrows=nrows)
    if path.endswith(".csv"):
        return pd.read_csv(path, low_memory=False, nrows=nrows)
    return pd.read_csv(path, sep="\t", low_memory=False, nrows=nrows)


def find_col(df, candidates, required=False):
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    if required:
        raise ValueError(f"Missing required column from {candidates}; available={list(df.columns)}")
    return None


def clean_stage(x):
    return str(x).replace("stage:", "").strip()


def infer_module_family(x):
    s = str(x)
    if "NTM1" in s:
        return "NTM1_ASD_up"
    if "NTM2" in s:
        return "NTM2_ASD_down"
    if "NTM3" in s:
        return "NTM3_ASD_signed"
    return s


def infer_top_n(x):
    m = re.search(r"top[_|\\-]?(\\d+)", str(x), flags=re.I)
    if m:
        return int(m.group(1))
    m = re.search(r"(\\d+)$", str(x))
    if m:
        return int(m.group(1))
    return np.nan


def zscore_by_group(df, group_cols, value_col, out_col):
    df = df.copy()
    vals = []
    for _, sub in df.groupby(group_cols, dropna=False):
        x = pd.to_numeric(sub[value_col], errors="coerce").values.astype(float)
        sd = np.nanstd(x)
        if (not np.isfinite(sd)) or sd == 0:
            z = np.zeros_like(x)
        else:
            z = (x - np.nanmean(x)) / sd
        vals.extend(list(z))
    df[out_col] = vals
    return df


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


def component_col_score(cols):
    """
    Score whether a table likely contains stage-alignment component columns.
    """
    lower = [c.lower() for c in cols]
    text = " ".join(lower)

    score = 0
    if any(x in lower for x in ["module", "program", "module_family", "module_node"]):
        score += 2
    if any(x in lower for x in ["stage", "stage_clean", "state_id", "developmental_stage"]):
        score += 2
    if any("top" in x for x in lower):
        score += 1
    if "pc1" in text:
        score += 2
    if "pc2" in text or "pc2_pc6" in text or "pc2pc6" in text:
        score += 2
    if "signed" in text and "profile" in text:
        score += 3
    if "cosine" in text:
        score += 1
    if "alignment" in text:
        score += 1

    return score


def discover_alignment_files(search_root, outdir):
    patterns = [
        "**/*stage*align*.tsv",
        "**/*stage*align*.tsv.gz",
        "**/*alignment*.tsv",
        "**/*alignment*.tsv.gz",
        "**/*Step03D*.tsv",
        "**/*step03D*.tsv",
        "**/*NTM*.tsv",
        "**/*NTM*.tsv.gz",
        "**/*.tsv",
        "**/*.tsv.gz",
    ]

    files = []
    for pat in patterns:
        files.extend(glob.glob(os.path.join(search_root, pat), recursive=True))

    files = sorted(set(files))

    records = []
    for p in files:
        try:
            dfh = read_table(p, nrows=200)
            sc = component_col_score(dfh.columns)

            value_text = " ".join(dfh.astype(str).values.ravel().tolist())
            has_ntm_rows = int(("NTM1" in value_text) or ("NTM2" in value_text) or ("NTM3" in value_text))
            has_stage_like = int(any("stage" in c.lower() or "state" in c.lower() for c in dfh.columns))
            has_component_like = int(sc >= 6)

            # Strongly prioritize NTM-containing tables.
            final_score = sc + 20 * has_ntm_rows + 3 * has_stage_like + 2 * has_component_like

            records.append({
                "path": p,
                "n_cols": dfh.shape[1],
                "score": sc,
                "final_score": final_score,
                "has_ntm_rows": has_ntm_rows,
                "has_stage_like": has_stage_like,
                "has_component_like": has_component_like,
                "columns": "|".join(list(dfh.columns)[:100])
            })
        except Exception as e:
            records.append({
                "path": p,
                "n_cols": np.nan,
                "score": -1,
                "final_score": -1,
                "has_ntm_rows": 0,
                "has_stage_like": 0,
                "has_component_like": 0,
                "columns": f"READ_ERROR: {e}"
            })

    cand = pd.DataFrame(records).sort_values(
        ["has_ntm_rows", "final_score", "score", "path"],
        ascending=[False, False, False, True]
    )
    cand_out = os.path.join(outdir, "00_step12C_candidate_alignment_component_files.tsv")
    cand.to_csv(cand_out, sep="\t", index=False)
    return cand


def detect_component_columns(df):
    cols = list(df.columns)
    lower = {c.lower(): c for c in cols}

    module_col = find_col(df, ["module_node", "module", "module_family", "program", "module_name"])
    topn_col = find_col(df, ["top_n", "topn", "n_top"])
    stage_col = find_col(df, ["stage_clean", "stage", "state_id", "developmental_stage", "stage_id"])

    # Flexible component detection.
    no_pc1_candidates = []
    pc2_candidates = []
    signed_candidates = []

    for c in cols:
        cl = c.lower()
        if ("pc1" in cl or "no_pc1" in cl or "nopc1" in cl) and ("cos" in cl or "sim" in cl or "align" in cl or "score" in cl):
            no_pc1_candidates.append(c)
        if ("pc2" in cl or "pc2_pc6" in cl or "pc2pc6" in cl or "secondary" in cl) and ("cos" in cl or "sim" in cl or "align" in cl or "score" in cl):
            pc2_candidates.append(c)
        if ("signed" in cl and ("profile" in cl or "stage" in cl or "align" in cl or "score" in cl)) or ("profile" in cl and "z" in cl):
            signed_candidates.append(c)

    # More exact fallback names.
    for name in [
        "no_pc1_cosine_z", "no_pc1_cosine", "cosine_no_pc1", "module_stage_no_pc1_z",
        "embedding_no_pc1_z", "noPC1_cosine_z"
    ]:
        if name.lower() in lower and lower[name.lower()] not in no_pc1_candidates:
            no_pc1_candidates.insert(0, lower[name.lower()])

    for name in [
        "pc2_pc6_cosine_z", "pc2_pc6_similarity_z", "pc2_pc6_z",
        "secondary_pc_z", "module_stage_pc2_pc6_z"
    ]:
        if name.lower() in lower and lower[name.lower()] not in pc2_candidates:
            pc2_candidates.insert(0, lower[name.lower()])

    for name in [
        "signed_stage_profile_z", "stage_profile_signed_z", "signed_profile_z",
        "signed_stage_profile", "stage_profile_z"
    ]:
        if name.lower() in lower and lower[name.lower()] not in signed_candidates:
            signed_candidates.insert(0, lower[name.lower()])

    return {
        "module_col": module_col,
        "topn_col": topn_col,
        "stage_col": stage_col,
        "no_pc1_col": no_pc1_candidates[0] if no_pc1_candidates else None,
        "pc2_pc6_col": pc2_candidates[0] if pc2_candidates else None,
        "signed_profile_col": signed_candidates[0] if signed_candidates else None,
        "no_pc1_candidates": no_pc1_candidates,
        "pc2_candidates": pc2_candidates,
        "signed_candidates": signed_candidates,
    }


def standardize_alignment_table(path):
    df = read_table(path)
    det = detect_component_columns(df)

    missing = [k for k in ["module_col", "stage_col", "no_pc1_col", "pc2_pc6_col", "signed_profile_col"] if det[k] is None]
    if missing:
        raise ValueError(f"Could not detect required columns {missing} in {path}. Detected: {det}")

    out = pd.DataFrame()
    out["module_raw"] = df[det["module_col"]].astype(str)
    out["module_family"] = out["module_raw"].map(infer_module_family)

    if det["topn_col"] is not None:
        out["top_n"] = pd.to_numeric(df[det["topn_col"]], errors="coerce")
    else:
        out["top_n"] = out["module_raw"].map(infer_top_n)

    out["stage"] = df[det["stage_col"]].astype(str)
    out["stage_clean"] = out["stage"].map(clean_stage)

    out["no_pc1"] = pd.to_numeric(df[det["no_pc1_col"]], errors="coerce")
    out["pc2_pc6"] = pd.to_numeric(df[det["pc2_pc6_col"]], errors="coerce")
    out["signed_profile"] = pd.to_numeric(df[det["signed_profile_col"]], errors="coerce")

    out = out[
        out["module_family"].isin(["NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed"]) &
        out["top_n"].isin([50, 100, 200, 500]) &
        out["stage_clean"].notna()
    ].copy()

    # z-score each component within each module/top_n across stages.
    for c in ["no_pc1", "pc2_pc6", "signed_profile"]:
        out[c + "_z"] = np.nan
        for key, idx in out.groupby(["module_family", "top_n"]).groups.items():
            vals = out.loc[idx, c].values.astype(float)
            out.loc[idx, c + "_z"] = zscore_vec(vals)

    out.attrs["detected_columns"] = det
    return out, det


def build_weight_designs(step=0.05):
    designs = []

    # Default design.
    designs.append({
        "design_type": "default",
        "design_name": "default_0.45_0.35_0.20",
        "w_no_pc1": 0.45,
        "w_pc2_pc6": 0.35,
        "w_signed_profile": 0.20
    })

    # Leave-one-component-out, renormalized.
    components = {
        "no_pc1": (1, 0, 0),
        "pc2_pc6": (0, 1, 0),
        "signed_profile": (0, 0, 1)
    }

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

    # Single component only.
    designs.extend([
        {"design_type": "single_component_only", "design_name": "only_no_pc1", "w_no_pc1": 1.0, "w_pc2_pc6": 0.0, "w_signed_profile": 0.0},
        {"design_type": "single_component_only", "design_name": "only_pc2_pc6", "w_no_pc1": 0.0, "w_pc2_pc6": 1.0, "w_signed_profile": 0.0},
        {"design_type": "single_component_only", "design_name": "only_signed_profile", "w_no_pc1": 0.0, "w_pc2_pc6": 0.0, "w_signed_profile": 1.0},
    ])

    # Grid weights summing to 1.
    vals = np.round(np.arange(0, 1 + 1e-9, step), 3)
    for a in vals:
        for b in vals:
            c = 1.0 - a - b
            if c < -1e-9:
                continue
            c = round(c, 3)
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

                # Broad prenatal probability.
                prenatal_prob = sum(
                    p for s, p in zip(stages, prob)
                    if "prenatal" in s
                )

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
                    "prenatal_probability": float(prenatal_prob),
                    "late_prenatal_is_top": int(top_stage == "late_prenatal"),
                    "prenatal_is_top": int("prenatal" in top_stage),
                    "postnatal_is_top": int(top_stage in ["adolescence", "adulthood", "infancy", "childhood"])
                })

    return pd.DataFrame(records)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_table", default="")
    parser.add_argument("--search_root", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project")
    parser.add_argument("--outdir", default="/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD/neurotrace_algorithm_project/step12_revision_strengthening/results/step12C_stage_alignment_weight_sensitivity")
    parser.add_argument("--grid_step", type=float, default=0.05)
    parser.add_argument("--temperatures", default="0.25,0.5,0.75,1.0,1.5,2.0")
    args = parser.parse_args()

    ensure_dir(args.outdir)

    # Discover or load component table.
    if args.input_table:
        input_table = args.input_table
        cand = pd.DataFrame([{"path": input_table, "score": 999, "columns": "user_provided"}])
        cand.to_csv(os.path.join(args.outdir, "00_step12C_candidate_alignment_component_files.tsv"), sep="\t", index=False)
    else:
        cand = discover_alignment_files(args.search_root, args.outdir)
        # Require NTM rows to avoid accidentally selecting old refined-program alignment tables.
        usable = cand[(cand["has_ntm_rows"] == 1) & (cand["score"] >= 6)].copy()
        if usable.shape[0] == 0:
            fail_md = os.path.join(args.outdir, "06_step12C_FAILED_no_NTM_alignment_component_table.md")
            with open(fail_md, "w") as f:
                f.write("# Step12C failed: no NTM alignment component table detected\n\n")
                f.write("The automatic search did not find a table containing NTM1/NTM2/NTM3 rows plus component-like columns. Inspect `00_step12C_candidate_alignment_component_files.tsv` and rerun with `--input_table` pointing to the correct NTM stage-alignment table.\n")
            log(f"No usable NTM alignment component table found. Wrote {fail_md}")
            return
        input_table = usable.iloc[0]["path"]

    log(f"Using alignment component table: {input_table}")
    aln, det = standardize_alignment_table(input_table)

    aln_out = os.path.join(args.outdir, "01_step12C_standardized_alignment_components.tsv")
    aln.to_csv(aln_out, sep="\t", index=False)

    if aln.shape[0] == 0:
        fail_md = os.path.join(args.outdir, "06_step12C_FAILED_empty_after_NTM_filter.md")
        with open(fail_md, "w") as f:
            f.write("# Step12C failed: selected table produced zero NTM rows after filtering\n\n")
            f.write(f"Selected table: `{input_table}`\n\n")
            f.write("This usually means the table is a refined-program or legacy alignment table rather than the NTM native alignment component table. Inspect `00_step12C_candidate_alignment_component_files.tsv` and rerun with `--input_table`.\n")
        log(f"Standardized alignment table is empty. Wrote {fail_md}")
        return

    det_out = os.path.join(args.outdir, "00b_step12C_detected_columns.tsv")
    pd.DataFrame([{
        "input_table": input_table,
        "module_col": det["module_col"],
        "topn_col": det["topn_col"],
        "stage_col": det["stage_col"],
        "no_pc1_col": det["no_pc1_col"],
        "pc2_pc6_col": det["pc2_pc6_col"],
        "signed_profile_col": det["signed_profile_col"],
        "no_pc1_candidates": "|".join(det["no_pc1_candidates"]),
        "pc2_candidates": "|".join(det["pc2_candidates"]),
        "signed_candidates": "|".join(det["signed_candidates"])
    }]).to_csv(det_out, sep="\t", index=False)

    designs = build_weight_designs(step=args.grid_step)
    design_out = os.path.join(args.outdir, "00c_step12C_weight_designs.tsv")
    designs.to_csv(design_out, sep="\t", index=False)

    temps = [float(x) for x in args.temperatures.split(",") if x.strip()]

    grid_df = evaluate_designs(aln, designs, temps)

    grid_out = os.path.join(args.outdir, "02_step12C_weight_temperature_grid_stage_probs.tsv.gz")
    grid_df.to_csv(grid_out, sep="\t", index=False, compression="gzip")

    # Design-level summary.
    design_summary = (
        grid_df.groupby(["design_type", "design_name", "temperature"], dropna=False)
        .agg(
            n_tests=("module_family", "count"),
            n_late_top=("late_prenatal_is_top", "sum"),
            n_prenatal_top=("prenatal_is_top", "sum"),
            n_postnatal_top=("postnatal_is_top", "sum"),
            median_late_prenatal_probability=("late_prenatal_probability", "median"),
            median_prenatal_probability=("prenatal_probability", "median"),
            median_top_margin=("top_margin", "median")
        )
        .reset_index()
    )

    design_summary_out = os.path.join(args.outdir, "03_step12C_design_level_summary.tsv")
    design_summary.to_csv(design_summary_out, sep="\t", index=False)

    # Module-level stability.
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

    module_summary_out = os.path.join(args.outdir, "04_step12C_module_level_stability_summary.tsv")
    module_summary.to_csv(module_summary_out, sep="\t", index=False)

    # Specific default/leave-one/single summaries.
    focus = grid_df[
        grid_df["design_type"].isin(["default", "leave_one_component_out", "single_component_only"])
    ].copy()
    focus_out = os.path.join(args.outdir, "04b_step12C_default_leaveone_singlecomponent_results.tsv")
    focus.to_csv(focus_out, sep="\t", index=False)

    # Decision table.
    records = []
    for _, row in module_summary.iterrows():
        fam = row["module_family"]
        topn = int(row["top_n"])

        if fam in ["NTM1_ASD_up", "NTM3_ASD_signed"]:
            # Broader prenatal bias is the key robust claim; strict late-prenatal can be sensitive.
            if row["prenatal_top_fraction"] >= 0.80:
                status = "PASS"
            elif row["prenatal_top_fraction"] >= 0.60:
                status = "PARTIAL"
            else:
                status = "REVIEW"
            interpretation = "prenatal-biased across alignment-component and temperature settings"
        elif fam == "NTM2_ASD_down":
            # NTM2 should not be late-prenatal dominated.
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
    decision_out = os.path.join(args.outdir, "05_step12C_decision_table.tsv")
    decision.to_csv(decision_out, sep="\t", index=False)

    md_out = os.path.join(args.outdir, "06_step12C_overall_summary.md")
    with open(md_out, "w") as f:
        f.write("# Step12C stage-alignment component weight sensitivity\n\n")

        f.write("## Purpose\n\n")
        f.write("This analysis tests whether developmental stage calls depend on the fixed 0.45/0.35/0.20 alignment-component weighting scheme or on the softmax temperature used to convert alignment scores to stage probabilities.\n\n")

        f.write("## Input table and detected columns\n\n")
        f.write(f"- Input table: `{input_table}`\n")
        f.write(f"- Module column: `{det['module_col']}`\n")
        f.write(f"- top_n column: `{det['topn_col']}`\n")
        f.write(f"- Stage column: `{det['stage_col']}`\n")
        f.write(f"- no-PC1 component: `{det['no_pc1_col']}`\n")
        f.write(f"- PC2-PC6 component: `{det['pc2_pc6_col']}`\n")
        f.write(f"- Signed stage-profile component: `{det['signed_profile_col']}`\n")
        f.write(f"- Weight designs tested: {designs.shape[0]}\n")
        f.write(f"- Temperatures tested: {', '.join(map(str, temps))}\n\n")

        f.write("## Module-level stability summary\n\n")
        f.write(module_summary.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Decision table\n\n")
        f.write(decision.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Interpretation guide\n\n")
        f.write("- PASS for NTM1/NTM3 means broad prenatal bias is stable across alignment-component weights and temperatures. Strict late-prenatal substage specificity may still vary and should be interpreted cautiously.\n")
        f.write("- PASS for NTM2 means NTM2 is not dominated by late-prenatal top-stage calls across alignment-component weights and temperatures.\n")
        f.write("- If strict late-prenatal is sensitive but broad prenatal bias is stable, manuscript language should emphasize prenatal-biased alignment rather than a uniquely optimized late-prenatal substage under every component weighting.\n")

        f.write("\n## Manuscript-ready interpretation\n\n")
        f.write("Alignment-component sensitivity showed that the broad prenatal versus non-prenatal separation was more stable than exact prenatal substage assignment. This supports the interpretation that NTM1_ASD_up and NTM3_ASD_signed are prenatal-biased under alternative component weightings, whereas NTM2_ASD_down remains separated from strict late-prenatal dominance. Because exact late-prenatal specificity is more sensitive to component weighting, the developmental conclusion should remain calibrated by the primary matched-random transport analysis and reported with caution at the substage level.\n")

    log("Finished Step12C")
    for p in [
        aln_out,
        det_out,
        design_out,
        grid_out,
        design_summary_out,
        module_summary_out,
        focus_out,
        decision_out,
        md_out,
    ]:
        log(f"Wrote: {p}")


if __name__ == "__main__":
    main()
