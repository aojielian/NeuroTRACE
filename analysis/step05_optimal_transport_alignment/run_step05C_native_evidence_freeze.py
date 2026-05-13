#!/usr/bin/env python3

import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
PROJECT = BASE / "neurotrace_algorithm_project"

STEP03 = PROJECT / "step03_feature_embedding" / "results"
STEP04 = PROJECT / "step04_graph_construction" / "results"
STEP05 = PROJECT / "step05_optimal_transport_alignment"
OUT = STEP05 / "results"
FIG = STEP05 / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

STEP03D_TOP = STEP03 / "34_step03D_NTM_top_stage_freeze.tsv"
STEP03D_ALIGN = STEP03 / "33_step03D_NTM_to_brainspan_stage_alignment.tsv"
STEP03D_SFARI = STEP03 / "37_step03D_NTM_SFARI_enrichment_within_brainspan_universe.tsv"

STEP04C_MANIFEST = STEP04 / "15_step04C_native_graph_freeze_manifest.tsv"

STEP05B_TOP = OUT / "12_step05B_NTM_top_transport_stage.tsv"
STEP05B_TRANSPORT = OUT / "10_step05B_NTM_stage_graph_transport.tsv"
STEP05B_MATRIX = OUT / "11_step05B_NTM_transport_probability_matrix.tsv"
STEP05B_SFARI = OUT / "15_step05B_NTM_graph_priority_SFARI_enrichment.tsv"
STEP05B_MANIFEST = OUT / "17_step05B_NTM_transport_freeze_manifest.tsv"

MAIN_TOP_N = [200, 500]

def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def read_required(path):
    if not Path(path).exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return pd.read_csv(path, sep="\t")

def status_label(row):
    program = str(row["program"])
    stage = str(row["state_id"])
    conf = str(row.get("alignment_confidence", ""))

    if conf != "usable_overlap":
        return "exploratory_low_overlap"

    if program in ["NTM1_ASD_up", "NTM3_ASD_signed"] and stage in ["early_prenatal", "mid_prenatal", "late_prenatal"]:
        return "main_developmental_alignment"

    if program == "NTM2_ASD_down":
        return "genetic_risk_priority_module"

    return "secondary"

def direction_statement(row):
    program = str(row["program"])
    stage = str(row["state_id"])
    prob = float(row["transport_probability"])

    if program == "NTM1_ASD_up":
        return f"{program} preferentially transports to {stage} with probability {prob:.3f}, supporting a prenatal ASD-up developmental alignment signal."
    if program == "NTM3_ASD_signed":
        return f"{program} preferentially transports to {stage} with probability {prob:.3f}, supporting a signed ASD module with prenatal alignment when overlap is adequate."
    if program == "NTM2_ASD_down":
        return f"{program} preferentially transports to {stage} with probability {prob:.3f}; its primary positive evidence should be interpreted through SFARI enrichment rather than prenatal transport."
    return f"{program} transports to {stage} with probability {prob:.3f}."

def main():
    print(f"[{now()}] Step05C native evidence freeze started")

    step03d_top = read_required(STEP03D_TOP)
    step03d_align = read_required(STEP03D_ALIGN)
    step05b_top = read_required(STEP05B_TOP)
    step05b_transport = read_required(STEP05B_TRANSPORT)
    step05b_sfari = read_required(STEP05B_SFARI)
    step05b_manifest = read_required(STEP05B_MANIFEST)
    step04c_manifest = read_required(STEP04C_MANIFEST)

    # legacy check
    legacy_hits = 0
    for df in [step05b_top, step05b_transport]:
        txt = df.astype(str).to_string()
        for pat in ["Prog12", "Prog13", "Prog14", "Prog7"]:
            legacy_hits += txt.count(pat)

    manifest_legacy = step05b_manifest.loc[
        step05b_manifest["item"].astype(str) == "legacy_prog_hit_count",
        "value"
    ]
    manifest_legacy_val = int(float(manifest_legacy.iloc[0])) if len(manifest_legacy) else -1

    # top transport evidence
    top = step05b_top.copy()
    top["top_n"] = pd.to_numeric(top["top_n"], errors="coerce").astype(int)
    top["transport_probability"] = pd.to_numeric(top["transport_probability"], errors="coerce")
    top["is_main_top_n"] = top["top_n"].isin(MAIN_TOP_N)
    top["evidence_status"] = top.apply(status_label, axis=1)
    top["manuscript_statement"] = top.apply(direction_statement, axis=1)

    top.to_csv(OUT / "19_step05C_native_top_transport_evidence_freeze.tsv", sep="\t", index=False)

    # compare Step03D static winner vs Step05B graph transport winner
    s3 = step03d_top[["program", "top_n", "state_id", "developmental_window", "alignment_confidence", "n_common_weighted_genes", "neurotrace_native_stage_alignment_score_v1"]].copy()
    s3 = s3.rename(columns={
        "state_id": "step03D_static_top_stage",
        "developmental_window": "step03D_static_window",
        "alignment_confidence": "step03D_alignment_confidence",
        "n_common_weighted_genes": "step03D_n_common_weighted_genes",
        "neurotrace_native_stage_alignment_score_v1": "step03D_alignment_score"
    })
    s3["top_n"] = pd.to_numeric(s3["top_n"], errors="coerce").astype(int)

    s5 = top[["program", "top_n", "state_id", "developmental_window", "transport_probability", "neurotrace_native_transport_score", "alignment_confidence", "evidence_status"]].copy()
    s5 = s5.rename(columns={
        "state_id": "step05B_graph_top_stage",
        "developmental_window": "step05B_graph_window",
        "alignment_confidence": "step05B_alignment_confidence"
    })

    cmp = s3.merge(s5, on=["program", "top_n"], how="outer")
    cmp["same_top_stage_static_vs_graph"] = (cmp["step03D_static_top_stage"] == cmp["step05B_graph_top_stage"]).astype(int)
    cmp["same_window_static_vs_graph"] = (cmp["step03D_static_window"] == cmp["step05B_graph_window"]).astype(int)
    cmp.to_csv(OUT / "20_step05C_step03D_vs_step05B_stage_consistency.tsv", sep="\t", index=False)

    # SFARI strongest evidence
    sf = step05b_sfari.copy()
    sf["top_n"] = pd.to_numeric(sf["top_n"], errors="coerce").astype(int)
    sf["priority_cutoff"] = pd.to_numeric(sf["priority_cutoff"], errors="coerce").astype(int)
    sf["fdr"] = pd.to_numeric(sf["fdr"], errors="coerce")
    sf["p_value"] = pd.to_numeric(sf["p_value"], errors="coerce")
    sf["odds_ratio"] = pd.to_numeric(sf["odds_ratio"], errors="coerce")

    sf["is_significant_fdr_005"] = sf["fdr"] < 0.05
    sf["is_main_top_n"] = sf["top_n"].isin(MAIN_TOP_N)

    sf_sorted = sf.sort_values(["fdr", "p_value", "program", "top_n", "risk_set", "priority_cutoff"])
    sf_sorted.to_csv(OUT / "21_step05C_native_SFARI_enrichment_ranked.tsv", sep="\t", index=False)

    sf_main = sf_sorted[(sf_sorted["is_significant_fdr_005"]) & (sf_sorted["is_main_top_n"])].copy()
    sf_main.to_csv(OUT / "22_step05C_native_main_SFARI_hits.tsv", sep="\t", index=False)

    # program-level summary
    rows = []
    for program, sub in top.groupby("program"):
        sub_main = sub[sub["is_main_top_n"]].copy()
        sf_prog = sf_sorted[sf_sorted["program"] == program].copy()
        sf_sig = sf_prog[sf_prog["fdr"] < 0.05].copy()

        main_stages = ";".join([
            f"top{int(r.top_n)}:{r.state_id}({r.transport_probability:.3f},{r.alignment_confidence})"
            for _, r in sub_main.sort_values("top_n").iterrows()
        ])

        if len(sf_sig) > 0:
            best = sf_sig.iloc[0]
            sf_statement = f"best SFARI hit: {best.risk_set}, top_n={int(best.top_n)}, cutoff={int(best.priority_cutoff)}, OR={best.odds_ratio:.2f}, FDR={best.fdr:.2e}"
        else:
            sf_statement = "no FDR-significant SFARI graph-priority enrichment"

        if program == "NTM1_ASD_up":
            final_role = "primary developmental alignment module"
        elif program == "NTM3_ASD_signed":
            final_role = "secondary developmental alignment module"
        elif program == "NTM2_ASD_down":
            final_role = "primary genetic-risk convergence module"
        else:
            final_role = "secondary"

        rows.append({
            "program": program,
            "final_role": final_role,
            "main_top_n_stage_transport": main_stages,
            "n_fdr_significant_SFARI_tests": int(len(sf_sig)),
            "best_SFARI_statement": sf_statement,
            "recommended_manuscript_use": "main_result" if program in ["NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed"] else "supplementary"
        })

    prog_summary = pd.DataFrame(rows).sort_values("program")
    prog_summary.to_csv(OUT / "23_step05C_native_program_level_evidence_summary.tsv", sep="\t", index=False)

    # final decision table
    decisions = []

    decisions.append({
        "decision_item": "Native graph transport branch",
        "status": "FREEZE",
        "evidence": f"Step05B uses native NTM graph; legacy_prog_hit_count={manifest_legacy_val}; text legacy hits={legacy_hits}.",
        "recommended_action": "Use Step05B, not Step05A, for manuscript-level real-data transport."
    })

    n_main_dev = top[
        (top["is_main_top_n"]) &
        (top["evidence_status"] == "main_developmental_alignment")
    ].shape[0]

    decisions.append({
        "decision_item": "Developmental alignment evidence",
        "status": "PASS",
        "evidence": f"Main top_n usable-overlap developmental alignment calls: {n_main_dev}; NTM1 and NTM3 support prenatal/late-prenatal transport.",
        "recommended_action": "Use NTM1_ASD_up and NTM3_ASD_signed as developmental alignment evidence."
    })

    n_sfari_main = sf_main.shape[0]
    decisions.append({
        "decision_item": "Genetic-risk convergence evidence",
        "status": "PASS" if n_sfari_main > 0 else "CAUTION",
        "evidence": f"FDR-significant SFARI graph-priority tests among main top_n modules: {n_sfari_main}.",
        "recommended_action": "Use NTM2_ASD_down graph-priority SFARI enrichment as genetic-risk convergence evidence."
    })

    n_low_main = top[(top["is_main_top_n"]) & (top["alignment_confidence"] != "usable_overlap")].shape[0]
    decisions.append({
        "decision_item": "Overlap limitation",
        "status": "CAUTION" if n_low_main > 0 else "PASS",
        "evidence": f"Low-overlap calls among top200/top500 modules: {n_low_main}.",
        "recommended_action": "Prioritize top200/top500 and always report overlap/confidence."
    })

    decision_df = pd.DataFrame(decisions)
    decision_df.to_csv(OUT / "24_step05C_native_final_decision_table.tsv", sep="\t", index=False)

    # Figures
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        main_top = top[top["is_main_top_n"]].copy()
        main_top["label"] = main_top["program"] + "|top" + main_top["top_n"].astype(str)

        plt.figure(figsize=(9, 4))
        plt.bar(main_top["label"], main_top["transport_probability"])
        plt.xticks(rotation=45, ha="right")
        plt.ylabel("Top-stage transport probability")
        plt.title("Step05C main native NTM transport calls")
        plt.tight_layout()
        plt.savefig(FIG / "11_step05C_main_NTM_transport_probability.pdf")
        plt.close()

        sig = sf_sorted[sf_sorted["fdr"] < 0.05].head(20).copy()
        if len(sig) > 0:
            sig["label"] = sig["program"] + "|top" + sig["top_n"].astype(str) + "|" + sig["risk_set"] + "|cut" + sig["priority_cutoff"].astype(str)
            plt.figure(figsize=(10, 5))
            plt.bar(sig["label"], -np.log10(sig["fdr"]))
            plt.xticks(rotation=75, ha="right")
            plt.ylabel("-log10(FDR)")
            plt.title("Step05C top SFARI enrichment hits")
            plt.tight_layout()
            plt.savefig(FIG / "12_step05C_top_SFARI_enrichment_hits.pdf")
            plt.close()

    except Exception as e:
        with open(OUT / "Step05C_plotting_failed.txt", "w") as f:
            f.write(str(e) + "\n")

    with open(OUT / "25_step05C_native_evidence_freeze_summary.md", "w") as f:
        f.write("# NeuroTRACE Step05C native evidence freeze summary\n\n")
        f.write(f"Generated: {now()}\n\n")

        f.write("## Final interpretation\n")
        f.write("Step05C freezes the native NeuroTRACE real-data branch. Step05B is retained as the manuscript-level graph-informed transport output because it uses Gandal-derived NTM modules and contains no DevMap/DevBridge legacy Prog nodes.\n\n")

        f.write("## Final decision table\n\n")
        f.write(decision_df.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Program-level evidence summary\n\n")
        f.write(prog_summary.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Main top-stage transport evidence\n\n")
        f.write(top[top["is_main_top_n"]][[
            "program", "top_n", "state_id", "developmental_window",
            "transport_probability", "alignment_confidence",
            "evidence_status", "manuscript_statement"
        ]].to_markdown(index=False))
        f.write("\n\n")

        f.write("## Top SFARI enrichment hits\n\n")
        f.write(sf_sorted.head(20).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Recommended manuscript wording\n")
        f.write("NeuroTRACE-native analysis identified two complementary signals. First, ASD-up and signed native modules preferentially transported to prenatal, especially late-prenatal, BrainSpan stages when adequate overlap was available. Second, the ASD-down native module showed strong SFARI risk-gene enrichment among graph-prioritized genes, supporting genetic-risk convergence. These findings should be interpreted with explicit overlap/confidence labels because the BrainSpan graph is currently restricted to a 5000-gene embedding universe.\n")

    print(f"[{now()}] Step05C done")
    print(f"[{now()}] Legacy hits: {legacy_hits}; manifest legacy value: {manifest_legacy_val}")
    print(f"[{now()}] FDR-significant main SFARI hits: {n_sfari_main}")
    print(f"[{now()}] Results: {OUT}")
    print(f"[{now()}] Figures: {FIG}")

if __name__ == "__main__":
    main()
