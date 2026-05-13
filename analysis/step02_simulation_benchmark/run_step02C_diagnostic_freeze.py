#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
STEP = BASE / "neurotrace_algorithm_project" / "step02_simulation_benchmark"
OUT = STEP / "results"
FIG = STEP / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SUMMARY = OUT / "13_step02B_summary_by_condition_method.tsv"
GLOBAL = OUT / "16_step02B_global_method_summary.tsv"
ABLATION = OUT / "15_step02B_ablation_component_loss.tsv"
DELTAS = OUT / "14_step02B_paired_deltas_full_vs_baselines.tsv"

PRIMARY_ALIGNMENT_METRICS = [
    "alignment_top1_mean",
    "alignment_top2_mean",
    "true_transport_mass_mean"
]

DECOY_METRICS = [
    "decoy_transport_mass_mean"
]

GENE_METRICS = [
    "gene_auprc_mean",
    "precision_at_50_mean",
    "precision_at_100_mean"
]

NEUROTRACE_METHODS = [
    "neurotrace_full",
    "neurotrace_no_embedding",
    "neurotrace_no_graph",
    "neurotrace_no_risk",
    "neurotrace_no_invariance"
]

BASELINES = [
    "mean_signature",
    "embedding_nn",
    "simple_ot",
    "graph_ot"
]

MAIN_CONDITIONS = [
    "weak_signal",
    "low_overlap",
    "composition_confounded",
    "high_dropout",
    "false_prior_stress"
]

STRESS_ONLY_CONDITIONS = [
    "combined_hard"
]


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def rank_within_condition(df, metric, higher_is_better=True):
    ascending = not higher_is_better
    out = df[["condition", "method", metric]].copy()
    out[f"{metric}_rank"] = out.groupby("condition")[metric].rank(
        method="min",
        ascending=ascending
    )
    return out


def add_ranks(summary):
    parts = []

    for metric in PRIMARY_ALIGNMENT_METRICS:
        parts.append(rank_within_condition(summary, metric, higher_is_better=True))

    for metric in GENE_METRICS:
        parts.append(rank_within_condition(summary, metric, higher_is_better=True))

    for metric in DECOY_METRICS:
        parts.append(rank_within_condition(summary, metric, higher_is_better=False))

    out = summary.copy()

    for p in parts:
        rank_cols = [c for c in p.columns if c.endswith("_rank")]
        out = out.merge(
            p[["condition", "method"] + rank_cols],
            on=["condition", "method"],
            how="left"
        )

    align_rank_cols = [f"{m}_rank" for m in PRIMARY_ALIGNMENT_METRICS]
    gene_rank_cols = [f"{m}_rank" for m in GENE_METRICS]
    decoy_rank_cols = [f"{m}_rank" for m in DECOY_METRICS]

    out["mean_alignment_rank"] = out[align_rank_cols].mean(axis=1)
    out["mean_gene_rank"] = out[gene_rank_cols].mean(axis=1)
    out["mean_decoy_rank"] = out[decoy_rank_cols].mean(axis=1)
    out["mean_overall_rank"] = out[
        ["mean_alignment_rank", "mean_gene_rank", "mean_decoy_rank"]
    ].mean(axis=1)

    out = out.sort_values(["condition", "mean_alignment_rank", "mean_gene_rank", "method"])
    return out


def condition_winners(ranked):
    rows = []

    for cond, sub in ranked.groupby("condition"):
        sub = sub.copy()

        align_winner = sub.sort_values(
            ["mean_alignment_rank", "alignment_top1_mean", "true_transport_mass_mean"],
            ascending=[True, False, False]
        ).iloc[0]

        gene_winner = sub.sort_values(
            ["mean_gene_rank", "gene_auprc_mean", "precision_at_50_mean"],
            ascending=[True, False, False]
        ).iloc[0]

        decoy_winner = sub.sort_values(
            ["mean_decoy_rank", "decoy_transport_mass_mean"],
            ascending=[True, True]
        ).iloc[0]

        overall_winner = sub.sort_values(
            ["mean_overall_rank", "mean_alignment_rank", "mean_gene_rank"],
            ascending=[True, True, True]
        ).iloc[0]

        full = sub[sub["method"] == "neurotrace_full"].iloc[0]
        best_baseline = sub[sub["method"].isin(BASELINES)].sort_values(
            ["mean_alignment_rank", "gene_auprc_mean"],
            ascending=[True, False]
        ).iloc[0]

        rows.append({
            "condition": cond,
            "alignment_winner": align_winner["method"],
            "gene_winner": gene_winner["method"],
            "decoy_lowest_method": decoy_winner["method"],
            "overall_winner": overall_winner["method"],
            "neurotrace_full_alignment_top1": full["alignment_top1_mean"],
            "neurotrace_full_true_transport_mass": full["true_transport_mass_mean"],
            "neurotrace_full_decoy_transport_mass": full["decoy_transport_mass_mean"],
            "neurotrace_full_gene_auprc": full["gene_auprc_mean"],
            "best_baseline_method_by_alignment": best_baseline["method"],
            "best_baseline_alignment_top1": best_baseline["alignment_top1_mean"],
            "delta_full_minus_best_baseline_alignment_top1": full["alignment_top1_mean"] - best_baseline["alignment_top1_mean"],
            "delta_full_minus_best_baseline_true_transport_mass": full["true_transport_mass_mean"] - best_baseline["true_transport_mass_mean"],
            "delta_full_minus_best_baseline_gene_auprc": full["gene_auprc_mean"] - best_baseline["gene_auprc_mean"],
            "condition_use_recommendation": "main_benchmark" if cond in MAIN_CONDITIONS else "stress_test_only"
        })

    return pd.DataFrame(rows)


def global_objective_summary(global_df):
    rows = []

    g = global_df.copy()

    align_cols = ["alignment_top1_mean", "true_transport_mass_mean"]
    gene_cols = ["gene_auprc_mean", "precision_at_50_mean"]
    decoy_cols = ["decoy_transport_mass_mean"]

    g["alignment_objective_rank"] = (
        g["alignment_top1_mean"].rank(method="min", ascending=False) +
        g["true_transport_mass_mean"].rank(method="min", ascending=False)
    ) / 2

    g["gene_objective_rank"] = (
        g["gene_auprc_mean"].rank(method="min", ascending=False) +
        g["precision_at_50_mean"].rank(method="min", ascending=False)
    ) / 2

    g["decoy_objective_rank"] = g["decoy_transport_mass_mean"].rank(
        method="min",
        ascending=True
    )

    g["overall_objective_rank"] = g[
        ["alignment_objective_rank", "gene_objective_rank", "decoy_objective_rank"]
    ].mean(axis=1)

    g = g.sort_values(["alignment_objective_rank", "gene_objective_rank", "method"])
    return g


def make_decision_table(condition_df, global_obj, ablation):
    full_global = global_obj[global_obj["method"] == "neurotrace_full"].iloc[0]

    n_main = condition_df["condition_use_recommendation"].eq("main_benchmark").sum()
    n_full_alignment_wins_main = (
        condition_df[
            condition_df["condition_use_recommendation"] == "main_benchmark"
        ]["alignment_winner"].eq("neurotrace_full").sum()
    )

    n_full_positive_vs_baseline_main = (
        condition_df[
            condition_df["condition_use_recommendation"] == "main_benchmark"
        ]["delta_full_minus_best_baseline_alignment_top1"].gt(0).sum()
    )

    # ablation evidence
    ab = ablation.copy()
    if len(ab) > 0:
        ab_summary = ab.groupby("method").agg(
            mean_loss_alignment=("loss_alignment_top1", "mean"),
            mean_loss_gene_auprc=("loss_gene_auprc", "mean"),
            mean_loss_transport=("loss_true_transport_mass", "mean"),
            mean_decoy_gain_when_removed=("gain_decoy_transport_mass_when_removed", "mean")
        ).reset_index()
    else:
        ab_summary = pd.DataFrame()

    rows = [
        {
            "decision_item": "Step02B developmental alignment benchmark",
            "status": "PASS" if full_global["alignment_objective_rank"] <= 1.5 else "CAUTION",
            "evidence": f"neurotrace_full global alignment_top1={full_global['alignment_top1_mean']:.3f}, true_transport_mass={full_global['true_transport_mass_mean']:.3f}",
            "recommendation": "Use as main simulation benchmark endpoint."
        },
        {
            "decision_item": "Step02B gene prioritization benchmark",
            "status": "CAUTION" if full_global["gene_objective_rank"] > 1.5 else "PASS",
            "evidence": f"neurotrace_full global gene_AUPRC={full_global['gene_auprc_mean']:.3f}; gene_objective_rank={full_global['gene_objective_rank']:.2f}",
            "recommendation": "Do not claim full NeuroTRACE is best for gene prioritization; present gene ranking as secondary or tune separately."
        },
        {
            "decision_item": "Condition-level main benchmark stability",
            "status": "PASS" if n_full_positive_vs_baseline_main >= 4 else "CAUTION",
            "evidence": f"full NeuroTRACE beats best baseline by alignment top1 in {n_full_positive_vs_baseline_main}/{n_main} main conditions",
            "recommendation": "Use main conditions separately; keep combined_hard as stress-only."
        },
        {
            "decision_item": "Ablation interpretability",
            "status": "PARTIAL",
            "evidence": "Component removal affects alignment and transport, but no_invariance improves gene AUPRC.",
            "recommendation": "Separate alignment objective from gene-prioritization objective in the manuscript."
        }
    ]

    return pd.DataFrame(rows), ab_summary


def make_plots(ranked, condition_df, global_obj):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Main condition alignment top1
        plot_df = ranked[ranked["condition"].isin(MAIN_CONDITIONS)].copy()
        pivot = plot_df.pivot(index="condition", columns="method", values="alignment_top1_mean")
        pivot = pivot.loc[MAIN_CONDITIONS]

        ax = pivot.plot(kind="bar", figsize=(14, 5), width=0.85)
        ax.set_ylabel("Top-1 alignment accuracy")
        ax.set_xlabel("Main hard simulation condition")
        ax.set_title("Step02C: developmental alignment benchmark")
        ax.legend(title="Method", bbox_to_anchor=(1.02, 1.0), loc="upper left")
        plt.tight_layout()
        plt.savefig(FIG / "13_step02C_main_conditions_alignment_top1.pdf")
        plt.close()

        # Gene AUPRC
        pivot2 = plot_df.pivot(index="condition", columns="method", values="gene_auprc_mean")
        pivot2 = pivot2.loc[MAIN_CONDITIONS]
        ax = pivot2.plot(kind="bar", figsize=(14, 5), width=0.85)
        ax.set_ylabel("Gene prioritization AUPRC")
        ax.set_xlabel("Main hard simulation condition")
        ax.set_title("Step02C: gene prioritization is a distinct objective")
        ax.legend(title="Method", bbox_to_anchor=(1.02, 1.0), loc="upper left")
        plt.tight_layout()
        plt.savefig(FIG / "14_step02C_main_conditions_gene_AUPRC.pdf")
        plt.close()

        # Objective ranks
        g = global_obj.sort_values("overall_objective_rank")
        ax = g.set_index("method")[[
            "alignment_objective_rank",
            "gene_objective_rank",
            "decoy_objective_rank"
        ]].plot(kind="bar", figsize=(12, 5), width=0.8)
        ax.set_ylabel("Rank; lower is better")
        ax.set_xlabel("Method")
        ax.set_title("Step02C: objective-specific method ranks")
        plt.tight_layout()
        plt.savefig(FIG / "15_step02C_objective_specific_ranks.pdf")
        plt.close()

    except Exception as e:
        with open(OUT / "Step02C_plotting_failed.txt", "w") as f:
            f.write(str(e) + "\n")


def main():
    print(f"[{now()}] Step02C diagnostic freeze started")

    for f in [SUMMARY, GLOBAL, ABLATION, DELTAS]:
        if not f.exists():
            raise FileNotFoundError(f"Missing required Step02B file: {f}")

    summary = pd.read_csv(SUMMARY, sep="\t")
    global_df = pd.read_csv(GLOBAL, sep="\t")
    ablation = pd.read_csv(ABLATION, sep="\t")

    ranked = add_ranks(summary)
    ranked.to_csv(OUT / "18_step02C_condition_method_metric_ranks.tsv", sep="\t", index=False)

    cond_winners = condition_winners(ranked)
    cond_winners.to_csv(OUT / "19_step02C_condition_objective_winners.tsv", sep="\t", index=False)

    global_obj = global_objective_summary(global_df)
    global_obj.to_csv(OUT / "20_step02C_global_objective_specific_ranks.tsv", sep="\t", index=False)

    decision, ab_summary = make_decision_table(cond_winners, global_obj, ablation)
    decision.to_csv(OUT / "21_step02C_freeze_decision_table.tsv", sep="\t", index=False)
    ab_summary.to_csv(OUT / "22_step02C_ablation_summary_by_component.tsv", sep="\t", index=False)

    make_plots(ranked, cond_winners, global_obj)

    with open(OUT / "23_step02C_diagnostic_freeze_summary.md", "w") as f:
        f.write("# NeuroTRACE Step02C diagnostic freeze summary\n\n")
        f.write(f"Generated: {now()}\n\n")

        f.write("## Core conclusion\n")
        f.write("Step02B supports NeuroTRACE-full primarily as a developmental alignment algorithm. ")
        f.write("It should not yet be claimed as the strongest gene-prioritization model because no_invariance and some baseline methods show higher gene-level AUPRC.\n\n")

        f.write("## Freeze decision table\n\n")
        f.write(decision.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Global objective-specific ranks\n\n")
        f.write(global_obj.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Condition-level winners\n\n")
        f.write(cond_winners.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Recommended manuscript wording\n")
        f.write("Use Step02B/Step02C to support the claim that NeuroTRACE improves disease-to-developmental state alignment under weak signal, low overlap, dropout, and composition-confounded settings. ")
        f.write("State explicitly that gene prioritization is a secondary output and may require objective-specific calibration.\n")

    print(f"[{now()}] Step02C done")
    print(f"[{now()}] Results written to: {OUT}")
    print(f"[{now()}] Figures written to: {FIG}")


if __name__ == "__main__":
    main()
