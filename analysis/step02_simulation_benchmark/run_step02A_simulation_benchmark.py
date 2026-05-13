#!/usr/bin/env python3

import os
import math
import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
STEP = BASE / "neurotrace_algorithm_project" / "step02_simulation_benchmark"
OUT = STEP / "results"
FIG = STEP / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SEED = 20260508
N_REP = 120
N_GENES = 1500
N_STATES = 8
N_COMMUNITIES = 30
N_COHORTS = 3
EMBED_DIM = 24

CONDITIONS = {
    "baseline": {
        "batch_sd": 0.25,
        "composition_sd": 0.20,
        "dropout_rate": 0.05,
        "false_prior_rate": 0.04,
        "cohort_noise_sd": 0.45
    },
    "high_batch": {
        "batch_sd": 0.85,
        "composition_sd": 0.20,
        "dropout_rate": 0.05,
        "false_prior_rate": 0.04,
        "cohort_noise_sd": 0.60
    },
    "composition_shift": {
        "batch_sd": 0.25,
        "composition_sd": 0.80,
        "dropout_rate": 0.05,
        "false_prior_rate": 0.04,
        "cohort_noise_sd": 0.55
    },
    "high_dropout": {
        "batch_sd": 0.25,
        "composition_sd": 0.20,
        "dropout_rate": 0.30,
        "false_prior_rate": 0.04,
        "cohort_noise_sd": 0.55
    },
    "false_prior": {
        "batch_sd": 0.25,
        "composition_sd": 0.20,
        "dropout_rate": 0.05,
        "false_prior_rate": 0.22,
        "cohort_noise_sd": 0.55
    },
    "combined_stress": {
        "batch_sd": 0.70,
        "composition_sd": 0.70,
        "dropout_rate": 0.25,
        "false_prior_rate": 0.18,
        "cohort_noise_sd": 0.70
    }
}

METHODS = [
    "random",
    "mean_signature",
    "simple_embedding_nn",
    "simple_ot",
    "graph_ot",
    "neurotrace_lite"
]


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def zscore(x, eps=1e-8):
    x = np.asarray(x, dtype=float)
    return (x - np.nanmean(x)) / (np.nanstd(x) + eps)


def softmax(x, temp=1.0):
    x = np.asarray(x, dtype=float) / max(temp, 1e-8)
    x = x - np.nanmax(x)
    ex = np.exp(x)
    return ex / (np.sum(ex) + 1e-12)


def corr(a, b):
    a = zscore(a)
    b = zscore(b)
    return float(np.nanmean(a * b))


def cosine(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.dot(a, b) / ((np.linalg.norm(a) + 1e-12) * (np.linalg.norm(b) + 1e-12)))


def rankdata_average(x):
    x = np.asarray(x)
    order = np.argsort(x)
    ranks = np.empty(len(x), dtype=float)
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and x[order[j + 1]] == x[order[i]]:
            j += 1
        avg = (i + j + 2) / 2.0
        ranks[order[i:j + 1]] = avg
        i = j + 1
    return ranks


def auroc_score(y_true, y_score):
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    pos = y_true == 1
    neg = y_true == 0
    n_pos = int(pos.sum())
    n_neg = int(neg.sum())
    if n_pos == 0 or n_neg == 0:
        return np.nan
    ranks = rankdata_average(y_score)
    rank_sum_pos = ranks[pos].sum()
    auc = (rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def average_precision_score_manual(y_true, y_score):
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    n_pos = int(y_true.sum())
    if n_pos == 0:
        return np.nan
    order = np.argsort(-y_score)
    y = y_true[order]
    tp = np.cumsum(y)
    fp = np.cumsum(1 - y)
    precision = tp / np.maximum(tp + fp, 1)
    recall = tp / n_pos
    recall_prev = np.concatenate([[0.0], recall[:-1]])
    ap = np.sum((recall - recall_prev) * precision)
    return float(ap)


def assign_communities(rng, n_genes, n_communities):
    return rng.integers(0, n_communities, size=n_genes)


def community_smooth(x, communities, alpha=0.35):
    x = np.asarray(x, dtype=float)
    out = x.copy()
    for c in np.unique(communities):
        idx = np.where(communities == c)[0]
        if len(idx) > 0:
            out[idx] = (1 - alpha) * x[idx] + alpha * np.mean(x[idx])
    return out


def simulate_reference(rng):
    genes = np.array([f"G{i:05d}" for i in range(N_GENES)])
    states = np.array([f"DevState_{i + 1}" for i in range(N_STATES)])
    communities = assign_communities(rng, N_GENES, N_COMMUNITIES)

    # community centers for embedding
    comm_centers = rng.normal(0, 1, size=(N_COMMUNITIES, EMBED_DIM))
    gene_embed = comm_centers[communities] + rng.normal(0, 0.45, size=(N_GENES, EMBED_DIM))

    dev = rng.normal(0, 0.12, size=(N_STATES, N_GENES))

    state_marker_sets = []
    for k in range(N_STATES):
        marker_communities = rng.choice(N_COMMUNITIES, size=4, replace=False)
        idx = np.where(np.isin(communities, marker_communities))[0]
        if len(idx) > 120:
            idx = rng.choice(idx, size=120, replace=False)
        sign = rng.choice([-1, 1], size=len(idx), p=[0.25, 0.75])
        dev[k, idx] += sign * rng.normal(1.8, 0.35, size=len(idx))
        state_marker_sets.append(set(idx.tolist()))

    dev = np.apply_along_axis(zscore, 1, dev)

    state_embed = []
    for k in range(N_STATES):
        w = np.abs(dev[k])
        top = np.argsort(-w)[:150]
        weights = w[top] / (w[top].sum() + 1e-12)
        state_embed.append((gene_embed[top] * weights[:, None]).sum(axis=0))
    state_embed = np.vstack(state_embed)

    return {
        "genes": genes,
        "states": states,
        "communities": communities,
        "gene_embed": gene_embed,
        "state_embed": state_embed,
        "dev": dev,
        "state_marker_sets": state_marker_sets
    }


def simulate_one_replicate(rng, condition_name, condition_params):
    ref = simulate_reference(rng)

    true_primary = int(rng.integers(0, N_STATES))
    remaining = [x for x in range(N_STATES) if x != true_primary]
    true_secondary = int(rng.choice(remaining))
    true_states = [true_primary, true_secondary]

    true_mix = 0.72 * ref["dev"][true_primary] + 0.38 * ref["dev"][true_secondary]
    true_mix = zscore(true_mix)

    true_gene_set = set()
    for s in true_states:
        top_idx = np.argsort(-np.abs(ref["dev"][s]))[:140]
        true_gene_set.update(top_idx.tolist())

    y_true_gene = np.zeros(N_GENES, dtype=int)
    y_true_gene[list(true_gene_set)] = 1

    # risk prior: enriched among true genes but with false positives
    risk_prior = np.zeros(N_GENES)
    true_idx = np.array(list(true_gene_set))
    non_true_idx = np.where(y_true_gene == 0)[0]
    n_true_prior = max(1, int(0.34 * len(true_idx)))
    n_false_prior = max(1, int(condition_params["false_prior_rate"] * len(non_true_idx)))

    risk_true = rng.choice(true_idx, size=min(n_true_prior, len(true_idx)), replace=False)
    risk_false = rng.choice(non_true_idx, size=min(n_false_prior, len(non_true_idx)), replace=False)
    risk_prior[risk_true] = 1.0
    risk_prior[risk_false] = 1.0

    # composition confound is partly related to support/neural background but not true disease
    comp_state = int(rng.choice([x for x in range(N_STATES) if x not in true_states]))
    comp_vector = zscore(ref["dev"][comp_state] + rng.normal(0, 0.3, size=N_GENES))

    cohort_effects = []
    cohort_vectors = []
    for c in range(N_COHORTS):
        batch_noise = rng.normal(0, condition_params["batch_sd"], size=N_GENES)
        cohort_noise = rng.normal(0, condition_params["cohort_noise_sd"], size=N_GENES)
        comp_weight = rng.normal(0.0, condition_params["composition_sd"])

        observed = true_mix + comp_weight * comp_vector + batch_noise + cohort_noise

        # dropout masks gene-level disease signal
        dropout_mask = rng.random(N_GENES) < condition_params["dropout_rate"]
        observed[dropout_mask] = observed[dropout_mask] * rng.normal(0.05, 0.03, size=dropout_mask.sum())

        observed = zscore(observed)
        cohort_vectors.append(observed)
        cohort_effects.append({
            "cohort": f"cohort_{c + 1}",
            "comp_weight": comp_weight
        })

    cohort_matrix = np.vstack(cohort_vectors)
    mean_disease = zscore(cohort_matrix.mean(axis=0))
    sd_disease = cohort_matrix.std(axis=0)
    invariance = 1.0 / (1.0 + sd_disease)
    invariance = zscore(invariance)

    return ref, {
        "condition": condition_name,
        "true_primary": true_primary,
        "true_secondary": true_secondary,
        "true_states": true_states,
        "true_mix": true_mix,
        "y_true_gene": y_true_gene,
        "risk_prior": risk_prior,
        "cohort_matrix": cohort_matrix,
        "mean_disease": mean_disease,
        "invariance": invariance,
        "composition_state": comp_state
    }


def disease_embedding_from_scores(scores, gene_embed, top_n=180):
    scores = np.asarray(scores)
    idx = np.argsort(-np.abs(scores))[:top_n]
    weights = np.abs(scores[idx])
    weights = weights / (weights.sum() + 1e-12)
    return (gene_embed[idx] * weights[:, None]).sum(axis=0)


def state_scores_mean_signature(disease, dev):
    return np.array([corr(disease, dev[k]) for k in range(dev.shape[0])])


def state_scores_embedding_nn(disease, gene_embed, state_embed):
    emb = disease_embedding_from_scores(disease, gene_embed)
    return np.array([cosine(emb, state_embed[k]) for k in range(state_embed.shape[0])])


def state_scores_simple_ot(disease, dev, gene_embed, state_embed):
    sig = state_scores_mean_signature(disease, dev)
    emb = state_scores_embedding_nn(disease, gene_embed, state_embed)
    combined = 0.60 * zscore(sig) + 0.40 * zscore(emb)
    return combined


def state_scores_graph_ot(disease, dev, gene_embed, state_embed, communities):
    sm = community_smooth(disease, communities, alpha=0.45)
    sig = state_scores_mean_signature(sm, dev)
    emb = state_scores_embedding_nn(sm, gene_embed, state_embed)
    return 0.60 * zscore(sig) + 0.40 * zscore(emb)


def state_scores_neurotrace_lite(disease, dev, gene_embed, state_embed, communities, risk_prior, invariance):
    # disease vector is stabilized by cross-cohort invariance and graph smoothing
    stable_disease = zscore(disease) * (1.0 + 0.25 * zscore(invariance))
    stable_disease = zscore(stable_disease)
    sm = community_smooth(stable_disease, communities, alpha=0.42)

    sig = state_scores_mean_signature(sm, dev)
    emb = state_scores_embedding_nn(sm, gene_embed, state_embed)

    risk_scores = []
    for k in range(dev.shape[0]):
        top = np.argsort(-np.abs(dev[k]))[:180]
        risk_scores.append(np.mean(risk_prior[top]))
    risk_scores = np.array(risk_scores)

    combined = (
        0.52 * zscore(sig) +
        0.28 * zscore(emb) +
        0.20 * zscore(risk_scores)
    )
    return combined


def gene_scores_by_method(method, sim, ref):
    disease = sim["mean_disease"]
    risk = sim["risk_prior"]
    inv = sim["invariance"]
    communities = ref["communities"]

    if method == "random":
        rng = np.random.default_rng(SEED + 999)
        return rng.random(len(disease))

    if method == "mean_signature":
        return np.abs(disease)

    if method == "simple_embedding_nn":
        return np.abs(disease)

    if method == "simple_ot":
        return np.abs(disease)

    if method == "graph_ot":
        return np.abs(community_smooth(disease, communities, alpha=0.45))

    if method == "neurotrace_lite":
        base = np.abs(community_smooth(disease, communities, alpha=0.42))
        score = 0.70 * zscore(base) + 0.18 * zscore(risk) + 0.12 * zscore(inv)
        return zscore(score)

    raise ValueError(method)


def state_scores_by_method(method, sim, ref):
    disease = sim["mean_disease"]

    if method == "random":
        return np.zeros(N_STATES)

    if method == "mean_signature":
        return state_scores_mean_signature(disease, ref["dev"])

    if method == "simple_embedding_nn":
        return state_scores_embedding_nn(disease, ref["gene_embed"], ref["state_embed"])

    if method == "simple_ot":
        return state_scores_simple_ot(disease, ref["dev"], ref["gene_embed"], ref["state_embed"])

    if method == "graph_ot":
        return state_scores_graph_ot(
            disease,
            ref["dev"],
            ref["gene_embed"],
            ref["state_embed"],
            ref["communities"]
        )

    if method == "neurotrace_lite":
        return state_scores_neurotrace_lite(
            disease,
            ref["dev"],
            ref["gene_embed"],
            ref["state_embed"],
            ref["communities"],
            sim["risk_prior"],
            sim["invariance"]
        )

    raise ValueError(method)


def evaluate_method(method, rep_id, condition_name, sim, ref):
    scores = state_scores_by_method(method, sim, ref)

    if method == "random":
        pred_state = int((rep_id + 3) % N_STATES)
    else:
        pred_state = int(np.argmax(scores))

    true_states = sim["true_states"]
    alignment_hit_top1 = int(pred_state in true_states)

    state_order = np.argsort(-scores)
    rank_primary = int(np.where(state_order == sim["true_primary"])[0][0] + 1)
    rank_secondary = int(np.where(state_order == sim["true_secondary"])[0][0] + 1)
    best_true_rank = min(rank_primary, rank_secondary)

    gene_score = gene_scores_by_method(method, sim, ref)
    y = sim["y_true_gene"]
    gene_auroc = auroc_score(y, gene_score)
    gene_auprc = average_precision_score_manual(y, gene_score)

    # transport mass on true states
    transport = softmax(scores, temp=0.75)
    true_transport_mass = float(transport[true_states].sum())

    return {
        "replicate": rep_id,
        "condition": condition_name,
        "method": method,
        "true_primary_state": ref["states"][sim["true_primary"]],
        "true_secondary_state": ref["states"][sim["true_secondary"]],
        "predicted_state": ref["states"][pred_state],
        "alignment_hit_top1": alignment_hit_top1,
        "best_true_state_rank": best_true_rank,
        "true_transport_mass": true_transport_mass,
        "gene_auroc": gene_auroc,
        "gene_auprc": gene_auprc
    }


def summarize(df):
    group_cols = ["condition", "method"]
    rows = []

    for (cond, method), sub in df.groupby(group_cols):
        row = {
            "condition": cond,
            "method": method,
            "n_replicates": len(sub),
            "alignment_accuracy_mean": sub["alignment_hit_top1"].mean(),
            "alignment_accuracy_sd": sub["alignment_hit_top1"].std(ddof=1),
            "best_true_rank_mean": sub["best_true_state_rank"].mean(),
            "true_transport_mass_mean": sub["true_transport_mass"].mean(),
            "true_transport_mass_sd": sub["true_transport_mass"].std(ddof=1),
            "gene_auroc_mean": sub["gene_auroc"].mean(),
            "gene_auroc_sd": sub["gene_auroc"].std(ddof=1),
            "gene_auprc_mean": sub["gene_auprc"].mean(),
            "gene_auprc_sd": sub["gene_auprc"].std(ddof=1)
        }
        rows.append(row)

    out = pd.DataFrame(rows)
    out = out.sort_values(["condition", "alignment_accuracy_mean", "gene_auprc_mean"], ascending=[True, False, False])
    return out


def make_ranking(summary_df):
    rows = []
    for cond, sub in summary_df.groupby("condition"):
        sub = sub.copy()
        sub["rank_alignment"] = sub["alignment_accuracy_mean"].rank(ascending=False, method="min")
        sub["rank_gene_auprc"] = sub["gene_auprc_mean"].rank(ascending=False, method="min")
        sub["rank_transport"] = sub["true_transport_mass_mean"].rank(ascending=False, method="min")
        sub["mean_rank"] = sub[["rank_alignment", "rank_gene_auprc", "rank_transport"]].mean(axis=1)
        rows.append(sub)
    rank = pd.concat(rows, ignore_index=True)
    rank = rank.sort_values(["condition", "mean_rank", "method"])
    return rank


def try_plot(summary_df):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        for metric, fname, ylabel in [
            ("alignment_accuracy_mean", "05_step02A_alignment_accuracy_by_condition.pdf", "Top-1 developmental alignment accuracy"),
            ("gene_auprc_mean", "06_step02A_gene_AUPRC_by_condition.pdf", "Gene prioritization AUPRC"),
            ("true_transport_mass_mean", "07_step02A_true_transport_mass_by_condition.pdf", "Transport mass assigned to true states")
        ]:
            pivot = summary_df.pivot(index="condition", columns="method", values=metric)
            pivot = pivot.loc[list(CONDITIONS.keys()), METHODS]

            ax = pivot.plot(kind="bar", figsize=(12, 5), width=0.82)
            ax.set_ylabel(ylabel)
            ax.set_xlabel("Simulation condition")
            ax.set_title(f"NeuroTRACE Step02A simulation benchmark: {ylabel}")
            ax.legend(title="Method", bbox_to_anchor=(1.02, 1.0), loc="upper left")
            plt.tight_layout()
            plt.savefig(FIG / fname)
            plt.close()

    except Exception as e:
        with open(OUT / "plotting_failed.txt", "w") as f:
            f.write(str(e) + "\n")


def main():
    print(f"[{now()}] NeuroTRACE Step02A simulation benchmark started")
    print(f"[{now()}] N_REP per condition: {N_REP}")
    print(f"[{now()}] N_GENES: {N_GENES}; N_STATES: {N_STATES}; N_COHORTS: {N_COHORTS}")

    rng = np.random.default_rng(SEED)

    config_rows = []
    for cond, pars in CONDITIONS.items():
        row = {"condition": cond}
        row.update(pars)
        config_rows.append(row)
    pd.DataFrame(config_rows).to_csv(OUT / "01_step02A_simulation_config.tsv", sep="\t", index=False)

    all_metrics = []

    rep_global = 0
    for cond_name, cond_params in CONDITIONS.items():
        print(f"[{now()}] Running condition: {cond_name}")
        for rep in range(1, N_REP + 1):
            rep_global += 1
            ref, sim = simulate_one_replicate(rng, cond_name, cond_params)

            for method in METHODS:
                all_metrics.append(
                    evaluate_method(
                        method=method,
                        rep_id=rep_global,
                        condition_name=cond_name,
                        sim=sim,
                        ref=ref
                    )
                )

    metrics = pd.DataFrame(all_metrics)
    metrics.to_csv(OUT / "02_step02A_replicate_metrics.tsv", sep="\t", index=False)

    summary = summarize(metrics)
    summary.to_csv(OUT / "03_step02A_summary_by_condition_method.tsv", sep="\t", index=False)

    ranking = make_ranking(summary)
    ranking.to_csv(OUT / "04_step02A_method_ranking.tsv", sep="\t", index=False)

    try_plot(summary)

    # compact global summary
    global_summary = (
        metrics
        .groupby("method")
        .agg(
            n=("replicate", "count"),
            alignment_accuracy_mean=("alignment_hit_top1", "mean"),
            best_true_rank_mean=("best_true_state_rank", "mean"),
            true_transport_mass_mean=("true_transport_mass", "mean"),
            gene_auroc_mean=("gene_auroc", "mean"),
            gene_auprc_mean=("gene_auprc", "mean")
        )
        .reset_index()
        .sort_values(["alignment_accuracy_mean", "gene_auprc_mean"], ascending=[False, False])
    )
    global_summary.to_csv(OUT / "05_step02A_global_method_summary.tsv", sep="\t", index=False)

    best_method = global_summary.iloc[0]["method"]

    with open(OUT / "06_step02A_simulation_benchmark_summary.md", "w") as f:
        f.write("# NeuroTRACE Step02A simulation benchmark summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Simulation design\n")
        f.write(f"- Replicates per condition: {N_REP}\n")
        f.write(f"- Genes: {N_GENES}\n")
        f.write(f"- Developmental states: {N_STATES}\n")
        f.write(f"- Adult cohorts per replicate: {N_COHORTS}\n")
        f.write(f"- Conditions: {', '.join(CONDITIONS.keys())}\n")
        f.write(f"- Methods: {', '.join(METHODS)}\n\n")
        f.write("## Main endpoints\n")
        f.write("- Top-1 disease-to-development alignment accuracy\n")
        f.write("- Rank of true developmental state\n")
        f.write("- Transport mass assigned to true developmental states\n")
        f.write("- Gene prioritization AUROC/AUPRC\n\n")
        f.write("## Global method summary\n\n")
        f.write(global_summary.to_markdown(index=False))
        f.write("\n\n")
        f.write(f"Best global method by alignment accuracy and AUPRC: **{best_method}**\n\n")
        f.write("## Interpretation template\n")
        f.write("If NeuroTRACE-lite ranks first across stress conditions, this supports the algorithmic claim that cross-cohort invariance, graph smoothing, embedding similarity, and risk-prior integration improve disease-to-development alignment beyond simple scoring or simple OT baselines.\n")

    print(f"[{now()}] Step02A done")
    print(f"[{now()}] Results: {OUT}")
    print(f"[{now()}] Figures: {FIG}")


if __name__ == "__main__":
    main()
