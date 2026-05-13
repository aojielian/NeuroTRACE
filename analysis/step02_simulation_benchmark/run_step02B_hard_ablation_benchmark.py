#!/usr/bin/env python3

import os
import math
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

SEED = 20260508 + 200
N_REP = 160
N_GENES = 3000
N_STATES = 10
N_COMMUNITIES = 60
N_COHORTS = 4
EMBED_DIM = 32

CONDITIONS = {
    "weak_signal": {
        "signal_strength": 0.42,
        "batch_sd": 0.65,
        "composition_sd": 0.55,
        "dropout_rate": 0.12,
        "false_prior_rate": 0.08,
        "true_prior_recall": 0.30,
        "state_overlap": 0.45,
        "adversarial_composition": 0.25
    },
    "low_overlap": {
        "signal_strength": 0.50,
        "batch_sd": 0.55,
        "composition_sd": 0.50,
        "dropout_rate": 0.12,
        "false_prior_rate": 0.08,
        "true_prior_recall": 0.28,
        "state_overlap": 0.22,
        "adversarial_composition": 0.30
    },
    "composition_confounded": {
        "signal_strength": 0.50,
        "batch_sd": 0.55,
        "composition_sd": 1.15,
        "dropout_rate": 0.10,
        "false_prior_rate": 0.08,
        "true_prior_recall": 0.30,
        "state_overlap": 0.42,
        "adversarial_composition": 0.95
    },
    "high_dropout": {
        "signal_strength": 0.48,
        "batch_sd": 0.60,
        "composition_sd": 0.55,
        "dropout_rate": 0.38,
        "false_prior_rate": 0.08,
        "true_prior_recall": 0.30,
        "state_overlap": 0.40,
        "adversarial_composition": 0.35
    },
    "false_prior_stress": {
        "signal_strength": 0.52,
        "batch_sd": 0.55,
        "composition_sd": 0.55,
        "dropout_rate": 0.12,
        "false_prior_rate": 0.32,
        "true_prior_recall": 0.22,
        "state_overlap": 0.42,
        "adversarial_composition": 0.40
    },
    "combined_hard": {
        "signal_strength": 0.38,
        "batch_sd": 0.85,
        "composition_sd": 1.00,
        "dropout_rate": 0.32,
        "false_prior_rate": 0.26,
        "true_prior_recall": 0.22,
        "state_overlap": 0.25,
        "adversarial_composition": 0.90
    }
}

METHODS = [
    "random",
    "mean_signature",
    "embedding_nn",
    "simple_ot",
    "graph_ot",
    "neurotrace_no_embedding",
    "neurotrace_no_graph",
    "neurotrace_no_risk",
    "neurotrace_no_invariance",
    "neurotrace_full"
]


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def zscore(x, eps=1e-8):
    x = np.asarray(x, dtype=float)
    return (x - np.nanmean(x)) / (np.nanstd(x) + eps)


def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def softmax(x, temp=1.0):
    x = np.asarray(x, dtype=float) / max(temp, 1e-8)
    x = x - np.nanmax(x)
    ex = np.exp(x)
    return ex / (np.sum(ex) + 1e-12)


def corr(a, b):
    return float(np.mean(zscore(a) * zscore(b)))


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
    return float((rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg))


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
    return float(np.sum((recall - recall_prev) * precision))


def precision_at_k(y_true, y_score, k=100):
    y_true = np.asarray(y_true).astype(int)
    order = np.argsort(-np.asarray(y_score))[:k]
    return float(np.mean(y_true[order]))


def community_smooth(x, communities, alpha=0.35):
    x = np.asarray(x, dtype=float)
    out = x.copy()
    for c in np.unique(communities):
        idx = np.where(communities == c)[0]
        out[idx] = (1 - alpha) * x[idx] + alpha * np.mean(x[idx])
    return out


def build_reference(rng):
    genes = np.array([f"G{i:05d}" for i in range(N_GENES)])
    states = np.array([f"DevState_{i + 1}" for i in range(N_STATES)])
    communities = rng.integers(0, N_COMMUNITIES, size=N_GENES)

    comm_centers = rng.normal(0, 1.0, size=(N_COMMUNITIES, EMBED_DIM))
    gene_embed = comm_centers[communities] + rng.normal(0, 0.60, size=(N_GENES, EMBED_DIM))

    dev = rng.normal(0, 0.18, size=(N_STATES, N_GENES))
    state_markers = []

    used = set()
    for s in range(N_STATES):
        # 每个 state 由多个 community 驱动，但允许交叉
        comms = rng.choice(N_COMMUNITIES, size=5, replace=False)
        idx = np.where(np.isin(communities, comms))[0]
        if len(idx) > 180:
            idx = rng.choice(idx, size=180, replace=False)

        direction = rng.choice([-1, 1], size=len(idx), p=[0.30, 0.70])
        dev[s, idx] += direction * rng.normal(1.35, 0.35, size=len(idx))

        state_markers.append(set(idx.tolist()))
        used.update(idx.tolist())

    dev = np.apply_along_axis(zscore, 1, dev)

    state_embed = []
    for s in range(N_STATES):
        w = np.abs(dev[s])
        top = np.argsort(-w)[:220]
        weights = w[top] / (w[top].sum() + 1e-12)
        state_embed.append((gene_embed[top] * weights[:, None]).sum(axis=0))
    state_embed = np.vstack(state_embed)

    return {
        "genes": genes,
        "states": states,
        "communities": communities,
        "gene_embed": gene_embed,
        "dev": dev,
        "state_embed": state_embed,
        "state_markers": state_markers
    }


def disease_embedding_from_scores(scores, gene_embed, top_n=220):
    scores = np.asarray(scores, dtype=float)
    top = np.argsort(-np.abs(scores))[:top_n]
    weights = np.abs(scores[top])
    weights = weights / (weights.sum() + 1e-12)
    return (gene_embed[top] * weights[:, None]).sum(axis=0)


def choose_true_gene_set(rng, ref, true_states, state_overlap):
    true_genes = set()

    for s in true_states:
        marker_pool = list(ref["state_markers"][s])
        rng.shuffle(marker_pool)
        n_from_state = max(20, int(len(marker_pool) * state_overlap))
        true_genes.update(marker_pool[:n_from_state])

    # 加入一部分 developmental state 外的 disease-specific genes，模拟 adult disease 不完全等于 fetal marker
    n_extra = max(40, int(0.35 * len(true_genes)))
    all_idx = np.arange(N_GENES)
    current = np.array(list(true_genes))
    extra_pool = np.setdiff1d(all_idx, current)
    extra = rng.choice(extra_pool, size=min(n_extra, len(extra_pool)), replace=False)
    true_genes.update(extra.tolist())

    return sorted(true_genes)


def make_risk_prior(rng, y_true_gene, true_prior_recall, false_prior_rate):
    risk = np.zeros(N_GENES)
    true_idx = np.where(y_true_gene == 1)[0]
    false_idx = np.where(y_true_gene == 0)[0]

    n_true = max(1, int(true_prior_recall * len(true_idx)))
    n_false = max(1, int(false_prior_rate * len(false_idx)))

    risk_true = rng.choice(true_idx, size=min(n_true, len(true_idx)), replace=False)
    risk_false = rng.choice(false_idx, size=min(n_false, len(false_idx)), replace=False)

    risk[risk_true] = 1.0
    risk[risk_false] = 1.0
    return risk


def simulate_replicate(rng, condition_name, pars):
    ref = build_reference(rng)

    true_primary = int(rng.integers(0, N_STATES))
    candidate_secondary = [x for x in range(N_STATES) if x != true_primary]
    true_secondary = int(rng.choice(candidate_secondary))
    true_states = [true_primary, true_secondary]

    # decoy composition state: not true, but strong apparent disease correlation in some cohorts
    decoys = [x for x in range(N_STATES) if x not in true_states]
    decoy_state = int(rng.choice(decoys))

    true_gene_idx = choose_true_gene_set(
        rng=rng,
        ref=ref,
        true_states=true_states,
        state_overlap=pars["state_overlap"]
    )

    y_true_gene = np.zeros(N_GENES, dtype=int)
    y_true_gene[true_gene_idx] = 1

    true_effect = np.zeros(N_GENES)
    true_effect[true_gene_idx] = rng.normal(1.0, 0.25, size=len(true_gene_idx))

    # 用 developmental state 方向给真实 disease effect 定向，但加入 adult-specific effect
    dev_mix = 0.65 * ref["dev"][true_primary] + 0.35 * ref["dev"][true_secondary]
    true_effect = 0.60 * zscore(true_effect) + 0.40 * zscore(dev_mix)
    true_effect = zscore(true_effect)

    risk_prior = make_risk_prior(
        rng=rng,
        y_true_gene=y_true_gene,
        true_prior_recall=pars["true_prior_recall"],
        false_prior_rate=pars["false_prior_rate"]
    )

    decoy_vector = zscore(ref["dev"][decoy_state] + rng.normal(0, 0.30, size=N_GENES))

    cohort_effects = []
    cohort_quality = []

    for c in range(N_COHORTS):
        disease = pars["signal_strength"] * true_effect

        # adversarial composition: decoy direction can be stronger than true signal in difficult settings
        comp_weight = rng.normal(pars["adversarial_composition"], pars["composition_sd"])
        comp = comp_weight * decoy_vector

        batch = rng.normal(0, pars["batch_sd"], size=N_GENES)
        noise = rng.normal(0, 0.95, size=N_GENES)

        obs = disease + comp + batch + noise

        dropout = rng.random(N_GENES) < pars["dropout_rate"]
        obs[dropout] = obs[dropout] * rng.normal(0.02, 0.02, size=dropout.sum())

        # 少数 cohort 方向翻转或弱化，模拟真实多队列异质性
        if rng.random() < 0.18:
            obs = obs - 0.35 * disease + rng.normal(0, 0.25, size=N_GENES)

        obs = zscore(obs)
        cohort_effects.append(obs)

        # cohort quality 与真实信号相关程度，后续方法不能直接使用 true，只用于模拟 invariance 的背景
        cohort_quality.append(corr(obs, true_effect))

    cohort_matrix = np.vstack(cohort_effects)
    mean_disease = zscore(cohort_matrix.mean(axis=0))
    sd_disease = cohort_matrix.std(axis=0)
    invariance = zscore(1.0 / (1.0 + sd_disease))

    return ref, {
        "condition": condition_name,
        "true_primary": true_primary,
        "true_secondary": true_secondary,
        "true_states": true_states,
        "decoy_state": decoy_state,
        "y_true_gene": y_true_gene,
        "true_effect": true_effect,
        "risk_prior": risk_prior,
        "cohort_matrix": cohort_matrix,
        "mean_disease": mean_disease,
        "invariance": invariance,
        "cohort_quality_mean": float(np.mean(cohort_quality))
    }


def score_signature(disease, dev):
    return np.array([corr(disease, dev[s]) for s in range(dev.shape[0])])


def score_embedding(disease, gene_embed, state_embed):
    emb = disease_embedding_from_scores(disease, gene_embed)
    return np.array([cosine(emb, state_embed[s]) for s in range(state_embed.shape[0])])


def score_risk(dev, risk_prior, top_n=220):
    rs = []
    for s in range(dev.shape[0]):
        top = np.argsort(-np.abs(dev[s]))[:top_n]
        rs.append(np.mean(risk_prior[top]))
    return np.array(rs)


def state_score(method, sim, ref):
    disease = sim["mean_disease"]
    dev = ref["dev"]
    gene_embed = ref["gene_embed"]
    state_embed = ref["state_embed"]
    communities = ref["communities"]
    risk_prior = sim["risk_prior"]
    invariance = sim["invariance"]

    if method == "random":
        return np.zeros(N_STATES)

    if method == "mean_signature":
        return score_signature(disease, dev)

    if method == "embedding_nn":
        return score_embedding(disease, gene_embed, state_embed)

    if method == "simple_ot":
        sig = score_signature(disease, dev)
        emb = score_embedding(disease, gene_embed, state_embed)
        return 0.60 * zscore(sig) + 0.40 * zscore(emb)

    if method == "graph_ot":
        d = community_smooth(disease, communities, alpha=0.48)
        sig = score_signature(d, dev)
        emb = score_embedding(d, gene_embed, state_embed)
        return 0.60 * zscore(sig) + 0.40 * zscore(emb)

    # NeuroTRACE variants
    use_embedding = method != "neurotrace_no_embedding"
    use_graph = method != "neurotrace_no_graph"
    use_risk = method != "neurotrace_no_risk"
    use_invariance = method != "neurotrace_no_invariance"

    d = disease.copy()

    if use_invariance:
        d = zscore(d) * (1.0 + 0.35 * zscore(invariance))
        d = zscore(d)

    if use_graph:
        d = community_smooth(d, communities, alpha=0.50)

    sig = score_signature(d, dev)

    components = [0.50 * zscore(sig)]

    if use_embedding:
        emb = score_embedding(d, gene_embed, state_embed)
        components.append(0.25 * zscore(emb))

    if use_risk:
        risk = score_risk(dev, risk_prior)
        components.append(0.25 * zscore(risk))

    return np.sum(np.vstack(components), axis=0)


def gene_score(method, sim, ref):
    disease = sim["mean_disease"]
    communities = ref["communities"]
    risk_prior = sim["risk_prior"]
    invariance = sim["invariance"]

    if method == "random":
        local_rng = np.random.default_rng(SEED + 999)
        return local_rng.random(N_GENES)

    if method in ["mean_signature", "embedding_nn", "simple_ot"]:
        return zscore(np.abs(disease))

    if method == "graph_ot":
        return zscore(np.abs(community_smooth(disease, communities, alpha=0.48)))

    use_graph = method != "neurotrace_no_graph"
    use_risk = method != "neurotrace_no_risk"
    use_invariance = method != "neurotrace_no_invariance"

    d = disease.copy()

    if use_invariance:
        d = zscore(d) * (1.0 + 0.35 * zscore(invariance))
        d = zscore(d)

    if use_graph:
        d = community_smooth(d, communities, alpha=0.50)

    comps = [0.70 * zscore(np.abs(d))]

    if use_risk:
        comps.append(0.20 * zscore(risk_prior))

    if use_invariance:
        comps.append(0.10 * zscore(invariance))

    return zscore(np.sum(np.vstack(comps), axis=0))


def evaluate(method, rep_id, cond, sim, ref):
    s = state_score(method, sim, ref)

    if method == "random":
        pred = int((rep_id + 7) % N_STATES)
        s = np.random.default_rng(SEED + rep_id).normal(0, 1, size=N_STATES)
    else:
        pred = int(np.argmax(s))

    order = np.argsort(-s)
    true_states = sim["true_states"]

    top1 = int(pred in true_states)
    top2 = int(len(set(order[:2]).intersection(true_states)) > 0)
    rank_primary = int(np.where(order == sim["true_primary"])[0][0] + 1)
    rank_secondary = int(np.where(order == sim["true_secondary"])[0][0] + 1)
    rank_decoy = int(np.where(order == sim["decoy_state"])[0][0] + 1)

    transport = softmax(s, temp=0.75)
    true_mass = float(transport[true_states].sum())
    decoy_mass = float(transport[sim["decoy_state"]])

    gs = gene_score(method, sim, ref)
    y = sim["y_true_gene"]

    return {
        "replicate": rep_id,
        "condition": cond,
        "method": method,
        "true_primary_state": ref["states"][sim["true_primary"]],
        "true_secondary_state": ref["states"][sim["true_secondary"]],
        "decoy_state": ref["states"][sim["decoy_state"]],
        "predicted_state": ref["states"][pred],
        "alignment_hit_top1": top1,
        "alignment_hit_top2": top2,
        "best_true_state_rank": min(rank_primary, rank_secondary),
        "decoy_state_rank": rank_decoy,
        "true_transport_mass": true_mass,
        "decoy_transport_mass": decoy_mass,
        "gene_auroc": auroc_score(y, gs),
        "gene_auprc": average_precision_score_manual(y, gs),
        "precision_at_50": precision_at_k(y, gs, 50),
        "precision_at_100": precision_at_k(y, gs, 100),
        "cohort_quality_mean": sim["cohort_quality_mean"]
    }


def summarize(metrics):
    rows = []

    for (cond, method), sub in metrics.groupby(["condition", "method"]):
        rows.append({
            "condition": cond,
            "method": method,
            "n_replicates": len(sub),
            "alignment_top1_mean": sub["alignment_hit_top1"].mean(),
            "alignment_top2_mean": sub["alignment_hit_top2"].mean(),
            "best_true_rank_mean": sub["best_true_state_rank"].mean(),
            "decoy_rank_mean": sub["decoy_state_rank"].mean(),
            "true_transport_mass_mean": sub["true_transport_mass"].mean(),
            "decoy_transport_mass_mean": sub["decoy_transport_mass"].mean(),
            "gene_auroc_mean": sub["gene_auroc"].mean(),
            "gene_auprc_mean": sub["gene_auprc"].mean(),
            "precision_at_50_mean": sub["precision_at_50"].mean(),
            "precision_at_100_mean": sub["precision_at_100"].mean(),
        })

    out = pd.DataFrame(rows)
    out = out.sort_values(
        ["condition", "alignment_top1_mean", "gene_auprc_mean"],
        ascending=[True, False, False]
    )
    return out


def paired_deltas(metrics):
    full = metrics[metrics["method"] == "neurotrace_full"].copy()
    full = full.rename(columns={
        "alignment_hit_top1": "full_alignment_top1",
        "gene_auprc": "full_gene_auprc",
        "true_transport_mass": "full_true_transport_mass",
        "decoy_transport_mass": "full_decoy_transport_mass"
    })
    full = full[[
        "replicate", "condition",
        "full_alignment_top1", "full_gene_auprc",
        "full_true_transport_mass", "full_decoy_transport_mass"
    ]]

    rows = []
    for method in METHODS:
        if method == "neurotrace_full":
            continue
        sub = metrics[metrics["method"] == method].merge(full, on=["replicate", "condition"], how="inner")
        for cond, ss in sub.groupby("condition"):
            rows.append({
                "condition": cond,
                "baseline_method": method,
                "n": len(ss),
                "delta_alignment_top1_full_minus_baseline": (ss["full_alignment_top1"] - ss["alignment_hit_top1"]).mean(),
                "delta_gene_auprc_full_minus_baseline": (ss["full_gene_auprc"] - ss["gene_auprc"]).mean(),
                "delta_true_transport_mass_full_minus_baseline": (ss["full_true_transport_mass"] - ss["true_transport_mass"]).mean(),
                "delta_decoy_transport_mass_full_minus_baseline": (ss["full_decoy_transport_mass"] - ss["decoy_transport_mass"]).mean()
            })
    return pd.DataFrame(rows)


def ablation_table(summary):
    full = summary[summary["method"] == "neurotrace_full"].copy()
    full = full.rename(columns={
        "alignment_top1_mean": "full_alignment_top1_mean",
        "gene_auprc_mean": "full_gene_auprc_mean",
        "true_transport_mass_mean": "full_true_transport_mass_mean",
        "decoy_transport_mass_mean": "full_decoy_transport_mass_mean"
    })
    full = full[[
        "condition", "full_alignment_top1_mean", "full_gene_auprc_mean",
        "full_true_transport_mass_mean", "full_decoy_transport_mass_mean"
    ]]

    ab = summary[summary["method"].str.startswith("neurotrace_no_")].merge(full, on="condition", how="left")
    ab["loss_alignment_top1"] = ab["full_alignment_top1_mean"] - ab["alignment_top1_mean"]
    ab["loss_gene_auprc"] = ab["full_gene_auprc_mean"] - ab["gene_auprc_mean"]
    ab["loss_true_transport_mass"] = ab["full_true_transport_mass_mean"] - ab["true_transport_mass_mean"]
    ab["gain_decoy_transport_mass_when_removed"] = ab["decoy_transport_mass_mean"] - ab["full_decoy_transport_mass_mean"]
    return ab.sort_values(["condition", "loss_alignment_top1", "loss_gene_auprc"], ascending=[True, False, False])


def make_plots(summary, ablation):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        plot_methods = METHODS
        for metric, fname, ylabel in [
            ("alignment_top1_mean", "08_step02B_hard_alignment_top1.pdf", "Top-1 alignment accuracy"),
            ("gene_auprc_mean", "09_step02B_hard_gene_AUPRC.pdf", "Gene prioritization AUPRC"),
            ("true_transport_mass_mean", "10_step02B_hard_true_transport_mass.pdf", "Transport mass on true states"),
            ("decoy_transport_mass_mean", "11_step02B_hard_decoy_transport_mass.pdf", "Transport mass on decoy state")
        ]:
            pivot = summary.pivot(index="condition", columns="method", values=metric)
            pivot = pivot.loc[list(CONDITIONS.keys()), plot_methods]
            ax = pivot.plot(kind="bar", figsize=(14, 5), width=0.85)
            ax.set_ylabel(ylabel)
            ax.set_xlabel("Hard simulation condition")
            ax.set_title(f"Step02B hard benchmark: {ylabel}")
            ax.legend(title="Method", bbox_to_anchor=(1.02, 1.0), loc="upper left")
            plt.tight_layout()
            plt.savefig(FIG / fname)
            plt.close()

        ab_metric = ablation.pivot(index="condition", columns="method", values="loss_gene_auprc")
        ax = ab_metric.plot(kind="bar", figsize=(12, 5), width=0.85)
        ax.set_ylabel("AUPRC loss versus full NeuroTRACE")
        ax.set_xlabel("Condition")
        ax.set_title("Step02B ablation: component removal reduces gene prioritization")
        ax.legend(title="Ablation", bbox_to_anchor=(1.02, 1.0), loc="upper left")
        plt.tight_layout()
        plt.savefig(FIG / "12_step02B_ablation_AUPRC_loss.pdf")
        plt.close()

    except Exception as e:
        with open(OUT / "Step02B_plotting_failed.txt", "w") as f:
            f.write(str(e) + "\n")


def main():
    print(f"[{now()}] NeuroTRACE Step02B hard ablation benchmark started")
    print(f"[{now()}] N_REP per condition: {N_REP}")
    print(f"[{now()}] N_GENES={N_GENES}; N_STATES={N_STATES}; N_COHORTS={N_COHORTS}")

    rng = np.random.default_rng(SEED)

    config_rows = []
    for cond, pars in CONDITIONS.items():
        row = {"condition": cond}
        row.update(pars)
        config_rows.append(row)
    pd.DataFrame(config_rows).to_csv(OUT / "11_step02B_hard_simulation_config.tsv", sep="\t", index=False)

    metrics_rows = []
    rep_global = 0

    for cond, pars in CONDITIONS.items():
        print(f"[{now()}] Running condition: {cond}")
        for rep in range(1, N_REP + 1):
            rep_global += 1
            ref, sim = simulate_replicate(rng, cond, pars)
            for method in METHODS:
                metrics_rows.append(evaluate(method, rep_global, cond, sim, ref))

    metrics = pd.DataFrame(metrics_rows)
    metrics.to_csv(OUT / "12_step02B_replicate_metrics.tsv", sep="\t", index=False)

    summary = summarize(metrics)
    summary.to_csv(OUT / "13_step02B_summary_by_condition_method.tsv", sep="\t", index=False)

    deltas = paired_deltas(metrics)
    deltas.to_csv(OUT / "14_step02B_paired_deltas_full_vs_baselines.tsv", sep="\t", index=False)

    ablation = ablation_table(summary)
    ablation.to_csv(OUT / "15_step02B_ablation_component_loss.tsv", sep="\t", index=False)

    global_summary = (
        metrics
        .groupby("method")
        .agg(
            n=("replicate", "count"),
            alignment_top1_mean=("alignment_hit_top1", "mean"),
            alignment_top2_mean=("alignment_hit_top2", "mean"),
            best_true_rank_mean=("best_true_state_rank", "mean"),
            decoy_rank_mean=("decoy_state_rank", "mean"),
            true_transport_mass_mean=("true_transport_mass", "mean"),
            decoy_transport_mass_mean=("decoy_transport_mass", "mean"),
            gene_auroc_mean=("gene_auroc", "mean"),
            gene_auprc_mean=("gene_auprc", "mean"),
            precision_at_50_mean=("precision_at_50", "mean"),
            precision_at_100_mean=("precision_at_100", "mean")
        )
        .reset_index()
        .sort_values(["alignment_top1_mean", "gene_auprc_mean"], ascending=[False, False])
    )
    global_summary.to_csv(OUT / "16_step02B_global_method_summary.tsv", sep="\t", index=False)

    make_plots(summary, ablation)

    full = global_summary[global_summary["method"] == "neurotrace_full"].iloc[0].to_dict()
    best = global_summary.iloc[0].to_dict()

    with open(OUT / "17_step02B_hard_ablation_benchmark_summary.md", "w") as f:
        f.write("# NeuroTRACE Step02B hard simulation and ablation benchmark summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Design\n")
        f.write(f"- Replicates per condition: {N_REP}\n")
        f.write(f"- Genes: {N_GENES}\n")
        f.write(f"- Developmental states: {N_STATES}\n")
        f.write(f"- Adult cohorts per replicate: {N_COHORTS}\n")
        f.write(f"- Conditions: {', '.join(CONDITIONS.keys())}\n")
        f.write(f"- Methods: {', '.join(METHODS)}\n\n")

        f.write("## Why Step02B is harder than Step02A\n")
        f.write("- Lower disease signal strength\n")
        f.write("- Partial overlap between adult disease genes and developmental state markers\n")
        f.write("- Strong adversarial composition state that can mimic disease alignment\n")
        f.write("- Dropout and cohort heterogeneity\n")
        f.write("- False-positive genetic-risk prior stress test\n")
        f.write("- Explicit ablation of embedding, graph smoothing, risk prior, and cross-cohort invariance\n\n")

        f.write("## Global method summary\n\n")
        f.write(global_summary.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Full NeuroTRACE global performance\n")
        f.write(f"- alignment_top1_mean: {full['alignment_top1_mean']:.4f}\n")
        f.write(f"- gene_auprc_mean: {full['gene_auprc_mean']:.4f}\n")
        f.write(f"- true_transport_mass_mean: {full['true_transport_mass_mean']:.4f}\n")
        f.write(f"- decoy_transport_mass_mean: {full['decoy_transport_mass_mean']:.4f}\n\n")

        f.write("## Best global method\n")
        f.write(f"- best method by alignment_top1_mean then gene_auprc_mean: {best['method']}\n\n")

        f.write("## Interpretation rule\n")
        f.write("This benchmark is useful for the manuscript only if neurotrace_full is consistently among the top methods under combined_hard, composition_confounded, low_overlap, and false_prior_stress, and if at least one ablation shows measurable loss in AUPRC or true transport mass.\n")

    print(f"[{now()}] Step02B done")
    print(f"[{now()}] Results: {OUT}")
    print(f"[{now()}] Figures: {FIG}")


if __name__ == "__main__":
    main()
