#!/usr/bin/env python3
"""final Part H+I: confirmatory held-out simulation with FRESH seeds for the
single predeclared dual-head method DUAL_HEAD_GENE_FIRST_PPR.

Design (protocol Parts H and I; no tuning):

- Two generative designs reused byte-identically from Stages 3/4 (archived
  step02A / step02B drivers, loaded via filesystem-path-only patch through the
  core runner module s3; baseline PPR baselines loaded from their scripts):
    step02A: 6 conditions x N_REP_CONFIRM = 200 reps, SEED_A fresh;
    step02B: 6 conditions x N_REP_CONFIRM = 200 reps, SEED_B fresh.
  N_REP_CONFIRM = 200 is recorded in the protocol and in table 08 BEFORE any
  result is inspected (protocol I2 minimum 160).
- Fresh seeds not used in Stages 1-4: CONFIRM_SEED_A = 202609041,
  CONFIRM_SEED_B = 202609042. Each design instantiates
  rng = np.random.default_rng(CONFIRM_SEED_x) and replays the archived main
  loop (conditions in dict order, N_REP_CONFIRM consecutive replicates per
  condition, one generator call per replicate:
  step02A simulate_one_replicate(rng, cond, pars) / step02B
  simulate_replicate(rng, cond, pars)), exactly as the archived drivers
  consume the stream. Nothing else draws from the shared rng; archived
  comparator methods that use fixed-seed fixed rngs (chance controls) are
  replayed verbatim with their own fixed seeds.

Method set (predeclared in protocol I3; protocol names are the canonical
method labels in all output tables):
    DUAL_HEAD_GENE_FIRST_PPR        stage head q_magnitude = q_pos + q_neg;
                                    gene head q_signed = q_pos - q_neg,
                                    priority |q_signed|, direction sign
    SIGNED_GENE_FIRST_PPR_SINGLE_HEAD   reference signed stage readout
                                        (projection of q_signed); gene head
                                        identical to the dual head by
                                        construction (same q_signed)
    UNSIGNED_GENE_PPR               q_pos + q_neg for both readouts
    DIRECT_NATIVE_PROJECTION        no diffusion; signed cosine of the
                                    observed mean_disease (reference rule)
    GENE_ONLY_PPR_BASELINE            core core GENE_ONLY_PPR (module node +
                                    module-gene + within-community edges;
                                    stage ranking downstream corr(qg, dev)),
                                    replayed unchanged via s3.evaluate_core
    HETERO_PPR_STANDARD             baseline genuine hetero PPR baseline
    STAGE_CORRELATION               baseline non-graph baseline
    GRAPH_OT / SIMPLE_OT / MEAN_SIGNATURE / RANDOM   archived comparators
    SIMULATION_PROXY     historical-only archived NeuroTRACE
                                    simulation proxy (step02A neurotrace_lite,
                                    step02B neurotrace_full); NOT the final
                                    method (protocol I3)

Stage head (protocol I4): alignment_hit_top1, best_true_state_rank,
true_stage_mass from the magnitude readout; affinity =
softmax(zscore(scores), temp 1.0) (family rule; identical to real-data Part D
and reference); compatible true_transport_mass = softmax(scores, 0.75) on true
states (archived convention).
Gene head (protocol I4): gene_AUROC, gene_AUPRC, precision_at_50 ranked by the
priority-magnitude readout |q_signed| for the signed methods (Part H frozen
rule; identical to the real-data 02 priority ranking rule), and
direction_recovery_accuracy = sign agreement with the hidden direction vector
stored by the truth generator (step02A sim['true_mix'], step02B
sim['true_effect']) over y_true genes (reference rule). reference scored sim gene
metrics on the raw signed vector; the frozen final method definition ranks
by |q_signed| (protocol Part H), and the confirmatory run applies that rule to
the two signed-family methods (DUAL_HEAD and SIGNED_GENE_FIRST_PPR_SINGLE_HEAD
- same q_signed by construction, so their gene metrics are identical by
construction). This is a final definitional freeze applied to the new
confirmatory tables only; no baseline-4 file is altered.

The two-channel PPR solves of DUAL_HEAD / SIGNED_GENE_FIRST_PPR_SINGLE_HEAD /
UNSIGNED_GENE_PPR are bit-identical (same gene-only graph, same signed module
seed list, same per-channel restart vectors); they are therefore solved once
per replicate and shared across the three variants, and each variant's
convergence diagnostic rows report that shared solve metadata.

PPR frozen settings: alpha 0.35, tol 1e-10, max_iter 120, gene-only network
(within-community edges scale 0.35; no stage/module/risk nodes; the frozen
real-data builder convention leaves zero-degree rows zero). Module seed list =
top-180 (step02A) / top-220 (step02B) |mean_disease| genes (core/4
module-size convention).

Paired condition-wise evaluation (protocol I5): every primary metric, every
comparator, per condition and pooled over the design: mean_delta_dualhead_
minus_baseline on identical replicates, paired n, n_positive, n_negative,
n_ties, Wilcoxon signed-rank test (continuous metrics) or two-sided exact
binomial sign test (alignment_hit_top1, whose deltas are trinary),
BH-FDR within the declared families (benchmark x metric over the per-condition
rows only; pooled rows are reported without an FDR and excluded from the BH
family). Metric availability mirrors reference: true_stage_mass is defined only
for affinity-bearing family methods (affinity = softmax(zscore(cosine),1.0));
direction_recovery_accuracy only for methods with a signed gene vector
(DUAL_HEAD, SIGNED_GENE_FIRST_PPR_SINGLE_HEAD, DIRECT_NATIVE_PROJECTION);
undefined metric-comparator pairs are recorded as absent (no invented
readouts for archived comparators).

Outputs (results_confirmatory_simulation/01..08):
  01 replicate metrics (long)
  02 summary by benchmark x condition x method (means/sds)
  03 global method summary by benchmark x method
  04 method ranking by metric (per condition and global; rank_scope column)
  05 paired deltas (I5 long form with tests and BH within declared families)
  06 primary comparison table (head-to-head of DUAL_HEAD vs each comparator:
     pooled mean deltas, pooled p, per-condition win consistency)
  07 convergence diagnostics (family solves + core GENE_ONLY_PPR solves +
     baseline hetero PPR runtime rows, benchmark-tagged)
  08 numeric summary md (numbers only; no manuscript claims)
"""
import importlib.util
import os
import pathlib
from datetime import datetime

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = pathlib.Path(__file__).resolve().parents[2]
INPUT_ROOT = pathlib.Path(os.environ.get("NEUROTRACE_INPUT_ROOT", ROOT / "data" / "upstream"))
STAGE = ROOT / "data" / "processed" / "results"
RES = STAGE / "confirmatory_simulation"
RES.mkdir(parents=True, exist_ok=True)
LOGS = ROOT / "data" / "processed" / "logs"
LOGS.mkdir(parents=True, exist_ok=True)
SCRATCH = ROOT / "data" / "processed" / "scratch"
SCRATCH.mkdir(parents=True, exist_ok=True)
S3 = (INPUT_ROOT / "helpers" / "confirmatory_simulation.py")

# ------------------------------------------------------------------ freeze
ALPHA = 0.35
TOL = 1e-10
MAX_ITER = 120
GG_SCALE = 0.35

CONFIRM_SEED_A = 202609041
CONFIRM_SEED_B = 202609042
N_REP_CONFIRM = 200          # recorded before inspecting any result (I2 min 160)
MODULE_TOP_N = {"step02A": 180, "step02B": 220}

DUAL = "DUAL_HEAD_GENE_FIRST_PPR"
SIGNED = "SIGNED_GENE_FIRST_PPR_SINGLE_HEAD"
UNSIGNED = "UNSIGNED_GENE_PPR"
DIRECT = "DIRECT_NATIVE_PROJECTION"
GENE_ONLY_S3 = "GENE_ONLY_PPR_BASELINE"
HETERO = "HETERO_PPR_STANDARD"
STAGE_CORR = "STAGE_CORRELATION"
GRAPH_OT = "GRAPH_OT"
SIMPLE_OT = "SIMPLE_OT"
MEAN_SIG = "MEAN_SIGNATURE"
RANDOM = "RANDOM"
PROXY = "SIMULATION_PROXY"

CONFIRM_METHODS = [DUAL, SIGNED, UNSIGNED, DIRECT, GENE_ONLY_S3, HETERO,
                   STAGE_CORR, GRAPH_OT, SIMPLE_OT, MEAN_SIG, RANDOM, PROXY]

# protocol I4 metric names -> fixed snake columns (documented in 08)
PRIMARY_METRICS = ["alignment_hit_top1", "best_true_state_rank",
                   "true_stage_mass", "gene_auroc", "gene_auprc",
                   "precision_at_50", "direction_recovery_accuracy"]
METRIC_DIRECTION = {  # True = higher is better; False = lower is better
    "alignment_hit_top1": True, "best_true_state_rank": False,
    "true_stage_mass": True, "gene_auroc": True, "gene_auprc": True,
    "precision_at_50": True, "direction_recovery_accuracy": True}
STAGE_METRICS = ["alignment_hit_top1", "best_true_state_rank",
                 "true_stage_mass"]
GENE_METRICS = ["gene_auroc", "gene_auprc", "precision_at_50",
                "direction_recovery_accuracy"]

# methods whose row defines each metric (reference conventions)
STAGE_MASS_METHODS = {DUAL, SIGNED, UNSIGNED, DIRECT}
SIGNED_GENE_METHODS = {DUAL, SIGNED, DIRECT}
NO_SIGNED_STAGE_GENE = {"UNSIGNED_GENE_PPR"}  # family methods w/o signed gene

# protocol name -> implemented comparator name inside the archived/Stage1/3
# evaluators (verbatim), per design
ARCH_IMPL = {
    GRAPH_OT: "graph_ot", SIMPLE_OT: "simple_ot", MEAN_SIG: "mean_signature",
    RANDOM: "random",
    PROXY: {"step02A": "neurotrace_lite", "step02B": "neurotrace_full"}}
S1_IMPL = {HETERO: "hetero_ppr_standard", STAGE_CORR: "stage_correlation"}
S1_PPR_METHODS = {HETERO}


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def load_plain(path, modname):
    spec = importlib.util.spec_from_file_location(modname, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------ s3 framework
s3 = load_plain(S3, "s3confirmview")
archA, archB = s3.archA, s3.archB
s1A, s1B = s3.s1A, s3.s1B
CORR = s3.CORR
evaluate_core = s3.evaluate_core
row_from_scores = s3.row_from_scores
softmax_impl = s3.softmax_impl
auroc_impl = s3.auroc_impl
avg_precision_impl = s3.avg_precision_impl
precision_at_k_impl = s3.precision_at_k_impl
_community_edges = s3._community_edges
ppr_solve = s3.ppr_solve


def zscore(x, eps=1e-8):
    x = np.asarray(x, dtype=float)
    return (x - np.nanmean(x)) / (np.nanstd(x) + eps)


def cosine_vs_rows(vec, rows):
    nv = float(np.linalg.norm(vec))
    out = np.zeros(len(rows))
    for si in range(len(rows)):
        nz = float(np.linalg.norm(rows[si]))
        out[si] = (float(np.dot(vec, rows[si]) / (nv * nz))
                   if nv > 0 and nz > 0 else np.nan)
    return out


def build_gene_only_P(communities, n_g):
    i, j, v = _community_edges(communities, n_g)
    W = sp.coo_matrix((v * GG_SCALE, (i, j)),
                      shape=(n_g, n_g)).tocsr()
    outdeg = np.asarray(W.sum(axis=1)).ravel()
    n_dangle = int((outdeg <= 0).sum())
    inv = np.zeros_like(outdeg)
    inv[outdeg > 0] = 1.0 / outdeg[outdeg > 0]
    P = sp.diags(inv).dot(W).tocsr()
    return P, n_dangle


def family_channel_seed(w_mod, gene_pos, sign_positive, n_genes):
    sel = w_mod > 0 if sign_positive else w_mod < 0
    if int(sel.sum()) == 0:
        return None, 0
    mag = np.abs(w_mod[sel])
    s = np.zeros(n_genes)
    s[gene_pos[sel]] = mag / mag.sum()
    return s, int(sel.sum())


def direction_recovery(res_dir, sim):
    tv = sim.get("true_mix", sim.get("true_effect"))
    if tv is None or res_dir is None:
        return np.nan, 0
    tv = np.asarray(tv, dtype=float)
    sv = np.asarray(res_dir, dtype=float)
    y = np.asarray(sim["y_true_gene"]).astype(int)
    mask = (y == 1) & (np.abs(tv) > 1e-8) & (np.abs(sv) > 1e-12)
    n = int(mask.sum())
    if n == 0:
        return np.nan, 0
    return float(np.mean(np.sign(sv[mask]) == np.sign(tv[mask]))), n


def family_evaluate(variant, sim, ref, n_states, module_top_n):
    """One replicate, one dual-head-family variant. Channel solves are shared
    across the three PPR variants inside one call: q_pos/q_neg solved once per
    replicate here (the variant decides the readouts)."""
    n_g = len(ref["communities"])
    d = np.asarray(sim["mean_disease"], dtype=float)
    module_idx = np.argsort(-np.abs(d))[:module_top_n]
    module_idx = module_idx[np.abs(d[module_idx]) > 0]
    base = dict(stage_present=1, downstream_map=1, n_dangle=0,
                n_solves=0, pos_present=0, neg_present=0, pos_n=0, neg_n=0,
                pos_iter=np.nan, neg_iter=np.nan,
                pos_delta=np.nan, neg_delta=np.nan,
                n_iter=0, delta=0.0, converged=1)
    if variant == DIRECT:
        cs = cosine_vs_rows(d, ref["dev"])
        aff = softmax_impl(zscore(cs), temp=1.0)
        base.update(scores=cs, aff=aff, qg=np.abs(d), dir_vec=d.copy())
        return base
    P, n_dangle = build_gene_only_P(ref["communities"], n_g)
    sp_, sn_ = (family_channel_seed(d[module_idx], module_idx, True, n_g),
                family_channel_seed(d[module_idx], module_idx, False, n_g))
    q_pos = np.zeros(n_g)
    q_neg = np.zeros(n_g)
    iters, deltas, counts = [], [], []
    for ch, s_ in (("pos", sp_[0]), ("neg", sn_[0])):
        if s_ is None:
            continue
        q, n_iter, delta = ppr_solve(P, s_)
        iters.append(n_iter)
        deltas.append(float(delta))
        counts.append(ch)
        if ch == "pos":
            q_pos = q
        else:
            q_neg = q
    if not counts:
        cs = np.zeros(n_states)
        base.update(scores=cs, aff=np.full(n_states, 1.0 / n_states),
                    stage_present=0, qg=np.zeros(n_g), dir_vec=None,
                    n_dangle=n_dangle)
        return base
    q_mag = q_pos + q_neg
    q_sig = q_pos - q_neg
    if variant == UNSIGNED:
        qg, dir_vec = q_mag.copy(), None
    else:
        qg = np.abs(q_sig) if variant in (DUAL, SIGNED) else q_mag.copy()
        dir_vec = q_sig.copy() if variant in (DUAL, SIGNED) else None
    stage_vec = q_mag if variant in (DUAL, UNSIGNED) else q_sig
    cs = cosine_vs_rows(stage_vec, ref["dev"])
    aff = softmax_impl(zscore(cs), temp=1.0)
    pos_present, neg_present = sp_[0] is not None, sn_[0] is not None
    pos_iter = iters[0] if counts[0] == "pos" else (
        iters[1] if len(iters) > 1 else np.nan)
    neg_iter = (iters[1] if len(iters) > 1 and counts[1] == "neg" else
                (iters[0] if counts[0] == "neg" else np.nan))
    pos_delta = deltas[0] if counts[0] == "pos" else (
        deltas[1] if len(deltas) > 1 else np.nan)
    neg_delta = (deltas[1] if len(deltas) > 1 and counts[1] == "neg" else
                 (deltas[0] if counts[0] == "neg" else np.nan))
    base.update(scores=cs, aff=aff, qg=qg, dir_vec=dir_vec,
                n_dangle=n_dangle, n_solves=len(counts),
                pos_present=int(pos_present), neg_present=int(neg_present),
                pos_n=sp_[1], neg_n=sn_[1],
                pos_iter=pos_iter, neg_iter=neg_iter,
                pos_delta=pos_delta, neg_delta=neg_delta,
                n_iter=max(iters), delta=max(deltas),
                converged=int(all([q_delta < TOL for q_delta in deltas])))
    return base


def family_row(variant, rep_id, cond_name, sim, ref, res, n_states):
    y = np.asarray(sim["y_true_gene"]).astype(int)
    cs = res["scores"]
    order = np.argsort(-cs)
    pred = int(np.argmax(cs))
    aff = res["aff"]
    true_mass = float(aff[list(sim["true_states"])].sum())
    rank_p = int(np.where(order == sim["true_primary"])[0][0] + 1)
    rank_s = int(np.where(order == sim["true_secondary"])[0][0] + 1)
    transport = softmax_impl(cs, temp=0.75)
    compat_mass = float(transport[list(sim["true_states"])].sum())
    dir_acc, dir_n = direction_recovery(res["dir_vec"], sim)
    return dict(
        replicate=rep_id, condition=cond_name, method=variant,
        implementation_method=variant,
        true_primary_state=ref["states"][sim["true_primary"]],
        true_secondary_state=ref["states"][sim["true_secondary"]],
        predicted_state=ref["states"][pred],
        alignment_hit_top1=int(pred in sim["true_states"]),
        alignment_hit_top2=int(len(set(order[:2].tolist()).intersection(
            sim["true_states"])) > 0),
        best_true_state_rank=int(min(rank_p, rank_s)),
        true_transport_mass=float(compat_mass),
        gene_auroc=auroc_impl(y, res["qg"]),
        gene_auprc=avg_precision_impl(y, res["qg"]),
        precision_at_50=precision_at_k_impl(y, res["qg"], 50),
        true_stage_mass=float(true_mass), true_stage_mass_defined=1,
        stage_signal_present=res["stage_present"], downstream_stage_map=1,
        n_iter=res["n_iter"], converged=float(res["converged"]),
        direction_recovery_accuracy=dir_acc, direction_recovery_n=dir_n,
        method_origin="dual_head_gene_first_ppr")


def family_diag(variant, rep_id, cond_name, res, benchmark):
    return dict(benchmark=benchmark, method=variant, condition=cond_name,
                replicate=rep_id, n_ppr_solves=res["n_solves"],
                n_iter=res["n_iter"], delta=res["delta"],
                converged=res["converged"], n_dangle=res["n_dangle"],
                pos_channel_present=res["pos_present"],
                neg_channel_present=res["neg_present"],
                pos_n_seed_genes=res["pos_n"], neg_n_seed_genes=res["neg_n"],
                pos_n_iter=res["pos_iter"], neg_n_iter=res["neg_iter"],
                pos_delta=res["pos_delta"], neg_delta=res["neg_delta"],
                stage_signal_present=res["stage_present"],
                downstream_stage_map=res["downstream_map"],
                method_origin="dual_head_gene_first_ppr")


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full(len(p), np.nan)
    finite = ~np.isnan(p)
    if finite.sum() == 0:
        return out
    pv = p[finite]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.minimum(adj, 1.0)
    out[finite] = adj[np.argsort(order)]
    return out


def wilcoxon_paired(d, rng_stat):
    """Two-sided paired test on delta vector d (d = dual - baseline).
    Returns (p, test_name)."""
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d)]
    if len(d) == 0:
        return np.nan, "no_pairs"
    nz = d[np.abs(d) > 1e-12]
    if len(nz) == 0:
        return np.nan, "all_ties_no_test"
    try:
        from scipy import stats
    except Exception:
        return np.nan, "scipy_unavailable"
    try:
        stat, p = stats.wilcoxon(d, zero_method="wilcox")
        if not np.isfinite(p):
            raise ValueError
        return float(p), "wilcoxon_signed_rank_two_sided"
    except ValueError:
        p = np.nan
    # degenerate / all-ties-after-zero-drop alternative: exact sign test
    npos = int((nz > 0).sum())
    p = 2.0 * float(min(
        _binom_tail(len(nz), npos),
        _binom_tail(len(nz), len(nz) - npos)))
    return min(p, 1.0), "exact_sign_test_two_sided"


def _binom_tail(n, k):
    from math import comb
    return sum(comb(n, i) for i in range(k, n + 1)) / (2.0 ** n)


def sign_test_paired(d):
    d = np.asarray(d, dtype=float)
    d = d[np.isfinite(d)]
    nz = d[np.abs(d) > 1e-12]
    if len(nz) == 0:
        return np.nan, "all_ties_no_test"
    npos = int((nz > 0).sum())
    p = 2.0 * min(_binom_tail(len(nz), npos),
                  _binom_tail(len(nz), len(nz) - npos))
    return min(p, 1.0), "exact_sign_test_two_sided"


# ------------------------------------------------------------ replicate loop
def run_design(arch, s1, name, module_top_n, seed):
    rng = np.random.default_rng(seed)
    out_rows, diag_rows = [], []
    rep_global = 0
    total = len(arch.CONDITIONS) * N_REP_CONFIRM
    for cond_name, pars in arch.CONDITIONS.items():
        for rep in range(1, N_REP_CONFIRM + 1):
            rep_global += 1
            if name == "step02A":
                ref, sim = arch.simulate_one_replicate(rng, cond_name, pars)
            else:
                ref, sim = arch.simulate_replicate(rng, cond_name, pars)

            # archived comparator methods (verbatim evaluators)
            for m_proto, m_arch in [
                    (GRAPH_OT, ARCH_IMPL[GRAPH_OT]),
                    (SIMPLE_OT, ARCH_IMPL[SIMPLE_OT]),
                    (MEAN_SIG, ARCH_IMPL[MEAN_SIG]),
                    (RANDOM, ARCH_IMPL[RANDOM])]:
                row = (arch.evaluate_method(m_arch, rep_global, cond_name,
                                            sim, ref)
                       if name == "step02A"
                       else arch.evaluate(m_arch, rep_global, cond_name,
                                          sim, ref))
                row = dict(row)
                row["method"] = m_proto
                row["implementation_method"] = m_arch
                row["method_origin"] = "archived"
                if name == "step02A":
                    gs = arch.gene_scores_by_method(m_arch, sim, ref)
                    row["precision_at_50"] = precision_at_k_impl(
                        sim["y_true_gene"], gs, 50)
                out_rows.append(row)
            # archived archived NeuroTRACE proxy (historical-only label)
            m_arch = ARCH_IMPL[PROXY][name]
            row = (arch.evaluate_method(m_arch, rep_global, cond_name, sim,
                                        ref)
                   if name == "step02A"
                   else arch.evaluate(m_arch, rep_global, cond_name, sim,
                                      ref))
            row = dict(row)
            row["method"] = PROXY
            row["implementation_method"] = m_arch
            row["method_origin"] = "archived_proxy_historical_only"
            if name == "step02A":
                gs = arch.gene_scores_by_method(m_arch, sim, ref)
                row["precision_at_50"] = precision_at_k_impl(
                    sim["y_true_gene"], gs, 50)
            out_rows.append(row)

            # baseline PPR baselines (verbatim)
            for m_proto, m_s1 in S1_IMPL.items():
                row = (s1.evaluate_method(m_s1, rep_global, cond_name, sim,
                                          ref)
                       if name == "step02A"
                       else s1.evaluate(m_s1, rep_global, cond_name, sim,
                                        ref))
                row = dict(row)
                row["method"] = m_proto
                row["implementation_method"] = m_s1
                row["method_origin"] = "stage1_ppr_baseline"
                if name == "step02A":
                    _, gs = s1.score_pair(m_s1, sim, ref, rep_global,
                                          cond_name)
                    row["precision_at_50"] = precision_at_k_impl(
                        sim["y_true_gene"], gs, 50)
                out_rows.append(row)

            # core core GENE_ONLY_PPR replayed unchanged (protocol
            # comparator gene_only_ppr_stage3)
            res = evaluate_core("GENE_ONLY_PPR", rep_global, cond_name, sim,
                                ref, arch.N_STATES, module_top_n)
            row = row_from_scores("GENE_ONLY_PPR", rep_global, cond_name,
                                  sim, ref, res, arch.N_STATES)
            row = dict(row)
            row["method"] = GENE_ONLY_S3
            row["implementation_method"] = "GENE_ONLY_PPR"
            row["method_origin"] = "neurotrace_ppr_core"
            out_rows.append(row)
            diag_rows.append(dict(
                benchmark=name, method=GENE_ONLY_S3, condition=cond_name,
                replicate=rep_global, n_ppr_solves=1, n_iter=res["n_iter"],
                delta=res["delta"], converged=res["converged"],
                n_dangle=res["n_dangle"],
                pos_channel_present=np.nan, neg_channel_present=np.nan,
                pos_n_seed_genes=np.nan, neg_n_seed_genes=np.nan,
                pos_n_iter=np.nan, neg_n_iter=np.nan,
                pos_delta=np.nan, neg_delta=np.nan,
                stage_signal_present=res["stage_present"],
                downstream_stage_map=res["downstream_map"],
                method_origin="neurotrace_ppr_core"))

            # final dual-head family (one shared two-channel solve per
            # replicate for the three PPR variants)
            for variant in (DUAL, SIGNED, UNSIGNED, DIRECT):
                res_f = family_evaluate(variant, sim, ref, arch.N_STATES,
                                        module_top_n)
                out_rows.append(family_row(variant, rep_global, cond_name,
                                           sim, ref, res_f, arch.N_STATES))
                diag_rows.append(family_diag(variant, rep_global, cond_name,
                                             res_f, name))
            if rep_global % 100 == 0:
                print(f"[{now()}] {name} replicate {rep_global}/{total} "
                      f"(seed {seed})", flush=True)
    return out_rows, diag_rows, rep_global


# ------------------------------------------------------------ post summary
def fmt_markdown(df, cols, float_fmt="{:.4f}"):
    lines = ["| " + " | ".join(cols) + " |",
             "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            cells.append(float_fmt.format(v) if isinstance(v, float)
                         else str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def summarize(m):
    rows = []
    for (bench, cond, method), sub in m.groupby(["benchmark", "condition",
                                                 "method"]):
        row = {"benchmark": bench, "condition": cond, "method": method,
               "n_replicates": len(sub)}
        for metric in PRIMARY_METRICS:
            v = pd.to_numeric(sub[metric], errors="coerce").dropna()
            if metric == "true_stage_mass":
                keep = sub["true_stage_mass_defined"].fillna(0).astype(bool)
                v = v[keep.loc[v.index].values] if len(v) else v
                row["true_stage_mass_defined_n"] = int(len(v))
            if metric == "direction_recovery_accuracy":
                row["n_direction_defined"] = int(len(v))
            row[f"{metric}_mean"] = float(v.mean()) if len(v) else np.nan
            row[f"{metric}_sd"] = (float(v.std(ddof=1)) if len(v) > 1
                                   else np.nan)
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["benchmark", "condition", "method"])


def global_summary(m):
    rows = []
    for (bench, method), sub in m.groupby(["benchmark", "method"]):
        row = {"benchmark": bench, "method": method, "n": len(sub)}
        for metric in PRIMARY_METRICS:
            v = pd.to_numeric(sub[metric], errors="coerce").dropna()
            if metric == "true_stage_mass":
                keep = sub["true_stage_mass_defined"].fillna(0).astype(bool)
                v = v[keep.loc[v.index].values] if len(v) else v
                row["true_stage_mass_defined_n"] = int(len(v))
            if metric == "direction_recovery_accuracy":
                row["n_direction_defined"] = int(len(v))
            row[f"{metric}_mean"] = float(v.mean()) if len(v) else np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["benchmark", "method"])


def ranking(m, s):
    """Per-metric ranks within benchmark x condition (rank_scope=condition,
    ranks on the per-condition method means) and within benchmark over pooled
    replicates (rank_scope=global, ranks on pooled means computed directly
    from replicate rows). Metric direction respected (lower-is-better metrics
    ranked ascending)."""
    rows = []
    for (bench, cond), sub in s.groupby(["benchmark", "condition"]):
        for metric in PRIMARY_METRICS:
            if metric in ("true_stage_mass",):
                sel = sub[sub["true_stage_mass_defined_n"] > 0]
            elif metric == "direction_recovery_accuracy":
                sel = sub[sub["n_direction_defined"] > 0]
            else:
                sel = sub[sub[f"{metric}_mean"].notna()]
            if len(sel) == 0:
                continue
            ascending = not METRIC_DIRECTION[metric]
            ranked = sel.copy()
            ranked["rank_in_condition"] = ranked[f"{metric}_mean"].rank(
                ascending=ascending, method="min")
            for _, r in ranked.iterrows():
                rows.append(dict(benchmark=bench, condition=cond, metric=metric,
                                 method=r["method"],
                                 metric_mean=float(r[f"{metric}_mean"]),
                                 n_methods_ranked=len(ranked),
                                 rank_in_condition=int(r["rank_in_condition"]),
                                 rank_scope="condition"))
    for bench, sub in m.groupby("benchmark"):
        for metric in PRIMARY_METRICS:
            vals = {}
            for method, ss in sub.groupby("method"):
                v = pd.to_numeric(ss[metric], errors="coerce")
                if metric == "true_stage_mass":
                    v = v[ss["true_stage_mass_defined"].fillna(0).astype(
                        bool).values]
                elif metric == "direction_recovery_accuracy":
                    v = v[v.notna()]
                if len(v) == 0:
                    continue
                vals[method] = float(v.mean())
            if len(vals) < 2:
                continue
            ascending = not METRIC_DIRECTION[metric]
            order = sorted(vals, key=lambda k: vals[k], reverse=not ascending)
            rank_map = {k: order.index(k) + 1 for k in order}
            for method, mean_v in vals.items():
                rows.append(dict(benchmark=bench,
                                 condition="ALL_CONDITIONS", metric=metric,
                                 method=method, metric_mean=mean_v,
                                 n_methods_ranked=len(vals),
                                 rank_in_condition=rank_map[method],
                                 rank_scope="global"))
    return pd.DataFrame(rows).sort_values(
        ["benchmark", "rank_scope", "condition", "metric", "method"])


def paired_tests(m):
    """I5: paired condition-wise deltas DUAL minus each comparator on
    identical replicates, per condition and pooled; tests + BH within the
    declared (benchmark, metric) family over per-condition rows."""
    rows = []
    dual = m[m["method"] == DUAL].copy()
    for bench in ["step02A", "step02B"]:
        d_all = dual[dual["benchmark"] == bench]
        base_methods = [x for x in CONFIRM_METHODS if x != DUAL]
        for bl in base_methods:
            b = m[(m["benchmark"] == bench) & (m["method"] == bl)]
            merged = d_all.merge(b, on=["condition", "replicate"],
                                 suffixes=("_dual", "_base"))
            for metric in PRIMARY_METRICS:
                mc = f"{metric}_base"
                if mc not in merged.columns:
                    continue
                md = pd.to_numeric(merged[f"{metric}_dual"],
                                   errors="coerce")
                mb = pd.to_numeric(merged[mc], errors="coerce")
                if metric == "true_stage_mass":
                    ok = md.notna() & mb.notna() & merged[
                        "true_stage_mass_defined_dual"].fillna(0).astype(
                        bool) & merged["true_stage_mass_defined_base"].fillna(
                        0).astype(bool)
                elif metric == "direction_recovery_accuracy":
                    ok = md.notna() & mb.notna()
                else:
                    ok = md.notna() & mb.notna()
                if not bool(ok.any()):
                    continue
                conditions = sorted(merged["condition"].unique())
                for cond in conditions + ["ALL_CONDITIONS"]:
                    sub = merged[ok]
                    if cond != "ALL_CONDITIONS":
                        sub = sub[sub["condition"] == cond]
                    d = (md[sub.index] - mb[sub.index]).values
                    n_pairs = int(len(d))
                    if n_pairs == 0:
                        continue
                    if metric == "alignment_hit_top1":
                        p, test_name = sign_test_paired(d)
                    else:
                        p, test_name = wilcoxon_paired(d, None)
                    rows.append(dict(
                        benchmark=bench, condition=cond,
                        baseline_method=bl, metric=metric,
                        n_pairs=n_pairs,
                        mean_delta_dualhead_minus_baseline=float(d.mean()),
                        sd_delta=float(d.std(ddof=1)) if len(d) > 1 else np.nan,
                        n_positive=int((d > 1e-12).sum()),
                        n_negative=int((d < -1e-12).sum()),
                        n_ties=int((np.abs(d) <= 1e-12).sum()),
                        test_name=test_name, p_value=p,
                        metric_direction=("higher_is_better" if
                                          METRIC_DIRECTION[metric]
                                          else "lower_is_better"),
                        bh_family="benchmark_x_metric_per_condition_only"))
    df = pd.DataFrame(rows)
    # BH within (benchmark, metric): per-condition rows only
    for (bench, metric), idx in df.groupby(["benchmark", "metric"]).groups.items():
        cond_idx = idx[df.loc[idx, "condition"] != "ALL_CONDITIONS"]
        if len(cond_idx) == 0:
            continue
        df.loc[cond_idx, "fdr_BH_within_family"] = bh_fdr(
            df.loc[cond_idx, "p_value"].values)
    return df.sort_values(["benchmark", "condition", "baseline_method",
                           "metric"])


def primary_comparison(m):
    """06: pooled (whole-benchmark) head-to-head of DUAL_HEAD vs each
    comparator on every primary metric plus per-condition win consistency."""
    rows = []
    dual = m[m["method"] == DUAL].copy()
    for bench in ["step02A", "step02B"]:
        d_all = dual[dual["benchmark"] == bench]
        merged_all = None
        for bl in [x for x in CONFIRM_METHODS if x != DUAL]:
            b = m[(m["benchmark"] == bench) & (m["method"] == bl)]
            merged = d_all.merge(b, on=["condition", "replicate"],
                                 suffixes=("_dual", "_base"))
            row = dict(benchmark=bench, baseline_method=bl,
                       n_replicates=len(merged))
            for metric in PRIMARY_METRICS:
                mc = f"{metric}_base"
                if mc not in merged.columns:
                    row[f"{metric}_comparator_defined"] = 0
                    continue
                md = pd.to_numeric(merged[f"{metric}_dual"], errors="coerce")
                mb = pd.to_numeric(merged[mc], errors="coerce")
                if metric == "true_stage_mass":
                    ok = md.notna() & mb.notna() & merged[
                        "true_stage_mass_defined_dual"].fillna(0).astype(
                        bool) & merged["true_stage_mass_defined_base"].fillna(
                        0).astype(bool)
                else:
                    ok = md.notna() & mb.notna()
                row[f"{metric}_comparator_defined"] = int(ok.sum())
                if not bool(ok.any()):
                    row[f"{metric}_dual_mean"] = np.nan
                    row[f"{metric}_baseline_mean"] = np.nan
                    row[f"{metric}_mean_delta"] = np.nan
                    row[f"{metric}_pooled_p"] = np.nan
                    row[f"{metric}_n_cond_dual_better"] = np.nan
                    continue
                d = (md[ok] - mb[ok]).values
                row[f"{metric}_dual_mean"] = float(md[ok].mean())
                row[f"{metric}_baseline_mean"] = float(mb[ok].mean())
                row[f"{metric}_mean_delta"] = float(d.mean())
                if metric == "alignment_hit_top1":
                    p, _ = sign_test_paired(d)
                else:
                    p, _ = wilcoxon_paired(d, None)
                row[f"{metric}_pooled_p"] = p
                wins = 0
                for cond in sorted(merged["condition"].unique()):
                    idx = ok & (merged["condition"] == cond)
                    if idx.sum() == 0:
                        continue
                    dmean = float((md[idx] - mb[idx]).mean())
                    if METRIC_DIRECTION[metric]:
                        wins += int(dmean > 1e-12)
                    else:
                        wins += int(dmean < -1e-12)
                row[f"{metric}_n_cond_dual_better"] = wins
            rows.append(row)
    return pd.DataFrame(rows).sort_values(["benchmark", "baseline_method"])


def main():
    t0 = datetime.now()
    print(f"[{now()}] final Part H+I confirmatory simulation start")
    print(f"[{now()}] N_REP_CONFIRM = {N_REP_CONFIRM} per condition "
          f"(recorded before results are inspected; protocol I2 minimum "
          f"160); CONFIRM_SEED_A = {CONFIRM_SEED_A}; "
          f"CONFIRM_SEED_B = {CONFIRM_SEED_B}", flush=True)
    print(f"[{now()}] methods ({len(CONFIRM_METHODS)}): "
          f"{', '.join(CONFIRM_METHODS)}", flush=True)
    for mod, names in [(archA, ["simulate_one_replicate", "evaluate_method",
                                "gene_scores_by_method"]),
                       (archB, ["simulate_replicate", "evaluate"]),
                       (s1A, ["evaluate_method", "score_pair"]),
                       (s1B, ["evaluate", "score_pair"])]:
        for n in names:
            assert hasattr(mod, n), (mod.__name__, n)

    rows_a, diag_a, n_rep_a = run_design(archA, s1A, "step02A",
                                         MODULE_TOP_N["step02A"],
                                         CONFIRM_SEED_A)
    print(f"[{now()}] design step02A complete: {n_rep_a} replicates, "
          f"{len(rows_a)} method rows", flush=True)
    for r in rows_a:
        r["benchmark"] = "step02A"
    pd.DataFrame(rows_a).to_csv(RES / "_checkpoint_designA_metrics.tsv",
                                sep="\t", index=False)

    rows_b, diag_b, n_rep_b = run_design(archB, s1B, "step02B",
                                         MODULE_TOP_N["step02B"],
                                         CONFIRM_SEED_B)
    print(f"[{now()}] design step02B complete: {n_rep_b} replicates, "
          f"{len(rows_b)} method rows", flush=True)
    for r in rows_b:
        r["benchmark"] = "step02B"
    pd.DataFrame(rows_b).to_csv(RES / "_checkpoint_designB_metrics.tsv",
                                sep="\t", index=False)

    m = pd.DataFrame(rows_a + rows_b)
    # normalize missing true_stage_mass defined flag (archived/stage1 rows do
    # not carry the key; core GENE_ONLY_PPR_BASELINE rows carry defined=0)
    if "true_stage_mass_defined" not in m.columns:
        m["true_stage_mass_defined"] = 0
    m.to_csv(RES / "01_confirmatory_replicate_metrics.tsv", sep="\t",
             index=False)
    (RES / "_checkpoint_designA_metrics.tsv").unlink(missing_ok=True)
    (RES / "_checkpoint_designB_metrics.tsv").unlink(missing_ok=True)
    print(f"[{now()}] 01 written: {len(m)} replicate-method rows "
          f"(2 designs x {n_rep_a}/{n_rep_b} reps x {len(CONFIRM_METHODS)} "
          f"methods)")

    s = summarize(m)
    s.to_csv(RES / "02_confirmatory_summary_by_condition_method.tsv",
             sep="\t", index=False)
    g = global_summary(m)
    g.to_csv(RES / "03_confirmatory_global_method_summary.tsv", sep="\t",
             index=False)
    rk = ranking(m, s)
    rk.to_csv(RES / "04_confirmatory_method_ranking_by_metric.tsv",
              sep="\t", index=False)
    pw = paired_tests(m)
    pw.to_csv(RES / "05_confirmatory_paired_deltas.tsv", sep="\t",
              index=False)
    pc = primary_comparison(m)
    pc.to_csv(RES / "06_confirmatory_primary_comparison_table.tsv",
              sep="\t", index=False)

    # ---- 07 convergence diagnostics (family + core + stage1 runtime rows)
    fam_diag = pd.DataFrame(diag_a + diag_b)
    s1_rows = []
    for df_, bname in [(pd.DataFrame(s1A.PPR_RUNTIME_ROWS), "step02A"),
                       (pd.DataFrame(s1B.PPR_RUNTIME_ROWS), "step02B")]:
        if len(df_):
            df_ = df_.copy()
            df_["benchmark"] = bname
            df_["method_origin"] = "stage1_ppr_baseline"
            s1_rows.append(df_)
    diag_all = pd.concat([fam_diag] + s1_rows, ignore_index=True) \
        if s1_rows else fam_diag
    diag_all.to_csv(RES / "07_confirmatory_convergence_diagnostics.tsv",
                    sep="\t", index=False)
    fam_conv = fam_diag[fam_diag["method_origin"] ==
                        "dual_head_gene_first_ppr"]
    conv_all = bool((fam_conv["converged"] == 1).all())
    print(f"[{now()}] family solve rows {len(fam_conv)}; "
          f"all converged = {conv_all}")

    # ---- 08 summary md (numbers only)
    gA = g[g["benchmark"] == "step02A"].sort_values("method")
    gB = g[g["benchmark"] == "step02B"].sort_values("method")
    pw_all = pw[pw["condition"] == "ALL_CONDITIONS"]
    pw_al = pw_all[pw_all["metric"].isin(
        ["alignment_hit_top1", "true_stage_mass", "gene_auroc",
         "gene_auprc", "precision_at_50"])]
    txt = ["# final Part H+I summary - confirmatory held-out simulation "
           "with fresh seeds",
           "",
           f"Generated {now()} from results_confirmatory_simulation/01-07. "
           "Numbers only; interpretation and manuscript claims are deferred "
           "to the analysis.",
           "",
           "## Design",
           "",
           "- N_REP_CONFIRM = 200 replicates per condition (recorded before "
           "inspecting any result; protocol I2 minimum 160).",
           "- Fresh seeds not used in Stages 1-4: CONFIRM_SEED_A = "
           f"{CONFIRM_SEED_A} (step02A: 6 conditions x 200 = "
           f"{n_rep_a} replicates), CONFIRM_SEED_B = {CONFIRM_SEED_B} "
           f"(step02B: 6 hard conditions x 200 = {n_rep_b} replicates).",
           "- step02A conditions: baseline, high_batch, composition_shift, "
           "high_dropout, false_prior, combined_stress; step02B conditions: "
           "weak_signal, low_overlap, composition_confounded, high_dropout, "
           "false_prior_stress, combined_hard. Conditions in the archived "
           "dict order; one generator call per replicate (verbatim archived "
           "main-loop consumption pattern, fresh rng per design).",
           "- Generative settings unchanged from Stages 3/4 (1500 genes / 8 "
           "states / 30 communities / 3 cohorts step02A; 3000 genes / 10 "
           "states / 60 communities / 4 cohorts step02B).",
           "- Module seed list size for the PPR family: top-180 "
           "(step02A) / top-220 (step02B) |mean_disease| genes (core/4 "
           "convention).",
           "- PPR frozen settings: alpha 0.35, tol 1e-10, max_iter 120; "
           "gene-only network (within-community edges scale 0.35, no "
           "stage/module/risk nodes); two-channel restart from the signed "
           "module weights (per-channel sum-1 magnitudes).",
           "- DUAL_HEAD_GENE_FIRST_PPR: stage head q_magnitude = q_pos + "
           "q_neg -> magnitude cosine vs dev -> affinity "
           "softmax(zscore, 1.0); gene head q_signed = q_pos - q_neg, gene "
           "priority |q_signed|, direction = sign(q_signed). One method with "
           "two task-specific outputs (protocol Part H).",
           "- The two-channel PPR solves are bit-identical for "
           "DUAL_HEAD / SIGNED_GENE_FIRST_PPR_SINGLE_HEAD / UNSIGNED_GENE_PPR "
           "(same graph and same channel seeds); solved once per replicate "
           "and shared (each variant's 07 rows report that shared metadata).",
           "- Gene-metric ranking freeze: gene_AUROC/AUPRC/precision_at_50 "
           "rank by the priority magnitude |q_signed| for the two signed "
           "methods (protocol Part H; the real-data 02 priority rule). Stage "
           "4 scored its simulation gene metrics on the raw signed vector; "
           "the final frozen definition applies |q_signed| to the new "
           "confirmatory tables only (no baseline-4 file altered).",
           "- Methods (protocol I3; canonical method labels): " +
           ", ".join(CONFIRM_METHODS) +
           ". SIMULATION_PROXY = archived NeuroTRACE simulation "
           "proxy replayed unchanged (step02A neurotrace_lite, step02B "
           "neurotrace_full; NOT the final method, historical context only).",
           "- All comparators are replayed with their verbatim Stage "
           "1/3/archived evaluators; no method consumes the shared simulation "
           "rng; test truth (true states, y_true_gene, stored direction "
           "vector) is used only in evaluation.",
           "- Metric mapping (snake_case columns of 01): gene_AUROC = "
           "gene_auroc, gene_AUPRC = gene_auprc. Metric availability mirrors "
           "reference: true_stage_mass is defined only for the affinity-bearing "
           "family methods (DUAL/SIGNED/UNSIGNED/DIRECT; affinity = "
           "softmax(zscore(cosine), 1.0)); direction_recovery_accuracy only "
           "for methods with a signed gene vector "
           "(DUAL/SIGNED/DIRECT); undefined metric-comparator pairs are "
           "absent rather than filled by invented readouts.",
           "",
           "## Convergence (07)",
           "",
           f"- final family solve rows: {len(fam_conv)} "
           f"(2 designs x {n_rep_a}+{n_rep_b} replicates x 4 family "
           f"methods); all converged = {conv_all} "
           "(channel solves: converged = both channels delta < 1e-10 within "
           "max_iter 120).",
           f"- GENE_ONLY_PPR_BASELINE core solve rows and baseline "
           f"HETERO_PPR_STANDARD runtime rows are in 07 with method_origin "
           "labels; total diagnostic rows "
           f"{len(diag_all)}.",
           "",
           "## Global method summary (03)",
           "",
           "### step02A (means over 1200 replicates/method; metric direction: "
           "alignment/best_true_stage_rank are hit-rate / lower-better rank)",
           "",
           fmt_markdown(gA, ["method", "alignment_hit_top1_mean",
                             "best_true_state_rank_mean",
                             "true_stage_mass_mean",
                             "gene_auroc_mean", "gene_auprc_mean",
                             "precision_at_50_mean",
                             "direction_recovery_accuracy_mean"]),
           "",
           "### step02B",
           "",
           fmt_markdown(gB, ["method", "alignment_hit_top1_mean",
                             "best_true_state_rank_mean",
                             "true_stage_mass_mean",
                             "gene_auroc_mean", "gene_auprc_mean",
                             "precision_at_50_mean",
                             "direction_recovery_accuracy_mean"]),
           "",
           "## Key paired deltas (05; DUAL_HEAD minus comparator on identical "
           "replicates, pooled over the design; p from the paired test of "
           "the metric class)",
           "",
           fmt_markdown(pw_al.sort_values(["benchmark", "metric",
                                           "baseline_method"]),
                        ["benchmark", "baseline_method", "metric", "n_pairs",
                         "mean_delta_dualhead_minus_baseline",
                         "n_positive", "n_negative", "n_ties", "p_value"]),
           "",
           "## Primary comparison table (06)",
           "",
           "- 06 reports, per benchmark and comparator, pooled mean deltas "
           "with pooled paired p and the number of the 6 conditions in which "
           "the per-condition mean delta favors DUAL_HEAD "
           "(n_cond_dual_better) for every primary metric.",
           "",
           "## Test and FDR conventions (05)",
           "",
           "- Paired deltas are computed on identical replicates "
           "(condition x replicate merge).",
           "- alignment_hit_top1: two-sided exact binomial sign test on "
           "non-tied pairs (deltas are trinary); continuous metrics: "
           "two-sided Wilcoxon signed-rank (zero_method wilcox), falling "
           "back to the exact sign test if degenerate.",
           "- BH-FDR computed within the declared families benchmark x "
           "metric over the per-condition rows only "
           "(fdr_BH_within_family); pooled ALL_CONDITIONS rows are reported "
           "for reference without an FDR and excluded from the BH family.",
           "- 09_dualhead_pareto_summary.tsv (Part J, domain-level "
           "better/tied/worse) is produced from 01-06 by the Part J step.",
           "",
           "## Caveats",
           "",
           "- DUAL_HEAD and UNSIGNED_GENE_PPR share the magnitude stage "
           "readout by construction (identical q_pos + q_neg stage head); "
           "DUAL_HEAD and SIGNED_GENE_FIRST_PPR_SINGLE_HEAD share the signed "
           "gene output by construction (identical |q_signed| gene head). "
           "Their paired comparisons on those shared readouts are expected "
           "ties at machine precision; the informative contrasts are the "
           "cross-head ones and the comparisons against the other methods.",
           "- true_stage_mass for the family uses affinity = "
           "softmax(zscore(cosine), 1.0) (real-data Part D readout); "
           "true_transport_mass (informational, present in 01) uses the "
           "archived softmax(scores, 0.75) convention.",
           "- direction_recovery_accuracy compares sign(q_signed) (or the "
           "direct native vector) with the truth generator's stored "
           "direction vector (step02A sim['true_mix'], step02B "
           "sim['true_effect']) over y_true genes; undefined for methods "
           "without a signed gene score.",
           ]
    (RES / "08_confirmatory_simulation_summary.md").write_text("\n".join(txt))
    print(f"[{now()}] 08 written")
    print("\n[post] global step02A:")
    print(gA[["method", "alignment_hit_top1_mean", "best_true_state_rank_mean",
              "true_stage_mass_mean", "gene_auroc_mean", "gene_auprc_mean",
              "precision_at_50_mean",
              "direction_recovery_accuracy_mean"]].to_string(index=False))
    print("\n[post] global step02B:")
    print(gB[["method", "alignment_hit_top1_mean", "best_true_state_rank_mean",
              "true_stage_mass_mean", "gene_auroc_mean", "gene_auprc_mean",
              "precision_at_50_mean",
              "direction_recovery_accuracy_mean"]].to_string(index=False))
    print("\n[post] paired deltas (ALL_CONDITIONS, stage metrics):")
    print(pw_all[pw_all["metric"].isin(["alignment_hit_top1",
                                        "best_true_state_rank",
                                        "true_stage_mass"])].to_string(
                                        index=False))
    print("\n[post] paired deltas (ALL_CONDITIONS, gene metrics):")
    print(pw_all[pw_all["metric"].isin(["gene_auroc", "gene_auprc",
                                        "precision_at_50",
                                        "direction_recovery_accuracy"
                                        ])].to_string(index=False))
    print(f"[{now()}] elapsed {(datetime.now() - t0).total_seconds():.0f} s")
    print("[post] files written 01-08")


if __name__ == "__main__":
    main()
