#!/usr/bin/env python3
"""final Parts B+C+D+F: dual-head gene-first PPR on real frozen graph.

Fixed final dual-head framework (protocol Part B, no alteration allowed):
  positive restart  = abs(native weight) for weight > 0, normalized to sum 1
  negative restart  = abs(native weight) for weight < 0, normalized to sum 1
  q_pos / q_neg     = PPR on the gene-only graph (official
                      gene_gene_embedding_knn edges only, alpha 0.35,
                      tol 1e-10, max_iter 120)
  q_signed          = q_pos - q_neg
  q_magnitude       = q_pos + q_neg

  Developmental localization head  : q_magnitude
      primary stage rank = descending magnitude_cosine vs the exact processed
      six-stage BrainSpan profiles (z-softmax affinity for visualization only)
      localization_contrast = mean(magnitude_cosine, pre) - mean(post)
  Gene-priority head              : abs(q_signed)
  Gene direction                  : sign(q_signed)
  External disease score weight   : q_signed
  Signed developmental orientation: q_signed only, labeled
      SIGNED_DEVELOPMENTAL_ORIENTATION_SECONDARY (never developmental
      timing/localization).

No PPR is re-solved here: the per-gene q_pos/q_neg vectors for every
program x top_n in {50,100,200,500} are read from the frozen upstream output
results_realdata/01_signed_ppr_gene_scores.tsv.gz (already prepared upstream).
This script performs the Part C algebraic/sign identity check, the Part D
magnitude-head developmental localization, and the Part F secondary signed
orientation table.

Outputs:
  method/method_definition.tsv   (Part B)
  method/method_definition.md       (Part B)
  method/channel_identity.tsv       (Part C)
  results_realdata/01_dualhead_magnitude_stage_localization_long.tsv        (D)
  results_realdata/02_dualhead_magnitude_stage_affinity_matrix.tsv          (D)
  results_realdata/03_dualhead_magnitude_top_stage_calls.tsv                (D)
  results_realdata/04_dualhead_magnitude_prenatal_postnatal_contrast.tsv    (D)
  results_realdata/05_dualhead_pos_neg_channel_stage_localization.tsv       (D)
  results_realdata/06_dualhead_vs_reference_signed_orientation.tsv             (D)
  results_realdata/10_secondary_signed_developmental_orientation.tsv        (F)
"""
import importlib.util
import os
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
INPUT_ROOT = pathlib.Path(os.environ.get("NEUROTRACE_INPUT_ROOT", ROOT / "data" / "upstream"))
GENE_INPUTS = INPUT_ROOT / "gene_first_ppr_inputs"
RESULTS_ROOT = ROOT / "data" / "processed" / "results"
S4_RES = GENE_INPUTS / "realdata"
RES = RESULTS_ROOT / "realdata"
METHOD_DIR = RESULTS_ROOT / "method"
RES.mkdir(parents=True, exist_ok=True)
METHOD_DIR.mkdir(parents=True, exist_ok=True)
FROZEN = INPUT_ROOT / "gene_graph" / "results"
STEP03 = INPUT_ROOT / "developmental_reference"

# Shared stage-order and affinity helpers supplied with prepared inputs.
spec = importlib.util.spec_from_file_location(
    "step05b_driver",
    INPUT_ROOT / "helpers" / "stage_profiles.py")
step05b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(step05b)
STAGE_ORDER = step05b.STAGE_ORDER
zscore = step05b.zscore
softmax = step05b.softmax

ALPHA = 0.35
TOL = 1e-10
MAX_ITER = 120
PROGRAMS = ["NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed"]
TOPS = [50, 100, 200, 500]


def clean_gene(x):
    if pd.isna(x):
        return ""
    return str(x).strip().upper()


def stage_cosines(vec, Z):
    """Raw cosine of vec vs each of the six stage profile rows (NA when the
    vector has no mass)."""
    nv = np.linalg.norm(vec)
    out = []
    for si in range(len(STAGE_ORDER)):
        nz = np.linalg.norm(Z[si])
        out.append(float(np.dot(vec, Z[si]) / (nv * nz))
                   if nv > 0 and nz > 0 else np.nan)
    return np.array(out)


def load_inputs():
    """Mirror of upstream load_inputs (gene universe order + stage matrix Z),
    without the (unneeded) transition-graph build: q vectors are read from the
    frozen upstream output file instead of re-solved."""
    nodes = pd.read_csv(FROZEN / "09_step04C_native_graph_nodes.tsv", sep="\t")
    prof = pd.read_csv(STEP03 / "04_step03A_brainspan_stage_gene_profiles.tsv.gz",
                       sep="\t")
    gene_df = nodes[nodes["node_type"] == "gene"].copy()
    gene_df["sym"] = gene_df["feature_id"].map(clean_gene)
    gene_df = gene_df.drop_duplicates("sym")
    universe = gene_df["sym"].astype(str).tolist()
    sym_to_idx = {s: i for i, s in enumerate(universe)}
    assert len(universe) == 5000, len(universe)
    pp = prof.copy()
    pp["sym"] = pp["gene_key"].map(clean_gene)
    Z = np.zeros((len(STAGE_ORDER), len(universe)))
    for si, st in enumerate(STAGE_ORDER):
        sub = pp[pp["state_id"] == st].set_index("sym")["stage_profile_z"]
        for g, v in sub.items():
            if g in sym_to_idx:
                Z[si, sym_to_idx[g]] = v
    assert (np.abs(Z.sum(axis=1)) > 0).all()
    return dict(universe=universe, sym_to_idx=sym_to_idx, Z=Z)


def module_block(ntm, program, top_n, sym_to_idx):
    """Signed weight vector + per-sign restart statistics over the universe."""
    rows = ntm[(ntm["program"] == program) & (ntm["top_n"] == top_n)]
    w = np.zeros(len(sym_to_idx))
    n_mapped = 0
    for _, r in rows.iterrows():
        g = clean_gene(r["gene_symbol_fixed"])
        if g in sym_to_idx:
            w[sym_to_idx[g]] = float(r["weight"])
            n_mapped += 1
    pos_genes = int((w > 0).sum())
    neg_genes = int((w < 0).sum())
    mass_pos = float(w[w > 0].sum())
    mass_neg = float(-w[w < 0].sum())
    return dict(w=w, n_module=len(rows), n_mapped=n_mapped,
                n_pos=pos_genes, n_neg=neg_genes,
                mass_pos=mass_pos, mass_neg=mass_neg)


def load_q_vectors(df1, universe):
    """Return {(program, top_n): (q_pos, q_neg)} aligned to `universe` order."""
    out = {}
    for (program, top_n), grp in df1.groupby(["program", "top_n"]):
        m = grp.set_index("feature_id").reindex(universe)
        if m["q_pos"].isna().any() or m["q_neg"].isna().any():
            raise RuntimeError(f"{program} top{top_n}: universe/gene mismatch")
        out[(program, int(top_n))] = (
            np.asarray(m["q_pos"], dtype=float),
            np.asarray(m["q_neg"], dtype=float))
    return out


def projection_identity(c_signed, c_mag, atol=1e-9):
    """Classify the signed vs magnitude six-cosine projections.

    SAME_AS_MAGNITUDE          : signed cosine vector == magnitude cosine vector
    EXACT_SIGN_INVERSE_OF_MAGNITUDE : signed == -magnitude (elementwise)
    NEITHER                    : mixed channel, neither identity holds
    """
    d_same = float(np.max(np.abs(np.asarray(c_signed) - np.asarray(c_mag))))
    d_inv = float(np.max(np.abs(np.asarray(c_signed) + np.asarray(c_mag))))
    if d_same <= atol:
        return "SAME_AS_MAGNITUDE", d_same, d_inv
    if d_inv <= atol:
        return "EXACT_SIGN_INVERSE_OF_MAGNITUDE", d_same, d_inv
    return "NEITHER", d_same, d_inv


def main():
    inp = load_inputs()
    universe, sym_to_idx, Z = (inp["universe"], inp["sym_to_idx"], inp["Z"])
    print(f"[final] universe {len(universe)}; Z rows {len(STAGE_ORDER)}; "
          f"Z cols {Z.shape[1]}")

    # ---------- frozen upstream channel vectors ----------
    df1 = pd.read_csv(S4_RES / "01_signed_ppr_gene_scores.tsv.gz", sep="\t",
                      compression="gzip")
    grp_sizes = df1.groupby(["program", "top_n"]).size()
    assert len(grp_sizes) == 12 and (grp_sizes == 5000).all(), \
        f"unexpected 01 layout {len(grp_sizes)}"
    # file-consistency: recomputed signed/magnitude == stored columns
    re_signed = (df1["q_pos"] - df1["q_neg"]).to_numpy()
    re_mag = (df1["q_pos"] + df1["q_neg"]).to_numpy()
    e_s = float(np.max(np.abs(re_signed - df1["signed_score"].to_numpy())))
    e_m = float(np.max(np.abs(re_mag - df1["magnitude_score"].to_numpy())))
    assert e_s < 1e-10 and e_m < 1e-10, (e_s, e_m)
    print(f"[final] 01 consistency: max|signed recomputed-stored|={e_s:.2e}, "
          f"max|magnitude recomputed-stored|={e_m:.2e}")

    qvec = load_q_vectors(df1, universe)
    ntm = pd.read_csv(STEP03 /
                      "24_step03C2_neurotrace_native_module_weights_symbol.tsv",
                      sep="\t")

    # ---------- PART B: frozen method definition (02 tsv, 03 md) ----------
    b_rows = [
        ("ppr_graph", "node_types", "gene only (no module/stage/risk nodes)"),
        ("ppr_graph", "edges",
         "official frozen gene_gene_embedding_knn edges only "
         "(no module-stage/stage-gene/risk-gene edges)"),
        ("ppr_graph", "source",
         "gene_graph/step04_graph_construction/results/"
         "10_step04C_native_graph_edges.tsv.gz (identical to upstream)"),
        ("ppr_graph", "alpha", "0.35"),
        ("ppr_graph", "tol", "1e-10"),
        ("ppr_graph", "max_iter", "120"),
        ("ppr_graph", "no_direct_stage_blend", "yes (projection is post-hoc "
         "cosine onto processed stage profiles only)"),
        ("channel", "positive_restart",
         "abs(native_weight) for weight > 0, normalized to sum 1"),
        ("channel", "negative_restart",
         "abs(native_weight) for weight < 0, normalized to sum 1"),
        ("channel", "q_pos", "PPR restart from the positive channel"),
        ("channel", "q_neg", "PPR restart from the negative channel"),
        ("head", "q_signed", "q_pos - q_neg"),
        ("head", "q_magnitude", "q_pos + q_neg"),
        ("head", "developmental_localization",
         "q_magnitude; primary stage rank = descending magnitude_cosine "
         "against the exact processed six-stage BrainSpan profiles"),
        ("head", "localization_contrast",
         "mean(magnitude_cosine early/mid/late prenatal) - "
         "mean(magnitude_cosine childhood/adolescence/adulthood)"),
        ("head", "gene_priority", "abs(q_signed)"),
        ("head", "gene_direction", "sign(q_signed)"),
        ("head", "external_disease_score_weight", "q_signed"),
        ("head", "signed_orientation_label",
         "SIGNED_DEVELOPMENTAL_ORIENTATION_SECONDARY (never developmental "
         "timing/localization)"),
        ("affinity", "transform",
         "softmax(zscore(magnitude_cosine across six stages)), temperature "
         "1.0, no fitted temperature; visualization only"),
        ("affinity", "profiles",
         "processed stage_profile_z matrix on the shared 5000-gene universe "
         "(04_step03A_brainspan_stage_gene_profiles.tsv.gz, exact upstream)"),
        ("nonnegotiable", "no_tuning",
         "freeze all definitions before confirmatory simulation; no parameter "
         "selection from confirmatory runs"),
        ("nonnegotiable", "no_architecture_change",
         "no re-derivation, no architecture alteration based on results"),
        ("nonnegotiable", "no_new_downloads",
         "no new biological data downloads during the final run"),
        ("nonnegotiable", "no_manuscript_rewrite",
         "no rewrite of the manuscript; frozen Figures 1-6 untouched"),
    ]
    df_b = pd.DataFrame(b_rows, columns=["section", "component",
                                         "frozen_definition"])
    df_b.to_csv(METHOD_DIR / "method_definition.tsv", sep="\t",
                index=False)

    md_b = f"""# final frozen method definition: dual-head gene-first PPR

Freeze date: 2026-09-04 (freeze precedes all confirmatory simulation; no
parameter or architecture selection is made from final results).

## 1. PPR graph and parameters (verbatim B1)

- Graph: the upstream gene-only PPR graph, exactly: official frozen
  gene_gene_embedding_knn edges only; gene nodes only.
- alpha = 0.35, tol = 1e-10, max_iter = 120.
- Excluded from the PPR graph: module nodes; stage nodes; risk nodes;
  module-stage edges; stage-gene edges; risk-gene edges.
- No direct stage blend: developmental readout is a post-hoc cosine
  projection onto the processed six-stage BrainSpan profiles; stage nodes
  never participate in the diffusion.

## 2. Restart channels (verbatim B2)

For every module/top_n:

    positive channel restart  = abs(native_weight) for weight > 0
    negative channel restart  = abs(native_weight) for weight < 0

Each non-empty channel is normalized independently to sum to 1, then:

    q_pos
    q_neg
    q_signed    = q_pos - q_neg
    q_magnitude = q_pos + q_neg

## 3. Final frozen outputs (verbatim B3)

| purpose | quantity |
|---|---|
| developmental localization head | q_magnitude (top stage = argmax magnitude_cosine) |
| gene-priority head | abs(q_signed) |
| gene direction | sign(q_signed) |
| external disease score weight | q_signed |
| secondary signed developmental orientation | q_signed, labeled SIGNED_DEVELOPMENTAL_ORIENTATION_SECONDARY |

The secondary signed orientation is never called "developmental timing" or
"developmental localization".

## 4. Provenance notes

- q_pos/q_neg per gene were solved in upstream (identical graph, alpha 0.35,
  tol 1e-10, max_iter 120; convergence diagnostics in the upstream file
  03_signed_ppr_convergence_diagnostics.tsv: every channel converged) and are
  stored frozen per (program, top_n, gene) in the upstream file
  results_realdata/01_signed_ppr_gene_scores.tsv.gz. final reads those
  frozen vectors (verified above: recomputed q_signed/q_magnitude match the
  stored signed_score/magnitude_score columns to < 1e-10); no re-solve, no
  tuning.
- Affinity = softmax(zscore(magnitude_cosine across the six stages)) with
  temperature 1.0, no fitted temperature, for visualization only.
- Companion identity table: method/channel_identity.tsv.
"""
    (METHOD_DIR / "method_definition.md").write_text(md_b)
    print(f"[final] Part B: 02 tsv rows {len(df_b)}; 03 md bytes "
          f"{len(md_b)}")

    # ---------- PART C: sign-channel identity check (04) ----------
    c_rows = []
    cos_cache = {}  # (program, top_n, kind) -> six-cosine array
    for program in PROGRAMS:
        for top_n in TOPS:
            mb = module_block(ntm, program, top_n, sym_to_idx)
            qp, qn = qvec[(program, top_n)]
            qs = qp - qn
            qm = qp + qn
            cos_s = stage_cosines(qs, Z)
            cos_m = stage_cosines(qm, Z)
            cos_cache[(program, top_n, "signed")] = cos_s
            cos_cache[(program, top_n, "magnitude")] = cos_m
            assert np.isfinite(cos_s).all() and np.isfinite(cos_m).all(), \
                f"{program} top{top_n}: non-finite cosine"
            if mb["n_pos"] > 0 and mb["n_neg"] == 0:
                chan_cls, exp_ident = "all_positive", "q_signed==q_magnitude"
            elif mb["n_neg"] > 0 and mb["n_pos"] == 0:
                chan_cls, exp_ident = "all_negative", "q_signed==-q_magnitude"
            else:
                chan_cls, exp_ident = "mixed", "neither_expected"
            err_same = float(np.max(np.abs(qs - qm)))
            err_inv = float(np.max(np.abs(qs + qm)))
            if err_same < 1e-10:
                best_fit = "+1 (q_signed == q_magnitude)"
            elif err_inv < 1e-10:
                best_fit = "-1 (q_signed == -q_magnitude)"
            else:
                best_fit = "neither"
            pcls, pe_same, pe_inv = projection_identity(cos_s, cos_m)
            # top-stage/contrast context from this script's projection
            order = sorted(range(6), key=lambda i: -cos_s[i])
            pre_s = float(cos_s[:3].mean())
            post_s = float(cos_s[3:].mean())
            order_m = sorted(range(6), key=lambda i: -cos_m[i])
            pre_m = float(cos_m[:3].mean())
            post_m = float(cos_m[3:].mean())
            c_rows.append(dict(
                program=program, top_n=top_n,
                n_positive_restart_genes=mb["n_pos"],
                n_negative_restart_genes=mb["n_neg"],
                positive_channel_present=int(mb["n_pos"] > 0),
                negative_channel_present=int(mb["n_neg"] > 0),
                channel_class=chan_cls,
                expected_vector_identity=exp_ident,
                max_abs_error_q_signed_minus_q_magnitude=err_same,
                max_abs_error_q_signed_plus_q_magnitude=err_inv,
                vector_best_fit_sign=best_fit,
                signed_projection_vs_magnitude_max_abs_error=pe_same,
                signed_projection_vs_negated_magnitude_max_abs_error=pe_inv,
                projection_identity_class=pcls,
                signed_top_stage=STAGE_ORDER[order[0]],
                signed_prenatal_minus_postnatal=pre_s - post_s,
                magnitude_top_stage=STAGE_ORDER[order_m[0]],
                magnitude_localization_contrast=pre_m - post_m))
    df_c = pd.DataFrame(c_rows)
    df_c.to_csv(METHOD_DIR / "channel_identity.tsv", sep="\t",
                index=False)
    print(f"[final] Part C: 04 rows {len(df_c)}")
    print(df_c[df_c.top_n.isin([200, 500])].to_string(index=False))

    # ---------- PART D: magnitude-head localization (01-06) ----------
    long_rows, wide_rows, top_rows, contrast_rows, chan_rows = [], [], [], [], []
    for program in PROGRAMS:
        for top_n in TOPS:
            cos_m = cos_cache[(program, top_n, "magnitude")]
            zc = zscore(cos_m)
            aff = softmax(zc, temp=1.0)
            order = sorted(range(6), key=lambda i: -cos_m[i])
            pre = float(cos_m[:3].mean())
            post = float(cos_m[3:].mean())
            contr = pre - post
            mb = module_block(ntm, program, top_n, sym_to_idx)
            top_rows.append(dict(
                program=program, top_n=top_n,
                top_stage=STAGE_ORDER[order[0]],
                top_stage_magnitude_cosine=float(cos_m[order[0]]),
                second_stage=STAGE_ORDER[order[1]],
                second_stage_magnitude_cosine=float(cos_m[order[1]]),
                top_stage_margin=float(cos_m[order[0]] - cos_m[order[1]]),
                localization_contrast=contr,
                prenatal_or_postnatal=("prenatal" if contr >= 0
                                       else "postnatal"),
                affinity_prenatal=float(aff[:3].sum()),
                affinity_postnatal=float(aff[3:].sum()),
                n_module_genes=mb["n_module"],
                n_mapped_to_gene_network=mb["n_mapped"]))
            contrast_rows.append(dict(
                program=program, top_n=top_n,
                prenatal_localization=pre, postnatal_localization=post,
                localization_contrast=contr,
                prenatal_or_postnatal=("prenatal" if contr >= 0
                                       else "postnatal")))
            for si, s in enumerate(STAGE_ORDER):
                long_rows.append(dict(
                    program=program, top_n=top_n, state_id=s,
                    stage_order=si + 1,
                    magnitude_cosine=float(cos_m[si]),
                    cosine_z=float(zc[si]),
                    stage_affinity_z_softmax=float(aff[si]),
                    stage_rank=int(order.index(si) + 1),
                    top_stage=STAGE_ORDER[order[0]],
                    localization_contrast=contr))
            # ---- channel-level (D3): q_pos and q_neg separately ----
            for ch_name, qch in (("positive", qvec[(program, top_n)][0]),
                                 ("negative", qvec[(program, top_n)][1])):
                present = float(np.linalg.norm(qch)) > 0
                if present:
                    cch = stage_cosines(qch, Z)
                    assert np.isfinite(cch).all()
                    och = sorted(range(6), key=lambda i: -cch[i])
                    ch_pre = float(cch[:3].mean())
                    ch_post = float(cch[3:].mean())
                    ch_contr = ch_pre - ch_post
                    chan_rows.append(dict(
                        program=program, top_n=top_n, channel=ch_name,
                        channel_present=1,
                        n_restart_genes=(mb["n_pos"] if ch_name == "positive"
                                         else mb["n_neg"]),
                        channel_mass_before_normalization=(
                            mb["mass_pos"] if ch_name == "positive"
                            else mb["mass_neg"]),
                        channel_mass_share_of_total_abs_weight=(
                            mb["mass_pos"] / (mb["mass_pos"] + mb["mass_neg"])
                            if (mb["mass_pos"] + mb["mass_neg"]) > 0
                            and ch_name == "positive"
                            else (mb["mass_neg"] / (mb["mass_pos"] + mb["mass_neg"])
                                  if (mb["mass_pos"] + mb["mass_neg"]) > 0
                                  else np.nan)),
                        top_stage=STAGE_ORDER[och[0]],
                        top_stage_channel_cosine=float(cch[och[0]]),
                        prenatal_localization=ch_pre,
                        postnatal_localization=ch_post,
                        localization_contrast=ch_contr,
                        prenatal_or_postnatal=("prenatal" if ch_contr >= 0
                                               else "postnatal")))
                else:
                    chan_rows.append(dict(
                        program=program, top_n=top_n, channel=ch_name,
                        channel_present=0, n_restart_genes=0,
                        channel_mass_before_normalization=0.0,
                        channel_mass_share_of_total_abs_weight=0.0,
                        top_stage="NA", top_stage_channel_cosine=np.nan,
                        prenatal_localization=np.nan,
                        postnatal_localization=np.nan,
                        localization_contrast=np.nan,
                        prenatal_or_postnatal="NA"))
    df01 = pd.DataFrame(long_rows)
    df01.to_csv(RES / "01_dualhead_magnitude_stage_localization_long.tsv",
                sep="\t", index=False)
    # 02 wide affinity/cosine matrix
    wide_cols = []
    for program in PROGRAMS:
        for top_n in TOPS:
            sub = [x for x in long_rows if x["program"] == program and
                   x["top_n"] == top_n]
            cosv = [next(x["magnitude_cosine"] for x in sub
                         if x["state_id"] == s) for s in STAGE_ORDER]
            affv = [next(x["stage_affinity_z_softmax"] for x in sub
                         if x["state_id"] == s) for s in STAGE_ORDER]
            row = dict(program=program, top_n=top_n)
            row.update({f"cosine_{s}": cosv[i] for i, s in enumerate(STAGE_ORDER)})
            row.update({f"affinity_{s}": affv[i] for i, s in enumerate(STAGE_ORDER)})
            wide_cols.append(row)
    df02 = pd.DataFrame(wide_cols)
    cols = ["program", "top_n"] + [f"cosine_{s}" for s in STAGE_ORDER] + \
           [f"affinity_{s}" for s in STAGE_ORDER]
    df02[cols].to_csv(RES / "02_dualhead_magnitude_stage_affinity_matrix.tsv",
                      sep="\t", index=False)
    df03 = pd.DataFrame(top_rows)
    df03.to_csv(RES / "03_dualhead_magnitude_top_stage_calls.tsv", sep="\t",
                index=False)
    df04 = pd.DataFrame(contrast_rows)
    df04.to_csv(RES / "04_dualhead_magnitude_prenatal_postnatal_contrast.tsv",
                sep="\t", index=False)
    df05 = pd.DataFrame(chan_rows)
    df05.to_csv(RES / "05_dualhead_pos_neg_channel_stage_localization.tsv",
                sep="\t", index=False)

    # ---- 06: magnitude localization vs Stage-4 signed orientation ----
    s4_12 = pd.read_csv(S4_RES / "12_gene_first_top_stage_calls.tsv", sep="\t")
    s4_13 = pd.read_csv(S4_RES / "13_gene_first_prenatal_postnatal_contrast.tsv",
                        sep="\t")
    m_rows = []
    for program in PROGRAMS:
        for top_n in TOPS:
            s4t = s4_12[(s4_12["program"] == program) &
                        (s4_12["top_n"] == top_n)].iloc[0]
            s4c = s4_13[(s4_13["program"] == program) &
                        (s4_13["top_n"] == top_n)].iloc[0]
            ct = df03[(df03["program"] == program) &
                      (df03["top_n"] == top_n)].iloc[0]
            cr = df04[(df04["program"] == program) &
                      (df04["top_n"] == top_n)].iloc[0]
            crec = df_c[(df_c["program"] == program) &
                        (df_c["top_n"] == top_n)].iloc[0]
            pcls = crec["projection_identity_class"]
            m_rows.append(dict(
                program=program, top_n=top_n,
                reference_signed_top_stage=str(s4t["top_stage"]),
                reference_signed_prenatal_minus_postnatal=float(
                    s4c["prenatal_minus_postnatal"]),
                reference_signed_window=str(s4c["prenatal_or_postnatal"]),
                reference_signed_source=str(S4_RES /
                                         "12_gene_first_top_stage_calls.tsv"),
                magnitude_top_stage=str(ct["top_stage"]),
                magnitude_localization_contrast=float(
                    cr["localization_contrast"]),
                magnitude_window=str(cr["prenatal_or_postnatal"]),
                relation_to_reference_signed_orientation=pcls,
                delta_contrast_magnitude_minus_signed=(
                    float(cr["localization_contrast"]) -
                    float(s4c["prenatal_minus_postnatal"]))))
    df06 = pd.DataFrame(m_rows)
    df06.to_csv(RES / "06_dualhead_vs_reference_signed_orientation.tsv", sep="\t",
                index=False)
    print(f"[final] Part D: 01 {len(df01)}, 02 {len(df02)}, "
          f"03 {len(df03)}, 04 {len(df04)}, 05 {len(df05)}, 06 {len(df06)}")
    print(df04[df04.top_n.isin([200, 500])].to_string(index=False))

    # ---------- PART F: 10 secondary signed developmental orientation -------
    f_rows = []
    for program in PROGRAMS:
        for top_n in TOPS:
            s4t = s4_12[(s4_12["program"] == program) &
                        (s4_12["top_n"] == top_n)].iloc[0]
            s4c = s4_13[(s4_13["program"] == program) &
                        (s4_13["top_n"] == top_n)].iloc[0]
            ct = df03[(df03["program"] == program) &
                      (df03["top_n"] == top_n)].iloc[0]
            crec = df_c[(df_c["program"] == program) &
                        (df_c["top_n"] == top_n)].iloc[0]
            f_rows.append(dict(
                program=program, top_n=top_n,
                signed_top_stage=str(s4t["top_stage"]),
                signed_prenatal_minus_postnatal=float(
                    s4c["prenatal_minus_postnatal"]),
                signed_window=str(s4c["prenatal_or_postnatal"]),
                magnitude_top_stage=str(ct["top_stage"]),
                magnitude_localization_contrast=float(
                    ct["localization_contrast"]),
                magnitude_window=str(ct["prenatal_or_postnatal"]),
                identity_class=str(crec["projection_identity_class"]),
                identity_note=(
                    "SIGNED_DEVELOPMENTAL_ORIENTATION_SECONDARY: signed "
                    "readout is retained only as a secondary descriptor and "
                    "must not be conflated with magnitude developmental "
                    "localization")))
    df10 = pd.DataFrame(f_rows)
    df10.to_csv(RES / "10_secondary_signed_developmental_orientation.tsv",
                sep="\t", index=False)
    print(f"[final] Part F: 10 rows {len(df10)}")
    print(df10[df10.top_n.isin([200, 500])].to_string(index=False))
    print("[final] done")


if __name__ == "__main__":
    main()
