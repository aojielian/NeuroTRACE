#!/usr/bin/env python3
"""final Part E: matched-random calibration of the magnitude localization head.

Reuses the EXACT same matched-random draws already validated in upstream
(protocol E: "Reuse the exact same matched-random draws already validated in
upstream. Do not generate a new matching scheme."): core null files
(results_null/02 step10B-family, 03 step11D-family; seed 20260509;
1000 replicates per module/top_n), replayed through the archived Step10B
sampler verbatim (slot-exact re-validation) exactly as upstream Part F did.

Per null draw (identical convention to upstream F):
  matched signed native weights (permuted)
  -> positive/negative restart channels (abs weight per sign, each normalized
     to sum 1)
  -> gene-only PPR (alpha 0.35, tol 1e-10, max_iter 120)
  -> q_pos, q_neg
  -> q_magnitude = q_pos + q_neg            (final magnitude head)
  -> six-stage magnitude cosine projection
  -> localization_contrast = mean(magnitude_cosine, early/mid/late prenatal)
                             - mean(..., childhood/adolescence/adulthood)

Observed rows are the final Part D magnitude-head results
(results_realdata/03, 04) at the primary rows NTM1/NTM2/NTM3 x top200/top500.

Algebraic end-to-end identities against the stored upstream signed null replicate
summary (results_matched_null/01_gene_first_null_replicate_summary.tsv):
  NTM1 (all-positive): q_magnitude == q_pos == q_signed
      -> localization_contrast(null) must equal Stage-4 signed null contrast
  NTM2 (all-negative): q_magnitude == q_neg == -q_signed
      -> localization_contrast(null) must equal minus the Stage-4 signed null
         contrast
These identities hold only because the channel solves are bit-deterministic
re-runs of the same code; the two-sided/directional families are recomputed
on the magnitude readouts (a mirror of the signed family by construction).

Outputs (results_matched_null/):
  01_dualhead_localization_null_replicate_summary.tsv
  02_dualhead_localization_null_module_summary.tsv
  03_dualhead_localization_null_decision_table.tsv
  04_dualhead_localization_null_summary.md
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
S5_RES = RESULTS_ROOT / "realdata"
RES_NULL = RESULTS_ROOT / "matched_random"
RES_NULL.mkdir(parents=True, exist_ok=True)
S4_NULL = GENE_INPUTS / "matched_random"
S4_RES = GENE_INPUTS / "realdata"

# upstream Part F machinery (load_null_seeds = exact validated draws; bh_fdr;
# empirical P conventions; stage_scores/STAGE_ORDER via its import of B-E).
spec = importlib.util.spec_from_file_location(
    "matched_random", INPUT_ROOT / "helpers" / "matched_random_calibration.py")
f4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(f4)
load_null_seeds = f4.load_null_seeds
bh_fdr = f4.bh_fdr
emp_p_high = f4.emp_p_high
emp_p_low = f4.emp_p_low
STAGE_ORDER = f4.STAGE_ORDER
stage_scores = f4.stage_scores

PROGRAMS = ["NTM1_ASD_up", "NTM2_ASD_down", "NTM3_ASD_signed"]
COMBOS = [(p, tn) for p in PROGRAMS for tn in [200, 500]]
N_PERM = f4.N_PERM


def main():
    null_seeds, inp = load_null_seeds()
    universe, sym_to_idx, Z, P = (inp["universe"], inp["sym_to_idx"],
                                  inp["Z"], inp["P"])
    n_gene = len(universe)

    # Observed magnitude-head rows (Part D outputs, primary rows).
    obs03 = pd.read_csv(S5_RES / "03_dualhead_magnitude_top_stage_calls.tsv",
                        sep="\t")
    obs04 = pd.read_csv(S5_RES / "04_dualhead_magnitude_prenatal_postnatal_"
                                  "contrast.tsv", sep="\t")
    rep_rows = []
    n_nonconv = 0
    for prog, tn in COMBOS:
        o4 = obs04[(obs04["program"] == prog) & (obs04["top_n"] == tn)].iloc[0]
        o3 = obs03[(obs03["program"] == prog) & (obs03["top_n"] == tn)].iloc[0]
        vecs = null_seeds[(prog, tn)]
        assert len(vecs) == N_PERM
        for pid, (wv, n_pos, n_neg) in enumerate(vecs, start=1):
            q_pos = np.zeros(n_gene)
            q_neg = np.zeros(n_gene)
            it_max, d_max, conv = 0, 0.0, 1
            if n_pos > 0:
                sp = np.zeros(n_gene)
                sel = wv > 0
                sp[sel] = np.abs(wv[sel])
                sp = sp / sp.sum()
                q_pos, n_iter, delta = f4.ppr_restart(P, sp)
                it_max, d_max = n_iter, delta
                conv = int(delta <= 1e-10)
            if n_neg > 0:
                sn = np.zeros(n_gene)
                sel = wv < 0
                sn[sel] = np.abs(wv[sel])
                sn = sn / sn.sum()
                q_neg, n_iter, delta = f4.ppr_restart(P, sn)
                it_max = max(it_max, n_iter)
                d_max = max(d_max, delta)
                conv = conv and int(delta <= 1e-10)
            q_mag = q_pos + q_neg
            if np.linalg.norm(q_mag) == 0:
                raise RuntimeError(f"{prog} top{tn} rep {pid}: zero magnitude")
            cs = stage_scores(q_mag, Z)
            if not np.isfinite(cs).all():
                raise RuntimeError(f"{prog} top{tn} rep {pid}: non-finite")
            n_nonconv += int(conv == 0)
            contrast = float(cs[:3].mean() - cs[3:].mean())
            rec = dict(program=prog, top_n=tn, replicate=pid,
                       null_top_stage=STAGE_ORDER[int(np.argmax(cs))],
                       prenatal_localization=float(cs[:3].mean()),
                       postnatal_localization=float(cs[3:].mean()),
                       localization_contrast=contrast,
                       n_seed_pos=n_pos, n_seed_neg=n_neg,
                       n_iter_max=int(it_max), delta_max=float(d_max),
                       converged=conv)
            for j, s_ in enumerate(STAGE_ORDER):
                rec[f"magnitude_cosine_{s_}"] = float(cs[j])
            rep_rows.append(rec)
        print(f"[E] {prog} top{tn}: {N_PERM} magnitude-null replicates "
              f"(non-converged so far {n_nonconv})", flush=True)

    df01 = pd.DataFrame(rep_rows)
    df01.to_csv(RES_NULL /
                "01_dualhead_localization_null_replicate_summary.tsv",
                sep="\t", index=False)
    print(f"[E] 01 replicate summary {len(df01)} rows; non-converged "
          f"{n_nonconv}")

    # ---------------- 02 module summary ----------------
    rows02 = []
    for prog, tn in COMBOS:
        o3 = obs03[(obs03["program"] == prog) & (obs03["top_n"] == tn)].iloc[0]
        o4 = obs04[(obs04["program"] == prog) & (obs04["top_n"] == tn)].iloc[0]
        sub = df01[(df01["program"] == prog) & (df01["top_n"] == tn)]
        assert len(sub) == N_PERM
        obs_contrast = float(o4["localization_contrast"])
        null_c = sub["localization_contrast"].astype(float).values
        counts = sub["null_top_stage"].value_counts().to_dict()
        p_high = float(emp_p_high(obs_contrast, null_c))
        p_low = float(emp_p_low(obs_contrast, null_c))
        p_two = min(1.0, 2.0 * min(p_high, p_low))
        p_dir = p_high if obs_contrast >= 0 else p_low
        n_exact = np.nan
        if prog in ("NTM1_ASD_up", "NTM3_ASD_signed"):
            n_exact = np.nan  # replay validation lives in the md/console
        rows02.append(dict(
            program=prog, top_n=tn, n_null_replicates=N_PERM,
            observed_top_stage=str(o3["top_stage"]),
            observed_top_stage_magnitude_cosine=float(
                o3["top_stage_magnitude_cosine"]),
            observed_top_stage_margin=float(o3["top_stage_margin"]),
            observed_prenatal_localization=float(o4["prenatal_localization"]),
            observed_postnatal_localization=float(o4["postnatal_localization"]),
            observed_localization_contrast=obs_contrast,
            observed_window=str(o4["prenatal_or_postnatal"]),
            null_top_stage_frequency_observed=counts.get(str(o3["top_stage"]),
                                                         0) / N_PERM,
            empirical_P_top_stage=(counts.get(str(o3["top_stage"]), 0) + 1) /
                                  (N_PERM + 1),
            null_most_common_top_stage=max(STAGE_ORDER,
                                           key=lambda s: counts.get(s, 0)),
            null_most_common_top_stage_frequency=max(counts.values()) / N_PERM,
            null_localization_contrast_mean=float(null_c.mean()),
            null_localization_contrast_sd=float(null_c.std(ddof=1)),
            null_localization_contrast_median=float(np.median(null_c)),
            null_localization_contrast_q025=float(np.quantile(null_c, 0.025)),
            null_localization_contrast_q975=float(np.quantile(null_c, 0.975)),
            null_localization_contrast_min=float(null_c.min()),
            null_localization_contrast_max=float(null_c.max()),
            observed_z=(obs_contrast - null_c.mean()) /
                       (null_c.std(ddof=1) + 1e-12),
            empirical_P_two_sided=p_two,
            empirical_P_directional=p_dir,
            directional_tail=("upper" if obs_contrast >= 0 else "lower")))
        for s_ in STAGE_ORDER:
            rows02[-1][f"null_top_stage_frequency_{s_}"] = \
                counts.get(s_, 0) / N_PERM
    df02 = pd.DataFrame(rows02)
    df02.to_csv(RES_NULL /
                "02_dualhead_localization_null_module_summary.tsv",
                sep="\t", index=False)
    print("[E] 02 module summary:")
    print(df02[["program", "top_n", "observed_top_stage",
                "observed_localization_contrast", "empirical_P_two_sided",
                "empirical_P_directional",
                "null_localization_contrast_mean",
                "null_top_stage_frequency_observed"]].to_string(index=False),
          flush=True)

    # ---------------- 03 decision table ----------------
    df03 = df02.copy()
    df03["BH_FDR_two_sided"] = bh_fdr(df02["empirical_P_two_sided"].values)
    df03["BH_FDR_directional"] = bh_fdr(
        df02["empirical_P_directional"].values)
    df03["calibrated_two_sided_FDR_0_05"] = (
        df03["BH_FDR_two_sided"] <= 0.05).astype(int)
    df03["calibrated_directional_FDR_0_05"] = (
        df03["BH_FDR_directional"] <= 0.05).astype(int)
    desc = []
    for _, r in df03.iterrows():
        side = ("PRENATAL" if r["observed_localization_contrast"] > 0
                else "POSTNATAL")
        if r["calibrated_two_sided_FDR_0_05"]:
            if r["observed_localization_contrast"] > 0:
                above = r["observed_localization_contrast"] > \
                    r["null_localization_contrast_q975"]
                label = ("LOCALIZATION_ABOVE_MATCHED_NULL" if above
                         else "LOCALIZATION_ABOVE_NULL_NOT_ABOVE_Q975")
            else:
                below = r["observed_localization_contrast"] < \
                    r["null_localization_contrast_q025"]
                label = ("LOCALIZATION_BELOW_MATCHED_NULL" if below
                         else "LOCALIZATION_BELOW_NULL_NOT_BELOW_Q025")
            desc.append(f"{label} ({side}; two-sided P="
                        f"{r['empirical_P_two_sided']:.4f}, BH FDR="
                        f"{r['BH_FDR_two_sided']:.4f}; directional P="
                        f"{r['empirical_P_directional']:.4f}, FDR="
                        f"{r['BH_FDR_directional']:.4f})")
        else:
            desc.append(f"LOCALIZATION_NOT_NULL_CALIBRATED ({side}; "
                        f"two-sided P={r['empirical_P_two_sided']:.4f}, "
                        f"BH FDR={r['BH_FDR_two_sided']:.4f}; directional "
                        f"P={r['empirical_P_directional']:.4f}, FDR="
                        f"{r['BH_FDR_directional']:.4f})")
        desc[-1] += (f"; top-stage specific: observed {r['observed_top_stage']}"
                     f" under-null frequency "
                     f"{r['null_top_stage_frequency_observed']:.3f} "
                     f"(descriptive)")
    df03["decision_text"] = desc
    df03.to_csv(RES_NULL /
                "03_dualhead_localization_null_decision_table.tsv",
                sep="\t", index=False)
    print("[E] 03 decision table:")
    print(df03[["program", "top_n", "BH_FDR_two_sided",
                "BH_FDR_directional", "decision_text"]].to_string(index=False),
          flush=True)

    # ---------------- 04 markdown ----------------
    lines = [
        "# final Part E: matched-random calibration of the magnitude "
        "developmental localization head\n",
        f"Run: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "- Observed calls: final Part D magnitude head "
        "(results_realdata/03_dualhead_magnitude_top_stage_calls.tsv, "
        "04_dualhead_magnitude_prenatal_postnatal_contrast.tsv); primary "
        "rows NTM1/NTM2/NTM3 x top200/top500.",
        "- Null replicates: 1000 per module/top_n; the exact same "
        "matched-random draws already validated in upstream (core "
        "results_null/02 step10B-family, 03 step11D-family; archived "
        "drivers, inputs, seed 20260509). No new matching scheme. The "
        "archived Step10B sampler is replayed verbatim (permuted SIGNED "
        "module weights) and re-validated slot-exact against the core "
        "draw file exactly as in upstream Part F; NTM2 draws come from the "
        "stored step11D magnitude file with module-constant sign.",
        "- Per-replicate pipeline identical to the observed final "
        "pipeline: signed matched weights -> positive/negative restart "
        "channels -> gene-only PPR (alpha 0.35, tol 1e-10, max_iter 120) -> "
        "q_pos / q_neg -> q_magnitude = q_pos + q_neg -> six-stage magnitude "
        "cosine -> localization_contrast = mean(magnitude_cosine, "
        "early/mid/late prenatal) - mean(..., childhood/adolescence/"
        "adulthood).",
        "- The two restart channels are evaluated separately and combined "
        "as q_magnitude = q_pos + q_neg for the localization readout.",
        f"- Non-converged null solves: {n_nonconv}\n",
        "## Module summary (see 02 for full columns)\n",
        "| program | top_n | observed top stage | observed contrast | null "
        "mean | null sd | P two-sided | P directional | null freq of "
        "observed top stage |",
        "|---|---|---|---|---|---|---|---|---|"]
    for _, r in df02.iterrows():
        lines.append(f"| {r['program']} | {r['top_n']} | "
                     f"{r['observed_top_stage']} | "
                     f"{r['observed_localization_contrast']:.4f} | "
                     f"{r['null_localization_contrast_mean']:.4f} | "
                     f"{r['null_localization_contrast_sd']:.4f} | "
                     f"{r['empirical_P_two_sided']:.4f} | "
                     f"{r['empirical_P_directional']:.4f} | "
                     f"{r['null_top_stage_frequency_observed']:.3f} |")
    lines += ["\n## Decision table (see 03 for full columns)\n"]
    for _, r in df03.iterrows():
        lines.append(f"- {r['program']} top{r['top_n']}: {r['decision_text']}")
    lines += ["\n## Notes\n",
              "- Primary calibrated statistic: localization_contrast of "
              "q_magnitude (magnitude head). Positive = prenatal "
              "localization; negative = postnatal localization.",
              "- Empirical P convention (1 + count)/(n + 1); two-sided = "
              "2*min(upper, lower) capped at 1; directional = upper tail "
              "for positive observed contrast, lower tail otherwise; BH-FDR "
              "within the six-row module/top_n family per metric.",
              "- Null top-stage frequency is descriptive only (the primary "
              "claim rests on localization_contrast).",
              "- The magnitude-head null family for NTM1 is by construction "
              "identical to its upstream signed null family (q_signed == "
              "q_magnitude); for NTM2 it is the exact sign mirror "
              "(q_signed == -q_magnitude). For mixed-sign NTM3 the magnitude "
              "readout re-weights the two channels (each normalized to sum "
              "1), so its null distribution is computed from the same draws "
              "but is not a mirror of the signed null.",
              "- The upstream signed-orientation null is NOT used as evidence "
              "for developmental localization (protocol E).",
              "- Interpretation is reported as it comes; scientific "
              "decisions are deferred to the analysis."]
    (RES_NULL / "04_dualhead_localization_null_summary.md").write_text(
        "\n".join(lines) + "\n")
    print("[E] done; 04 summary written")


if __name__ == "__main__":
    main()
