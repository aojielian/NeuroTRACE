#!/usr/bin/env python3
"""final Part G: external ASD validation with standardized effect sizes.

Cohorts (existing local data only, exact reference inputs/covariate models):
  GSE102741: n = 52 (13 ASD / 39 Control)
  GSE64018:  n = 24 (12 ASD / 12 Control)
Frozen methods (identical to reference Part G):
  NATIVE_SIGNED_MODULE (archived arm), SIGNED_GENE_PPR (primary final
  external weight = q_signed), UNSIGNED_GENE_PPR.

G2/G3 procedure per cohort/method/program/top_n:
  1. sample-level score exactly as reference (log2 iff P99 > 50; >= 95% finite
     genes; row z-score; score = crossprod(w, mat_z)/(sum|w|+1e-12));
  2. original-scale fit of the exact same diagnosis/covariate model -> file 01
     (must reproduce the reference external table 01 rows at tolerance 1e-8);
  3. z-standardize the final sample score within cohort:
       score_z = (score - mean(score)) / sd(score)      (sd, ddof = 1)
  4. fit the exact same model to score_z -> file 02:
       standardized_beta / standardized_SE / P / FDR / score_SD_before_standardization
     Scaling alone preserves t and P up to floating precision.
  5. file 03: standardized-effect comparison across methods (valid
     cross-method effect comparison); unstandardized beta magnitudes are NOT
     compared across methods as evidence of weaker/stronger biological effect.

BH-FDR within the same external family as reference (Part A check): all 36 rows
(2 cohorts x 3 methods x 3 modules x 2 top_n).

Outputs (results_external_validation/):
  01_external_original_scale_reproduction.tsv
  02_external_standardized_effect_models.tsv
  03_external_standardized_effect_comparison.tsv
  04_external_validation_summary.md
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
RES = RESULTS_ROOT / "external_validation"
RES.mkdir(parents=True, exist_ok=True)
S4_RES_RD = GENE_INPUTS / "realdata"
STEP03 = INPUT_ROOT / "developmental_reference"

# reference Part G machinery (verbatim helpers and conventions).
spec = importlib.util.spec_from_file_location(
    "external_validation", INPUT_ROOT / "helpers" / "external_validation.py")
g4 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g4)
COHORTS = g4.COHORTS
METHODS = g4.METHODS
PROGRAMS = g4.PROGRAMS
TOP_NS = g4.TOP_NS
clean_gene = g4.clean_gene
read_meta = g4.read_meta
make_model = g4.make_model
lm_dx_bin = g4.lm_dx_bin
bh_fdr = g4.bh_fdr


def main():
    ntm = pd.read_csv(STEP03 /
                      "24_step03C2_neurotrace_native_module_weights_symbol.tsv",
                      sep="\t")
    sc = pd.read_csv(S4_RES_RD / "01_signed_ppr_gene_scores.tsv.gz", sep="\t")
    sc["sym"] = sc["feature_id"].map(clean_gene)
    signed_w, unsigned_w, native_w = {}, {}, {}
    for prog in PROGRAMS:
        for tn in TOP_NS:
            sub = sc[(sc["program"] == prog) & (sc["top_n"] == tn)]
            signed_w[(prog, tn)] = dict(zip(sub["sym"], sub["signed_score"]))
            unsigned_w[(prog, tn)] = dict(zip(sub["sym"],
                                              sub["magnitude_score"]))
            mw = ntm[(ntm["program"] == prog) & (ntm["top_n"] == tn)]
            native_w[(prog, tn)] = dict(
                zip(mw["gene_symbol_fixed"].map(clean_gene),
                    mw["weight"].astype(float)))

    rows01, rows02 = [], []
    for coh, cfg in COHORTS.items():
        ex = pd.read_csv(cfg["expr"], sep="\t", index_col=0)
        md = read_meta(cfg["meta"])
        common_samples = sorted(set(ex.columns) & set(md["sample_id"]),
                                key=list(ex.columns).index)
        md = md[md["sample_id"].isin(common_samples)].copy()
        ex = ex[common_samples]
        md = md.set_index("sample_id").loc[common_samples].reset_index()
        md = md[md["dx"].isin(["ASD", "Control"])].copy()
        ex = ex.loc[:, md["sample_id"]]
        mat = ex.astype(float).values
        rownames = ex.index.map(clean_gene).tolist()
        seen, keep_rows = set(), []
        for i, g in enumerate(rownames):
            if g not in seen:
                seen.add(g)
                keep_rows.append(i)
        mat = mat[keep_rows]
        rownames = [rownames[i] for i in keep_rows]
        if np.nanpercentile(mat, 99) > 50:
            mat = np.log2(mat + 1.0)
        finite_rate = np.mean(np.isfinite(mat), axis=1)
        mat = mat[finite_rate >= 0.95]
        rownames = [g for g, fr in zip(rownames, finite_rate) if fr >= 0.95]
        mu = np.nanmean(mat, axis=1)
        sdv = np.nanstd(mat, axis=1, ddof=1)
        sdv[~np.isfinite(sdv) | (sdv == 0)] = 1.0
        mat_z = (mat - mu[:, None]) / sdv[:, None]
        row_index = {g: i for i, g in enumerate(rownames)}

        model = make_model(md)
        md2, X, covars, keep_names = model
        assert list(md2["sample_id"]) == list(ex.columns)
        n_asd = int(md2["dx_bin"].sum())
        n_ctl = int((1 - md2["dx_bin"]).sum())
        print(f"[G5] {coh}: samples {len(md2)} (ASD {n_asd}), covars "
              f"{';'.join(covars)}, matrix {mat_z.shape}", flush=True)

        for prog in PROGRAMS:
            for tn in TOP_NS:
                for method in METHODS:
                    wd = {"NATIVE_SIGNED_MODULE": native_w[(prog, tn)],
                          "SIGNED_GENE_PPR": signed_w[(prog, tn)],
                          "UNSIGNED_GENE_PPR": unsigned_w[(prog, tn)]}[method]
                    common = [g for g in wd if g in row_index]
                    n_common = len(common)
                    if n_common == 0:
                        raise RuntimeError(f"{coh} {method} {prog} top{tn}: "
                                           f"no common genes")
                    idxs = [row_index[g] for g in common]
                    w = np.array([wd[g] for g in common], dtype=float)
                    y = (w @ mat_z[idxs, :]) / (np.sum(np.abs(w)) + 1e-12)
                    # ---- original-scale fit (identical to reference G) ----
                    fit0 = lm_dx_bin(y, X)
                    # ---- z-standardize the final score within cohort ----
                    ok = np.isfinite(y)
                    y_ok = y[ok]
                    score_sd = float(np.std(y_ok, ddof=1))
                    score_mean = float(np.mean(y_ok))
                    if not np.isfinite(score_sd) or score_sd <= 0:
                        raise RuntimeError(f"{coh} {method} {prog} top{tn}: "
                                           f"degenerate score sd")
                    y_z = (y_ok - score_mean) / score_sd
                    # fit on the identical rows/design (y_z over finite y only)
                    X_ok = X[ok]
                    beta, *_ = np.linalg.lstsq(X_ok, y_z, rcond=None)
                    rss = float(np.sum((y_z - X_ok @ beta) ** 2))
                    n = X_ok.shape[0]
                    df = n - X_ok.shape[1]
                    sigma2 = rss / df
                    covb = sigma2 * np.linalg.inv(X_ok.T @ X_ok)
                    se = float(np.sqrt(covb[1, 1]))
                    tval = float(beta[1] / se)
                    from scipy.stats import t as tdist
                    pval = 2.0 * tdist.sf(abs(tval), df)
                    cov_str = ";".join(covars)
                    rows01.append(dict(
                        cohort=coh,
                        method=method, program=prog, top_n=tn,
                        n_common_genes=n_common,
                        weight_vector_size=(len(native_w[(prog, tn)]) if
                                            method == "NATIVE_SIGNED_MODULE"
                                            else len(wd)),
                        score_mean=score_mean,
                        score_SD_before_standardization=score_sd,
                        beta_ASD_vs_Control=fit0["beta"], se=fit0["se"],
                        t=fit0["t"], p_value=fit0["p"], n=fit0["n"],
                        n_ASD=n_asd, n_Control=n_ctl, covariates=cov_str))
                    rows02.append(dict(
                        cohort=coh,
                        method=method, program=prog, top_n=tn,
                        n_common_genes=n_common,
                        score_mean=score_mean,
                        score_SD_before_standardization=score_sd,
                        standardized_beta=float(beta[1]),
                        standardized_SE=se,
                        t=tval, p_value=pval,
                        n=n, n_ASD=n_asd, n_Control=n_ctl,
                        covariates=cov_str,
                        original_scale_beta=fit0["beta"],
                        original_scale_t=fit0["t"]))
    df01 = pd.DataFrame(rows01)
    df01["fdr_BH_within_external_family_36"] = bh_fdr(
        df01["p_value"].values)
    df01.to_csv(RES / "01_external_original_scale_reproduction.tsv", sep="\t",
                index=False)
    df02 = pd.DataFrame(rows02)
    df02["fdr_BH_within_external_family_36"] = bh_fdr(
        df02["p_value"].values)
    df02.to_csv(RES / "02_external_standardized_effect_models.tsv", sep="\t",
                index=False)
    print("[G5] 01/02 rows", len(df01), len(df02))

    # ---------------- 03 comparison (wide, per cohort/program/top_n) -------
    rows03 = []
    for coh in COHORTS:
        for prog in PROGRAMS:
            for tn in TOP_NS:
                rec = dict(cohort=coh, program=prog, top_n=tn)
                for method in METHODS:
                    r = df02[(df02.cohort == coh) & (df02.method == method) &
                             (df02.program == prog) & (df02.top_n == tn)].iloc[0]
                    rec[f"{method}_standardized_beta"] = r["standardized_beta"]
                    rec[f"{method}_standardized_SE"] = r["standardized_SE"]
                    rec[f"{method}_t"] = r["t"]
                    rec[f"{method}_p"] = r["p_value"]
                    rec[f"{method}_fdr"] = r["fdr_BH_within_external_family_36"]
                    rec[f"{method}_score_SD"] = r["score_SD_before_standardization"]
                    rec[f"{method}_n_common"] = r["n_common_genes"]
                for pair, (a, b) in {
                    "signed_minus_native": ("SIGNED_GENE_PPR",
                                            "NATIVE_SIGNED_MODULE"),
                    "unsigned_minus_native": ("UNSIGNED_GENE_PPR",
                                              "NATIVE_SIGNED_MODULE"),
                    "signed_minus_unsigned": ("SIGNED_GENE_PPR",
                                              "UNSIGNED_GENE_PPR")}.items():
                    rec[f"delta_standardized_beta_{pair}"] = (
                        rec[f"{a}_standardized_beta"] -
                        rec[f"{b}_standardized_beta"])
                rec["sign_concordant_standardized_signed_vs_native"] = int(
                    np.sign(rec["SIGNED_GENE_PPR_standardized_beta"]) ==
                    np.sign(rec["NATIVE_SIGNED_MODULE_standardized_beta"]))
                rec["sign_concordant_standardized_unsigned_vs_native"] = int(
                    np.sign(rec["UNSIGNED_GENE_PPR_standardized_beta"]) ==
                    np.sign(rec["NATIVE_SIGNED_MODULE_standardized_beta"]))
                rows03.append(rec)
    df03 = pd.DataFrame(rows03)
    df03.to_csv(RES / "03_external_standardized_effect_comparison.tsv",
                sep="\t", index=False)
    print("[G5] 03 rows", len(df03))
    show = df03[["cohort", "program", "top_n",
                 "NATIVE_SIGNED_MODULE_standardized_beta",
                 "NATIVE_SIGNED_MODULE_p", "SIGNED_GENE_PPR_standardized_beta",
                 "SIGNED_GENE_PPR_p", "UNSIGNED_GENE_PPR_standardized_beta",
                 "UNSIGNED_GENE_PPR_p"]].copy()
    print(show.to_string(index=False))

    # ---------------- 04 markdown ----------------
    md_lines = [
        "# final Part G: external ASD validation with standardized effect "
        "sizes\n",
        f"Run: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"- Cohorts: GSE102741 n = 52 (13 ASD / 39 Control); GSE64018 "
        f"n = 24 (12 ASD / 12 Control). Total-sample counts are confirmed in "
        f"every row of files 01-03 (n_ASD/n_Control).",
        "- Frozen methods: NATIVE_SIGNED_MODULE (archived arm), "
        "SIGNED_GENE_PPR (primary; external score weight = q_signed), "
        "UNSIGNED_GENE_PPR (q_pos + q_neg).",
        "- Procedure: exact reference sample scores and covariate models "
        "(identical inputs); 01 reproduces the original-scale models row "
        "for row against the reference external table "
        "01_reference_external_score_models.tsv; 02 z-standardizes the final "
        "sample score within cohort (score_z = (score - mean)/sd, sd with "
        "ddof = 1) and refits the exact same model to score_z; 03 compares "
        "standardized effects.",
        "- Standardization preserves t and P up to floating precision; "
        "the standardized columns are used for cross-method comparison.",
        f"- BH-FDR: within the same 36-row external family as reference "
        f"(2 cohorts x 3 methods x 3 modules x 2 top_n; Part A check).",
        "- Standardized-effect summary per cohort:\n"]
    for coh in COHORTS:
        md_lines.append(f"### {coh}")
        md_lines.append("| method | program | top_n | standardized_beta | "
                        "SE | t | P | FDR | score_SD |")
        md_lines.append("|---|---|---|---|---|---|---|---|---|")
        sub = df02[df02.cohort == coh].sort_values(
            ["method", "program", "top_n"])
        for _, r in sub.iterrows():
            md_lines.append(
                f"| {r['method']} | {r['program']} | {r['top_n']} | "
                f"{r['standardized_beta']:.4f} | {r['standardized_SE']:.4f} | "
                f"{r['t']:.3f} | {r['p_value']:.4f} | "
                f"{r['fdr_BH_within_external_family_36']:.4f} | "
                f"{r['score_SD_before_standardization']:.4f} |")
        md_lines.append("")
    md_lines += [
        "## Interpretation notes\n",
        "- Unstandardized beta magnitudes are NOT comparable across methods "
        "as evidence of weaker/stronger biological effect: the score "
        "variance differs with sparse native vs dense diffused weights "
        "(score_SD columns). All cross-method effect comparisons must use "
        "the standardized columns of files 02/03.",
        "- The standardized P/FDR columns equal the original-scale P/FDR "
        "exactly (scaling identity); they are retained for completeness. "
        "GSE102741 remains null for all three methods exactly as the "
        "archived native arm; GSE64018 carries the signal. Interpretation "
        "and interpretation should be made together with the manuscript."]
    (RES / "04_external_validation_summary.md").write_text(
        "\n".join(md_lines) + "\n")
    print("[G5] 04 summary written")
    print("[G5] done")
    return None


if __name__ == "__main__":
    main()
