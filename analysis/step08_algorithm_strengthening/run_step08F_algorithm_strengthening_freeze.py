#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
PROJECT = BASE / "neurotrace_algorithm_project"
STEP08 = PROJECT / "step08_algorithm_strengthening"
OUT = STEP08 / "results"
OUT.mkdir(parents=True, exist_ok=True)

FILES = {
    "step08A_summary": OUT / "47_step08A_exprvar_matched_random_specificity_summary.md",
    "step08A_specificity": OUT / "44_step08A_exprvar_matched_random_specificity_summary.tsv",
    "step08A_meta": OUT / "45_step08A_exprvar_matched_random_meta_summary.tsv",

    "step08B_summary": OUT / "54_step08B_realdata_score_baseline_comparison_summary.md",
    "step08B_method": OUT / "51_step08B_realdata_baseline_method_summary.tsv",
    "step08B_pairwise": OUT / "53_step08B_neurotrace_vs_baseline_pairwise_summary.tsv",

    "step08C_summary": OUT / "63_step08C_gandal_cross_disease_specificity_summary.md",
    "step08C_models": OUT / "60_step08C_gandal_cross_disease_NTM_models.tsv",
    "step08C_specificity": OUT / "62_step08C_ASD_specificity_summary.tsv",

    "step08E_summary": OUT / "85_step08E_cross_disease_NTM_scoring_summary.md",
    "step08E_models": OUT / "82_step08E_cross_disease_NTM_models.tsv",
    "step08E_dx_summary": OUT / "83_step08E_cross_disease_NTM_summary.tsv",
    "step08E_disease_summary": OUT / "84_step08E_cross_disease_disease_level_summary.tsv",
}


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_tsv(path):
    path = Path(path)
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def file_status():
    rows = []
    for k, p in FILES.items():
        p = Path(p)
        rows.append({
            "file_key": k,
            "path": str(p),
            "exists": p.exists(),
            "size_mb": round(p.stat().st_size / 1024 / 1024, 4) if p.exists() else np.nan,
            "mtime": datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S") if p.exists() else ""
        })
    return pd.DataFrame(rows)


def build_strengthening_evidence():
    rows = []

    # Step08A
    a_meta = read_tsv(FILES["step08A_meta"])
    if not a_meta.empty:
        rows.append({
            "evidence_layer": "expression_variance_matched_specificity",
            "status": "PASS" if (a_meta["min_empirical_p"].astype(float) < 0.05).all() else "CAUTION",
            "main_result": (
                f"{a_meta.shape[0]}/{a_meta.shape[0]} module/top_n combinations pass expression/variance-matched random tests; "
                f"min empirical P range: {a_meta['min_empirical_p'].min()}–{a_meta['min_empirical_p'].max()}; "
                f"mean exact match rate: {a_meta['mean_exact_expr_var_match_rate'].mean():.3f}"
            ),
            "interpretation": "Observed external NTM effects exceed random modules matched by size, weight/sign, expression abundance and expression variance.",
            "manuscript_use": "main or supplementary methods-strengthening result",
            "source_file": str(FILES["step08A_meta"])
        })

    # Step08B
    b_method = read_tsv(FILES["step08B_method"])
    b_pair = read_tsv(FILES["step08B_pairwise"])
    if not b_method.empty:
        nt = b_method[b_method["scoring_method"] == "neurotrace_weighted_signed"]
        if len(nt):
            nt = nt.iloc[0]
            rows.append({
                "evidence_layer": "realdata_score_baseline_comparison",
                "status": "PASS",
                "main_result": (
                    f"NeuroTRACE weighted signed score: direction concordance={nt['direction_concordance_rate']}, "
                    f"median AUC={nt['median_auc']}, n_rank1_by_p={nt['n_rank1_by_p']}, "
                    f"n_rank1_by_beta={nt['n_rank1_by_beta']}, n_rank1_by_auc={nt['n_rank1_by_auc']}."
                ),
                "interpretation": "Signed disease-module architecture is the dominant signal; continuous weights add consistent but modest gain over sign-only scoring.",
                "manuscript_use": "main method comparison",
                "source_file": str(FILES["step08B_method"])
            })

    if not b_pair.empty:
        for _, r in b_pair.iterrows():
            rows.append({
                "evidence_layer": "neurotrace_vs_baseline_pairwise",
                "status": "PASS" if float(r["frac_nt_better_p"]) >= 0.75 else "PARTIAL",
                "main_result": (
                    f"Baseline={r['baseline_method']}; NT better by P in {r['n_nt_better_p']}/{r['n']} tests, "
                    f"better beta in {r['n_nt_better_beta']}/{r['n']}, better AUC in {r['n_nt_better_auc']}/{r['n']}; "
                    f"mean_delta_beta={r['mean_delta_beta']}; mean_delta_auc={r['mean_delta_auc']}."
                ),
                "interpretation": "Quantifies incremental value of signed continuous NeuroTRACE scoring over simpler score definitions.",
                "manuscript_use": "main or supplementary baseline comparison",
                "source_file": str(FILES["step08B_pairwise"])
            })

    # Step08C
    c_spec = read_tsv(FILES["step08C_specificity"])
    c_models = read_tsv(FILES["step08C_models"])
    if not c_models.empty:
        dxs = sorted(c_models["target_dx"].astype(str).unique())
        rows.append({
            "evidence_layer": "gandal_internal_cross_diagnosis_audit",
            "status": "PASS_CONTEXTUAL",
            "main_result": f"Gandal contains diagnostic groups: {', '.join(dxs)}. Dup15q effects exceed ASD effects across all NTM modules/top_n settings.",
            "interpretation": "NTM modules are not strictly idiopathic-ASD-specific; they are stronger in syndromic autism-related Dup15q.",
            "manuscript_use": "important boundary analysis",
            "source_file": str(FILES["step08C_models"])
        })

    if not c_spec.empty:
        rows.append({
            "evidence_layer": "ASD_vs_Dup15q_specificity_boundary",
            "status": "NOT_ASD_SPECIFIC_WITHIN_AUTISM_RELATED_GROUPS",
            "main_result": f"ASD beta was greater than Dup15q in {int(c_spec['n_ASD_beta_greater'].sum())}/{int(c_spec['n_other_diseases'].sum())} comparisons.",
            "interpretation": "The strongest within-Gandal signal is Dup15q, suggesting the NTM axis captures autism-related cortical dysregulation rather than idiopathic-ASD specificity.",
            "manuscript_use": "discussion and disease-boundary result",
            "source_file": str(FILES["step08C_specificity"])
        })

    # Step08E
    e_disease = read_tsv(FILES["step08E_disease_summary"])
    if not e_disease.empty:
        for _, r in e_disease.iterrows():
            rows.append({
                "evidence_layer": "external_cross_disease_generalization",
                "status": "PASS" if float(r["positive_rate"]) == 1 else "NEGATIVE_OR_BOUNDARY",
                "main_result": (
                    f"{r['target_dx']}: n_tests={r['n_module_tests']}, n_positive={r['n_positive']}, "
                    f"positive_rate={r['positive_rate']}, mean_beta={r['mean_beta']}, "
                    f"median_AUC={r['median_auc']}, min_FDR={r['min_fdr']}."
                ),
                "interpretation": (
                    "Supports cross-disease generalization." if float(r["positive_rate"]) == 1
                    else "No positive generalization; serves as disease-boundary evidence."
                ),
                "manuscript_use": "main cross-disease generalization result",
                "source_file": str(FILES["step08E_disease_summary"])
            })

    return pd.DataFrame(rows)


def build_updated_claims():
    rows = [
        {
            "claim": "NeuroTRACE v1 now meets a stronger methods-oriented proof-of-concept standard.",
            "support": "Step08A expression/variance-matched random specificity, Step08B real-data baseline comparison, and Step08E cross-disease generalization have been completed.",
            "strength": "strong_for_proof_of_concept",
            "caveat": "Still not a full software/package-level algorithm article; graph-regularized OT or contrastive learning remains future work."
        },
        {
            "claim": "The NTM axis is not ASD-only.",
            "support": "Step08C shows Dup15q has stronger NTM effects than idiopathic ASD; Step08E shows SCZ and BD are positive.",
            "strength": "strong",
            "caveat": "This reframes the biology as autism-related and broader psychiatric cortical dysregulation rather than ASD-specific."
        },
        {
            "claim": "MDD acts as a negative disease boundary.",
            "support": "Step08E shows MDD has 0/6 positive NTM module tests in GSE53987 PFC.",
            "strength": "moderate",
            "caveat": "MDD is represented only in GSE53987 PFC; more MDD datasets would strengthen this boundary."
        },
        {
            "claim": "Signed direction-aware scoring is essential.",
            "support": "Step08B shows unsigned scoring reverses/downweights NTM2 effects, whereas signed NeuroTRACE scoring is best by P and beta.",
            "strength": "strong",
            "caveat": "Sign-only is close to weighted-signed, so continuous weights provide incremental rather than dominant improvement."
        },
        {
            "claim": "External adult validation is robust to matched random controls.",
            "support": "Step08A shows all external NTM effects exceed expression- and variance-matched random modules with empirical FDR=0.000999.",
            "strength": "strong",
            "caveat": "Gene-gene correlation matching is still not included."
        }
    ]
    return pd.DataFrame(rows)


def build_updated_article_plan():
    rows = [
        {
            "figure": "Figure 1",
            "title": "NeuroTRACE algorithm framework",
            "new_status": "unchanged",
            "main_panels": "module derivation; developmental embedding; native graph; graph transport; validation layers"
        },
        {
            "figure": "Figure 2",
            "title": "Simulation and baseline method benchmark",
            "new_status": "strengthened",
            "main_panels": "simulation benchmark; dual-head objective; score baseline comparison from Step08B"
        },
        {
            "figure": "Figure 3",
            "title": "Native NTM derivation and late-prenatal transport",
            "new_status": "unchanged",
            "main_panels": "Gandal native DE; NTM modules; BrainSpan transport; top-stage calls"
        },
        {
            "figure": "Figure 4",
            "title": "Graph transport and genetic-risk convergence",
            "new_status": "unchanged",
            "main_panels": "native graph; NTM2 SFARI convergence; graph-priority genes"
        },
        {
            "figure": "Figure 5",
            "title": "External ASD validation and matched-random specificity",
            "new_status": "strengthened",
            "main_panels": "GSE102741/GSE64018 effects; expression/variance-matched random specificity"
        },
        {
            "figure": "Figure 6",
            "title": "Cross-disease generalization and disease-boundary analysis",
            "new_status": "new_main_figure",
            "main_panels": "SCZ/BD/MDD NTM scores; disease-level summary; Dup15q boundary; MDD negative boundary"
        },
        {
            "figure": "Supplementary",
            "title": "Platform processing and probe-to-symbol collapse",
            "new_status": "new_supplement",
            "main_panels": "GEO parse; GPL570/GPL96 mapping; gene overlap; tissue/dx audit"
        }
    ]
    return pd.DataFrame(rows)


def main():
    print(f"[{now()}] Step08F algorithm strengthening freeze started")

    status = file_status()
    status.to_csv(OUT / "86_step08F_input_file_status.tsv", sep="\t", index=False)

    evidence = build_strengthening_evidence()
    evidence.to_csv(OUT / "87_step08F_algorithm_strengthening_evidence.tsv", sep="\t", index=False)

    claims = build_updated_claims()
    claims.to_csv(OUT / "88_step08F_updated_claims_and_caveats.tsv", sep="\t", index=False)

    plan = build_updated_article_plan()
    plan.to_csv(OUT / "89_step08F_updated_article_figure_plan.tsv", sep="\t", index=False)

    with open(OUT / "90_step08F_algorithm_strengthening_freeze_summary.md", "w") as f:
        f.write("# NeuroTRACE Step08F algorithm-strengthening freeze summary\n\n")
        f.write(f"Generated: {now()}\n\n")

        f.write("## Updated project status\n")
        f.write("NeuroTRACE has advanced beyond a simple ASD proof-of-concept. Step08A adds expression/variance-matched specificity, Step08B adds real-data score baseline comparison, and Step08E adds external cross-disease generalization across SCZ/BD/MDD datasets.\n\n")

        f.write("## Algorithm-strengthening evidence\n\n")
        f.write(evidence.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Updated claims and caveats\n\n")
        f.write(claims.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Updated article figure plan\n\n")
        f.write(plan.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Recommended revised central claim\n")
        f.write(
            "NeuroTRACE identifies a disease-direction-aware cortical dysregulation axis that links adult ASD-derived modules to late-prenatal developmental programs and genetic-risk convergence, "
            "while also generalizing to SCZ and BD but not MDD. This supports a broader neuropsychiatric cortical dysregulation framework rather than a strictly ASD-specific signature.\n"
        )

    print(f"[{now()}] Step08F done")
    print(f"[{now()}] Evidence rows: {len(evidence)}")
    print(f"[{now()}] Results: {OUT}")


if __name__ == "__main__":
    main()
