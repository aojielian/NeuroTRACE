#!/usr/bin/env python3
"""final Part J: required confirmatory interpretation table
results_confirmatory_simulation/09_dualhead_pareto_summary.tsv

For each comparator method the two major domains are reported separately
(protocol Part J):

    developmental_alignment:   alignment_hit_top1, true_stage_mass
    gene_prioritization:       gene_AUPRC, precision_at_50

and whether DUAL_HEAD_GENE_FIRST_PPR is better / tied / worse than the
comparator on each domain. No global-comparison claim is made.

Verdict rule (declared here, deterministic, per (benchmark, comparator,
domain), all read from 05_confirmatory_paired_deltas.tsv, condition-level
rows only):

1. Cell = one (benchmark, condition, comparator, metric) paired comparison,
   metric direction as in 05 (all four Part J metrics are higher-is-better).
   Cells exist only for metric-comparator pairs the comparator defines
   (reference availability conventions recorded in 05; e.g. true_stage_mass is
   defined only for the affinity-bearing family comparators, and 05 contains
   no cell otherwise).
2. nominal_cell = BETTER/WORSE if |mean_delta_dualhead_minus_baseline| >
   1e-12 (sign with respect to metric direction), else TIED.
3. significant = fdr_BH_within_family < 0.05 (BH families are
   benchmark x metric over per-condition rows, declared in 05). If the cell
   has no p (all-ties / no-test) it cannot be significant.
4. cell_verdict = nominal_cell if significant else TIED.
5. domain verdict over its metric cells (up to 2 metrics x 6 conditions,
   whatever exists): BETTER if n_better_cells > n_worse_cells, WORSE if
   n_worse_cells > n_better_cells, else TIED.

Both the per-domain verdict and the full cell-level evidence are recorded in
this table (one row per benchmark x comparator x domain plus supporting
count columns). DUAL_HEAD vs comparators with identical shared readouts
(UNSIGNED_GENE_PPR stage head; SIGNED_GENE_FIRST_PPR_SINGLE_HEAD gene head)
produce all-tied cells by construction (documented in 08).
"""
import pathlib

import numpy as np
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[2]
STAGE = ROOT / "data" / "processed" / "results"
RES = STAGE / "results_confirmatory_simulation"
TIE_EPS = 1e-12

DOMAIN_METRICS = {
    "developmental_alignment": ["alignment_hit_top1", "true_stage_mass"],
    "gene_prioritization": ["gene_auprc", "precision_at_50"]}
BENCHES = ["step02A", "step02B"]
COMPARATORS = [
    "SIGNED_GENE_FIRST_PPR_SINGLE_HEAD", "UNSIGNED_GENE_PPR",
    "DIRECT_NATIVE_PROJECTION", "GENE_ONLY_PPR_BASELINE",
    "HETERO_PPR_STANDARD", "STAGE_CORRELATION", "GRAPH_OT", "SIMPLE_OT",
    "MEAN_SIGNATURE", "RANDOM", "SIMULATION_PROXY"]


def main():
    pw = pd.read_csv(RES / "05_confirmatory_paired_deltas.tsv", sep="\t")
    rows = []
    for bench in BENCHES:
        b = pw[(pw["benchmark"] == bench) &
               (pw["condition"] != "ALL_CONDITIONS")]
        for comparator in COMPARATORS:
            sub = b[(b["baseline_method"] == comparator)]
            for domain, metrics in DOMAIN_METRICS.items():
                cells = sub[sub["metric"].isin(metrics)]
                n_total = 0
                n_better = 0
                n_worse = 0
                n_tied = 0
                cell_notes = []
                for metric in metrics:
                    mc = cells[cells["metric"] == metric]
                    nb = nw = nt = 0
                    for _, r in mc.iterrows():
                        d = float(r["mean_delta_dualhead_minus_baseline"])
                        q = r.get("fdr_BH_within_family")
                        if not np.isfinite(d) or abs(d) <= TIE_EPS:
                            nt += 1
                            continue
                        nominal = ("BETTER" if d > 0 else "WORSE")
                        sig = (isinstance(q, float) and np.isfinite(q)
                               and q < 0.05)
                        if sig:
                            if nominal == "BETTER":
                                nb += 1
                            else:
                                nw += 1
                        else:
                            nt += 1
                    n_total += nb + nw + nt
                    n_better += nb
                    n_worse += nw
                    n_tied += nt
                    if nb or nw:
                        cell_notes.append(
                            f"{metric}: better {nb} / worse {nw} / tied {nt}")
                    else:
                        cell_notes.append(f"{metric}: all tied ({nt})")
                if n_better > n_worse:
                    verdict = "BETTER"
                elif n_worse > n_better:
                    verdict = "WORSE"
                else:
                    verdict = "TIED"
                rows.append(dict(
                    benchmark=bench, comparator=comparator, domain=domain,
                    domain_metrics=",".join(DOMAIN_METRICS[domain]),
                    n_cells_classified=n_total,
                    n_cells_dual_better=n_better,
                    n_cells_comparator_better=n_worse,
                    n_cells_tied=n_tied,
                    domain_verdict_dual_vs_comparator=verdict,
                    cell_rule=("cell = condition-level paired mean delta; "
                               "significant iff BH-FDR<0.05 within the "
                               "benchmark x metric family (05); "
                               "nominal-tie threshold |delta|<=1e-12"),
                    domain_rule=("BETTER if n_better_cells > "
                                 "n_worse_cells; WORSE if n_worse_cells > "
                                 "n_better_cells; else TIED"),
                    cells_detail="; ".join(cell_notes),
                    global_superiority_claim="NO - domain-level verdicts "
                                             "only (protocol Part J)"))
    out = pd.DataFrame(rows)
    out = out.sort_values(["benchmark", "domain", "comparator"])
    out.to_csv(RES / "09_dualhead_pareto_summary.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    print("[post] 09_dualhead_pareto_summary.tsv written")


if __name__ == "__main__":
    main()
