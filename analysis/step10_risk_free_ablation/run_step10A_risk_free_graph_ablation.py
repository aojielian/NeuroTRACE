#!/usr/bin/env python3

import math
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
PROJECT = BASE / "neurotrace_algorithm_project"

STEP04 = PROJECT / "step04_graph_construction" / "results"
STEP05 = PROJECT / "step05_optimal_transport_alignment" / "results"
STEP10 = PROJECT / "step10_risk_free_ablation"
OUT = STEP10 / "results"
FIG = STEP10 / "figures"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

NODES_FILE = STEP04 / "09_step04C_native_graph_nodes.tsv"
EDGES_FILE = STEP04 / "10_step04C_native_graph_edges.tsv.gz"
RISK_INCLUDED_TRANSPORT = STEP05 / "10_step05B_NTM_stage_graph_transport.tsv"
RISK_INCLUDED_TOP_TRANSPORT = STEP05 / "12_step05B_NTM_top_transport_stage.tsv"
RISK_INCLUDED_SFARI = STEP05 / "15_step05B_NTM_graph_priority_SFARI_enrichment.tsv"

ALPHA_RESTART = 0.35
MAX_ITER = 120
TOL = 1e-10
TRANSPORT_TEMP = 0.75
TOP_GENES_EXPORT = 500

EDGE_TYPE_SCALE = {
    "gene_gene_embedding_knn": 0.35,
    "module_gene_weight": 1.00,
    "stage_gene_profile": 0.80,
    "module_stage_native_alignment": 1.25,
    # risk_gene_prior intentionally excluded in risk-free graph
}

STAGE_ORDER = [
    "early_prenatal",
    "mid_prenatal",
    "late_prenatal",
    "childhood",
    "adolescence",
    "adulthood",
]

RISK_SETS = [
    "is_SFARI_all",
    "is_SFARI_S_plus_1",
    "is_SFARI_S_plus_1_2",
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


def parse_module_node(x):
    s = str(x)
    body = s.replace("module:", "")
    if "|top" in body:
        program, top_n = body.split("|top", 1)
    else:
        program, top_n = body, ""
    return program, top_n, body


def parse_stage_node(x):
    return str(x).replace("stage:", "")


def hypergeom_sf(k, K, n, N):
    if k <= 0:
        return 1.0
    upper = min(K, n)
    if k > upper:
        return 0.0

    def log_choose(a, b):
        if b < 0 or b > a:
            return -np.inf
        return math.lgamma(a + 1) - math.lgamma(b + 1) - math.lgamma(a - b + 1)

    logs = []
    denom = log_choose(N, n)
    for i in range(k, upper + 1):
        logs.append(log_choose(K, i) + log_choose(N - K, n - i) - denom)

    m = max(logs)
    return float(np.exp(m) * np.sum(np.exp(np.array(logs) - m)))


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    out = np.empty(n, dtype=float)
    prev = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        val = p[idx] * n / rank
        prev = min(prev, val)
        out[idx] = min(prev, 1.0)
    return out


def build_transition_matrix(nodes, edges):
    print(f"[{now()}] Building risk-free weighted transition matrix")

    node_ids = nodes["node_id"].astype(str).tolist()
    node_to_idx = {n: i for i, n in enumerate(node_ids)}
    n = len(node_ids)

    rows, cols, vals = [], [], []

    for _, e in edges.iterrows():
        src = str(e["source"])
        tgt = str(e["target"])
        et = str(e["edge_type"])

        if et == "risk_gene_prior":
            continue

        if src not in node_to_idx or tgt not in node_to_idx:
            continue

        base_w = float(e["weight"])
        scale = EDGE_TYPE_SCALE.get(et, 1.0)
        w = base_w * scale

        if not np.isfinite(w) or w <= 0:
            continue

        i = node_to_idx[src]
        j = node_to_idx[tgt]

        # Undirected PPR graph, same as Step05B except risk edges removed
        rows.extend([i, j])
        cols.extend([j, i])
        vals.extend([w, w])

    try:
        from scipy import sparse

        A = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=np.float64)
        row_sum = np.asarray(A.sum(axis=1)).ravel()
        inv = np.zeros_like(row_sum)
        inv[row_sum > 0] = 1.0 / row_sum[row_sum > 0]
        P = sparse.diags(inv) @ A
        return P, node_to_idx, node_ids, "scipy_sparse", len(vals) // 2

    except Exception:
        A = np.zeros((n, n), dtype=np.float64)
        for r, c, v in zip(rows, cols, vals):
            A[r, c] += v
        row_sum = A.sum(axis=1, keepdims=True)
        row_sum[row_sum == 0] = 1.0
        P = A / row_sum
        return P, node_to_idx, node_ids, "dense_numpy", len(vals) // 2


def personalized_pagerank(P, seed_idx, alpha=0.35, max_iter=120, tol=1e-10):
    n = P.shape[0]
    s = np.zeros(n, dtype=np.float64)
    s[seed_idx] = 1.0
    x = s.copy()

    sparse_like = not isinstance(P, np.ndarray)

    for it in range(max_iter):
        if sparse_like:
            x_new = alpha * s + (1.0 - alpha) * (P.T @ x)
        else:
            x_new = alpha * s + (1.0 - alpha) * (P.T.dot(x))

        delta = np.sum(np.abs(x_new - x))
        x = x_new
        if delta < tol:
            break

    return x, it + 1, delta


def safe_read(path):
    path = Path(path)
    if not path.exists():
        return pd.DataFrame()
    if str(path).endswith(".gz"):
        return pd.read_csv(path, sep="\t", compression="gzip")
    return pd.read_csv(path, sep="\t")


def main():
    print(f"[{now()}] Step10A risk-free graph ablation started")

    for f in [NODES_FILE, EDGES_FILE]:
        if not Path(f).exists():
            raise FileNotFoundError(f"Missing required input: {f}")

    nodes = pd.read_csv(NODES_FILE, sep="\t")
    edges = pd.read_csv(EDGES_FILE, sep="\t", compression="gzip")

    legacy_hit = 0
    for col in ["node_id", "label"]:
        if col in nodes.columns:
            legacy_hit += nodes[col].astype(str).str.contains("Prog12|Prog13|Prog14|Prog7", regex=True).sum()

    for col in ["source", "target"]:
        if col in edges.columns:
            legacy_hit += edges[col].astype(str).str.contains("Prog12|Prog13|Prog14|Prog7", regex=True).sum()

    if legacy_hit > 0:
        raise RuntimeError(f"Legacy Prog terms detected in native graph: {legacy_hit}")

    n_edges_original = edges.shape[0]
    n_risk_edges_removed = int((edges["edge_type"] == "risk_gene_prior").sum())
    edges_no_risk = edges[edges["edge_type"] != "risk_gene_prior"].copy()

    print(f"[{now()}] Original edges: {n_edges_original}")
    print(f"[{now()}] Removed risk_gene_prior edges: {n_risk_edges_removed}")
    print(f"[{now()}] Risk-free edges: {edges_no_risk.shape[0]}")

    P, node_to_idx, node_ids, matrix_backend, n_graph_edges_used = build_transition_matrix(nodes, edges_no_risk)

    module_nodes = nodes[nodes["node_type"] == "native_module"]["node_id"].astype(str).tolist()
    stage_nodes = nodes[nodes["node_type"] == "developmental_stage"]["node_id"].astype(str).tolist()
    gene_nodes = nodes[nodes["node_type"] == "gene"]["node_id"].astype(str).tolist()

    print(f"[{now()}] Matrix backend: {matrix_backend}")
    print(f"[{now()}] Native module nodes: {len(module_nodes)}")
    print(f"[{now()}] Stage nodes: {len(stage_nodes)}")
    print(f"[{now()}] Gene nodes: {len(gene_nodes)}")

    # Direct native alignment remains available through module-stage edges
    direct_map = {}
    conf_map = {}

    ms_edges = edges_no_risk[edges_no_risk["edge_type"] == "module_stage_native_alignment"].copy()
    for _, r in ms_edges.iterrows():
        direct_map[(str(r["source"]), str(r["target"]))] = float(r["signed_weight"])
        conf_map[(str(r["source"]), str(r["target"]))] = str(r.get("alignment_confidence", ""))

    gene_meta = nodes[nodes["node_type"] == "gene"].copy().set_index("node_id")

    transport_rows = []
    gene_priority_rows = []
    runtime_rows = []

    for mnode in module_nodes:
        midx = node_to_idx[mnode]
        program, top_n, module_id = parse_module_node(mnode)

        score, n_iter, delta = personalized_pagerank(
            P,
            seed_idx=midx,
            alpha=ALPHA_RESTART,
            max_iter=MAX_ITER,
            tol=TOL,
        )

        runtime_rows.append({
            "module_node": mnode,
            "program": program,
            "top_n": top_n,
            "module_id": module_id,
            "n_iter": n_iter,
            "final_delta": delta,
        })

        # Module-stage transport
        stage_records = []
        for snode in stage_nodes:
            sidx = node_to_idx[snode]
            state = parse_stage_node(snode)
            stage_records.append({
                "module_node": mnode,
                "program": program,
                "top_n": top_n,
                "module_id": module_id,
                "stage_node": snode,
                "state_id": state,
                "stage_order": STAGE_ORDER.index(state) + 1 if state in STAGE_ORDER else 999,
                "graph_diffusion_score": float(score[sidx]),
                "direct_native_alignment_score": direct_map.get((mnode, snode), np.nan),
                "alignment_confidence": conf_map.get((mnode, snode), ""),
            })

        sdf = pd.DataFrame(stage_records)
        sdf["graph_diffusion_z"] = zscore(sdf["graph_diffusion_score"])
        sdf["direct_native_z"] = zscore(sdf["direct_native_alignment_score"])

        sdf["neurotrace_no_risk_transport_score"] = (
            0.55 * sdf["direct_native_z"] +
            0.45 * sdf["graph_diffusion_z"]
        )

        sdf["transport_probability"] = softmax(
            sdf["neurotrace_no_risk_transport_score"].values,
            temp=TRANSPORT_TEMP,
        )
        sdf["transport_rank"] = sdf["transport_probability"].rank(ascending=False, method="min")
        transport_rows.append(sdf)

        # Module-gene graph priority
        gene_records = []
        for gnode in gene_nodes:
            gidx = node_to_idx[gnode]
            meta = gene_meta.loc[gnode]
            gene_records.append({
                "module_node": mnode,
                "program": program,
                "top_n": top_n,
                "module_id": module_id,
                "gene_node": gnode,
                "gene": str(meta["label"]),
                "graph_priority_score": float(score[gidx]),
                "is_SFARI_all": int(meta.get("is_SFARI_all", 0)),
                "is_SFARI_S_plus_1": int(meta.get("is_SFARI_S_plus_1", 0)),
                "is_SFARI_S_plus_1_2": int(meta.get("is_SFARI_S_plus_1_2", 0)),
            })

        gdf = pd.DataFrame(gene_records)
        gdf["priority_rank"] = gdf["graph_priority_score"].rank(ascending=False, method="first")
        gdf = gdf.sort_values("priority_rank")
        gene_priority_rows.append(gdf)

    transport = pd.concat(transport_rows, ignore_index=True)
    transport = transport.sort_values(["program", "top_n", "transport_rank", "stage_order"])
    transport.to_csv(OUT / "91_step10A_no_risk_NTM_stage_graph_transport.tsv", sep="\t", index=False)

    top_stage = transport[transport["transport_rank"] == 1].copy()

    def window_label(stage):
        if stage in ["early_prenatal", "mid_prenatal", "late_prenatal"]:
            return "prenatal"
        if stage in ["childhood", "adolescence"]:
            return "postnatal_developmental"
        if stage == "adulthood":
            return "adult_like"
        return "unknown"

    top_stage["developmental_window"] = top_stage["state_id"].map(window_label)
    top_stage = top_stage.sort_values(["program", "top_n"])
    top_stage.to_csv(OUT / "92_step10A_no_risk_top_transport_stage.tsv", sep="\t", index=False)

    gene_priority = pd.concat(gene_priority_rows, ignore_index=True)
    gene_priority.to_csv(OUT / "93_step10A_no_risk_gene_graph_priority.tsv.gz", sep="\t", index=False, compression="gzip")

    top500 = gene_priority[gene_priority["priority_rank"] <= TOP_GENES_EXPORT].copy()
    top500.to_csv(OUT / "94_step10A_no_risk_top500_graph_priority_genes.tsv.gz", sep="\t", index=False, compression="gzip")

    # Post hoc SFARI enrichment after no-risk diffusion
    universe_meta = gene_priority.drop_duplicates("gene_node")
    universe_n = universe_meta.shape[0]

    risk_rows = []
    for risk_col in RISK_SETS:
        K = int(universe_meta[risk_col].sum())
        for (program, top_n), sub in gene_priority.groupby(["program", "top_n"]):
            sub = sub.sort_values("priority_rank")
            for cutoff in [50, 100, 200, 500]:
                top = sub.head(cutoff)
                k = int(top[risk_col].sum())
                n = int(top.shape[0])
                N = int(universe_n)

                pval = hypergeom_sf(k, K, n, N)

                a = k
                b = n - k
                c = K - k
                d = N - K - b
                OR = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

                risk_rows.append({
                    "program": program,
                    "top_n": top_n,
                    "risk_set": risk_col,
                    "priority_cutoff": cutoff,
                    "overlap_n": k,
                    "top_n_genes": n,
                    "risk_genes_in_universe": K,
                    "universe_n": N,
                    "odds_ratio": OR,
                    "p_value": pval,
                    "analysis": "no_risk_graph_posthoc",
                })

    risk_enrich = pd.DataFrame(risk_rows)
    risk_enrich["fdr"] = bh_fdr(risk_enrich["p_value"].values)
    risk_enrich = risk_enrich.sort_values(["fdr", "p_value", "program", "top_n", "risk_set", "priority_cutoff"])
    risk_enrich.to_csv(OUT / "95_step10A_no_risk_posthoc_SFARI_enrichment.tsv", sep="\t", index=False)

    # Compare transport with risk-included Step05B
    risk_transport = safe_read(RISK_INCLUDED_TRANSPORT)
    if not risk_transport.empty:
        rt = risk_transport[[
            "program", "top_n", "state_id",
            "transport_probability",
            "neurotrace_native_transport_score",
            "graph_diffusion_score",
            "direct_native_alignment_score",
            "alignment_confidence",
            "transport_rank",
        ]].copy()
        rt = rt.rename(columns={
            "transport_probability": "risk_included_transport_probability",
            "neurotrace_native_transport_score": "risk_included_transport_score",
            "graph_diffusion_score": "risk_included_graph_diffusion_score",
            "direct_native_alignment_score": "risk_included_direct_alignment_score",
            "alignment_confidence": "risk_included_alignment_confidence",
            "transport_rank": "risk_included_transport_rank",
        })

        nr = transport[[
            "program", "top_n", "state_id",
            "transport_probability",
            "neurotrace_no_risk_transport_score",
            "graph_diffusion_score",
            "direct_native_alignment_score",
            "alignment_confidence",
            "transport_rank",
        ]].copy()
        nr = nr.rename(columns={
            "transport_probability": "no_risk_transport_probability",
            "neurotrace_no_risk_transport_score": "no_risk_transport_score",
            "graph_diffusion_score": "no_risk_graph_diffusion_score",
            "direct_native_alignment_score": "no_risk_direct_alignment_score",
            "alignment_confidence": "no_risk_alignment_confidence",
            "transport_rank": "no_risk_transport_rank",
        })

        # Ensure merge keys have identical types across Step05B and Step10A outputs
        for df in (rt, nr):
            df["program"] = df["program"].astype(str)
            df["top_n"] = df["top_n"].astype(str)
            df["state_id"] = df["state_id"].astype(str)

        comp = rt.merge(nr, on=["program", "top_n", "state_id"], how="outer")
        comp["delta_no_risk_minus_risk_transport_probability"] = (
            comp["no_risk_transport_probability"] - comp["risk_included_transport_probability"]
        )
        comp.to_csv(OUT / "96_step10A_risk_vs_no_risk_transport_comparison.tsv", sep="\t", index=False)

    # Compare SFARI enrichment with risk-included Step05B
    risk_sfari = safe_read(RISK_INCLUDED_SFARI)
    if not risk_sfari.empty:
        risk_sfari = risk_sfari.copy()
        risk_sfari["analysis"] = "risk_included_graph"
        risk_sfari = risk_sfari[[
            "program", "top_n", "risk_set", "priority_cutoff",
            "overlap_n", "top_n_genes", "risk_genes_in_universe",
            "universe_n", "odds_ratio", "p_value", "fdr", "analysis"
        ]]

        no_risk_sfari = risk_enrich[[
            "program", "top_n", "risk_set", "priority_cutoff",
            "overlap_n", "top_n_genes", "risk_genes_in_universe",
            "universe_n", "odds_ratio", "p_value", "fdr", "analysis"
        ]]

        sfari_long = pd.concat([risk_sfari, no_risk_sfari], ignore_index=True)
        sfari_long.to_csv(OUT / "97_step10A_risk_vs_no_risk_SFARI_enrichment_long.tsv", sep="\t", index=False)

        # Ensure merge keys have identical types across Step05B and Step10A enrichment outputs
        for df in (risk_sfari, no_risk_sfari):
            df["program"] = df["program"].astype(str)
            df["top_n"] = df["top_n"].astype(str)
            df["risk_set"] = df["risk_set"].astype(str)
            df["priority_cutoff"] = df["priority_cutoff"].astype(str)

        wide = risk_sfari.merge(
            no_risk_sfari,
            on=["program", "top_n", "risk_set", "priority_cutoff"],
            how="outer",
            suffixes=("_risk", "_no_risk")
        )

        wide["delta_log10FDR_no_risk_minus_risk"] = (
            -np.log10(np.maximum(wide["fdr_no_risk"].astype(float), 1e-300)) -
            -np.log10(np.maximum(wide["fdr_risk"].astype(float), 1e-300))
        )
        wide["delta_OR_no_risk_minus_risk"] = wide["odds_ratio_no_risk"] - wide["odds_ratio_risk"]

        wide.to_csv(OUT / "98_step10A_risk_vs_no_risk_SFARI_enrichment_comparison.tsv", sep="\t", index=False)

    pd.DataFrame(runtime_rows).to_csv(OUT / "99_step10A_no_risk_diffusion_runtime.tsv", sep="\t", index=False)

    manifest_rows = [
        {"item": "n_nodes", "value": nodes.shape[0]},
        {"item": "n_edges_original", "value": n_edges_original},
        {"item": "n_risk_gene_prior_edges_removed", "value": n_risk_edges_removed},
        {"item": "n_edges_no_risk", "value": edges_no_risk.shape[0]},
        {"item": "n_graph_edges_used_undirected_unique", "value": n_graph_edges_used},
        {"item": "n_module_nodes", "value": len(module_nodes)},
        {"item": "n_stage_nodes", "value": len(stage_nodes)},
        {"item": "n_gene_nodes", "value": len(gene_nodes)},
        {"item": "matrix_backend", "value": matrix_backend},
        {"item": "alpha_restart", "value": ALPHA_RESTART},
        {"item": "transport_temp", "value": TRANSPORT_TEMP},
        {"item": "legacy_prog_hit_count", "value": int(legacy_hit)},
    ]
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(OUT / "100_step10A_no_risk_manifest.tsv", sep="\t", index=False)

    # Summary interpretation
    top_sig = risk_enrich.head(20)

    n_no_risk_main_sig = risk_enrich[
        (risk_enrich["fdr"] < 0.05) &
        (risk_enrich["top_n"].astype(int).isin([200, 500]))
    ].shape[0]

    with open(OUT / "101_step10A_risk_free_graph_ablation_summary.md", "w") as f:
        f.write("# NeuroTRACE Step10A risk-free graph ablation summary\n\n")
        f.write(f"Generated: {now()}\n\n")

        f.write("## Purpose\n")
        f.write("This analysis removes all risk-gene prior edges from the native heterogeneous graph, reruns graph diffusion, and then evaluates SFARI enrichment only as a post hoc endpoint. This directly tests whether SFARI enrichment among graph-prioritized genes is driven by circular inclusion of SFARI risk nodes.\n\n")

        f.write("## Graph ablation\n")
        f.write(f"- Original edges: {n_edges_original}\n")
        f.write(f"- Removed risk_gene_prior edges: {n_risk_edges_removed}\n")
        f.write(f"- Risk-free edges: {edges_no_risk.shape[0]}\n")
        f.write(f"- Matrix backend: {matrix_backend}\n")
        f.write(f"- Legacy Prog hit count: {legacy_hit}\n\n")

        f.write("## No-risk top transport stage\n\n")
        f.write(top_stage[[
            "program", "top_n", "state_id", "developmental_window",
            "transport_probability",
            "neurotrace_no_risk_transport_score",
            "direct_native_alignment_score",
            "graph_diffusion_score",
            "alignment_confidence"
        ]].to_markdown(index=False))
        f.write("\n\n")

        f.write("## No-risk post hoc SFARI enrichment: top rows\n\n")
        f.write(top_sig.to_markdown(index=False))
        f.write("\n\n")

        f.write("## Main top_n no-risk significant SFARI tests\n")
        f.write(f"- Number of FDR<0.05 post hoc SFARI tests among top200/top500 modules: {n_no_risk_main_sig}\n\n")

        f.write("## Interpretation\n")
        if n_no_risk_main_sig > 0:
            f.write("SFARI enrichment persists after removing risk-gene prior edges from graph diffusion. This supports the interpretation that genetic-risk convergence is not solely an artifact of SFARI priors embedded in the graph.\n")
        else:
            f.write("SFARI enrichment does not persist after removing risk-gene prior edges. Genetic-risk convergence should be interpreted as prior-guided ranking rather than independent validation.\n")

    print(f"[{now()}] Step10A done")
    print(f"[{now()}] No-risk main significant SFARI tests: {n_no_risk_main_sig}")
    print(f"[{now()}] Results: {OUT}")


if __name__ == "__main__":
    main()
