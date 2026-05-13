#!/usr/bin/env python3

import math
import gzip
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
PROJECT = BASE / "neurotrace_algorithm_project"

STEP04 = PROJECT / "step04_graph_construction" / "results"
STEP05 = PROJECT / "step05_optimal_transport_alignment"
OUT = STEP05 / "results"
FIG = STEP05 / "figures"

OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

NODES_FILE = STEP04 / "01_step04A_graph_nodes.tsv"
EDGES_FILE = STEP04 / "02_step04A_graph_edges.tsv.gz"
PROGRAM_STAGE_EDGES_FILE = STEP04 / "06_step04A_program_stage_alignment_edges.tsv"

ALPHA_RESTART = 0.35
MAX_ITER = 120
TOL = 1e-10
TRANSPORT_TEMP = 0.75
TOP_GENES_EXPORT = 500

EDGE_TYPE_SCALE = {
    "gene_gene_embedding_knn": 0.35,
    "program_gene_weight": 1.00,
    "stage_gene_profile": 0.80,
    "program_stage_refined_alignment": 1.25,
    "risk_gene_prior": 0.30
}

STAGE_ORDER = [
    "early_prenatal",
    "mid_prenatal",
    "late_prenatal",
    "childhood",
    "adolescence",
    "adulthood"
]

RISK_COLS = [
    "is_SFARI_all",
    "is_SFARI_S_plus_1",
    "is_SFARI_S_plus_1_2"
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


def parse_program_node(x):
    s = str(x)
    if not s.startswith("program:"):
        return "", ""
    body = s.replace("program:", "")
    if "|top" in body:
        program, topn = body.split("|top", 1)
    else:
        program, topn = body, ""
    return program, topn


def parse_stage_node(x):
    s = str(x)
    return s.replace("stage:", "") if s.startswith("stage:") else ""


def hypergeom_sf(k, K, n, N):
    """
    P[X >= k] for X ~ Hypergeom(N population, K success, n draws).
    Uses log-combination to avoid scipy dependency.
    """
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


def benjamini_hochberg(pvals):
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
    print(f"[{now()}] Building weighted adjacency")

    node_ids = nodes["node_id"].astype(str).tolist()
    node_to_idx = {n: i for i, n in enumerate(node_ids)}
    n = len(node_ids)

    rows = []
    cols = []
    vals = []

    for _, e in edges.iterrows():
        src = str(e["source"])
        tgt = str(e["target"])
        et = str(e["edge_type"])
        if src not in node_to_idx or tgt not in node_to_idx:
            continue

        base_w = float(e["weight"])
        scale = EDGE_TYPE_SCALE.get(et, 1.0)
        w = base_w * scale
        if not np.isfinite(w) or w <= 0:
            continue

        i = node_to_idx[src]
        j = node_to_idx[tgt]

        # undirected graph for first graph-diffusion MVP
        rows.extend([i, j])
        cols.extend([j, i])
        vals.extend([w, w])

    try:
        from scipy import sparse
        A = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=np.float64)
        row_sum = np.asarray(A.sum(axis=1)).ravel()
        inv = np.zeros_like(row_sum)
        inv[row_sum > 0] = 1.0 / row_sum[row_sum > 0]
        Dinv = sparse.diags(inv)
        P = Dinv @ A
        return P, node_to_idx, node_ids, "scipy_sparse"
    except Exception:
        A = np.zeros((n, n), dtype=np.float64)
        for r, c, v in zip(rows, cols, vals):
            A[r, c] += v
        row_sum = A.sum(axis=1, keepdims=True)
        row_sum[row_sum == 0] = 1.0
        P = A / row_sum
        return P, node_to_idx, node_ids, "dense_numpy"


def personalized_pagerank(P, seed_idx, alpha=0.35, max_iter=120, tol=1e-10):
    n = P.shape[0]
    s = np.zeros(n, dtype=np.float64)
    s[seed_idx] = 1.0
    x = s.copy()

    use_sparse = hasattr(P, "T") and not isinstance(P, np.ndarray)

    for it in range(max_iter):
        if use_sparse:
            x_new = alpha * s + (1.0 - alpha) * (P.T @ x)
        else:
            x_new = alpha * s + (1.0 - alpha) * (P.T.dot(x))

        delta = np.sum(np.abs(x_new - x))
        x = x_new
        if delta < tol:
            break

    return x, it + 1, delta


def main():
    print(f"[{now()}] Step05A graph-informed transport started")

    for f in [NODES_FILE, EDGES_FILE, PROGRAM_STAGE_EDGES_FILE]:
        if not f.exists():
            raise FileNotFoundError(f"Missing required input: {f}")

    nodes = pd.read_csv(NODES_FILE, sep="\t")
    edges = pd.read_csv(EDGES_FILE, sep="\t")
    ps_edges = pd.read_csv(PROGRAM_STAGE_EDGES_FILE, sep="\t")

    P, node_to_idx, node_ids, matrix_backend = build_transition_matrix(nodes, edges)

    program_nodes = nodes[nodes["node_type"] == "program"]["node_id"].astype(str).tolist()
    stage_nodes = nodes[nodes["node_type"] == "developmental_stage"]["node_id"].astype(str).tolist()
    gene_nodes = nodes[nodes["node_type"] == "gene"]["node_id"].astype(str).tolist()

    print(f"[{now()}] Matrix backend: {matrix_backend}")
    print(f"[{now()}] Program nodes: {len(program_nodes)}")
    print(f"[{now()}] Stage nodes: {len(stage_nodes)}")
    print(f"[{now()}] Gene nodes: {len(gene_nodes)}")

    # direct Step03B score from program-stage edge signed_weight
    direct_map = {}
    for _, r in ps_edges.iterrows():
        direct_map[(str(r["source"]), str(r["target"]))] = float(r["signed_weight"])

    transport_rows = []
    gene_priority_rows = []
    runtime_rows = []

    gene_meta = nodes[nodes["node_type"] == "gene"].copy()
    gene_meta = gene_meta.set_index("node_id")

    for pnode in program_nodes:
        pidx = node_to_idx[pnode]
        program, top_n = parse_program_node(pnode)

        score, n_iter, delta = personalized_pagerank(
            P,
            seed_idx=pidx,
            alpha=ALPHA_RESTART,
            max_iter=MAX_ITER,
            tol=TOL
        )

        runtime_rows.append({
            "program_node": pnode,
            "program": program,
            "top_n": top_n,
            "n_iter": n_iter,
            "final_delta": delta
        })

        # Stage transport
        stage_raw = []
        for snode in stage_nodes:
            sidx = node_to_idx[snode]
            state = parse_stage_node(snode)
            stage_raw.append({
                "program_node": pnode,
                "program": program,
                "top_n": top_n,
                "stage_node": snode,
                "state_id": state,
                "stage_order": STAGE_ORDER.index(state) + 1 if state in STAGE_ORDER else 999,
                "graph_diffusion_score": float(score[sidx]),
                "direct_refined_alignment_score": direct_map.get((pnode, snode), np.nan)
            })

        sdf = pd.DataFrame(stage_raw)
        sdf["graph_diffusion_z"] = zscore(sdf["graph_diffusion_score"].values)
        sdf["direct_refined_z"] = zscore(sdf["direct_refined_alignment_score"].values)
        sdf["neurotrace_transport_score"] = (
            0.55 * sdf["direct_refined_z"] +
            0.45 * sdf["graph_diffusion_z"]
        )
        sdf["transport_probability"] = softmax(sdf["neurotrace_transport_score"].values, temp=TRANSPORT_TEMP)
        sdf["transport_rank"] = sdf["transport_probability"].rank(ascending=False, method="min")
        transport_rows.append(sdf)

        # Gene priority
        gene_records = []
        for gnode in gene_nodes:
            gidx = node_to_idx[gnode]
            meta = gene_meta.loc[gnode]
            gene_records.append({
                "program_node": pnode,
                "program": program,
                "top_n": top_n,
                "gene_node": gnode,
                "gene": str(meta["label"]),
                "graph_priority_score": float(score[gidx]),
                "is_SFARI_all": int(meta.get("is_SFARI_all", 0)),
                "is_SFARI_S_plus_1": int(meta.get("is_SFARI_S_plus_1", 0)),
                "is_SFARI_S_plus_1_2": int(meta.get("is_SFARI_S_plus_1_2", 0))
            })

        gdf = pd.DataFrame(gene_records)
        gdf["priority_rank"] = gdf["graph_priority_score"].rank(ascending=False, method="first")
        gdf = gdf.sort_values("priority_rank")
        gene_priority_rows.append(gdf)

    transport = pd.concat(transport_rows, ignore_index=True)
    transport = transport.sort_values(["program", "top_n", "transport_rank", "stage_order"])
    transport.to_csv(OUT / "01_step05A_program_stage_graph_transport.tsv", sep="\t", index=False)

    wide = transport.pivot_table(
        index=["program", "top_n"],
        columns="state_id",
        values="transport_probability",
        aggfunc="mean"
    ).reset_index()
    wide.to_csv(OUT / "02_step05A_transport_probability_matrix_wide.tsv", sep="\t", index=False)

    top_stage = transport[transport["transport_rank"] == 1].copy()
    top_stage = top_stage.sort_values(["program", "top_n"])
    top_stage.to_csv(OUT / "03_step05A_top_transport_stage_per_program.tsv", sep="\t", index=False)

    gene_priority = pd.concat(gene_priority_rows, ignore_index=True)
    gene_priority.to_csv(OUT / "04_step05A_program_gene_graph_priority.tsv.gz", sep="\t", index=False, compression="gzip")

    top_genes = gene_priority[gene_priority["priority_rank"] <= TOP_GENES_EXPORT].copy()
    top_genes.to_csv(OUT / "05_step05A_top500_graph_priority_genes_per_program.tsv.gz", sep="\t", index=False, compression="gzip")

    # Risk enrichment among top N genes
    risk_rows = []
    universe_n = gene_priority[["gene_node"]].drop_duplicates().shape[0]

    risk_sets = [
        "is_SFARI_all",
        "is_SFARI_S_plus_1",
        "is_SFARI_S_plus_1_2"
    ]

    universe_meta = gene_priority.drop_duplicates("gene_node")
    for risk_col in risk_sets:
        K = int(universe_meta[risk_col].sum())
        for (program, top_n), sub in gene_priority.groupby(["program", "top_n"]):
            sub = sub.sort_values("priority_rank")
            for cutoff in [50, 100, 200, 500]:
                top = sub.head(cutoff)
                k = int(top[risk_col].sum())
                n = int(top.shape[0])
                N = int(universe_n)

                pval = hypergeom_sf(k, K, n, N)

                # odds ratio with 0.5 correction
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
                    "p_value": pval
                })

    risk_enrich = pd.DataFrame(risk_rows)
    risk_enrich["fdr"] = benjamini_hochberg(risk_enrich["p_value"].values)
    risk_enrich = risk_enrich.sort_values(["fdr", "p_value", "program", "top_n"])
    risk_enrich.to_csv(OUT / "06_step05A_graph_priority_risk_enrichment.tsv", sep="\t", index=False)

    pd.DataFrame(runtime_rows).to_csv(OUT / "07_step05A_diffusion_runtime.tsv", sep="\t", index=False)

    # Manifest
    manifest_rows = [
        {"item": "n_program_nodes", "value": len(program_nodes)},
        {"item": "n_stage_nodes", "value": len(stage_nodes)},
        {"item": "n_gene_nodes", "value": len(gene_nodes)},
        {"item": "n_transport_rows", "value": int(transport.shape[0])},
        {"item": "n_gene_priority_rows", "value": int(gene_priority.shape[0])},
        {"item": "matrix_backend", "value": matrix_backend},
        {"item": "alpha_restart", "value": ALPHA_RESTART},
        {"item": "max_iter", "value": MAX_ITER},
        {"item": "transport_temp", "value": TRANSPORT_TEMP}
    ]
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(OUT / "08_step05A_transport_freeze_manifest.tsv", sep="\t", index=False)

    # Plots
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        heat = transport.pivot_table(
            index=["program", "top_n"],
            columns="state_id",
            values="transport_probability",
            aggfunc="mean"
        )
        ordered_cols = [s for s in STAGE_ORDER if s in heat.columns]
        heat = heat[ordered_cols]

        plt.figure(figsize=(max(8, 0.8 * heat.shape[1]), max(5, 0.35 * heat.shape[0])))
        plt.imshow(heat.values, aspect="auto")
        plt.xticks(range(heat.shape[1]), heat.columns, rotation=45, ha="right")
        plt.yticks(range(heat.shape[0]), [f"{a}|top{b}" for a, b in heat.index])
        plt.colorbar(label="Transport probability")
        plt.title("Step05A graph-informed program-to-stage transport")
        plt.tight_layout()
        plt.savefig(FIG / "01_step05A_program_stage_transport_heatmap.pdf")
        plt.close()

        best = top_stage.copy()
        labels = best["program"].astype(str) + "|top" + best["top_n"].astype(str)
        plt.figure(figsize=(10, 4))
        plt.bar(labels, best["transport_probability"])
        plt.xticks(rotation=45, ha="right")
        plt.ylabel("Top-stage transport probability")
        plt.title("Step05A top transport probability per program")
        plt.tight_layout()
        plt.savefig(FIG / "02_step05A_top_stage_transport_probability.pdf")
        plt.close()

    except Exception as e:
        with open(OUT / "Step05A_plotting_failed.txt", "w") as f:
            f.write(str(e) + "\n")

    with open(OUT / "09_step05A_graph_informed_transport_summary.md", "w") as f:
        f.write("# NeuroTRACE Step05A graph-informed transport summary\n\n")
        f.write(f"Generated: {now()}\n\n")

        f.write("## Inputs\n")
        f.write(f"- Nodes: `{NODES_FILE}`\n")
        f.write(f"- Edges: `{EDGES_FILE}`\n")
        f.write(f"- Program-stage refined edges: `{PROGRAM_STAGE_EDGES_FILE}`\n\n")

        f.write("## Configuration\n")
        f.write(f"- Restart alpha: {ALPHA_RESTART}\n")
        f.write(f"- Max iterations: {MAX_ITER}\n")
        f.write(f"- Transport softmax temperature: {TRANSPORT_TEMP}\n")
        f.write(f"- Matrix backend: {matrix_backend}\n\n")

        f.write("## Graph-diffusion outputs\n")
        f.write(f"- Program nodes: {len(program_nodes)}\n")
        f.write(f"- Stage nodes: {len(stage_nodes)}\n")
        f.write(f"- Gene nodes: {len(gene_nodes)}\n")
        f.write(f"- Transport rows: {transport.shape[0]}\n")
        f.write(f"- Gene-priority rows: {gene_priority.shape[0]}\n\n")

        f.write("## Top transport stage per program\n\n")
        f.write(top_stage[[
            "program", "top_n", "state_id",
            "transport_probability",
            "neurotrace_transport_score",
            "direct_refined_alignment_score",
            "graph_diffusion_score"
        ]].to_markdown(index=False))
        f.write("\n\n")

        f.write("## Interpretation\n")
        f.write("Step05A upgrades the Step03B static program-stage alignment into a graph-informed transport map by combining direct refined alignment with graph diffusion through program-gene, stage-gene, gene-gene, and risk-gene edges. This is the first operational NeuroTRACE transport output on real data.\n")

    print(f"[{now()}] Step05A done")
    print(f"[{now()}] Program nodes: {len(program_nodes)}")
    print(f"[{now()}] Stage nodes: {len(stage_nodes)}")
    print(f"[{now()}] Gene nodes: {len(gene_nodes)}")
    print(f"[{now()}] Transport rows: {transport.shape[0]}")
    print(f"[{now()}] Gene-priority rows: {gene_priority.shape[0]}")
    print(f"[{now()}] Results: {OUT}")
    print(f"[{now()}] Figures: {FIG}")


if __name__ == "__main__":
    main()
