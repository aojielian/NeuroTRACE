#!/usr/bin/env python3

import os
import re
import csv
import gzip
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
PROJECT = BASE / "neurotrace_algorithm_project"
STEP = PROJECT / "step06_cross_cohort_validation"
OUT = STEP / "results"
OUT.mkdir(parents=True, exist_ok=True)

STEP01_MANIFEST = PROJECT / "step01_resource_manifest/results/01_neurotrace_step01_resource_manifest.tsv"

SCAN_ROOTS = [
    Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD"),
    Path("/gpfs/hpc/home/lijc/lianaoj/autism_scRNA"),
    Path("/gpfs/hpc/home/lijc/lianaoj/AAAY_data"),
]

DATASETS = ["GSE102741", "GSE64018"]

TEXT_EXT = (
    ".tsv", ".tsv.gz", ".csv", ".csv.gz",
    ".txt", ".txt.gz"
)

OBJECT_EXT = (
    ".rds", ".RDS", ".rdata", ".RData"
)

EXCLUDE_PATTERNS = [
    "/logs/", "/outs/", "/errs/", "/figures/",
    "/scripts/", "/tmp/", "/slurm",
    "module_and_devbridge_scores",
    "deconvolution_fractions",
    "selected_deconv_marker_genes",
    "step38_canonical_module_comparison/results_step38B",
    "step37_reference_deconvolution/results",
    "step42_incremental_value",
    "step43_exploratory",
    "neurotrace_algorithm_project/step02",
    "neurotrace_algorithm_project/step03",
    "neurotrace_algorithm_project/step04",
    "neurotrace_algorithm_project/step05",
]

POSITIVE_EXPR_TERMS = [
    "expr", "expression", "normalized", "normalize",
    "count", "counts", "matrix", "rsem", "fpkm", "tpm", "cpm",
    "voom", "logcpm", "gene"
]

POSITIVE_META_TERMS = [
    "meta", "metadata", "pheno", "phenotype",
    "clinical", "sample", "samples", "subject",
    "covariate", "covariates", "diagnosis"
]

METADATA_COL_TERMS = [
    "diagnosis", "dx", "condition", "group",
    "case", "control", "asd",
    "sample", "sample_id", "subject", "donor",
    "age", "sex", "gender", "rin", "pmi", "batch", "region"
]

EXPRESSION_FIRST_COL_TERMS = [
    "gene", "gene_symbol", "symbol", "feature", "feature_id",
    "ensembl", "id", "rowname", "rownames"
]


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def open_text(path):
    p = str(path)
    if p.endswith(".gz"):
        return gzip.open(p, "rt", errors="replace")
    return open(p, "rt", errors="replace")


def read_header(path):
    try:
        with open_text(path) as f:
            line = f.readline().rstrip("\n\r")
        return line
    except Exception:
        return ""


def split_header(header):
    if "\t" in header:
        return header.split("\t"), "\t"
    if "," in header:
        return header.split(","), ","
    return header.split(), None


def quick_count_rows(path, max_lines=200000):
    try:
        n = 0
        with open_text(path) as f:
            for _ in f:
                n += 1
                if n >= max_lines:
                    return f">={max_lines}"
        return max(0, n - 1)
    except Exception:
        return ""


def safe_stat(path):
    try:
        st = Path(path).stat()
        return {
            "exists": True,
            "readable": os.access(path, os.R_OK),
            "size_mb": round(st.st_size / 1024 / 1024, 3),
            "mtime": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S")
        }
    except Exception:
        return {
            "exists": False,
            "readable": False,
            "size_mb": "",
            "mtime": ""
        }


def read_manifest(path):
    if not Path(path).exists():
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(rows, path, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def is_excluded(path):
    low = str(path).lower()
    return any(p.lower() in low for p in EXCLUDE_PATTERNS)


def dataset_hit(path, dataset):
    low = str(path).lower()
    d = dataset.lower()
    return d in low or d.replace("gse", "gse_") in low


def role_guess(path, header_cols):
    low = str(path).lower()
    name = Path(path).name.lower()
    cols_low = [c.lower() for c in header_cols]

    expr_name_score = sum(t in low for t in POSITIVE_EXPR_TERMS)
    meta_name_score = sum(t in low for t in POSITIVE_META_TERMS)

    n_cols = len(header_cols)

    metadata_col_score = sum(any(t in c for t in METADATA_COL_TERMS) for c in cols_low)
    first_col = cols_low[0] if cols_low else ""
    first_col_expr = any(t in first_col for t in EXPRESSION_FIRST_COL_TERMS)

    # expression matrix usually has many columns and gene-like first column
    if n_cols >= 20 and (first_col_expr or expr_name_score >= 1) and metadata_col_score <= 5:
        return "expression_matrix_candidate"

    # metadata has clinical/dx/sample terms and not too many sample columns
    if metadata_col_score >= 2 or meta_name_score >= 1:
        return "metadata_candidate"

    if n_cols >= 20 and expr_name_score >= 1:
        return "expression_matrix_candidate"

    if str(path).endswith(OBJECT_EXT):
        return "serialized_object_candidate"

    return "uncertain_candidate"


def candidate_priority(path, dataset, role, n_cols, n_rows):
    low = str(path).lower()
    score = 0

    if dataset.lower() in low:
        score += 30

    if role == "expression_matrix_candidate":
        score += 30
    elif role == "metadata_candidate":
        score += 25
    elif role == "serialized_object_candidate":
        score += 20

    if any(t in low for t in ["normalized", "expression", "expr", "matrix", "count", "counts"]):
        score += 10

    if any(t in low for t in ["meta", "metadata", "pheno", "sample", "clinical", "covariate"]):
        score += 10

    if "raw" in low:
        score += 2
    if "old" in low or "backup" in low or "test" in low:
        score -= 10

    try:
        nc = int(n_cols)
        if nc >= 20:
            score += min(20, nc // 50)
    except Exception:
        pass

    if isinstance(n_rows, int) and n_rows >= 1000:
        score += 10

    return score


def scan_filesystem_candidates():
    rows = []

    for root in SCAN_ROOTS:
        if not root.exists():
            continue

        for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
            # skip noise directories early
            lowdir = str(dirpath).lower()
            if any(x in lowdir for x in ["/logs", "/outs", "/errs", "/figures", "/tmp"]):
                continue

            for fn in filenames:
                p = Path(dirpath) / fn
                plow = str(p).lower()

                if is_excluded(p):
                    continue

                if not (plow.endswith(TEXT_EXT) or plow.endswith(tuple(e.lower() for e in OBJECT_EXT))):
                    continue

                matched_datasets = [d for d in DATASETS if dataset_hit(p, d)]
                if not matched_datasets:
                    continue

                st = safe_stat(p)

                header = ""
                cols = []
                n_cols = ""
                n_rows = ""
                header_status = "not_text_checked"

                if plow.endswith(TEXT_EXT):
                    header = read_header(p)
                    cols, _ = split_header(header)
                    cols = [c.strip() for c in cols]
                    n_cols = len(cols)
                    header_status = "ok" if header else "empty_or_unreadable"

                    # only count rows for reasonably small tabular files
                    try:
                        if float(st["size_mb"]) <= 500:
                            n_rows = quick_count_rows(p)
                    except Exception:
                        pass

                for d in matched_datasets:
                    role = role_guess(p, cols)
                    rows.append({
                        "dataset": d,
                        "path": str(p),
                        "file_name": p.name,
                        "parent_dir": str(p.parent),
                        "role_guess": role,
                        "priority_score": candidate_priority(p, d, role, n_cols, n_rows),
                        "exists": st["exists"],
                        "readable": st["readable"],
                        "size_mb": st["size_mb"],
                        "mtime": st["mtime"],
                        "n_columns": n_cols,
                        "n_data_rows_minus_header": n_rows,
                        "header_status": header_status,
                        "columns_preview": "|".join(cols[:80]),
                        "source": "filesystem_scan"
                    })

    return rows


def scan_manifest_candidates(manifest):
    rows = []

    for r in manifest:
        p = r.get("path", "")
        if not p:
            continue
        if is_excluded(p):
            continue

        matched_datasets = [d for d in DATASETS if dataset_hit(p, d)]
        if not matched_datasets:
            continue

        plow = p.lower()
        if not (plow.endswith(TEXT_EXT) or plow.endswith(tuple(e.lower() for e in OBJECT_EXT))):
            continue

        cols = []
        detected = r.get("detected_columns", "")
        if detected:
            cols = detected.split("|")

        n_cols = r.get("n_detected_columns", "")
        try:
            n_cols_int = int(n_cols)
        except Exception:
            n_cols_int = len(cols)

        st = safe_stat(p)

        for d in matched_datasets:
            role = role_guess(p, cols)
            rows.append({
                "dataset": d,
                "path": p,
                "file_name": Path(p).name,
                "parent_dir": str(Path(p).parent),
                "role_guess": role,
                "priority_score": candidate_priority(p, d, role, n_cols_int, ""),
                "exists": st["exists"],
                "readable": st["readable"],
                "size_mb": st["size_mb"],
                "mtime": st["mtime"],
                "n_columns": n_cols,
                "n_data_rows_minus_header": "",
                "header_status": r.get("header_read_status", ""),
                "columns_preview": detected,
                "source": "step01_manifest"
            })

    return rows


def choose_shortlist(rows):
    out = []
    df_by_key = {}

    # de-duplicate by dataset + path
    for r in rows:
        key = (r["dataset"], r["path"])
        if key not in df_by_key or int(r["priority_score"]) > int(df_by_key[key]["priority_score"]):
            df_by_key[key] = r

    rows = list(df_by_key.values())

    for dataset in DATASETS:
        sub = [r for r in rows if r["dataset"] == dataset]

        expr = [r for r in sub if r["role_guess"] == "expression_matrix_candidate"]
        meta = [r for r in sub if r["role_guess"] == "metadata_candidate"]
        obj = [r for r in sub if r["role_guess"] == "serialized_object_candidate"]

        expr = sorted(expr, key=lambda x: (-int(x["priority_score"]), -float(x["size_mb"] or 0), x["path"]))[:10]
        meta = sorted(meta, key=lambda x: (-int(x["priority_score"]), -float(x["size_mb"] or 0), x["path"]))[:10]
        obj = sorted(obj, key=lambda x: (-int(x["priority_score"]), -float(x["size_mb"] or 0), x["path"]))[:10]

        for role_group, items in [
            ("top_expression_candidates", expr),
            ("top_metadata_candidates", meta),
            ("top_serialized_object_candidates", obj)
        ]:
            for rank, r in enumerate(items, start=1):
                rr = dict(r)
                rr["shortlist_group"] = role_group
                rr["rank_in_group"] = rank
                out.append(rr)

    return out


def main():
    print(f"[{now()}] Step06A locating external adult validation inputs started")

    manifest = read_manifest(STEP01_MANIFEST)
    rows = []
    rows.extend(scan_manifest_candidates(manifest))
    rows.extend(scan_filesystem_candidates())

    # de-duplicate exact duplicate rows
    uniq = {}
    for r in rows:
        key = (r["dataset"], r["path"], r["role_guess"])
        if key not in uniq or int(r["priority_score"]) > int(uniq[key]["priority_score"]):
            uniq[key] = r
    rows = list(uniq.values())

    rows = sorted(rows, key=lambda x: (x["dataset"], -int(x["priority_score"]), x["role_guess"], x["path"]))

    fields = [
        "dataset", "path", "file_name", "parent_dir",
        "role_guess", "priority_score",
        "exists", "readable", "size_mb", "mtime",
        "n_columns", "n_data_rows_minus_header",
        "header_status", "columns_preview", "source"
    ]

    write_tsv(rows, OUT / "26_step06A_external_adult_input_candidates.tsv", fields)

    shortlist = choose_shortlist(rows)
    short_fields = fields + ["shortlist_group", "rank_in_group"]
    write_tsv(shortlist, OUT / "27_step06A_external_adult_input_shortlist.tsv", short_fields)

    # Dataset-level summary
    summary = []
    for d in DATASETS:
        sub = [r for r in rows if r["dataset"] == d]
        summary.append({
            "dataset": d,
            "n_total_candidates": len(sub),
            "n_expression_candidates": sum(r["role_guess"] == "expression_matrix_candidate" for r in sub),
            "n_metadata_candidates": sum(r["role_guess"] == "metadata_candidate" for r in sub),
            "n_serialized_object_candidates": sum(r["role_guess"] == "serialized_object_candidate" for r in sub),
            "best_expression_path": next((r["path"] for r in sorted([x for x in sub if x["role_guess"] == "expression_matrix_candidate"], key=lambda x: -int(x["priority_score"]))), ""),
            "best_metadata_path": next((r["path"] for r in sorted([x for x in sub if x["role_guess"] == "metadata_candidate"], key=lambda x: -int(x["priority_score"]))), ""),
        })

    write_tsv(
        summary,
        OUT / "28_step06A_external_adult_input_summary.tsv",
        [
            "dataset",
            "n_total_candidates",
            "n_expression_candidates",
            "n_metadata_candidates",
            "n_serialized_object_candidates",
            "best_expression_path",
            "best_metadata_path"
        ]
    )

    with open(OUT / "29_step06A_external_adult_input_locator_summary.md", "w") as f:
        f.write("# NeuroTRACE Step06A external adult input locator summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Purpose\n")
        f.write("Step06A locates candidate expression and metadata inputs for GSE102741 and GSE64018. It intentionally excludes DevMap/DevBridge score tables, deconvolution outputs and NeuroTRACE downstream outputs. These candidates will be manually reviewed before module-score validation.\n\n")
        f.write("## Dataset-level summary\n\n")

        for row in summary:
            f.write(f"### {row['dataset']}\n")
            f.write(f"- Total candidates: {row['n_total_candidates']}\n")
            f.write(f"- Expression candidates: {row['n_expression_candidates']}\n")
            f.write(f"- Metadata candidates: {row['n_metadata_candidates']}\n")
            f.write(f"- Serialized object candidates: {row['n_serialized_object_candidates']}\n")
            f.write(f"- Best expression path: `{row['best_expression_path']}`\n")
            f.write(f"- Best metadata path: `{row['best_metadata_path']}`\n\n")

        f.write("## Main outputs\n")
        f.write("- `26_step06A_external_adult_input_candidates.tsv`\n")
        f.write("- `27_step06A_external_adult_input_shortlist.tsv`\n")
        f.write("- `28_step06A_external_adult_input_summary.tsv`\n\n")
        f.write("## Next step\n")
        f.write("After checking the shortlist, freeze exact expression and metadata paths for GSE102741/GSE64018, then run Step06B to score NTM1/NTM2/NTM3 in external adult cohorts.\n")

    print(f"[{now()}] Step06A done")
    print(f"[{now()}] Candidates: {len(rows)}")
    print(f"[{now()}] Shortlist rows: {len(shortlist)}")
    print(f"[{now()}] Results: {OUT}")


if __name__ == "__main__":
    main()
