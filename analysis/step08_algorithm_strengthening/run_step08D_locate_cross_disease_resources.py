#!/usr/bin/env python3

import os
import csv
import gzip
from pathlib import Path
from datetime import datetime

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
PROJECT = BASE / "neurotrace_algorithm_project"
STEP = PROJECT / "step08_algorithm_strengthening"
OUT = STEP / "results"
OUT.mkdir(parents=True, exist_ok=True)

SCAN_ROOTS = [
    Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD"),
    Path("/gpfs/hpc/home/lijc/lianaoj/autism_scRNA"),
    Path("/gpfs/hpc/home/lijc/lianaoj/AAAY_data"),
    Path("/gpfs/hpc/home/lijc/lianaoj"),
]

DATASET_TERMS = [
    "SCZ", "schiz", "schizophrenia",
    "BD", "bipolar",
    "MDD", "depression",
    "psychencode", "commonmind", "CMC",
    "Gandal", "BrainGVEX",
    "ROSMAP", "AD", "Alzheimer",
    "epilepsy"
]

FILE_EXT = (
    ".tsv", ".tsv.gz", ".csv", ".csv.gz",
    ".txt", ".txt.gz",
    ".rds", ".RDS", ".RData", ".rdata",
    ".mtx", ".mtx.gz", ".h5ad", ".h5seurat", ".h5"
)

EXCLUDE_DIRS = [
    "/tmp/", "/logs/", "/outs/", "/errs/", "/figures/",
    "/.git/", "/conda/", "/mambaforge/",
    "/neurotrace_algorithm_project/step02",
    "/neurotrace_algorithm_project/step03",
    "/neurotrace_algorithm_project/step04",
    "/neurotrace_algorithm_project/step05",
    "/neurotrace_algorithm_project/step06",
    "/neurotrace_algorithm_project/step07",
]

EXPR_TERMS = [
    "expr", "expression", "count", "counts", "matrix",
    "rsem", "fpkm", "tpm", "cpm", "normalized", "voom"
]

META_TERMS = [
    "meta", "metadata", "pheno", "phenotype",
    "clinical", "sample", "subject", "diagnosis", "covariate"
]


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def safe_stat(path):
    try:
        st = path.stat()
        return {
            "size_mb": round(st.st_size / 1024 / 1024, 3),
            "mtime": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
            "readable": os.access(path, os.R_OK)
        }
    except Exception:
        return {"size_mb": "", "mtime": "", "readable": False}


def is_excluded(path):
    low = str(path).lower()
    return any(x.lower() in low for x in EXCLUDE_DIRS)


def term_hits(path):
    low = str(path).lower()
    hits = []
    for t in DATASET_TERMS:
        if t.lower() in low:
            hits.append(t)
    return sorted(set(hits))


def read_header(path):
    try:
        p = str(path)
        if p.endswith(".gz"):
            with gzip.open(p, "rt", errors="replace") as f:
                return f.readline().rstrip("\n\r")
        else:
            with open(p, "rt", errors="replace") as f:
                return f.readline().rstrip("\n\r")
    except Exception:
        return ""


def split_header(h):
    if "\t" in h:
        return h.split("\t")
    if "," in h:
        return h.split(",")
    return h.split()


def role_guess(path, cols):
    low = str(path).lower()
    cols_low = [c.lower() for c in cols]
    ncol = len(cols)

    expr_score = sum(t in low for t in EXPR_TERMS)
    meta_score = sum(t in low for t in META_TERMS)
    dx_col_score = sum(any(x in c for x in ["diagnosis", "dx", "condition", "disease", "group", "case", "control", "scz", "bd"]) for c in cols_low)

    if ncol >= 20 and expr_score > 0:
        return "expression_candidate"
    if dx_col_score >= 2 or meta_score > 0:
        return "metadata_candidate"
    if str(path).lower().endswith((".rds", ".rdata")):
        return "serialized_object_candidate"
    return "uncertain"


def priority(path, role, hits, ncol, size_mb):
    low = str(path).lower()
    score = 0
    score += 20 * len(hits)
    if role == "expression_candidate":
        score += 30
    elif role == "metadata_candidate":
        score += 25
    elif role == "serialized_object_candidate":
        score += 20

    if any(x in low for x in ["commonmind", "psychencode", "gandal", "brainseq", "scz", "schizophrenia", "bipolar"]):
        score += 20

    if any(x in low for x in ["test", "backup", "old", "tmp"]):
        score -= 10

    try:
        if float(size_mb) > 1:
            score += 5
        if float(size_mb) > 50:
            score += 10
    except Exception:
        pass

    try:
        if int(ncol) >= 20:
            score += 5
    except Exception:
        pass

    return score


def write_tsv(rows, path, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    print(f"[{now()}] Step08D locate cross-disease resources started")

    rows = []

    for root in SCAN_ROOTS:
        if not root.exists():
            continue

        print(f"[{now()}] Scanning: {root}")

        for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
            dp = Path(dirpath)
            if is_excluded(dp):
                continue

            for fn in filenames:
                p = dp / fn
                low = str(p).lower()

                if not low.endswith(tuple(e.lower() for e in FILE_EXT)):
                    continue
                if is_excluded(p):
                    continue

                hits = term_hits(p)
                if not hits:
                    continue

                st = safe_stat(p)
                header = ""
                cols = []
                if low.endswith((".tsv", ".tsv.gz", ".csv", ".csv.gz", ".txt", ".txt.gz")):
                    header = read_header(p)
                    cols = split_header(header)

                role = role_guess(p, cols)
                pr = priority(p, role, hits, len(cols), st["size_mb"])

                rows.append({
                    "path": str(p),
                    "file_name": p.name,
                    "parent_dir": str(p.parent),
                    "term_hits": ";".join(hits),
                    "role_guess": role,
                    "priority_score": pr,
                    "size_mb": st["size_mb"],
                    "mtime": st["mtime"],
                    "readable": st["readable"],
                    "n_columns": len(cols) if cols else "",
                    "columns_preview": "|".join(cols[:80]) if cols else "",
                })

    # de-duplicate
    dedup = {}
    for r in rows:
        key = r["path"]
        if key not in dedup or r["priority_score"] > dedup[key]["priority_score"]:
            dedup[key] = r
    rows = list(dedup.values())
    rows.sort(key=lambda x: (-int(x["priority_score"]), str(x["path"])))

    fields = [
        "path", "file_name", "parent_dir", "term_hits",
        "role_guess", "priority_score", "size_mb", "mtime",
        "readable", "n_columns", "columns_preview"
    ]
    write_tsv(rows, OUT / "64_step08D_cross_disease_resource_candidates.tsv", fields)

    shortlist = rows[:100]
    write_tsv(shortlist, OUT / "65_step08D_cross_disease_resource_shortlist_top100.tsv", fields)

    summary = {}
    for r in rows:
        for h in r["term_hits"].split(";"):
            if not h:
                continue
            summary.setdefault(h, {"term": h, "n_files": 0, "n_expression": 0, "n_metadata": 0, "n_serialized": 0})
            summary[h]["n_files"] += 1
            if r["role_guess"] == "expression_candidate":
                summary[h]["n_expression"] += 1
            if r["role_guess"] == "metadata_candidate":
                summary[h]["n_metadata"] += 1
            if r["role_guess"] == "serialized_object_candidate":
                summary[h]["n_serialized"] += 1

    summary_rows = list(summary.values())
    summary_rows.sort(key=lambda x: (-x["n_files"], x["term"]))
    write_tsv(summary_rows, OUT / "66_step08D_cross_disease_resource_summary.tsv",
              ["term", "n_files", "n_expression", "n_metadata", "n_serialized"])

    with open(OUT / "67_step08D_cross_disease_resource_locator_summary.md", "w") as f:
        f.write("# NeuroTRACE Step08D cross-disease resource locator summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Purpose\n")
        f.write("Locate local candidate resources for cross-disease generalization/specificity testing, especially SCZ/BD/MDD adult cortical expression resources.\n\n")
        f.write("## Counts\n")
        f.write(f"- Total candidate files: {len(rows)}\n")
        f.write(f"- Shortlist rows: {len(shortlist)}\n\n")
        f.write("## Term summary\n\n")
        for r in summary_rows:
            f.write(f"- {r['term']}: files={r['n_files']}, expression={r['n_expression']}, metadata={r['n_metadata']}, serialized={r['n_serialized']}\n")
        f.write("\n## Main outputs\n")
        f.write("- `64_step08D_cross_disease_resource_candidates.tsv`\n")
        f.write("- `65_step08D_cross_disease_resource_shortlist_top100.tsv`\n")
        f.write("- `66_step08D_cross_disease_resource_summary.tsv`\n")

    print(f"[{now()}] Step08D done")
    print(f"[{now()}] Candidate files: {len(rows)}")
    print(f"[{now()}] Results: {OUT}")


if __name__ == "__main__":
    main()
