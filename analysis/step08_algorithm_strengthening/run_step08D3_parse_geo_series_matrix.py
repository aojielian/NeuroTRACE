#!/usr/bin/env python3

import gzip
import csv
import re
from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
GEO_DIR = BASE / "neurotrace_algorithm_project" / "cross_disease_GEO"
PARSED = GEO_DIR / "parsed"

STEP = BASE / "neurotrace_algorithm_project" / "step08_algorithm_strengthening"
OUT = STEP / "results"

PARSED.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

GSES = ["GSE53987", "GSE12649", "GSE21138"]


def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def open_maybe_gz(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def split_geo_line(line):
    # GEO series matrix uses quoted tab-delimited fields
    return next(csv.reader([line.rstrip("\n\r")], delimiter="\t"))


def clean_field(x):
    if x is None:
        return ""
    x = str(x)
    x = x.strip()
    if len(x) >= 2 and x[0] == '"' and x[-1] == '"':
        x = x[1:-1]
    return x.strip()


def sanitize_colname(x):
    x = clean_field(x)
    x = re.sub(r"[^A-Za-z0-9_.:-]+", "_", x)
    x = re.sub(r"_+", "_", x)
    x = x.strip("_")
    return x if x else "unknown"


def parse_series_matrix(path):
    meta = {}
    table_lines = []
    in_table = False

    with open_maybe_gz(path) as f:
        for line in f:
            line = line.rstrip("\n\r")

            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue

            if in_table:
                table_lines.append(line)
                continue

            if line.startswith("!"):
                fields = split_geo_line(line)
                key = clean_field(fields[0])
                vals = [clean_field(x) for x in fields[1:]]
                meta.setdefault(key, []).append(vals)

    if len(table_lines) == 0:
        raise RuntimeError(f"No series_matrix_table found in {path}")

    header = [clean_field(x) for x in split_geo_line(table_lines[0])]
    rows = []
    for ln in table_lines[1:]:
        vals = [clean_field(x) for x in split_geo_line(ln)]
        if len(vals) < len(header):
            vals += [""] * (len(header) - len(vals))
        elif len(vals) > len(header):
            vals = vals[:len(header)]
        rows.append(vals)

    expr = pd.DataFrame(rows, columns=header)

    return meta, expr


def build_sample_metadata(meta, sample_ids):
    sample_ids = list(sample_ids)
    n = len(sample_ids)
    md = pd.DataFrame({"sample_id": sample_ids})

    sample_keys = [k for k in meta.keys() if k.startswith("!Sample_")]

    for key in sample_keys:
        values_lists = meta[key]
        # Usually one row per key with n values
        for idx, vals in enumerate(values_lists):
            if len(vals) != n:
                continue

            base_col = sanitize_colname(key.replace("!Sample_", ""))
            col = base_col if idx == 0 else f"{base_col}_{idx+1}"
            md[col] = vals

    # Parse characteristics columns into key-value fields
    char_cols = [c for c in md.columns if c.lower().startswith("characteristics")]
    parsed = {}

    for c in char_cols:
        for i, val in enumerate(md[c].astype(str).tolist()):
            # Typical format: "diagnosis: Control"
            if ":" in val:
                k, v = val.split(":", 1)
                k = sanitize_colname(k.lower())
                v = v.strip()
                parsed.setdefault(k, [""] * n)
                if not parsed[k][i]:
                    parsed[k][i] = v

    for k, vals in parsed.items():
        # avoid overwriting existing sample columns
        col = k
        if col in md.columns:
            col = f"characteristics_{k}"
        md[col] = vals

    return md


def infer_dx_columns(md):
    candidates = []
    for c in md.columns:
        low = c.lower()
        if any(x in low for x in ["diagnosis", "disease", "condition", "group", "phenotype", "status"]):
            candidates.append(c)

    value_hits = []
    for c in md.columns:
        vals = md[c].astype(str).str.lower().unique().tolist()
        joined = " | ".join(vals[:200])
        if any(x in joined for x in ["schiz", "bipolar", "depression", "mdd", "control", "healthy", "autism", "asd"]):
            value_hits.append(c)

    return sorted(set(candidates + value_hits))


def infer_tissue_columns(md):
    candidates = []
    for c in md.columns:
        low = c.lower()
        if any(x in low for x in ["tissue", "brain", "region", "area", "cortex", "source_name"]):
            candidates.append(c)

    return sorted(set(candidates))


def expression_audit(expr):
    first_col = expr.columns[0]
    sample_cols = list(expr.columns[1:])
    n_rows = expr.shape[0]
    n_samples = len(sample_cols)

    ids = expr[first_col].astype(str).head(200).tolist()

    n_gene_symbol_like = sum(bool(re.match(r"^[A-Z][A-Z0-9.-]{1,15}$", x)) for x in ids)
    n_ensembl_like = sum(x.startswith("ENSG") for x in ids)
    n_affy_like = sum(bool(re.match(r"^[0-9A-Za-z]+(_[a-z]+)?_at$", x)) or "_s_at" in x or "_x_at" in x for x in ids)
    n_illumina_like = sum(x.startswith("ILMN_") for x in ids)

    return {
        "feature_id_col": first_col,
        "n_features": n_rows,
        "n_samples": n_samples,
        "n_gene_symbol_like_first200": n_gene_symbol_like,
        "n_ensembl_like_first200": n_ensembl_like,
        "n_affy_like_first200": n_affy_like,
        "n_illumina_like_first200": n_illumina_like,
        "first_10_feature_ids": ";".join(ids[:10])
    }


def normalize_expression_table(expr):
    first_col = expr.columns[0]
    expr = expr.rename(columns={first_col: "feature_id"})
    sample_cols = [c for c in expr.columns if c != "feature_id"]

    out = expr.copy()
    out["feature_id"] = out["feature_id"].astype(str)

    for c in sample_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce")

    # Drop empty feature rows
    out = out[out["feature_id"].astype(str).str.len() > 0].copy()
    out = out.drop_duplicates("feature_id")

    return out


def write_tsv(df, path, compression=None):
    df.to_csv(path, sep="\t", index=False, compression=compression)


def main():
    print(f"[{now()}] Step08D3 parse GEO series matrix started")

    manifest_rows = []
    sample_summary_rows = []
    dx_audit_rows = []
    tissue_audit_rows = []
    platform_rows = []

    for gse in GSES:
        in_path = GEO_DIR / f"{gse}_series_matrix.txt.gz"

        if not in_path.exists():
            # also support nested matrix path, just in case
            nested = GEO_DIR / gse / "matrix" / f"{gse}_series_matrix.txt.gz"
            if nested.exists():
                in_path = nested

        if not in_path.exists():
            manifest_rows.append({
                "dataset": gse,
                "status": "missing",
                "input_path": str(in_path),
                "message": "series matrix not found"
            })
            continue

        print(f"[{now()}] Parsing {gse}: {in_path}")

        meta, expr_raw = parse_series_matrix(in_path)
        expr = normalize_expression_table(expr_raw)

        sample_cols = [c for c in expr.columns if c != "feature_id"]
        sample_md = build_sample_metadata(meta, sample_cols)

        dataset_dir = PARSED / gse
        dataset_dir.mkdir(parents=True, exist_ok=True)

        expr_path = dataset_dir / f"{gse}.expression_probe.tsv.gz"
        sample_path = dataset_dir / f"{gse}.samples.tsv"

        write_tsv(expr, expr_path, compression="gzip")
        write_tsv(sample_md, sample_path)

        # series/platform metadata
        platforms = []
        for k in ["!Series_platform_id", "!Series_platform_taxid"]:
            if k in meta:
                platforms.append(f"{k}=" + "|".join([";".join(x) for x in meta[k]]))
        platform_info = "; ".join(platforms)

        audit = expression_audit(expr)
        platform_rows.append({
            "dataset": gse,
            "input_path": str(in_path),
            "expression_path": str(expr_path),
            "sample_path": str(sample_path),
            "platform_info": platform_info,
            **audit
        })

        dx_cols = infer_dx_columns(sample_md)
        tissue_cols = infer_tissue_columns(sample_md)

        # Diagnosis audit
        if dx_cols:
            for c in dx_cols:
                vc = sample_md[c].astype(str).value_counts(dropna=False).reset_index()
                vc.columns = ["value", "n"]
                for _, r in vc.head(100).iterrows():
                    dx_audit_rows.append({
                        "dataset": gse,
                        "candidate_dx_col": c,
                        "value": r["value"],
                        "n": int(r["n"])
                    })
        else:
            dx_audit_rows.append({
                "dataset": gse,
                "candidate_dx_col": "",
                "value": "NO_DX_COLUMN_DETECTED",
                "n": 0
            })

        # Tissue audit
        if tissue_cols:
            for c in tissue_cols:
                vc = sample_md[c].astype(str).value_counts(dropna=False).reset_index()
                vc.columns = ["value", "n"]
                for _, r in vc.head(100).iterrows():
                    tissue_audit_rows.append({
                        "dataset": gse,
                        "candidate_tissue_col": c,
                        "value": r["value"],
                        "n": int(r["n"])
                    })
        else:
            tissue_audit_rows.append({
                "dataset": gse,
                "candidate_tissue_col": "",
                "value": "NO_TISSUE_COLUMN_DETECTED",
                "n": 0
            })

        sample_summary_rows.append({
            "dataset": gse,
            "n_samples": sample_md.shape[0],
            "n_sample_metadata_cols": sample_md.shape[1],
            "candidate_dx_cols": ";".join(dx_cols),
            "candidate_tissue_cols": ";".join(tissue_cols)
        })

        manifest_rows.append({
            "dataset": gse,
            "status": "parsed",
            "input_path": str(in_path),
            "expression_path": str(expr_path),
            "sample_path": str(sample_path),
            "n_features": expr.shape[0],
            "n_samples": len(sample_cols),
            "n_sample_metadata_cols": sample_md.shape[1]
        })

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(OUT / "69_step08D3_geo_series_parse_manifest.tsv", sep="\t", index=False)

    pd.DataFrame(platform_rows).to_csv(OUT / "70_step08D3_expression_platform_audit.tsv", sep="\t", index=False)
    pd.DataFrame(sample_summary_rows).to_csv(OUT / "71_step08D3_sample_metadata_summary.tsv", sep="\t", index=False)
    pd.DataFrame(dx_audit_rows).to_csv(OUT / "72_step08D3_diagnosis_column_audit.tsv", sep="\t", index=False)
    pd.DataFrame(tissue_audit_rows).to_csv(OUT / "73_step08D3_tissue_column_audit.tsv", sep="\t", index=False)

    with open(OUT / "74_step08D3_geo_series_parse_summary.md", "w") as f:
        f.write("# NeuroTRACE Step08D3 GEO series matrix parse summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Purpose\n")
        f.write("Parse locally uploaded GEO series matrix files for cross-disease adult cortex validation.\n\n")
        f.write("## Parse manifest\n\n")
        f.write(manifest.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Platform / expression audit\n\n")
        if platform_rows:
            f.write(pd.DataFrame(platform_rows).to_markdown(index=False))
        else:
            f.write("No platform rows.\n")
        f.write("\n\n")
        f.write("## Sample metadata summary\n\n")
        if sample_summary_rows:
            f.write(pd.DataFrame(sample_summary_rows).to_markdown(index=False))
        else:
            f.write("No sample metadata rows.\n")
        f.write("\n\n")
        f.write("## Interpretation\n")
        f.write("If feature IDs are probe IDs rather than gene symbols, the next step must map probes to gene symbols using platform annotation before NTM scoring.\n")

    print(f"[{now()}] Step08D3 done")
    print(f"[{now()}] Results: {OUT}")
    print(f"[{now()}] Parsed data: {PARSED}")


if __name__ == "__main__":
    main()
