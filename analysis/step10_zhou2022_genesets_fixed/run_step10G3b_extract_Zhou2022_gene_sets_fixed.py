#!/usr/bin/env python3

import re
import math
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd

BASE = Path("/gpfs/hpc/home/lijc/lianaoj/DevMap_ASD")
STEP = BASE / "neurotrace_algorithm_project" / "step10_zhou2022_genesets_fixed"
OUT = STEP / "results"
GS = STEP / "gene_sets"

OUT.mkdir(parents=True, exist_ok=True)
GS.mkdir(parents=True, exist_ok=True)

XLSX = BASE / "41588_2022_1148_MOESM4_ESM.xlsx"

HEADER_ROWS = {
    "Table S2": 1,
    "Table S3": 2,
    "Table S5": 2,
    "Table S6": 2,
    "Table S7": 2,
    "Table S9": 2,
    "Table S13": 1,
    "Table S17": 1,
}

def now():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def clean_gene(x):
    if x is None or pd.isna(x):
        return ""
    g = str(x).strip().upper()
    g = g.replace('"', '').replace("'", "")
    g = g.split(";")[0].split(",")[0].split("|")[0].strip()
    if g in {"", "NA", "NAN", "NULL", "---", ".", "GENE", "GENES", "HGNC", "SYMBOL", "GENE_SYMBOL"}:
        return ""
    if g.startswith("ENSG"):
        return ""
    if re.match(r"^[A-Z0-9][A-Z0-9._-]{1,30}$", g):
        return g
    return ""

def norm_cols(df):
    new = []
    seen = {}
    for c in df.columns:
        x = str(c).strip()
        x = re.sub(r"\s+", "_", x)
        x = re.sub(r"[^\w.\-/:>=<]+", "_", x)
        x = x.strip("_")
        if not x:
            x = "unnamed"
        if x in seen:
            seen[x] += 1
            x = f"{x}.{seen[x]}"
        else:
            seen[x] = 0
        new.append(x)
    df.columns = new
    return df

def read_sheet(sheet):
    df = pd.read_excel(XLSX, sheet_name=sheet, header=HEADER_ROWS[sheet], engine="openpyxl")
    df = norm_cols(df)
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    return df

def genes_from(df, gene_col):
    if gene_col not in df.columns:
        return []
    genes = [clean_gene(x) for x in df[gene_col].tolist()]
    return sorted(set([g for g in genes if g]))

def num(df, col):
    if col is None or col not in df.columns:
        return pd.Series([np.nan] * len(df))
    return pd.to_numeric(df[col], errors="coerce")

def flag_x(df, col):
    if col is None or col not in df.columns:
        return pd.Series([False] * len(df))
    s = df[col].astype(str).str.lower().str.strip()
    return s.isin(["x", "1", "true", "yes", "y"])

def contains(df, col, pats):
    if col is None or col not in df.columns:
        return pd.Series([False] * len(df))
    s = df[col].astype(str).str.lower()
    m = pd.Series([False] * len(df))
    for p in pats:
        m = m | s.str.contains(p.lower(), regex=False, na=False)
    return m

def chi2_sf_1df(x):
    # survival function for chi-square df=1 = erfc(sqrt(x/2))
    x = pd.to_numeric(x, errors="coerce")
    return x.apply(lambda v: math.erfc(math.sqrt(v / 2.0)) if pd.notna(v) and v >= 0 else np.nan)

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    idxs = np.where(ok)[0]
    if len(idxs) == 0:
        return out
    pv = p[idxs]
    order = np.argsort(pv)
    tmp = np.empty(len(pv))
    prev = 1.0
    n = len(pv)
    for rank in range(n, 0, -1):
        ii = order[rank - 1]
        val = pv[ii] * n / rank
        prev = min(prev, val)
        tmp[ii] = min(prev, 1.0)
    out[idxs] = tmp
    return out

records = []

def write_set(name, genes, desc, source, rule):
    genes = sorted(set([g for g in genes if clean_gene(g)]))
    path = GS / f"{name}.txt"
    with open(path, "w") as f:
        for g in genes:
            f.write(g + "\n")
    records.append({
        "geneset": name,
        "n_genes": len(genes),
        "path": str(path),
        "description": desc,
        "source_table": source,
        "filter_rule": rule,
        "status": "ok" if genes else "empty",
    })
    return genes

def write_table(df, name):
    df.to_csv(OUT / name, sep="\t", index=False)

def main():
    print(f"[{now()}] Step10G3b fixed extraction started")
    if not XLSX.exists():
        raise FileNotFoundError(XLSX)

    sheet_audit = []

    # ---------- Table S2 ----------
    s2 = read_sheet("Table S2")
    sheet_audit.append({"sheet": "Table S2", "n_rows": len(s2), "n_cols": s2.shape[1], "columns": "|".join(s2.columns)})
    write_set("Zhou_Known_ASD_NDD_618", genes_from(s2, "HGNC"),
              "Known dominant or X-linked ASD/NDD genes.", "Table S2", "all HGNC")
    write_set("Zhou_Known_ASD_NDD_high_pLI", genes_from(s2[num(s2, "ExACpLI") >= 0.9], "HGNC"),
              "Known ASD/NDD genes with ExAC pLI >= 0.9.", "Table S2", "ExACpLI >= 0.9")
    write_set("Zhou_Known_ASD_NDD_LOEUF_top", genes_from(s2[num(s2, "LOEUF_decile") <= 2], "HGNC"),
              "Known ASD/NDD genes with LOEUF decile <= 2.", "Table S2", "LOEUF_decile <= 2")
    write_set("Zhou_Known_ASD_NDD_GeneScore_1_or_1S",
              genes_from(s2[s2["GeneScore"].astype(str).str.upper().isin(["1", "1.0", "1S"])], "HGNC"),
              "Known ASD/NDD genes with GeneScore 1 or 1S.", "Table S2", "GeneScore in 1/1S")

    # ---------- Table S3 ----------
    s3 = read_sheet("Table S3")
    sheet_audit.append({"sheet": "Table S3", "n_rows": len(s3), "n_cols": s3.shape[1], "columns": "|".join(s3.columns)})
    s3_p = "The_final_DenovoWEST_p-value_as_the_minimum_of_pAllEnrich_and_pMisComb"
    s3_final_p = num(s3, s3_p)
    s3_fdr = bh_fdr(s3_final_p.values)
    s3["computed_BH_FDR_from_final_DenovoWEST_p"] = s3_fdr
    write_table(s3, "01_step10G3b_TableS3_with_computed_FDR.tsv")

    write_set("Zhou_DeNovoWEST_stage1_p_lt_0.001", genes_from(s3[s3_final_p < 0.001], "Gene_symbol"),
              "Stage 1 DeNovoWEST genes with final p < 0.001.", "Table S3", f"{s3_p} < 0.001")
    write_set("Zhou_DeNovoWEST_stage1_FDR_lt_0.1", genes_from(s3[s3["computed_BH_FDR_from_final_DenovoWEST_p"] <= 0.1], "Gene_symbol"),
              "Stage 1 DeNovoWEST genes with BH-FDR <= 0.1 computed from final p.", "Table S3", "computed BH-FDR <= 0.1")
    write_set("Zhou_DeNovoWEST_stage1_exomewide", genes_from(s3[s3_final_p < 2.5e-6], "Gene_symbol"),
              "Stage 1 DeNovoWEST exome-wide proxy final p < 2.5e-6.", "Table S3", f"{s3_p} < 2.5e-6")
    lof = num(s3, "AutismMerged_LoF")
    dmis = num(s3, "Total_number_of_de_novo_Dmis_variants_REVEL>=0.5_in_discovery_cohort")
    write_set("Zhou_DeNovoWEST_stage1_LoF_or_Dmis_driven", genes_from(s3[(lof > 0) | (dmis > 0)], "Gene_symbol"),
              "Stage 1 genes with at least one LoF or damaging missense de novo event.", "Table S3", "AutismMerged_LoF > 0 OR Dmis_REVEL>=0.5 > 0")

    # ---------- Table S5 ----------
    s5 = read_sheet("Table S5")
    sheet_audit.append({"sheet": "Table S5", "n_rows": len(s5), "n_cols": s5.shape[1], "columns": "|".join(s5.columns)})
    write_set("Zhou_constrained_background_5754", genes_from(s5, "HGNC"),
              "Constrained background genes used in Zhou gene-set membership table.", "Table S5", "all HGNC")
    write_set("Zhou_ASD_deNovo_enrich_P001", genes_from(s5[flag_x(s5, "ASD_de_novo_enrich_P<0.01")], "HGNC"),
              "ASD de novo enrichment P<0.01 membership from Table S5.", "Table S5", "ASD_de_novo_enrich_P<0.01 == x")
    write_set("Zhou_ASD_deNovo_enrich_P005", genes_from(s5[flag_x(s5, "ASD_de_novo_enrich_P<0.05")], "HGNC"),
              "ASD de novo enrichment P<0.05 membership from Table S5.", "Table S5", "ASD_de_novo_enrich_P<0.05 == x")
    ddg_mask = flag_x(s5, "DDG2P_de_novo_enrich_P<0.01") | flag_x(s5, "DDG2P_de_novo_enrich_P<0.05")
    write_set("Zhou_DDG2P_deNovo_enrich", genes_from(s5[ddg_mask], "HGNC"),
              "DDG2P de novo enrichment membership from Table S5.", "Table S5", "DDG2P P<0.01 or P<0.05 == x")
    write_set("Zhou_SCZ_case_control", genes_from(s5[flag_x(s5, "SCZ_case-control_P<0.05")], "HGNC"),
              "SCZ case-control P<0.05 membership from Table S5.", "Table S5", "SCZ_case-control_P<0.05 == x")
    loeuf_mask = flag_x(s5, "gnomAD_LOEUF_0~10%") | flag_x(s5, "gnomAD_LOEUF_10~20%")
    write_set("Zhou_LOEUF_top", genes_from(s5[loeuf_mask], "HGNC"),
              "gnomAD LOEUF top 20 percent membership from Table S5.", "Table S5", "LOEUF 0-10 or 10-20 == x")
    exac_cols = [c for c in s5.columns if "ExAC" in c or "pLI" in c or "sHet" in c]
    exac_mask = pd.Series([False] * len(s5))
    for c in exac_cols:
        exac_mask = exac_mask | flag_x(s5, c)
    write_set("Zhou_ExACpLI", genes_from(s5[exac_mask], "HGNC"),
              "ExAC pLI/sHet constrained membership from Table S5.", "Table S5", "any ExAC/sHet/pLI membership == x")

    # ---------- Table S6 ----------
    s6 = read_sheet("Table S6")
    sheet_audit.append({"sheet": "Table S6", "n_rows": len(s6), "n_cols": s6.shape[1], "columns": "|".join(s6.columns)})
    prioritized = flag_x(s6, "Prioritized")
    write_set("Zhou_TDT_prioritized_260", genes_from(s6[prioritized], "Gene_symbol"),
              "Prioritized rare inherited LoF TDT genes from Table S6.", "Table S6", "Prioritized == x")
    cyto = s6["CytoBand"].astype(str).str.upper()
    autosomal = prioritized & ~cyto.str.startswith("X")
    write_set("Zhou_TDT_prioritized_autosomal", genes_from(s6[autosomal], "Gene_symbol"),
              "Autosomal prioritized TDT genes from Table S6.", "Table S6", "Prioritized == x and CytoBand not X")

    tdt_cols = [c for c in s6.columns if c.endswith("TDTStat")]
    pmat = pd.DataFrame({c: chi2_sf_1df(s6[c]) for c in tdt_cols})
    s6["min_TDT_p_from_TDTStat"] = pmat.min(axis=1)
    s6["computed_TDT_BH_FDR"] = bh_fdr(s6["min_TDT_p_from_TDTStat"].values)
    write_table(s6, "02_step10G3b_TableS6_with_TDT_p_FDR.tsv")
    write_set("Zhou_TDT_nominal_p_lt_0.05", genes_from(s6[s6["min_TDT_p_from_TDTStat"] < 0.05], "Gene_symbol"),
              "TDT genes with chi-square-derived min p < 0.05.", "Table S6", "min chi-square p from TDTStat columns < 0.05")
    write_set("Zhou_TDT_FDR_lt_0.1", genes_from(s6[s6["computed_TDT_BH_FDR"] <= 0.1], "Gene_symbol"),
              "TDT genes with BH-FDR <= 0.1 from chi-square-derived p.", "Table S6", "computed_TDT_BH_FDR <= 0.1")

    trans_mask = pd.Series([False] * len(s6))
    for prefix in ["AllLoF,pExt>=0.9", "AllLoF,pExt>=0.1", "AllLoF", "URLoF,pExt>=0.9", "URLoF,pExt>=0.1", "URLoF"]:
        tc = f"{prefix}Trans"
        nc = f"{prefix}NonTrans"
        if tc in s6.columns and nc in s6.columns:
            trans_mask = trans_mask | (num(s6, tc) > num(s6, nc))
    write_set("Zhou_TDT_overtransmitted_HC_LoF", genes_from(s6[prioritized & trans_mask], "Gene_symbol"),
              "Prioritized genes with over-transmitted LoF evidence.", "Table S6", "Prioritized == x and Trans > NonTrans in LoF filters")

    # ---------- Table S7 ----------
    s7 = read_sheet("Table S7")
    sheet_audit.append({"sheet": "Table S7", "n_rows": len(s7), "n_cols": s7.shape[1], "columns": "|".join(s7.columns)})
    write_set("Zhou_Selected_404_combined", genes_from(s7, "HGNC"),
              "404 selected genes in combined Stage 1+2 de novo analysis.", "Table S7", "all HGNC")
    # Table title states 159 de novo + 245 inherited/TDT; no explicit source column, so split by row order.
    write_set("Zhou_Selected_159_deNovo", genes_from(s7.iloc[:159], "HGNC"),
              "First 159 selected genes, corresponding to top de novo-enriched set described in Table S7 title.", "Table S7", "first 159 rows")
    write_set("Zhou_Selected_245_inherited_TDT", genes_from(s7.iloc[159:], "HGNC"),
              "Remaining 245 selected genes, corresponding to inherited TDT-prioritized set described in Table S7 title.", "Table S7", "rows 160 onward")

    p7 = num(s7, "pDenovoWEST_Meta")
    s7["computed_BH_FDR_from_pDenovoWEST_Meta"] = bh_fdr(p7.values)
    write_table(s7, "03_step10G3b_TableS7_with_computed_FDR.tsv")
    write_set("Zhou_Combined_DeNovoWEST_FDR_lt_0.1", genes_from(s7[s7["computed_BH_FDR_from_pDenovoWEST_Meta"] <= 0.1], "HGNC"),
              "Combined Stage 1+2 DeNovoWEST genes with BH-FDR <= 0.1.", "Table S7", "computed BH-FDR <= 0.1 from pDenovoWEST_Meta")
    write_set("Zhou_Combined_DeNovoWEST_exomewide", genes_from(s7[p7 < 2.5e-6], "HGNC"),
              "Combined Stage 1+2 DeNovoWEST exome-wide proxy p < 2.5e-6.", "Table S7", "pDenovoWEST_Meta < 2.5e-6")

    # ---------- Table S9 ----------
    s9 = read_sheet("Table S9")
    sheet_audit.append({"sheet": "Table S9", "n_rows": len(s9), "n_cols": s9.shape[1], "columns": "|".join(s9.columns)})
    write_set("Zhou_Meta_391_all_selected", genes_from(s9, "HGNC"),
              "391 selected autosomal genes in combined meta-analysis.", "Table S9", "all HGNC")
    p9 = num(s9, "Maximum_combined_p-value_of_DNVs_TDT_case-control_and_VS_TopMed")
    if p9.isna().all():
        p9 = num(s9, "Maximum_combined_p-value_of_DNVs_TDT_case-control_VS_gnomADexomeNonNeuro_and_VS_TopMed")
    if p9.isna().all():
        p9 = num(s9, "Combined_p-value_from_TDT_case-control_and_DNVs")
    s9["computed_Meta_BH_FDR"] = bh_fdr(p9.values)
    write_table(s9, "04_step10G3b_TableS9_with_computed_FDR.tsv")

    write_set("Zhou_Meta_significant_FDR_lt_0.1", genes_from(s9[s9["computed_Meta_BH_FDR"] <= 0.1], "HGNC"),
              "Meta-analysis genes with BH-FDR <= 0.1.", "Table S9", "computed Meta BH-FDR <= 0.1")
    write_set("Zhou_Meta_exomewide", genes_from(s9[contains(s9, "Study-wide_significance_based_on_5_754_constraint_genes_p<8.69E-06", ["yes"])], "HGNC"),
              "Study-wide significant meta-analysis genes.", "Table S9", "Study-wide significance == Yes")
    write_set("Zhou_Meta_known", genes_from(s9[flag_x(s9, "Known")], "HGNC"),
              "Known meta-analysis genes.", "Table S9", "Known == x")
    write_set("Zhou_Meta_novel", genes_from(s9[~flag_x(s9, "Known")], "HGNC"),
              "Novel meta-analysis genes.", "Table S9", "Known != x")
    write_set("Zhou_Meta_inherited_plus_deNovo", genes_from(s9[s9["Ascertainment"].astype(str).str.lower().isin(["denovo", "tdt"])], "HGNC"),
              "Meta-analysis genes ascertained by de novo or TDT evidence.", "Table S9", "Ascertainment in Denovo/TDT")

    # ---------- Table S17 ----------
    s17 = read_sheet("Table S17")
    sheet_audit.append({"sheet": "Table S17", "n_rows": len(s17), "n_cols": s17.shape[1], "columns": "|".join(s17.columns)})
    write_set("Zhou_dnLoF_enriched_96_constrained", genes_from(s17, "HGNC"),
              "96 de novo LoF enriched constrained genes from Table S17.", "Table S17", "all HGNC")
    write_set("Zhou_dnLoF_enriched_ASDdnSignif", genes_from(s17[flag_x(s17, "ASDdnSignif")], "HGNC"),
              "Table S17 genes with ASDdnSignif flag.", "Table S17", "ASDdnSignif == x")
    sfari_mask = s17["SFARICategory"].notna() & ~s17["SFARICategory"].astype(str).isin([".", "NA", "nan"])
    write_set("Zhou_dnLoF_enriched_SFARI_category", genes_from(s17[sfari_mask], "HGNC"),
              "Table S17 genes with SFARI category annotation.", "Table S17", "SFARICategory non-empty")

    # ---------- Table S13 ----------
    s13 = read_sheet("Table S13")
    sheet_audit.append({"sheet": "Table S13", "n_rows": len(s13), "n_cols": s13.shape[1], "columns": "|".join(s13.columns)})
    arch_cols = {
        "A1_neurotransmission": "A1:_neurotransmission",
        "A2_chromatin_modification": "A2:_chromatin_modification",
        "A3_RNA_processing": "A3:_RNA_processing",
        "A4_vesicle_mediated_transport": "A4:_vesicle_mediated_transport",
        "A5_MAPK_signaling_migration": "A5:_MAPK_signaling_and_migration",
        "A6_cytoskeleton_mitosis": "A6:_cytoskeleton_and_mitosis",
    }
    arch_long = []
    for name, col in arch_cols.items():
        vals = num(s13, col)
        if vals.notna().sum() == 0:
            write_set(f"Zhou_Archetype_{name}", [], f"Archetype {name}; column not detected.", "Table S13", f"missing {col}")
            continue
        thr = vals.quantile(0.9)
        write_set(f"Zhou_Archetype_{name}", genes_from(s13[vals >= thr], "gene"),
                  f"Top-decile genes for archetype {name}.", "Table S13", f"{col} >= 90th percentile {thr}")
        tmp = pd.DataFrame({
            "gene": s13["gene"].map(clean_gene),
            "archetype": name,
            "score": vals
        })
        arch_long.append(tmp.dropna())
    if arch_long:
        pd.concat(arch_long, ignore_index=True).to_csv(OUT / "04_step10G3b_Zhou_archetype_gene_scores_long.tsv", sep="\t", index=False)

    # ---------- write audit ----------
    inv = pd.DataFrame(records).sort_values(["status", "geneset"])
    inv.to_csv(OUT / "05_step10G3b_Zhou2022_gene_set_inventory_fixed.tsv", sep="\t", index=False)
    pd.DataFrame(sheet_audit).to_csv(OUT / "01_step10G3b_sheet_audit_fixed.tsv", sep="\t", index=False)

    with open(OUT / "06_step10G3b_Zhou2022_gene_sets_fixed.gmt", "w") as gmt:
        for r in records:
            if r["status"] != "ok" or r["n_genes"] == 0:
                continue
            with open(r["path"]) as f:
                genes = [line.strip() for line in f if line.strip()]
            gmt.write(r["geneset"] + "\t" + r["description"] + "\t" + "\t".join(genes) + "\n")

    failed = inv[inv["n_genes"] == 0].copy()
    failed.to_csv(OUT / "07_step10G3b_empty_or_failed_sets_fixed.tsv", sep="\t", index=False)

    with open(OUT / "08_step10G3b_Zhou2022_extraction_fixed_summary.md", "w") as f:
        f.write("# NeuroTRACE Step10G3b Zhou 2022 fixed gene-set extraction summary\n\n")
        f.write(f"Generated: {now()}\n\n")
        f.write("## Purpose\n")
        f.write("Re-extract Zhou 2022 gene sets using fixed sheet-specific header rows and exact column names.\n\n")
        f.write("## Inventory\n\n")
        f.write(inv.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Empty sets\n\n")
        f.write(failed.to_markdown(index=False))
        f.write("\n\n")
        f.write("## Interpretation\n")
        f.write("Use this fixed extraction instead of the previous auto-detected extraction. Previous output used incorrect header rows for multiple sheets and should not be used.\n")

    print(f"[{now()}] Step10G3b fixed extraction done")
    print(f"[{now()}] Non-empty sets: {(inv['n_genes'] > 0).sum()} / {inv.shape[0]}")
    print(f"[{now()}] Results: {OUT}")
    print(f"[{now()}] Gene sets: {GS}")

if __name__ == "__main__":
    main()
