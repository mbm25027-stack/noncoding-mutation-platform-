#!/usr/bin/env python
"""
25_pcawg_expression_test.py

Expression test for non-coding hotspots: within each cohort, is the gene's
expression different in donors carrying a mutation in its hotspot tile(s)
versus donors without one? Rank test + burden-adjusted rank regression,
Benjamini-Hochberg across every gene tested.

Run from the project folder:
    python 25_pcawg_expression_test.py

Edit the CONFIG block first. Nothing here has been run against your data.
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ----------------------------- CONFIG ---------------------------------------
ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"

SNV_FILE = RESULTS / "24_cohort_snvs.csv.gz"   # the cache stage C reported using
HOTSPOT_FILE = None        # None = auto-find a results/*hotspot*.csv from stage A
EXPR_FILE = ROOT / "data" / "tophat_star_fpkm_uq.v2_aliquot_gl.donor.log.gz"  # Xena matrix; if missing, data/ is searched
PROBEMAP_FILE = ROOT / "data" / "gencode.v19.annotation.gene.probemap"        # Ensembl ID -> gene symbol
MAP_FILE = None            # optional 2-column file: expression sample id, donor id (same ids as SNVs)

TILE = 50                  # tile width in bp (BCL2/BCL6 tile numbers imply 50)
COHORTS = ["Lymph-BNHL", "Lymph-CLL", "Skin-Melanoma"]
MIN_MUT = 3                # mutated donors with expression needed to test a gene
MIN_WT = 5                 # unmutated donors with expression needed
# TERT is excluded by request. Immunoglobulin / TCR segments are excluded too:
# their "expression" is clonal rearrangement, not transcription of the gene.
EXCLUDE_RE = re.compile(r"^(TERT$|IG[HKL][VDJC]|TR[ABDG][VDJC])")
PROMOTER_EFFECTS = {"5'Flank", "5'UTR"}   # used only in gene-level fallback mode
OUT_FILE = RESULTS / "25_expression_tests.csv"
# ----------------------------------------------------------------------------


def pick(cols, candidates, what):
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    sys.exit(f"Cannot find a {what} column. Columns present: {list(cols)}")


def norm_chr(s):
    return s.astype(str).str.replace("^chr", "", regex=True)


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def load_snvs():
    head = pd.read_csv(SNV_FILE, nrows=5)
    keep = [c for c in head.columns
            if c.lower() in {"hist", "cohort", "histology", "sample", "donor", "donor_id",
                             "icgc_donor_id", "tumor_sample_barcode", "chr", "chrom",
                             "chromosome", "pos", "start", "start_position",
                             "gene", "effect"}]
    chr_cols = {c: str for c in keep if c.lower() in {"chr", "chrom", "chromosome"}}
    snv = pd.read_csv(SNV_FILE, usecols=keep, dtype=chr_cols)
    c_h = pick(snv.columns, ["hist", "cohort", "histology"], "cohort")
    c_d = pick(snv.columns, ["Sample", "donor", "donor_id", "icgc_donor_id", "Tumor_Sample_Barcode"], "donor")
    c_c = pick(snv.columns, ["chr", "chrom", "Chromosome"], "chromosome")
    c_p = pick(snv.columns, ["pos", "start", "Start_position", "Start_Position"], "position")
    snv = snv.rename(columns={c_h: "hist", c_d: "donor", c_c: "chr", c_p: "pos"})
    snv["chr"] = norm_chr(snv["chr"])
    snv["pos"] = snv["pos"].astype(np.int64)
    snv["tile"] = snv["pos"] // TILE
    return snv[[c for c in ["hist", "donor", "chr", "pos", "tile", "gene", "effect"]
                if c in snv.columns]]


def load_hotspots():
    need = {"cohort": ["cohort", "hist"], "chr": ["chr", "chrom"],
            "tile": ["tile"], "gene": ["gene"]}

    def usable(p):
        try:
            cols = {c.lower() for c in pd.read_csv(p, nrows=2).columns}
        except Exception:
            return False
        return all(any(a in cols for a in alts) for alts in need.values())

    path = HOTSPOT_FILE
    if path is None:
        cands = sorted(RESULTS.glob("*.csv*"))
        good = [p for p in cands if usable(p)]
        if not good:
            print("No file in results/ has cohort, chr, tile and gene columns, so there is no "
                  "stage A hotspot table on disk.\nFalling back to GENE-LEVEL mode: a donor counts as "
                  f"mutated if it carries a SNV with effect in {sorted(PROMOTER_EFFECTS)} in the gene.")
            return None
        path = good[0]
    print(f"Using hotspot table: {path}")
    hs = pd.read_csv(path)
    c_h = pick(hs.columns, ["cohort", "hist"], "cohort")
    c_c = pick(hs.columns, ["chr", "chrom"], "chromosome")
    c_t = pick(hs.columns, ["tile"], "tile")
    c_g = pick(hs.columns, ["gene"], "gene")
    hs = hs.rename(columns={c_h: "hist", c_c: "chr", c_t: "tile", c_g: "gene"})
    hs["chr"] = norm_chr(hs["chr"])
    hs = hs.dropna(subset=["gene"])
    hs = hs[~hs["gene"].astype(str).str.match(EXCLUDE_RE)]
    return hs[["hist", "chr", "tile", "gene"]].drop_duplicates()


def load_expression(snv_donors):
    path = Path(EXPR_FILE)
    if not path.exists():
        alt = sorted((ROOT / "data").glob("*fpkm_uq*donor*")) if (ROOT / "data").exists() else []
        alt = [p for p in alt if "probemap" not in p.name.lower()]
        if not alt:
            sys.exit(f"Expression matrix not found: {EXPR_FILE}\n"
                     "Put the Xena expression matrix (not the probemap) in the data folder, "
                     "or point EXPR_FILE at it.")
        path = alt[0]
        print(f"Using expression matrix: {path}")
    ex = pd.read_csv(path, sep="\t", index_col=0)
    ex = ex[~ex.index.isna()]
    if str(ex.index[0]).startswith("ENSG"):
        pm_path = next((p for p in [Path(PROBEMAP_FILE),
                                    ROOT / "gencode.v19.annotation.gene.probemap"] if p.exists()), None)
        if pm_path is None:
            sys.exit("Expression rows are Ensembl IDs, so the probemap is needed to get gene symbols: "
                     f"{PROBEMAP_FILE} not found.")
        pm = pd.read_csv(pm_path, sep="\t").iloc[:, :2]
        pm.columns = ["id", "gene"]
        lookup = pm.drop_duplicates("id").set_index("id")["gene"]
        sym = ex.index.map(lookup)
        keep = pd.notna(sym)
        print(f"Probemap: {int(keep.sum()):,} of {len(keep):,} Ensembl rows mapped to gene symbols")
        ex = ex.loc[keep]
        ex.index = np.asarray(sym)[keep]
    ex = ex.groupby(level=0).mean(numeric_only=True)
    if MAP_FILE:
        mp = pd.read_csv(MAP_FILE, sep=None, engine="python")
        mp = mp.iloc[:, :2]
        mp.columns = ["sample", "donor"]
        mp = mp.drop_duplicates("sample").set_index("sample")["donor"]
        ex = ex.loc[:, ex.columns.isin(mp.index)]
        ex.columns = mp.loc[ex.columns].values
        ex = ex.T.groupby(level=0).mean().T       # one column per donor
    overlap = len(set(ex.columns) & set(snv_donors))
    print(f"Expression: {ex.shape[0]:,} genes x {ex.shape[1]:,} samples; "
          f"{overlap} sample ids match SNV donor ids")
    if overlap == 0:
        sys.exit("No expression sample ids match the SNV donor ids. PCAWG expression uses aliquot ids, "
                 "so supply MAP_FILE (sample id -> donor id).")
    # Raw counts/TPM get log-transformed; values that already look logged are kept.
    if np.nanmax(ex.values) > 100:
        ex = np.log2(ex + 1)
    return ex


def adjusted_rank_test(y, mut, burden):
    """OLS of rank(expression) on mutated + log10(burden); returns beta, p for 'mutated'."""
    r = stats.rankdata(y)
    X = np.column_stack([np.ones(len(r)), mut.astype(float), np.log10(burden + 1.0)])
    n, k = X.shape
    if n <= k:
        return np.nan, np.nan
    beta, *_ = np.linalg.lstsq(X, r, rcond=None)
    resid = r - X @ beta
    s2 = resid @ resid / (n - k)
    try:
        cov = s2 * np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return np.nan, np.nan
    t = beta[1] / np.sqrt(cov[1, 1])
    return beta[1] / n, 2 * stats.t.sf(abs(t), n - k)   # beta in rank units / n = fraction of rank range


def main():
    snv = load_snvs()
    hs = load_hotspots()
    ex = load_expression(snv["donor"].unique())
    burden = snv.groupby("donor").size()

    # Which donors carry a SNV inside each hotspot tile, per cohort.
    if hs is not None:
        hit = snv.drop(columns=[c for c in ["gene"] if c in snv.columns]).merge(
            hs, on=["hist", "chr", "tile"], how="inner")
    else:
        if "gene" not in snv.columns or "effect" not in snv.columns:
            sys.exit("Gene-level fallback needs 'gene' and 'effect' columns in the SNV file.")
        hit = snv[snv["effect"].isin(PROMOTER_EFFECTS) & snv["gene"].notna()].copy()
        hit = hit[~hit["gene"].astype(str).str.match(EXCLUDE_RE)]
    hit = hit[["hist", "gene", "donor"]].drop_duplicates()

    rows = []
    for cohort in COHORTS:
        donors = set(snv.loc[snv["hist"] == cohort, "donor"]) & set(ex.columns)
        print(f"[{cohort}] donors with expression: {len(donors)}")
        if not donors:
            continue
        for gene, g in hit[hit["hist"] == cohort].groupby("gene"):
            if gene not in ex.index:
                continue
            mut_d = sorted(set(g["donor"]) & donors)
            wt_d = sorted(donors - set(mut_d))
            if len(mut_d) < MIN_MUT or len(wt_d) < MIN_WT:
                continue
            a = ex.loc[gene, mut_d].astype(float).values
            b = ex.loc[gene, wt_d].astype(float).values
            ok_a, ok_b = ~np.isnan(a), ~np.isnan(b)
            a, b = a[ok_a], b[ok_b]
            if len(a) < MIN_MUT or len(b) < MIN_WT:
                continue
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            y = np.concatenate([a, b])
            m = np.concatenate([np.ones(len(a)), np.zeros(len(b))])
            bur = burden.reindex(
                [d for d, o in zip(mut_d, ok_a) if o] + [d for d, o in zip(wt_d, ok_b) if o]
            ).fillna(0).values
            adj_beta, adj_p = adjusted_rank_test(y, m, bur)
            rows.append(dict(
                cohort=cohort, gene=gene, n_mut=len(a), n_wt=len(b),
                median_mut=np.median(a), median_wt=np.median(b),
                median_diff=np.median(a) - np.median(b),
                rank_biserial=2 * u / (len(a) * len(b)) - 1,
                p_mwu=p, adj_effect=adj_beta, p_adj=adj_p,
            ))

    if not rows:
        sys.exit("No gene had enough mutated and unmutated donors with expression. "
                 "Lower MIN_MUT/MIN_WT or check the donor-id mapping.")

    res = pd.DataFrame(rows)
    res["q_mwu"] = bh(res["p_mwu"].values)             # BH across ALL tests, all cohorts
    ok = res["p_adj"].notna()
    res["q_adj"] = np.nan
    res.loc[ok, "q_adj"] = bh(res.loc[ok, "p_adj"].values)
    res = res.sort_values("p_mwu")
    RESULTS.mkdir(exist_ok=True)
    res.to_csv(OUT_FILE, index=False)

    print(f"\nGenes tested: {len(res)}  (BH applied across all of them)")
    print(res.groupby("cohort").agg(tested=("gene", "size"),
                                    q_mwu_lt_005=("q_mwu", lambda s: int((s < 0.05).sum())),
                                    q_adj_lt_005=("q_adj", lambda s: int((s < 0.05).sum()))))
    print("\nTop 15 by raw p:")
    print(res.head(15)[["cohort", "gene", "n_mut", "n_wt", "median_diff",
                        "p_mwu", "q_mwu", "q_adj"]].to_string(index=False))
    print(f"\nFull table: {OUT_FILE}")
    print("NOTE: low power. Check n_mut before reading anything into a q-value.")


if __name__ == "__main__":
    main()