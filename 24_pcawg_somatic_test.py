"""
24_pcawg_somatic_test.py  --  PCAWG somatic test of NC-MutScan calls (TERT excluded)

Usage (from D:\\noncoding mutation finding surver, venv active):
    python -m pip install pyliftover          # once, before stage B/C; needs internet for the chain file
    python 24_pcawg_somatic_test.py A         # full run of stage A (first run builds the cache, takes minutes)
    python 24_pcawg_somatic_test.py B         # then B, C, D

Stages (each writes its own files to results\\tables):
  A  NB hotspot calling (non-coding SNVs only) in native hg19 coordinates. No liftover needed.
  B  Liftover to GRCh38 + engine scoring of hotspot / recurrent / singleton mutations,
     then CMH enrichment of engine calls (strata = cohort x trinucleotide context).
     NOTE: descriptive only. UV-driven hotspots at ETS/CTCF sites are expected to score as
     "motif-disrupting" without any selection.
  C  AID/SHM mask: SHM-target genes vs other lymphoma genes vs UV loci vs random windows.
  D  Expression: mutated vs WT OLS (CNV, purity, log10 TMB), within cohort, BH-corrected.
"""
import argparse
import importlib.util
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import fisher_exact, nbinom
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

# ============================== CONFIG (edit here) ===========================
BASE = Path(r"D:\noncoding mutation finding surver")
MAF = BASE / "October_2016_whitelist_2583.snv_mnv_indel.maf.xena.nonUS"
DONOR_FILE = BASE / "pcawg_specimen_histology_August2016_v9_donor"   # must hold DO ids + histology_abbreviation
EXT = BASE / "results" / "external_pcawg"
FASTA = BASE / "GRCh38.d1.vd1.fa" / "GRCh38.d1.vd1.fa"
CHAIN = BASE / "hg19ToHg38.over.chain.gz"      # optional local chain; else pyliftover downloads it
ENGINE_FILE = BASE / "nc_mutscan_research_suite.py"
ENGINE_FUNC = "analyze_variant_delta"          # <-- SET to the engine's WT-vs-MUT scoring function
GENE_MAP = BASE / "gene_symbol_to_ensg.tsv"    # optional: two columns, symbol <tab> ENSG (no header)
OUT = BASE / "results" / "tables"
CACHE = BASE / "results" / "24_cohort_snvs.csv.gz"
OUT.mkdir(parents=True, exist_ok=True)

COHORTS = ["Skin-Melanoma", "Lymph-BNHL", "Lymph-CLL"]
LYMPH = {"Lymph-BNHL", "Lymph-CLL"}
TILE, GENOME_TILES, FLANK = 50, 56_912_480, 25
FDR, MIN_K, MAX_HOT, CTRL_X, SEED = 0.05, 3, 20000, 5, 1
TERT_HG19 = ("5", 1_293_000, 1_297_000)        # TERT promoter region, excluded everywhere
# Known aSHM targets in B-cell lymphoma: verify this list against the paper you cite.
SHM_GENES = ["BCL6", "BCL2", "PAX5", "MYC", "PIM1", "CD79B", "CD79A", "RHOH", "SOCS1",
             "CIITA", "BTG2", "CD83", "BACH2", "SGK1", "IGLL5"]
UV_GENES = ["PLEKHS1", "WDR74", "DPH3", "RPL13A"]
ENSG = {"BCL6": "ENSG00000113916", "BCL2": "ENSG00000171791", "PAX5": "ENSG00000196092",
        "MYC": "ENSG00000136997", "PIM1": "ENSG00000137193", "CD79B": "ENSG00000007312",
        "WDR74": "ENSG00000133316"}             # VERIFY against the expression file
# A label is "coding" if it matches CODING_PAT and does NOT match NONCODING_PAT.
CODING_PAT = r"missense|nonsense|silent|synonymous|frame|stop|start_lost|start_gained|coding_sequence|cds"
NONCODING_PAT = r"non[_ -]?coding|utr|intron|intergenic|upstream|downstream|promoter|flank|regulatory"
COMP = str.maketrans("ACGT", "TGCA")
# Mirrors the engine exactly (nc_mutscan_research_suite.py lines 366-368): non-overlapping
# re.findall of WRCY and RGYW on the MUT sequence, mask fires when the sum is >= 6.
WRCY_RE = re.compile(r"[AT][AG]C[CT]")
RGYW_RE = re.compile(r"[AG]G[CT][AT]")
MASK_MIN = 6
# =============================================================================


# ------------------------------- shared helpers ------------------------------
def noncoding(df):
    eff = df["effect"].fillna("")
    coding = eff.str.contains(CODING_PAT, case=False, regex=True) & ~eff.str.contains(NONCODING_PAT, case=False, regex=True)
    return ~coding


def load_hist():
    h = pd.read_csv(DONOR_FILE, sep="\t", usecols=["icgc_specimen_id", "histology_abbreviation"], dtype=str)
    h = h.dropna(subset=["icgc_specimen_id", "histology_abbreviation"])
    h["icgc_specimen_id"] = h["icgc_specimen_id"].str.strip()
    h["histology_abbreviation"] = h["histology_abbreviation"].str.strip()
    h = h[h["histology_abbreviation"] != ""]
    return dict(zip(h["icgc_specimen_id"], h["histology_abbreviation"]))


def load_snvs():
    if CACHE.exists():
        print(f"Using cached SNVs: {CACHE}")
        df = pd.read_csv(CACHE, dtype=str)
        df["start"] = df["start"].astype(np.int64)
        return df
    hist = load_hist()
    print(f"Donor histology map: {len(hist)} ids; histologies: {sorted(set(hist.values()))}")
    if sum(k.startswith("DO") for k in hist) < 0.5 * len(hist):
        sys.exit("DONOR_FILE keys are not DO ids (looks like the specimen file). Point DONOR_FILE at the donor file.")
    missing = [c for c in COHORTS if c not in set(hist.values())]
    if missing:
        sys.exit(f"Cohort names not found in histology map: {missing}. Edit COHORTS.")
    parts = []
    rd = pd.read_csv(MAF, sep="\t", dtype=str, chunksize=1_000_000,
                     usecols=["Sample", "chr", "start", "reference", "alt", "gene", "effect"])
    for i, ch in enumerate(rd):
        ch["hist"] = ch["Sample"].map(hist)
        ch = ch[ch["hist"].isin(COHORTS)]
        ch = ch[ch["reference"].isin(list("ACGT")) & ch["alt"].isin(list("ACGT"))]
        if len(ch):
            parts.append(ch)
        print(f"  chunk {i + 1} read, kept so far: {sum(len(p) for p in parts):,}")
    if not parts:
        sys.exit("No SNVs matched the three cohorts. Check that MAF 'Sample' ids match DONOR_FILE keys.")
    df = pd.concat(parts, ignore_index=True)
    df["chr"] = df["chr"].str.replace("chr", "", regex=False)
    df["start"] = df["start"].astype(np.int64)
    is_tert = ((df["chr"] == TERT_HG19[0]) & df["start"].between(TERT_HG19[1], TERT_HG19[2])) | (df["gene"] == "TERT")
    print(f"Dropped {int(is_tert.sum())} TERT-region SNVs")
    df = df[~is_tert].reset_index(drop=True)
    print(df.groupby("hist")["Sample"].agg(["size", "nunique"]).rename(columns={"size": "SNVs", "nunique": "donors"}))
    print("Top effect labels (check the coding filter):")
    print(df["effect"].value_counts().head(15))
    print(f"Classified non-coding: {noncoding(df).mean():.1%} of SNVs")
    df.to_csv(CACHE, index=False)
    return df


def bh_sorted(p, m):
    """BH q-values for the tested subset, with m = total number of tests (all tiles)."""
    o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * m / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1.0)
    return out


def load_engine():
    if not ENGINE_FILE.exists():
        sys.exit(f"Set ENGINE_FILE: {ENGINE_FILE} not found")
    spec = importlib.util.spec_from_file_location("ncms_engine", ENGINE_FILE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)   # runs engine top-level code: the engine needs an if __name__ == "__main__" guard
    if not hasattr(mod, ENGINE_FUNC):
        funcs = [n for n in dir(mod) if callable(getattr(mod, n)) and not n.startswith("_")]
        sys.exit(f"Engine has no function '{ENGINE_FUNC}'. Available: {funcs}. Set ENGINE_FUNC.")
    return mod


def get_lifter():
    from pyliftover import LiftOver
    return LiftOver(str(CHAIN)) if CHAIN.exists() else LiftOver("hg19", "hg38")


def lift_fetch(lo, fa, chrom, pos, ref, alt):
    """hg19 1-based pos -> hg38 window. Reference-base check validates the hg19 assumption."""
    res = lo.convert_coordinate("chr" + chrom, int(pos) - 1)
    if not res or len(res) != 1:
        return None, "nolift"
    c, p0, strand, _ = res[0]
    if p0 < FLANK:
        return None, "edge"
    try:
        seq = fa[c][p0 - FLANK:p0 + FLANK + 1].seq.upper()
    except Exception:
        return None, "fetch"
    if len(seq) != 2 * FLANK + 1:
        return None, "edge"
    if strand == "-":
        ref, alt = ref.translate(COMP), alt.translate(COMP)
    if seq[FLANK] != ref:
        return None, "refmismatch"
    return {"wt": seq, "mut": seq[:FLANK] + alt + seq[FLANK + 1:], "alt38": alt}, "ok"


def ctx96(wt, alt):
    tri = wt[FLANK - 1:FLANK + 2]
    if tri[1] in "AG":
        tri, alt = tri.translate(COMP)[::-1], alt.translate(COMP)
    return f"{tri}>{alt}"


def score_pair(eng, wt, mut):
    """Differential (WT vs MUT) path with neutral labels, so headers cannot leak the class."""
    r = getattr(eng, ENGINE_FUNC)(wt, mut, flank_bp=FLANK, locus="NA", title="NA")
    if not isinstance(r, dict):
        sys.exit(f"{ENGINE_FUNC} returned {type(r).__name__}, expected a dict. Send me the engine function so I can map the keys.")
    n_mask = len(WRCY_RE.findall(mut)) + len(RGYW_RE.findall(mut))
    return r.get("classification"), r.get("delta_pwm"), r.get("primary_tf"), n_mask


def score_table(rows, eng, lo, fa, label):
    out, status = [], {}
    for i, r in enumerate(rows.itertuples(index=False)):
        info, st = lift_fetch(lo, fa, r.chr, r.start, r.reference, r.alt)
        status[st] = status.get(st, 0) + 1
        if info is None:
            continue
        cls, dp, tf, nm = score_pair(eng, info["wt"], info["mut"])
        out.append({**r._asdict(), "classification": cls, "delta_pwm": dp, "primary_tf": tf,
                    "n_wrcy_rgyw": nm, "mask_trigger": nm >= MASK_MIN, "ctx96": ctx96(info["wt"], info["alt38"])})
        if (i + 1) % 5000 == 0:
            print(f"  [{label}] scored {i + 1}/{len(rows)}")
    tot = sum(status.values()) or 1
    print(f"  [{label}] liftover/ref status: {status}  (ok = {status.get('ok', 0) / tot:.1%})")
    if status.get("ok", 0) / tot < 0.9:
        print("  WARNING: <90% reference-base match. The MAF is probably not hg19, or the chain is wrong. Stop and check.")
    return pd.DataFrame(out)


# ------------------------------------ Stage A --------------------------------
def stage_A(df):
    all_tiles, tagged = [], []
    for cohort, d in df.groupby("hist"):
        d = d[noncoding(d)].copy()                       # hotspots are called on non-coding SNVs only
        d["tile"] = d["start"] // TILE
        dt = d.drop_duplicates(["Sample", "chr", "tile"]).groupby(["chr", "tile"]).size().rename("k").reset_index()
        M = GENOME_TILES
        mean = dt["k"].sum() / M
        var = (dt["k"] ** 2).sum() / M - mean ** 2
        r = min(1e6, mean ** 2 / max(var - mean, 1e-12))
        # local rate: shrunken mean of the surrounding 1 Mb (excluding the tile itself)
        dt["bin"] = dt["tile"] * TILE // 1_000_000
        binsum = dt.groupby(["chr", "bin"])["k"].transform("sum")
        W, per_bin = 2000, 1_000_000 // TILE
        dt["lam0"] = (binsum - dt["k"] + W * mean) / (per_bin - 1 + W)
        t = dt[dt["k"] >= 2].copy()
        t["p"] = np.maximum(nbinom.sf(t["k"] - 1, r, r / (r + t["lam0"])), 1e-300)
        t["q"] = bh_sorted(t["p"].to_numpy(), M)
        t["cohort"] = cohort
        t["hotspot"] = (t["q"] < FDR) & (t["k"] >= MIN_K)
        print(f"\n[{cohort}] donors={d['Sample'].nunique()} tiles with k>=1: {len(dt):,}  mean={mean:.4f}  r={r:.3f}")
        for kk in (2, 3, 4, 5):   # calibration: observed vs expected recurrent tiles under the global NB
            exp = M * nbinom.sf(kk - 1, r, r / (r + mean))
            print(f"   tiles with k>={kk}: observed {int((dt['k'] >= kk).sum()):,}  expected {exp:,.0f}")
        print(f"   hotspots (q<{FDR}, k>={MIN_K}): {int(t['hotspot'].sum())}")
        g = d.merge(t[["chr", "tile"]], on=["chr", "tile"]).groupby(["chr", "tile"])["gene"].agg(
            lambda s: s.dropna().mode().iat[0] if len(s.dropna()) else "")
        t = t.merge(g.rename("gene").reset_index(), on=["chr", "tile"], how="left")
        all_tiles.append(t)
        m = d.merge(dt[["chr", "tile", "k"]], on=["chr", "tile"], how="left").merge(
            t[["chr", "tile", "q", "hotspot"]], on=["chr", "tile"], how="left")
        m["hotspot"] = m["hotspot"].fillna(False).astype(bool)
        m["group"] = np.where(m["hotspot"], "hotspot", np.where(m["k"] >= 2, "recurrent", "singleton"))
        sing = m[m["group"] == "singleton"]
        sing = sing.sample(min(len(sing), 100_000), random_state=SEED)
        tagged.append(pd.concat([m[m["group"] != "singleton"], sing]))
    tiles = pd.concat(all_tiles, ignore_index=True)
    tiles.to_csv(OUT / "24a_hotspot_tiles.csv", index=False)
    pd.concat(tagged, ignore_index=True).to_csv(OUT / "24a_tagged_mutations.csv.gz", index=False)
    top = tiles[tiles["hotspot"]].sort_values("p").head(20)
    print("\nTop hotspots:\n", top[["cohort", "chr", "tile", "k", "lam0", "p", "q", "gene"]].to_string(index=False))


# ------------------------------------ Stage B --------------------------------
def cmh(sub, outcome):
    """CMH across strata (cohort x trinucleotide context): hotspot-type vs control."""
    tabs = []
    for _, s in sub.groupby("stratum"):
        a, b = s[s["is_case"]], s[~s["is_case"]]
        if len(a) == 0 or len(b) == 0:
            continue
        tabs.append(np.array([[a[outcome].sum(), (~a[outcome]).sum()], [b[outcome].sum(), (~b[outcome]).sum()]]))
    if not tabs:
        return None
    st = sm.stats.StratifiedTable(tabs)
    lo, hi = st.oddsratio_pooled_confint()
    return {"OR_MH": st.oddsratio_pooled, "CI_lo": lo, "CI_hi": hi, "p_CMH": st.test_null_odds(correction=True).pvalue,
            "n_strata": len(tabs)}


def stage_B():
    tag = pd.read_csv(OUT / "24a_tagged_mutations.csv.gz", dtype={"chr": str, "Sample": str, "gene": str})
    hot, rec, sing = (tag[tag["group"] == g] for g in ("hotspot", "recurrent", "singleton"))
    if hot.empty:
        sys.exit("No hotspots from stage A; nothing to score.")
    hot = hot.sample(min(len(hot), MAX_HOT), random_state=SEED)
    rec = rec.sample(min(len(rec), MAX_HOT), random_state=SEED)
    sing = sing.sample(min(len(sing), CTRL_X * len(hot)), random_state=SEED)
    from pyfaidx import Fasta
    fa, lo, eng = Fasta(str(FASTA)), get_lifter(), load_engine()
    sc = pd.concat([score_table(x, eng, lo, fa, n) for x, n in ((hot, "hotspot"), (rec, "recurrent"), (sing, "singleton"))],
                   ignore_index=True)
    sc["driver"] = sc["classification"].eq("Driver Variant")
    sc["functional"] = sc["classification"].isin(["Driver Variant", "Loss-of-Function Variant"])
    sc["stratum"] = sc["hist"] + "|" + sc["ctx96"]
    sc.to_csv(OUT / "24b_scored_calls.csv.gz", index=False)
    print("\nClassification by group:\n", pd.crosstab(sc["group"], sc["classification"]))
    rows = []
    for case in ("hotspot", "recurrent"):
        sub = sc[sc["group"].isin([case, "singleton"])].copy()
        sub["is_case"] = sub["group"] == case
        for outcome in ("driver", "functional"):
            res = cmh(sub, outcome)
            if res:
                rows.append({"case": case, "control": "singleton", "outcome": outcome,
                             "case_rate": sub[sub.is_case][outcome].mean(), "control_rate": sub[~sub.is_case][outcome].mean(),
                             "n_case": int(sub.is_case.sum()), "n_control": int((~sub.is_case).sum()), **res})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "24b_enrichment.csv", index=False)
    print("\nEnrichment (CMH, matched on cohort x trinucleotide context):\n", res.to_string(index=False))
    print("\nNOTE: descriptive only. UV hotspots at ETS/CTCF sites score as motif-disrupting without any selection.")


# ------------------------------------ Stage C --------------------------------
def summarize_group(name, s):
    n = len(s)
    ci = proportion_confint(int(s["mask_trigger"].sum()), n, method="wilson") if n else (np.nan, np.nan)
    return {"group": name, "n": n,
            "mask_trigger_rate": s["mask_trigger"].mean(), "trigger_CI_lo": ci[0], "trigger_CI_hi": ci[1],
            "labelled_passenger": s["classification"].eq("AID/SHM Passenger").mean(),
            "driver_rate": s["classification"].eq("Driver Variant").mean(),
            "LoF_rate": s["classification"].eq("Loss-of-Function Variant").mean(),
            "mask_missed_functional": (s["mask_trigger"] & s["classification"].isin(
                ["Driver Variant", "Loss-of-Function Variant"])).mean()}


def stage_C(df):
    from pyfaidx import Fasta
    fa, lo, eng = Fasta(str(FASTA)), get_lifter(), load_engine()
    rng = np.random.default_rng(SEED)
    nc = df[noncoding(df)]
    lym = nc[nc["hist"].isin(LYMPH)]
    shm = lym[lym["gene"].isin(SHM_GENES)]
    oth = lym[~lym["gene"].isin(SHM_GENES) & lym["gene"].notna()]
    uv = nc[(nc["hist"] == "Skin-Melanoma") & nc["gene"].isin(UV_GENES)]
    print(f"Group sizes before sampling: shm={len(shm)} lymph_other={len(oth)} uv={len(uv)}")
    n = min(5000, len(shm)) or 1
    groups = {"shm_targets": shm.sample(min(len(shm), 5000), random_state=SEED),
              "lymph_other_genes": oth.sample(min(len(oth), n), random_state=SEED),
              "uv_loci": uv.sample(min(len(uv), 5000), random_state=SEED)}
    scored = []
    for name, rows in groups.items():
        s = score_table(rows, eng, lo, fa, name)
        s["group"] = name
        scored.append(s)
    # random hg38 windows with random substitutions: false-flag rate of the mask
    chroms = [f"chr{i}" for i in range(1, 23)]
    sizes = np.array([len(fa[c]) for c in chroms], dtype=float)
    rnd = []
    while len(rnd) < 5000:
        c = rng.choice(chroms, p=sizes / sizes.sum())
        p0 = int(rng.integers(FLANK, len(fa[c]) - FLANK - 1))
        seq = fa[c][p0 - FLANK:p0 + FLANK + 1].seq.upper()
        if len(seq) != 2 * FLANK + 1 or set(seq) - set("ACGT"):
            continue
        alt = rng.choice([b for b in "ACGT" if b != seq[FLANK]])
        mut = seq[:FLANK] + alt + seq[FLANK + 1:]
        cls, dp, tf, nm = score_pair(eng, seq, mut)
        rnd.append({"classification": cls, "delta_pwm": dp, "primary_tf": tf, "n_wrcy_rgyw": nm,
                    "mask_trigger": nm >= MASK_MIN, "group": "random_windows"})
    scored.append(pd.DataFrame(rnd))
    sc = pd.concat([s for s in scored if len(s)], ignore_index=True)
    sc.to_csv(OUT / "24c_mask_scored.csv.gz", index=False)
    summ = pd.DataFrame([summarize_group(g, s) for g, s in sc.groupby("group")])
    summ.to_csv(OUT / "24c_mask_summary.csv", index=False)
    print("\nMask summary:\n", summ.to_string(index=False))
    tests = []
    for a, b in (("shm_targets", "lymph_other_genes"), ("shm_targets", "random_windows"), ("uv_loci", "random_windows")):
        for col in ("mask_trigger", "passenger"):
            sa, sb = sc[sc["group"] == a], sc[sc["group"] == b]
            fa_ = sa["mask_trigger"] if col == "mask_trigger" else sa["classification"].eq("AID/SHM Passenger")
            fb_ = sb["mask_trigger"] if col == "mask_trigger" else sb["classification"].eq("AID/SHM Passenger")
            if len(sa) and len(sb):
                t = [[int(fa_.sum()), int((~fa_).sum())], [int(fb_.sum()), int((~fb_).sum())]]
                orr, p = fisher_exact(t)
                tests.append({"group_a": a, "group_b": b, "outcome": col, "rate_a": fa_.mean(), "rate_b": fb_.mean(), "OR": orr, "p": p})
    tdf = pd.DataFrame(tests)
    tdf.to_csv(OUT / "24c_mask_tests.csv", index=False)
    print("\nFisher tests:\n", tdf.to_string(index=False))
    print("\nNOTE: in the engine the mask is evaluated AFTER the splice/TF checks, so it can never turn a Driver/LoF "
          "call into a Passenger. 'mask_missed_functional' counts mutations whose window triggers the mask but were "
          "already called Driver/LoF.")


# ------------------------------------ Stage D --------------------------------
def read_rows(path, bases, transform=None):
    """Stream a gene x sample table, keep rows whose base id/symbol is in `bases`."""
    out = {}
    for ch in pd.read_csv(path, sep="\t", index_col=0, chunksize=3000):
        base = ch.index.to_series().astype(str).str.split("|").str[0].str.split(".").str[0]
        keep = base.isin(bases).to_numpy()
        for b, (_, row) in zip(base[keep], ch[keep].iterrows()):
            v = pd.to_numeric(row, errors="coerce")
            out[b] = transform(v) if transform else v
    return out


def stage_D(df):
    ens = dict(ENSG)
    if GENE_MAP.exists():
        gm = pd.read_csv(GENE_MAP, sep="\t", header=None, dtype=str)
        ens.update(dict(zip(gm[0], gm[1])))
    else:
        print("WARNING: GENE_MAP missing; only the hard-coded ENSG ids resolve, other genes are skipped unless the "
              "expression file is keyed by symbol.")
    hot = pd.read_csv(OUT / "24a_hotspot_tiles.csv", dtype={"chr": str})
    hot = hot[hot["hotspot"] & hot["gene"].notna() & (hot["gene"] != "")]
    d = df.copy()
    d["tile"] = d["start"] // TILE
    nc_d = d[noncoding(d)]
    tests = []
    hd = nc_d.merge(hot[["cohort", "chr", "tile", "gene"]].rename(columns={"gene": "tgene", "cohort": "hist"}),
                    on=["hist", "chr", "tile"])
    for (h, c, t, g), s in hd.groupby(["hist", "chr", "tile", "tgene"]):
        tests.append({"cohort": h, "type": "hotspot_tile", "gene": g, "chr": c, "pos": int(t * TILE), "mut": set(s["Sample"])})
    nc = nc_d[nc_d["hist"].isin(LYMPH) & nc_d["gene"].isin(SHM_GENES)]
    for (h, g), s in nc.groupby(["hist", "gene"]):
        tests.append({"cohort": h, "type": "shm_gene", "gene": g, "chr": s["chr"].iloc[0], "pos": int(s["start"].median()),
                      "mut": set(s["Sample"])})
    tests = [t for t in tests if len(t["mut"]) >= 3]
    print(f"{len(tests)} candidate tests with >=3 mutated donors")

    syms = {t["gene"] for t in tests}
    keys = {ens.get(s, s) for s in syms} | syms
    expr = read_rows(EXT / "pcawg_rnaseq_fpkm_uq.tsv.gz", keys, lambda v: np.log2(v + 1.0))
    expr_ids = set(next(iter(expr.values())).index) if expr else set()
    tmb = d.groupby("Sample").size()   # genome-wide SNV count per donor (cohort donors only)
    cohort_of = d.drop_duplicates("Sample").set_index("Sample")["hist"]
    print(f"Expression genes found: {len(expr)}; donors in RNA-seq: {len(expr_ids)}; overlap with MAF donors: {len(expr_ids & set(tmb.index))}")
    if not expr_ids & set(tmb.index):
        sys.exit("No donor-ID overlap between RNA-seq columns and MAF Sample. Print the first RNA-seq header ids and fix the mapping.")

    pur = pd.read_csv(EXT / "pcawg_purity_ploidy.tsv.gz", sep="\t")
    pcol = next(c for c in pur.columns if "purity" in c.lower())
    purity = pd.Series(pd.to_numeric(pur[pcol], errors="coerce").values, index=pur.iloc[:, 0].astype(str))
    print(f"Purity: {len(purity)} rows, overlap with MAF donors: {len(set(purity.index) & set(tmb.index))}")

    hdr = pd.read_csv(EXT / "pcawg_consensus_cnv.tsv.gz", sep="\t", nrows=0).columns
    low = [c.lower() for c in hdr]
    long_fmt = any(c in low for c in ("chr", "chrom", "chromosome")) and "start" in low and "end" in low
    if long_fmt:
        cn = pd.read_csv(EXT / "pcawg_consensus_cnv.tsv.gz", sep="\t")
        cn.columns = [c.lower() for c in cn.columns]
        cc = next(c for c in ("chr", "chrom", "chromosome") if c in cn.columns)
        sc_ = next(c for c in ("sampleid", "samplename", "sample", "donor", "icgc_donor_id") if c in cn.columns)
        vc = next(c for c in ("total_cn", "value", "copy_number", "segment_mean") if c in cn.columns)
        cn[cc] = cn[cc].astype(str).str.replace("chr", "", regex=False)
        print(f"CNV long format: sample={sc_}, value={vc}; donor overlap {len(set(cn[sc_].astype(str)) & set(tmb.index))}")
        cnv_rows = None
    else:
        cnv_rows = read_rows(EXT / "pcawg_consensus_cnv.tsv.gz", keys)
        print(f"CNV wide format: {len(cnv_rows)} genes found")

    out = []
    for t in tests:
        y_all = expr.get(ens.get(t["gene"], t["gene"]))
        if y_all is None:
            y_all = expr.get(t["gene"])
        if y_all is None:
            continue
        donors = [x for x in y_all.index if x in cohort_of.index and cohort_of[x] == t["cohort"] and x in tmb.index]
        y = y_all.loc[donors].astype(float)
        if long_fmt:
            seg = cn[(cn[cc] == t["chr"]) & (cn["start"] <= t["pos"]) & (cn["end"] >= t["pos"])]
            cnv = seg.drop_duplicates(sc_).set_index(sc_)[vc].reindex(donors)
        else:
            row = cnv_rows.get(ens.get(t["gene"], t["gene"]), cnv_rows.get(t["gene"]))
            cnv = row.reindex(donors) if row is not None else pd.Series(np.nan, index=donors)
        pu = purity.reindex(donors)
        X = pd.DataFrame({"const": 1.0, "mut": [1.0 if x in t["mut"] else 0.0 for x in donors],
                          "cnv": cnv.values, "purity": pu.values, "log10_tmb": np.log10(tmb.reindex(donors).values + 1.0)},
                         index=donors)
        n_imp = {"cnv": int(X["cnv"].isna().sum()), "purity": int(X["purity"].isna().sum())}
        for c in ("cnv", "purity"):
            if X[c].isna().all():
                X = X.drop(columns=c)          # covariate unavailable for this cohort
            else:
                X[c] = X[c].fillna(X[c].median())
        n_mut = int(X["mut"].sum())
        if n_mut < 3 or len(X) - n_mut < 10:
            continue
        ok = y.notna().to_numpy()
        fit = sm.OLS(y.values[ok], X.values[ok]).fit(cov_type="HC3")
        i = list(X.columns).index("mut")
        out.append({"cohort": t["cohort"], "type": t["type"], "gene": t["gene"], "locus_hg19": f"chr{t['chr']}:{t['pos']}",
                    "n_wt": len(X) - n_mut, "n_mut": n_mut, "mean_wt": y[X["mut"] == 0].mean(), "mean_mut": y[X["mut"] == 1].mean(),
                    "beta_mut": fit.params[i], "se": fit.bse[i], "p": fit.pvalues[i],
                    "n_cnv_imputed": n_imp["cnv"], "n_purity_imputed": n_imp["purity"], "covariates": ",".join(X.columns[2:])})
    res = pd.DataFrame(out)
    if res.empty:
        sys.exit("No testable gene/cohort pairs (check ID mapping printouts above).")
    res["q_BH_all"] = multipletests(res["p"], method="fdr_bh")[1]
    res["q_BH_within_type"] = res.groupby("type")["p"].transform(lambda s: multipletests(s, method="fdr_bh")[1])
    res.sort_values("p").to_csv(OUT / "24d_expression_eqtl.csv", index=False)
    print(res.sort_values("p").head(25).to_string(index=False))


# ------------------------------------ main -----------------------------------
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["A", "B", "C", "D", "all"])
    a = ap.parse_args()
    data = load_snvs() if a.stage in ("A", "C", "D", "all") else None
    if a.stage in ("A", "all"):
        stage_A(data)
    if a.stage in ("B", "all"):
        stage_B()
    if a.stage in ("C", "all"):
        stage_C(data)
    if a.stage in ("D", "all"):
        stage_D(data)