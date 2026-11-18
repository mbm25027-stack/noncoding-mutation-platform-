#!/usr/bin/env python
"""
26_mpra_expanded_test.py

Expanded MPRA benchmark of motif-based variant scoring on saturation-mutagenesis
MPRA data (Kircher et al. 2019; portal: https://kircherlab.bihealth.org/satMutMPRA/).

What it does
  1. Scores every MPRA variant against the FULL JASPAR CORE vertebrate collection
     (data/JASPAR*) and, if present, HOCOMOCO v12 (data/H12CORE*) PFMs.
     dPWM = best-window log2-odds(alt) - best-window log2-odds(ref), both strands.
  2. Keeps only TFs expressed (TPM >= EXPR_TPM) in the cell line each element was assayed in.
  3. Direction-aware scoring: the sign of dPWM is compared with the sign of the MPRA
     effect, with Spearman correlation as the rank statistic.
  4. Leave-one-element-out (LOEO): motif weights (signed) are learned on all other
     elements and applied to the held-out one. Elements that share genomic sequence
     (e.g. TERT-*, SORT1*, PKLR-24h/48h, LDLR/LDLR.2) are held out TOGETHER.
  5. Controls: no expression filter, and a permutation null for the filter.
  6. cCRE prior: does "variant lies in an ENCODE cCRE" (PLS / ELS / CA-CTCF) enrich for
     MPRA-active variants (LR+, group-bootstrap CI)? Compared with the LOEO motif call
     and their combination.

Inputs (all under data/):
  mpra/<element>.tsv[.gz]   one file per element from the portal (GRCh38 download).
                            Columns: Chromosome, Position, Ref, Alt, Value, P-Value
  JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar   (already in your data/)
  H12CORE_jaspar_format.txt (optional; HOCOMOCO v12 CORE, JASPAR/PFM format)
  cell-line expression table (genes x cell lines, TPM), e.g. CCLE or ENCODE;
                            auto-found if the file name contains ccle / tpm / cellline
  GRCh38-cCREs.PLS.bed, GRCh38-cCREs.ELS.bed, GRCh38-cCREs.CA-CTCF.bed (already in your data/)
  gencode.v19.annotation.gene.probemap (only used if the expression rows are ENSG ids)

Run:  python 26_mpra_expanded_test.py
"""
import re
import sys
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ----------------------------- CONFIG ---------------------------------------
ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
MPRA_DIR = DATA / "mpra"
EXPR_FILE = None                    # None = auto-find in data/
GENE_MAP_FILE = None                # None = first data/*probemap*
CCRE_FILES = {"PLS": DATA / "GRCh38-cCREs.PLS.bed",
              "ELS": DATA / "GRCh38-cCREs.ELS.bed",
              "CTCF": DATA / "GRCh38-cCREs.CA-CTCF.bed"}

EXPR_TPM = 1.0                      # a TF is "expressed" at TPM >= this
MIN_REL = 0.80                      # a motif hit needs relative score >= this (0..1) in ref or alt
PSEUDOCOUNT = 0.8
TOP_K = 50                          # motifs kept per fold (largest |training correlation|)
MIN_TRAIN_N = 30                    # variants with nonzero dPWM needed to estimate a motif weight
MIN_TEST = 50                       # variants needed to evaluate a held-out element
SIG_P = 0.05                        # variants called significant for the sign test
ACTIVE_Q = 0.05                     # BH q (pooled) defining "MPRA-active" for LR+
CALL_FRAC = 0.10                    # motif call = top CALL_FRAC of |S| within each element
N_PERM = 20                         # permutations for the expression-filter null
N_BOOT = 1000
SEED = 1
# Elements to drop (lower-case names). Empty = keep all. Add the TERT-* names to match
# the earlier hotspot analysis, where TERT was excluded.
EXCLUDE_ELEMENTS = set()

# Cell line each element was assayed in (from the paper's element tables). Values are
# alias lists matched against expression-table column names (case/punctuation-insensitive;
# CCLE-style NAME_TISSUE columns match on the part before the first underscore).
# Elements with None, or with no matching column, are skipped.
ELEMENT_CELL = {
    "f9": ["HEPG2"], "ldlr": ["HEPG2"], "ldlr.2": ["HEPG2"],
    "sort1": ["HEPG2"], "sort1-flip": ["HEPG2"], "sort1.2": ["HEPG2"],
    "pklr-24h": ["K562"], "pklr-48h": ["K562"],
    "hnf4a": ["HEK293T", "293T", "HEK293"], "msmb": ["HEK293T", "293T", "HEK293"],
    "mycrs6983267": ["HEK293T", "293T", "HEK293"], "tert-hek": ["HEK293T", "293T", "HEK293"],
    "foxe1": ["HELA"],
    "gp1ba": ["HEL9217", "HEL"], "hbb": ["HEL9217", "HEL"], "hbg1": ["HEL9217", "HEL"],
    "bcl11a": ["HEL9217", "HEL"],
    "irf4": ["SKMEL28"], "irf6": ["HACAT"],
    "mycrs11986220": ["LNCAP", "LNCAPCLONEFGC"],
    "ret": ["NEURO2A", "N2A"], "uc88": ["NEURO2A", "N2A"],
    "tcf7l2": ["MIN6"], "zfand3": ["MIN6"],
    "zrsh-13": ["NIH3T3", "3T3"], "zrsh-13h2": ["NIH3T3", "3T3"],
    "tert-gbm": ["SF7996"], "tert-gaa": None, "tert-gsc": None,
}
# ----------------------------------------------------------------------------

BASES = "ACGT"
CODE = {b: i for i, b in enumerate(BASES)}


def norm_name(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def pick(cols, candidates, what):
    low = {re.sub(r"[^a-z0-9]", "", c.lower()): c for c in cols}
    for cand in candidates:
        k = re.sub(r"[^a-z0-9]", "", cand.lower())
        if k in low:
            return low[k]
    sys.exit(f"Cannot find a {what} column. Columns present: {list(cols)}")


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


# ----------------------------- motifs ---------------------------------------
def tf_genes(name):
    parts = name.split(None, 1)
    tf = parts[1] if len(parts) == 2 else parts[0].split(".")[0]
    tf = re.sub(r"\(var\.\d+\)", "", tf)
    return tuple(g.strip().upper() for g in tf.split("::") if g.strip())


def parse_pfm_file(path):
    """JASPAR/HOCOMOCO PFM text: '>ID NAME' then four rows A,C,G,T."""
    mats, name, rows = [], None, []

    def flush():
        if name is not None and len(rows) == 4 and len({len(r) for r in rows}) == 1:
            mats.append((name, np.array(rows, float).T))

    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                name, rows = line[1:].strip(), []
            else:
                body = line[1:] if line[0] in "ACGTacgt" else line
                nums = re.findall(r"-?\d+\.?\d*(?:[eE]-?\d+)?", body.replace("[", " ").replace("]", " "))
                if nums:
                    rows.append([float(x) for x in nums])
    flush()
    return mats


def load_motifs():
    def looks_like_pfm(p):
        try:
            if p.stat().st_size > 50_000_000:
                return False
            opener = gzip.open if p.suffix.lower() == ".gz" else open
            with opener(p, "rt", errors="ignore") as fh:
                head = fh.read(4000)
            return bool(re.match(r"\s*>", head)) and bool(re.search(r"^\s*A\s*\[?\s*[\d.]", head, re.M))
        except Exception:
            return False

    # any file under data/ (one level down too) whose content is a PFM list: JASPAR, HOCOMOCO, ...
    pool = [p for p in list(DATA.iterdir()) + list(DATA.glob("*/*")) if p.is_file()
            and p.suffix.lower() not in {".bed", ".zip", ".tsv", ".csv", ".png", ".pdf", ".py"}]
    files = sorted(p for p in pool if looks_like_pfm(p))
    if not files:
        listing = sorted(p.name for p in DATA.iterdir()) if DATA.exists() else ["(data/ folder not found)"]
        sys.exit(f"No motif (PFM) files found under {DATA}.\nFiles there: {listing}\n"
                 "Need a JASPAR-format PFM text file (starts with '>MA0001.1 NAME' then rows 'A [ ... ]'). "
                 "If your JASPAR file is a .zip, unzip it first.")
    mats = []
    per_dir = {}
    for f in files:
        got = parse_pfm_file(f)
        per_dir[f.parent.name] = per_dir.get(f.parent.name, 0) + len(got)
        mats += got
    for d, n in per_dir.items():
        print(f"Motifs: {n} matrices from folder '{d}'")
    out = []
    for name, counts in mats:
        if counts.shape[0] < 5:
            continue
        col = counts.sum(1, keepdims=True)
        col[col == 0] = 1.0
        p = (counts + PSEUDOCOUNT * 0.25) / (col + PSEUDOCOUNT)
        M = np.log2(p / 0.25)
        lo, hi = M.min(1).sum(), M.max(1).sum()
        out.append(dict(name=name, M=M, lo=lo, hi=hi, genes=tf_genes(name)))
    print(f"Total matrices scored: {len(out)}")
    return out


def matrix_deltas(seq, mat):
    """dPWM [L,4] for every substitution: best window (both strands) alt minus ref, 0 if no motif hit."""
    L, M = len(seq), mat["M"]
    w = M.shape[0]
    if L < w:
        return np.zeros((L, 4), np.float32)
    nW = L - w + 1
    M5 = np.hstack([M, np.zeros((w, 1))])
    idx = np.arange(nW)[:, None] + np.arange(w)[None, :]
    base = seq[idx]
    best_ref = np.full(L, -np.inf)
    best_alt = np.full((L, 4), -np.inf)
    for sm in (M5, M5[::-1][:, [3, 2, 1, 0, 4]]):
        contrib = sm[np.arange(w)[None, :], base]
        ref = contrib.sum(1)
        for j in range(w):
            sl = slice(j, j + nW)
            best_ref[sl] = np.maximum(best_ref[sl], ref)
            alt = ref[:, None] - contrib[:, j][:, None] + sm[j, :4][None, :]
            best_alt[sl] = np.maximum(best_alt[sl], alt)
    span = max(mat["hi"] - mat["lo"], 1e-9)
    rel = (np.maximum(best_alt, best_ref[:, None]) - mat["lo"]) / span
    d = best_alt - best_ref[:, None]
    d[rel < MIN_REL] = 0.0
    d[~np.isfinite(d)] = 0.0
    return d.astype(np.float32)


# ----------------------------- MPRA data ------------------------------------
def read_table(p):
    return pd.read_csv(p, sep=None, engine="python", compression="infer")


def load_mpra():
    files = sorted(list(MPRA_DIR.glob("*.tsv")) + list(MPRA_DIR.glob("*.tsv.gz")) +
                   list(MPRA_DIR.glob("*.csv")) + list(MPRA_DIR.glob("*.csv.gz")))
    if not files:
        sys.exit(f"No MPRA tables in {MPRA_DIR}.\nDownload the per-element GRCh38 files from "
                 "https://kircherlab.bihealth.org/satMutMPRA/ and save them there (e.g. F9.tsv, IRF4.tsv).")
    els = {}
    pairs = []
    for f in files:
        d0 = read_table(f)
        ecol = next((c for c in d0.columns if str(c).lower() == "element"), None)
        if ecol is not None:          # one combined file (e.g. GRCh38_ALL.tsv): split by its Element column
            pairs += [(str(n), g.drop(columns=[ecol])) for n, g in d0.groupby(ecol)]
        else:
            pairs.append((re.sub(r"\.(tsv|csv)(\.gz)?$", "", f.name, flags=re.I), d0))
    for name, d in pairs:
        if name.lower() in EXCLUDE_ELEMENTS:
            continue
        c = dict(chr=pick(d.columns, ["Chromosome", "chrom", "chr"], "chromosome"),
                 pos=pick(d.columns, ["Position", "pos"], "position"),
                 ref=pick(d.columns, ["Ref"], "ref"), alt=pick(d.columns, ["Alt"], "alt"),
                 val=pick(d.columns, ["Value", "effect", "log2FC"], "value"),
                 p=pick(d.columns, ["P-Value", "pvalue", "p"], "p-value"))
        d = d.rename(columns={v: k for k, v in c.items()})[list(c)]
        d["chr"] = "chr" + d["chr"].astype(str).str.replace("^chr", "", regex=True)
        d["ref"], d["alt"] = d["ref"].astype(str).str.upper(), d["alt"].astype(str).str.upper()
        d = d[d["ref"].isin(list(BASES)) & d["alt"].isin(list(BASES))].dropna(subset=["val", "p"])
        d["pos"] = d["pos"].astype(np.int64)
        d = d.drop_duplicates(["pos", "alt"])
        if len(d) < 20:
            print(f"  skip {name}: only {len(d)} usable variants")
            continue
        els[name] = d.reset_index(drop=True)
    print(f"MPRA elements loaded: {len(els)}; variants: {sum(len(d) for d in els.values()):,}")
    return els


def element_sequence(d):
    """Reconstruct the reference sequence from the Ref column (saturation mutagenesis covers every base)."""
    lo, hi = int(d["pos"].min()), int(d["pos"].max())
    seq = np.full(hi - lo + 1, 4, np.int8)
    ref = d.drop_duplicates("pos")
    seq[(ref["pos"] - lo).values] = [CODE[b] for b in ref["ref"]]
    return lo, seq


def group_elements(meta):
    names = list(meta)
    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if meta[a][0] == meta[b][0] and meta[a][1] <= meta[b][2] and meta[b][1] <= meta[a][2]:
                parent[find(a)] = find(b)
    return {n: find(n) for n in names}


# ----------------------------- expression -----------------------------------
def find_expr_file():
    if EXPR_FILE is not None:
        return Path(EXPR_FILE)
    c = [p for p in DATA.iterdir() if p.is_file() and p.suffix.lower() not in {".bed", ".py", ".png", ".pdf"}
         and re.search(r"ccle|tpm|cell[_ -]?line", p.name, re.I)]
    if not c:
        sys.exit("No cell-line expression table found in data/.\nSave a genes x cell lines TPM table there "
                 "(file name containing 'ccle', 'tpm' or 'cellline'), e.g. CCLE or ENCODE RNA-seq TPM.\n"
                 "Or set EXPR_FILE at the top of this script.")
    best = max(c, key=lambda p: p.stat().st_size)
    print(f"Using expression table: {best.name}")
    return best


def load_expression(tf_gene_set):
    raw = pd.read_csv(find_expr_file(), sep=None, engine="python", compression="infer")
    lc = {c.lower().replace("_", " "): c for c in raw.columns}
    if "cell line" in lc:                      # long format, e.g. Human Protein Atlas rna_celline.tsv
        gcol = lc.get("gene name") or lc.get("gene") or raw.columns[0]
        vcol = next((lc[k] for k in ("tpm", "ptpm", "ntpm") if k in lc), None)
        if vcol is None:
            sys.exit(f"Long-format expression table has no TPM column. Columns: {list(raw.columns)}")
        print(f"Expression table is long format; using '{gcol}' x '{lc['cell line']}' values from '{vcol}'")
        ex = raw.pivot_table(index=gcol, columns=lc["cell line"], values=vcol, aggfunc="max")
    else:
        ex = raw.set_index(raw.columns[0])
    ex.index = ex.index.astype(str)
    if ex.index[0].startswith("ENSG"):
        gm_path = GENE_MAP_FILE or next(iter(sorted(DATA.glob("*probemap*"))), None)
        if gm_path is None:
            sys.exit("Expression rows are Ensembl ids but no data/*probemap* file was found.")
        gm = pd.read_csv(gm_path, sep="\t").iloc[:, :2]
        gm.columns = ["id", "gene"]
        gm["id"] = gm["id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        gm = gm.drop_duplicates("id").set_index("id")["gene"]
        ex.index = pd.Index(ex.index.str.replace(r"\.\d+$", "", regex=True)).map(gm)
        ex = ex[ex.index.notna()]
    # Rows should be genes. If the columns look like genes instead, transpose.
    if len(set(ex.index.astype(str).str.upper()) & tf_gene_set) < len(set(ex.columns.astype(str).str.upper()) & tf_gene_set):
        ex = ex.T
    ex.index = ex.index.astype(str).str.upper()
    ex = ex.apply(pd.to_numeric, errors="coerce")
    ex = ex.groupby(level=0).max()
    if np.nanmax(ex.values) < 25:              # values look like log2(TPM+1) (DepMap style): back-transform
        print("Expression values look log2(TPM+1); converting to TPM.")
        ex = 2.0 ** ex - 1.0
    print(f"Expression: {ex.shape[0]:,} genes x {ex.shape[1]} cell lines; "
          f"{len(set(ex.index) & tf_gene_set)} motif TFs present")
    return ex


def match_cell_column(aliases, columns):
    if not aliases:
        return None
    al = {norm_name(a) for a in aliases}
    for c in columns:
        if norm_name(c) in al or norm_name(str(c).split("_")[0]) in al:
            return c
    return None


# ----------------------------- cCRE overlap ---------------------------------
def ccre_overlap(chrom, pos):
    out = {}
    for cls, path in CCRE_FILES.items():
        if not Path(path).exists():
            print(f"  cCRE file missing: {path}")
            continue
        b = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
        hit = np.zeros(len(chrom), bool)
        for c, g in b.groupby("c"):
            g = g.sort_values("s")
            sel = np.where(chrom == c)[0]
            if not len(sel):
                continue
            s, e = g["s"].values, g["e"].values
            i = np.searchsorted(s, pos[sel] - 1, side="right") - 1
            ok = i >= 0
            hit[sel[ok]] = (pos[sel][ok] - 1) < e[i[ok]]
        out[cls] = hit
    return out


# ----------------------------- LOEO engine ----------------------------------
def group_stats(F, z, groups, ug):
    st = {}
    for g in ug:
        r = groups == g
        X = F[r]
        nz = (X != 0)
        y = z[r]
        st[g] = dict(n=nz.sum(0).astype(float), Sx=X.sum(0, dtype=np.float64),
                     Sxx=(X.astype(np.float64) ** 2).sum(0), Sy=nz.T.astype(np.float64) @ y,
                     Syy=nz.T.astype(np.float64) @ (y * y), Sxy=X.T.astype(np.float64) @ y)
    tot = {k: sum(s[k] for s in st.values()) for k in next(iter(st.values()))}
    return st, tot


def weights_from(stats_):
    n, Sx, Sxx, Sy, Syy, Sxy = (stats_[k] for k in ("n", "Sx", "Sxx", "Sy", "Syy", "Sxy"))
    with np.errstate(all="ignore"):
        r = (n * Sxy - Sx * Sy) / np.sqrt((n * Sxx - Sx ** 2) * (n * Syy - Sy ** 2))
    r[~np.isfinite(r)] = 0.0
    r[n < MIN_TRAIN_N] = 0.0
    return r


def loeo_scores(F, z, groups, usable):
    """Signed motif score S for each variant, with weights learned excluding the variant's own group."""
    S = np.full(len(F), np.nan)
    ug = [g for g in pd.unique(groups[usable])]
    if len(ug) < 2:
        return S
    st, tot = group_stats(F[usable], z[usable], groups[usable], ug)
    for g in ug:
        te = (groups == g) & usable
        if te.sum() < MIN_TEST:
            continue
        r = weights_from({k: tot[k] - st[g][k] for k in tot})
        sel = np.argsort(-np.abs(r))[:TOP_K]
        sel = sel[r[sel] != 0]
        if len(sel):
            S[te] = F[te][:, sel].astype(np.float64) @ r[sel]
    return S


def evaluate(S, V, label):
    rows = []
    for el, d in V.groupby("element"):
        m = np.isfinite(S[d.index])
        if m.sum() < MIN_TEST:
            continue
        s, v = S[d.index][m], d["val"].values[m]
        rho = stats.spearmanr(s, v)[0] if np.std(s) > 0 else np.nan
        called = (s != 0) & (d["p"].values[m] < SIG_P)
        agree = int((np.sign(s[called]) == np.sign(v[called])).sum())
        rows.append(dict(cond=label, element=el, group=d["group"].iloc[0], n=int(m.sum()), rho=rho,
                         n_called_sig=int(called.sum()), n_sign_agree=agree))
    return pd.DataFrame(rows)


def summarise(tab, rng):
    t = tab.dropna(subset=["rho"])
    if t.empty:
        return dict(mean_rho=np.nan, ci_lo=np.nan, ci_hi=np.nan, sign_conc=np.nan, sign_n=0, sign_p=np.nan)
    gm = t.groupby("group").agg(rho=("rho", "mean"), a=("n_sign_agree", "sum"), n=("n_called_sig", "sum"))
    boots = []
    for _ in range(N_BOOT):
        pick_ = rng.integers(0, len(gm), len(gm))
        b = gm.iloc[pick_]
        boots.append((b["rho"].mean(), b["a"].sum() / b["n"].sum() if b["n"].sum() else np.nan))
    boots = np.array(boots, float)
    a, n = int(gm["a"].sum()), int(gm["n"].sum())
    return dict(mean_rho=gm["rho"].mean(), ci_lo=np.nanpercentile(boots[:, 0], 2.5), ci_hi=np.nanpercentile(boots[:, 0], 97.5),
                sign_conc=a / n if n else np.nan, sign_n=n,
                sign_p=stats.binomtest(a, n, 0.5).pvalue if n else np.nan,
                sign_ci_lo=np.nanpercentile(boots[:, 1], 2.5), sign_ci_hi=np.nanpercentile(boots[:, 1], 97.5))


def lr_plus(call, active, groups, rng):
    def one(c, a):
        tp, fn = (c & a).sum(), (~c & a).sum()
        fp, tn = (c & ~a).sum(), (~c & ~a).sum()
        sens = tp / (tp + fn) if tp + fn else np.nan
        fpr = fp / (fp + tn) if fp + tn else np.nan
        return sens / fpr if fpr and np.isfinite(fpr) and fpr > 0 else np.nan, sens, fpr
    est = one(call, active)
    ug = pd.unique(groups)
    idx = {g: np.where(groups == g)[0] for g in ug}
    bs = []
    for _ in range(N_BOOT):
        ii = np.concatenate([idx[g] for g in rng.choice(ug, len(ug))])
        bs.append(one(call[ii], active[ii])[0])
    bs = np.array(bs, float)
    return est, np.nanpercentile(bs, 2.5), np.nanpercentile(bs, 97.5)


# ----------------------------- main -----------------------------------------
def main():
    rng = np.random.default_rng(SEED)
    RESULTS.mkdir(exist_ok=True)
    mats = load_motifs()
    nM = len(mats)
    tf_set = {g for m in mats for g in m["genes"]}
    els = load_mpra()
    ex = load_expression(tf_set)

    # per-element features
    frames, feats, meta, masks = [], [], {}, {}
    for name, d in els.items():
        key = name.lower()
        aliases = ELEMENT_CELL.get(key)
        if aliases is None and re.sub(r"[-.][^-.]*$", "", key) in ELEMENT_CELL:
            aliases = ELEMENT_CELL[re.sub(r"[-.][^-.]*$", "", key)]
        col = match_cell_column(aliases, ex.columns)
        lo, seq = element_sequence(d)
        meta[name] = (d["chr"].iloc[0], int(d["pos"].min()), int(d["pos"].max()))
        if col is None:
            print(f"  {name}: no expression column for {aliases} -> excluded from expression-filtered analysis")
            masks[name] = None
        else:
            expd = set(ex.index[ex[col] >= EXPR_TPM])
            masks[name] = np.array([all(g in expd for g in m["genes"]) for m in mats])
            print(f"  {name}: cell line column '{col}', {masks[name].sum()} of {nM} motifs expressed")
        D = np.stack([matrix_deltas(seq, m) for m in mats], axis=-1)       # [L,4,nM]
        pi = (d["pos"].values - lo)
        ai = np.array([CODE[a] for a in d["alt"]])
        feats.append(D[pi, ai, :])
        d = d.copy()
        d["element"] = name
        frames.append(d)
    V = pd.concat(frames, ignore_index=True)
    F = np.vstack(feats).astype(np.float32)
    gmap = group_elements(meta)
    V["group"] = V["element"].map(gmap)
    V["q"] = bh(V["p"].values)
    V["active"] = V["q"] < ACTIVE_Q
    # within-element robust z of the effect (winsorised) used only for training weights
    z = np.zeros(len(V))
    for el, d in V.groupby("element"):
        v = d["val"].values
        mad = 1.4826 * np.median(np.abs(v - np.median(v)))
        s = mad if mad > 0 else (v.std() or 1.0)
        z[d.index] = np.clip((v - np.median(v)) / s, -4, 4)
    groups = V["group"].values
    elem_mask = np.zeros((len(V), nM), bool)
    usable = np.zeros(len(V), bool)
    for name in els:
        r = (V["element"] == name).values
        if masks[name] is not None:
            elem_mask[r] = masks[name]
            usable[r] = True
    if usable.sum() < MIN_TEST:
        sys.exit("No element has a matching expression column, so nothing can be tested. "
                 "Check ELEMENT_CELL aliases against your expression table's column names.")
    print(f"\nVariants: {len(V):,} total, {int(usable.sum()):,} in elements with matched expression, "
          f"{len(pd.unique(groups[usable]))} independent sequence groups")

    # ---- conditions ----
    F_expr = F * elem_mask
    S_expr = loeo_scores(F_expr, z, groups, usable)
    S_all = loeo_scores(F, z, groups, usable)
    S_naive = np.where(usable, F_expr.sum(1), np.nan)       # unsupervised: every expressed TF assumed an activator

    tabs = {"expressed_LOEO": evaluate(S_expr, V, "expressed_LOEO"),
            "all_TFs_LOEO": evaluate(S_all, V, "all_TFs_LOEO"),
            "expressed_unsupervised_sum": evaluate(S_naive, V, "expressed_unsupervised_sum")}
    summ = {k: summarise(t, rng) for k, t in tabs.items()}

    # permutation null for the expression filter (same number of expressed TFs per element)
    null_rho, null_sc = [], []
    for _ in range(N_PERM):
        pm = np.zeros_like(elem_mask)
        for name in els:
            if masks[name] is None:
                continue
            r = (V["element"] == name).values
            pm[r] = rng.permutation(masks[name])
        Sp = loeo_scores(F * pm, z, groups, usable)
        sp = summarise(evaluate(Sp, V, "perm"), rng)
        null_rho.append(sp["mean_rho"])
        null_sc.append(sp["sign_conc"])
    null_rho, null_sc = np.array(null_rho, float), np.array(null_sc, float)
    obs = summ["expressed_LOEO"]
    p_rho = (1 + np.nansum(null_rho >= obs["mean_rho"])) / (1 + np.sum(np.isfinite(null_rho)))
    p_sc = (1 + np.nansum(null_sc >= obs["sign_conc"])) / (1 + np.sum(np.isfinite(null_sc)))

    # ---- cCRE prior + LR+ ----
    cc = ccre_overlap(V["chr"].values, V["pos"].values)
    in_ccre = np.zeros(len(V), bool)
    for h in cc.values():
        in_ccre |= h
    U = usable & np.isfinite(S_expr)
    Vu, act = V[U].reset_index(drop=True), V.loc[U, "active"].values
    gu, Su = groups[U], S_expr[U]
    call_reg = in_ccre[U]
    thr = Vu.assign(a=np.abs(Su)).groupby("element")["a"].transform(lambda x: x.quantile(1 - CALL_FRAC)).values
    call_motif = (np.abs(Su) >= thr) & (Su != 0)
    val_sign = np.sign(Vu["val"].values)
    # directed call: a hit only if the predicted direction matches the observed direction of an active variant
    lr_rows = []

    def add(label, call, active_def, note=""):
        (lr, sens, fpr), lo, hi = lr_plus(call, active_def, gu, rng)
        lr_rows.append(dict(call=label, call_rate=call.mean(), sensitivity=sens, FPR=fpr, LR_plus=lr,
                            CI_lo=lo, CI_hi=hi, note=note))

    add("registry: in any cCRE", call_reg, act)
    for cls, h in cc.items():
        add(f"registry: {cls}", h[U], act)
    add(f"motif LOEO: top {int(CALL_FRAC*100)}% |S|", call_motif, act)
    add("registry AND motif", call_reg & call_motif, act)
    # a true positive needs the predicted sign to match the observed sign; wrong-sign hits on active variants count as misses
    call_dir = call_motif & ((np.sign(Su) == val_sign) | ~act)
    add("motif LOEO, direction-matched", call_dir, act,
        "true positive = active AND predicted sign equals observed sign")

    lr = pd.DataFrame(lr_rows)
    frac_ccre = float(call_reg.mean())

    # ---- outputs ----
    by_el = pd.concat(tabs.values(), ignore_index=True)
    by_el.to_csv(RESULTS / "26_mpra_loeo_by_element.csv", index=False)
    lr.to_csv(RESULTS / "26_mpra_lr_plus.csv", index=False)
    srows = [dict(cond=k, **v) for k, v in summ.items()]
    srows.append(dict(cond="permuted_expression_null", mean_rho=np.nanmean(null_rho), sign_conc=np.nanmean(null_sc),
                      note=f"perm p(rho)={p_rho:.3f}, perm p(sign)={p_sc:.3f}, n_perm={N_PERM}"))
    pd.DataFrame(srows).to_csv(RESULTS / "26_mpra_summary.csv", index=False)
    r_all = weights_from(group_stats(F_expr[usable], z[usable], groups[usable], list(pd.unique(groups[usable])))[1])
    w = pd.DataFrame(dict(motif=[m["name"] for m in mats], r_in_sample=r_all)).assign(abs_r=lambda t: t.r_in_sample.abs())
    w.sort_values("abs_r", ascending=False).head(100).to_csv(RESULTS / "26_mpra_top_motifs_insample.csv", index=False)

    pd.set_option("display.width", 200)
    print("\n=== Direction-aware motif scoring, leave-one-element-out (group-level 95% bootstrap CI) ===")
    for k, v in summ.items():
        print(f"{k:30s} mean Spearman rho = {v['mean_rho']:+.3f} [{v['ci_lo']:+.3f},{v['ci_hi']:+.3f}]   "
              f"sign agreement = {v['sign_conc']:.3f} (n={v['sign_n']}, binomial p={v['sign_p']:.2g})")
    print(f"permuted-expression null: rho {np.nanmean(null_rho):+.3f}, sign {np.nanmean(null_sc):.3f}; "
          f"observed vs null: p(rho)={p_rho:.3f}, p(sign)={p_sc:.3f}")
    print("\n=== LR+ for predicting MPRA-active variants (q<%.2f) ===" % ACTIVE_Q)
    print(lr.drop(columns="note").round(3).to_string(index=False))
    print(f"\nShare of tested variants inside any cCRE: {frac_ccre:.3f}")
    if frac_ccre > 0.9 or frac_ccre < 0.1:
        print("WARNING: the cCRE call is nearly constant across these variants (the MPRA elements were chosen as "
              "regulatory regions), so the variant-level LR+ for the registry is close to 1 by construction and "
              "cannot show whether the registry is limiting.")
    print(f"\nFiles written to {RESULTS}: 26_mpra_summary.csv, 26_mpra_loeo_by_element.csv, "
          "26_mpra_lr_plus.csv, 26_mpra_top_motifs_insample.csv")


if __name__ == "__main__":
    main()