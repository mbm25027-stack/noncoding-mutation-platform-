#!/usr/bin/env python3
"""
nc_junction_control.py - distance-to-exon-junction control for the FROZEN
NC-MutScan v3.3 engine. No engine parameter is touched: it only reads the
cached scores.tsv written by heldout_full_pool_eval.py.

Question answered: does the engine add anything beyond "how close is the
variant to an exon-intron junction"?

Steps
  1. Distance d from each variant to the nearest annotated splice junction
     (GTF). Convention: intronic side, first intronic base = 1 (so +/-1,2 are
     the canonical splice positions); exonic side, last/first exonic base = 1.
  2. Strata: d 1-2, 3-20, 21-50, >50 bp (all, intronic only, exonic only).
  3. Distance-matched benign controls: benign variants are re-weighted so
     their distance distribution (fine distance bins x exonic/intronic side)
     equals the pathogenic one. A proximity-only scorer is then at chance.
  4. Proximity vs proximity + engine features:
       - likelihood-ratio test of nested logistic models (in-sample)
       - gene-grouped K-fold cross-fitted AUROC / AUPRC and delta, with a
         gene-clustered bootstrap CI
       - sensitivity of a proximity-only scorer at the engine's own FPR
  Dependencies: numpy + standard library only.

Pools: the +/-1,2 stratum only exists in a pool that keeps canonical splice
variants, so use the splice-inclusive pool for the full analysis (default
--exclude-mc) and the clean pool (add the two splice terms to --exclude-mc)
for the strictly noncoding analysis.

Example (splice-inclusive pool, chr12-22):
  python nc_junction_control.py --pool check_blind\\pool_all.vcf ^
     --scores eval_full_pool_with_splice\\scores.tsv --gtf gencode.gtf.gz ^
     --out eval_junction_splice
"""
import argparse, gzip, json, math, os, sys, time
from collections import Counter, defaultdict
import numpy as np

BIG = 10 ** 9
PATH_SIG = {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"}
BEN_SIG = {"benign", "likely_benign", "benign/likely_benign"}
COARSE = ["d1-2", "d3-20", "d21-50", "d>50"]
LIMITS = np.array([1, 2, 3, 4, 5, 8, 12, 20, 30, 50, 100, 300, 1000, 5000])
NBIN = len(LIMITS) + 1
NONE_CELL = 2 * NBIN
NCELL = NONE_CELL + 1


# ------------------------------------------------------------------ inputs
def load_pool(path, exclude_mc):
    rows, skipped, excl = [], Counter(), set(exclude_mc)
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.rstrip("\r\n").split("\t")
            info = dict(kv.split("=", 1) for kv in t[7].split(";") if "=" in kv)
            sig = info.get("CLNSIG", "").lower()
            if sig in PATH_SIG:
                y = 1
            elif sig in BEN_SIG:
                y = 0
            else:
                skipped["not_clean_label"] += 1
                continue
            terms = {x.split("|")[-1] for x in info.get("MC", "").split(",") if x}
            if terms & excl:
                skipped["excluded_mc"] += 1
                continue
            gi = info.get("GENEINFO", "")
            gene = gi.split("|")[0].split(":")[0] if gi else ""
            rows.append((t[0].replace("chr", ""), int(t[1]), t[3].upper(), t[4].upper(), y, gene))
    return rows, skipped


def load_scores(path):
    d = {}
    with open(path) as f:
        next(f, None)
        for line in f:
            k, s, c, st = line.rstrip("\n").split("\t")
            if st == "ok":
                d[k] = (float(s), int(c))
    return d


def load_cuts(path, tag=""):
    """Splice-junction 'cuts' per contig. Each intron [i1,i2] gives a cut at
    i1 (type 0: intron to the right) and at i2+1 (type 1: intron to the left)."""
    op = gzip.open if path.endswith(".gz") else open
    tx = defaultdict(list)
    with op(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line[0] == "#":
                continue
            t = line.split("\t")
            if len(t) < 9 or t[2] != "exon":
                continue
            a = t[8]
            if tag and tag not in a:
                continue
            i = a.find('transcript_id "')
            if i < 0:
                continue
            j = a.find('"', i + 15)
            tx[(t[0].replace("chr", ""), a[i + 15:j])].append((int(t[3]), int(t[4])))
    introns = defaultdict(set)
    for (ch, _), ex in tx.items():
        if len(ex) < 2:
            continue
        ex.sort()
        for (s1, e1), (s2, e2) in zip(ex[:-1], ex[1:]):
            if s2 - e1 > 1:
                introns[ch].add((e1 + 1, s2 - 1))
    cuts = {}
    for ch, S in introns.items():
        arr = sorted({(i1, 0) for i1, i2 in S} | {(i2 + 1, 1) for i1, i2 in S})
        cuts[ch] = (np.array([a for a, _ in arr], dtype=np.int64),
                    np.array([b for _, b in arr], dtype=np.int8))
    return cuts, len(tx), sum(len(v) for v in introns.values())


# --------------------------------------------------------------- distances
def nearest_cut(ch_arr, p, cuts):
    d = np.full(len(p), BIG, dtype=np.int64)
    side = np.full(len(p), -1, dtype=np.int8)       # 1 intronic, 0 exonic, -1 none
    for ch in np.unique(ch_arr):
        if ch not in cuts:
            continue
        m = ch_arr == ch
        c, t = cuts[ch]
        pp = p[m]
        idx = np.searchsorted(c, pp)
        bd = np.full(pp.shape, BIG, dtype=np.int64)
        bs = np.full(pp.shape, -1, dtype=np.int8)
        for off in (-2, -1, 0, 1):
            k = np.clip(idx + off, 0, len(c) - 1)
            cc, tt = c[k], t[k]
            isL = tt == 0
            intr = np.where(isL, pp >= cc, pp < cc)
            dd = np.where(isL, np.where(pp >= cc, pp - cc + 1, cc - pp),
                          np.where(pp < cc, cc - pp, pp - cc + 1))
            better = dd < bd
            bd = np.where(better, dd, bd)
            bs = np.where(better, intr.astype(np.int8), bs)
        d[m], side[m] = bd, bs
    return d, side


def junction_distance(chrom, pos, reflen, cuts):
    d1, s1 = nearest_cut(chrom, pos, cuts)
    d2, s2 = nearest_cut(chrom, pos + reflen - 1, cuts)
    use2 = d2 < d1
    return np.where(use2, d2, d1), np.where(use2, s2, s1)


def cell_ids(d, side):
    b = np.searchsorted(LIMITS, d, side="left")
    cell = side.astype(np.int64) * NBIN + b
    cell[side < 0] = NONE_CELL
    return cell


def coarse_bin(d):
    c = np.full(len(d), 3)
    c[d <= 50] = 2
    c[d <= 20] = 1
    c[d <= 2] = 0
    return c


# ----------------------------------------------------------------- models
def sigmoid(z):
    return 1.0 / (1.0 + np.exp(-np.clip(z, -35, 35)))


def loglik(X, y, b):
    eta = X @ b
    return float(np.sum(y * eta - np.logaddexp(0, eta)))


def fit_logit(X, y, ridge=1e-3, iters=100):
    """Ridge-stabilised logistic regression by DAMPED Newton (step halving on the
    penalised log-likelihood). Plain Newton diverges when cells differ a lot in
    prevalence, so the line search is required."""
    n, k = X.shape
    pen = np.full(k, ridge)
    pen[0] = 0.0
    b = np.zeros(k)
    m = min(max(y.mean(), 1e-6), 1 - 1e-6)
    b[0] = math.log(m / (1 - m))
    obj = lambda bb: loglik(X, y, bb) - 0.5 * float(np.sum(pen * bb * bb))
    cur = obj(b)
    for _ in range(iters):
        mu = sigmoid(X @ b)
        W = np.maximum(mu * (1 - mu), 1e-9)
        g = X.T @ (y - mu) - pen * b
        H = (X * W[:, None]).T @ X + np.diag(pen)
        step = np.linalg.solve(H, g)
        t, ok = 1.0, False
        for _h in range(40):
            new = obj(b + t * step)
            if np.isfinite(new) and new >= cur - 1e-9:
                ok = True
                break
            t *= 0.5
        if not ok:
            break
        b = b + t * step
        gain, cur = new - cur, new
        if np.max(np.abs(t * step)) < 1e-7 or abs(gain) < 1e-9:
            break
    return b


def chi2_sf(x, k):
    if x <= 0:
        return 1.0
    a, xx = k / 2.0, x / 2.0
    if xx < a + 1:
        term = 1.0 / a
        s, n = term, a
        for _ in range(2000):
            n += 1
            term *= xx / n
            s += term
            if abs(term) < abs(s) * 1e-15:
                break
        return max(0.0, 1.0 - s * math.exp(-xx + a * math.log(xx) - math.lgamma(a)))
    tiny = 1e-300
    b = xx + 1 - a
    c, d = 1 / tiny, 1 / b
    h = d
    for i in range(1, 2000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-15:
            break
    return h * math.exp(-xx + a * math.log(xx) - math.lgamma(a))


def build_X(cell, F, cells_used, ref):
    X0 = np.column_stack([np.ones(len(cell))] + [(cell == c).astype(float) for c in cells_used if c != ref])
    return X0, np.column_stack([X0, F])


def engine_features(score, call, sd):
    s = np.maximum(score, 0.0)
    return np.column_stack([(score > 0).astype(float), np.log1p(s) / sd, call.astype(float)])


def lrt(cell, F, y):
    cells_used = np.unique(cell)
    ref = cells_used[np.argmax([np.sum(cell == c) for c in cells_used])]
    X0, X1 = build_X(cell, F, cells_used, ref)
    b0, b1 = fit_logit(X0, y), fit_logit(X1, y)
    stat = 2 * (loglik(X1, y, b1) - loglik(X0, y, b0))
    df = F.shape[1]
    return stat, df, chi2_sf(stat, df), b1[-df:]


# ---------------------------------------------------------------- metrics
class Ranked:
    def __init__(self, score, y):
        o = np.argsort(-score, kind="stable")
        self.o = o
        s = score[o]
        self.y = y[o].astype(float)
        self.ends = np.r_[np.nonzero(np.diff(s))[0], len(s) - 1]

    def _curve(self, w):
        w = w[self.o]
        return np.cumsum(w * self.y)[self.ends], np.cumsum(w * (1 - self.y))[self.ends]

    def auc_ap(self, w):
        ctp, cfp = self._curve(w)
        P, N = ctp[-1], cfp[-1]
        if P <= 0 or N <= 0:
            return float("nan"), float("nan")
        rec, den = ctp / P, ctp + cfp
        prec = np.divide(ctp, den, out=np.zeros_like(den), where=den > 0)
        ap = float(np.sum(np.diff(np.r_[0.0, rec]) * prec))
        tpr, fpr = np.r_[0.0, rec], np.r_[0.0, cfp / N]
        return float(np.sum(np.diff(fpr) * (tpr[1:] + tpr[:-1]) / 2)), ap

    def sens_at_fpr(self, w, f):
        ctp, cfp = self._curve(w)
        P, N = ctp[-1], cfp[-1]
        if P <= 0 or N <= 0 or not np.isfinite(f):
            return float("nan")
        return float(np.interp(f, np.r_[0.0, cfp / N], np.r_[0.0, ctp / P]))


def call_metrics(w, y, call):
    tp, fp = np.sum(w * y * call), np.sum(w * (1 - y) * call)
    P, N = np.sum(w * y), np.sum(w * (1 - y))
    nan = float("nan")
    return dict(sens=tp / P if P > 0 else nan, fpr=fp / N if N > 0 else nan,
                ppv=tp / (tp + fp) if tp + fp > 0 else nan)


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--scores", required=True, help="scores.tsv for the same pool (heldout_full_pool_eval.py output)")
    ap.add_argument("--gtf", required=True, help="GENCODE-style GTF (.gtf or .gtf.gz), GRCh38")
    ap.add_argument("--out", default="eval_junction_control")
    ap.add_argument("--exclude-mc", default="nonsense,inframe_indel,no_sequence_alteration")
    ap.add_argument("--tag", default="", help='keep only transcripts whose attributes contain this text, e.g. MANE_Select')
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=12345)
    a = ap.parse_args()
    t0 = time.time()

    rows, skipped = load_pool(a.pool, [x for x in a.exclude_mc.split(",") if x])
    scores = load_scores(a.scores)
    keep, sc, cl, missing = [], [], [], 0
    for r in rows:
        s = scores.get(f"{r[0]}:{r[1]}:{r[2]}:{r[3]}")
        if s is None:
            missing += 1
            continue
        keep.append(r)
        sc.append(s[0])
        cl.append(s[1])
    chrom = np.array([r[0] for r in keep])
    pos = np.array([r[1] for r in keep], dtype=np.int64)
    reflen = np.array([len(r[2]) for r in keep], dtype=np.int64)
    y = np.array([r[4] for r in keep], dtype=float)
    score, call = np.array(sc), np.array(cl, dtype=float)
    n = len(y)
    print(f"Pool: {n:,} variants | P {int(y.sum()):,} | B {int(n - y.sum()):,} | "
          f"skipped {dict(skipped)} | no score {missing:,}")
    if missing:
        print("WARNING: some pool variants have no ok score in --scores; use the scores.tsv of the same pool.")

    cuts, ntx, nintron = load_cuts(a.gtf, a.tag)
    print(f"GTF: {ntx:,} transcripts, {nintron:,} unique introns" + (f" (tag filter '{a.tag}')" if a.tag else ""))
    d, side = junction_distance(chrom, pos, reflen, cuts)
    cell, cb = cell_ids(d, side), coarse_bin(d)
    print(f"Variants with no junction on their contig: {int(np.sum(side < 0)):,}")

    names = [r[5] if r[5] else f"nogene_{r[0]}_{r[1] // 1_000_000}" for r in keep]
    uniq = {}
    cidx = np.array([uniq.setdefault(nm, len(uniq)) for nm in names])
    nc = len(uniq)

    # --- distance-matched benign weights
    cP = np.bincount(cell[y == 1], minlength=NCELL).astype(float)
    cB = np.bincount(cell[y == 0], minlength=NCELL).astype(float)
    ratio = np.divide(cP, cB, out=np.zeros_like(cP), where=cB > 0)
    wm = np.where(y == 1, (cB[cell] > 0).astype(float), ratio[cell])
    p_kept = float(np.sum(wm[y == 1]))
    print(f"Distance matching: {int(p_kept):,} of {int(y.sum()):,} pathogenic variants have benign controls in their "
          f"distance cell; matched benign weight {np.sum(wm[y == 0]):,.0f}")

    # --- models
    sd = float(np.std(np.log1p(np.maximum(score, 0.0)))) or 1.0
    F = engine_features(score, call, sd)
    cells_used = np.unique(cell)
    ref = cells_used[np.argmax([np.sum(cell == c) for c in cells_used])]
    X0, X1 = build_X(cell, F, cells_used, ref)
    stat, df, pval, coefs = lrt(cell, F, y)
    print(f"LRT proximity vs proximity+engine: chi2={stat:.2f}, df={df}, p={pval:.3g}")

    rng = np.random.default_rng(a.seed)
    fold = (rng.permutation(nc) % a.folds)[cidx]
    eta0, eta1 = np.zeros(n), np.zeros(n)
    for k in range(a.folds):
        tr, te = fold != k, fold == k
        eta0[te] = X0[te] @ fit_logit(X0[tr], y[tr])
        eta1[te] = X1[te] @ fit_logit(X1[tr], y[tr])
    print(f"Cross-fitted predictions done ({a.folds} gene-grouped folds)")

    R_eng, R_prox = Ranked(score, y), Ranked(-d.astype(float), y)
    R_m0, R_m1 = Ranked(eta0, y), Ranked(eta1, y)

    groups = {}
    for ci, nm in enumerate(COARSE):
        for lab, m in (("all", cb == ci), ("intronic", (cb == ci) & (side == 1)), ("exonic", (cb == ci) & (side == 0))):
            idx = np.nonzero(m)[0]
            if len(idx) and 0 < y[idx].sum() < len(idx):
                groups[f"{nm} {lab}"] = (idx, Ranked(score[idx], y[idx]))

    def evaluate(cw):
        out = {}
        for setting, wmv in (("unmatched", None), ("matched", wm)):
            w = cw if wmv is None else cw * wmv
            cm = call_metrics(w, y, call)
            ea, ep = R_eng.auc_ap(w)
            pa, pp = R_prox.auc_ap(w)
            ps = R_prox.sens_at_fpr(w, cm["fpr"])
            a0, p0 = R_m0.auc_ap(w)
            a1, p1 = R_m1.auc_ap(w)
            for k, v in (("engine_call_sens", cm["sens"]), ("engine_call_fpr", cm["fpr"]), ("engine_call_ppv", cm["ppv"]),
                         ("engine_auroc", ea), ("engine_auprc", ep), ("proximity_auroc", pa), ("proximity_auprc", pp),
                         ("proximity_sens_at_engine_fpr", ps), ("sens_engine_minus_proximity", cm["sens"] - ps),
                         ("model_prox_auroc_oof", a0), ("model_prox_plus_engine_auroc_oof", a1), ("delta_auroc", a1 - a0),
                         ("model_prox_auprc_oof", p0), ("model_prox_plus_engine_auprc_oof", p1), ("delta_auprc", p1 - p0)):
                out[f"{setting}|{k}"] = v
            for gname, (idx, Rg) in groups.items():
                wg = w[idx]
                cg = call_metrics(wg, y[idx], call[idx])
                out[f"{setting}|G|{gname}|sens"] = cg["sens"]
                out[f"{setting}|G|{gname}|fpr"] = cg["fpr"]
                out[f"{setting}|G|{gname}|ppv"] = cg["ppv"]
                out[f"{setting}|G|{gname}|auroc"] = Rg.auc_ap(wg)[0]
        return out

    point = evaluate(np.ones(n))
    draws = defaultdict(list)
    for b in range(a.n_boot):
        cnt = np.bincount(rng.integers(0, nc, nc), minlength=nc)
        for k, v in evaluate(cnt[cidx].astype(float)).items():
            draws[k].append(v)
        if (b + 1) % 200 == 0:
            print(f"  bootstrap {b + 1}/{a.n_boot}  ({time.time() - t0:,.0f}s)")
    ci = {k: [float(x) for x in np.nanpercentile(v, [2.5, 97.5])] for k, v in draws.items()} if a.n_boot else {}
    frac_pos = {k: float(np.nanmean(np.array(draws[k]) > 0)) for k in ("matched|delta_auroc", "unmatched|delta_auroc",
                                                                       "matched|delta_auprc", "unmatched|delta_auprc")} if a.n_boot else {}

    # --- stratified LRT (in-sample)
    strat_lrt = {}
    for ci_, nm in enumerate(COARSE):
        idx = np.nonzero(cb == ci_)[0]
        if len(idx) and y[idx].sum() >= 30 and (1 - y[idx]).sum() >= 30:
            st, dff, pv, co = lrt(cell[idx], F[idx], y[idx])
            strat_lrt[nm] = dict(chi2=st, df=dff, p=pv, n_pathogenic=int(y[idx].sum()), n_benign=int((1 - y[idx]).sum()))

    # --- write outputs
    os.makedirs(a.out, exist_ok=True)
    rep = dict(pool=a.pool, scores=a.scores, gtf=a.gtf, tag=a.tag, exclude_mc=a.exclude_mc, n=n,
               n_pathogenic=int(y.sum()), n_benign=int(n - y.sum()), n_clusters=nc, n_boot=a.n_boot,
               folds=a.folds, seed=a.seed, distance_matching_pathogenic_kept=p_kept,
               lrt_overall=dict(chi2=stat, df=df, p=pval,
                                coef_log_odds=dict(zip(["S_raw>0", "log1p(S_raw)/sd", "frozen_call"], map(float, coefs)))),
               lrt_by_distance_bin=strat_lrt,
               metrics={k: dict(estimate=point[k], ci95=ci.get(k)) for k in point},
               fraction_bootstrap_draws_delta_positive=frac_pos)
    json.dump(rep, open(os.path.join(a.out, "report.json"), "w"), indent=2)

    with open(os.path.join(a.out, "variant_distance.tsv"), "w") as o:
        o.write("key\tlabel\tgene\tside\tdist\tcoarse_bin\tcell\tS_raw\tcall\tmatched_weight\n")
        for i, r in enumerate(keep):
            o.write(f"{r[0]}:{r[1]}:{r[2]}:{r[3]}\t{int(y[i])}\t{r[5]}\t{['exonic', 'intronic', 'none'][int(side[i])] if side[i] >= 0 else 'none'}\t"
                    f"{int(d[i]) if d[i] < BIG else ''}\t{COARSE[cb[i]]}\t{cell[i]}\t{score[i]}\t{int(call[i])}\t{wm[i]:.6g}\n")

    def f(k, pct=True):
        v, c = point[k], ci.get(k)
        s = f"{100 * v:.2f}%" if pct else f"{v:.4f}"
        if c:
            s += f"  [{100 * c[0]:.2f}%, {100 * c[1]:.2f}%]" if pct else f"  [{c[0]:.4f}, {c[1]:.4f}]"
        return s

    L = [f"JUNCTION-DISTANCE CONTROL  pool={a.pool}  n={n:,} (P {int(y.sum()):,} / B {int(n - y.sum()):,})",
         f"gene-clustered bootstrap: {a.n_boot} draws over {nc:,} clusters; cross-fitting: {a.folds} gene-grouped folds",
         f"distance convention: intronic first base = 1; exonic last/first base = 1; GTF tag filter: '{a.tag or 'none'}'", ""]
    for setting, title in (("unmatched", "A. ALL BENIGN (raw ClinVar-pool distances)"),
                           ("matched", "B. DISTANCE-MATCHED BENIGN CONTROLS (benign re-weighted to the pathogenic distance profile)")):
        L += [title,
              f"  engine frozen call: sensitivity {f(setting + '|engine_call_sens')} | FPR {f(setting + '|engine_call_fpr')} | PPV {f(setting + '|engine_call_ppv')}",
              f"  proximity-only sensitivity at the engine's FPR: {f(setting + '|proximity_sens_at_engine_fpr')}",
              f"  engine minus proximity (sensitivity at matched FPR): {f(setting + '|sens_engine_minus_proximity')}",
              f"  AUROC  engine S_raw {f(setting + '|engine_auroc', False)} | raw proximity {f(setting + '|proximity_auroc', False)}",
              f"  AUPRC  engine S_raw {f(setting + '|engine_auprc', False)} | raw proximity {f(setting + '|proximity_auprc', False)}",
              f"  cross-fitted AUROC  proximity model {f(setting + '|model_prox_auroc_oof', False)} | + engine {f(setting + '|model_prox_plus_engine_auroc_oof', False)}",
              f"  delta AUROC (+engine - proximity) {f(setting + '|delta_auroc', False)}",
              f"  cross-fitted AUPRC  proximity model {f(setting + '|model_prox_auprc_oof', False)} | + engine {f(setting + '|model_prox_plus_engine_auprc_oof', False)}",
              f"  delta AUPRC {f(setting + '|delta_auprc', False)}", ""]
    L += [f"C. LIKELIHOOD-RATIO TEST (in-sample, nested logistic models; engine adds {df} terms)",
          f"  overall: chi2 = {stat:.2f}, df = {df}, p = {pval:.3g}   (variants within a gene are not independent; rely on the bootstrap deltas)"]
    for nm, v in strat_lrt.items():
        L.append(f"  {nm:8s}: chi2 = {v['chi2']:.2f}, p = {v['p']:.3g}  (P {v['n_pathogenic']:,} / B {v['n_benign']:,})")
    L += ["", "D. STRATIFIED BY DISTANCE (engine frozen call; AUROC of S_raw within the stratum)"]
    for setting in ("unmatched", "matched"):
        L.append(f"  [{setting}]  group | P / B (raw n) | sens | FPR | PPV | AUROC")
        for gname, (idx, _) in groups.items():
            k = f"{setting}|G|{gname}|"
            L.append(f"    {gname:16s} {int(y[idx].sum()):>6,} / {int(len(idx) - y[idx].sum()):>8,} | {f(k + 'sens')} | {f(k + 'fpr')} | {f(k + 'ppv')} | {f(k + 'auroc', False)}")
    L += ["", "Fraction of bootstrap draws with positive delta: " + json.dumps(frac_pos),
          "",
          "NOTES",
          "  - In section B the benign set is re-weighted so pathogenic prevalence is 50%; PPV there is PPV at 50% prevalence,",
          "    not at the pool prevalence. Compare sensitivity, FPR and AUROC across A and B; do not quote B's PPV as a real-world PPV.",
          "  - The delta-AUROC/AUPRC intervals resample genes but keep the cross-fitted predictions fixed, so they do not include",
          "    model-fitting uncertainty. In a null test (engine carrying only distance information) the delta was +0.002 to +0.003",
          "    AUROC with an interval that barely excluded 0. Treat deltas below about 0.01 as null-level and read them with the LRT.",
          "  - The LRT treats variants as independent; variants within a gene are correlated, so its p-values are optimistic.",
          "  - Engine features in the models: S_raw>0, log1p(S_raw) and the frozen call. Other engine outputs (delta_pwm, direction, ...)",
          "    are not in scores.tsv, so they are not tested here."]
    txt = "\n".join(L)
    open(os.path.join(a.out, "report.txt"), "w").write(txt + "\n")
    print("\n" + txt)
    print(f"\nOutputs in {a.out}\\  (report.txt, report.json, variant_distance.tsv)  [{time.time() - t0:,.0f}s]")


if __name__ == "__main__":
    main()