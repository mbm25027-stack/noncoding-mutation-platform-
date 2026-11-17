#!/usr/bin/env python3
"""
nc_spliceai_complementarity.py - does SpliceAI add to / overlap with the FROZEN
NC-MutScan v3.3 call?  Reads cached engine scores (scores.tsv) and a per-variant
SpliceAI score file; no rescoring, no engine change. numpy + standard library.

  A. Coverage: which variants SpliceAI has no score for (SNV/indel, label), and
     how NC-MutScan does on exactly those. Unscorable variants stay in every
     denominator and count as SpliceAI-negative.
  B. 2x2 discordance tables (NC call x SpliceAI call at 0.2 / 0.5 / 0.8), shown
     separately for pathogenic and benign variants, with exact McNemar tests on
     the discordant pairs (NC-only vs SpliceAI-only), and gene-clustered
     bootstrap CIs for the paired sensitivity difference.
  C. Matched specificity: the rule "NC-MutScan frozen call OR SpliceAI >= t".
     For each target FPR (a multiple of the NC call's own FPR) t is chosen as the
     lowest threshold keeping the union's FPR within the target; SpliceAI alone
     is given the SAME FPR by its own threshold. Union gain = union sensitivity
     minus SpliceAI-alone sensitivity at equal FPR (bootstrap re-selects t in
     every draw, so selection uncertainty is included).

SpliceAI file: tab-separated, columns  key<TAB>ds_max  (key = chrom:pos:ref:alt,
chrom without "chr"; ds_max = max of the four delta scores). Make it with
extract_spliceai_scores.py. Variants absent from the file are "unscorable".

Example:
  python nc_spliceai_complementarity.py --pool check_blind_v2\\pool_all.vcf ^
     --scores eval_check_v2\\scores.tsv --spliceai spliceai_clean.tsv ^
     --exclude-mc none --out eval_spliceai_clean
"""
import argparse, json, math, os, time
from collections import Counter, defaultdict
import numpy as np

PATH_SIG = {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"}
BEN_SIG = {"benign", "likely_benign", "benign/likely_benign"}
THRS = [0.2, 0.5, 0.8]
FACTORS = [1.0, 1.25, 1.5, 2.0, 3.0]
NAN = float("nan")


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
            rows.append((t[0].replace("chr", ""), int(t[1]), t[3].upper(), t[4].upper(), y, gene, "|".join(sorted(terms))))
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


def load_spliceai(path):
    d = {}
    with open(path) as f:
        for line in f:
            t = line.rstrip("\r\n").split("\t")
            if len(t) < 2 or t[0] == "key":
                continue
            try:
                v = float(t[1])
            except ValueError:
                continue
            if t[0] not in d or v > d[t[0]]:
                d[t[0]] = v
    return d


def mcnemar_exact(b, c):
    n = int(b + c)
    if n == 0:
        return 1.0
    k = int(min(b, c))
    lp = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * math.log(2) for i in range(k + 1)]
    m = max(lp)
    return min(1.0, 2 * math.exp(m) * sum(math.exp(x - m) for x in lp))


# ----------------------------------------------------- matched-FPR solver
class Solver:
    """Union rule  call OR ds>=t  (or SpliceAI alone), t chosen under an FPR cap."""

    def __init__(self, y, c, ds):
        self.y1, self.y0 = y == 1, y == 0
        sc = ds >= 0
        self.P_c1 = np.nonzero(self.y1 & (c == 1))[0]
        self.B_c1 = np.nonzero(self.y0 & (c == 1))[0]

        def srt(mask):
            idx = np.nonzero(mask)[0]
            idx = idx[np.argsort(-ds[idx], kind="stable")]
            v = ds[idx]
            ends = np.r_[np.nonzero(np.diff(v))[0], len(v) - 1] if len(v) else np.array([], dtype=int)
            return idx, v, ends

        self.U = (srt(self.y1 & (c == 0) & sc), srt(self.y0 & (c == 0) & sc))
        self.S = (srt(self.y1 & sc), srt(self.y0 & sc))

    def solve(self, w, target, union):
        wy1, wy0 = w[self.y1].sum(), w[self.y0].sum()
        if wy1 <= 0 or wy0 <= 0:
            return dict(t=NAN, fpr=NAN, sens=NAN)
        (Pi, Pv, _), (Bi, Bv, Be) = self.U if union else self.S
        baseP = w[self.P_c1].sum() if union else 0.0
        baseB = w[self.B_c1].sum() if union else 0.0
        cap = target * wy0
        if baseB > cap * (1 + 1e-9) + 1e-12:
            return dict(t=NAN, fpr=NAN, sens=NAN)
        cwp = np.cumsum(w[Pi]) if len(Pi) else np.array([])

        def p_gt(v):   # pathogenic weight with ds >  v
            m = int(np.searchsorted(-Pv, -v, side="left"))
            return float(cwp[m - 1]) if m > 0 else 0.0

        def p_ge(v):   # pathogenic weight with ds >= v
            m = int(np.searchsorted(-Pv, -v, side="right"))
            return float(cwp[m - 1]) if m > 0 else 0.0

        if len(Bi) == 0:    # no scorable benign competitors: every scorable pathogenic is free
            return dict(t=0.0, fpr=baseB / wy0, sens=(baseP + (cwp[-1] if len(cwp) else 0.0)) / wy1, phi=1.0)
        cw = np.cumsum(w[Bi])
        tot = baseB + cw[Be]                       # non-decreasing over tie-group ends
        nfeas = int(np.searchsorted(tot, cap * (1 + 1e-9) + 1e-12, side="right"))
        if nfeas == len(Be):                       # every benign group fits: include all scorable
            return dict(t=0.0, fpr=(baseB + cw[-1]) / wy0, sens=(baseP + (cwp[-1] if len(cwp) else 0.0)) / wy1, phi=1.0)
        # groups before nfeas are fully included; the next tie group (value vinf) is included
        # fractionally (standard ROC interpolation across ties) so the FPR matches the target
        vinf = float(Bv[Be[nfeas]])
        bin_ = float(tot[nfeas - 1]) if nfeas > 0 else baseB
        gB = float(cw[Be[nfeas]] - (cw[Be[nfeas - 1]] if nfeas > 0 else 0.0))
        phi = min(1.0, max(0.0, (cap - bin_) / gB)) if gB > 0 else 0.0
        sensP = baseP + p_gt(vinf) + phi * (p_ge(vinf) - p_gt(vinf))
        return dict(t=vinf, fpr=(bin_ + phi * gB) / wy0, sens=sensP / wy1, phi=phi)


# ------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--scores", required=True, help="scores.tsv of the same pool (heldout_full_pool_eval.py)")
    ap.add_argument("--spliceai", required=True, help="TSV key<TAB>ds_max (see extract_spliceai_scores.py)")
    ap.add_argument("--exclude-mc", default="nonsense,inframe_indel,no_sequence_alteration")
    ap.add_argument("--out", default="eval_spliceai")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--snv-only", action="store_true",
                    help="drop all indels from the pool (the AbSplice/Zenodo file has no indel scores, so indels are unscorable)")
    a = ap.parse_args()
    t0 = time.time()

    rows, skipped = load_pool(a.pool, [x for x in a.exclude_mc.split(",") if x])
    if a.snv_only:
        n0 = len(rows)
        rows = [r for r in rows if len(r[2]) == 1 and len(r[3]) == 1]
        skipped["indel_dropped_snv_only"] = n0 - len(rows)
    sc = load_scores(a.scores)
    sp = load_spliceai(a.spliceai)
    keep, call, ds, srw, miss = [], [], [], [], 0
    for r in rows:
        k = f"{r[0]}:{r[1]}:{r[2]}:{r[3]}"
        s = sc.get(k)
        if s is None:
            miss += 1
            continue
        keep.append(r)
        call.append(s[1])
        srw.append(s[0])
        ds.append(sp.get(k, -1.0))
    y = np.array([r[4] for r in keep])
    c = np.array(call, dtype=int)
    ds = np.array(ds, dtype=float)
    has = ds >= 0
    snv = np.array([len(r[2]) == 1 and len(r[3]) == 1 for r in keep])
    n = len(y)
    P, B = int((y == 1).sum()), int((y == 0).sum())
    print(f"Pool: {n:,} | P {P:,} | B {B:,} | skipped {dict(skipped)} | no engine score {miss:,}")
    print(f"SpliceAI file: {len(sp):,} keys; pool variants with a score: {int(has.sum()):,} ({has.mean():.1%})")
    if has.sum() == 0:
        raise SystemExit("No pool variant found in the SpliceAI file - check the key format (chrom:pos:ref:alt, no 'chr').")

    # ---- pathogenic variants with their discordance category (for a supplementary table)
    os.makedirs(a.out, exist_ok=True)
    cats = Counter()
    with open(os.path.join(a.out, "pathogenic_discordance.tsv"), "w") as o:
        o.write("key\tgene\tmc\tindel\tnc_call\tS_raw\tspliceai_ds\tscorable\tcat_0.2\tcat_0.5\tcat_0.8\n")
        for i, r in enumerate(keep):
            if r[4] != 1:
                continue
            cs = []
            for th in (0.2, 0.5, 0.8):
                sa = ds[i] >= th
                cs.append("both" if call[i] and sa else "NC_only" if call[i] else "SpliceAI_only" if sa else "neither")
            if cs[0] == "NC_only":
                cats[(r[6] or "NA")] += 1
            o.write(f"{r[0]}:{r[1]}:{r[2]}:{r[3]}\t{r[5]}\t{r[6]}\t{int(not (len(r[2]) == 1 and len(r[3]) == 1))}\t{call[i]}\t{srw[i]}\t"
                    f"{ds[i] if ds[i] >= 0 else 'NA'}\t{int(ds[i] >= 0)}\t" + "\t".join(cs) + "\n")
    print("NC-only pathogenic variants at SpliceAI>=0.2 by MC term (top 8): " + "; ".join(f"{k}: {v}" for k, v in cats.most_common(8)))

    names = [r[5] if r[5] else f"nogene_{r[0]}_{r[1] // 1_000_000}" for r in keep]
    uniq = {}
    cidx = np.array([uniq.setdefault(nm, len(uniq)) for nm in names])
    nc = len(uniq)
    solver = Solver(y, c, ds)
    P1, B0 = y == 1, y == 0

    def wsum(w, m):
        return float(w[m].sum())

    def evaluate(w):
        o = {}
        wP, wB = wsum(w, P1), wsum(w, B0)
        for thr in THRS:
            s = ds >= thr
            cell = {}
            for lab, m, den in (("P", P1, wP), ("B", B0, wB)):
                cell[lab] = (wsum(w, m & (c == 1) & s), wsum(w, m & (c == 1) & ~s),
                             wsum(w, m & (c == 0) & s), wsum(w, m & (c == 0) & ~s))
            aP, bP, cP, dP = cell["P"]
            aB, bB, cB, dB = cell["B"]
            tg = f"thr{thr}"
            o[f"{tg}|sens_NC"] = (aP + bP) / wP
            o[f"{tg}|sens_SpliceAI"] = (aP + cP) / wP
            o[f"{tg}|sens_union"] = (aP + bP + cP) / wP
            o[f"{tg}|sens_diff_NC_minus_SpliceAI"] = (bP - cP) / wP
            o[f"{tg}|union_gain_over_SpliceAI"] = bP / wP
            o[f"{tg}|union_gain_over_NC"] = cP / wP
            o[f"{tg}|fpr_NC"] = (aB + bB) / wB
            o[f"{tg}|fpr_SpliceAI"] = (aB + cB) / wB
            o[f"{tg}|fpr_union"] = (aB + bB + cB) / wB
            o[f"{tg}|ppv_union"] = (aP + bP + cP) / (aP + bP + cP + aB + bB + cB) if (aP + bP + cP + aB + bB + cB) > 0 else NAN
        fN = (wsum(w, B0 & (c == 1))) / wB
        o["fpr_NC_call"] = fN
        for fac in FACTORS:
            tgt = fac * fN
            su = solver.solve(w, tgt, True)
            ss = solver.solve(w, tgt, False)
            g = f"F{fac}"
            o[f"{g}|target_fpr"] = tgt
            o[f"{g}|sens_SpliceAI_alone"] = ss["sens"]
            o[f"{g}|fpr_SpliceAI_alone"] = ss["fpr"]
            o[f"{g}|thr_SpliceAI_alone"] = ss["t"]
            o[f"{g}|sens_union"] = su["sens"]
            o[f"{g}|fpr_union"] = su["fpr"]
            o[f"{g}|thr_union"] = su["t"]
            o[f"{g}|union_gain_vs_SpliceAI"] = su["sens"] - ss["sens"]
        # NC alone on the unscorable / scorable pathogenic variants
        for nm, m in (("unscorable", ~has), ("scorable", has)):
            d_ = wsum(w, P1 & m)
            o[f"NC_sens_among_{nm}_pathogenic"] = wsum(w, P1 & m & (c == 1)) / d_ if d_ > 0 else NAN
        return o

    one = np.ones(n)
    point = evaluate(one)

    rng = np.random.default_rng(a.seed)
    draws = defaultdict(list)
    for b in range(a.n_boot):
        cnt = np.bincount(rng.integers(0, nc, nc), minlength=nc)
        for k, v in evaluate(cnt[cidx].astype(float)).items():
            draws[k].append(v)
        if (b + 1) % 250 == 0:
            print(f"  bootstrap {b + 1}/{a.n_boot}  ({time.time() - t0:,.0f}s)")
    ci = {}
    for k, v in draws.items():
        if "|thr_" in k:
            continue
        arr = np.array(v, dtype=float)
        ci[k] = [float(x) for x in np.nanpercentile(arr, [2.5, 97.5])] if np.isfinite(arr).any() else [NAN, NAN]

    # counts-based tables and exact McNemar (point data)
    tables = {}
    for thr in THRS:
        s = ds >= thr
        t_ = {}
        for lab, m in (("pathogenic", P1), ("benign", B0)):
            ab = int((m & (c == 1) & s).sum()); bb = int((m & (c == 1) & ~s).sum())
            cc = int((m & (c == 0) & s).sum()); dd = int((m & (c == 0) & ~s).sum())
            t_[lab] = dict(both=ab, NC_only=bb, SpliceAI_only=cc, neither=dd, mcnemar_exact_p=mcnemar_exact(bb, cc))
        tables[str(thr)] = t_

    cov = {}
    for lab, m in (("pathogenic", P1), ("benign", B0)):
        for vt, mv in (("SNV", snv), ("indel", ~snv)):
            q = m & mv
            cov[f"{lab}_{vt}"] = dict(n=int(q.sum()), scorable=int((q & has).sum()))

    os.makedirs(a.out, exist_ok=True)
    rep = dict(pool=a.pool, scores=a.scores, spliceai=a.spliceai, exclude_mc=a.exclude_mc, n=n, n_pathogenic=P,
               n_benign=B, n_clusters=nc, n_boot=a.n_boot, seed=a.seed, coverage=cov, tables=tables,
               metrics={k: dict(estimate=v, ci95=ci.get(k)) for k, v in point.items()})
    json.dump(rep, open(os.path.join(a.out, "report.json"), "w"), indent=2)

    def f(k, pct=True):
        v, q = point[k], ci.get(k)
        if v is None or not np.isfinite(v):
            return "n/a"
        s = f"{100 * v:.2f}%" if pct else f"{v:.4g}"
        if q and np.isfinite(q[0]):
            s += f" [{100 * q[0]:.2f}%, {100 * q[1]:.2f}%]" if pct else f" [{q[0]:.4g}, {q[1]:.4g}]"
        return s

    L = [f"SPLICEAI COMPLEMENTARITY  pool={a.pool}  n={n:,} (P {P:,} / B {B:,})  clusters {nc:,}  boot {a.n_boot}",
         f"SpliceAI file: {a.spliceai}; unscorable = absent from that file (kept in all denominators, counted SpliceAI-negative)", "",
         "A. COVERAGE (variants with a SpliceAI score)"]
    for k, v in cov.items():
        L.append(f"  {k:20s} {v['scorable']:>8,} / {v['n']:>8,}  ({(v['scorable'] / v['n']) if v['n'] else float('nan'):.1%})")
    L.append(f"  NC-MutScan frozen call on UNSCORABLE pathogenic variants: {f('NC_sens_among_unscorable_pathogenic')}")
    L.append(f"  NC-MutScan frozen call on SCORABLE pathogenic variants:   {f('NC_sens_among_scorable_pathogenic')}")
    L.append(f"  NC-MutScan frozen call FPR (all benign): {f('fpr_NC_call')}")
    L += ["", "B. 2x2 DISCORDANCE (NC frozen call x SpliceAI >= threshold); exact McNemar on NC-only vs SpliceAI-only"]
    for thr in THRS:
        tg = f"thr{thr}"
        L.append(f"  SpliceAI >= {thr}")
        for lab in ("pathogenic", "benign"):
            t_ = tables[str(thr)][lab]
            L.append(f"    {lab:10s} both {t_['both']:>7,} | NC only {t_['NC_only']:>7,} | SpliceAI only {t_['SpliceAI_only']:>7,} | "
                     f"neither {t_['neither']:>9,} | McNemar p = " + (f"{t_['mcnemar_exact_p']:.3g}" if t_['mcnemar_exact_p'] > 1e-300 else "< 1e-300"))
        L += [f"    sensitivity: NC {f(tg + '|sens_NC')} | SpliceAI {f(tg + '|sens_SpliceAI')} | union {f(tg + '|sens_union')}",
              f"    FPR:         NC {f(tg + '|fpr_NC')} | SpliceAI {f(tg + '|fpr_SpliceAI')} | union {f(tg + '|fpr_union')}",
              f"    paired sensitivity difference NC - SpliceAI: {f(tg + '|sens_diff_NC_minus_SpliceAI')}",
              f"    union gain over SpliceAI (NC-only pathogenic / P): {f(tg + '|union_gain_over_SpliceAI')}",
              f"    union gain over NC (SpliceAI-only pathogenic / P): {f(tg + '|union_gain_over_NC')}"]
    L += ["", "C. MATCHED SPECIFICITY: union = NC frozen call OR SpliceAI >= t, vs SpliceAI alone at the SAME FPR",
          "   (target FPR = factor x the NC call's own FPR; bootstrap re-selects the thresholds each draw)"]
    for fac in FACTORS:
        g = f"F{fac}"
        L += [f"  factor {fac}: target FPR {f(g + '|target_fpr')}",
              f"    SpliceAI alone: sens {f(g + '|sens_SpliceAI_alone')} at FPR {f(g + '|fpr_SpliceAI_alone')} (t = {point[g + '|thr_SpliceAI_alone']:.3g})",
              f"    union:          sens {f(g + '|sens_union')} at FPR {f(g + '|fpr_union')} (t = {point[g + '|thr_union']:.3g})",
              f"    union gain over SpliceAI alone: {f(g + '|union_gain_vs_SpliceAI')}"]
    L += ["", "NOTES",
          "  - 'sens ... (t = x)': variants scoring above x are included; the tie group at exactly x is included fractionally so the FPR",
          "    equals the target (standard ROC interpolation; SpliceAI scores have many ties at 0). t = 0 means every scorable variant fits.",
          "  - factor 1.0: the union may not exceed the NC call's own FPR, so SpliceAI can only add pathogenic variants scoring above every",
          "    benign non-NC variant; read it also as NC call vs SpliceAI alone at the NC call's FPR (SpliceAI-alone line).",
          "    The union threshold there hangs on the single highest-scoring benign variant, so its interval is wide; lean on factors >= 1.25.",
          "  - If 'union gain' intervals include 0 and the SpliceAI-only column in B is small, the union adds nothing.",
          "  - McNemar p-values treat variants as independent; variants in one gene are correlated, so lean on the clustered intervals.",
          "  - The NC call is the FROZEN call; it is not re-tuned. Only the SpliceAI threshold is selected."]
    txt = "\n".join(L)
    open(os.path.join(a.out, "report.txt"), "w").write(txt + "\n")
    print("\n" + txt)
    print(f"\nOutputs in {a.out}\\  (report.txt, report.json, pathogenic_discordance.tsv)  [{time.time() - t0:,.0f}s]")


if __name__ == "__main__":
    main()