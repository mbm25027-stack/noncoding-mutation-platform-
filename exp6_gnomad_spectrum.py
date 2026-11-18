#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
exp6_gnomad_spectrum.py -- Experiment 6: allele-frequency spectrum of NC-MutScan v3.3 calls (gnomAD)
=====================================================================================================
Label-free, population-level test.  If NC-MutScan calls flag variants that disrupt functional
signals, purifying selection should make calls

    * DEPLETED among common variants (AF >= 1%), and
    * ENRICHED among singletons (AC = 1),

after holding constant the factors that change both mutability and allele frequency:
distance to the exon edge (or TSS), CpG status / transition-vs-transversion, and local GC.

Design
------
1. SAMPLE  stream gnomAD genomes sites VCF(s); keep PASS, biallelic SNVs, adequate coverage (AN),
           outside every CDS, inside one of three noncoding domains (utr / splice +-100 bp / promoter);
           reservoir-sample up to --per-class variants from each frequency class.
2. SCORE   build +-25 bp WT/MUT windows from the GRCh38 FASTA, call your frozen engine
           (my_engine.py, the same adapter as experiment 1), derive CpG / GC.
3. ANALYSE exact stratification (= exact matching) on domain x distance-bin x CpG/ts-tv x GC-bin.
           Per frequency class: observed vs expected calls (O/E; expected from stratum-specific
           call rates), Mantel-Haenszel odds ratio vs the common class, and a stratified
           log-odds trend per frequency-class step.  All CIs are block-bootstrap (100 kb blocks).
           A dose-response over score thresholds and per-module / per-domain / per-distance
           breakdowns are included.  Sampling per class is random within class, so
           class comparisons inside strata are unbiased.

Run (PowerShell, project folder, venv active; start with the two smallest chromosomes):

  python exp6_gnomad_spectrum.py --gnomad gnomad.genomes.v4.1.sites.chr21.vcf.bgz gnomad.genomes.v4.1.sites.chr22.vcf.bgz ^
         --gtf gencode.v36.annotation.gtf --fasta GRCh38.d1.vd1.fa --build-hash b276d0d44045

Stages can be run separately (--stage sample | score | analyse) and are restartable from the files in --out.
The gnomAD files are big (10-100 GB per chromosome); streaming is decompression-bound, so use the
smallest chromosomes first.  chr12-22 are held-out chromosomes for the frozen build, a natural choice.
"""
from __future__ import annotations

import argparse
import csv
import json
import random
import re
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

from nc_exp_common import (Fasta, chrom_order, chrom_set, load_engine, log, make_window, normalise_result,
                           norm_chrom, open_text, write_csv)
from nc_genome_features import GenomeFeatures

CLASSES = ["singleton", "ultra_rare", "rare", "low_freq", "common", "very_common"]
CLASS_DEF = ["AC = 1", "AC >= 2, AF < 1e-4", "1e-4 <= AF < 1e-3", "1e-3 <= AF < common cut",
             "common cut <= AF < 5%", "AF >= 5%"]
REF_CLASS = CLASSES.index("common")
_AC = re.compile(r"(?:^|;)AC=(\d+)")
_AN = re.compile(r"(?:^|;)AN=(\d+)")
_AF = re.compile(r"(?:^|;)AF=([0-9.eE+-]+)")


def af_class(ac, af, common):
    if ac == 1:
        return 0
    if af < 1e-4:
        return 1
    if af < 1e-3:
        return 2
    if af < common:
        return 3
    if af < 0.05:
        return 4
    return 5


# ===========================================================================
# stage 1: stream + reservoir sample
# ===========================================================================
def stage_sample(args, out):
    chroms = chrom_set(args.chroms)
    gf = GenomeFeatures(args.gtf, chroms, out / "gtf_features.npz" if not args.no_cache else None,
                        args.splice_window, args.prom_up, args.prom_down)
    rng = random.Random(args.seed)
    res = {k: [] for k in range(len(CLASSES))}
    seen = Counter()
    stats = Counter()

    def consider(rec, thr):
        chrom, pos, ref, alt, ac, an, af, dom, dist = rec
        if an < thr:
            stats["low_an"] += 1
            return
        k = af_class(ac, af, args.common_af)
        seen[k] += 1
        r = res[k]
        item = (chrom, pos, ref, alt, ac, an, af, CLASSES[k], dom, dist)
        if len(r) < args.per_class:
            r.append(item)
        else:
            j = rng.randrange(seen[k])
            if j < args.per_class:
                r[j] = item

    for path in args.gnomad:
        with open_text(path) as _fh:
            first = _fh.readline()
        if not first.startswith("##fileformat"):
            raise SystemExit(f"{path} does not look like a VCF (first line: {first[:80]!r}). If it is a few KB "
                             "it may be an error page from a failed download: delete it and download again.")
        log(f"streaming {path}")
        buf, thr, t0 = [], None, time.time()
        with open_text(path) as fh:
            for ln, line in enumerate(fh, 1):
                if args.max_lines and ln > args.max_lines:
                    break
                if line[0] == "#":
                    continue
                if ln % 5_000_000 == 0:
                    log(f"  {ln:,} lines  ({ln / (time.time() - t0):,.0f}/s)  candidates {stats['cand']:,}")
                f = line.split("\t", 8)
                if len(f) < 8 or len(f[3]) != 1 or len(f[4]) != 1 or f[6] != "PASS":
                    continue
                chrom = norm_chrom(f[0])
                pos = int(f[1])
                if chroms and chrom not in chroms or not gf.prefilter(chrom, pos):
                    continue
                cl = gf.classify(chrom, pos)
                if cl is None:
                    stats["coding_or_outside"] += 1
                    continue
                info = f[7]
                m_ac, m_an = _AC.search(info), _AN.search(info)
                if not m_ac or not m_an:
                    continue
                ac, an = int(m_ac.group(1)), int(m_an.group(1))
                if ac == 0 or an == 0:
                    continue
                m_af = _AF.search(info)
                af = float(m_af.group(1)) if m_af else ac / an
                stats["cand"] += 1
                rec = (chrom, pos, f[3].upper(), f[4].upper(), ac, an, af, cl[0], cl[1])
                if thr is None:
                    buf.append(rec)
                    if len(buf) >= 5000:
                        thr = args.min_an_frac * float(np.percentile([b[5] for b in buf], 90))
                        log(f"  AN threshold for this file: {thr:,.0f}")
                        for b in buf:
                            consider(b, thr)
                        buf = []
                else:
                    consider(rec, thr)
        if buf:
            thr = args.min_an_frac * float(np.percentile([b[5] for b in buf], 90))
            for b in buf:
                consider(b, thr)
    rows = [dict(chrom=i[0], pos=i[1], ref=i[2], alt=i[3], ac=i[4], an=i[5], af=f"{i[6]:.6g}",
                 afclass=i[7], domain=i[8], dist=i[9]) for k in res for i in res[k]]
    rows.sort(key=lambda r: (chrom_order(r["chrom"]), r["pos"]))
    write_csv(out / "exp6_sampled.tsv", rows)
    meta = dict(seen={CLASSES[k]: seen[k] for k in range(len(CLASSES))},
                kept={CLASSES[k]: len(res[k]) for k in range(len(CLASSES))}, stats=dict(stats),
                files=[str(p) for p in args.gnomad], per_class=args.per_class, common_af=args.common_af)
    (out / "exp6_sampling.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    log(f"sampled {len(rows):,} variants; seen per class {meta['seen']}")
    if not rows:
        raise SystemExit("No variants sampled - check chromosome naming, --chroms and the gnomAD file.")


# ===========================================================================
# stage 2: score
# ===========================================================================
def stage_score(args, out):
    fa = Fasta(args.fasta)
    engine = load_engine(args.engine)
    with open(out / "exp6_sampled.tsv", newline="", encoding="utf-8") as fh:
        sampled = list(csv.DictReader(fh))
    rows, bad, t0 = [], 0, time.time()
    for i, r in enumerate(sampled):
        pos = int(r["pos"])
        w = make_window(fa, r["chrom"], pos, r["ref"], r["alt"], args.flank)
        if w is None:
            bad += 1
            continue
        wt, mut, c = w
        meta = dict(chrom=r["chrom"], pos=pos, ref=r["ref"], alt=r["alt"], gene="", center=c)
        try:
            res = normalise_result(engine(wt, mut, meta))
        except NotImplementedError:
            raise SystemExit("my_engine.py is not wired yet (see my_engine_template.py).")
        prev, nxt = (wt[c - 1] if c > 0 else ""), (wt[c + 1] if c + 1 < len(wt) else "")
        cpg = (r["ref"] == "C" and nxt == "G") or (r["ref"] == "G" and prev == "C")
        ts = {r["ref"], r["alt"]} in ({"A", "G"}, {"C", "T"})
        rows.append({**r, **res, "cpg": ("cpg_" if cpg else "non_") + ("ts" if ts else "tv"),
                     "gc": round((wt.count("G") + wt.count("C")) / len(wt), 4)})
        if (i + 1) % 50000 == 0:
            log(f"  scored {i + 1:,}/{len(sampled):,} ({(i + 1) / (time.time() - t0):,.0f}/s)")
    feat = sorted({k for r in rows for k in r if k.startswith("feat_")})
    cols = ["chrom", "pos", "ref", "alt", "ac", "an", "af", "afclass", "domain", "dist", "cpg", "gc",
            "s_raw", "s_unified", "call", "module"] + feat
    write_csv(out / "exp6_per_variant.csv", rows, cols)
    log(f"scored {len(rows):,} variants ({bad:,} reference mismatches skipped)")


# ===========================================================================
# stage 3: stratified analysis
# ===========================================================================
def dist_bin(domain, d):
    if domain == "utr":
        return "utr"
    if domain == "splice":
        return "s1-2" if d <= 2 else "s3-5" if d <= 5 else "s6-10" if d <= 10 else "s11-20" if d <= 20 \
            else "s21-50" if d <= 50 else "s51+"
    return "p0-50" if d <= 50 else "p51-150" if d <= 150 else "p151-400" if d <= 400 else "p401+"


def gc_bin(g):
    return "gc1" if g < 0.35 else "gc2" if g < 0.50 else "gc3" if g < 0.65 else "gc4"


class Cells:
    """Per-block counts of (stratum, class) cells; bootstrap = multiply by block multiplicities."""

    def __init__(self, strat, cls, blk, K):
        self.sidx, self.S = np.unique(strat, return_inverse=True)[1], len(np.unique(strat))
        self.bidx, self.nb = np.unique(blk, return_inverse=True)[1], len(np.unique(blk))
        self.K = K
        self.cell = self.sidx * K + np.asarray(cls)
        flat = self.bidx * (self.S * K) + self.cell
        self.tot = np.bincount(flat, minlength=self.nb * self.S * K).reshape(self.nb, self.S * K).astype(float)

    def cal(self, y):
        flat = self.bidx * (self.S * self.K) + self.cell
        return np.bincount(flat, weights=np.asarray(y, float), minlength=self.nb * self.S * self.K) \
            .reshape(self.nb, self.S * self.K)


def point_stats(tot, cal, ref):
    S, K = tot.shape
    n_s, c_s = tot.sum(1), cal.sum(1)
    p_s = np.divide(c_s, n_s, out=np.zeros(S), where=n_s > 0)
    exp, obs = (tot * p_s[:, None]).sum(0), cal.sum(0)
    oe = np.divide(obs, exp, out=np.full(K, np.nan), where=exp > 0)
    mh = np.full(K, np.nan)
    c, d = cal[:, ref], tot[:, ref] - cal[:, ref]
    for k in range(K):
        if k == ref:
            continue
        a, b = cal[:, k], tot[:, k] - cal[:, k]
        n = tot[:, k] + tot[:, ref]
        ok = n > 0
        num, den = (a * d / np.where(ok, n, 1))[ok].sum(), (b * c / np.where(ok, n, 1))[ok].sum()
        if num > 0 and den > 0:
            mh[k] = num / den
    return oe, mh, obs, exp


def trend_slope(tot, cal, xs):
    """Stratified binomial log-odds slope per class step (stratum fixed effects, IRLS)."""
    S, K = tot.shape
    ns, cs = tot.sum(1), cal.sum(1)
    keep = (ns > 0) & (cs > 0) & (cs < ns) & ((tot > 0).sum(1) >= 2)
    if keep.sum() < 2:
        return np.nan
    T, C = tot[keep], cal[keep]
    si, ki = np.nonzero(T > 0)
    n, y = T[si, ki], C[si, ki]
    X = np.zeros((len(si), keep.sum() + 1))
    X[np.arange(len(si)), si] = 1.0
    X[:, -1] = np.asarray(xs, float)[ki] - np.mean(xs)
    beta = np.zeros(X.shape[1])
    pr = (cs[keep] / ns[keep])
    beta[:-1] = np.log(pr / (1 - pr))
    for _ in range(30):
        eta = np.clip(X @ beta, -30, 30)
        p = 1 / (1 + np.exp(-eta))
        g = X.T @ (y - n * p)
        H = X.T @ (X * (n * p * (1 - p))[:, None]) + 1e-8 * np.eye(len(beta))
        step = np.linalg.solve(H, g)
        beta += step
        if np.abs(step).max() < 1e-9:
            break
    return float(beta[-1])


def analyse_one(cells, y, B, Bt, rng, ref, do_trend=True):
    cal = cells.cal(y)
    S, K = cells.S, cells.K
    tot0, cal0 = cells.tot.sum(0).reshape(S, K), cal.sum(0).reshape(S, K)
    oe, mh, obs, exp = point_stats(tot0, cal0, ref)
    xs = np.arange(K)
    slope = trend_slope(tot0, cal0, xs) if do_trend else np.nan
    bo, bm, bs = [], [], []
    for b in range(B):
        w = np.bincount(rng.integers(0, cells.nb, cells.nb), minlength=cells.nb).astype(float)
        t, c = (w @ cells.tot).reshape(S, K), (w @ cal).reshape(S, K)
        o, m, _, _ = point_stats(t, c, ref)
        bo.append(o)
        bm.append(m)
        if do_trend and b < Bt:
            bs.append(trend_slope(t, c, xs))
    q = lambda a: np.nanpercentile(np.array(a), [2.5, 97.5], axis=0)
    res = dict(n=tot0.sum(0).tolist(), obs=obs.tolist(), exp=exp.tolist(), oe=oe.tolist(),
               oe_ci=q(bo).tolist(), mh_or=mh.tolist(), mh_ci=q(bm).tolist())
    if do_trend and bs:
        bs = np.array(bs, float)
        bs = bs[~np.isnan(bs)]
        res["slope"] = dict(est=slope, lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)),
                            p=float(min(1.0, 2 * min((bs >= 0).mean(), (bs <= 0).mean()))), n_boot=len(bs))
    return res


def stage_analyse(args, out):
    import warnings
    warnings.simplefilter("ignore", RuntimeWarning)      # all-NaN slices for empty cells are expected
    rng = np.random.default_rng(args.seed)
    with open(out / "exp6_per_variant.csv", newline="", encoding="utf-8") as fh:
        R = list(csv.DictReader(fh))
    n = len(R)
    log(f"analysing {n:,} scored variants")
    cls = np.array([CLASSES.index(r["afclass"]) for r in R])
    dom = np.array([r["domain"] for r in R])
    dist = np.array([int(r["dist"]) for r in R])
    s_raw = np.array([float(r["s_raw"]) for r in R])
    call = np.array([int(float(r["call"])) for r in R])
    module = np.array([r["module"] for r in R])
    strat = np.array([f"{dist_bin(r['domain'], int(r['dist']))}|{r['cpg']}|{gc_bin(float(r['gc']))}" for r in R])
    blk = np.array([f"{r['chrom']}:{int(r['pos']) // args.block_size}" for r in R])
    K = len(CLASSES)
    results = {}
    n_ref = int((cls == REF_CLASS).sum())
    if n_ref < 500:
        log(f"WARNING: only {n_ref} variants in the reference (common) class; intervals will be wide - "
            "add chromosomes (more --gnomad files) or lower --common-af")

    def run(name, mask, y, trend=True):
        if mask.sum() < 200 or y[mask].sum() < 5:
            log(f"  skip {name} (too few variants/calls)")
            return
        cells = Cells(strat[mask], cls[mask], blk[mask], K)
        results[name] = analyse_one(cells, y[mask], args.boot, args.boot_trend, rng, REF_CLASS, trend)
        r = results[name]
        log(f"  {name:32s} O/E singleton {r['oe'][0]:.2f}  common {r['oe'][REF_CLASS]:.2f}"
            + (f"  slope {r['slope']['est']:+.3f} (p={r['slope']['p']:.3f})" if "slope" in r else ""))

    allm = np.ones(n, bool)
    for nm, m in (("all", allm), ("splice", dom == "splice"), ("utr", dom == "utr"), ("promoter", dom == "promoter"),
                  ("canonical_1-2", (dom == "splice") & (dist <= 2))):
        run(f"call | {nm}", m, call)
    for lo, hi, nm in ((3, 5, "3-5"), (6, 10, "6-10"), (11, 20, "11-20"), (21, 50, "21-50"), (51, 100, "51-100")):
        run(f"call | splice d{nm}", (dom == "splice") & (dist >= lo) & (dist <= hi), call, trend=False)
    for thr in (2.8, 4.2, 6.0, 8.0):
        run(f"s_raw>={thr:g} | all", allm, (s_raw >= thr).astype(int))
    for mod in sorted({m for m in module if m and m != "none"}):
        run(f"call & module={mod} | all", allm, ((call == 1) & (module == mod)).astype(int))

    meta = json.loads((out / "exp6_sampling.json").read_text()) if (out / "exp6_sampling.json").exists() else {}
    (out / "exp6_summary.json").write_text(json.dumps(dict(results=results, sampling=meta, n=n, classes=CLASSES,
                                                           build_hash=args.build_hash), indent=2, default=float))
    write_tables(results, out)
    write_report(results, meta, args, out, n)
    if not args.no_figure:
        try:
            make_figure(results, out)
        except Exception as e:                                         # figure is optional
            log(f"figure skipped: {type(e).__name__}: {e}")


def write_tables(results, out):
    rows = []
    for name, r in results.items():
        for k, c in enumerate(CLASSES):
            rows.append(dict(analysis=name, afclass=c, n=int(r["n"][k]), observed=r["obs"][k], expected=r["exp"][k],
                             oe=r["oe"][k], oe_lo=r["oe_ci"][0][k], oe_hi=r["oe_ci"][1][k],
                             mh_or_vs_common=r["mh_or"][k], mh_lo=r["mh_ci"][0][k], mh_hi=r["mh_ci"][1][k]))
    write_csv(out / "exp6_class_table.csv", rows)
    write_csv(out / "exp6_trend.csv", [dict(analysis=n, **r["slope"]) for n, r in results.items() if "slope" in r])


def write_report(results, meta, args, out, n):
    f2 = lambda x: "NA" if x is None or np.isnan(x) else f"{x:.2f}"
    L = ["# Experiment 6: gnomAD allele-frequency spectrum of NC-MutScan v3.3 calls\n",
         f"- Build hash recorded: `{args.build_hash or 'not given'}`; engine `{args.engine}`; flank \u00b1{args.flank} bp; seed {args.seed}.",
         f"- gnomAD files: {', '.join(Path(p).name for p in meta.get('files', []))}; common cut-off AF \u2265 {args.common_af:g}; "
         f"PASS biallelic SNVs, AN \u2265 {args.min_an_frac:.0%} of the 90th-percentile AN; CDS of any isoform excluded.",
         f"- Domains: utr (exon, non-CDS), splice (intronic 1\u2013{args.splice_window} bp from an exon edge), promoter "
         f"(TSS \u2212{args.prom_up}/+{args.prom_down}). Scored sample: {n:,} variants.",
         f"- Strata: domain-distance bin \u00d7 CpG/ts-tv (4) \u00d7 GC bin (4). Intervals: block bootstrap, "
         f"{args.block_size // 1000} kb blocks, {args.boot} replicates ({args.boot_trend} for the slope).\n",
         "## Variants per frequency class\n", "| Class | Definition | Seen in domains | Sampled |", "|---|---|---|---|"]
    for c, d in zip(CLASSES, CLASS_DEF):
        L.append(f"| {c} | {d} | {meta.get('seen', {}).get(c, 'NA')} | {meta.get('kept', {}).get(c, 'NA')} |")
    L.append("\n*O/E = observed calls / calls expected from the stratum-specific call rate (1 = no frequency effect). "
             "MH OR = Mantel\u2013Haenszel odds ratio of a call versus the **common** class. Slope = stratified log-odds "
             "change of a call per frequency-class step (negative = calls fall as frequency rises).*\n")
    for name, r in results.items():
        L.append(f"## {name}\n")
        s = r.get("slope")
        if s:
            L.append(f"Trend: slope {s['est']:+.3f} per class step (95% CI {s['lo']:+.3f} to {s['hi']:+.3f}; "
                     f"bootstrap two-sided P = {s['p']:.3f}).\n")
        L.append("| Class | n | Calls | Expected | O/E (95% CI) | MH OR vs common (95% CI) |\n|---|---|---|---|---|---|")
        for k, c in enumerate(CLASSES):
            mh = "reference" if k == REF_CLASS else f"{f2(r['mh_or'][k])} ({f2(r['mh_ci'][0][k])}\u2013{f2(r['mh_ci'][1][k])})"
            L.append(f"| {c} | {int(r['n'][k]):,} | {int(r['obs'][k]):,} | {r['exp'][k]:.1f} | {f2(r['oe'][k])} "
                     f"({f2(r['oe_ci'][0][k])}\u2013{f2(r['oe_ci'][1][k])}) | {mh} |")
        L.append("")
    L.append("## How to read this\n\n- Classes are sampled at different rates, so O/E is relative to the sampled class mix: "
             "use its pattern across classes, and take effect sizes from the MH OR and the slope.\n- Expected pattern if calls mark functional disruption: O/E > 1 for singletons, falling "
             "below 1 for the common classes, strongest for canonical \u00b11\u20132 bp splice sites and at the higher score "
             "thresholds (dose-response).\n- No frequency gradient is itself informative: it limits the claim that calls mark "
             "selected positions in that domain (for example promoter calls from generic motifs).\n- Remaining confounders: CpG "
             "methylation level, local recombination/background selection and gene-level constraint are not matched; report them "
             "as limitations.")
    (out / "exp6_report.md").write_text("\n".join(L), encoding="utf-8")


def make_figure(results, out):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                         "font.size": 7, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})
    col = {"call | all": "#222222", "call | splice": "#0072B2", "call | utr": "#009E73", "call | promoter": "#E69F00",
           "call | canonical_1-2": "#D55E00"}
    fig, ax = plt.subplots(1, 3, figsize=(183 / 25.4, 62 / 25.4), gridspec_kw=dict(wspace=0.45))
    x = np.arange(len(CLASSES))
    for name, c in col.items():
        if name in results:
            r = results[name]
            oe, lo, hi = np.array(r["oe"]), np.array(r["oe_ci"][0]), np.array(r["oe_ci"][1])
            ax[0].errorbar(x + 0.06 * list(col).index(name), oe, yerr=[oe - lo, hi - oe], fmt="o-", ms=3, lw=0.9,
                           capsize=1.5, color=c, label=name.replace("call | ", ""))
    ax[0].axhline(1, color="grey", lw=0.6, ls="--")
    ax[0].set_xticks(x)
    ax[0].set_xticklabels([c.replace("_", "\n") for c in CLASSES], fontsize=5.5)
    ax[0].set_ylabel("Observed / expected calls")
    ax[0].set_title("Calls by frequency class", fontsize=7)
    ax[0].legend(fontsize=5.5, frameon=False)
    bins = [("splice d3-5", "3\u20135"), ("splice d6-10", "6\u201310"), ("splice d11-20", "11\u201320"),
            ("splice d21-50", "21\u201350"), ("splice d51-100", "51\u2013100")]
    labs, vals, los, his = ["1\u20132"], [], [], []
    names = ["call | canonical_1-2"] + [f"call | {b}" for b, _ in bins]
    labs += [b for _, b in bins]
    xs = []
    for i, nm in enumerate(names):
        if nm in results:
            r = results[nm]
            v, lo, hi = r["mh_or"][0], r["mh_ci"][0][0], r["mh_ci"][1][0]
            xs.append(i)
            vals.append(v)
            los.append(v - lo)
            his.append(hi - v)
    ax[1].errorbar(xs, vals, yerr=[los, his], fmt="o", color="#0072B2", capsize=1.5, ms=3)
    ax[1].axhline(1, color="grey", lw=0.6, ls="--")
    ax[1].set_yscale("log")
    ax[1].set_xticks(range(len(labs)))
    ax[1].set_xticklabels(labs, fontsize=6)
    ax[1].set_xlabel("Distance to exon edge (bp)")
    ax[1].set_ylabel("MH OR, singleton vs common")
    ax[1].set_title("Splice domain", fontsize=7)
    th = ["call | all", "s_raw>=2.8 | all", "s_raw>=4.2 | all", "s_raw>=6 | all", "s_raw>=8 | all"]
    tl = ["call", "\u22652.8", "\u22654.2", "\u22656", "\u22658"]
    for k, c, lab in ((0, "#D55E00", "singleton"), (REF_CLASS, "#0072B2", "common")):
        xx = [i for i, t in enumerate(th) if t in results]
        v = [results[th[i]]["oe"][k] for i in xx]
        lo = [results[th[i]]["oe"][k] - results[th[i]]["oe_ci"][0][k] for i in xx]
        hi = [results[th[i]]["oe_ci"][1][k] - results[th[i]]["oe"][k] for i in xx]
        ax[2].errorbar(xx, v, yerr=[lo, hi], fmt="o-", color=c, capsize=1.5, ms=3, label=lab)
    ax[2].axhline(1, color="grey", lw=0.6, ls="--")
    ax[2].set_xticks(range(len(tl)))
    ax[2].set_xticklabels(tl, fontsize=6)
    ax[2].set_xlabel("S_raw threshold")
    ax[2].set_ylabel("Observed / expected")
    ax[2].set_title("Dose-response", fontsize=7)
    ax[2].legend(fontsize=5.5, frameon=False)
    for a, t in zip(ax, "abc"):
        a.text(-0.2, 1.08, t, transform=a.transAxes, fontsize=9, fontweight="bold")
    fig.savefig(out / "exp6_figure.png", dpi=600, bbox_inches="tight")
    fig.savefig(out / "exp6_figure.pdf", bbox_inches="tight")
    plt.close(fig)


# ===========================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gnomad", nargs="+", help="gnomAD genomes sites VCF(.bgz/.gz), one or more chromosomes")
    ap.add_argument("--gtf", help="GENCODE GTF (file or folder)")
    ap.add_argument("--fasta", help="GRCh38 FASTA (file or folder)")
    ap.add_argument("--engine", default="my_engine:score")
    ap.add_argument("--chroms", default="all", help="restrict to these chromosomes, e.g. 21-22 or 12-22")
    ap.add_argument("--flank", type=int, default=25)
    ap.add_argument("--splice-window", type=int, default=100)
    ap.add_argument("--prom-up", type=int, default=1000)
    ap.add_argument("--prom-down", type=int, default=100)
    ap.add_argument("--common-af", type=float, default=0.01)
    ap.add_argument("--per-class", type=int, default=150000)
    ap.add_argument("--min-an-frac", type=float, default=0.8)
    ap.add_argument("--stage", choices=["all", "sample", "score", "analyse"], default="all")
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--boot-trend", type=int, default=200)
    ap.add_argument("--block-size", type=int, default=100000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--build-hash", default="")
    ap.add_argument("--max-lines", type=int, help="testing: stop after N lines per gnomAD file")
    ap.add_argument("--no-cache", action="store_true")
    ap.add_argument("--no-figure", action="store_true")
    ap.add_argument("--out", default="exp6_out")
    args = ap.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.stage in ("all", "sample"):
        if not (args.gnomad and args.gtf):
            ap.error("--gnomad and --gtf are required for the sample stage")
        stage_sample(args, out)
    if args.stage in ("all", "score"):
        if not args.fasta:
            ap.error("--fasta is required for the score stage")
        stage_score(args, out)
    if args.stage in ("all", "analyse"):
        stage_analyse(args, out)
    log(f"done -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
