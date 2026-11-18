#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
exp1_full_pool_eval.py -- Experiment 1: full-pool held-out evaluation of NC-MutScan v3.3
=========================================================================================
Scores EVERY strictly-noncoding ClinVar variant on the held-out chromosomes (default chr12-22;
about 3,472 pathogenic and 193,490 benign) with the frozen build, changing no parameter, and
reports sensitivity / specificity / PPV / LR+ / AUROC / AUPRC with gene-clustered bootstrap CIs.

Typical run (PowerShell, from the project folder, venv active):

  1) wire the engine once:     copy my_engine_template.py my_engine.py   (then edit it)
  2) smoke test + sanity check on the 150 locked blind variants (expects TP37 FP1 TN74 FN38):
       python exp1_full_pool_eval.py --pool locked_blind_chr12_22.vcf --fasta GRCh38.d1.vd1.fa ^
              --sanity-only --sanity-vcf locked_blind_chr12_22.vcf
  3) the real run with the pool your own build_holdout.py produced (VCF or TSV):
       python exp1_full_pool_eval.py --pool <pool file> --fasta GRCh38.d1.vd1.fa ^
              --sanity-vcf locked_blind_chr12_22.vcf --build-hash b276d0d44045
     or, if you do not have a pool file, an APPROXIMATE one built from ClinVar:
       python exp1_full_pool_eval.py --build-pool clinvar_GRCh38.vcf.gz --save-pool pool_chr12_22.tsv ^
              --fasta GRCh38.d1.vd1.fa --sanity-vcf locked_blind_chr12_22.vcf
     (compare the printed counts with 3,472 / 193,490; use --min-stars to tune the review filter)

Pool file formats: ClinVar-style VCF (label from CLNSIG) or TSV/CSV: chrom,pos,ref,alt,label[,gene]
(label 1 = pathogenic / likely pathogenic, 0 = benign / likely benign).

Outputs (folder --out, default exp1_out): exp1_per_variant.csv (input of experiment 2),
exp1_report.md, exp1_summary.json, exp1_per_chrom.csv, exp1_false_positives.csv.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from nc_exp_common import (Fasta, RankPrep, build_pool_from_clinvar, chrom_order, chrom_set, cluster_bootstrap,
                           cluster_index, confusion, fmt_ci, fmt_pct, load_engine, load_pool, log, lr_plus,
                           make_window, normalise_result, pct_ci, ppv_at_prior, sens_at_spec, wilson, write_csv)

PRIORS = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10]


# ---------------------------------------------------------------------------
def score_variants(rows, fa, engine, flank, limit=None, tag=""):
    out, skipped_ref, skipped_err, t0 = [], 0, 0, time.time()
    for i, r in enumerate(rows):
        if limit and i >= limit:
            break
        if not fa.has(r["chrom"]):
            skipped_ref += 1
            continue
        w = make_window(fa, r["chrom"], r["pos"], r["ref"], r["alt"], flank)
        if w is None:
            skipped_ref += 1
            continue
        wt, mut, center = w
        meta = dict(chrom=r["chrom"], pos=r["pos"], ref=r["ref"], alt=r["alt"], gene=r["gene"], center=center)
        try:
            res = normalise_result(engine(wt, mut, meta))
        except NotImplementedError:
            raise SystemExit("my_engine.py is not wired yet: edit score() (see my_engine_template.py).")
        except Exception as e:                         # keep going, but count and show the first few
            skipped_err += 1
            if skipped_err <= 3:
                log(f"  engine error on {r['id']}: {type(e).__name__}: {e}")
            continue
        out.append({**{k: r[k] for k in ("id", "chrom", "pos", "ref", "alt", "label", "gene")}, **res})
        if (i + 1) % 20000 == 0:
            rate = (i + 1) / (time.time() - t0)
            log(f"  {tag}scored {i + 1:,}/{len(rows):,}  ({rate:,.0f} variants/s)")
    return out, skipped_ref, skipped_err


def sanity_check(args, fa, engine):
    rows = load_pool(args.sanity_vcf)
    log(f"sanity set: {len(rows)} variants from {args.sanity_vcf}")
    scored, sref, serr = score_variants(rows, fa, engine, args.flank, tag="sanity ")
    y = np.array([r["label"] for r in scored])
    c = np.array([r["call"] for r in scored])
    tp, fp, tn, fn = (int(v) for v in confusion(y, c))
    exp = tuple(int(v) for v in args.sanity_expect.split(","))
    log(f"sanity result  TP={tp} FP={fp} TN={tn} FN={fn}   (expected {exp}; "
        f"skipped: ref-mismatch/missing={sref}, errors={serr})")
    if (tp, fp, tn, fn) != exp:
        msg = ("SANITY CHECK FAILED: the adapter does not reproduce the locked blind result.\n"
               "  Likely causes: (1) --flank differs from the flank your evaluation used (try 25, 50, 100);\n"
               "  (2) my_engine.score() maps the wrong keys / threshold; (3) a different build than v3.3;\n"
               "  (4) reference allele mismatches (check the skipped count above).")
        if args.ignore_sanity:
            log("WARNING: " + msg)
        else:
            raise SystemExit(msg + "\n  Re-run with --ignore-sanity to continue anyway.")
    else:
        log("sanity check PASSED (matches the manuscript blind-set confusion matrix)")
    return dict(tp=tp, fp=fp, tn=tn, fn=fn, expected=list(exp), skipped_ref=sref, skipped_err=serr)


# ---------------------------------------------------------------------------
def analyse(rows, args, out):
    rng = np.random.default_rng(args.seed)
    y = np.array([r["label"] for r in rows], int)
    call = np.array([r["call"] for r in rows], int)
    s_raw = np.array([r["s_raw"] for r in rows], float)
    s_uni = np.array([r["s_unified"] for r in rows], float)
    gene = [r["gene"] or f"{r['chrom']}:{r['pos'] // 1_000_000}Mb" for r in rows]
    cl_idx, n_cl = cluster_index(gene)
    n = len(y)
    tp, fp, tn, fn = (int(v) for v in confusion(y, call))
    sens, spec = tp / (tp + fn), tn / (tn + fp)
    prev = (tp + fn) / n
    ppv_nat = tp / (tp + fp) if tp + fp else float("nan")
    lr, lr_lo, lr_hi = lr_plus(tp, fn, fp, tn)
    S = dict(n=n, n_path=tp + fn, n_benign=tn + fp, n_gene_clusters=int(n_cl), prevalence=prev,
             tp=tp, fp=fp, tn=tn, fn=fn, sens=sens, spec=spec, ppv_natural=ppv_nat,
             sens_ci_wilson=wilson(tp, tp + fn), spec_ci_wilson=wilson(tn, tn + fp), lr_plus=[lr, lr_lo, lr_hi],
             ppv_at_prior={str(p): ppv_at_prior(sens, spec, p) for p in PRIORS})

    # ---- gene-clustered bootstrap
    log(f"gene-cluster bootstrap ({args.boot} reps for counts, {args.boot_auc} for AUC/AUPRC; {n_cl:,} clusters) ...")

    def counts_fn(w):
        a, b, c, d = confusion(y, call, w)
        se, sp = a / (a + d), c / (c + b)
        return se, sp, (a / (a + b) if a + b else np.nan), (se / (1 - sp) if sp < 1 else np.nan)

    bc = cluster_bootstrap(cl_idx, n_cl, args.boot, rng, counts_fn)
    S["sens_ci_boot"], S["spec_ci_boot"], S["ppv_natural_ci_boot"], S["lr_plus_ci_boot"] = (
        list(pct_ci(bc[:, i])) for i in range(4))
    S["ppv_at_prior_ci_boot"] = {}
    for p in PRIORS:
        v = bc[:, 0] * p / (bc[:, 0] * p + (1 - bc[:, 1]) * (1 - p))
        S["ppv_at_prior_ci_boot"][str(p)] = list(pct_ci(v))

    prep = {"s_raw": RankPrep(s_raw), "s_unified": RankPrep(s_uni)}
    scores = {"s_raw": s_raw, "s_unified": s_uni}
    S["auc"], S["ap"] = {}, {}
    for k, pr in prep.items():
        a, ap = pr.metrics(y)
        bb = cluster_bootstrap(cl_idx, n_cl, args.boot_auc, rng, lambda w, pr=pr: pr.metrics(y, w))
        S["auc"][k] = [a, *pct_ci(bb[:, 0])]
        S["ap"][k] = [ap, *pct_ci(bb[:, 1])]

    # ---- optional comparator scores (id + numeric columns)
    S["extra"] = {}
    if args.extra_scores:
        ids = [r["id"] for r in rows]
        with open(args.extra_scores, newline="", encoding="utf-8-sig") as fh:
            head = fh.readline()
            fh.seek(0)
            ex = {d["id"]: d for d in csv.DictReader(fh, delimiter="\t" if "\t" in head else ",")}
        cols = [c for c in next(iter(ex.values())).keys() if c != "id"]
        for c in cols:
            idx = [i for i, k in enumerate(ids) if k in ex and ex[k][c] not in ("", "NA", "nan", None)]
            if len(idx) < 50:
                continue
            sc = np.array([float(ex[ids[i]][c]) for i in idx])
            yy, nc = y[idx], s_raw[idx]
            a, ap = RankPrep(sc).metrics(yy)
            a0, ap0 = RankPrep(nc).metrics(yy)
            S["extra"][c] = dict(n=len(idx), auc=a, ap=ap, nc_auc_same_variants=a0, nc_ap_same_variants=ap0,
                                 sens_at_nc_spec=sens_at_spec(sc, yy, spec),
                                 nc_sens=float(np.mean(call[idx][yy == 1])) if (yy == 1).any() else np.nan)

    # ---- per chromosome
    per_chr = []
    for ch in sorted({r["chrom"] for r in rows}, key=chrom_order):
        m = np.array([r["chrom"] == ch for r in rows])
        a, b, c, d = (int(v) for v in confusion(y[m], call[m]))
        per_chr.append(dict(chrom=ch, n_path=a + d, n_benign=b + c, tp=a, fp=b, tn=c, fn=d,
                            sens=a / (a + d) if a + d else np.nan, spec=c / (b + c) if b + c else np.nan,
                            spec_lo=wilson(c, b + c)[0], spec_hi=wilson(c, b + c)[1]))
    write_csv(out / "exp1_per_chrom.csv", per_chr)

    # ---- which module makes the calls?
    mods = {}
    if any(r["module"] for r in rows):
        for r in rows:
            if r["call"]:
                key = r["module"] or "unspecified"
                d = mods.setdefault(key, [0, 0])
                d[0 if r["label"] == 1 else 1] += 1
        S["calls_by_module"] = {k: dict(true_pos=v[0], false_pos=v[1],
                                        precision=v[0] / (v[0] + v[1])) for k, v in mods.items()}
    write_csv(out / "exp1_false_positives.csv",
              [r for r in rows if r["label"] == 0 and r["call"]],
              ["id", "chrom", "pos", "ref", "alt", "gene", "s_raw", "s_unified", "module"])
    return S, per_chr


def write_report(S, per_chr, args, out, sanity, pool_info):
    ci = lambda k: S[k]
    L = []
    A = L.append
    A("# Experiment 1: full-pool held-out evaluation (NC-MutScan v3.3, frozen build)\n")
    A(f"- Build hash recorded: `{args.build_hash or 'not given'}`; engine `{args.engine}`; flank \u00b1{args.flank} bp; "
      f"seed {args.seed}; chromosomes `{args.chroms}`; no parameter was changed.")
    A(f"- Pool: **{S['n_path']:,} pathogenic / {S['n_benign']:,} benign** (prevalence {fmt_pct(S['prevalence'], 3)}), "
      f"{S['n_gene_clusters']:,} gene clusters. Skipped before scoring (reference mismatch / missing contig): "
      f"{pool_info['skipped_ref']:,}; engine errors: {pool_info['skipped_err']:,}.")
    if sanity:
        A(f"- Sanity check on the locked blind set: TP {sanity['tp']} FP {sanity['fp']} TN {sanity['tn']} FN {sanity['fn']} "
          f"(expected {tuple(sanity['expected'])}).")
    A("\n## Confusion matrix and rates\n")
    A(f"| | Called (rule-in) | Not called |\n|---|---|---|\n| Pathogenic | TP {S['tp']:,} | FN {S['fn']:,} |\n"
      f"| Benign | FP {S['fp']:,} | TN {S['tn']:,} |\n")
    sb, pb = S["sens_ci_boot"], S["spec_ci_boot"]
    A("| Metric | Estimate | Wilson 95% CI | Gene-clustered bootstrap 95% CI |\n|---|---|---|---|")
    A(f"| Sensitivity | {fmt_pct(S['sens'])} | {fmt_pct(S['sens_ci_wilson'][0])}\u2013{fmt_pct(S['sens_ci_wilson'][1])} "
      f"| {fmt_pct(sb[0])}\u2013{fmt_pct(sb[1])} |")
    A(f"| Specificity | {fmt_pct(S['spec'], 3)} | {fmt_pct(S['spec_ci_wilson'][0], 3)}\u2013{fmt_pct(S['spec_ci_wilson'][1], 3)} "
      f"| {fmt_pct(pb[0], 3)}\u2013{fmt_pct(pb[1], 3)} |")
    lr = S["lr_plus"]
    A(f"| LR+ | {lr[0]:.2f} | {lr[1]:.2f}\u2013{lr[2]:.2f} (Katz) | {S['lr_plus_ci_boot'][0]:.2f}\u2013{S['lr_plus_ci_boot'][1]:.2f} |")
    A(f"| PPV at pool prevalence ({fmt_pct(S['prevalence'], 2)}) | {fmt_pct(S['ppv_natural'])} | | "
      f"{fmt_pct(S['ppv_natural_ci_boot'][0])}\u2013{fmt_pct(S['ppv_natural_ci_boot'][1])} |")
    for p in PRIORS:
        c = S["ppv_at_prior_ci_boot"][str(p)]
        A(f"| PPV at prior {100 * p:g}% | {fmt_pct(S['ppv_at_prior'][str(p)])} | | {fmt_pct(c[0])}\u2013{fmt_pct(c[1])} |")
    A("\n*Earlier estimates for comparison: blind set (n=150) specificity 98.67% (92.8\u201399.8), sensitivity 49.33% "
      "(38.3\u201360.4); 23-chromosome LOCO (n=300) specificity 94.00% (89.0\u201396.8), sensitivity 44.67%.*\n")
    A("## Threshold-free discrimination (gene-clustered bootstrap 95% CI)\n")
    A(f"AUPRC chance level = prevalence = {S['prevalence']:.4f}.\n")
    A("| Score | ROC-AUC | AUPRC |\n|---|---|---|")
    for k in ("s_raw", "s_unified"):
        a, p = S["auc"][k], S["ap"][k]
        A(f"| {k} | {a[0]:.3f} ({a[1]:.3f}\u2013{a[2]:.3f}) | {p[0]:.4f} ({p[1]:.4f}\u2013{p[2]:.4f}) |")
    if S.get("extra"):
        A("\n## Comparator scores (same variants only)\n")
        A("| Score | n | ROC-AUC | AUPRC | Sens. at NC-MutScan specificity | NC-MutScan sens. | NC S_raw ROC-AUC on same variants |\n|---|---|---|---|---|---|---|")
        for k, v in S["extra"].items():
            A(f"| {k} | {v['n']:,} | {v['auc']:.3f} | {v['ap']:.4f} | {fmt_pct(v['sens_at_nc_spec'])} | "
              f"{fmt_pct(v['nc_sens'])} | {v['nc_auc_same_variants']:.3f} |")
    if S.get("calls_by_module"):
        A("\n## Calls by module\n\n| Module | True positives | False positives | Precision |\n|---|---|---|---|")
        for k, v in S["calls_by_module"].items():
            A(f"| {k} | {v['true_pos']:,} | {v['false_pos']:,} | {fmt_pct(v['precision'])} |")
    A("\n## Per chromosome\n\n| Chr | Pathogenic | Benign | Sensitivity | Specificity (Wilson) |\n|---|---|---|---|---|")
    for r in per_chr:
        A(f"| {r['chrom']} | {r['n_path']:,} | {r['n_benign']:,} | {fmt_pct(r['sens'])} | "
          f"{fmt_pct(r['spec'], 3)} ({fmt_pct(r['spec_lo'], 3)}\u2013{fmt_pct(r['spec_hi'], 3)}) |")
    A("\n*Note: ClinVar-wide prevalence is not the clinical prevalence; use the prior-adjusted PPV rows for "
      "screening settings. Variants cluster within genes, so the clustered intervals are the ones to report.*")
    (out / "exp1_report.md").write_text("\n".join(L), encoding="utf-8")


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--pool", help="pool VCF/TSV (your own build_holdout.py output)")
    g.add_argument("--build-pool", metavar="CLINVAR_VCF", help="build an APPROXIMATE noncoding pool from ClinVar")
    ap.add_argument("--save-pool", help="write the built pool to this TSV")
    ap.add_argument("--min-stars", type=int, default=0, help="min ClinVar review stars for --build-pool (0 = none)")
    ap.add_argument("--fasta", required=True, help="GRCh38 FASTA (file or the folder containing it)")
    ap.add_argument("--engine", default="my_engine:score")
    ap.add_argument("--chroms", default="12-22")
    ap.add_argument("--flank", type=int, default=25)
    ap.add_argument("--max-len", type=int, default=15)
    ap.add_argument("--sanity-vcf", help="locked_blind_chr12_22.vcf")
    ap.add_argument("--sanity-expect", default="37,1,74,38", help="TP,FP,TN,FN expected on the sanity set")
    ap.add_argument("--sanity-only", action="store_true")
    ap.add_argument("--ignore-sanity", action="store_true")
    ap.add_argument("--extra-scores", help="TSV/CSV: id (chrom:pos:ref:alt) + numeric score columns, e.g. spliceai")
    ap.add_argument("--build-hash", default="")
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--boot-auc", type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--limit", type=int, help="score only the first N variants (testing)")
    ap.add_argument("--out", default="exp1_out")
    args = ap.parse_args(argv)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    chroms = chrom_set(args.chroms)

    fa = Fasta(args.fasta)
    engine = load_engine(args.engine)
    sanity = sanity_check(args, fa, engine) if args.sanity_vcf else None
    if args.sanity_only:
        return 0

    if args.build_pool:
        log("building pool from ClinVar (strictly noncoding, <=%d bp, chroms %s) ..." % (args.max_len, args.chroms))
        pool = build_pool_from_clinvar(args.build_pool, chroms, args.max_len, args.min_stars)
        if args.save_pool:
            write_csv(args.save_pool, [{k: r[k] for k in ("chrom", "pos", "ref", "alt", "label", "gene")} for r in pool],
                      ["chrom", "pos", "ref", "alt", "label", "gene"])
    elif args.pool:
        pool = load_pool(args.pool, chroms, args.max_len)
    else:
        ap.error("give --pool or --build-pool")
    npath = sum(r["label"] == 1 for r in pool)
    log(f"pool: {npath:,} pathogenic, {len(pool) - npath:,} benign (manuscript reference for chr12-22: 3,472 / 193,490)")

    t0 = time.time()
    rows, sref, serr = score_variants(pool, fa, engine, args.flank, args.limit)
    log(f"scored {len(rows):,} variants in {time.time() - t0:,.0f} s (ref mismatch/missing {sref:,}, errors {serr:,})")
    if not rows:
        raise SystemExit("No variants were scored.")
    feat_cols = sorted({k for r in rows for k in r if k.startswith("feat_")})
    cols = ["id", "chrom", "pos", "ref", "alt", "label", "gene", "s_raw", "s_unified", "call", "module"] + feat_cols
    write_csv(out / "exp1_per_variant.csv", rows, cols)

    S, per_chr = analyse(rows, args, out)
    S["run"] = dict(build_hash=args.build_hash, engine=args.engine, flank=args.flank, chroms=args.chroms,
                    seed=args.seed, sanity=sanity, skipped_ref=sref, skipped_err=serr)
    (out / "exp1_summary.json").write_text(json.dumps(S, indent=2, default=float), encoding="utf-8")
    write_report(S, per_chr, args, out, sanity, dict(skipped_ref=sref, skipped_err=serr))
    log(f"done -> {out}/exp1_report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
