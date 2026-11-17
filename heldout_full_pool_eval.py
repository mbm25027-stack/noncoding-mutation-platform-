#!/usr/bin/env python3
"""
Full-pool held-out evaluation of the FROZEN NC-MutScan v3.3 engine on the
chr12-22 ClinVar pool (check_blind/pool_all.vcf).  No parameters are changed.

Dependencies: numpy + standard library only. Run from the folder with server.py.

Pool definition is explicit and recorded in report.json:
  --min-stars     minimum ClinVar review-status stars (default 0)
  --exclude-mc    molecular-consequence terms removed from the pool
                  (default removes coding/canonical-splice classes that leaked
                  into pool_all.vcf: nonsense, splice donor/acceptor,
                  inframe_indel, no_sequence_alteration)

Smoke test:  python heldout_full_pool_eval.py --pool check_blind\\pool_all.vcf --fasta $fa --smoke 200
Full run:    python heldout_full_pool_eval.py --pool check_blind\\pool_all.vcf --fasta $fa --out eval_full_pool
Scores are cached in <out>/scores.tsv, so re-running with the same --out and
pool definition skips scoring and only redoes the statistics.
"""
import argparse, json, math, os, random, sys, time
from collections import Counter
import numpy as np

FLANK_BP = 25
DEFAULT_EXCLUDE_MC = "nonsense,splice_donor_variant,splice_acceptor_variant,inframe_indel,no_sequence_alteration"

STARS = {"reviewed_by_expert_panel": 3,
         "criteria_provided,_multiple_submitters,_no_conflicts": 2,
         "criteria_provided,_single_submitter": 1}


# ---------------------------------------------------------------- FASTA (.fai)
class Fasta:
    def __init__(self, path):
        self.path = path
        fai = path + ".fai"
        if not os.path.exists(fai):
            print(f"Building index {fai} (one-time, may take a few minutes)...")
            self._build(fai)
        self.idx = {}
        with open(fai) as f:
            for line in f:
                n, L, off, lb, lw = line.rstrip("\n").split("\t")[:5]
                self.idx[n] = (int(L), int(off), int(lb), int(lw))
        self.fh = open(path, "rb")

    def _build(self, fai_path):
        ents, name, offset = [], None, 0
        length = seq_off = lb = lw = 0
        with open(self.path, "rb") as f:
            for line in f:
                if line.startswith(b">"):
                    if name:
                        ents.append((name, length, seq_off, lb, lw))
                    name = line[1:].split()[0].decode()
                    length, seq_off, lb, lw = 0, offset + len(line), None, None
                else:
                    s = line.rstrip(b"\r\n")
                    if lb is None:
                        lb, lw = len(s), len(line)
                    length += len(s)
                offset += len(line)
            if name:
                ents.append((name, length, seq_off, lb, lw))
        with open(fai_path, "w") as o:
            for e in ents:
                o.write("\t".join(map(str, e)) + "\n")

    def contig(self, chrom):
        for n in (f"chr{chrom}", chrom):
            if n in self.idx:
                return n
        return None

    def fetch(self, name, start, end):  # 0-based, half-open
        L, off, lb, lw = self.idx[name]
        if start < 0 or end > L:
            return None
        first, last = start // lb, (end - 1) // lb
        b0 = off + first * lw + start % lb
        b1 = off + last * lw + (end - 1) % lb + 1
        self.fh.seek(b0)
        return self.fh.read(b1 - b0).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


# ---------------------------------------------------------------- pool
def load_pool(path, min_stars=0, exclude_mc=()):
    rows, skipped, excl = [], Counter(), set(exclude_mc)
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.rstrip("\r\n").split("\t")
            chrom, pos, ref, alt, info_s = t[0].replace("chr", ""), int(t[1]), t[3].upper(), t[4].upper(), t[7]
            info = dict(kv.split("=", 1) for kv in info_s.split(";") if "=" in kv)
            sig = info.get("CLNSIG", "").lower()
            if "pathogenic" in sig and "benign" not in sig:
                lab = 1
            elif "benign" in sig and "pathogenic" not in sig:
                lab = 0
            else:
                skipped["unclassifiable_CLNSIG"] += 1
                continue
            if STARS.get(info.get("CLNREVSTAT", ""), 0) < min_stars:
                skipped["review_status"] += 1
                continue
            terms = {x.split("|")[-1] for x in info.get("MC", "").split(",") if x}
            if terms & excl:
                skipped["excluded_mc"] += 1
                continue
            gi = info.get("GENEINFO", "")
            gene = gi.split("|")[0].split(":")[0] if gi else ""
            rows.append((chrom, pos, ref, alt, lab, gene))
    return rows, skipped


# ---------------------------------------------------------------- frozen engine
def is_frozen_driver_call(res):   # verbatim from build_dev_flanks_and_loco.py
    cls_txt = str(res.get("classification", "")).lower()
    if any(neg in cls_txt for neg in ("benign", "neutral", "polymorphism", "passenger", "sub-threshold", "subthreshold", "non-driver", "compet")):
        return 0
    if any(pos in cls_txt for pos in ("driver", "pathogenic", "loss", "gain", "splice", "polya", "disrupt")):
        return 1
    return 0

def load_engine(server_dir):
    sys.path.insert(0, os.path.abspath(server_dir))
    import server
    def score(wt, mut):
        res = server.run_v33_prediction({"wt_seq": wt, "mut_seq": mut, "flank_bp": FLANK_BP})
        if not res:
            return float("nan"), 0, None
        if "S_raw" in res:
            s = float(res["S_raw"])
        elif "score" in res:
            s = float(res["score"])
        else:
            raise KeyError(f"No S_raw/score in engine output; keys = {sorted(res)}")
        return s, is_frozen_driver_call(res), res
    return score

def build_flanks(fa, chrom, pos, ref, alt):
    name = fa.contig(chrom)
    if name is None:
        return None, "no_contig"
    start0 = pos - 1 - FLANK_BP
    end0 = pos - 1 + len(ref) + FLANK_BP
    dna = fa.fetch(name, start0, end0)
    if dna is None or len(dna) != end0 - start0:
        return None, "flank_out_of_range"
    if dna[FLANK_BP:FLANK_BP + len(ref)] != ref:
        return None, "ref_mismatch"
    return (dna, dna[:FLANK_BP] + alt + dna[FLANK_BP + len(ref):]), None


# ---------------------------------------------------------------- metrics
class Ranked:
    """Pre-sorted by score so weighted (bootstrap) metrics are one vectorised pass."""
    def __init__(self, score, y, call):
        o = np.argsort(-score, kind="stable")
        self.order = o
        s = score[o]
        self.y = y[o].astype(float)
        self.pred = call[o].astype(bool)
        self.ends = np.r_[np.nonzero(np.diff(s))[0], len(s) - 1]   # tie-group ends

    def metrics(self, w):
        w = w[self.order]
        ctp = np.cumsum(w * self.y)[self.ends]
        cfp = np.cumsum(w * (1 - self.y))[self.ends]
        P, N = ctp[-1], cfp[-1]

        den = ctp + cfp
        rec = ctp / P
        # groups with zero total weight so far have undefined precision, but
        # their recall increment is 0, so they contribute nothing to AP
        prec = np.divide(ctp, den, out=np.zeros_like(den), where=den > 0)
        ap = float(np.sum(np.diff(np.r_[0.0, rec]) * prec))

        tpr, fpr = np.r_[0.0, ctp / P], np.r_[0.0, cfp / N]
        auroc = float(np.sum(np.diff(fpr) * (tpr[1:] + tpr[:-1]) / 2))

        tp = float(np.sum(w * self.y * self.pred))
        fp = float(np.sum(w * (1 - self.y) * self.pred))
        fn, tn = P - tp, N - fp
        return dict(auprc=ap, auroc=auroc, prevalence=P / (P + N),
                    sensitivity=tp / P, specificity=tn / N,
                    ppv=tp / (tp + fp) if tp + fp > 0 else float("nan"),
                    npv=tn / (tn + fn) if tn + fn > 0 else float("nan"),
                    tp=tp, fp=fp, fn=fn, tn=tn)

def _binom_cdf(k, n, p):
    if p <= 0: return 1.0
    if p >= 1: return 1.0 if k >= n else 0.0
    lp, lq, lg = math.log(p), math.log1p(-p), math.lgamma
    terms = [lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq for i in range(k + 1)]
    m = max(terms)
    return min(1.0, math.exp(m) * sum(math.exp(t - m) for t in terms))

def clopper_pearson(k, n, alpha=0.05):
    def bis(f, lo, hi):
        for _ in range(80):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if f(mid) else (lo, mid)
        return (lo + hi) / 2
    lo = 0.0 if k == 0 else bis(lambda p: 1 - _binom_cdf(k - 1, n, p) < alpha / 2, 0.0, 1.0)
    hi = 1.0 if k == n else bis(lambda p: _binom_cdf(k, n, p) > alpha / 2, 0.0, 1.0)
    return lo, hi


# ---------------------------------------------------------------- scoring loop
def score_all(rows, fa, score_fn, ckpt):
    done = {}
    if os.path.exists(ckpt):
        with open(ckpt) as f:
            next(f, None)
            for line in f:
                k, s, c, why = line.rstrip("\n").split("\t")
                done[k] = (float(s), int(c), why)
        print(f"Resuming: {len(done):,} variants already in {ckpt}")
    out = open(ckpt, "a")
    if not done and os.path.getsize(ckpt) == 0:
        out.write("key\tscore\tcall\tstatus\n")
    res, t0, n_new = [], time.time(), 0
    for i, (chrom, pos, ref, alt, lab, gene) in enumerate(rows):
        k = f"{chrom}:{pos}:{ref}:{alt}"
        if k not in done:
            fl, why = build_flanks(fa, chrom, pos, ref, alt)
            if fl is None:
                rec = (float("nan"), 0, why)
            else:
                try:
                    s, c, _ = score_fn(*fl)
                    rec = (s, c, "ok" if not math.isnan(s) else "engine_none")
                except Exception as e:
                    rec = (float("nan"), 0, "engine_error:" + type(e).__name__)
                    if n_new < 5:
                        print(f"  engine error at {k}: {e}", file=sys.stderr)
            out.write(f"{k}\t{rec[0]}\t{rec[1]}\t{rec[2]}\n")
            done[k] = rec
            n_new += 1
            if n_new % 5000 == 0:
                out.flush()
                print(f"  {i+1:,}/{len(rows):,}  {n_new/(time.time()-t0):,.0f} variants/s")
        res.append(done[k])
    out.close()
    return res


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--server-dir", default=".")
    ap.add_argument("--out", default="eval_full_pool")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--smoke", type=int, default=0, help="score N random variants and exit")
    ap.add_argument("--min-stars", type=int, default=0)
    ap.add_argument("--exclude-mc", default=DEFAULT_EXCLUDE_MC,
                    help="comma-separated molecular-consequence terms to remove from the pool")
    a = ap.parse_args()

    exclude_mc = [x for x in a.exclude_mc.split(",") if x]
    rows, skipped = load_pool(a.pool, a.min_stars, exclude_mc)
    labels = np.array([r[4] for r in rows])
    n_p, n_b = int((labels == 1).sum()), int((labels == 0).sum())
    print(f"Pool (min_stars={a.min_stars}, excluded MC={exclude_mc}):")
    print(f"  {len(rows):,} variants | pathogenic {n_p:,} | benign {n_b:,} | "
          f"ClinVar-pool prevalence {n_p/len(rows):.2%} | removed {dict(skipped)}")
    fa = Fasta(a.fasta)
    score_fn = load_engine(a.server_dir)

    if a.smoke:
        random.seed(1)
        t0, stat = time.time(), Counter()
        for r in random.sample(rows, a.smoke):
            fl, why = build_flanks(fa, r[0], r[1], r[2], r[3])
            if fl is None:
                stat[why] += 1
                continue
            s, c, res = score_fn(*fl)
            if stat["ok"] == 0:
                print("Engine output keys:", sorted(res) if res else res)
                print("Example:", {k: res[k] for k in res if k in ("classification", "S_raw", "score", "S_unified", "is_driver")})
            stat["ok"] += 1
            stat[f"call={c}"] += 1
        dt = time.time() - t0
        print(dict(stat), f"| {a.smoke/dt:,.0f} variants/s")
        return

    os.makedirs(a.out, exist_ok=True)
    res = score_all(rows, fa, score_fn, os.path.join(a.out, "scores.tsv"))
    status = Counter(r[2] for r in res)
    print("Scoring status:", dict(status))
    keep = np.array([r[2] == "ok" for r in res])
    score = np.array([r[0] for r in res])[keep]
    call = np.array([r[1] for r in res])[keep]
    y = labels[keep]
    kept_rows = [r for r, k in zip(rows, keep) if k]
    print(f"Evaluated {len(y):,} variants ({int((y==1).sum()):,} P / {int((y==0).sum()):,} B); "
          f"excluded {int((~keep).sum()):,} (see status counts; report these).")

    zero_p = float(np.mean(score[y == 1] == 0.0))
    zero_b = float(np.mean(score[y == 0] == 0.0))
    print(f"Fraction with S_raw == 0 (tied at the bottom): pathogenic {zero_p:.1%}, benign {zero_b:.1%}")

    names = [r[5] if r[5] else f"nogene_{r[0]}_{r[1] // 1_000_000}" for r in kept_rows]
    uniq = {}
    cidx = np.array([uniq.setdefault(n, len(uniq)) for n in names])
    nc = len(uniq)
    print(f"Clusters (genes / 1Mb bins for gene-less variants): {nc:,}")

    ranked = Ranked(score, y, call)
    point = ranked.metrics(np.ones(len(y)))

    rng = np.random.default_rng(a.seed)
    draws = []
    for b in range(a.n_boot):
        cnt = np.bincount(rng.integers(0, nc, nc), minlength=nc)
        draws.append(ranked.metrics(cnt[cidx].astype(float)))
        if (b + 1) % 500 == 0:
            print(f"  bootstrap {b+1}/{a.n_boot}")
    boot = {k: np.array([d[k] for d in draws]) for k in draws[0]}
    n_nan = {k: int(np.isnan(boot[k]).sum()) for k in
             ("auprc", "auroc", "sensitivity", "specificity", "ppv", "npv")}
    if any(n_nan.values()):
        print("NOTE: bootstrap draws with undefined values (dropped from CIs):", n_nan)

    rep = {"pool_definition": {"source_file": a.pool, "min_stars": a.min_stars,
                               "excluded_mc_terms": exclude_mc, "removed_counts": dict(skipped),
                               "n_pathogenic_in_pool": n_p, "n_benign_in_pool": n_b},
           "n_pathogenic_evaluated": int((y == 1).sum()), "n_benign_evaluated": int((y == 0).sum()),
           "scoring_status": dict(status), "n_clusters": nc, "n_boot": a.n_boot,
           "n_bootstrap_draws_with_nan": n_nan,
           "pool_prevalence": point["prevalence"], "score_used_for_ranking": "S_raw",
           "call_used_for_threshold_metrics": "is_frozen_driver_call",
           "frac_S_raw_zero": {"pathogenic": zero_p, "benign": zero_b}, "metrics": {}}
    for k in ["auprc", "auroc", "sensitivity", "specificity", "ppv", "npv"]:
        lo, hi = np.nanpercentile(boot[k], [2.5, 97.5])
        rep["metrics"][k] = {"estimate": point[k], "ci95_gene_clustered": [float(lo), float(hi)]}
    rep["confusion"] = {k: int(point[k]) for k in ("tp", "fp", "fn", "tn")}
    nben = int(point["fp"] + point["tn"])
    cp = clopper_pearson(int(point["fp"]), nben)
    rep["fpr_point"] = point["fp"] / nben
    rep["fpr_clopper_pearson_95_(ignores_clustering)"] = list(cp)
    json.dump(rep, open(os.path.join(a.out, "report.json"), "w"), indent=2)
    np.savez(os.path.join(a.out, "bootstrap_draws.npz"), **boot)

    print(f"\n=== FULL-POOL HELD-OUT (chr12-22), ClinVar-pool prevalence {point['prevalence']:.2%} ===")
    for k, v in rep["metrics"].items():
        print(f"{k:12s} {v['estimate']:.4f}  95% CI [{v['ci95_gene_clustered'][0]:.4f}, {v['ci95_gene_clustered'][1]:.4f}]")
    print(f"AUPRC of a random scorer = prevalence = {point['prevalence']:.4f}")
    print("TP/FP/FN/TN:", rep["confusion"])
    print(f"FPR {rep['fpr_point']:.5f}  exact 95% CI [{cp[0]:.5f}, {cp[1]:.5f}]")
    print(f"Outputs in {a.out}/")


if __name__ == "__main__":
    main()