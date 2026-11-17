#!/usr/bin/env python3
"""
Full-pool held-out evaluation (chr12-22 ClinVar pool) -- NO parameter changes.

Steps
  1. Load pool (3,472 pathogenic + 193,490 benign expected).
  2. Score every variant with the frozen scorer (resumable checkpoint file).
  3. Report AUROC, AUPRC, and sensitivity / specificity / PPV at the frozen
     threshold, at natural prevalence.
  4. Gene-clustered bootstrap CIs (resample genes, keep all their variants).
  5. Exact (Clopper-Pearson) CI for specificity / FPR.

Usage
  python full_pool_heldout_eval.py --pool holdout_v34/pool.tsv \
      --scorer server:score_variant --threshold 0.5 --out eval_full_pool

Pool file (TSV with header, or .vcf/.vcf.gz from ClinVar):
  TSV columns: chrom  pos  ref  alt  label  [gene]      (label: 1=pathogenic, 0=benign)
  VCF: label derived from CLNSIG, gene from GENEINFO.
  IMPORTANT: use the SAME pool-building filters as build_holdout.py. Prefer
  pointing --pool at the file that script already produced.

Scorer adapter (--scorer module:function)
  function(chrom: str, pos: int, ref: str, alt: str) -> float   (higher = more pathogenic)
  Wrap whatever server.py uses (genome FASTA fetch, engines) in such a function.
"""
import argparse, gzip, importlib, json, os, sys, time
import numpy as np
import pandas as pd

# ----------------------------------------------------------------- pool loading
def _gene_from_geneinfo(s):
    return s.split("|")[0].split(":")[0] if s else ""

def load_pool_vcf(path, chroms):
    opener = gzip.open if path.endswith(".gz") else open
    rows = []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            chrom = c[0].replace("chr", "")
            if chrom not in chroms:
                continue
            info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
            sig = info.get("CLNSIG", "")
            if "Conflicting" in sig or "Uncertain" in sig:
                continue
            if "athogenic" in sig and "Benign" not in sig:
                lab = 1
            elif "enign" in sig and "athogenic" not in sig:
                lab = 0
            else:
                continue
            if len(c[3]) != 1 or len(c[4]) != 1:      # SNVs only (match your pool!)
                continue
            rows.append((chrom, int(c[1]), c[3], c[4], lab,
                         _gene_from_geneinfo(info.get("GENEINFO", ""))))
    return pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt", "label", "gene"])

def load_pool(path):
    chroms = {str(i) for i in range(12, 23)}
    if path.endswith((".vcf", ".vcf.gz")):
        df = load_pool_vcf(path, chroms)
    else:
        df = pd.read_csv(path, sep="\t", dtype={"chrom": str})
        df["chrom"] = df["chrom"].str.replace("chr", "", regex=False)
        if "gene" not in df.columns:
            df["gene"] = ""
    df["gene"] = df["gene"].fillna("").astype(str)
    return df.reset_index(drop=True)

# ----------------------------------------------------------------- scoring
def load_scorer(spec):
    mod, fn = spec.split(":")
    sys.path.insert(0, os.getcwd())
    return getattr(importlib.import_module(mod), fn)

def score_pool(df, scorer, ckpt, add_chr):
    key = lambda r: f"{r.chrom}:{r.pos}:{r.ref}:{r.alt}"
    done = {}
    if os.path.exists(ckpt):
        old = pd.read_csv(ckpt, sep="\t", dtype={"key": str})
        done = dict(zip(old["key"], old["score"]))
        print(f"Resuming: {len(done):,} scores already in {ckpt}")
    new = open(ckpt, "a")
    if not done:
        new.write("key\tscore\n")
    t0, n_new, scores = time.time(), 0, np.empty(len(df))
    for i, r in enumerate(df.itertuples(index=False)):
        k = key(r)
        if k in done:
            scores[i] = done[k]
            continue
        chrom = ("chr" + r.chrom) if add_chr else r.chrom
        try:
            s = float(scorer(chrom, int(r.pos), r.ref, r.alt))
        except Exception as e:                    # record, do not silently drop
            print(f"  scorer failed at {k}: {e}", file=sys.stderr)
            s = np.nan
        scores[i] = s
        new.write(f"{k}\t{s}\n")
        n_new += 1
        if n_new % 5000 == 0:
            new.flush()
            rate = n_new / (time.time() - t0)
            print(f"  {i+1:,}/{len(df):,}  {rate:,.0f} variants/s")
    new.close()
    return scores

# ----------------------------------------------------------------- metrics
class Ranked:
    """Pre-sorted structure so weighted (bootstrap) metrics are O(n) vectorised."""
    def __init__(self, score, y, thr):
        o = np.argsort(-score, kind="stable")
        self.order = o
        s = score[o]
        self.y = y[o].astype(float)
        self.ends = np.r_[np.nonzero(np.diff(s))[0], len(s) - 1]   # tie-group ends
        self.pred = s >= thr

    def metrics(self, w):
        w = w[self.order]
        ctp = np.cumsum(w * self.y)[self.ends]
        cfp = np.cumsum(w * (1 - self.y))[self.ends]
        P, N = ctp[-1], cfp[-1]
        rec = ctp / P
        prec = ctp / (ctp + cfp)
        ap = np.sum(np.diff(np.r_[0.0, rec]) * prec)
        auroc = np.trapz(np.r_[0.0, ctp / P], np.r_[0.0, cfp / N])
        tp = np.sum(w * self.y * self.pred); fp = np.sum(w * (1 - self.y) * self.pred)
        fn = P - tp; tn = N - fp
        return dict(
            auprc=ap, auroc=auroc, prevalence=P / (P + N),
            sensitivity=tp / P, specificity=tn / N,
            ppv=tp / (tp + fp) if tp + fp > 0 else np.nan,
            npv=tn / (tn + fn) if tn + fn > 0 else np.nan,
            tp=tp, fp=fp, fn=fn, tn=tn)

def clopper_pearson(k, n, alpha=0.05):
    from scipy.stats import beta
    lo = 0.0 if k == 0 else beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - alpha / 2, k + 1, n - k)
    return lo, hi

def cluster_bootstrap(ranked, cluster_idx, n_clusters, n_boot, seed):
    rng = np.random.default_rng(seed)
    out = []
    for b in range(n_boot):
        counts = np.bincount(rng.integers(0, n_clusters, n_clusters), minlength=n_clusters)
        m = ranked.metrics(counts[cluster_idx].astype(float))
        out.append(m)
        if (b + 1) % 500 == 0:
            print(f"  bootstrap {b+1}/{n_boot}")
    return pd.DataFrame(out)

# ----------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)
    ap.add_argument("--scorer", required=True, help="module:function")
    ap.add_argument("--threshold", type=float, required=True,
                    help="FROZEN decision threshold (score >= thr => pathogenic)")
    ap.add_argument("--out", default="eval_full_pool")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--add-chr", action="store_true", help="pass 'chr12' instead of '12'")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    df = load_pool(a.pool)
    n_p, n_b = int((df.label == 1).sum()), int((df.label == 0).sum())
    print(f"Pool: {len(df):,} variants | pathogenic {n_p:,} | benign {n_b:,}")
    if (n_p, n_b) != (3472, 193490):
        print("WARNING: pool differs from the expected 3,472 / 193,490 -- check filters.")

    scorer = load_scorer(a.scorer)
    df["score"] = score_pool(df, scorer, os.path.join(a.out, "scores.tsv"), a.add_chr)
    n_fail = int(df.score.isna().sum())
    if n_fail:
        print(f"WARNING: {n_fail:,} variants failed scoring and are EXCLUDED; report this.")
    df = df.dropna(subset=["score"]).reset_index(drop=True)
    df.to_csv(os.path.join(a.out, "scored_pool.tsv"), sep="\t", index=False)

    # clusters: gene symbol; variants without a gene -> 1 Mb genomic bin
    cl = np.where(df.gene != "", df.gene,
                  "nogene_" + df.chrom + "_" + (df.pos // 1_000_000).astype(str))
    cluster_idx, uniq = pd.factorize(cl)
    n_clusters = len(uniq)
    print(f"Clusters (genes/bins): {n_clusters:,}")

    y = df.label.to_numpy()
    ranked = Ranked(df.score.to_numpy(), y, a.threshold)
    point = ranked.metrics(np.ones(len(df)))

    try:                                           # sanity check vs sklearn
        from sklearn.metrics import average_precision_score, roc_auc_score
        assert abs(point["auprc"] - average_precision_score(y, df.score)) < 1e-6
        assert abs(point["auroc"] - roc_auc_score(y, df.score)) < 1e-6
    except ImportError:
        pass

    print(f"Gene-clustered bootstrap ({a.n_boot} reps)...")
    boot = cluster_bootstrap(ranked, cluster_idx, n_clusters, a.n_boot, a.seed)
    boot.to_csv(os.path.join(a.out, "bootstrap_draws.csv"), index=False)

    report = {"n_pathogenic": n_p, "n_benign": n_b, "n_scored": len(df),
              "n_scoring_failures": n_fail, "n_clusters": int(n_clusters),
              "threshold": a.threshold, "n_boot": a.n_boot, "metrics": {}}
    for k in ["auprc", "auroc", "sensitivity", "specificity", "ppv", "npv"]:
        lo, hi = np.nanpercentile(boot[k], [2.5, 97.5])
        report["metrics"][k] = {"estimate": float(point[k]),
                                "ci95_gene_clustered": [float(lo), float(hi)]}
    report["prevalence"] = float(point["prevalence"])
    report["auprc_baseline_(=prevalence)"] = float(point["prevalence"])
    report["confusion"] = {k: int(point[k]) for k in ["tp", "fp", "fn", "tn"]}

    # exact binomial CI on specificity (ignores clustering; the bootstrap is primary)
    lo, hi = clopper_pearson(int(point["fp"]), n_b)
    report["false_positive_rate_clopper_pearson_95"] = [lo, hi]
    report["fpr_point"] = float(point["fp"] / n_b)

    with open(os.path.join(a.out, "report.json"), "w") as f:
        json.dump(report, f, indent=2)

    print("\n=== FULL-POOL HELD-OUT RESULTS (chr12-22, natural prevalence "
          f"{point['prevalence']:.4%}) ===")
    for k, v in report["metrics"].items():
        print(f"{k:12s} {v['estimate']:.4f}   95% CI [{v['ci95_gene_clustered'][0]:.4f}, "
              f"{v['ci95_gene_clustered'][1]:.4f}]")
    print(f"AUPRC baseline (random) = {point['prevalence']:.4f}")
    print(f"TP={report['confusion']['tp']} FP={report['confusion']['fp']} "
          f"FN={report['confusion']['fn']} TN={report['confusion']['tn']}")
    print(f"FPR {report['fpr_point']:.5f}, exact 95% CI [{lo:.5f}, {hi:.5f}]")
    print(f"Outputs in {a.out}/")

if __name__ == "__main__":
    main()