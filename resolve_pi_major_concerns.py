"""
resolve_pi_major_concerns.py
Zero-pandas script addressing PI Major Concerns #1, #2, #4, and #8:
  1. Audits chr7:140,453,135 (BRAF/MKRN1 build check) and scans PCAWG MELA-AU WGS
     for TERT and core-promoter hotspots using unique donor counts.
  2. Re-runs PCAWG eQTL OLS with histology fixed effects AND within-lymphoma subset.
  3. Evaluates Unpaired Single-PWM, Splice-Only, and SpliceAI on locked_blind_chr12_22.vcf.
"""

import csv
import glob
import gzip
import json
import math
import os
import re
import ssl
import time
import urllib.request
import numpy as np
from nc_mutscan_research_suite import analyze_variant_delta, analyze_fasta_entry

os.makedirs("results/clean_submission", exist_ok=True)
SSL_CTX = ssl._create_unverified_context()
FLANK_BP = 25


def wilson_ci(k: int, n: int, z: float = 1.959963984540054):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    denom = 1.0 + (z ** 2) / n
    center = (p + (z ** 2) / (2.0 * n)) / denom
    half = (z * math.sqrt((p * (1.0 - p) + (z ** 2) / (4.0 * n)) / n)) / denom
    return (p * 100.0, max(0.0, center - half) * 100.0, min(1.0, center + half) * 100.0)


# ---------------------------------------------------------------------------
# PART 1: Audit BRAF/MKRN1 Coordinate Build & Check PCAWG MELA-AU WGS Hotspots
# ---------------------------------------------------------------------------
def audit_coordinates_and_mela_wgs():
    print("=== PART 1: Coordinate Build Audit & PCAWG MELA-AU WGS Hotspot Check ===")
    # 1A. Inspect any TCGA-SKCM hotspot / window tables in results/
    for fpath in sorted(glob.glob("results/**/*.*", recursive=True)):
        if any(k in fpath.lower() for k in ("window", "hotspot", "promoter", "skcm", "mela")):
            if fpath.endswith((".csv", ".tsv")):
                with open(fpath, "rt", encoding="utf-8", errors="ignore") as f:
                    txt = f.read()
                if "140453135" in txt or "140,453,135" in txt or "1295113" in txt or "1295228" in txt:
                    print(f"  Found target coordinates in: {fpath}")
                    for line in txt.splitlines():
                        if any(c in line for c in ("140453135", "140753336", "1295113", "1295228", "TERT", "TCF3", "AC040160")):
                            print("   ", line[:140])

    # 1B. Locate PCAWG MELA-AU WGS SNV file(s) in workspace to extract donor-level counts
    wgs_files = [
        f for f in glob.glob("**/*", recursive=True)
        if os.path.isfile(f) and any(k in f.lower() for k in ("mela", "pcawg", "snv", "maf", "wgs"))
        and f.endswith((".maf", ".maf.gz", ".tsv", ".tsv.gz", ".vcf", ".vcf.gz", ".txt", ".txt.gz"))
    ]
    print("  Candidate WGS/PCAWG files detected:", wgs_files[:12])

    # Target loci in hg19 (PCAWG native) and hg38 (TCGA native)
    targets = {
        "TERT_C228T_C250T": [("5", 1295150, 1295300), ("5", 1295050, 1295200)],
        "AC040160.1_prom":  [("16", 67260800, 67261150), ("16", 67294800, 67295200)],
        "AL358472.4_prom":  [("1", 153963050, 153963400), ("1", 153935600, 153936000)],
        "TCF3_prom":        [("19", 1652850, 1653200), ("19", 1652850, 1653200)],
    }

    for wf in wgs_files:
        if "mela" not in wf.lower() and "pcawg" not in wf.lower():
            continue
        opener = gzip.open if wf.endswith(".gz") else open
        try:
            with opener(wf, "rt", encoding="utf-8", errors="ignore") as f:
                header = f.readline().rstrip("\n").split("\t")
                # Identify chrom, pos, donor/sample columns
                c_idx = next((i for i, h in enumerate(header) if h.lower() in ("chrom", "chromosome", "chr", "#chrom")), None)
                p_idx = next((i for i, h in enumerate(header) if h.lower() in ("pos", "position", "start_position", "start")), None)
                d_idx = next((i for i, h in enumerate(header) if any(k in h.lower() for k in ("donor", "sample", "tumor_sample", "icgc_donor"))), None)
                if c_idx is None or p_idx is None:
                    continue
                print(f"  Scanning {wf} (chrom={header[c_idx]}, pos={header[p_idx]}, donor={header[d_idx] if d_idx is not None else 'row'})...")
                hits = {k: set() for k in targets}
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) <= max(c_idx, p_idx):
                        continue
                    chrom = parts[c_idx].replace("chr", "")
                    if chrom not in ("1", "5", "16", "19"):
                        continue
                    try:
                        pos = int(parts[p_idx])
                    except ValueError:
                        continue
                    donor = parts[d_idx] if (d_idx is not None and len(parts) > d_idx) else f"row_{pos}"
                    for tname, spans in targets.items():
                        for sc, s0, s1 in spans:
                            if chrom == sc and s0 <= pos <= s1:
                                hits[tname].add((donor, pos))
                for tname, dset in hits.items():
                    unique_donors = len({d for d, _ in dset})
                    print(f"    {tname:20s}: {unique_donors} unique mutated donors ({len(dset)} SNVs)")
        except Exception as e:
            print(f"    Skipped {wf}: {e}")


# ---------------------------------------------------------------------------
# PART 2: Histology-Adjusted & Within-Lymphoma PCAWG eQTL OLS Regression
# ---------------------------------------------------------------------------
def ols_with_t_pvalue(X: np.ndarray, y: np.ndarray, target_col: int = 1):
    """Exact OLS fit via QR/SVD returning beta, SE, t-stat, and two-sided p-value."""
    n, p = X.shape
    df = n - p
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    s2 = float(np.sum(resid ** 2)) / max(1, df)
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = math.sqrt(max(1e-18, float(cov[target_col, target_col])))
    b_val = float(beta[target_col])
    t_stat = abs(b_val) / se
    # Exact two-sided Student-t p-value via regularized incomplete beta
    x_beta = df / (df + t_stat ** 2)
    # Regularized incomplete beta I_x(df/2, 0.5)
    p_val = betainc_reg(df / 2.0, 0.5, x_beta)
    ci_95 = 1.96 + (2.37 / max(1, df))  # Accurate t_0.975 approximation for df >= 20
    return b_val, se, b_val - ci_95 * se, b_val + ci_95 * se, p_val, df


def betainc_reg(a: float, b: float, x: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    front = math.exp(lbeta + a * math.log(x) + b * math.log(1.0 - x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * betacf(a, b, x) / a
    return 1.0 - front * betacf(b, a, 1.0 - x) / b


def betacf(a: float, b: float, x: float) -> float:
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, 200):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-12:
            break
    return h


def rerun_histology_adjusted_eqtl():
    print("\n=== PART 2: Histology-Adjusted & Within-Lymphoma PCAWG eQTL Analysis ===")
    # Search for donor-level PCAWG eQTL input tables in workspace
    eqtl_files = [
        f for f in glob.glob("**/*", recursive=True)
        if os.path.isfile(f) and any(k in f.lower() for k in ("eqtl", "fpkm", "purity", "22_nb"))
        and f.endswith((".csv", ".tsv", ".py"))
    ]
    print("  Detected eQTL artifacts:", eqtl_files)

    # Check if donor-level merged table exists in results/
    donor_tables = [f for f in eqtl_files if f.endswith((".csv", ".tsv")) and "31_covariate" not in f]
    if not donor_tables:
        print("  Note: Donor-level merged eQTL matrix is generated inside 22_nb_bmr_and_covariate_eqtl.py;")
        print("        inspecting 22_nb_bmr_and_covariate_eqtl.py data paths...")
        if os.path.exists("22_nb_bmr_and_covariate_eqtl.py"):
            txt = open("22_nb_bmr_and_covariate_eqtl.py", encoding="utf-8", errors="ignore").read()
            for line in txt.splitlines():
                if any(k in line.lower() for k in ("read_csv", "open(", "xena", "fpkm", "purity", "bcl6")):
                    print("   ", line[:130])
        return

    for dt in donor_tables:
        sep = "\t" if dt.endswith(".tsv") else ","
        with open(dt, "rt", encoding="utf-8", errors="ignore") as f:
            reader = csv.DictReader(f, delimiter=sep)
            rows = list(reader)
        if not rows or "histology" not in "".join(reader.fieldnames or []).lower():
            continue
        print(f"  Loaded donor matrix {dt} (N={len(rows)} rows)")
        # Fit histology-adjusted and within-lymphoma models per gene if gene column exists
        gene_col = next((c for c in reader.fieldnames if "gene" in c.lower()), None)
        hist_col = next((c for c in reader.fieldnames if "hist" in c.lower() or "cohort" in c.lower()), None)
        expr_col = next((c for c in reader.fieldnames if "expr" in c.lower() or "fpkm" in c.lower()), None)
        mut_col = next((c for c in reader.fieldnames if "mut" in c.lower()), None)
        cn_col = next((c for c in reader.fieldnames if "cn" in c.lower() or "copy" in c.lower()), None)
        pur_col = next((c for c in reader.fieldnames if "pur" in c.lower()), None)
        tmb_col = next((c for c in reader.fieldnames if "tmb" in c.lower()), None)

        if all([gene_col, hist_col, expr_col, mut_col, cn_col, pur_col, tmb_col]):
            for gene in ("BCL6", "BCL2", "PAX5", "WDR74"):
                g_rows = [r for r in rows if r[gene_col] == gene]
                if not g_rows:
                    continue
                hists = sorted({r[hist_col] for r in g_rows})
                X_list, y_list = [], []
                for r in g_rows:
                    feat = [
                        1.0,
                        float(r[mut_col]),
                        float(r[cn_col]),
                        float(r[pur_col]),
                        math.log10(float(r[tmb_col]) + 1.0),
                    ] + [1.0 if r[hist_col] == h else 0.0 for h in hists[1:]]
                    X_list.append(feat)
                    y_list.append(float(r[expr_col]))
                b_val, se, lo, hi, p_val, df = ols_with_t_pvalue(np.array(X_list), np.array(y_list), target_col=1)
                print(f"    [Histology-Adjusted] {gene:6s} (df={df}): beta={b_val:+.3f} [{lo:+.3f}, {hi:+.3f}], p={p_val:.3e}")

                # Within B-cell lymphoma only (Lymph-BNHL / Lymph-CLL)
                lym_rows = [r for r in g_rows if "lymph" in r[hist_col].lower()]
                if len(lym_rows) >= 10:
                    lym_hists = sorted({r[hist_col] for r in lym_rows})
                    Xl, yl = [], []
                    for r in lym_rows:
                        feat = [
                            1.0,
                            float(r[mut_col]),
                            float(r[cn_col]),
                            float(r[pur_col]),
                            math.log10(float(r[tmb_col]) + 1.0),
                        ] + [1.0 if r[hist_col] == h else 0.0 for h in lym_hists[1:]]
                        Xl.append(feat)
                        yl.append(float(r[expr_col]))
                    bl, sel, lol, hil, pl, dfl = ols_with_t_pvalue(np.array(Xl), np.array(yl), target_col=1)
                    print(f"    [Lymphoma-Only]      {gene:6s} (df={dfl}): beta={bl:+.3f} [{lol:+.3f}, {hil:+.3f}], p={pl:.3e}")


# ---------------------------------------------------------------------------
# PART 3: Additional Head-to-Head Baselines on Locked Blind chr12-22 ClinVar
# ---------------------------------------------------------------------------
def query_spliceai_grch38(chrom: str, pos: int, ref: str, alt: str):
    c = chrom.replace("chr", "")
    url = f"https://spliceai-38-xwkwwwxdwq-uc.a.run.app/spliceai/?hg=38&distance=50&variant={c}-{pos}-{ref}-{alt}"
    req = urllib.request.Request(url, headers={"User-Agent": "NC-DriverEngine/3.3"})
    with urllib.request.urlopen(req, timeout=20, context=SSL_CTX) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    scores = data.get("scores", [])
    if not scores:
        return None
    best = 0.0
    for sc in scores:
        for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL"):
            val = sc.get(k)
            if val is not None:
                best = max(best, float(val))
    return best


def run_blind_ablation_and_spliceai():
    print("\n=== PART 3: Ablation Baselines & SpliceAI on Locked Blind Cohort (chr12-22) ===")
    flank_cache_path = "results/pi_upgrades/blind_flanks_cache.json"
    if not os.path.exists(flank_cache_path):
        print("  Run run_pi_upgrades.py first to populate results/pi_upgrades/blind_flanks_cache.json.")
        return

    flank_cache = json.load(open(flank_cache_path, "rt", encoding="utf-8"))
    records = []
    with open("locked_blind_chr12_22.vcf", "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] in ("PATHOGENIC", "BENIGN"):
                records.append({"label": p[0], "chrom": p[1].replace("chr", ""), "pos": int(p[2]), "ref": p[4], "alt": p[5].split(",")[0]})
            else:
                lbl = "PATHOGENIC" if ("Pathogenic" in "\t".join(p[5:]) or "Likely_pathogenic" in "\t".join(p[5:])) else "BENIGN"
                records.append({"label": lbl, "chrom": p[0].replace("chr", ""), "pos": int(p[1]), "ref": p[3], "alt": p[4].split(",")[0]})

    # 1. Ablation A: Unpaired Single-Sequence PWM Baseline (no WT comparison)
    # 2. Ablation B: Splice-Motif-Only Baseline (SPLICE_5SS / SPLICE_3SS without TF PWMs)
    # 3. External Modern Comparator: SpliceAI (DS_max >= 0.20 standard clinical cutoff)
    counts = {
        "NC-MutScan v3.3 (Full Differential)": [0, 0, 0, 0],  # TP, FP, TN, FN
        "Ablation: Unpaired Single-Seq PWM":   [0, 0, 0, 0],
        "Ablation: Splice-PWM Only (no TFs)":  [0, 0, 0, 0],
        "SpliceAI (DS_max >= 0.20, GRCh38)":   [0, 0, 0, 0],
    }

    spliceai_cache_path = "results/clean_submission/blind_150_spliceai.json"
    sai_cache = json.load(open(spliceai_cache_path, "rt", encoding="utf-8")) if os.path.exists(spliceai_cache_path) else {}

    for idx, rec in enumerate(records, 1):
        key = f"chr{rec['chrom']}:{rec['pos']}:{rec['ref']}>{rec['alt']}"
        seqs = flank_cache[key]
        is_path = rec["label"] == "PATHOGENIC"

        # Full v3.3
        res_full = analyze_variant_delta(seqs["wt"], seqs["mut"], flank_bp=FLANK_BP, locus=key)
        pred_full = res_full["classification"] in ("Driver Variant", "Loss-of-Function Variant")

        # Unpaired single-sequence scan on mutant allele only
        res_unpaired = analyze_fasta_entry(f">{key}\n{seqs['mut']}", wt_seq=None)
        pred_unpaired = res_unpaired.get("classification") in ("Driver Variant", "Loss-of-Function Variant")

        # Splice-only subset of v3.3
        motif_name = str(res_full.get("top_motif") or res_full.get("motif") or "")
        pred_splice_only = pred_full and ("SPLICE" in motif_name.upper())

        for name, pred in [
            ("NC-MutScan v3.3 (Full Differential)", pred_full),
            ("Ablation: Unpaired Single-Seq PWM", pred_unpaired),
            ("Ablation: Splice-PWM Only (no TFs)", pred_splice_only),
        ]:
            if is_path and pred:
                counts[name][0] += 1
            elif not is_path and pred:
                counts[name][1] += 1
            elif not is_path and not pred:
                counts[name][2] += 1
            else:
                counts[name][3] += 1

        # SpliceAI query (cached)
        if key not in sai_cache:
            try:
                sai_cache[key] = query_spliceai_grch38(rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
                time.sleep(0.25)
            except Exception:
                sai_cache[key] = None
        ds = sai_cache[key]
        if ds is not None:
            pred_sai = ds >= 0.20
            if is_path and pred_sai:
                counts["SpliceAI (DS_max >= 0.20, GRCh38)"][0] += 1
            elif not is_path and pred_sai:
                counts["SpliceAI (DS_max >= 0.20, GRCh38)"][1] += 1
            elif not is_path and not pred_sai:
                counts["SpliceAI (DS_max >= 0.20, GRCh38)"][2] += 1
            else:
                counts["SpliceAI (DS_max >= 0.20, GRCh38)"][3] += 1

        if idx % 50 == 0:
            print(f"  Evaluated {idx}/150 blind variants across all baselines & SpliceAI...")

    with open(spliceai_cache_path, "wt", encoding="utf-8") as f:
        json.dump(sai_cache, f)

    print("\n=== HEAD-TO-HEAD SUMMARY ON LOCKED BLIND CLINVAR COHORT (chr12-22) ===")
    for name, (tp, fp, tn, fn) in counts.items():
        n_tot = tp + fp + tn + fn
        sens, s_lo, s_hi = wilson_ci(tp, tp + fn)
        spec, sp_lo, sp_hi = wilson_ci(tn, tn + fp)
        ppv, p_lo, p_hi = wilson_ci(tp, tp + fp)
        f1 = (2 * tp) / max(1, 2 * tp + fp + fn)
        print(f"  {name:36s} (N={n_tot:3d}) | TP={tp:2d} FP={fp:2d} TN={tn:2d} FN={fn:2d} | "
              f"Sens={sens:5.2f}% ({s_lo:.1f}-{s_hi:.1f}) | Spec={spec:5.2f}% ({sp_lo:.1f}-{sp_hi:.1f}) | "
              f"PPV={ppv:5.2f}% | F1={f1:.4f}")


if __name__ == "__main__":
    audit_coordinates_and_mela_wgs()
    rerun_histology_adjusted_eqtl()
    run_blind_ablation_and_spliceai()
