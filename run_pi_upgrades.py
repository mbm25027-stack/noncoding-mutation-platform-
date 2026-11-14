"""
run_pi_upgrades.py (Zero-pandas version for Windows Application Control)
Executes PI-requested analytical upgrades (Items 1, 4, 5) and regenerates
Figure 13B (forest plot with 95% CIs) and Figure 22B (legend clearance).
"""

import csv
import json
import math
import os
import re
import ssl
import time
import urllib.request
import numpy as np
import matplotlib.pyplot as plt
from nc_mutscan_research_suite import analyze_variant_delta

os.makedirs("results/pi_upgrades", exist_ok=True)
SSL_CTX = ssl._create_unverified_context()
FLANK_BP = 25
CACHE_FILE = "results/pi_upgrades/blind_flanks_cache.json"


def wilson_ci(k: int, n: int, z: float = 1.959963984540054):
    if n == 0:
        return (0.0, 0.0, 0.0)
    p = k / n
    denom = 1.0 + (z ** 2) / n
    center = (p + (z ** 2) / (2.0 * n)) / denom
    half = (z * math.sqrt((p * (1.0 - p) + (z ** 2) / (4.0 * n)) / n)) / denom
    return (p * 100.0, max(0.0, center - half) * 100.0, min(1.0, center + half) * 100.0)


def load_blind_records(vcf_path="locked_blind_chr12_22.vcf"):
    records = []
    with open(vcf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] in ("PATHOGENIC", "BENIGN"):
                label = parts[0]
                chrom = parts[1].replace("chr", "").strip()
                pos = int(parts[2])
                var_id = parts[3].strip()
                ref = parts[4].strip()
                alt = parts[5].split(",")[0].strip()
                extra = "\t".join(parts[6:]) if len(parts) > 6 else ""
            else:
                chrom = parts[0].replace("chr", "").strip()
                pos = int(parts[1])
                var_id = parts[2].strip()
                ref = parts[3].strip()
                alt = parts[4].split(",")[0].strip()
                extra = "\t".join(parts[5:]) if len(parts) > 5 else ""
                label = "PATHOGENIC" if "Pathogenic" in extra or "Likely_pathogenic" in extra else "BENIGN"
            records.append({
                "label": label,
                "chrom": chrom,
                "pos": pos,
                "var_id": var_id,
                "ref": ref,
                "alt": alt,
                "extra": extra,
            })
    return records


def get_flanks_cached(records):
    cache = {}
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, "rt", encoding="utf-8") as f:
            cache = json.load(f)

    updated = False
    for rec in records:
        key = f"chr{rec['chrom']}:{rec['pos']}:{rec['ref']}>{rec['alt']}"
        if key not in cache:
            start_0 = max(0, rec["pos"] - 1 - FLANK_BP)
            end_0 = rec["pos"] - 1 + len(rec["ref"]) + FLANK_BP
            url = (
                f"https://api.genome.ucsc.edu/getData/sequence?"
                f"genome=hg38;chrom=chr{rec['chrom']};start={start_0};end={end_0}"
            )
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=15, context=SSL_CTX) as resp:
                dna = json.loads(resp.read().decode("utf-8")).get("dna", "").upper()
            offset = (rec["pos"] - 1) - start_0
            mut_dna = dna[:offset] + rec["alt"].upper() + dna[offset + len(rec["ref"]):]
            cache[key] = {"wt": dna, "mut": mut_dna}
            updated = True
            time.sleep(0.05)

    if updated:
        with open(CACHE_FILE, "wt", encoding="utf-8") as f:
            json.dump(cache, f)
    return cache


def fetch_clinvar_mc(var_id: str, chrom: str, pos: int) -> str:
    uid = var_id if var_id.isdigit() else ""
    if not uid:
        esearch = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
            f"db=clinvar&term={chrom}[CHR]+AND+{pos}[CHRPOS38]&retmode=json"
        )
        req = urllib.request.Request(esearch, headers={"User-Agent": "NC-DriverEngine/3.3"})
        with urllib.request.urlopen(req, timeout=15, context=SSL_CTX) as resp:
            ids = json.loads(resp.read().decode("utf-8")).get("esearchresult", {}).get("idlist", [])
        if ids:
            uid = ids[0]
        time.sleep(0.35)

    if not uid:
        return ""

    esum = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={uid}&retmode=json"
    req = urllib.request.Request(esum, headers={"User-Agent": "NC-DriverEngine/3.3"})
    with urllib.request.urlopen(req, timeout=15, context=SSL_CTX) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    rec = data.get("result", {}).get(str(uid), {})
    mol_cons = rec.get("molecular_consequence_list", [])
    time.sleep(0.35)
    return "|".join(mol_cons) if mol_cons else ""


def classify_clinvar_mc(mc_raw: str) -> str:
    mc = (mc_raw or "").lower()
    if any(k in mc for k in ("splice_donor", "splice_acceptor", "splice_region", "splice")):
        return "Splice-region / canonical splice (below PWM cutoff)"
    if "intron" in mc:
        return "Deep-intronic (>20 bp from exon junction)"
    if any(k in mc for k in ("5_prime_utr", "upstream", "promoter", "tf_binding", "regulatory")):
        return "Promoter / 5'-UTR / TF motif outside 9 core PWMs"
    if "3_prime_utr" in mc:
        return "3'-UTR (non-polyadenylation signal)"
    if any(k in mc for k in ("non_coding", "ncrna", "lncrna", "mirna")):
        return "Non-coding RNA (ncRNA / miRNA)"
    return "Other / intergenic regulatory"


def run_item1_fn_breakdown(records, flank_cache):
    print("Running Item 1: Annotating blind false negatives by ClinVar mechanism...")
    fn_rows = []
    counts = {}

    for rec in records:
        if rec["label"] != "PATHOGENIC":
            continue
        key = f"chr{rec['chrom']}:{rec['pos']}:{rec['ref']}>{rec['alt']}"
        seqs = flank_cache[key]
        res = analyze_variant_delta(seqs["wt"], seqs["mut"], flank_bp=FLANK_BP, locus=f"chr{rec['chrom']}:{rec['pos']}")
        is_driver = res["classification"] in ("Driver Variant", "Loss-of-Function Variant")

        if not is_driver:
            mc_match = re.search(r"MC=([^;\t]+)", rec["extra"])
            mc = mc_match.group(1) if mc_match else ""
            if not mc:
                try:
                    mc = fetch_clinvar_mc(rec["var_id"], rec["chrom"], rec["pos"])
                except Exception:
                    mc = ""
            cat = classify_clinvar_mc(mc)
            counts[cat] = counts.get(cat, 0) + 1
            fn_rows.append({
                "chrom": f"chr{rec['chrom']}",
                "pos": rec["pos"],
                "var_id": rec["var_id"],
                "ref": rec["ref"],
                "alt": rec["alt"],
                "clinvar_mc": mc or "UNANNOTATED",
                "mechanism_category": cat,
                "v33_delta_pwm": res.get("delta_pwm", 0.0),
            })

    with open("results/pi_upgrades/blind_38_FN_detailed.csv", "wt", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(fn_rows[0].keys()))
        writer.writeheader()
        writer.writerows(fn_rows)

    total_fn = len(fn_rows)
    summary_rows = []
    for cat, cnt in sorted(counts.items(), key=lambda x: -x[1]):
        summary_rows.append({
            "Mechanism Category": cat,
            "Missed Pathogenic Variants (n)": cnt,
            "Share of FNs (%)": f"{cnt / total_fn * 100.0:.1f}%",
        })

    with open("results/pi_upgrades/Table16_blind_FN_mechanisms.csv", "wt", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Mechanism Category", "Missed Pathogenic Variants (n)", "Share of FNs (%)"])
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"\n=== ITEM 1: Table 16 — Mechanism Breakdown of {total_fn} Blind False Negatives ===")
    for r in summary_rows:
        print(f"  {r['Mechanism Category']:54s} | n = {r['Missed Pathogenic Variants (n)']:2d} ({r['Share of FNs (%)']})")


def query_cadd_v16(chrom: str, pos: int, ref: str, alt: str):
    url = f"https://cadd.gs.washington.edu/api/v1.0/GRCh38-v1.6/{chrom}:{pos}_{ref}_{alt}"
    req = urllib.request.Request(url, headers={"User-Agent": "NC-DriverEngine/3.3"})
    with urllib.request.urlopen(req, timeout=20, context=SSL_CTX) as resp:
        rows = json.loads(resp.read().decode("utf-8"))
    if not rows:
        return None
    return float(rows[0].get("PHRED"))


def run_item5_cadd_head_to_head(records, phred_cutoff=15.0):
    print("\nRunning Item 5: Head-to-head CADD v1.6 query on locked blind cohort...")
    scored_rows = []
    tp = fn = fp = tn = missing = 0

    for idx, rec in enumerate(records, 1):
        phred = None
        try:
            phred = query_cadd_v16(rec["chrom"], rec["pos"], rec["ref"], rec["alt"])
            time.sleep(0.12)
        except Exception:
            phred = None

        scored_rows.append({
            "chrom": f"chr{rec['chrom']}",
            "pos": rec["pos"],
            "ref": rec["ref"],
            "alt": rec["alt"],
            "label": rec["label"],
            "cadd_v16_phred": "" if phred is None else phred,
        })

        if phred is None:
            missing += 1
            continue

        pred_pos = phred >= phred_cutoff
        if rec["label"] == "PATHOGENIC":
            if pred_pos:
                tp += 1
            else:
                fn += 1
        else:
            if pred_pos:
                fp += 1
            else:
                tn += 1

        if idx % 50 == 0:
            print(f"  Queried CADD v1.6 for {idx}/{len(records)} variants...")

    with open("results/pi_upgrades/blind_150_cadd_v16_scores.csv", "wt", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["chrom", "pos", "ref", "alt", "label", "cadd_v16_phred"])
        writer.writeheader()
        writer.writerows(scored_rows)

    n_scored = tp + fn + fp + tn
    sens, s_lo, s_hi = wilson_ci(tp, tp + fn)
    spec, sp_lo, sp_hi = wilson_ci(tn, tn + fp)
    ppv, p_lo, p_hi = wilson_ci(tp, tp + fp)
    f1 = (2 * tp) / (2 * tp + fp + fn) if (2 * tp + fp + fn) else 0.0

    print("\n=== ITEM 5: Head-to-Head on Locked Blind ClinVar Cohort (chr12-22) ===")
    print(f"  NC-MutScan v3.3 (N=150)          : TP=37, FP=1, TN=74, FN=38 | Sens=49.33% | Spec=98.67% | PPV=97.37% | F1=0.6549")
    print(f"  CADD v1.6 >= {phred_cutoff:.0f} (N={n_scored}, miss={missing}): TP={tp}, FP={fp}, TN={tn}, FN={fn} | "
          f"Sens={sens:.2f}% ({s_lo:.1f}-{s_hi:.1f}) | Spec={spec:.2f}% ({sp_lo:.1f}-{sp_hi:.1f}) | "
          f"PPV={ppv:.2f}% ({p_lo:.1f}-{p_hi:.1f}) | F1={f1:.4f}")


def redraw_figures():
    # Figure 13B: True Forest Plot with 95% CIs
    genes = ["BCL6", "BCL2", "PAX5", "WDR74\n(neg. control)"]
    betas = [0.903, 1.322, 3.197, -0.158]
    ci_low = [0.653, 0.720, 1.011, -0.690]
    ci_high = [1.153, 1.924, 5.383, 0.374]
    q_labels = [
        "β = +0.903 [0.653, 1.153]; q = 1.68 × 10⁻¹¹",
        "β = +1.322 [0.720, 1.924]; q = 4.95 × 10⁻⁵",
        "β = +3.197 [1.011, 5.383]; q = 5.72 × 10⁻³",
        "β = −0.158 [−0.690, 0.374]; q = 0.559",
    ]
    colors = ["#b84a27", "#b84a27", "#b84a27", "#6c757d"]

    fig, ax = plt.subplots(figsize=(8.2, 3.4), dpi=300)
    y_pos = np.arange(len(genes))[::-1]
    ax.axvline(0, color="#222222", linewidth=1.0, zorder=1)
    for y, b, lo, hi, lbl, col in zip(y_pos, betas, ci_low, ci_high, q_labels, colors):
        ax.errorbar(b, y, xerr=[[b - lo], [hi - b]], fmt="o", color=col, ecolor=col, elinewidth=2.2, capsize=4, markersize=7, zorder=3)
        ax.text(hi + 0.14, y, lbl, va="center", fontsize=8.8, color="#111111")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(genes, fontsize=9.5, fontstyle="italic")
    ax.set_xlim(-1.0, 8.8)
    ax.set_xlabel(r"Covariate-adjusted $\beta_{\mathrm{mut}}$ ($\log_2$ FPKM-UQ + 1, 95% CI)", fontsize=9.5)
    ax.set_title("B  Covariate-adjusted eQTL effects with 95% CIs (matched PCAWG WGS + RNA-seq)", loc="left", fontweight="bold", fontsize=10.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig("results/pi_upgrades/Figure_13B_Forest_CI.png", dpi=300)
    plt.close(fig)

    # Figure 22B: Clean Legend Headroom
    metrics = ["Sensitivity\n(TPR)", "Specificity\n(TNR)", "Precision\n(PPV)", "F1-Score\n(×100)"]
    v31_dev = [72.0, 21.3, 47.8, 57.5]
    v33_dev = [66.7, 89.3, 86.2, 75.2]
    v33_blind = [49.3, 98.7, 97.4, 65.5]
    x = np.arange(len(metrics))
    w = 0.26

    fig, ax = plt.subplots(figsize=(6.2, 4.3), dpi=300)
    b1 = ax.bar(x - w, v31_dev, width=w, label="v3.1 Single-Seq (Dev)", color="#9aa0a6")
    b2 = ax.bar(x, v33_dev, width=w, label="v3.3 ΔPWM (Dev chr1–11)", color="#2b6cb0")
    b3 = ax.bar(x + w, v33_blind, width=w, label="v3.3 ΔPWM (Blind chr12–22)", color="#1b7f46")
    for bars in (b1, b2, b3):
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2.0, h + 1.5, f"{h:.1f}", ha="center", va="bottom", fontsize=7.8, fontweight="bold")
    ax.set_ylim(0, 126)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics, fontsize=9)
    ax.set_ylabel("Percentage / Score (%)", fontsize=9.5)
    ax.set_title("B. Differential ΔPWM Specificity Gain", fontweight="bold", fontsize=10.5)
    ax.legend(loc="upper right", frameon=True, fontsize=8.2)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig("results/pi_upgrades/Figure_22B_Fixed_Legend.png", dpi=300)
    plt.close(fig)
    print("Saved Figure_13B_Forest_CI.png and Figure_22B_Fixed_Legend.png in results/pi_upgrades/")


if __name__ == "__main__":
    redraw_figures()
    recs = load_blind_records("locked_blind_chr12_22.vcf")
    flanks = get_flanks_cached(recs)
    run_item1_fn_breakdown(recs, flanks)
    run_item5_cadd_head_to_head(recs, phred_cutoff=15.0)
