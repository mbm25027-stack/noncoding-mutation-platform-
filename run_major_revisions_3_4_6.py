"""
run_major_revisions_3_4_6.py
Executes Major Revision Items 3, 4, and 6 on local project files:
  - Item 3: Mechanistic breakdown of all 38 False Negatives on locked_blind_chr12_22.vcf
            + side-by-side Full N=150 vs SpliceAI-Matched N=139 vs Unmatched N=11 metrics.
  - Item 4: 23-fold Leave-One-Chromosome-Out (LOCO) cross-validation across all N=1,105
            variants (dev_chr1_11.vcf + locked_blind_chr12_22.vcf), both fixed v3.3 and
            nested fold-wise TF threshold tuning on training chromosomes only.
  - Item 6: SHA-256 checksums and packaged release bundle in results/release_bundle/.
"""

import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil
from collections import Counter


def norm_key(s):
    s = str(s).strip().lower()
    if s.startswith("chr"):
        s = s[3:]
    return re.sub(r"[^a-z0-9]+", ":", s).strip(":")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def load_vcf(vcf_path, split_name):
    rows = []
    with open(vcf_path, "r", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            ln = ln.rstrip("\r\n")
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split("\t")
            if parts[0].upper() in ("PATHOGENIC", "BENIGN", "LIKELY_PATHOGENIC", "LIKELY_BENIGN"):
                label = 1 if "PATHOGENIC" in parts[0].upper() else 0
                chrom, pos, vid, ref, alt = parts[1], int(parts[2]), parts[3], parts[4], parts[5]
                info = parts[7] if len(parts) > 7 else ""
            else:
                chrom, pos, vid, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
                info = parts[7] if len(parts) > 7 else ""
                label = 1 if "pathogenic" in info.lower() else 0
            c_clean = chrom if str(chrom).lower().startswith("chr") else f"chr{chrom}"
            rows.append({
                "split": split_name,
                "chrom": c_clean,
                "pos": pos,
                "vid": str(vid).strip(),
                "ref": ref,
                "alt": alt,
                "info": info,
                "nk_full": norm_key(f"{chrom}:{pos}:{ref}:{alt}"),
                "nk_cp": norm_key(f"{chrom}:{pos}"),
                "label": label,
            })
    return rows


def is_frozen_driver_call(res):
    cls_txt = str(res.get("classification", "")).lower()
    if any(neg in cls_txt for neg in ("benign", "neutral", "polymorphism", "passenger", "sub-threshold", "subthreshold", "non-driver", "compet")):
        return 0
    if any(pos in cls_txt for pos in ("driver", "pathogenic", "loss", "gain", "splice", "polya", "disrupt")):
        return 1
    return 0


def wilson_ci(k, n, z=1.959963984540054):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1.0 + z * z / n
    center = p + z * z / (2.0 * n)
    rad = z * math.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n)
    return max(0.0, (center - rad) / denom) * 100.0, min(1.0, (center + rad) / denom) * 100.0


def metrics_dict(y_true, y_pred):
    tp = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 1)
    fp = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 1)
    tn = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 0)
    fn = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 0)
    sens = 100.0 * tp / (tp + fn) if (tp + fn) else 0.0
    spec = 100.0 * tn / (tn + fp) if (tn + fp) else 0.0
    ppv = 100.0 * tp / (tp + fp) if (tp + fp) else 0.0
    f1 = (2.0 * tp) / (2.0 * tp + fp + fn) if (2.0 * tp + fp + fn) else 0.0
    s_lo, s_hi = wilson_ci(tp, tp + fn)
    sp_lo, sp_hi = wilson_ci(tn, tn + fp)
    return {
        "N": len(y_true), "TP": tp, "FP": fp, "TN": tn, "FN": fn,
        "Sens": f"{sens:.2f}% ({s_lo:.1f}-{s_hi:.1f}%)",
        "Spec": f"{spec:.2f}% ({sp_lo:.1f}-{sp_hi:.1f}%)",
        "PPV": f"{ppv:.2f}%",
        "F1": f"{f1:.3f}",
    }


def classify_fn_mechanism(v, res, dist_bp, sai_ds):
    info_l = v["info"].lower()
    wt_seq = str(res.get("wt_seq", "") or "").upper()
    raw_dpwm = abs(float(res.get("delta_pwm", 0.0) or 0.0))
    cls_txt = str(res.get("classification", "")).lower()

    if "aid" in cls_txt or "passenger" in cls_txt:
        return "AID/SHM passenger-masked regulatory locus"
    if dist_bp <= 2:
        return "Canonical ±1–2 bp splice site (WT score or |ΔPWM| < τ_splice)"
    if 3 <= dist_bp <= 17:
        return "Extended splice region (±3–17 bp; sub-threshold |ΔPWM_splice|)"
    if 18 <= dist_bp <= 45 and ("intron" in info_l or sai_ds is not None):
        if re.search(r"[CT]T[AG][AC][CT]", wt_seq):
            return "Putative intronic branch-point window (−18 to −45 bp; YUNAY motif)"
        return "Proximal intronic splicing regulatory region (18–45 bp)"
    if "utr" in info_l or "3_prime" in info_l or "5_prime" in info_l:
        return "5′/3′-UTR regulatory variant (non-AATAAA/ATTAAA or unmodeled RBP/miRNA)"
    if "promoter" in info_l or "upstream" in info_l or "enhancer" in info_l or "tf" in info_l:
        return f"Promoter/enhancer TFBS outside 9-TF library or sub-4.2 bits (max |ΔPWM|={raw_dpwm:.2f})"
    if dist_bp > 45:
        if sai_ds is not None and sai_ds >= 0.20:
            return "Deep-intronic pseudoexon / cryptic splice variant (>45 bp; SpliceAI DS≥0.20)"
        return "Deep-intronic or distal cis-regulatory variant (>45 bp; outside local PWM scope)"
    return "Other noncoding regulatory mechanism"


def main():
    os.makedirs(os.path.join("results", "major_revisions"), exist_ok=True)
    os.makedirs(os.path.join("results", "release_bundle"), exist_ok=True)

    # ITEM 6: SHA-256 Checksums & Release Bundle
    release_files = [
        "nc_mutscan_research_suite.py",
        "dev_chr1_11.vcf",
        "locked_blind_chr12_22.vcf",
    ]
    print("=== ITEM 6: CODE & VCF SHA-256 RELEASE MANIFEST ===")
    manifest_lines = []
    for rf in release_files:
        if os.path.isfile(rf):
            digest = sha256_file(rf)
            shutil.copy2(rf, os.path.join("results", "release_bundle", os.path.basename(rf)))
            line = f"{digest}  {os.path.basename(rf)}"
            manifest_lines.append(line)
            print(" ", line)
    with open(os.path.join("results", "release_bundle", "SHA256SUMS.txt"), "w", encoding="utf-8") as mf:
        mf.write("\n".join(manifest_lines) + "\n")

    spec = importlib.util.spec_from_file_location("frozen_mutscan", "nc_mutscan_research_suite.py")
    ms_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ms_mod)

    flanks = {}
    for cache_candidate in [
        os.path.join("results", "pi_upgrades", "blind_flanks_cache.json"),
        os.path.join("results", "pi_upgrades", "dev_flanks_cache.json"),
        os.path.join("results", "clean_submission", "all_flanks_cache.json"),
    ]:
        if os.path.isfile(cache_candidate):
            with open(cache_candidate, "r", encoding="utf-8") as ff:
                raw_c = json.load(ff)
            for k, val in raw_c.items():
                flanks[str(k).strip()] = val
                flanks[norm_key(k)] = val

    with open(os.path.join("results", "pi_upgrades", "blind_150_exon_distances.json"), "r", encoding="utf-8") as ef:
        exon_dists = json.load(ef)

    with open(os.path.join("results", "clean_submission", "blind_150_spliceai.json"), "r", encoding="utf-8") as jf:
        sdata = json.load(jf)
    spliceai_scores = {}
    for k, v in sdata.items():
        if isinstance(v, dict):
            vals = [float(x) for x in v.values() if isinstance(x, (int, float)) and not isinstance(x, bool)]
            if vals:
                spliceai_scores[norm_key(k)] = max(vals)
                spliceai_scores[str(k).strip()] = max(vals)
        elif isinstance(v, (int, float)):
            spliceai_scores[norm_key(k)] = float(v)

    # ITEM 3: 38 False Negatives + Full N=150 vs Matched N=139
    blind_vars = load_vcf("locked_blind_chr12_22.vcf", "blind_chr12_22")
    y_full, p_full = [], []
    y_match, p_match = [], []
    y_unmatch, p_unmatch = [], []
    fn_records = []

    for v in blind_vars:
        entry = flanks.get(v["nk_full"], flanks.get(v["vid"], flanks.get(v["nk_cp"])))
        res = ms_mod.analyze_variant_delta(wt_seq=entry["wt"], mut_seq=entry["mut"], flank_bp=25)
        res["wt_seq"] = entry["wt"]
        pred = is_frozen_driver_call(res)
        ds = spliceai_scores.get(v["nk_full"], spliceai_scores.get(v["vid"], spliceai_scores.get(v["nk_cp"])))
        dist_bp = int(exon_dists.get(v["nk_cp"], 500))

        y_full.append(v["label"])
        p_full.append(pred)
        if ds is not None:
            y_match.append(v["label"])
            p_match.append(pred)
        else:
            y_unmatch.append(v["label"])
            p_unmatch.append(pred)

        if v["label"] == 1 and pred == 0:
            mech = classify_fn_mechanism(v, res, dist_bp, ds)
            fn_records.append({
                "chrom": v["chrom"],
                "pos": v["pos"],
                "vid": v["vid"],
                "ref": v["ref"],
                "alt": v["alt"],
                "exon_dist_bp": dist_bp,
                "spliceai_ds_max": f"{ds:.3f}" if ds is not None else "NA (unmatched)",
                "max_raw_dpwm": f"{abs(float(res.get('delta_pwm', 0.0) or 0.0)):.2f}",
                "mechanism": mech,
                "info": v["info"],
            })

    print("\n=== ITEM 3A: LOCKED BLIND COHORT FULL (N=150) VS MATCHED (N=139) VS UNMATCHED (N=11) ===")
    for cohort_lbl, yt, yp in [
        ("Full Locked Blind Cohort (N=150)", y_full, p_full),
        ("SpliceAI-Matched Subset (N=139)", y_match, p_match),
        ("Unmatched Non-Transcript/Distal Subset (N=11)", y_unmatch, p_unmatch),
    ]:
        m = metrics_dict(yt, yp)
        print(f"  {cohort_lbl}: TP={m['TP']}, FP={m['FP']}, TN={m['TN']}, FN={m['FN']} | Sens={m['Sens']} | Spec={m['Spec']} | PPV={m['PPV']} | F1={m['F1']}")

    print(f"\n=== ITEM 3B: MECHANISTIC BREAKDOWN OF ALL {len(fn_records)} FALSE NEGATIVES ===")
    mech_counts = Counter(r["mechanism"] for r in fn_records)
    for mech, cnt in mech_counts.most_common():
        pct = 100.0 * cnt / len(fn_records)
        sai_hits = sum(1 for r in fn_records if r["mechanism"] == mech and r["spliceai_ds_max"] != "NA (unmatched)" and float(r["spliceai_ds_max"]) >= 0.20)
        print(f"  - {mech}: {cnt}/{len(fn_records)} ({pct:.1f}%) [SpliceAI DS>=0.20 in {sai_hits}/{cnt}]")

    fn_csv = os.path.join("results", "major_revisions", "blind_38_false_negatives_mechanisms.csv")
    with open(fn_csv, "w", newline="", encoding="utf-8") as cf:
        w = csv.DictWriter(cf, fieldnames=list(fn_records[0].keys()))
        w.writeheader()
        w.writerows(fn_records)
    print(f"Saved detailed 38-FN table -> {fn_csv}")

    # ITEM 4: 23-Fold Leave-One-Chromosome-Out (LOCO) Cross-Validation
    dev_vars = load_vcf("dev_chr1_11.vcf", "dev_chr1_11_X")
    all_vars = dev_vars + blind_vars
    evaluated_all = []
    missing_flanks = 0
    for v in all_vars:
        entry = flanks.get(v["nk_full"], flanks.get(v["vid"], flanks.get(v["nk_cp"])))
        if entry is None:
            missing_flanks += 1
            continue
        res = ms_mod.analyze_variant_delta(wt_seq=entry["wt"], mut_seq=entry["mut"], flank_bp=25)
        pred_v33 = is_frozen_driver_call(res)
        cls_txt = str(res.get("classification", "")).lower()
        raw_dpwm = abs(float(res.get("delta_pwm", 0.0) or 0.0))
        is_non_tf_hit = 1 if (pred_v33 == 1 and any(k in cls_txt for k in ("splice", "polya"))) else 0
        is_masked = 1 if any(k in cls_txt for k in ("aid", "passenger", "compet")) else 0
        evaluated_all.append({
            "chrom": v["chrom"],
            "label": v["label"],
            "pred_v33": pred_v33,
            "is_non_tf_hit": is_non_tf_hit,
            "is_masked": is_masked,
            "raw_dpwm": raw_dpwm,
        })

    chroms = sorted(set(x["chrom"] for x in evaluated_all), key=lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, c))
    print(f"\n=== ITEM 4: LEAVE-ONE-CHROMOSOME-OUT (LOCO) CROSS-VALIDATION ({len(chroms)} CHROMOSOMES, N={len(evaluated_all)}) ===")
    if missing_flanks:
        print(f"  (Note: {missing_flanks} dev variants lacked cached flanks; build dev_flanks_cache.json to evaluate all 1,105)")

    oof_fixed = metrics_dict([x["label"] for x in evaluated_all], [x["pred_v33"] for x in evaluated_all])
    print(f"  Pooled LOCO (Fixed v3.3 rules): TP={oof_fixed['TP']}, FP={oof_fixed['FP']}, TN={oof_fixed['TN']}, FN={oof_fixed['FN']} | Sens={oof_fixed['Sens']} | Spec={oof_fixed['Spec']} | PPV={oof_fixed['PPV']}")

    tau_grid = [round(3.5 + 0.1 * i, 1) for i in range(16)]
    oof_tuned_preds = []
    oof_tuned_labels = []
    chosen_taus = []
    for held_chr in chroms:
        train = [x for x in evaluated_all if x["chrom"] != held_chr]
        test = [x for x in evaluated_all if x["chrom"] == held_chr]
        best_tau = 4.2
        best_obj = -1.0
        for tau in tau_grid:
            yt_tr = [x["label"] for x in train]
            yp_tr = [1 if (x["is_non_tf_hit"] or (not x["is_masked"] and x["raw_dpwm"] >= tau)) else 0 for x in train]
            tp = sum(1 for a, b in zip(yt_tr, yp_tr) if a == 1 and b == 1)
            fp = sum(1 for a, b in zip(yt_tr, yp_tr) if a == 0 and b == 1)
            fn = sum(1 for a, b in zip(yt_tr, yp_tr) if a == 1 and b == 0)
            f05 = (1.25 * tp) / (1.25 * tp + 0.25 * fn + fp) if (tp + fp + fn) else 0.0
            if f05 > best_obj:
                best_obj = f05
                best_tau = tau
        chosen_taus.append(best_tau)
        for x in test:
            pred_te = 1 if (x["is_non_tf_hit"] or (not x["is_masked"] and x["raw_dpwm"] >= best_tau)) else 0
            oof_tuned_labels.append(x["label"])
            oof_tuned_preds.append(pred_te)

    oof_nested = metrics_dict(oof_tuned_labels, oof_tuned_preds)
    print(f"  Nested Fold-Tuned LOCO (tau_TF tuned on 22 train chroms; mode={Counter(chosen_taus).most_common(1)[0][0]} bits): TP={oof_nested['TP']}, FP={oof_nested['FP']}, TN={oof_nested['TN']}, FN={oof_nested['FN']} | Sens={oof_nested['Sens']} | Spec={oof_nested['Spec']} | PPV={oof_nested['PPV']}")


if __name__ == "__main__":
    main()
