"""
build_dev_flanks_and_loco.py
Fetches hg38 ±25 bp flanking sequences for dev_chr1_11.vcf (cached to
results/major_revisions/dev_flanks_cache.json) and runs true 23-chromosome
Leave-One-Chromosome-Out (LOCO) cross-validation across all N=300 variants
(dev_chr1_11.vcf + locked_blind_chr12_22.vcf).
"""

import importlib.util
import json
import math
import os
import re
import time
import urllib.request
from collections import Counter


def norm_key(s):
    s = str(s).strip().lower()
    if s.startswith("chr"):
        s = s[3:]
    return re.sub(r"[^a-z0-9]+", ":", s).strip(":")


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
                "ref": ref.upper(),
                "alt": alt.upper(),
                "nk_full": norm_key(f"{chrom}:{pos}:{ref}:{alt}"),
                "nk_cp": norm_key(f"{chrom}:{pos}"),
                "label": label,
            })
    return rows


def fetch_hg38_flank(chrom, pos, ref, alt, flank_bp=25):
    start0 = pos - 1 - flank_bp
    end1 = pos - 1 + len(ref) + flank_bp
    ucsc_url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start0};end={end1}"
    try:
        req = urllib.request.Request(ucsc_url, headers={"User-Agent": "NC-MutScan/3.3"})
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        dna = str(data.get("dna", "")).upper()
        if len(dna) == (end1 - start0) and dna[flank_bp : flank_bp + len(ref)] == ref:
            mut_dna = dna[:flank_bp] + alt + dna[flank_bp + len(ref) :]
            return {"wt": dna, "mut": mut_dna}
    except Exception:
        pass

    c_num = chrom[3:] if chrom.startswith("chr") else chrom
    ens_url = f"https://rest.ensembl.org/sequence/region/human/{c_num}:{pos - flank_bp}..{pos + len(ref) - 1 + flank_bp}:1?content-type=application/json"
    req = urllib.request.Request(ens_url, headers={"User-Agent": "NC-MutScan/3.3"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    dna = str(data.get("seq", "")).upper()
    assert dna[flank_bp : flank_bp + len(ref)] == ref, f"Ref mismatch at {chrom}:{pos}: expected {ref}, got {dna[flank_bp:flank_bp+len(ref)]}"
    mut_dna = dna[:flank_bp] + alt + dna[flank_bp + len(ref) :]
    return {"wt": dna, "mut": mut_dna}


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


def main():
    os.makedirs(os.path.join("results", "major_revisions"), exist_ok=True)
    dev_cache_path = os.path.join("results", "major_revisions", "dev_flanks_cache.json")

    with open(os.path.join("results", "pi_upgrades", "blind_flanks_cache.json"), "r", encoding="utf-8") as ff:
        raw_blind = json.load(ff)
    flanks = {}
    for k, val in raw_blind.items():
        flanks[str(k).strip()] = val
        flanks[norm_key(k)] = val

    if os.path.isfile(dev_cache_path):
        with open(dev_cache_path, "r", encoding="utf-8") as df:
            raw_dev = json.load(df)
        for k, val in raw_dev.items():
            flanks[str(k).strip()] = val
            flanks[norm_key(k)] = val
    else:
        raw_dev = {}

    dev_vars = load_vcf("dev_chr1_11.vcf", "dev_chr1_11_X")
    blind_vars = load_vcf("locked_blind_chr12_22.vcf", "blind_chr12_22")

    fetched = 0
    for v in dev_vars:
        if v["nk_full"] not in flanks and v["nk_cp"] not in flanks:
            entry = fetch_hg38_flank(v["chrom"], v["pos"], v["ref"], v["alt"], flank_bp=25)
            raw_dev[v["nk_full"]] = entry
            flanks[v["nk_full"]] = entry
            flanks[v["nk_cp"]] = entry
            fetched += 1
            time.sleep(0.08)
    if fetched:
        with open(dev_cache_path, "w", encoding="utf-8") as df:
            json.dump(raw_dev, df, indent=2)
        print(f"Cached {fetched} hg38 flanking sequences -> {dev_cache_path}")

    spec = importlib.util.spec_from_file_location("frozen_mutscan", "nc_mutscan_research_suite.py")
    ms_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ms_mod)

    evaluated_all = []
    for v in dev_vars + blind_vars:
        entry = flanks.get(v["nk_full"], flanks.get(v["vid"], flanks.get(v["nk_cp"])))
        res = ms_mod.analyze_variant_delta(wt_seq=entry["wt"], mut_seq=entry["mut"], flank_bp=25)
        pred_v33 = is_frozen_driver_call(res)
        cls_txt = str(res.get("classification", "")).lower()
        raw_dpwm = abs(float(res.get("delta_pwm", 0.0) or 0.0))
        is_non_tf_hit = 1 if (pred_v33 == 1 and any(k in cls_txt for k in ("splice", "polya"))) else 0
        is_masked = 1 if any(k in cls_txt for k in ("aid", "passenger", "compet")) else 0
        evaluated_all.append({
            "split": v["split"],
            "chrom": v["chrom"],
            "label": v["label"],
            "pred_v33": pred_v33,
            "is_non_tf_hit": is_non_tf_hit,
            "is_masked": is_masked,
            "raw_dpwm": raw_dpwm,
        })

    chroms = sorted(set(x["chrom"] for x in evaluated_all), key=lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, c))
    dev_only = [x for x in evaluated_all if x["split"] == "dev_chr1_11_X"]
    m_dev = metrics_dict([x["label"] for x in dev_only], [x["pred_v33"] for x in dev_only])
    m_all = metrics_dict([x["label"] for x in evaluated_all], [x["pred_v33"] for x in evaluated_all])

    print(f"\n=== FULL 23-CHROMOSOME LOCO CROSS-VALIDATION ({len(chroms)} CHROMOSOMES, N={len(evaluated_all)}) ===")
    print(f"  Dev Set Only (chr1-11,X; N={m_dev['N']}): TP={m_dev['TP']}, FP={m_dev['FP']}, TN={m_dev['TN']}, FN={m_dev['FN']} | Sens={m_dev['Sens']} | Spec={m_dev['Spec']} | PPV={m_dev['PPV']} | F1={m_dev['F1']}")
    print(f"  Pooled 23-Chr LOCO (Fixed v3.3; N={m_all['N']}): TP={m_all['TP']}, FP={m_all['FP']}, TN={m_all['TN']}, FN={m_all['FN']} | Sens={m_all['Sens']} | Spec={m_all['Spec']} | PPV={m_all['PPV']} | F1={m_all['F1']}")

    tau_grid = [round(3.5 + 0.1 * i, 1) for i in range(16)]
    oof_tuned_preds, oof_tuned_labels, chosen_taus = [], [], []
    for held_chr in chroms:
        train = [x for x in evaluated_all if x["chrom"] != held_chr]
        test = [x for x in evaluated_all if x["chrom"] == held_chr]
        best_tau, best_obj = 4.2, -1.0
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
            oof_tuned_labels.append(x["label"])
            oof_tuned_preds.append(1 if (x["is_non_tf_hit"] or (not x["is_masked"] and x["raw_dpwm"] >= best_tau)) else 0)

    m_nest = metrics_dict(oof_tuned_labels, oof_tuned_preds)
    print(f"  Nested Fold-Tuned 23-Chr LOCO (mode tau_TF={Counter(chosen_taus).most_common(1)[0][0]} bits): TP={m_nest['TP']}, FP={m_nest['FP']}, TN={m_nest['TN']}, FN={m_nest['FN']} | Sens={m_nest['Sens']} | Spec={m_nest['Spec']} | PPV={m_nest['PPV']} | F1={m_nest['F1']}")


if __name__ == "__main__":
    main()
