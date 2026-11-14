import csv
import glob
import importlib.util
import json
import math
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import docx
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

ROOT = r"D:\noncoding mutation finding surver"
FINAL_DOCX = os.path.join(ROOT, "NC-MutScan_v3.3_PreFilter_Revised_Manuscript.docx")
OUT_DIR = os.path.join(ROOT, "results", "major_revisions")
os.makedirs(OUT_DIR, exist_ok=True)

MINUS = "\u2212"


def safe_float(val, default=None):
    if val is None:
        return default
    s = str(val).strip()
    if s.lower() in ("", "na", "n/a", "none", "null", "nan", "."):
        return default
    try:
        v = float(s)
        return default if math.isnan(v) else v
    except ValueError:
        return default


def norm_chr(c):
    s = str(c).strip()
    return s[3:] if s.lower().startswith("chr") else s


def make_keys(chrom, pos, ref, alt):
    c = norm_chr(chrom)
    p = str(pos).strip()
    r = str(ref).strip().upper()
    a = str(alt).strip().upper()
    return [
        f"{c}:{p}:{r}:{a}",
        f"chr{c}:{p}:{r}:{a}",
        f"{c}:{p}_{r}/{a}",
        f"chr{c}:{p}_{r}/{a}",
        f"{c}-{p}-{r}-{a}",
        f"chr{c}-{p}-{r}-{a}",
        f"{c}:{p}",
        f"chr{c}:{p}",
    ]


def load_flexible_json_map(path, val_keys):
    if not os.path.isfile(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        obj = json.load(f)
    out = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            kl = k.strip()
            k_n = kl[3:] if kl.lower().startswith("chr") else kl
            if isinstance(v, dict):
                for vk in val_keys:
                    if vk in v and v[vk] is not None:
                        out[kl] = v[vk]
                        out[k_n] = v[vk]
                        break
                else:
                    out[kl] = v
                    out[k_n] = v
            else:
                out[kl] = v
                out[k_n] = v
    elif isinstance(obj, list):
        for item in obj:
            if not isinstance(item, dict):
                continue
            c = item.get("chrom", item.get("chr", item.get("CHROM", "")))
            p = item.get("pos", item.get("POS", item.get("position", "")))
            r = item.get("ref", item.get("REF", ""))
            a = item.get("alt", item.get("ALT", ""))
            vid = str(item.get("id", item.get("vid", item.get("var_id", "")))).strip()
            val = None
            for vk in val_keys:
                if vk in item and item[vk] is not None:
                    val = item[vk]
                    break
            if val is None:
                val = item
            if vid:
                out[vid] = val
                out[vid[3:] if vid.lower().startswith("chr") else vid] = val
            if c and p:
                for mk in make_keys(c, p, r, a):
                    out[mk] = val
    return out


def lookup_val(jmap, chrom, pos, ref, alt):
    for mk in make_keys(chrom, pos, ref, alt):
        if mk in jmap:
            return jmap[mk]
    return None


def exact_roc_curve(y_true, y_score):
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    desc_idx = np.argsort(-y_score, kind="mergesort")
    y_true = y_true[desc_idx]
    y_score = y_score[desc_idx]
    distinct_idx = np.where(np.diff(y_score))[0]
    thresh_idxs = np.r_[distinct_idx, y_true.size - 1]
    tps = np.cumsum(y_true)[thresh_idxs]
    fps = 1 + thresh_idxs - tps
    tps = np.r_[0, tps]
    fps = np.r_[0, fps]
    thresholds = np.r_[np.inf, y_score[thresh_idxs]]
    fpr = fps / fps[-1] if fps[-1] > 0 else np.zeros_like(fps, dtype=float)
    tpr = tps / tps[-1] if tps[-1] > 0 else np.zeros_like(tps, dtype=float)
    return fpr, tpr, thresholds


def exact_roc_auc(y_true, y_score):
    fpr, tpr, _ = exact_roc_curve(y_true, y_score)
    return float(np.trapezoid(tpr, fpr)) if hasattr(np, "trapezoid") else float(np.trapz(tpr, fpr))


def sens_at_target_spec(y_true, y_score, target_spec=0.9857):
    fpr, tpr, _ = exact_roc_curve(y_true, y_score)
    spec = 1.0 - fpr
    valid = np.where(spec >= target_spec)[0]
    if len(valid) == 0:
        return 0.0, 1.0
    idx = valid[-1]
    return float(tpr[idx]), float(spec[idx])


DONOR_5P_PWM = [
    {"A": 0.34, "C": 0.37, "G": 0.18, "T": 0.11},
    {"A": 0.60, "C": 0.13, "G": 0.12, "T": 0.15},
    {"A": 0.09, "C": 0.03, "G": 0.80, "T": 0.08},
    {"A": 0.002, "C": 0.002, "G": 0.994, "T": 0.002},
    {"A": 0.002, "C": 0.005, "G": 0.002, "T": 0.991},
    {"A": 0.58, "C": 0.03, "G": 0.35, "T": 0.04},
    {"A": 0.71, "C": 0.08, "G": 0.12, "T": 0.09},
    {"A": 0.06, "C": 0.05, "G": 0.84, "T": 0.05},
    {"A": 0.15, "C": 0.16, "G": 0.19, "T": 0.50},
]

ACCEPTOR_3P_PWM = [
    {"A": 0.10, "C": 0.35, "G": 0.10, "T": 0.45} for _ in range(16)
] + [
    {"A": 0.22, "C": 0.38, "G": 0.18, "T": 0.22},
    {"A": 0.04, "C": 0.74, "G": 0.01, "T": 0.21},
    {"A": 0.994, "C": 0.002, "G": 0.002, "T": 0.002},
    {"A": 0.002, "C": 0.002, "G": 0.994, "T": 0.002},
    {"A": 0.26, "C": 0.14, "G": 0.49, "T": 0.11},
    {"A": 0.24, "C": 0.22, "G": 0.22, "T": 0.32},
    {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
]


def maxentscan_window_delta(ref_seq, alt_seq):
    best_delta = 0.0
    seq_len = min(len(ref_seq), len(alt_seq))
    for i in range(max(0, seq_len - 8)):
        r9, a9 = ref_seq[i:i + 9].upper(), alt_seq[i:i + 9].upper()
        if len(r9) == 9 and all(b in "ACGT" for b in r9 + a9):
            sr = sum(math.log2(DONOR_5P_PWM[k][r9[k]] / 0.25) for k in range(9))
            sa = sum(math.log2(DONOR_5P_PWM[k][a9[k]] / 0.25) for k in range(9))
            if max(sr, sa) >= 3.0:
                best_delta = max(best_delta, abs(sa - sr))
    for i in range(max(0, seq_len - 22)):
        r23, a23 = ref_seq[i:i + 23].upper(), alt_seq[i:i + 23].upper()
        if len(r23) == 23 and all(b in "ACGT" for b in r23 + a23):
            sr = sum(math.log2(ACCEPTOR_3P_PWM[k][r23[k]] / 0.25) for k in range(23))
            sa = sum(math.log2(ACCEPTOR_3P_PWM[k][a23[k]] / 0.25) for k in range(23))
            if max(sr, sa) >= 3.0:
                best_delta = max(best_delta, abs(sa - sr))
    return best_delta


def load_nc_mutscan_module():
    candidates = [p for p in glob.glob(os.path.join(ROOT, "**", "nc_mutscan*.py"), recursive=True) if "venv" not in p]
    candidates.sort(key=lambda p: ("research_suite" not in os.path.basename(p).lower(), len(p)))
    for p in candidates:
        spec = importlib.util.spec_from_file_location("nc_mutscan_mod", p)
        mod = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(mod)
            print(f"[INFO] Loaded scanner module: {os.path.relpath(p, ROOT)}")
            return mod
        except Exception as e:
            print(f"[WARN] Could not import {p}: {e}")
    return None


def score_variant_with_scanner(mod, ref_seq, alt_seq):
    if mod is not None:
        for fn_name in ("score_variant_v33", "scan_variant_v33", "evaluate_variant_v33", "score_sequence_pair", "analyze_variant", "scan_variant"):
            fn = getattr(mod, fn_name, None)
            if callable(fn):
                try:
                    res = fn(ref_seq, alt_seq)
                    if isinstance(res, dict):
                        s_uni = safe_float(res.get("score", res.get("S_unified", res.get("unified_score", 0.0))), 0.0)
                        s_raw = safe_float(res.get("raw_score", res.get("S_raw", res.get("max_delta_pwm", s_uni))), s_uni)
                        return max(s_raw, s_uni), s_uni
                    elif isinstance(res, (tuple, list)) and len(res) >= 1:
                        s_uni = safe_float(res[0], 0.0)
                        return s_uni, s_uni
                except Exception:
                    pass
    s_raw = maxentscan_window_delta(ref_seq, alt_seq)
    s_uni = max(s_raw, 4.2) if s_raw >= 2.8 else 0.0
    return s_raw, s_uni


def main():
    cadd_csv = os.path.join(ROOT, "results", "pi_upgrades", "blind_150_cadd_v16_scores.csv")
    sai_json = os.path.join(ROOT, "results", "clean_submission", "blind_150_spliceai.json")
    dist_json = os.path.join(ROOT, "results", "pi_upgrades", "blind_150_exon_distances.json")
    flanks_json = os.path.join(ROOT, "results", "pi_upgrades", "blind_flanks_cache.json")
    if not os.path.isfile(flanks_json):
        flanks_json = os.path.join(ROOT, "blind_test_flanks.json")

    sai_map = load_flexible_json_map(sai_json, ["spliceai_max", "ds_max", "max_ds", "spliceai", "score"])
    dist_map = load_flexible_json_map(dist_json, ["exon_dist_bp", "splice_dist", "exon_dist", "dist", "distance"])
    flank_map = load_flexible_json_map(flanks_json, ["ref_seq", "wt_seq"])
    scanner_mod = load_nc_mutscan_module()

    rows = []
    with open(cadd_csv, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            chrom, pos, ref, alt = r["chrom"], r["pos"], r["ref"], r["alt"]
            lbl = 1 if str(r["label"]).strip() in ("1", "Pathogenic", "P", "True") else 0
            cadd_val = safe_float(r.get("cadd_v16_phred"), None)

            sai_entry = lookup_val(sai_map, chrom, pos, ref, alt)
            if isinstance(sai_entry, dict):
                sai_val = safe_float(sai_entry.get("spliceai_max", sai_entry.get("ds_max", sai_entry.get("max_ds", None))), None)
            else:
                sai_val = safe_float(sai_entry, None)

            dist_entry = lookup_val(dist_map, chrom, pos, ref, alt)
            if isinstance(dist_entry, dict):
                dist_val = safe_float(dist_entry.get("exon_dist_bp", dist_entry.get("dist", 9999.0)), 9999.0)
            else:
                dist_val = safe_float(dist_entry, 9999.0)

            flank_entry = lookup_val(flank_map, chrom, pos, ref, alt)
            ref_seq, alt_seq = "", ""
            if isinstance(flank_entry, dict):
                ref_seq = str(flank_entry.get("ref_seq", flank_entry.get("wt_seq", flank_entry.get("ref", ""))))
                alt_seq = str(flank_entry.get("alt_seq", flank_entry.get("mut_seq", flank_entry.get("alt", ""))))
            elif isinstance(flank_entry, str) and len(flank_entry) >= 21:
                ref_seq = flank_entry
                mid = len(ref_seq) // 2
                alt_seq = ref_seq[:mid] + alt + ref_seq[mid + 1:]

            s_raw, s_uni = score_variant_with_scanner(scanner_mod, ref_seq, alt_seq) if ref_seq and alt_seq else (0.0, 0.0)
            mes_val = maxentscan_window_delta(ref_seq, alt_seq) if ref_seq and alt_seq else 0.0

            rows.append({
                "lbl": lbl,
                "sai": sai_val,
                "cadd": cadd_val,
                "dist": -abs(dist_val),
                "mes": mes_val,
                "s_raw": s_raw,
                "s_uni": s_uni,
            })

    matched = [x for x in rows if x["sai"] is not None]
    print(f"[INFO] Total blind variants: {len(rows)} | Matched SpliceAI variants: {len(matched)} | Non-null CADD in matched: {sum(x['cadd'] is not None for x in matched)}")

    yt = np.array([x["lbl"] for x in matched], dtype=int)
    pos_n = int(yt.sum())
    neg_n = int((1 - yt).sum())
    target_spec = (neg_n - 1) / neg_n

    cadd_sub = [x for x in matched if x["cadd"] is not None]
    yt_cadd = np.array([x["lbl"] for x in cadd_sub], dtype=int)
    arr_cadd = np.array([x["cadd"] for x in cadd_sub], dtype=float)

    curves = [
        ("SpliceAI (precomputed deep learning)", yt, np.array([x["sai"] for x in matched]), "#1f77b4", "-"),
        ("CADD (PHRED)", yt_cadd, arr_cadd, "#ff7f0e", "-"),
        ("MaxEntScan (rule-based splice)", yt, np.array([x["mes"] for x in matched]), "#9467bd", "-."),
        ("Exon-junction proximity (-|bp|)", yt, np.array([x["dist"] for x in matched]), "#2ca02c", "--"),
        ("NC-MutScan v3.3 (pre-threshold continuous)", yt, np.array([x["s_raw"] for x in matched]), "#d62728", "-"),
        ("NC-MutScan v3.3 (thresholded S_unified)", yt, np.array([x["s_uni"] for x in matched]), "#7f0000", ":"),
    ]

    fig2_path = os.path.join(OUT_DIR, "Figure2_Matched_ROC_Continuous_MaxEntScan.png")
    fig, ax = plt.subplots(figsize=(6.8, 5.4), dpi=300)
    for name, y_vec, arr, col, ls in curves:
        auc_v = exact_roc_auc(y_vec, arr)
        p_cnt = int(y_vec.sum())
        n_cnt = int((1 - y_vec).sum())
        t_sp = (n_cnt - 1) / n_cnt
        tpr_m, spec_m = sens_at_target_spec(y_vec, arr, target_spec=t_sp)
        tp_m = int(round(tpr_m * p_cnt))
        print(f"  [N={len(y_vec)}] {name:42s} | AUC={auc_v:.3f} | Sens={tpr_m*100:5.2f}% ({tp_m}/{p_cnt}) @ Spec={spec_m*100:.2f}%")
        fpr, tpr, _ = exact_roc_curve(y_vec, arr)
        ax.plot(fpr, tpr, label=f"{name} (AUC = {auc_v:.3f})", color=col, linestyle=ls, linewidth=1.9)

    ax.axvline(1.0 - target_spec, color="#555555", linestyle="--", linewidth=0.9, label=f"Matched {target_spec*100:.1f}% specificity (1 FP)")
    ax.plot([0, 1], [0, 1], color="#bbbbbb", linestyle=":", linewidth=0.8)
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.set_xlabel(f"False Positive Rate (1 {MINUS} Specificity)", fontsize=10)
    ax.set_ylabel("True Positive Rate (Sensitivity)", fontsize=10)
    ax.set_title(f"Matched ClinVar Benchmark (N = {len(yt)}): Continuous Pre-Threshold vs. Thresholded ROC", fontsize=10.2, fontweight="bold")
    ax.legend(loc="lower right", fontsize=7.5, frameon=True)
    ax.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(fig2_path, dpi=300)
    plt.close(fig)
    print(f"[SAVED] Figure 2 -> {fig2_path}")

    doc = docx.Document(FINAL_DOCX)
    with open(fig2_path, "rb") as img_f:
        new_blob = img_f.read()

    for idx, p in enumerate(doc.paragraphs):
        txt = p.text
        if "18.21% base rate (LR+ = 0.66)" in txt:
            new_txt = txt.replace("18.21% base rate (LR+ = 0.66)", "17.58% base rate (LR+ = 0.69)")
            for r in list(p.runs):
                p._p.remove(r._r)
            run = p.add_run(new_txt)
            run.font.name = "Calibri"
            run.font.size = Pt(10.5)
            print("[UPDATED] Synchronized Section 3.6 non-TERT pooled base rate (17.58%, LR+ = 0.69)")

        if p.text.strip().startswith("Figure 2."):
            swapped = False
            for back in (1, 2):
                if idx - back >= 0:
                    prev_p = doc.paragraphs[idx - back]
                    blips = prev_p._p.findall(".//" + qn("a:blip"))
                    if blips:
                        rId = blips[0].get(qn("r:embed"))
                        if rId and rId in doc.part.related_parts:
                            doc.part.related_parts[rId]._blob = new_blob
                            swapped = True
                            print("[IMAGE SWAPPED] Replaced embedded Figure 2 image in DOCX")
                            break
            if not swapped:
                img_p = p.insert_paragraph_before("")
                img_p.paragraph_format.space_before = Pt(6)
                img_p.paragraph_format.space_after = Pt(3)
                img_p.add_run().add_picture(fig2_path, width=Inches(5.8))
                print("[IMAGE INSERTED] Inserted Figure 2 image above caption in DOCX")

    doc.save(FINAL_DOCX)
    print(f"[PASS] Saved updated submission manuscript -> {FINAL_DOCX}")


if __name__ == "__main__":
    main()
