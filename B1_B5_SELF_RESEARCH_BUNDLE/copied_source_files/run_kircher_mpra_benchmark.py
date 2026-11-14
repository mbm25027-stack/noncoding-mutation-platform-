"""
run_kircher_mpra_benchmark.py
Evaluates frozen NC-MutScan v3.3 across all GRCh38 regulatory elements in the cached
Kircher et al. (2019, Nat. Commun. 10:3583) dataset (elements.tsv.gz).
"""

import csv
import gzip
import importlib.util
import json
import math
import os
import urllib.request
from collections import defaultdict

csv.field_size_limit(50_000_000)


def fetch_seq_window(genome, chrom, start1, end1):
    c = chrom if str(chrom).lower().startswith("chr") else f"chr{chrom}"
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome={genome};chrom={c};start={start1 - 1};end={end1}"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "NC-MutScan-Benchmark/3.3"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
        dna = str(data.get("dna", "")).upper()
        if len(dna) == (end1 - start1 + 1):
            return dna
    except Exception:
        pass

    c_num = c[3:] if c.startswith("chr") else c
    ens_url = f"https://rest.ensembl.org/sequence/region/human/{c_num}:{start1}..{end1}:1?content-type=application/json"
    req = urllib.request.Request(ens_url, headers={"User-Agent": "NC-MutScan-Benchmark/3.3"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    return str(data.get("seq", "")).upper()


def is_driver_call(res):
    cls_txt = str(res.get("classification", "")).lower()
    if any(neg in cls_txt for neg in ("benign", "neutral", "polymorphism", "passenger", "sub-threshold", "subthreshold", "non-driver", "compet")):
        return 0
    return 1


def wilson_ci(k, n, z=1.959963984540054):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1.0 + z * z / n
    center = p + z * z / (2.0 * n)
    rad = z * math.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n)
    return max(0.0, (center - rad) / denom) * 100.0, min(1.0, (center + rad) / denom) * 100.0


def main():
    local_gz = os.path.join("results", "major_revisions", "kircher_elements_cache", "elements.tsv.gz")
    spec = importlib.util.spec_from_file_location("frozen_mutscan", "nc_mutscan_research_suite.py")
    ms_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ms_mod)

    grouped_snvs = defaultdict(list)
    with gzip.open(local_gz, "rt", encoding="utf-8", errors="replace") as gz_in:
        reader = csv.DictReader(gz_in, delimiter="\t")
        for r in reader:
            rel = str(r.get("Release", "")).strip()
            if rel and rel.upper() != "GRCH38":
                continue
            elem = str(r.get("Element", "")).strip()
            chrom = str(r.get("Chrom", "")).strip()
            pos_s = str(r.get("Pos", "")).strip()
            ref = str(r.get("Ref", "")).strip().upper()
            alt = str(r.get("Alt", "")).strip().upper()
            bc_s = str(r.get("Barcodes", "0")).strip()
            coef_s = str(r.get("Coefficient", "")).strip()
            pval_s = str(r.get("pValue", "1.0")).strip()

            if not elem or not chrom or not pos_s.isdigit():
                continue
            if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
                continue
            try:
                if int(float(bc_s or 0)) < 10:
                    continue
                coef_f = float(coef_s)
                pval_f = float(pval_s) if pval_s not in ("", "NA", "None", ".") else 1.0
            except ValueError:
                continue

            c_clean = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
            grouped_snvs[(elem, c_clean)].append({
                "chrom": c_clean,
                "pos": int(pos_s),
                "ref": ref,
                "alt": alt,
                "val": coef_f,
                "pval": pval_f,
            })

    print(f"=== ITEM 2: KIRCHER ET AL. (2019) SATURATION MUTAGENESIS MPRA ({len(grouped_snvs)} ELEMENTS) ===")
    summary_rows = []
    tot_snvs = tot_tp = tot_fp = tot_tn = tot_fn = 0
    tot_act_tp = tot_act_tot = 0

    for (elem_name, chrom), snv_rows in sorted(grouped_snvs.items()):
        # Deduplicate physical SNVs within each element across cell-line assays, retaining peak |log2| effect
        dedup = {}
        for x in snv_rows:
            k = (x["pos"], x["ref"], x["alt"])
            if k not in dedup or abs(x["val"]) > abs(dedup[k]["val"]):
                dedup[k] = x
        uniq_snvs = list(dedup.values())

        min_pos = min(x["pos"] for x in uniq_snvs)
        max_pos = max(x["pos"] for x in uniq_snvs)
        if (max_pos - min_pos) > 50000:
            continue
        win_start = min_pos - 25
        win_end = max_pos + 25

        genome_used = "hg38"
        ref_seq = fetch_seq_window("hg38", chrom, win_start, win_end)
        matches = sum(1 for x in uniq_snvs[:30] if len(ref_seq) > (x["pos"] - win_start) and ref_seq[x["pos"] - win_start] == x["ref"])
        if matches < min(10, len(uniq_snvs[:30])):
            genome_used = "hg19"
            ref_seq = fetch_seq_window("hg19", chrom, win_start, win_end)

        tp = fp = tn = fn = n_eval = 0
        act_tp = act_fn = 0

        for x in uniq_snvs:
            idx = x["pos"] - win_start
            if idx < 25 or idx + 26 > len(ref_seq) or ref_seq[idx] != x["ref"]:
                continue

            # High-confidence functional vs. neutral tiers (Barcodes >= 10):
            # Functional regulatory variant: |Coefficient| >= 0.30 and pValue < 1e-5
            # Neutral control variant:       |Coefficient| <  0.10 and pValue > 0.10
            if abs(x["val"]) >= 0.30 and x["pval"] < 1e-5:
                y_true = 1
            elif abs(x["val"]) < 0.10 and x["pval"] > 0.10:
                y_true = 0
            else:
                continue

            wt_flank = ref_seq[idx - 25 : idx + 26]
            mut_flank = wt_flank[:25] + x["alt"] + wt_flank[26:]
            res = ms_mod.analyze_variant_delta(wt_seq=wt_flank, mut_seq=mut_flank, flank_bp=25)
            y_pred = is_driver_call(res)

            n_eval += 1
            if y_true == 1 and y_pred == 1:
                tp += 1
            elif y_true == 0 and y_pred == 1:
                fp += 1
            elif y_true == 0 and y_pred == 0:
                tn += 1
            else:
                fn += 1

            if x["val"] >= 0.30 and x["pval"] < 1e-5:
                if y_pred == 1:
                    act_tp += 1
                else:
                    act_fn += 1

        if n_eval == 0:
            continue

        spec_pct = 100.0 * tn / (tn + fp) if (tn + fp) else 0.0
        sens_pct = 100.0 * tp / (tp + fn) if (tp + fn) else 0.0
        ppv_pct = 100.0 * tp / (tp + fp) if (tp + fp) else 0.0
        tot_snvs += len(uniq_snvs)
        tot_tp += tp; tot_fp += fp; tot_tn += tn; tot_fn += fn
        tot_act_tp += act_tp; tot_act_tot += (act_tp + act_fn)

        summary_rows.append({
            "Element": elem_name,
            "Chrom": chrom,
            "Build": genome_used,
            "Unique_QC_SNVs": len(uniq_snvs),
            "Evaluated_Tiered_SNVs": n_eval,
            "TP": tp, "FP": fp, "TN": tn, "FN": fn,
            "Activating_TP_of_Total": f"{act_tp}/{act_tp + act_fn}",
            "Sensitivity_%": f"{sens_pct:.2f}",
            "Specificity_%": f"{spec_pct:.2f}",
            "Precision_PPV_%": f"{ppv_pct:.2f}",
        })
        print(f"  {elem_name:22s} ({chrom:5s}, {genome_used}) | QC={len(uniq_snvs):4d} | Eval={n_eval:4d} | TP={tp:2d}, FP={fp:2d}, TN={tn:4d}, FN={fn:3d} | Sens={sens_pct:5.2f}% | Spec={spec_pct:6.2f}% | PPV={ppv_pct:6.2f}%", flush=True)

    out_csv = os.path.join("results", "major_revisions", "kircher_mpra_multi_element_benchmark.csv")
    if summary_rows:
        with open(out_csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
            w.writeheader()
            w.writerows(summary_rows)

    tot_spec = 100.0 * tot_tn / (tot_tn + tot_fp) if (tot_tn + tot_fp) else 0.0
    tot_sens = 100.0 * tot_tp / (tot_tp + tot_fn) if (tot_tp + tot_fn) else 0.0
    tot_ppv = 100.0 * tot_tp / (tot_tp + tot_fp) if (tot_tp + tot_fp) else 0.0
    sp_lo, sp_hi = wilson_ci(tot_tn, tot_tn + tot_fp)
    print(f"\n  POOLED ACROSS {len(summary_rows)} ELEMENTS ({tot_snvs} Unique QC SNVs): TP={tot_tp}, FP={tot_fp}, TN={tot_tn}, FN={tot_fn}")
    print(f"  Pooled Sensitivity = {tot_sens:.2f}% | Pooled Specificity = {tot_spec:.2f}% ({sp_lo:.1f}-{sp_hi:.1f}%) | Pooled PPV = {tot_ppv:.2f}%")
    print(f"  Strong Activating Alleles (Coefficient>=+0.30, pValue<1e-5): {tot_act_tp}/{tot_act_tot} recovered")
    print(f"  Saved table -> {out_csv}")


if __name__ == "__main__":
    main()
