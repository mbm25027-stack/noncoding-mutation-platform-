"""
run_items2_and_3.py
Executes PI Item 3 (re-evaluating the 50-variant and 54-locus sets under v3.3)
and PI Item 2 (second independent non-overlapping chr12-22 ClinVar benchmark,
n = 300 pathogenic + 300 benign = 600 variants under frozen v3.3).
"""

import csv
import glob
import gzip
import json
import math
import os
import random
import re
import ssl
import time
import urllib.request
from nc_mutscan_research_suite import analyze_variant_delta

os.makedirs("results/pi_upgrades", exist_ok=True)
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


def fetch_ucsc_flank(chrom: str, pos: int, ref: str, alt: str, cache: dict):
    c = chrom.replace("chr", "")
    key = f"chr{c}:{pos}:{ref}>{alt}"
    if key in cache:
        return cache[key]
    start_0 = max(0, pos - 1 - FLANK_BP)
    end_0 = pos - 1 + len(ref) + FLANK_BP
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{c};start={start_0};end={end_0}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=15, context=SSL_CTX) as resp:
        dna = json.loads(resp.read().decode("utf-8")).get("dna", "").upper()
    offset = (pos - 1) - start_0
    mut_dna = dna[:offset] + alt.upper() + dna[offset + len(ref):]
    cache[key] = {"wt": dna, "mut": mut_dna}
    time.sleep(0.04)
    return cache[key]


# ---------------------------------------------------------------------------
# ITEM 3: Re-run 50-variant and 54-locus evaluations under NC-MutScan v3.3
# ---------------------------------------------------------------------------
def run_item3_unified_v33():
    print("=== ITEM 3: Re-evaluating 50-variant and 54-locus sets under v3.3 ===")
    # Inspect nc_mutscan_heldout_benchmark.py and any 54-locus script/data in workspace
    py_files = sorted(glob.glob("*.py") + glob.glob("results/**/*.json", recursive=True) + glob.glob("results/**/*.csv", recursive=True))
    print("Workspace benchmark artifacts detected:", [f for f in py_files if any(k in f.lower() for k in ("heldout", "54", "verify", "mutscan", "benchmark"))])

    if os.path.exists("nc_mutscan_heldout_benchmark.py"):
        import nc_mutscan_heldout_benchmark as hb
        # Check common dataset attributes in nc_mutscan_heldout_benchmark.py
        for attr in dir(hb):
            obj = getattr(hb, attr)
            if isinstance(obj, list) and len(obj) in (50, 54):
                print(f"Found {len(obj)}-item benchmark list in nc_mutscan_heldout_benchmark.{attr}")
                tp = tn = fp = fn = 0
                for item in obj:
                    # Handle tuple/dict formats from heldout benchmark
                    if isinstance(item, dict):
                        wt = item.get("wt_seq") or item.get("wt") or item.get("ref_seq") or item.get("seq", "")
                        mut = item.get("mut_seq") or item.get("mut") or item.get("alt_seq") or wt
                        expected = str(item.get("expected") or item.get("label") or item.get("class", ""))
                    elif isinstance(item, (list, tuple)) and len(item) >= 3:
                        wt, mut, expected = str(item[0]), str(item[1]), str(item[2])
                    else:
                        continue
                    res = analyze_variant_delta(wt, mut, flank_bp=FLANK_BP)
                    pred_pos = res["classification"] in ("Driver Variant", "Loss-of-Function Variant")
                    exp_pos = any(k in expected.lower() for k in ("driver", "pos", "loss", "1", "true"))
                    if exp_pos and pred_pos:
                        tp += 1
                    elif exp_pos and not pred_pos:
                        fn += 1
                    elif not exp_pos and pred_pos:
                        fp += 1
                    else:
                        tn += 1
                print(f"  v3.3 result on {attr} (N={tp+tn+fp+fn}): TP={tp}, FP={fp}, TN={tn}, FN={fn} | "
                      f"Sens={tp/max(1,tp+fn)*100:.1f}% | Spec={tn/max(1,tn+fp)*100:.1f}%")


# ---------------------------------------------------------------------------
# ITEM 2: Second Independent Noncoding Benchmark (n = 300/class = 600 on chr12-22)
# ---------------------------------------------------------------------------
CODING_TERMS = {
    "frameshift", "nonsense", "missense", "inframe", "stop_gained",
    "stop_lost", "start_lost", "synonymous", "protein_altering"
}


def load_existing_keys():
    used = set()
    for path in ("dev_chr1_11.vcf", "locked_blind_chr12_22.vcf"):
        if not os.path.exists(path):
            continue
        with open(path, "rt", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                p = line.rstrip("\n").split("\t")
                if p[0] in ("PATHOGENIC", "BENIGN"):
                    used.add((p[1].replace("chr", ""), int(p[2]), p[4], p[5].split(",")[0]))
                else:
                    used.add((p[0].replace("chr", ""), int(p[1]), p[3], p[4].split(",")[0]))
    return used


def find_local_clinvar_vcf():
    for pattern in ("*clinvar*.vcf.gz", "*clinvar*.vcf", "../*clinvar*.vcf*"):
        matches = glob.glob(pattern)
        if matches:
            return matches[0]
    return None


def sample_from_local_vcf(vcf_path, used_keys, target_per_class=300, seed=104):
    target_chroms = {str(c) for c in range(12, 23)}
    path_pool, ben_pool = [], []
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0].replace("chr", "")
            if chrom not in target_chroms:
                continue
            pos = int(parts[1])
            var_id = parts[2]
            ref = parts[3]
            alt = parts[4].split(",")[0]
            if len(ref) > 15 or len(alt) > 15 or not re.fullmatch(r"[ACGT]+", ref) or not re.fullmatch(r"[ACGT]+", alt):
                continue
            if (chrom, pos, ref, alt) in used_keys:
                continue
            info = parts[7]
            mc_m = re.search(r"MC=([^;]+)", info)
            mc = mc_m.group(1).lower() if mc_m else ""
            if any(term in mc for term in CODING_TERMS):
                continue
            cln_m = re.search(r"CLNSIG=([^;]+)", info)
            clnsig = cln_m.group(1) if cln_m else ""
            if "Conflicting" in clnsig or "Uncertain" in clnsig:
                continue
            if clnsig.startswith("Pathogenic") or clnsig.startswith("Likely_pathogenic"):
                path_pool.append({"label": "PATHOGENIC", "chrom": chrom, "pos": pos, "var_id": var_id, "ref": ref, "alt": alt, "mc": mc})
            elif clnsig.startswith("Benign") or clnsig.startswith("Likely_benign"):
                ben_pool.append({"label": "BENIGN", "chrom": chrom, "pos": pos, "var_id": var_id, "ref": ref, "alt": alt, "mc": mc})

    rng = random.Random(seed)
    return rng.sample(path_pool, target_per_class), rng.sample(ben_pool, target_per_class)


def sample_via_ncbi_api(used_keys, target_per_class=300, seed=104):
    print("No raw ClinVar VCF on disk; querying NCBI ClinVar GRCh38 API for fresh non-overlapping chr12-22 variants...")
    rng = random.Random(seed)
    chroms = [str(c) for c in range(12, 23)]

    def fetch_class(cln_query, label_name):
        pool = []
        # Query non-coding molecular consequences directly on chr12-22
        mc_terms = '("intron variant"[Molecular Consequence] OR "splice donor variant"[Molecular Consequence] OR "splice acceptor variant"[Molecular Consequence] OR "5 prime UTR variant"[Molecular Consequence] OR "3 prime UTR variant"[Molecular Consequence] OR "non-coding transcript variant"[Molecular Consequence])'
        chr_terms = "(" + " OR ".join(f"{c}[CHR]" for c in chroms) + ")"
        term = f'{cln_query} AND {mc_terms} AND {chr_terms} AND "single nucleotide variant"[Type of Variation]'
        encoded = urllib.parse.quote(term)
        esearch = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={encoded}&retmax=2500&retmode=json"
        req = urllib.request.Request(esearch, headers={"User-Agent": "NC-DriverEngine/3.3"})
        with urllib.request.urlopen(req, timeout=25, context=SSL_CTX) as resp:
            ids = json.loads(resp.read().decode("utf-8")).get("esearchresult", {}).get("idlist", [])
        rng.shuffle(ids)

        for batch_start in range(0, len(ids), 150):
            if len(pool) >= target_per_class + 40:
                break
            batch = ids[batch_start:batch_start + 150]
            esum = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={','.join(batch)}&retmode=json"
            req = urllib.request.Request(esum, headers={"User-Agent": "NC-DriverEngine/3.3"})
            with urllib.request.urlopen(req, timeout=25, context=SSL_CTX) as resp:
                res = json.loads(resp.read().decode("utf-8")).get("result", {})
            for uid in batch:
                rec = res.get(uid, {})
                mcs = [m.lower() for m in rec.get("molecular_consequence_list", [])]
                if not mcs or any(any(ct in m for ct in CODING_TERMS) for m in mcs):
                    continue
                for loc in rec.get("variation_set", [{}])[0].get("variation_loc", []):
                    if loc.get("assembly_name") == "GRCh38" and loc.get("chr") in set(chroms):
                        pos = int(loc.get("start", 0))
                        spdi = rec.get("canonical_spdi", "")
                        spdi_parts = spdi.split(":")
                        if len(spdi_parts) == 4 and len(spdi_parts[2]) == 1 and len(spdi_parts[3]) == 1:
                            ref, alt = spdi_parts[2], spdi_parts[3]
                            if (loc["chr"], pos, ref, alt) not in used_keys:
                                pool.append({
                                    "label": label_name,
                                    "chrom": loc["chr"],
                                    "pos": pos,
                                    "var_id": uid,
                                    "ref": ref,
                                    "alt": alt,
                                    "mc": "|".join(mcs),
                                })
                                used_keys.add((loc["chr"], pos, ref, alt))
                        break
            time.sleep(0.4)
        return pool[:target_per_class]

    import urllib.parse
    path_sample = fetch_class('("clinsig pathogenic"[Properties] OR "clinsig likely pathogenic"[Properties])', "PATHOGENIC")
    ben_sample = fetch_class('("clinsig benign"[Properties] OR "clinsig likely benign"[Properties])', "BENIGN")
    return path_sample, ben_sample


def run_item2_second_benchmark():
    print("\n=== ITEM 2: Second Independent Noncoding Benchmark (chr12-22, n=300/class) ===")
    used_keys = load_existing_keys()
    local_vcf = find_local_clinvar_vcf()
    if local_vcf:
        print(f"Found local ClinVar VCF: {local_vcf}")
        path_sample, ben_sample = sample_from_local_vcf(local_vcf, used_keys, target_per_class=300, seed=104)
    else:
        path_sample, ben_sample = sample_via_ncbi_api(used_keys, target_per_class=300, seed=104)

    cohort = path_sample + ben_sample
    print(f"Constructed independent cohort: {len(path_sample)} pathogenic + {len(ben_sample)} benign = {len(cohort)} variants")

    cache_path = "results/pi_upgrades/independent_n600_flanks.json"
    cache = json.load(open(cache_path, "rt", encoding="utf-8")) if os.path.exists(cache_path) else {}

    tp = fp = tn = fn = 0
    out_rows = []
    for idx, rec in enumerate(cohort, 1):
        seqs = fetch_ucsc_flank(rec["chrom"], rec["pos"], rec["ref"], rec["alt"], cache)
        res = analyze_variant_delta(seqs["wt"], seqs["mut"], flank_bp=FLANK_BP, locus=f"chr{rec['chrom']}:{rec['pos']}")
        pred_pos = res["classification"] in ("Driver Variant", "Loss-of-Function Variant")

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

        out_rows.append({
            "label": rec["label"],
            "chrom": f"chr{rec['chrom']}",
            "pos": rec["pos"],
            "var_id": rec["var_id"],
            "ref": rec["ref"],
            "alt": rec["alt"],
            "mc": rec["mc"],
            "v33_call": res["classification"],
            "delta_pwm": res.get("delta_pwm", 0.0),
        })
        if idx % 100 == 0:
            print(f"  Scored {idx}/{len(cohort)} independent variants...")

    with open(cache_path, "wt", encoding="utf-8") as f:
        json.dump(cache, f)

    with open("results/pi_upgrades/independent_blind_n600_chr12_22.tsv", "wt", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(out_rows)

    sens, s_lo, s_hi = wilson_ci(tp, tp + fn)
    spec, sp_lo, sp_hi = wilson_ci(tn, tn + fp)
    ppv, p_lo, p_hi = wilson_ci(tp, tp + fp)
    acc, a_lo, a_hi = wilson_ci(tp + tn, len(cohort))
    fpr, f_lo, f_hi = wilson_ci(fp, tn + fp)
    f1 = (2 * tp) / (2 * tp + fp + fn) if (2 * tp + fp + fn) else 0.0

    # Realistic 2% prevalence PPV
    prev = 0.02
    tpr_f = tp / max(1, tp + fn)
    fpr_f = fp / max(1, tn + fp)
    ppv_2pct = (tpr_f * prev) / max(1e-12, (tpr_f * prev + fpr_f * (1.0 - prev))) * 100.0

    print("\n=== ITEM 2 RESULTS: Independent Held-Out ClinVar Cohort (chr12-22, N=600) ===")
    print(f"  TP / FP / TN / FN      : {tp} / {fp} / {tn} / {fn}")
    print(f"  Sensitivity (95% CI)   : {sens:.2f}% ({s_lo:.1f}–{s_hi:.1f})")
    print(f"  Specificity (95% CI)   : {spec:.2f}% ({sp_lo:.1f}–{sp_hi:.1f})")
    print(f"  Precision 1:1 (95% CI) : {ppv:.2f}% ({p_lo:.1f}–{p_hi:.1f})")
    print(f"  PPV at 2% prevalence   : {ppv_2pct:.2f}%")
    print(f"  Accuracy (95% CI)      : {acc:.2f}% ({a_lo:.1f}–{a_hi:.1f})")
    print(f"  FPR (95% CI)           : {fpr:.2f}% ({f_lo:.1f}–{f_hi:.1f})")
    print(f"  F1 Score               : {f1:.4f}")


if __name__ == "__main__":
    run_item3_unified_v33()
    run_item2_second_benchmark()
