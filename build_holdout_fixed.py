#!/usr/bin/env python3
"""
build_holdout_fixed.py - rebuild the noncoding ClinVar pool and balanced hold-out.

Fixes the CODING_MC vocabulary bug: ClinVar's MC field writes stop-gained as
"nonsense" and inframe variants as "inframe_indel", which the old set never
matched. This version excludes both, plus no_sequence_alteration, and (by
default) the canonical splice donor/acceptor classes.

Writes to a NEW folder (default check_blind_v2), never over check_blind/check_dev.

Examples (run from the project folder, one line each):
  python build_holdout_fixed.py --clinvar clinvar.vcf --chroms 12-22 --out check_blind_v2
  python build_holdout_fixed.py --clinvar clinvar.vcf --chroms 1-11,X --out check_dev_v2
  python build_holdout_fixed.py --clinvar clinvar.vcf --chroms 12-22 --splice include --out check_blind_v2_splice

Dependencies: standard library only.
"""
import argparse, json, os, random, sys
from collections import Counter

# ClinVar writes some consequences under different names than Sequence Ontology
# terms, so both vocabularies are listed.
CODING_MC = {
    "missense_variant", "nonsense_variant", "synonymous_variant",
    "frameshift_variant", "inframe_insertion", "inframe_deletion",
    "stop_lost", "stop_retained_variant", "start_lost",
    "initiator_codon_variant", "coding_sequence_variant",
    "protein_altering_variant", "incomplete_terminal_codon_variant",
    "stop_gained",
    "nonsense", "inframe_indel", "missense", "synonymous", "frameshift",
}
NO_CHANGE_MC = {"no_sequence_alteration"}
SPLICE_MC = {"splice_donor_variant", "splice_acceptor_variant"}

STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
}


def parse_chroms(spec):
    out = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, b = part.split("-", 1)
            out.update(str(i) for i in range(int(a), int(b) + 1))
        else:
            out.add(part.upper() if part.isalpha() else part)
    return out


def mc_terms(info):
    # MC looks like "SO:0001583|missense_variant,SO:0001627|intron_variant"
    return {x.split("|")[-1] for x in info.get("MC", "").split(",") if x}


PATH_SIG = {"pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic"}
BENIGN_SIG = {"benign", "likely_benign", "benign/likely_benign"}


def label(info):
    # exact matches only: a substring test wrongly treats
    # "Conflicting_classifications_of_pathogenicity" as pathogenic
    sig = info.get("CLNSIG", "").lower()
    if sig in PATH_SIG:
        return "P"
    if sig in BENIGN_SIG:
        return "B"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinvar", required=True, help="ClinVar VCF (GRCh38)")
    ap.add_argument("--chroms", required=True, help='e.g. "12-22" or "1-11,X"')
    ap.add_argument("--out", default="check_blind_v2")
    ap.add_argument("--min-stars", type=int, default=0)
    ap.add_argument("--max-len", type=int, default=15)
    ap.add_argument("--snv-only", action="store_true")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--n-per-class", type=int, default=75)
    ap.add_argument("--splice", choices=["exclude", "include"], default="exclude",
                    help="exclude (default): drop canonical splice donor/acceptor. "
                         "include: keep them (label the result as noncoding+splice).")
    a = ap.parse_args()

    if os.path.abspath(a.out) in (os.path.abspath("check_blind"), os.path.abspath("check_dev")):
        sys.exit("Refusing to overwrite check_blind / check_dev. Choose a new --out folder.")

    chroms = parse_chroms(a.chroms)
    excluded = set(CODING_MC) | NO_CHANGE_MC
    if a.splice == "exclude":
        excluded |= SPLICE_MC

    header, pool_p, pool_b = [], [], []
    dropped, file_date = Counter(), ""
    with open(a.clinvar, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#"):
                header.append(line)
                if line.startswith("##fileDate="):
                    file_date = line.strip().split("=", 1)[1]
                continue
            if not line.strip():
                continue
            t = line.rstrip("\r\n").split("\t")
            chrom = t[0].replace("chr", "")
            if chrom not in chroms:
                dropped["chrom"] += 1
                continue
            ref, alt, info_s = t[3].upper(), t[4].upper(), t[7]
            info = dict(kv.split("=", 1) for kv in info_s.split(";") if "=" in kv)
            lab = label(info)
            if lab is None:
                dropped["not_P_or_B"] += 1
                continue
            if STARS.get(info.get("CLNREVSTAT", ""), 0) < a.min_stars:
                dropped["review_status"] += 1
                continue
            terms = mc_terms(info)
            if not terms:
                dropped["no_molecular_consequence"] += 1
                continue
            if terms & excluded:
                dropped["excluded_consequence"] += 1
                continue
            if "," in alt or max(len(ref), len(alt)) > a.max_len:
                dropped["too_long_or_multiallelic"] += 1
                continue
            if a.snv_only and not (len(ref) == 1 and len(alt) == 1):
                dropped["not_snv"] += 1
                continue
            (pool_p if lab == "P" else pool_b).append((line if line.endswith("\n") else line + "\n", terms))

    # guard: nothing excluded may remain in the pool
    leaked = Counter()
    for _, terms in pool_p + pool_b:
        for x in terms & excluded:
            leaked[x] += 1
    assert not leaked, f"excluded consequence leaked into pool: {dict(leaked)}"

    os.makedirs(a.out, exist_ok=True)
    meta = [h for h in header if h.startswith("##")]
    chrom_line = [h for h in header if h.startswith("#CHROM")]

    def write_vcf(path, rows):
        with open(path, "w", encoding="utf-8", newline="\n") as o:
            o.writelines(meta + chrom_line)
            o.writelines(r[0] for r in rows)

    write_vcf(os.path.join(a.out, "pool_all.vcf"), pool_p + pool_b)  # pathogenic first

    rng = random.Random(a.seed)
    n = min(a.n_per_class, len(pool_p), len(pool_b))
    hold = rng.sample(pool_p, n) + rng.sample(pool_b, n)
    write_vcf(os.path.join(a.out, "holdout_balanced.vcf"), hold)

    mc_counts = {"P": Counter(), "B": Counter()}
    for lab, pool in (("P", pool_p), ("B", pool_b)):
        for _, terms in pool:
            for x in terms:
                mc_counts[lab][x] += 1

    manifest = {
        "clinvar_file": os.path.basename(a.clinvar), "clinvar_fileDate": file_date,
        "chroms": a.chroms, "min_stars": a.min_stars, "max_len": a.max_len,
        "snv_only": a.snv_only, "seed": a.seed, "n_per_class": n,
        "splice_scope": a.splice, "excluded_mc_terms": sorted(excluded),
        "pool_counts_after_filters": {"pathogenic_pool": len(pool_p), "benign_pool": len(pool_b)},
        "dropped_by_reason": dict(dropped),
        "mc_counts_in_pool": {k: dict(v.most_common()) for k, v in mc_counts.items()},
    }
    with open(os.path.join(a.out, "manifest.json"), "w") as o:
        json.dump(manifest, o, indent=2)

    print(f"Pool ({a.chroms}, splice={a.splice}, min_stars={a.min_stars}): "
          f"pathogenic {len(pool_p):,} | benign {len(pool_b):,}")
    print("Dropped:", dict(dropped))
    print(f"Wrote {a.out}\\pool_all.vcf, holdout_balanced.vcf ({n} P + {n} B), manifest.json")
    print("Check: for chr12-22 with defaults, expect about 1,847 P / 193,519 B "
          "(your clean pool from the evaluation).")


if __name__ == "__main__":
    main()