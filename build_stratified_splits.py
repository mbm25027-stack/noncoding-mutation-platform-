import random

INPUT_VCF = "clinvar.vcf"
DEV_OUT = "dev_chr1_11.vcf"
BLIND_OUT = "locked_blind_chr12_22.vcf"

DEV_CHROMS = {str(i) for i in range(1, 12)}
BLIND_CHROMS = {str(i) for i in range(12, 23)}

NC_TERMS = (
    "upstream_gene_variant",
    "5_prime_utr_variant",
    "3_prime_utr_variant",
    "regulatory_region_variant",
    "non-coding_transcript_variant",
    "non_coding_transcript_variant",
    "intron_variant",
)

CODING_EXCLUSIONS = (
    "frameshift",
    "nonsense",
    "missense",
    "inframe",
    "stop_gained",
    "stop_lost",
    "start_lost",
    "synonymous",
    "protein_altering",
)

def parse_info(info_str: str) -> dict:
    fields = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            fields[k] = v
    return fields

headers = []
buckets = {"dev_path": [], "dev_benign": [], "blind_path": [], "blind_benign": []}

with open(INPUT_VCF, "rt", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"):
            headers.append(line)
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue

        chrom = parts[0].replace("chr", "").strip()
        ref, alt, info_str = parts[3], parts[4], parts[7]
        if len(ref) > 15 or len(alt) > 15:
            continue

        info = parse_info(info_str)
        clnsig = info.get("CLNSIG", "").lower()
        mc = info.get("MC", "").lower()

        # Require non-coding annotation AND strictly exclude any coding transcript overlap
        if not any(t in mc for t in NC_TERMS):
            continue
        if any(bad in mc for bad in CODING_EXCLUSIONS):
            continue

        is_path = clnsig in ("pathogenic", "likely_pathogenic", "pathogenic/likely_pathogenic")
        is_benign = clnsig in ("benign", "likely_benign", "benign/likely_benign")

        if chrom in DEV_CHROMS:
            if is_path: buckets["dev_path"].append(line)
            elif is_benign: buckets["dev_benign"].append(line)
        elif chrom in BLIND_CHROMS:
            if is_path: buckets["blind_path"].append(line)
            elif is_benign: buckets["blind_benign"].append(line)

random.seed(42)
MAX_PER_CLASS = 75

def write_split(filename, path_list, benign_list):
    sampled_p = random.sample(path_list, min(MAX_PER_CLASS, len(path_list)))
    sampled_b = random.sample(benign_list, min(MAX_PER_CLASS, len(benign_list)))
    with open(filename, "w", encoding="utf-8") as out:
        out.writelines(headers)
        for row in sampled_p:
            out.write("PATHOGENIC\t" + row)
        for row in sampled_b:
            out.write("BENIGN\t" + row)
    print(f"Saved {filename}: {len(sampled_p)} Pathogenic + {len(sampled_b)} Benign (Pool: {len(path_list)}P / {len(benign_list)}B)")

write_split(DEV_OUT, buckets["dev_path"], buckets["dev_benign"])
write_split(BLIND_OUT, buckets["blind_path"], buckets["blind_benign"])
