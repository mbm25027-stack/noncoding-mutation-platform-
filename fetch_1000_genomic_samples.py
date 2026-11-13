import urllib.request
import json
import time
import random

# Core benchmark driver hotspots with reference and mutant alleles
hotspots = [
    ("chr5", 1295113, "TERT_promoter", "C>T"),
    ("chr10", 113824420, "PLEKHS1_promoter", "C>T"),
    ("chr11", 111957575, "SDHD_promoter", "C>T"),
    ("chr1", 47210120, "TAL1_enhancer", "A>G"),
    ("chr18", 63318000, "BCL2_SHM_Region", "C>T"),
]

loci = list(hotspots)

# Generate remaining genomic positions across chromosomes 1-22
random.seed(42)  # Fixed seed for reproducibility
for i in range(995):
    chrom = f"chr{random.randint(1, 22)}"
    pos = random.randint(1000000, 50000000)
    loci.append((chrom, pos, f"Intergenic_Locus_{i+1}", ""))

print(f"Fetching 1,000 real hg38 sequence snippets...")

def fetch_seq(chrom, pos, window=50):
    chr_clean = chrom.replace("chr", "")
    start = max(1, pos - window // 2)
    end = pos + window // 2
    url = f"https://rest.ensembl.org/sequence/region/human/{chr_clean}:{start}..{end}?content-type=application/json"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req) as resp:
            return json.loads(resp.read().decode()).get("seq", "")
    except Exception:
        return None

fasta_entries = []
for idx, (chrom, pos, label, var) in enumerate(loci):
    seq = fetch_seq(chrom, pos)
    if seq:
        var_str = f"|{var}" if var else ""
        header = f">{chrom}:{pos}|{label}{var_str}"
        fasta_entries.append(f"{header}\n{seq}")
    if (idx + 1) % 100 == 0 or idx == len(loci) - 1:
        print(f"Progress: {idx + 1} / 1000 sequences fetched")
    time.sleep(0.02)  # Respect API rate limits

with open("benchmark_1000.fasta", "w") as f:
    f.write("\n".join(fasta_entries))

print(f"\nSuccess: Saved {len(fasta_entries)} FASTA sequences to 'benchmark_1000.fasta'")