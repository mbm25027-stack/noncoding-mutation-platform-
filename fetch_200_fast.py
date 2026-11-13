import urllib.request
import json
import time
import random
import concurrent.futures

hotspots = [
    ("chr5", 1295113, "TERT_promoter", "C>T"),
    ("chr10", 113824420, "PLEKHS1_promoter", "C>T"),
    ("chr11", 111957575, "SDHD_promoter", "C>T"),
    ("chr1", 47210120, "TAL1_enhancer", "A>G"),
    ("chr18", 63318000, "BCL2_SHM_Region", "C>T"),
]

loci = list(hotspots)
random.seed(42)
for i in range(195):
    chrom = f"chr{random.randint(1, 22)}"
    pos = random.randint(1000000, 50000000)
    loci.append((chrom, pos, f"Intergenic_Locus_{i+1}", ""))

print("Fetching 200 real hg38 sequence snippets in parallel...")

def fetch_seq(item):
    chrom, pos, label, var = item
    chr_clean = chrom.replace("chr", "")
    start = max(1, pos - 25)
    end = pos + 25
    url = f"https://rest.ensembl.org/sequence/region/human/{chr_clean}:{start}..{end}?content-type=application/json"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json", "User-Agent": "NC-MutScan/3.1"})
    try:
        with urllib.request.urlopen(req, timeout=5) as resp:
            seq = json.loads(resp.read().decode()).get("seq", "")
            if seq:
                var_str = f"|{var}" if var else ""
                return f">{chrom}:{pos}|{label}{var_str}\n{seq}"
    except Exception:
        return None

t0 = time.time()
with concurrent.futures.ThreadPoolExecutor(max_workers=25) as executor:
    results = [r for r in executor.map(fetch_seq, loci) if r]

with open("benchmark_200.fasta", "w") as f:
    f.write("\n".join(results))

print(f"Done in {time.time() - t0:.2f}s! Saved {len(results)} sequences to 'benchmark_200.fasta'")