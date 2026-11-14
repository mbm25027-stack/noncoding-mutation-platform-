import gzip
import json
import urllib.request
import urllib.error
import time

VCF_FILE = "clinvar_GRCh38.vcf.gz"
SERVER_URL = "http://localhost:8000/api/analyze"
MAX_SAMPLES = 1000  # 500 Pathogenic + 500 Benign real loci

def parse_clinvar_noncoding(vcf_path, max_items=1000):
    """Parses real ClinVar VCF for noncoding pathogenic and benign variants."""
    print(f"[+] Parsing real ClinVar variants from {vcf_path}...")
    
    pathogenic = []
    benign = []
    
    # Common noncoding consequence terms in ClinVar
    noncoding_terms = [
        "intron_variant", "upstream_gene_variant", "downstream_gene_variant",
        "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_transcript_variant",
        "intergenic_variant"
    ]

    with gzip.open(vcf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            
            chrom = parts[0]
            pos = parts[1]
            info = parts[7]

            # Determine Clinical Significance
            clnsig = ""
            for item in info.split(";"):
                if item.startswith("CLNSIG="):
                    clnsig = item.split("=")[1]
                    break
            
            # Check if variant is noncoding
            is_noncoding = any(term in info for term in noncoding_terms)
            if not is_noncoding:
                continue

            record = {
                "chrom": chrom if chrom.startswith("chr") else f"chr{chrom}",
                "pos": pos,
                "ref": parts[3],
                "alt": parts[4]
            }

            if ("Pathogenic" in clnsig or "Likely_pathogenic" in clnsig) and len(pathogenic) < max_items // 2:
                record["ground_truth"] = 1
                pathogenic.append(record)
            elif ("Benign" in clnsig or "Likely_benign" in clnsig) and len(benign) < max_items // 2:
                record["ground_truth"] = 0
                benign.append(record)

            if len(pathogenic) >= max_items // 2 and len(benign) >= max_items // 2:
                break

    combined = pathogenic + benign
    print(f"[✔] Extracted {len(pathogenic)} Pathogenic & {len(benign)} Benign real noncoding loci.")
    return combined


def query_nc_mutscan_server(chrom, pos):
    """Sends POST request to live server.py at http://localhost:8000/api/analyze."""
    payload = json.dumps({"chrom": chrom, "pos": pos}).encode("utf-8")
    req = urllib.request.Request(
        SERVER_URL, 
        data=payload, 
        headers={"Content-Type": "application/json"}
    )
    
    try:
        with urllib.request.urlopen(req) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            return data
    except urllib.error.URLError as e:
        print(f"[!] Server connection error: {e}. Is server.py running on http://localhost:8000?")
        return None


def run_empirical_benchmark():
    dataset = parse_clinvar_noncoding(VCF_FILE, max_items=MAX_SAMPLES)
    if not dataset:
        print("[!] No noncoding variants extracted. Check VCF file.")
        return

    print("\n[+] Querying NC-MutScan server for empirical validation...")
    
    y_true = []
    y_pred = []
    scores = []
    
    start_time = time.time()
    
    for i, item in enumerate(dataset):
        res = query_nc_mutscan_server(item["chrom"], item["pos"])
        if not res or "error" in res:
            continue
        
        gt = item["ground_truth"]
        is_driver = 1 if res.get("is_driver", False) else 0
        score = res.get("score", 0.0)
        
        y_true.append(gt)
        y_pred.append(is_driver)
        scores.append(score)
        
        if (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{len(dataset)} loci...")

    elapsed = time.time() - start_time
    
    # Calculate Empirical Metrics
    tp = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 1)
    tn = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 0)
    fp = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 0 and yp == 1)
    fn = sum(1 for yt, yp in zip(y_true, y_pred) if yt == 1 and yp == 0)

    total = len(y_true)
    acc = (tp + tn) / total if total > 0 else 0
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0
    f1 = (2 * prec * sens) / (prec + sens) if (prec + sens) > 0 else 0

    print("\n==========================================================================")
    print("      REAL CLINVAR EMPIRICAL BENCHMARK RESULTS (NC-MutScan v3.1)")
    print("==========================================================================")
    print(f" Total Real ClinVar Loci Tested : {total}")
    print(f" Execution Time                 : {elapsed:.2f} seconds ({total/elapsed:.1f} req/sec)")
    print(f" Empirical Accuracy             : {acc * 100:.2f}%")
    print(f" Empirical Sensitivity (Recall) : {sens * 100:.2f}%")
    print(f" Empirical Specificity           : {spec * 100:.2f}%")
    print(f" Empirical Precision             : {prec * 100:.2f}%")
    print(f" Empirical F1-Score             : {f1:.4f}")
    print(f" True Positives (TP): {tp} | False Positives (FP): {fp}")
    print(f" True Negatives (TN): {tn} | False Negatives (FN): {fn}")
    print("==========================================================================\n")


if __name__ == "__main__":
    run_empirical_benchmark()