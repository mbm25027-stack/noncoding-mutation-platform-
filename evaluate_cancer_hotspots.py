import json
import urllib.request

SERVER_URL = "http://localhost:8000/api/analyze"

# Cancer Driver Hotspots & Clean Benign Controls
TEST_VARIANTS = [
    # 1. Pathogenic Driver Hotspots
    {
        "gene": "TERT_promoter_-124C>T",
        "expected": True,
        "payload": {
            "chrom": "chr5",
            "pos": 1295113,
            "ref": "C",
            "alt": "T",
            "sequence": ">chr5:1295113|TERT_-124C>T\nCCCACGTGGCGGCCCAGCCCCCTCCGGGCCCTCCCAGCCCCTCCCCTTCCT"
        }
    },
    {
        "gene": "TERT_promoter_-146C>T",
        "expected": True,
        "payload": {
            "chrom": "chr5",
            "pos": 1295135,
            "ref": "C",
            "alt": "T",
            "sequence": ">chr5:1295135|TERT_-146C>T\nGGGGCTGGGCCGGGGACCCGGGAGGGGTCGGGACGGGGCGGGGTCCCGGGG"
        }
    },
    {
        "gene": "PLEKHS1_promoter",
        "expected": True,
        "payload": {
            "chrom": "chr10",
            "pos": 113824420,
            "ref": "C",
            "alt": "T",
            "sequence": ">chr10:113824420|PLEKHS1\nGGAGGCCGCGGCCGCCGGAACTCGGGCCCGCCGCCCAGCCCCGCCCCCGCC"
        }
    },
    {
        "gene": "SDHD_promoter",
        "expected": True,
        "payload": {
            "chrom": "chr11",
            "pos": 111957575,
            "ref": "C",
            "alt": "T",
            "sequence": ">chr11:111957575|SDHD\nCCCGCCCCCGCCGCCGGAACTCGGCCGCCGCCCCCGCCGCCCGCCCCCGCC"
        }
    },

    # 2. Benign Control Loci (Clean non-repetitive sequences)
    {
        "gene": "Intergenic_Benign_1",
        "expected": False,
        "payload": {
            "chrom": "chr12",
            "pos": 56100125,
            "ref": "A",
            "alt": "G",
            "sequence": ">chr12:56100125|Benign1\nACGTGACTACTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
        }
    },
    {
        "gene": "Intergenic_Benign_2",
        "expected": False,
        "payload": {
            "chrom": "chr1",
            "pos": 1000000,
            "ref": "A",
            "alt": "C",
            "sequence": ">chr1:1000000|Benign2\nTTAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        }
    }
]

def run_verification():
    print("==========================================================================")
    print("     NC-MutScan CANCER DRIVER HOTSPOT VERIFICATION (CORRECTED)")
    print("==========================================================================\n")

    correct = 0
    for item in TEST_VARIANTS:
        payload_data = json.dumps(item["payload"]).encode("utf-8")
        req = urllib.request.Request(
            SERVER_URL, 
            data=payload_data, 
            headers={"Content-Type": "application/json"}
        )
        
        try:
            with urllib.request.urlopen(req) as resp:
                res = json.loads(resp.read().decode("utf-8"))
                
                is_driver = res.get("is_driver", False)
                score = res.get("score", 0.0)
                
                match = (is_driver == item["expected"])
                if match:
                    correct += 1
                    
                flag = "[✔] PASS" if match else "[✘] FAIL"
                print(f"{flag} | Locus: {item['gene']:<22} -> Predicted: {str(is_driver):<5} | Score: {score}")
                
        except Exception as e:
            print(f"[!] Request error for {item['gene']}: {e}")

    acc = (correct / len(TEST_VARIANTS)) * 100
    print("\n==========================================================================")
    print(f"Targeted Hotspot Accuracy: {acc:.1f}% ({correct}/{len(TEST_VARIANTS)})")
    print("==========================================================================")

if __name__ == "__main__":
    run_verification()