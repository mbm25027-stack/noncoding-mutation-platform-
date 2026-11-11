import json
import random
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

SERVER_URL = "http://localhost:8000/api/analyze"

def generate_500_hard_sequences():
    random.seed(42)  # Fixed seed for reproducibility
    dataset = []

    # 1. True Drivers (100 seqs)
    driver_templates = [
        ("TERT_PROMOTER_ETS", "CCCCTTCCGGGAACCCCCAGCGGTC", True),
        ("TERT_ALT_ETS", "AGCCGGAAAGGGGCTGGGAGGAGGG", True),
        ("CTCF_INSULATOR", "ATCGCCACCAGGGGGCGCTAGGCTA", True),
        ("PLEKHS1_HOTSPOT", "GCTAGCCAGGAAGTGCTAGCAGGAA", True),
        ("WDR74_PROMOTER", "CGCGGAAGGCCCCTTTTGGGAAGGC", True)
    ]
    for i in range(100):
        tmpl = driver_templates[i % len(driver_templates)]
        flank_l = "".join(random.choices("ACGT", k=5))
        flank_r = "".join(random.choices("ACGT", k=5))
        seq = flank_l + tmpl[1] + flank_r
        dataset.append({
            "id": f"CAT1_DRIVER_{i+1:03d}",
            "category": "True Positive Drivers",
            "sequence": seq,
            "expected_is_driver": True
        })

    # 2. High GC-Rich Traps (100 seqs)
    for i in range(100):
        seq_len = random.randint(30, 60)
        seq = "".join(random.choices(["G", "C"], weights=[0.5, 0.5], k=seq_len))
        dataset.append({
            "id": f"CAT2_GC_TRAP_{i+1:03d}",
            "category": "High GC Decoys",
            "sequence": seq,
            "expected_is_driver": False
        })

    # 3. Near-Match Mismatches (100 seqs)
    near_matches = [
        "CCCCTTCCGGAGACCCCCAGCGGTC",  # CCGGAA -> CCGGAG (Invalid ETS)
        "AGCCAGAAAGGGGCTGGGAGGAGGG",  # CCGGAA -> CCAGAA
        "ATCGCCACCAGGGGACGCTAGGCTA",  # GGGGG -> GGGAC
        "CGCGGAACGCCCCTTTTGGGAAGGC"   # GGAAGG -> GGAACG
    ]
    for i in range(100):
        base = near_matches[i % len(near_matches)]
        dataset.append({
            "id": f"CAT3_NEAR_MATCH_{i+1:03d}",
            "category": "Near-Match Mismatches",
            "sequence": base,
            "expected_is_driver": False
        })

    # 4. Low-Complexity & Repeats (100 seqs)
    repeat_patterns = ["A", "T", "AC", "GTC", "ATGC"]
    for i in range(100):
        pat = repeat_patterns[i % len(repeat_patterns)]
        seq = (pat * (40 // len(pat)))[:40]
        dataset.append({
            "id": f"CAT4_REPEAT_{i+1:03d}",
            "category": "Low-Complexity Repeats",
            "sequence": seq,
            "expected_is_driver": False
        })

    # 5. Background Noise (100 seqs)
    for i in range(100):
        seq_len = random.randint(30, 80)
        seq = "".join(random.choices("ACGT", weights=[0.3, 0.2, 0.2, 0.3], k=seq_len))
        dataset.append({
            "id": f"CAT5_NOISE_{i+1:03d}",
            "category": "Genomic Background Noise",
            "sequence": seq,
            "expected_is_driver": False
        })

    return dataset

def test_single_sequence(item):
    payload = json.dumps({"sequence": item["sequence"]}).encode("utf-8")
    req = urllib.request.Request(
        SERVER_URL,
        data=payload,
        headers={"Content-Type": "application/json"},
        method="POST"
    )
    
    start_time = time.time()
    try:
        with urllib.request.urlopen(req, timeout=5) as resp:
            elapsed = (time.time() - start_time) * 1000  # ms
            res = json.loads(resp.read().decode("utf-8"))
            pred_driver = res.get("is_driver", False)
            correct = (pred_driver == item["expected_is_driver"])
            return {
                "id": item["id"],
                "category": item["category"],
                "expected": item["expected_is_driver"],
                "predicted": pred_driver,
                "correct": correct,
                "latency_ms": elapsed,
                "status_code": resp.status
            }
    except Exception as e:
        return {
            "id": item["id"],
            "category": item["category"],
            "expected": item["expected_is_driver"],
            "predicted": None,
            "correct": False,
            "error": str(e),
            "latency_ms": 0
        }

def run_benchmark():
    print("==================================================")
    print(" GENERATING 500 NON-BIASED HARD TEST SEQUENCES")
    print("==================================================")
    dataset = generate_500_hard_sequences()
    print(f"[+] Created dataset with {len(dataset)} sequences across 5 categories.")

    print("\n[+] Executing concurrent API benchmark against http://localhost:8000 ...")
    results = []
    start_all = time.time()
    
    # 10 concurrent threads to run stress test quickly
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(test_single_sequence, item) for item in dataset]
        for future in as_completed(futures):
            results.append(future.result())

    total_time = time.time() - start_all

    # Calculate metrics
    tp = sum(1 for r in results if r["expected"] and r["predicted"])
    tn = sum(1 for r in results if not r["expected"] and r["predicted"] is False)
    fp = sum(1 for r in results if not r["expected"] and r["predicted"] is True)
    fn = sum(1 for r in results if r["expected"] and r["predicted"] is False)
    
    total_valid = tp + tn + fp + fn
    acc = (tp + tn) / total_valid if total_valid > 0 else 0
    avg_latency = sum(r["latency_ms"] for r in results) / len(results)

    print("\n==================================================")
    print(" 500-SEQUENCE HARD BENCHMARK RESULTS")
    print("==================================================")
    print(f" Total Executed   : {len(results)} sequences")
    print(f" Total Time       : {total_time:.2f} s ({len(results)/total_time:.1f} req/sec)")
    print(f" Avg Latency      : {avg_latency:.2f} ms / request")
    print(f" Accuracy         : {acc * 100:.2f}%")
    print(f" True Positives   : {tp} / 100")
    print(f" True Negatives   : {tn} / 400")
    print(f" False Positives  : {fp} (Decoy/Noise over-calls)")
    print(f" False Negatives  : {fn} (Missed drivers)")
    print("==================================================")

    # Save detailed output
    output_filename = "benchmark_500_results.json"
    with open(output_filename, "w", encoding="utf-8") as f:
        json.dump({
            "metrics": {
                "total": len(results),
                "accuracy": acc,
                "tp": tp, "tn": tn, "fp": fp, "fn": fn,
                "avg_latency_ms": avg_latency
            },
            "detailed_results": results
        }, f, indent=2)
    print(f"\n[+] Saved complete audit results to '{output_filename}'")

if __name__ == "__main__":
    run_benchmark()