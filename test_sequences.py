import json
import urllib.request

BASE_URL = "http://localhost:8000/api/analyze"

test_cases = [
    {
        "name": "1. TERT Promoter C228T (De Novo ETS Motif)",
        "sequence": "CCCCTTCCGGGAACCCCCAGCGGTC"
    },
    {
        "name": "2. CTCF Insulator Boundary Disruption",
        "sequence": "ATCGCCACCAGGGGGCGCTAGGCTA"
    },
    {
        "name": "3. FOXA1 Enhancer Motif",
        "sequence": "GTAATGTTTACACTGTAAAGGCAGC"
    },
    {
        "name": "4. High-GC Benign Sequence (Negative Control)",
        "sequence": "GCGCGCGCCGCGCGCGGCCCCGCGC"
    },
    {
        "name": "5. Poly-A Noncoding (Negative Control)",
        "sequence": "AAAAAAAAAAAAAAAAAAAAAAAAA"
    }
]

def run_sequence_suite():
    print("==================================================")
    print(" RUNNING MULTI-SCENARIO SEQUENCE MOTIF SUITE")
    print("==================================================")
    
    for case in test_cases:
        payload = json.dumps({"sequence": case["sequence"]}).encode("utf-8")
        req = urllib.request.Request(
            BASE_URL,
            data=payload,
            headers={"Content-Type": "application/json"},
            method="POST"
        )
        
        print(f"\n---> {case['name']}")
        try:
            with urllib.request.urlopen(req) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                status = data.get("status", "Unknown")
                is_driver = data.get("is_driver", False)
                motifs = len(data.get("detected_motifs", []))
                print(f"Status: {status}")
                print(f"Is Driver: {is_driver} | Motifs Detected: {motifs}")
        except Exception as e:
            print(f"Request Error: {e}")

if __name__ == "__main__":
    run_sequence_suite()