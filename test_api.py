import json
import urllib.request

BASE_URL = "http://localhost:8000"

def send_request(test_title, endpoint, method="GET", payload=None):
    url = f"{BASE_URL}{endpoint}"
    data = json.dumps(payload).encode("utf-8") if payload else None
    headers = {"Content-Type": "application/json"} if payload else {}
    
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    print(f"\n==================================================")
    print(f" {test_title}")
    print(f"==================================================")
    try:
        with urllib.request.urlopen(req) as resp:
            print(f"[+] {method} {endpoint} -> Status {resp.status}")
            res_json = json.loads(resp.read().decode("utf-8"))
            print(json.dumps(res_json, indent=2))
    except Exception as e:
        print(f"[-] {method} {endpoint} -> Failed: {e}")

if __name__ == "__main__":
    print("Starting NC-MUTSCAN 4-Part Full API Test Harness...")

    # 1. System Health Check Endpoint
    send_request(
        "TEST 1: Health & System Status Check",
        "/api/health",
        method="GET"
    )

    # 2. Structural Variant Analysis Endpoint
    send_request(
        "TEST 2: Structural Variant Engine (/api/analyze_sv)",
        "/api/analyze_sv",
        method="POST",
        payload={
            "sv_type": "TRANSLOCATION",
            "chrom_a": "chr5",
            "pos_a": 1295228,
            "chrom_b": "chr8",
            "pos_b": 128748315
        }
    )

    # 3. Locus Coordinate Analysis Endpoint (Mode B)
    send_request(
        "TEST 3: Coordinate Analysis Engine (/api/analyze - Mode B)",
        "/api/analyze",
        method="POST",
        payload={
            "chrom": "chr5",
            "pos": 1295228,
            "gene": "TERT"
        }
    )

    # 4. Raw Sequence Alignment & Motif Analysis Endpoint (Mode A)
    send_request(
        "TEST 4: Raw Sequence Alignment Engine (/api/analyze - Mode A)",
        "/api/analyze",
        method="POST",
        payload={
            "sequence": "AGGCCGGAAGGGGCTGGGAGGAGGGGCGCG"
        }
    )