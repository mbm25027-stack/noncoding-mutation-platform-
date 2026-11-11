import json
import urllib.request

def test_structural_variant_api():
    url = "http://localhost:8000/api/analyze_sv"
    
    payload = {
        "sv_type": "TRANSLOCATION",
        "chrom_a": "chr5",
        "pos_a": 1295228,
        "chrom_b": "chr8",
        "pos_b": 128748315
    }
    
    headers = {"Content-Type": "application/json"}
    json_bytes = json.dumps(payload).encode("utf-8")
    
    req = urllib.request.Request(
        url=url,
        data=json_bytes,
        headers=headers,
        method="POST"
    )
    
    try:
        with urllib.request.urlopen(req) as response:
            res_body = response.read().decode("utf-8")
            print("==================================================")
            print(" SV ENGINE RESPONSE RECEIVED")
            print("==================================================")
            print(res_body)
    except Exception as err:
        print(f"API Request Error: {err}")

if __name__ == "__main__":
    test_structural_variant_api()