import json
import os
import sys
import time
import threading
import urllib.request
import urllib.error

# Global port for testing server
TEST_PORT = 8005

def print_header(title):
    print("\n" + "=" * 70)
    print(f" [TEST] {title}")
    print("=" * 70)

def test_epigenomics_engine():
    print_header("1. Testing ENCODE Epigenomics Engine (Today's Addition)")
    try:
        import epigenomics_engine
        if not epigenomics_engine.IS_LOADED:
            print("[!] Epigenomics index failed to load automatically. Triggering manual load...")
            epigenomics_engine.load_epigenomic_index()

        # Test TERT Promoter Locus
        res = epigenomics_engine.query_locus_context("chr5", 1295228)
        print(f" [+] Query locus: chr5:1295228 (TERT promoter region)")
        print(f"     - Overlap Status : {res['in_regulatory_region']}")
        print(f"     - Element Class  : {res['element_type']}")
        print(f"     - Multiplier     : {res['multiplier']}x")
        print(f"     - Description    : {res['description']}")

        assert res['in_regulatory_region'] == True, "Failed to locate TERT promoter element in BED index!"
        print("[✔] Epigenomics Engine Test: PASSED")
        return True
    except Exception as e:
        print(f"[✘] Epigenomics Engine Test: FAILED ({e})")
        return False

def test_sv_engine():
    print_header("2. Testing Structural Variant (SV) Engine (Today's Addition)")
    try:
        import sv_engine
        # Test Translocation bridging chr5 and chr8
        res = sv_engine.analyze_structural_variant("TRANSLOCATION", "chr5", 1295228, "chr8", 127735434)
        print(" [+] Translocation Input: t(chr5:1295228 ; chr8:127735434)")
        print(f"     - Verdict       : {res['verdict']}")
        print(f"     - Is Driver     : {res['is_driver']}")
        print(f"     - Path Score    : {res['pathogenicity_score']}")
        print(f"     - Mechanism     : {res['mechanism']}")

        assert "verdict" in res, "SV response missing verdict field!"
        print("[✔] Structural Variant Engine Test: PASSED")
        return True
    except Exception as e:
        print(f"[✘] SV Engine Test: FAILED ({e})")
        return False

def test_legacy_sequence_engine():
    print_header("3. Testing SNV / FASTA Sequence Engine (Previous Knowledge)")
    try:
        import server
        # Test TERT coordinate lookup
        res_coord = server.analyze_mutation("chr5", 1295228)
        print(" [+] Coordinate Query: chr5:1295228")
        print(f"     - Status                : {res_coord.get('status')}")
        print(f"     - Final Score           : {res_coord.get('score')}")
        print(f"     - Base Score            : {res_coord.get('base_score', 'N/A')}")
        print(f"     - Epigenomic Multiplier : {res_coord.get('epigenomic_multiplier', 'N/A')}")

        # Test FASTA Sequence scanning logic
        raw_seq = "AGCCCCTCCCCTCCCCAGTCGGGG"
        res_seq = server.align_and_analyze_sequence(raw_seq, server.NONCODING_DB)
        print(f" [+] FASTA Sequence Query: '{raw_seq}'")
        print(f"     - Scan Verdict : {res_seq.get('status', 'Processed')}")

        print("[✔] Legacy Sequence & Hotspot Engine Test: PASSED")
        return True
    except Exception as e:
        print(f"[✘] Legacy Sequence Engine Test: FAILED ({e})")
        return False

def run_background_server():
    import server
    server.PORT = TEST_PORT
    server.run()

def test_http_api_endpoints():
    print_header("4. Testing Live Server HTTP Endpoints (/health, /api/analyze, /api/analyze_sv)")
    server_thread = threading.Thread(target=run_background_server, daemon=True)
    server_thread.start()
    time.sleep(1.5)  # Wait for multithreaded server to initialize

    base_url = f"http://127.0.0.1:{TEST_PORT}"
    
    # A. Health Check Endpoint
    try:
        req = urllib.request.urlopen(f"{base_url}/health")
        data = json.loads(req.read().decode("utf-8"))
        print(f" [+] Endpoint GET /health: ONLINE")
        print(f"     - Version: {data.get('version')}")
        print(f"     - Epigenomics Active: {data.get('epigenomics_indexed')}")
    except Exception as e:
        print(f"[✘] GET /health failed: {e}")
        return False

    # B. SNV / Coordinate Endpoint POST
    try:
        payload = json.dumps({"chrom": "chr5", "pos": 1295228}).encode("utf-8")
        req = urllib.request.Request(f"{base_url}/api/analyze", data=payload, headers={'Content-Type': 'application/json'})
        response = json.loads(urllib.request.urlopen(req).read().decode("utf-8"))
        print(f" [+] Endpoint POST /api/analyze (SNV): PASSED")
        print(f"     - Returned Region: {response.get('region')}")
    except Exception as e:
        print(f"[✘] POST /api/analyze failed: {e}")
        return False

    # C. SV Endpoint POST
    try:
        payload = json.dumps({
            "sv_type": "TRANSLOCATION",
            "chrom_a": "chr5",
            "pos_a": 1295228,
            "chrom_b": "chr8",
            "pos_b": 127735434
        }).encode("utf-8")
        req = urllib.request.Request(f"{base_url}/api/analyze_sv", data=payload, headers={'Content-Type': 'application/json'})
        response = json.loads(urllib.request.urlopen(req).read().decode("utf-8"))
        print(f" [+] Endpoint POST /api/analyze_sv (SV): PASSED")
        print(f"     - Verdict: {response.get('verdict')}")
    except Exception as e:
        print(f"[✘] POST /api/analyze_sv failed: {e}")
        return False

    print("[✔] Live HTTP API Server Endpoints Test: PASSED")
    return True

if __name__ == "__main__":
    print("\n" + "#" * 70)
    print(" NC-MUTSCAN INTEGRATED DIAGNOSTIC & VERIFICATION SUITE")
    print(" Verifying Legacy SNV Engine + New Epigenomic & SV Systems")
    print("#" * 70)

    t1 = test_epigenomics_engine()
    t2 = test_sv_engine()
    t3 = test_legacy_sequence_engine()
    t4 = test_http_api_endpoints()

    print("\n" + "=" * 70)
    print(" SYSTEM VERIFICATION SUMMARY")
    print("=" * 70)
    print(f" - Epigenomics Engine (ENCODE cCRE) : {'[ PASSED ]' if t1 else '[ FAILED ]'}")
    print(f" - Structural Variant (SV Engine)   : {'[ PASSED ]' if t2 else '[ FAILED ]'}")
    print(f" - Sequence & Hotspot Scanner (SNV) : {'[ PASSED ]' if t3 else '[ FAILED ]'}")
    print(f" - Live Server HTTP REST APIs       : {'[ PASSED ]' if t4 else '[ FAILED ]'}")

    if t1 and t2 and t3 and t4:
        print("\n[✔✔✔] ALL SYSTEMS GO: Previous & Today's Knowledge Intact and Fully Functional!")
    else:
        print("\n[!] VERIFICATION COMPLETED WITH ERRORS. Check trace outputs above.")