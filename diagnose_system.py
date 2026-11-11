import time
import json
import os
import sys
import urllib.request
import urllib.error

print("==================================================")
print("     NC-MUTSCAN DEEP SYSTEM DIAGNOSTIC AUDIT     ")
print("==================================================")

# 1. Inspect benchmark dataset payload format
print("\n[1] INSPECTING BENCHMARK DATASET STRUCTURE...")
if os.path.exists("benchmark_500_results.json"):
    try:
        with open("benchmark_500_results.json", "r", encoding="utf-8") as f:
            data = json.load(f)
            if isinstance(data, list) and len(data) > 0:
                sample = data[0]
                print(f" [+] Sample Benchmark Record Keys : {list(sample.keys())}")
                payload_key = "request_payload" if "request_payload" in sample else ("payload" if "payload" in sample else None)
                if payload_key:
                    print(f" [+] Benchmark Payload Structure : {sample[payload_key]}")
            elif isinstance(data, dict):
                print(f" [+] Summary Keys : {list(data.keys())}")
    except Exception as e:
        print(f" [!] Error reading benchmark_500_results.json: {e}")
elif os.path.exists("benchmark_500_hard_cases.py"):
    print(" [+] Found benchmark_500_hard_cases.py script.")

# 2. Test Live HTTP Server Endpoint Latency
print("\n[2] TESTING LIVE API ENDPOINT (http://localhost:8000/api/analyze)...")

def test_live_endpoint(name, payload):
    url = "http://localhost:8000/api/analyze"
    req_data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=req_data,
        headers={"Content-Type": "application/json"},
        method="POST"
    )
    t0 = time.perf_counter()
    try:
        with urllib.request.urlopen(req, timeout=10) as resp:
            body = resp.read().decode("utf-8")
            t1 = time.perf_counter()
            ms = (t1 - t0) * 1000
            res_json = json.loads(body)
            print(f"  -> {name}: {ms:.2f} ms")
            print(f"     Status Message : {res_json.get('status', 'N/A')}")
            print(f"     Is Driver      : {res_json.get('is_driver', 'N/A')}")
            return ms, res_json
    except Exception as e:
        print(f"  -> {name} FAILED: {e}")
        return None, None

seq_ms, _ = test_live_endpoint("Sequence Payload Test", {"sequence": "CCCCTTCCGGGAACCCCCAGCGGTC"})
coord_ms, _ = test_live_endpoint("Coordinate Payload Test", {"chrom": "chr5", "pos": 1295228})

# 3. Direct Subsystem Micro-Profiling
print("\n[3] DIRECT SUBSYSTEM MICRO-PROFILING...")

# Subsystem A: Sequence Aligner Engine
try:
    import sequence_aligner_engine
    t0 = time.perf_counter()
    _ = sequence_aligner_engine.align_and_analyze_sequence("CCCCTTCCGGGAACCCCCAGCGGTC")
    t1 = time.perf_counter()
    print(f" [+] sequence_aligner_engine : {((t1-t0)*1000):.4f} ms")
except Exception as e:
    print(f" [!] sequence_aligner_engine profile error: {e}")

# Subsystem B: Epigenomics Engine
try:
    import epigenomics_engine
    if hasattr(epigenomics_engine, "query_locus_context"):
        t0 = time.perf_counter()
        _ = epigenomics_engine.query_locus_context("chr5", 1295228)
        t1 = time.perf_counter()
        print(f" [+] epigenomics_engine.query_locus_context : {((t1-t0)*1000):.4f} ms")
except Exception as e:
    print(f" [!] epigenomics_engine profile error: {e}")

# Subsystem C: Server Coordinate Mutation Engine
try:
    import server
    if hasattr(server, "analyze_mutation"):
        t0 = time.perf_counter()
        _ = server.analyze_mutation("chr5", 1295228)
        t1 = time.perf_counter()
        print(f" [+] server.analyze_mutation : {((t1-t0)*1000):.4f} ms")
except Exception as e:
    print(f" [!] server.analyze_mutation profile error: {e}")

print("\n==================================================")
print("                DIAGNOSTIC VERDICT                ")
print("==================================================")
if seq_ms is not None and coord_ms is not None:
    if seq_ms < 50 and coord_ms > 1000:
        print(" [!] BOTTLENECK IDENTIFIED: Benchmark is testing Coordinate payloads ('chrom'/'pos'),")
        print("     which hit server.analyze_mutation() and an unindexed epigenomics search.")
    elif seq_ms > 1000 and coord_ms > 1000:
        print(" [!] BOTTLENECK IDENTIFIED: Global server middleware or unindexed ENCODE lookups")
        print("     are blocking both sequence and coordinate endpoints.")
    elif seq_ms < 50 and coord_ms < 50:
        print(" [✓] Both endpoints are optimized! Restart benchmark_500_hard_cases.py.")