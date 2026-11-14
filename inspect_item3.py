import inspect, nc_mutscan_heldout_benchmark as hb

print("Sample HELD_OUT_DATASET entries:")
for x in hb.HELD_OUT_DATASET[:3]:
    print(" ", x)

print("\nSearching benchmark scripts for 54-locus set...")
for fname in ("run_complete_benchmark.py", "benchmark_nc_mutscan.py", "verify_full_system.py"):
    try:
        txt = open(fname, encoding="utf-8").read()
        if "54" in txt or "intergenic" in txt.lower() or "wrcy" in txt.lower():
            print(f"  Match found in {fname} (lines={len(txt.splitlines())})")
    except Exception as e:
        print(f"  {fname}: {e}")
