import inspect, nc_mutscan_research_suite as suite

print("=== Functions in nc_mutscan_research_suite.py ===")
for name, obj in inspect.getmembers(suite, inspect.isfunction):
    print(" ", name, inspect.signature(obj))

print("\n=== run_complete_benchmark.py ===")
print(open("run_complete_benchmark.py", encoding="utf-8").read()[:2500])

print("\n=== benchmark_nc_mutscan.py (first 120 lines) ===")
print("\n".join(open("benchmark_nc_mutscan.py", encoding="utf-8").read().splitlines()[:120]))
