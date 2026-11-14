import hashlib, json, os, socket, subprocess, sys, time, urllib.request

EXPECTED_HASH = "b276d0d44045"
test_port = 18082
env = os.environ.copy()
env["PORT"] = str(test_port)
proc = subprocess.Popen([sys.executable, "server.py"], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

def post_json(url, payload):
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    t0 = time.perf_counter()
    with urllib.request.urlopen(req, timeout=10) as r:
        return r.status, json.loads(r.read().decode("utf-8")), (time.perf_counter() - t0) * 1000.0

try:
    base = f"http://127.0.0.1:{test_port}"
    for _ in range(30):
        try:
            with socket.create_connection(("127.0.0.1", test_port), timeout=0.3):
                break
        except OSError:
            time.sleep(0.15)

    cases = [
        (
            "TC-01 TERT -124 C>T (GABPA/ETS Gain)",
            "AAAAAAAAAAAAAAAAAAAAAAAAACCCCGAGGGAAAAAAAAAAAAAAAAA",
            "AAAAAAAAAAAAAAAAAAAAAAAAACCCGGAAGAAAAAAAAAAAAAAAAAA",
            "Driver Variant",
        ),
        (
            "TC-06 TP53 Promoter Sp1 Disruption",
            "TTTTTTTTTTTTTTTTTTTTTTTTTGCGCGCGGAAAAAAAAAAAAAAAAAA",
            "TTTTTTTTTTTTTTTTTTTTTTTTTGTGTGCGGAAAAAAAAAAAAAAAAAA",
            "Loss-of-Function Variant",
        ),
        (
            "U1 5' Splice Donor Loss (GTAAGT -> ATAAGT)",
            "TGCATGCATGCATGCATGCATCAGGTAAGTTGCATGCATGCATGCATGCAT",
            "TGCATGCATGCATGCATGCATCAGATAAGTTGCATGCATGCATGCATGCAT",
            "Loss-of-Function Variant",
        ),
        (
            "3'-UTR AATAAA PolyA Signal Loss",
            "ATATATATATATATATATATAAATAAAATATATATATATATATATATATAT",
            "ATATATATATATATATATATAAGAAAAATATATATATATATATATATATAT",
            "Loss-of-Function Variant",
        ),
        (
            "TC-03/12 Clustered AID/SHM WRCY/RGYW",
            "AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
            "AGCTAGCTAGCTAGCTAGCTAGCTGGCTAGCTAGCTAGCTAGCTAGCTAGC",
            "AID/SHM Passenger",
        ),
        (
            "TC-04 IUPAC Degenerate Base Sanitization",
            "ATGCATGCATNRYATGCATGCATGCATGCATGCATGCATGCATGCATGCAT",
            "ATGCATGCATNRYATGCATGCATGCGTGCATGCATGCATGCATGCATGCAT",
            "Ambiguous IUPAC Sequence",
        ),
        (
            "TC-10 Short Sequence (<6 bp) Bypass",
            "ACGTA",
            "ACTTA",
            "Sub-Threshold Short Sequence",
        ),
    ]

    all_ok = True
    for label, wt, mut, expected_cls in cases:
        st, res, ms = post_json(f"{base}/predict", {"wt_seq": wt, "mut_seq": mut, "title": label})
        got_cls = res.get("classification")
        ok = (st == 200 and got_cls == expected_cls and res.get("build_sha256") == EXPECTED_HASH)
        all_ok = all_ok and ok
        print(
            f"[{'PASS' if ok else 'FAIL'}] {label:<42} -> {got_cls} | dPWM={res.get('delta_pwm')} | S_raw={res.get('S_raw')} | S_unified={res.get('S_unified')} ({ms:.1f} ms)"
        )

    core_sha = hashlib.sha256(open("nc_mutscan_research_suite.py", "rb").read()).hexdigest()[:12]
    print(f"\nCore SHA-256: {core_sha} | All 7/7 Passed: {all_ok}")
finally:
    proc.terminate()
    proc.wait(timeout=5)