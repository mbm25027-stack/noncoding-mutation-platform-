import re
import time

def parse_fasta(file_path):
    with open(file_path, 'r') as f:
        content = f.read()
    blocks = content.strip().split('>')
    records = []
    for b in blocks:
        if not b.strip():
            continue
        lines = b.strip().split('\n')
        header = lines[0]
        seq = ''.join(lines[1:]).upper().replace(' ', '').replace('\r', '')
        records.append((header, seq))
    return records

def analyze_sequence(header, seq):
    if len(seq) < 6 or re.search(r'[^ACGT]', seq):
        return 'Sanitized / Short'
    
    # AID/SHM WRCY signature detection
    wrcy_count = len(re.findall(r'[AT][AG]C[CT]', seq))
    
    # Oncogenic motif scanning (ETS, MYB, GABP)
    has_driver_motif = bool(re.search(r'CCGGAA|AACGGT|GGGCGG', seq))
    is_target_locus = any(kw in header.lower() for kw in ['promoter', 'enhancer', 'tert', 'plekhs1', 'sdhd', 'tal1'])
    
    if wrcy_count >= 3 and not ('tert' in header.lower() or 'plekhs1' in header.lower()):
        return 'AID/SHM Passenger'
    elif has_driver_motif and is_target_locus:
        return 'Driver Variant'
    else:
        return 'Benign / Neutral'

print("=================================================")
print(" Running Benchmark Evaluation on benchmark_200.fasta")
print("=================================================\n")

try:
    records = parse_fasta('benchmark_200.fasta')
    stats = {
        'total': len(records),
        'driver': 0,
        'passenger': 0,
        'benign': 0,
        'sanitized': 0,
        'errors': 0
    }

    t0 = time.time()

    for header, seq in records:
        try:
            res = analyze_sequence(header, seq)
            if res == 'Driver Variant':
                stats['driver'] += 1
            elif res == 'AID/SHM Passenger':
                stats['passenger'] += 1
            elif res == 'Benign / Neutral':
                stats['benign'] += 1
            else:
                stats['sanitized'] += 1
        except Exception:
            stats['errors'] += 1

    total_time_ms = (time.time() - t0) * 1000

    print(f"Execution Time:     {total_time_ms:.2f} ms")
    print(f"Latency:            {(total_time_ms / stats['total']):.3f} ms / sample")
    print(f"Throughput:         {int((stats['total'] / max(total_time_ms / 1000, 0.001)))} seq/sec\n")
    print("--- Summary Metrics ---")
    print(f"Driver Variants:    {stats['driver']}")
    print(f"AID/SHM Passengers: {stats['passenger']}")
    print(f"Benign / Neutral:   {stats['benign']}")
    print(f"Sanitized / Short:  {stats['sanitized']}")
    print(f"Errors / Failures:  {stats['errors']}")
    print("=================================================")

except FileNotFoundError:
    print("Error: 'benchmark_200.fasta' not found. Ensure fetch_200_fast.py has finished successfully.")