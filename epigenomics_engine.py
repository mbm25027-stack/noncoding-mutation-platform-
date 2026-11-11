import os
import bisect

COMBINED_BED_PATH = os.path.join("data", "encode_cCREs_GRCh38_combined.bed")

# Global in-memory index organized by chromosome
GENOME_INDEX = {}
IS_LOADED = False

def load_epigenomic_index():
    """Loads and indexes the sorted BED file into memory for O(log N) binary search."""
    global GENOME_INDEX, IS_LOADED
    if IS_LOADED:
        return

    print("[+] Loading ENCODE Epigenomic Index into memory...")
    if not os.path.exists(COMBINED_BED_PATH):
        print(f"[!] Warning: {COMBINED_BED_PATH} not found. Running without chromatin context.")
        return

    count = 0
    with open(COMBINED_BED_PATH, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, label = parts[0], int(parts[1]), int(parts[2]), parts[3]

            if chrom not in GENOME_INDEX:
                GENOME_INDEX[chrom] = {"starts": [], "ends": [], "labels": []}

            GENOME_INDEX[chrom]["starts"].append(start)
            GENOME_INDEX[chrom]["ends"].append(end)
            GENOME_INDEX[chrom]["labels"].append(label)
            count += 1

    IS_LOADED = True
    print(f"[✔] Epigenomic Index ready: {count:,} regulatory elements loaded across {len(GENOME_INDEX)} chromosomes.")

def query_locus_context(chrom, start, end=None):
    """
    Queries whether a coordinate (or range) falls inside an ENCODE cCRE element.
    Returns: dict with overlap status, element type, and pathogenicity multiplier.
    """
    if not IS_LOADED:
        load_epigenomic_index()

    if end is None:
        end = start + 1

    # Standardize chrom string
    if not chrom.startswith("chr"):
        chrom = f"chr{chrom}"

    if chrom not in GENOME_INDEX:
        return {
            "in_regulatory_region": False,
            "element_type": "Intergenic / Unknown",
            "multiplier": 0.3,
            "description": "Locus outside known ENCODE regulatory elements"
        }

    starts = GENOME_INDEX[chrom]["starts"]
    ends = GENOME_INDEX[chrom]["ends"]
    labels = GENOME_INDEX[chrom]["labels"]

    # Binary search for fast range overlap
    idx = bisect.bisect_right(starts, end) - 1

    matches = []
    while idx >= 0:
        if ends[idx] > start and starts[idx] < end:
            matches.append(labels[idx])
        elif ends[idx] <= start:
            # Past the potential overlap window
            break
        idx -= 1

    if matches:
        element_type = matches[0]
        
        # Calculate impact multiplier based on element class
        multiplier = 1.0
        if "PLS" in element_type:
            desc = "Promoter-like signature (High transcriptional disruption risk)"
            multiplier = 1.2
        elif "ELS" in element_type:
            desc = "Enhancer-like signature (Distal gene expression disruption risk)"
            multiplier = 1.0
        elif "CTCF" in element_type:
            desc = "CTCF Insulator / TAD Boundary Anchor (3D loop breakdown risk)"
            multiplier = 1.15
        else:
            desc = f"Active Regulatory Element ({element_type})"

        return {
            "in_regulatory_region": True,
            "element_type": element_type,
            "multiplier": multiplier,
            "description": desc
        }

    return {
        "in_regulatory_region": False,
        "element_type": "Closed Chromatin / Heterochromatin",
        "multiplier": 0.4,
        "description": "Locus in closed/inactive chromatin region"
    }

# Pre-load index on module import
load_epigenomic_index()