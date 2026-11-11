import json
import os

BIN_SIZE = 10000  # 10kb genomic bins
EPI_INDEX = {}    # Key: (chrom, bin_id) -> List of elements
IS_LOADED = False

def build_epigenomic_index(elements_list):
    global EPI_INDEX, IS_LOADED
    EPI_INDEX.clear()
    
    for item in elements_list:
        chrom = item.get("chrom", "")
        start = item.get("start", 0)
        end = item.get("end", 0)
        
        start_bin = start // BIN_SIZE
        end_bin = end // BIN_SIZE
        
        for b in range(start_bin, end_bin + 1):
            key = (chrom, b)
            if key not in EPI_INDEX:
                EPI_INDEX[key] = []
            EPI_INDEX[key].append(item)
            
    IS_LOADED = True
    print(f"[✓] Epigenomic Index indexed into {len(EPI_INDEX):,} genomic bins.")

def query_locus_context(chrom: str, pos: int) -> dict:
    if not IS_LOADED or not chrom or not pos:
        return {"multiplier": 1.0, "encode_type": "None"}
    
    bin_id = pos // BIN_SIZE
    candidate_elements = EPI_INDEX.get((chrom, bin_id), [])
    
    # Fast scan of only 5-10 elements in the candidate bin instead of 1.89M
    for elem in candidate_elements:
        if elem["start"] <= pos <= elem["end"]:
            c_type = elem.get("type", "pCRE")
            mult = 1.5 if c_type in ["PLS", "CTCF"] else 1.2
            return {
                "multiplier": mult,
                "encode_type": c_type,
                "accession": elem.get("accession", "N/A")
            }
            
    return {"multiplier": 1.0, "encode_type": "None"}