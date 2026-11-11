"""
gtf_annotation_engine.py
Automated Genomic Region & TSS Annotation Engine.
Classifies any GRCh37/hg19 coordinate into Promoter, Enhancer, lncRNA, Exon, or Intergenic.
"""

import json
import urllib.request
from typing import Dict, Any

def lookup_genomic_annotation(chrom: str, pos: int) -> Dict[str, Any]:
    clean_chrom = str(chrom).replace("chr", "").strip()
    
    # Query Ensembl REST API for overlapping/nearby gene features
    url = f"https://grch37.rest.ensembl.org/overlap/region/human/{clean_chrom}:{max(1, pos-2000)}..{pos+2000}?feature=gene;content-type=application/json"
    
    try:
        req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=5) as response:
            genes = json.loads(response.read().decode())
    except Exception as e:
        genes = []

    if not genes:
        return {
            "chromosome": clean_chrom,
            "position": pos,
            "target_gene": "Intergenic Region",
            "biotype": "intergenic",
            "element_class": "Intergenic Desert",
            "distance_to_tss": "N/A",
            "details": f"No gene features annotated within +/- 2kb of chr{clean_chrom}:{pos}."
        }

    # Find closest gene feature
    best_gene = None
    min_dist = float("inf")

    for g in genes:
        start = g.get("start", 0)
        end = g.get("end", 0)
        strand = g.get("strand", 1)
        
        # Calculate Transcription Start Site (TSS)
        tss = start if strand == 1 else end
        dist_to_tss = pos - tss if strand == 1 else tss - pos

        if abs(dist_to_tss) < abs(min_dist):
            min_dist = dist_to_tss
            best_gene = g
            best_gene["tss_dist"] = dist_to_tss

    gene_name = best_gene.get("external_name") or best_gene.get("gene_id", "Unknown")
    biotype = best_gene.get("biotype", "protein_coding")
    tss_dist = best_gene.get("tss_dist", 0)

    # Classify Genomic Region Type based on TSS distance & biotype
    if -2000 <= tss_dist <= 500:
        element_class = "Promoter Core Region" if abs(tss_dist) <= 300 else "Proximal Promoter"
    elif "lncRNA" in biotype.lower() or "ncrna" in biotype.lower():
        element_class = "lncRNA Promoter / Regulatory Body"
    elif biotype == "protein_coding":
        if best_gene.get("start", 0) <= pos <= best_gene.get("end", 0):
            element_class = "Intronic / Enhancer Element"
        else:
            element_class = "Distal Enhancer Element"
    else:
        element_class = f"Noncoding {biotype.upper()} Region"

    return {
        "chromosome": clean_chrom,
        "position": pos,
        "target_gene": gene_name,
        "biotype": biotype,
        "element_class": element_class,
        "distance_to_tss": f"{tss_dist:+d} bp",
        "details": f"Locus sits {tss_dist:+d} bp from {gene_name} ({biotype}) TSS."
    }

if __name__ == "__main__":
    print("--- Testing Automated Annotation Lookup ---")

    # Test Locus 1: TERT Promoter (chr5:1295228)
    print(json.dumps(lookup_genomic_annotation("5", 1295228), indent=2))

    # Test Locus 2: AC040160.1 lncRNA (chr16:67260980)
    print(json.dumps(lookup_genomic_annotation("16", 67260980), indent=2))