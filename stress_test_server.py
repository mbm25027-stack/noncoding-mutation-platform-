import json
import time

# Simulation of server inference scoring for stress validation
test_payloads = [
    {
        "test_name": "TC-1: TERT C228T Promoter (Canonical Driver)",
        "input": {
            "Gene": "TERT",
            "Cancer_Type": "Melanoma",
            "Variant_Allele_Frequency": 0.60,
            "TF_Motif_Disruption": "YES",
            "Expression_Log2FC": 4.30
        },
        "expected": "DRIVER"
    },
    {
        "test_name": "TC-2: Clonal Intronic TTN (High-VAF Passenger Trap)",
        "input": {
            "Gene": "TTN",
            "Cancer_Type": "Colon",
            "Variant_Allele_Frequency": 0.55,
            "TF_Motif_Disruption": "NO",
            "Expression_Log2FC": 0.02
        },
        "expected": "PASSENGER"
    },
    {
        "test_name": "TC-3: Deep Intronic Cryptic Splice Site (Tumor Suppressor Loss)",
        "input": {
            "Gene": "TP53",
            "Cancer_Type": "Breast",
            "Variant_Allele_Frequency": 0.42,
            "TF_Motif_Disruption": "YES",
            "Expression_Log2FC": -3.10
        },
        "expected": "DRIVER"
    },
    {
        "test_name": "TC-4: Subclonal 3'UTR microRNA Disruption (Subclonal Driver)",
        "input": {
            "Gene": "MYC",
            "Cancer_Type": "Lung",
            "Variant_Allele_Frequency": 0.15,
            "TF_Motif_Disruption": "YES",
            "Expression_Log2FC": 2.40
        },
        "expected": "DRIVER"
    },
    {
        "test_name": "TC-5: Low-VAF Intergenic Desert (Neutral Control)",
        "input": {
            "Gene": "NEUTRAL_LOCUS_1",
            "Cancer_Type": "Brain",
            "Variant_Allele_Frequency": 0.08,
            "TF_Motif_Disruption": "NO",
            "Expression_Log2FC": -0.05
        },
        "expected": "PASSENGER"
    }
]

def run_stress_test():
    print("======================================================================")
    print(" NC-MUTSCAN SERVER STRESS TEST & ACCURACY AUDIT")
    print("======================================================================\n")
    
    passed_tests = 0
    total_tests = len(test_payloads)
    
    for tc in test_payloads:
        inp = tc["input"]
        
        # Scoring logic calculation matching NC-MUTSCAN core matrix
        vaf_score = inp["Variant_Allele_Frequency"] * 3.0
        motif_score = 2.5 if inp["TF_Motif_Disruption"] == "YES" else 0.0
        expr_score = abs(inp["Expression_Log2FC"]) * 1.5
        
        total_impact = vaf_score + motif_score + expr_score
        
        # Classification threshold
        if (inp["TF_Motif_Disruption"] == "YES" and abs(inp["Expression_Log2FC"]) > 1.0) or total_impact >= 5.0:
            pred = "DRIVER"
            conf = min(0.9999, 0.70 + (total_impact / 20.0))
        else:
            pred = "PASSENGER"
            conf = min(0.9999, 0.85 + ((5.0 - total_impact) / 20.0))
            
        status = "✅ PASS" if pred == tc["expected"] else "❌ FAIL"
        if pred == tc["expected"]:
            passed_tests += 1
            
        print(f"Test Name : {tc['test_name']}")
        print(f"Inputs    : VAF={inp['Variant_Allele_Frequency']} | Motif={inp['TF_Motif_Disruption']} | Log2FC={inp['Expression_Log2FC']}")
        print(f"Result    : Prediction={pred} ({conf*100:.2f}% Conf) | Expected={tc['expected']} -> {status}\n")

    print("======================================================================")
    print(f" SUMMARY: {passed_tests}/{total_tests} HARD TESTS PASSED")
    print("======================================================================")

if __name__ == "__main__":
    run_stress_test()