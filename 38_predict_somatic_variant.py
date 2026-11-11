"""
38_predict_somatic_variant.py
-----------------------------
Inference script for the Noncoding Mutation Finding Server.
Evaluates novel variants across Melanoma, Breast, Liver, Colon, Brain, Blood, and Lung.
"""

import json
import math
import sys

def sigmoid(z):
    return 1.0 / (1.0 + math.exp(-z))

def predict_variant(hugo_symbol, cancer_type, vaf, tf_effect, log2fc):
    abs_fc = abs(log2fc)
    is_clonal = 1.0 if vaf > 0.25 else 0.0
    
    # Model scoring function tuned during Pan-Cancer training
    score = (3.5 * vaf) + (1.2 * tf_effect) + (0.8 * abs_fc) + (2.0 * is_clonal) - 2.5
    probability = sigmoid(score)
    
    label = "DRIVER" if probability >= 0.5 else "PASSENGER"
    confidence = round(probability * 100 if label == "DRIVER" else (1 - probability) * 100, 2)
    
    return {
        "Gene": hugo_symbol,
        "Cancer_Type": cancer_type,
        "Variant_Allele_Frequency": vaf,
        "TF_Motif_Disruption": "YES" if tf_effect == 1 else "NO",
        "Expression_Log2FC": log2fc,
        "Prediction": label,
        "Confidence_Score": f"{confidence}%",
        "Driver_Probability": round(probability, 4)
    }

def main():
    # Test cases representing unknown patient samples
    test_samples = [
        {"gene": "TERT", "cancer": "Melanoma", "vaf": 0.60, "tf": 1, "log2fc": 4.3},
        {"gene": "APC", "cancer": "Colon", "vaf": 0.52, "tf": 0, "log2fc": -2.9},
        {"gene": "PAX5", "cancer": "Blood", "vaf": 0.44, "tf": 1, "log2fc": 1.5},
        {"gene": "UNKNOWN_GENE", "cancer": "Lung", "vaf": 0.04, "tf": 0, "log2fc": 0.02}
    ]
    
    print("=== NONCODING MUTATION SERVER INFERENCE TEST ===")
    for sample in test_samples:
        res = predict_variant(sample["gene"], sample["cancer"], sample["vaf"], sample["tf"], sample["log2fc"])
        print(json.dumps(res, indent=2))

if __name__ == "__main__":
    main()