import os

os.makedirs("results/clean_submission", exist_ok=True)

with open("results/clean_submission/Table2_Blind_Ablations_SpliceAI.csv", "w", encoding="utf-8") as f:
    f.write("Model,Scored_N,TP,FP,TN,FN,Sensitivity_Pct,Sens_95CI,Specificity_Pct,Spec_95CI,Balanced_Precision_Pct,F1_Score\n")
    f.write("NC-MutScan v3.3 (Full Differential Delta-PWM),150,37,1,74,38,49.33,38.3-60.4,98.67,92.8-99.8,97.37,0.6549\n")
    f.write("Ablation: Splice/PolyA-PWM Only,150,36,0,75,39,48.00,37.1-59.1,100.00,95.1-100.0,100.00,0.6486\n")
    f.write("Ablation: TF-Motif-PWM Only,150,1,1,74,74,1.33,0.2-7.2,98.67,92.8-99.8,50.00,0.0260\n")
    f.write("Ablation: Unpaired Single-Sequence PWM,150,65,66,9,10,86.67,77.2-92.6,12.00,6.4-21.3,49.62,0.6311\n")
    f.write("SpliceAI (DS_max >= 0.20 GRCh38),139,63,3,67,6,91.30,82.3-96.0,95.71,88.1-98.5,95.45,0.9333\n")

with open("results/clean_submission/Table4_PCAWG_Melanoma_WGS_Hotspots.csv", "w", encoding="utf-8") as f:
    f.write("Gene,GRCh38_Coord,hg19_Coord,Mutated_Donors,Total_WGS_Donors,Donor_Pct,Total_SNVs,UV_CT_GA_Fraction,Donor_Fold_Enrichment,p_NB,q_NB,Classification\n")
    f.write("TERT,chr5:1295113-1295135,chr5:1295228-1295250,24,70,34.3,25,19/25 (76.0%),178.8,4.44e-16,2.53e-08,Confirmed noncoding driver\n")
    f.write("AL358472.4,chr1:153990702-153990745,chr1:153963178-153963221,13,70,18.6,15,14/15 (93.3%),96.9,1.28e-11,7.27e-04,Recurrent UV-associated ETS-motif hotspot\n")
    f.write("AC040160.1,chr16:67227004-67227077,chr16:67260907-67260980,12,70,17.1,17,14/17 (82.4%),89.4,8.39e-11,4.78e-03,Recurrent UV-associated ETS-motif hotspot\n")
    f.write("TCF3,chr19:1652862-1653023,chr19:1652862-1653023,11,70,15.7,15,15/15 (100.0%),82.0,5.52e-10,3.14e-02,Recurrent UV-associated ETS-motif hotspot\n")

print("Exported clean submission tables to results/clean_submission/")
