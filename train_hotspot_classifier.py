"""
train_hotspot_classifier.py
Trains a Random Forest Machine Learning model to classify cancer types 
based purely on noncoding mutation genomic features.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, accuracy_score

RESULTS_DIR = Path(r"D:\noncoding mutation finding surver\results\tables")
MATRIX_FILE = RESULTS_DIR / "ml_hotspot_training_matrix.csv"

def train_model():
    if not MATRIX_FILE.exists():
        print(f"Error: {MATRIX_FILE.name} not found.")
        return

    # 1. Load the data
    print("Loading ML feature matrix...")
    df = pd.read_csv(MATRIX_FILE)

    # 2. Define Features (X) and Target (y)
    # We are asking the model to predict the Cancer Type based ONLY on genomic topology and motifs
    feature_cols = [
        "Mutated_Donors", 
        "Fold_Enrichment", 
        "Absolute_Distance_to_TSS", 
        "Is_Promoter", 
        "Is_Enhancer", 
        "Has_TF_Impact", 
        "Motif_Gain", 
        "Motif_Loss"
    ]
    
    X = df[feature_cols]
    y = df["Cancer_Type"]

    # 3. Split into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)

    # 4. Train the Random Forest
    print("\nTraining Random Forest Classifier...")
    clf = RandomForestClassifier(n_estimators=100, random_state=42, max_depth=5)
    clf.fit(X_train, y_train)

    # 5. Evaluate the model
    y_pred = clf.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    
    print("\n" + "=" * 60)
    print(f" MODEL EVALUATION (ACCURACY: {accuracy * 100:.1f}%)")
    print("=" * 60)
    print(classification_report(y_test, y_pred))

    # 6. Extract Feature Importances
    print("=" * 60)
    print(" WHAT DID THE MODEL LEARN? (FEATURE IMPORTANCE)")
    print("=" * 60)
    
    importances = clf.feature_importances_
    feat_imp = pd.DataFrame({
        "Feature": feature_cols,
        "Importance": importances
    }).sort_values(by="Importance", ascending=False)
    
    for _, row in feat_imp.iterrows():
        print(f"{row['Feature']:>25}: {row['Importance']:.4f}")
    print("=" * 60)

if __name__ == "__main__":
    train_model()