###### Serial-section TLS composition and subtype stability ######

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


###### input and output ######

repo = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description="Serial-section TLS subtype analysis")
parser.add_argument("--data-dir", type=Path, default=repo / "data/10-serial")
parser.add_argument("--output-dir", type=Path, default=repo / "result/Figure2_serial")
parser.add_argument("--validate-model", action="store_true", help="Run slide-level LOSO and refit the classifier")
args = parser.parse_args()
data_dir = args.data_dir
output_dir = args.output_dir
output_dir.mkdir(parents=True, exist_ok=True)

fovs = pd.read_csv(data_dir / "TLS/FOV_manifest.csv", dtype={"slide": str})
cells = pd.read_csv(data_dir / "cell_types/serial_cell_metadata.csv.gz", dtype={"slide": str})
with open(data_dir / "model/model_parameters.json", encoding="utf-8") as handle:
    model_parameters = json.load(handle)
classifier = joblib.load(data_dir / "model/serial_TLS_classifier.joblib")

tls_names = [f"TLS{number}" for number in range(1, 8)]
if len(fovs) != 119 or fovs.duplicated(["TLS", "section_index"]).any():
    raise ValueError("Expected 119 unique TLS-section FOVs.")
for tls_name in tls_names:
    sections = fovs.loc[fovs["TLS"] == tls_name, "section_index"]
    if sorted(sections) != list(range(1, 18)):
        raise ValueError(f"{tls_name} must contain all 17 sections.")
if cells["cell_id"].duplicated().any():
    raise ValueError("Duplicate cell identifiers.")
if set(cells["fov"]) != set(fovs["fov"]):
    raise ValueError("Cell metadata and FOV manifest disagree.")


###### composition of each TLS from its member cells ######

tls_cells = cells.loc[cells["in_tls"]].copy()
counts = pd.crosstab([tls_cells["TLS"], tls_cells["section_index"]], tls_cells["cell_type"])
cell_types = [
    "B", "B-CD21", "B-GC", "PC", "CD3T", "CD4T", "CD8T", "Bile_Duct",
    "DC", "Endothelial", "Fibroblast", "HEV", "Lymphatic", "Macrophage", "Neutrophil",
]
if not set(counts.columns).issubset(cell_types):
    raise ValueError("Unexpected cell type within the TLS membership.")
counts = counts.reindex(columns=cell_types, fill_value=0)
fractions = counts.div(counts.sum(axis=1), axis=0)
fractions["n_member_cells"] = counts.sum(axis=1)
serial = fovs.merge(fractions.reset_index(), on=["TLS", "section_index"], validate="one_to_one")
if not serial["n_member_cells"].eq(serial["final_tls_cell_count"]).all():
    raise ValueError("TLS membership and cell counts disagree.")
serial = serial.sort_values(["TLS", "section_index"]).reset_index(drop=True)
serial.to_csv(output_dir / "serial_TLS_composition_recomputed.csv", index=False)

# B-lineage includes B, CD21-positive B, germinal-centre B and plasma cells.
features = pd.DataFrame({
    "B": 100 * (serial["B"] + serial["B-CD21"] + serial["B-GC"] + serial["PC"]),
    "CD4T": 100 * serial["CD4T"],
    "CD8T": 100 * serial["CD8T"],
    "Cholangiocyte_like": 100 * serial["Bile_Duct"],
    "DC": 100 * serial["DC"],
    "Endothelial": 100 * serial["Endothelial"],
    "Fibroblast": 100 * serial["Fibroblast"],
    "HEV": 100 * serial["HEV"],
    "Lymphatic": 100 * serial["Lymphatic"],
    "Macrophage": 100 * serial["Macrophage"],
    "Neutrophil": 100 * serial["Neutrophil"],
})
feature_names = model_parameters["features"]
features = features[feature_names]
if list(classifier.feature_names_in_) != feature_names:
    raise ValueError("Classifier feature order differs from the metadata.")


###### subtype probabilities ######

tnc_column = list(classifier.classes_).index("TNC-TLS")
serial["prob_TNC_logistic"] = classifier.predict_proba(features)[:, tnc_column]
serial["pred_subtype_logistic"] = np.where(
    serial["prob_TNC_logistic"] >= model_parameters["probability_cutoff"], "TNC-TLS", "BF-TLS"
)
serial.to_csv(output_dir / "serial_TLS_predictions_TLS.csv", index=False)

feature_table = pd.concat([serial[["TLS", "fov", "section_number"]], features], axis=1)
feature_table.to_csv(output_dir / "serial_TLS_classifier_features.csv", index=False)
probability_matrix = serial.pivot(index="TLS", columns="section_number", values="prob_TNC_logistic")
probability_matrix = probability_matrix.reindex(tls_names)
probability_matrix.to_csv(output_dir / "serial_TLS_probability_matrix.csv")


###### stability over the 17 sections ######

stability_rows = []
for tls_name in tls_names:
    one_tls = serial.loc[serial["TLS"] == tls_name].sort_values("section_index")
    labels = one_tls["pred_subtype_logistic"].to_numpy()
    subtype_counts = one_tls["pred_subtype_logistic"].value_counts()
    stability_rows.append({
        "TLS": tls_name,
        "n_sections": len(one_tls),
        "n_BF_TLS": int((labels == "BF-TLS").sum()),
        "n_TNC_TLS": int((labels == "TNC-TLS").sum()),
        "majority_subtype": subtype_counts.index[0],
        "majority_fraction": subtype_counts.iloc[0] / len(one_tls),
        "perfectly_stable": len(subtype_counts) == 1,
        "switch_count": int((labels[1:] != labels[:-1]).sum()),
        "median_probability_TNC": one_tls["prob_TNC_logistic"].median(),
    })
stability = pd.DataFrame(stability_rows)
stability.to_csv(output_dir / "serial_TLS_z_axis_stability_summary.csv", index=False)
print(stability.to_string(index=False))


###### optional classifier training and slide-level LOSO ######

def make_classifier():
    return Pipeline([
        ("scale", StandardScaler()),
        ("clf", LogisticRegression(
            penalty="l2", C=model_parameters["C"],
            class_weight=model_parameters["class_weight"],
            solver=model_parameters["solver"], max_iter=model_parameters["max_iter"],
            random_state=model_parameters["random_seed"],
        )),
    ])


if args.validate_model:
    training = pd.read_csv(data_dir / "model/training_TLS_composition.csv", dtype={"slide": str})
    training_columns = [
        "pct_B", "pct_CD4T", "pct_CD8T", "pct_Cholangiocyte", "pct_DC", "pct_Endothelial",
        "pct_Fibroblast", "pct_HEV", "pct_Lymphatic", "pct_Macrophage", "pct_Neutrophil",
    ]
    x_train = training[training_columns].copy()
    x_train.columns = feature_names
    y_train = training["tls_subtype"]
    if len(training) != 667 or training["slide"].nunique() != 29:
        raise ValueError("Expected 667 training TLSs from 29 slides.")
    if x_train.isna().any(axis=None):
        raise ValueError("Missing training cell percentages.")

    cv_probability = np.full(len(training), np.nan)
    for slide in sorted(training["slide"].unique()):
        test_rows = training["slide"] == slide
        fold_model = make_classifier()
        fold_model.fit(x_train.loc[~test_rows], y_train.loc[~test_rows])
        fold_tnc_column = list(fold_model.classes_).index("TNC-TLS")
        cv_probability[test_rows] = fold_model.predict_proba(x_train.loc[test_rows])[:, fold_tnc_column]

    cv_prediction = np.where(cv_probability >= 0.5, "TNC-TLS", "BF-TLS")
    cv_table = training[["tls_id", "slide", "fov", "tls_subtype"]].copy()
    cv_table["prob_TNC_TLS"] = cv_probability
    cv_table["predicted_TLS_subtype"] = cv_prediction
    cv_table.to_csv(output_dir / "serial_TLS_classifier_LOSO_predictions.csv", index=False)

    metrics = pd.DataFrame([{
        "n_training_TLS": len(training),
        "n_training_slides": training["slide"].nunique(),
        "accuracy": accuracy_score(y_train, cv_prediction),
        "balanced_accuracy": balanced_accuracy_score(y_train, cv_prediction),
        "roc_auc_TNC_TLS": roc_auc_score(y_train.eq("TNC-TLS"), cv_probability),
    }])
    metrics.to_csv(output_dir / "serial_TLS_classifier_LOSO_metrics.csv", index=False)

    fitted_model = make_classifier().fit(x_train, y_train)
    coefficient_error = fitted_model.named_steps["clf"].coef_ - classifier.named_steps["clf"].coef_
    if np.max(np.abs(coefficient_error)) > 1e-9:
        raise ValueError("Training does not reproduce the supplied classifier coefficients.")
    joblib.dump(fitted_model, output_dir / "serial_TLS_classifier_refitted.joblib")
    print(metrics.to_string(index=False))
