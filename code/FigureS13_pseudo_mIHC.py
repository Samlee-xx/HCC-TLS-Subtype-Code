###### Figure S13E: pseudo-mIHC cell types and TLS predictions ######

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.metrics import roc_auc_score


###### input and output ######

repo = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description="Pseudo-mIHC TLS classification")
parser.add_argument("--data-dir", type=Path, default=repo / "data/09-TLSProfiler")
parser.add_argument("--output-dir", type=Path, default=repo / "result/FigureS13/pseudo_mIHC")
args = parser.parse_args()

input_dir = args.data_dir / "pseudo_mIHC"
output_dir = args.output_dir
output_dir.mkdir(parents=True, exist_ok=True)

tls = pd.read_csv(input_dir / "tls_manifest.csv", dtype={"slide": str})
cells = pd.read_csv(input_dir / "tls_cell_marker_folds.csv.gz")
with open(input_dir / "cell_rule.json", encoding="utf-8-sig") as handle:
    cell_rule = json.load(handle)
with open(args.data_dir / "model/tlsprofiler_log1p_BTN.json", encoding="utf-8-sig") as handle:
    model = json.load(handle)

if tls.duplicated(["cohort", "tls_id"]).any():
    raise ValueError("Duplicate TLS identifiers.")
if cells.duplicated(["cohort", "tls_id", "cell_row_id"]).any():
    raise ValueError("Duplicate cell membership within a TLS.")

cell_tls = cells[["cohort", "tls_id"]].drop_duplicates()
membership = tls[["cohort", "tls_id"]].merge(
    cell_tls, on=["cohort", "tls_id"], how="outer", indicator=True
)
if not membership["_merge"].eq("both").all():
    raise ValueError("Cell membership and TLS manifest have different TLS sets.")


###### marker thresholds ######

marker_columns = ["fold_CD45", "fold_CD20", "fold_CD3", "fold_CD66B"]
if not np.isfinite(cells[marker_columns].to_numpy()).all():
    raise ValueError("Non-finite marker measurements.")

rule_values = np.array(list(cell_rule.values()), dtype=float)
if not np.isfinite(rule_values).all() or (rule_values <= 0).any():
    raise ValueError("Cell-rule parameters must be finite and positive.")

cd45_positive = cells["fold_CD45"] >= cell_rule["CD45_multiplier"]
cd20_positive = cells["fold_CD20"] >= cell_rule["CD20_multiplier"]
cd3_positive = cells["fold_CD3"] >= cell_rule["CD3_multiplier"]
cd66b_positive = cells["fold_CD66B"] >= cell_rule["CD66B_multiplier"]

cd20_signal = cells["fold_CD20"] / cell_rule["CD20_multiplier"]
cd3_signal = cells["fold_CD3"] / cell_rule["CD3_multiplier"]
cd66b_signal = cells["fold_CD66B"] / cell_rule["CD66B_multiplier"]


###### hierarchical cell annotation ######

# CD66B-positive immune cells are assigned first, using the dominance margin.
no_b_or_t_marker = ~cd20_positive & ~cd3_positive
largest_bt_signal = np.maximum(cd20_signal, cd3_signal)
cd66b_dominant = cd66b_signal >= cell_rule["CD66B_dominance_margin"] * largest_bt_signal
neutrophil = cd45_positive & cd66b_positive & (no_b_or_t_marker | cd66b_dominant)

remaining_immune = cd45_positive & ~neutrophil
b_dominant = cd20_signal >= cell_rule["CD20_over_CD3_dominance_scale"] * cd3_signal
b_cell = remaining_immune & cd20_positive & (~cd3_positive | b_dominant)
t_cell = remaining_immune & cd3_positive & (~cd20_positive | ~b_dominant)


###### cell counts and percentages within each TLS ######

cell_counts = cells[["cohort", "tls_id"]].copy()
cell_counts["denominator_n"] = 1
cell_counts["n_immune"] = cd45_positive.astype(int)
cell_counts["n_B"] = b_cell.astype(int)
cell_counts["n_T"] = t_cell.astype(int)
cell_counts["n_Neu"] = neutrophil.astype(int)
tls_counts = cell_counts.groupby(["cohort", "tls_id"], sort=False).sum()
tls_counts = tls_counts.reset_index()

tls = tls.drop(columns="union_wkt")
tls = tls.merge(tls_counts, on=["cohort", "tls_id"], validate="one_to_one")

# All cells within the TLS contribute to the denominator.
tls["B_pct"] = 100 * tls["n_B"] / tls["denominator_n"]
tls["T_pct"] = 100 * tls["n_T"] / tls["denominator_n"]
tls["Neu_pct"] = 100 * tls["n_Neu"] / tls["denominator_n"]
tls["immune_pct"] = 100 * tls["n_immune"] / tls["denominator_n"]


###### TLSProfiler probabilities ######

if model["transform"] != "log1p" or model["input_scale"] != "percentage_0_to_100":
    raise ValueError("TLSProfiler requires log1p-transformed cell percentages.")
if model["positive_class"] != "TNC-TLS":
    raise ValueError("TLSProfiler positive class must be TNC-TLS.")

coefficients = model["coefficients"]
tls["locked_lr_logit"] = (
    model["intercept"]
    + coefficients["B_pct"] * np.log1p(tls["B_pct"])
    + coefficients["T_pct"] * np.log1p(tls["T_pct"])
    + coefficients["Neutrophil_pct"] * np.log1p(tls["Neu_pct"])
)
tls["locked_lr_prob_TNC"] = expit(tls["locked_lr_logit"])
tls["predicted_subtype"] = np.where(
    tls["locked_lr_prob_TNC"] >= model["probability_cutoff"], "TNC-TLS", "BF-TLS"
)


###### AUC for each cohort ######

if not tls["label_subtype"].isin(["BF_TLS", "TNC_TLS"]).all():
    raise ValueError("Unexpected reference TLS subtype.")
tls["target_TNC"] = tls["label_subtype"].eq("TNC_TLS").astype(int)

auc_rows = []
for cohort, cohort_tls in tls.groupby("cohort", sort=True):
    auc_value = roc_auc_score(cohort_tls["target_TNC"], cohort_tls["locked_lr_prob_TNC"])
    auc_rows.append({
        "cohort": cohort,
        "n_tls": len(cohort_tls),
        "n_slides": cohort_tls["slide"].nunique(),
        "n_BF_TLS": int(cohort_tls["target_TNC"].eq(0).sum()),
        "n_TNC_TLS": int(cohort_tls["target_TNC"].eq(1).sum()),
        "auc": auc_value,
    })

auc_summary = pd.DataFrame(auc_rows)
tls.to_csv(output_dir / "tls_predictions.csv", index=False)
auc_summary.to_csv(output_dir / "auc_summary.csv", index=False)
print(auc_summary.to_string(index=False))
