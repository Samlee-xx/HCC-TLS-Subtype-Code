###### Figure S13J: Cohort6 TLS subtype comparison ######

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import roc_auc_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import RobustScaler


###### input and output ######

repo = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description="Cohort6 TLS subtype comparison")
parser.add_argument("--data-dir", type=Path, default=repo / "data/09-TLSProfiler")
parser.add_argument("--output-dir", type=Path, default=repo / "result/FigureS13/Cohort6")
parser.add_argument("--cell-table", type=Path)
parser.add_argument("--membership", type=Path)
args = parser.parse_args()

input_dir = args.data_dir / "Cohort6_validation"
output_dir = args.output_dir
output_dir.mkdir(parents=True, exist_ok=True)

tls = pd.read_csv(input_dir / "tls_composition.csv", dtype={"slide": str})
samples = pd.read_csv(input_dir / "sample_manifest.csv", dtype={"slide": str})
with open(input_dir / "reference_parameters.json", encoding="utf-8-sig") as handle:
    parameters = json.load(handle)
with open(args.data_dir / "model/tlsprofiler_log1p_BTN.json", encoding="utf-8-sig") as handle:
    model = json.load(handle)

if tls["tls_id"].duplicated().any() or samples["slide"].duplicated().any():
    raise ValueError("Duplicate TLS or sample identifiers.")
if len(tls) != 184 or len(samples) != 103 or tls["slide"].nunique() != 52:
    raise ValueError("Expected 184 TLSs from 52 TLS-positive samples within 103 samples.")

tls_per_slide = tls.groupby("slide").size()
sample_counts = samples["slide"].map(tls_per_slide).fillna(0).astype(int)
if not sample_counts.eq(samples["n_analysis_tls"]).all():
    raise ValueError("TLS counts disagree with the sample manifest.")
if not tls["slide"].isin(samples["slide"]).all():
    raise ValueError("TLS table contains a sample absent from the manifest.")
if (args.cell_table is None) != (args.membership is None):
    raise ValueError("Supply both --cell-table and --membership for cell clustering.")

seed = parameters["random_seed"]


###### optional cell clustering from the external cell and membership tables ######

if args.cell_table is not None:
    markers = parameters["markers"]
    cell_columns = ["fov", "label", "cell_size"] + markers
    cells = pd.read_csv(args.cell_table, usecols=cell_columns, dtype={"fov": str})

    keep_sample = cells["fov"].isin(samples["slide"])
    large_enough = cells["cell_size"] > parameters["cell_size_min_exclusive"]
    small_enough = cells["cell_size"] < parameters["cell_size_max_exclusive"]
    cells = cells.loc[keep_sample & large_enough & small_enough].copy()
    cells = cells.dropna(subset=markers)
    cells["label"] = pd.to_numeric(cells["label"], errors="raise").astype(np.int64)

    if cells.duplicated(["fov", "label"]).any():
        raise ValueError("Duplicate cell identifiers.")
    if not np.isfinite(cells[markers].to_numpy()).all():
        raise ValueError("Non-finite marker measurements.")
    if set(cells["fov"]) != set(samples["slide"]):
        raise ValueError("Every sample must supply eligible cells.")

    membership = pd.read_csv(
        args.membership, usecols=["slide", "label", "tls_id"], dtype={"slide": str}
    )
    membership = membership.loc[membership["tls_id"].isin(tls["tls_id"])].copy()
    membership["label"] = pd.to_numeric(membership["label"], errors="raise").astype(np.int64)
    if membership.duplicated(["tls_id", "slide", "label"]).any():
        raise ValueError("Duplicate cell membership within a TLS.")

    membership_counts = membership.groupby(["slide", "tls_id"]).size()
    membership_counts = membership_counts.rename("membership_n")
    denominators = tls.set_index(["slide", "tls_id"])[["denominator_n"]]
    denominators = denominators.join(membership_counts)
    if not denominators["denominator_n"].eq(denominators["membership_n"]).all():
        raise ValueError("Cell membership and TLS denominators disagree.")

    tls_cells = membership.merge(
        cells[["fov", "label"] + markers],
        left_on=["slide", "label"], right_on=["fov", "label"],
        how="left", validate="many_to_one",
    )
    if tls_cells[markers].isna().any(axis=None):
        raise ValueError("TLS cells lack eligible marker measurements.")

    # Sample up to 2,000 cells per slide for fitting.
    sample_parts = []
    for slide, slide_cells in cells.groupby("fov", sort=True):
        n_cells = min(len(slide_cells), parameters["sample_per_slide"])
        sample_parts.append(slide_cells.sample(n=n_cells, random_state=seed, replace=False))
    sampled = pd.concat(sample_parts, ignore_index=True)
    sampled = sampled.sample(frac=1, random_state=seed).reset_index(drop=True)

    # Cluster the four marker measurements after robust scaling.
    scaler = RobustScaler(quantile_range=tuple(parameters["scaler_quantiles"]))
    sample_values = sampled[markers].to_numpy(float)
    tls_values = tls_cells[markers].to_numpy(float)
    sample_scaled = scaler.fit_transform(sample_values).astype(np.float32)
    tls_scaled = scaler.transform(tls_values).astype(np.float32)

    clustering = MiniBatchKMeans(
        n_clusters=parameters["n_clusters"],
        init="k-means++",
        n_init=parameters["kmeans_n_init"],
        max_iter=parameters["kmeans_max_iter"],
        batch_size=parameters["kmeans_batch_size"],
        reassignment_ratio=parameters["kmeans_reassignment_ratio"],
        random_state=seed,
    )
    sampled["reference_cluster"] = clustering.fit_predict(sample_scaled)
    tls_cells["reference_cluster"] = clustering.predict(tls_scaled)

    profile = sampled.groupby("reference_cluster")[markers].median()
    profile["sampled_cells"] = sampled.groupby("reference_cluster").size()
    profile = profile.reset_index()

    # The component with the highest mean defines the high-expression state.
    for marker in markers:
        mixture = GaussianMixture(
            n_components=parameters["gmm_components"][marker],
            covariance_type="full",
            n_init=parameters["gmm_n_init"],
            max_iter=parameters["gmm_max_iter"],
            random_state=seed,
        )
        mixture.fit(sampled[[marker]].to_numpy(float))
        high_component = np.argmax(mixture.means_.ravel())
        marker_probability = mixture.predict_proba(profile[[marker]].to_numpy(float))
        profile["prob_high_" + marker] = marker_probability[:, high_component]

    # CD45 separates immune cells; CD20, CD3 and CD66B assign the lineage.
    lineage_columns = ["prob_high_CD20", "prob_high_CD3", "prob_high_CD66B"]
    lineage_probability = profile[lineage_columns].to_numpy()
    lineage_names = np.array(["B_cell", "T_cell", "Neutrophil"])
    winning_lineage = lineage_names[lineage_probability.argmax(axis=1)]
    immune = profile["prob_high_CD45"] >= parameters["immune_posterior_cutoff"]
    lineage_positive = lineage_probability.max(axis=1) >= parameters["lineage_posterior_cutoff"]

    profile["reference_cell_type"] = "Other_cell"
    profile.loc[immune, "reference_cell_type"] = "Other_immune"
    assigned = immune & lineage_positive
    profile.loc[assigned, "reference_cell_type"] = winning_lineage[assigned]

    cluster_types = profile.set_index("reference_cluster")["reference_cell_type"]
    tls_cells["reference_cell_type"] = tls_cells["reference_cluster"].map(cluster_types)
    counts = tls_cells.groupby(["slide", "tls_id", "reference_cell_type"]).size()
    counts = counts.unstack(fill_value=0)
    cell_types = ["B_cell", "T_cell", "Neutrophil", "Other_immune", "Other_cell"]
    counts = counts.reindex(columns=cell_types, fill_value=0)
    counts["reference_denominator_n"] = counts.sum(axis=1)
    counts["reference_B_pct"] = 100 * counts["B_cell"] / counts["reference_denominator_n"]
    counts["reference_T_pct"] = 100 * counts["T_cell"] / counts["reference_denominator_n"]
    counts["reference_Neu_pct"] = 100 * counts["Neutrophil"] / counts["reference_denominator_n"]
    counts = counts.reset_index()

    automated_columns = [
        "slide", "tls_id", "denominator_n",
        "automated_B_pct", "automated_T_pct", "automated_Neu_pct",
    ]
    tls = counts.merge(tls[automated_columns], on=["slide", "tls_id"], validate="one_to_one")
    profile.to_csv(output_dir / "cluster_annotation.csv", index=False)


###### TLSProfiler probabilities for the two cell-typing methods ######

if len(tls) != 184 or not tls["reference_denominator_n"].eq(tls["denominator_n"]).all():
    raise ValueError("The two methods must use the same TLSs and cell denominators.")

reference_columns = ["reference_B_pct", "reference_T_pct", "reference_Neu_pct"]
automated_columns = ["automated_B_pct", "automated_T_pct", "automated_Neu_pct"]
percentages = tls[reference_columns + automated_columns].to_numpy()
if not np.isfinite(percentages).all() or (percentages < 0).any() or (percentages > 100).any():
    raise ValueError("Cell percentages must be finite and within 0-100.")
if (tls[reference_columns].sum(axis=1) > 100 + 1e-8).any():
    raise ValueError("Reference cell percentages exceed 100.")
if (tls[automated_columns].sum(axis=1) > 100 + 1e-8).any():
    raise ValueError("Automated cell percentages exceed 100.")
if model["transform"] != "log1p" or model["input_scale"] != "percentage_0_to_100":
    raise ValueError("TLSProfiler requires log1p-transformed cell percentages.")
if model["positive_class"] != "TNC-TLS":
    raise ValueError("TLSProfiler positive class must be TNC-TLS.")

coefficients = model["coefficients"]
reference_logit = (
    model["intercept"]
    + coefficients["B_pct"] * np.log1p(tls["reference_B_pct"])
    + coefficients["T_pct"] * np.log1p(tls["reference_T_pct"])
    + coefficients["Neutrophil_pct"] * np.log1p(tls["reference_Neu_pct"])
)
tls["reference_prob_TNC"] = expit(reference_logit)
tls["reference_tls_subtype"] = np.where(
    tls["reference_prob_TNC"] >= model["probability_cutoff"], "TNC-TLS", "BF-TLS"
)

automated_logit = (
    model["intercept"]
    + coefficients["B_pct"] * np.log1p(tls["automated_B_pct"])
    + coefficients["T_pct"] * np.log1p(tls["automated_T_pct"])
    + coefficients["Neutrophil_pct"] * np.log1p(tls["automated_Neu_pct"])
)
tls["automated_prob_TNC"] = expit(automated_logit)
tls["automated_tls_subtype"] = np.where(
    tls["automated_prob_TNC"] >= model["probability_cutoff"], "TNC-TLS", "BF-TLS"
)


###### AUC and slide-bootstrap confidence interval ######

target = tls["reference_tls_subtype"].eq("TNC-TLS").to_numpy(dtype=int)
probability = tls["automated_prob_TNC"].to_numpy()
auc_value = roc_auc_score(target, probability)

slides = tls["slide"].unique()
slide_rows = {}
for slide in slides:
    slide_rows[slide] = np.flatnonzero(tls["slide"].eq(slide))

rng = np.random.default_rng(seed)
bootstrap_auc = []
for iteration in range(parameters["bootstrap_replicates"]):
    sampled_slides = rng.choice(slides, size=len(slides), replace=True)
    sampled_parts = [slide_rows[slide] for slide in sampled_slides]
    sampled_rows = np.concatenate(sampled_parts)
    bootstrap_target = target[sampled_rows]
    if np.unique(bootstrap_target).size < 2:
        continue
    bootstrap_auc.append(roc_auc_score(bootstrap_target, probability[sampled_rows]))

if not bootstrap_auc:
    raise ValueError("No evaluable slide-bootstrap samples.")
ci_low, ci_high = np.quantile(bootstrap_auc, [0.025, 0.975])

summary = pd.DataFrame([{
    "n_tls": len(tls),
    "n_tls_positive_samples": len(slides),
    "n_reference_BF": int((target == 0).sum()),
    "n_reference_TNC": int((target == 1).sum()),
    "auc": auc_value,
    "CI_low": ci_low,
    "CI_high": ci_high,
    "n_valid_bootstrap": len(bootstrap_auc),
}])
tls.to_csv(output_dir / "tls_comparison.csv", index=False)
summary.to_csv(output_dir / "roc_summary.csv", index=False)
print(summary.to_string(index=False))
