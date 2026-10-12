"""Reproduce manuscript TLS regions and export DBSCAN core/border QC.

DBSCAN density parameters
and post-clustering inclusion rules are distinct. Shared-border diagnostics
do not split clusters or establish whether a region is a single anatomical TLS.
"""

from pathlib import Path
import argparse
import json

parser = argparse.ArgumentParser(description="Reproduce final TLS DBSCAN and export candidate QC.")
parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parent / 'data')
parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent / 'result/Figure1_DBSCAN')
args = parser.parse_args()
table_dir = args.output_dir / '01-tables'
table_dir.mkdir(parents=True, exist_ok=True)

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull, Delaunay, cKDTree
from shapely.geometry import LineString
from shapely.ops import polygonize, unary_union
from sklearn.cluster import DBSCAN

# ---------- settings ----------
FINAL_CELLTYPE_COL = "celltype_l2_after_DC_refine"

DBSCAN_EPS = 80
# A core point requires 50 points, including itself, within eps.
DBSCAN_MIN_SAMPLES = 50
# This separate final inclusion rule is strict: a retained region has >=51 cells.
TLS_MIN_CELLS_EXCLUSIVE = 50
ALPHA_RADIUS = DBSCAN_EPS * 1.5

LYMPHOID_TYPES = ["B", "CD4T", "CD8T", "T-B"]
MYELOID_TYPES = ["Macrophage", "DC", "Neutrophil"]
STROMA_TYPES = ["Fibroblast", "Endothelial", "HEV", "Lymphatic"]
OTHER_TYPES = ["Cholangiocyte", "Unassigned"]
TLS_INPUT_TYPES = LYMPHOID_TYPES + MYELOID_TYPES + STROMA_TYPES + OTHER_TYPES


# ---------- read final cell metadata ----------
cells = pd.read_csv(
    args.data_dir / "02-CODEX_Cohort2/cell_types/final_codex_celltyped_obs.csv",
    usecols=["cellid", "slide", "fov", "label", "centroid-0", "centroid-1", FINAL_CELLTYPE_COL],
    dtype={"cellid": str, "slide": str, "fov": str},
    low_memory=False,
)
cells = cells.rename(
    columns={
        "centroid-0": "centroid.0",
        "centroid-1": "centroid.1",
        FINAL_CELLTYPE_COL: "celltype_l2",
    }
)

cells["label"] = pd.to_numeric(cells["label"], errors="coerce")
cells["centroid.0"] = pd.to_numeric(cells["centroid.0"], errors="coerce")
cells["centroid.1"] = pd.to_numeric(cells["centroid.1"], errors="coerce")
cells["celltype_l2"] = cells["celltype_l2"].astype(str)


# ---------- DBSCAN input cells ----------
dbscan_cells = cells[cells["celltype_l2"].isin(TLS_INPUT_TYPES)].copy()
dbscan_cells = dbscan_cells[~dbscan_cells["fov"].str.endswith("B", na=False)].copy()
dbscan_cells = dbscan_cells.dropna(subset=["label", "centroid.0", "centroid.1"]).copy()
dbscan_cells["label"] = dbscan_cells["label"].astype(int)
dbscan_cells["dbscan_label"] = -1
dbscan_cells["is_dbscan_core"] = False
dbscan_cells["is_shared_border"] = False
shared_border_rows = []
dbscan_cells["tls_id"] = "Noise"
dbscan_cells["tls_number"] = pd.NA


# ---------- run DBSCAN within each FOV ----------
for (slide, fov), one_fov in dbscan_cells.groupby(["slide", "fov"], sort=False):
    if one_fov.shape[0] < DBSCAN_MIN_SAMPLES:
        continue

    coordinates = one_fov[["centroid.0", "centroid.1"]].to_numpy(dtype=float)
    model = DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN_SAMPLES).fit(coordinates)
    labels = model.labels_
    core = np.zeros(len(one_fov), dtype=bool)
    core[model.core_sample_indices_] = True
    dbscan_cells.loc[one_fov.index, "dbscan_label"] = labels
    dbscan_cells.loc[one_fov.index, "is_dbscan_core"] = core

    # A non-core cell can neighbor core points from several distinct clusters.
    # Its single assigned label depends on the input order; keep that order.
    core_idx = np.flatnonzero(core)
    border_idx = np.flatnonzero((labels >= 0) & ~core)
    if len(core_idx) and len(border_idx):
        core_tree = cKDTree(coordinates[core_idx])
        neighbors = core_tree.query_ball_point(coordinates[border_idx], r=DBSCAN_EPS)
        for i, adjacent_core in zip(border_idx, neighbors):
            adjacent_labels = np.unique(labels[core_idx[np.asarray(adjacent_core, dtype=int)]])
            if len(adjacent_labels) > 1:
                dbscan_cells.loc[one_fov.index[i], "is_shared_border"] = True
                shared_border_rows.append({
                    "cellid": one_fov.iloc[i]["cellid"], "slide": slide, "fov": fov,
                    "assigned_dbscan_label": int(labels[i]),
                    "adjacent_core_cluster_labels": ";".join(map(str, adjacent_labels.tolist())),
                })


# ---------- name TLS clusters ----------
for (slide, fov), one_fov in dbscan_cells[dbscan_cells["dbscan_label"] >= 0].groupby(["slide", "fov"], sort=False):
    labels_by_size = one_fov["dbscan_label"].value_counts().sort_values(ascending=False).index.tolist()

    for tls_number, dbscan_label in enumerate(labels_by_size, start=1):
        idx = one_fov.index[one_fov["dbscan_label"] == dbscan_label]
        tls_id = f"{slide}_{fov}_TLS{tls_number}"
        dbscan_cells.loc[idx, "tls_id"] = tls_id
        dbscan_cells.loc[idx, "tls_number"] = tls_number

dbscan_cells["tls_number"] = dbscan_cells["tls_number"].astype("Int64")
tls_cells_all = dbscan_cells[dbscan_cells["dbscan_label"] >= 0].copy()


# ---------- TLS composition ----------
composition = (
    tls_cells_all.groupby(["tls_id", "slide", "fov", "tls_number", "dbscan_label"])
    .size()
    .reset_index(name="n_cells")
)

celltype_count = (
    tls_cells_all.groupby(["tls_id", "celltype_l2"])
    .size()
    .reset_index(name="n")
    .pivot_table(index="tls_id", columns="celltype_l2", values="n", fill_value=0, aggfunc="sum")
    .reset_index()
)

composition = composition.merge(celltype_count, on="tls_id", how="left")

for celltype in TLS_INPUT_TYPES:
    if celltype not in composition.columns:
        composition[celltype] = 0
    composition[celltype] = composition[celltype].astype(int)
    composition[f"pct_{celltype}"] = composition[celltype] / composition["n_cells"] * 100

composition["n_lymphoid"] = composition[LYMPHOID_TYPES].sum(axis=1)
composition["n_myeloid"] = composition[MYELOID_TYPES].sum(axis=1)
composition["n_stroma"] = composition[STROMA_TYPES].sum(axis=1)
composition["n_T_combined"] = composition[["CD4T", "CD8T", "T-B"]].sum(axis=1)

composition["pct_lymphoid"] = composition["n_lymphoid"] / composition["n_cells"] * 100
composition["pct_myeloid"] = composition["n_myeloid"] / composition["n_cells"] * 100
composition["pct_stroma"] = composition["n_stroma"] / composition["n_cells"] * 100
composition["pct_T_combined"] = composition["n_T_combined"] / composition["n_cells"] * 100
composition["pct_T_Chol_Neu"] = (
    composition["CD4T"]
    + composition["CD8T"]
    + composition["T-B"]
    + composition["Cholangiocyte"]
    + composition["Neutrophil"]
) / composition["n_cells"] * 100


# ---------- final TLS inclusion ----------
keep_tls = (
    (composition["n_cells"] > TLS_MIN_CELLS_EXCLUSIVE)
    & (composition["pct_Unassigned"] <= 30)
    & (composition["pct_lymphoid"] >= 10)
    & (composition["pct_stroma"] <= 60)
)

# ---------- candidate QC before final inclusion ----------
reference = pd.read_csv(args.data_dir / "02-CODEX_Cohort2/TLS/tls_subtype_assignment_final.csv")
candidate_qc = composition.copy()
core_count = tls_cells_all.groupby("tls_id")["is_dbscan_core"].sum()
shared_count = tls_cells_all.groupby("tls_id")["is_shared_border"].sum()
candidate_qc["n_core"] = candidate_qc["tls_id"].map(core_count).astype(int)
candidate_qc["n_border"] = candidate_qc["n_cells"] - candidate_qc["n_core"]
candidate_qc["n_assigned_shared_border"] = candidate_qc["tls_id"].map(shared_count).astype(int)
candidate_qc["pass_size_gt50"] = candidate_qc["n_cells"] > TLS_MIN_CELLS_EXCLUSIVE
candidate_qc["pass_Unassigned_le30pct"] = candidate_qc["pct_Unassigned"] <= 30
candidate_qc["pass_lymphoid_ge10pct"] = candidate_qc["pct_lymphoid"] >= 10
candidate_qc["pass_stroma_le60pct"] = candidate_qc["pct_stroma"] <= 60
candidate_qc["passes_DBSCAN_QC"] = keep_tls
candidate_qc["in_manuscript_cohort"] = candidate_qc["tls_id"].isin(reference["tls_id"])
rules = ["pass_size_gt50", "pass_Unassigned_le30pct", "pass_lymphoid_ge10pct", "pass_stroma_le60pct"]
candidate_qc["failed_rules"] = candidate_qc.apply(
    lambda row: ";".join(rule for rule in rules if not row[rule]), axis=1
)
candidate_qc.to_csv(table_dir / "tls_candidate_qc_final_rules.csv", index=False)
pd.DataFrame(shared_border_rows, columns=["cellid", "slide", "fov", "assigned_dbscan_label", "adjacent_core_cluster_labels"]).to_csv(
    table_dir / "dbscan_shared_border_cells.csv", index=False
)
tls_cells_all = tls_cells_all.drop(columns=["is_dbscan_core", "is_shared_border"])

composition = composition[keep_tls].copy()
final_tls_ids = composition["tls_id"].tolist()
tls_cells = tls_cells_all[tls_cells_all["tls_id"].isin(final_tls_ids)].copy()


# The DBSCAN output also contains regions outside the manuscript TA-TLS cohort.
if not set(reference["tls_id"]).issubset(set(composition["tls_id"])):
    raise ValueError("Manuscript TLS IDs are missing from DBSCAN output; inspect candidate QC.")
counts = composition.set_index("tls_id")["n_cells"].reindex(reference["tls_id"])
if not counts.eq(reference.set_index("tls_id")["n_cells"].reindex(counts.index)).all():
    raise ValueError("Recomputed per-TLS cell counts differ from the archived cohort.")

# ---------- TLS polygons ----------
polygon_rows = []

for tls_id, one_tls in tls_cells.groupby("tls_id", sort=False):
    points = one_tls[["centroid.0", "centroid.1"]].to_numpy(dtype=float)
    polygon_wkt = ""
    boundary_method = "not_available"
    area_pixel2 = np.nan

    if points.shape[0] >= 3:
        hull = ConvexHull(points)
        hull_points = points[hull.vertices]
        hull_points = np.vstack([hull_points, hull_points[0]])
        polygon_wkt = "POLYGON (({}))".format(", ".join([f"{x} {y}" for x, y in hull_points]))
        boundary_method = "convex_hull"
        area_pixel2 = float(hull.volume)

    if points.shape[0] >= 4:
        tri = Delaunay(points)
        edges = set()

        for simplex in tri.simplices:
            pa = points[simplex[0]]
            pb = points[simplex[1]]
            pc = points[simplex[2]]

            a = np.linalg.norm(pb - pc)
            b = np.linalg.norm(pa - pc)
            c = np.linalg.norm(pa - pb)
            s = (a + b + c) / 2
            triangle_area = max(s * (s - a) * (s - b) * (s - c), 0) ** 0.5
            if triangle_area == 0:
                continue

            circum_radius = a * b * c / (4 * triangle_area)
            if circum_radius <= ALPHA_RADIUS:
                edges.add(tuple(sorted((simplex[0], simplex[1]))))
                edges.add(tuple(sorted((simplex[1], simplex[2]))))
                edges.add(tuple(sorted((simplex[2], simplex[0]))))

        if len(edges) > 0:
            lines = [LineString([points[i], points[j]]) for i, j in edges]
            alpha_polygon = unary_union(list(polygonize(unary_union(lines))))
            if (not alpha_polygon.is_empty) and alpha_polygon.area > 0:
                polygon_wkt = alpha_polygon.wkt
                boundary_method = "alpha_shape"
                area_pixel2 = float(alpha_polygon.area)

    first = one_tls.iloc[0]
    polygon_rows.append(
        {
            "tls_id": tls_id,
            "slide": first["slide"],
            "fov": first["fov"],
            "tls_number": first["tls_number"],
            "boundary_method": boundary_method,
            "alpha_radius": ALPHA_RADIUS,
            "tls_area_pixel2": area_pixel2,
            "polygon_wkt": polygon_wkt,
        }
    )

polygon = pd.DataFrame(polygon_rows)


# ---------- save final DBSCAN TLS outputs ----------
Path(table_dir).mkdir(
    parents=True,
    exist_ok=True,
)

tls_cells.to_csv(
    table_dir / 'dbscan_tls_cells_final.csv.gz',
    index=False,
    compression="gzip",
)
composition.to_csv(
    table_dir / 'tls_cell_composition_final.csv',
    index=False,
)
polygon.to_csv(
    table_dir / 'tls_polygon_final.csv',
    index=False,
)


# ---------- report ----------
print("DBSCAN regions after QC (all tissue regions):", composition["tls_id"].nunique())
print("DBSCAN cells after QC:", tls_cells.shape[0])
print("DBSCAN polygons:", polygon["tls_id"].nunique())
print("Manuscript TLS subset:", len(reference))
print("Saved TLS tables:", table_dir)
parameters = {
    "eps_pixels": DBSCAN_EPS,
    "min_samples_including_self": DBSCAN_MIN_SAMPLES,
    "final_cluster_size_rule": "n_cells > 50",
    "candidate_regions": len(candidate_qc),
    "DBSCAN_regions_after_QC_all_tissue_regions": len(composition),
    "DBSCAN_cells_after_QC_all_tissue_regions": len(tls_cells),
    "manuscript_TLS": len(reference),
    "manuscript_cells": int(counts.sum()),
    "shared_border_cells": len(shared_border_rows),
    "morphological_fusion_review": "not_performed_by_this_algorithm",
}
(args.output_dir / "DBSCAN_parameters_and_counts.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
