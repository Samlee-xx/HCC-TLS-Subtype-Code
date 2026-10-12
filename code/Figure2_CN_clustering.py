"""Radius-based CN clustering with K=9.

Fit on 696 candidate regions;
use Figure2_CN.py for manuscript statistics on the final 667 TLSs.
"""

import os
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
from sklearn.cluster import MiniBatchKMeans

parser = argparse.ArgumentParser(description="Fit the archived K=9 CN model on 696 candidate regions.")
parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent / "data")
parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "result/Figure2_CN/upstream_recomputed")
args = parser.parse_args()
upstream_dir = args.data_dir / "02-CODEX_Cohort2/CN/upstream_696"
cn_save_dir = str(args.output_dir)
os.makedirs(cn_save_dir, exist_ok=True)
cn_tls_cells_path = upstream_dir / "DBSCAN_TLS_cells_filtered_with_L2_after_DC_refine.csv.gz"
cn_tls_assignment_path = upstream_dir / "tls_assignment.csv"

CN_RADIUS = 100
CN_FINAL_K = 9
CN_RANDOM_STATE = 20260618
CN_BATCH_SIZE = 20000

cn_cell_type_order = [
    "B",
    "CD4T",
    "CD8T",
    "T-B",
    "Cholangiocyte",
    "DC",
    "Endothelial",
    "Fibroblast",
    "HEV",
    "Lymphatic",
    "Macrophage",
    "Neutrophil",
    "Unassigned",
]

print("CN output folder:", cn_save_dir)
print("CN radius:", CN_RADIUS)


# %%
# 10.1 Read rank02 TLS assignment and cell table

cn_tls_assignment = pd.read_csv(cn_tls_assignment_path)
cn_tls_assignment = cn_tls_assignment[[
    "tls_id",
    "slide",
    "fov",
    "slide_fov_key",
    "fov_region",
    "cluster_raw",
    "tls_subtype_rank02",
]].copy()
cn_tls_assignment["tls_id"] = cn_tls_assignment["tls_id"].astype(str)

cn_tls_ids = set(cn_tls_assignment["tls_id"].tolist())

cn_use_cols = [
    "cellid",
    "slide",
    "fov",
    "label",
    "centroid.0",
    "centroid.1",
    "tls_id",
    "celltype_l2",
]

cn_cells_list = []
for one_chunk in pd.read_csv(
    cn_tls_cells_path,
    usecols=cn_use_cols,
    chunksize=300000,
    dtype={
        "cellid": str,
        "slide": str,
        "fov": str,
        "label": str,
        "tls_id": str,
        "celltype_l2": str,
    },
):
    one_chunk["tls_id"] = one_chunk["tls_id"].astype(str)
    one_chunk = one_chunk[one_chunk["tls_id"].isin(cn_tls_ids)].copy()

    if one_chunk.shape[0] > 0:
        cn_cells_list.append(one_chunk)

cn_cells = pd.concat(cn_cells_list, ignore_index=True)
del cn_cells_list

cn_cells["cellid"] = cn_cells["cellid"].astype(str)
cn_cells["slide"] = cn_cells["slide"].astype(str)
cn_cells["fov"] = cn_cells["fov"].astype(str)
cn_cells["tls_id"] = cn_cells["tls_id"].astype(str)
cn_cells["celltype_l2"] = cn_cells["celltype_l2"].astype(str)
cn_cells["centroid.0"] = pd.to_numeric(cn_cells["centroid.0"], errors="coerce")
cn_cells["centroid.1"] = pd.to_numeric(cn_cells["centroid.1"], errors="coerce")

cn_cells_before_coord_filter = cn_cells.shape[0]
cn_cells = cn_cells.dropna(subset=["centroid.0", "centroid.1"]).copy()

cn_cells = cn_cells.merge(
    cn_tls_assignment[[
        "tls_id",
        "slide_fov_key",
        "fov_region",
        "tls_subtype_rank02",
    ]],
    on="tls_id",
    how="inner",
)

cn_cell_type_order = [
    x for x in cn_cell_type_order
    if x in cn_cells["celltype_l2"].unique()
]
cn_cell_type_to_code = {
    one_type: i
    for i, one_type in enumerate(cn_cell_type_order)
}
cn_cells["celltype_code"] = cn_cells["celltype_l2"].map(cn_cell_type_to_code)
cn_cells = cn_cells.dropna(subset=["celltype_code"]).copy()
cn_cells["celltype_code"] = cn_cells["celltype_code"].astype(int)

cn_cell_inclusion = pd.DataFrame({
    "item": [
        "rank02_TLS",
        "cells_before_coordinate_filter",
        "cells_after_coordinate_filter_and_TLS_join",
        "TLS_after_cell_join",
        "radius_pixels",
        "cell_types_used",
    ],
    "value": [
        int(cn_tls_assignment.shape[0]),
        int(cn_cells_before_coord_filter),
        int(cn_cells.shape[0]),
        int(cn_cells["tls_id"].nunique()),
        int(CN_RADIUS),
        ";".join(cn_cell_type_order),
    ],
})
cn_cell_inclusion.to_csv(
    os.path.join(cn_save_dir, "cell_inclusion.csv"),
    index=False,
    encoding="utf-8-sig",
)

print("Cells used for CN:", cn_cells.shape)
print("TLS used for CN:", cn_cells["tls_id"].nunique())
print(cn_cells["tls_subtype_rank02"].value_counts())
print(cn_cells["celltype_l2"].value_counts())


# %%
# 10.2 Calculate radius-based neighbor composition per cell

cn_neighborhood_records = []
cn_cells_no_neighbor = []

cn_tls_sizes = (
    cn_cells
    .groupby("tls_id")
    .size()
    .sort_values(ascending=False)
)

for cn_i, (one_tls_id, one_tls_cells) in enumerate(
    cn_cells.groupby("tls_id", sort=False),
    start=1,
):
    if cn_i == 1 or cn_i % 50 == 0 or cn_i == cn_cells["tls_id"].nunique():
        print(
            "CN neighbor calculation:",
            cn_i,
            "/",
            cn_cells["tls_id"].nunique(),
            one_tls_id,
            one_tls_cells.shape[0],
        )

    one_tls_cells = one_tls_cells.reset_index(drop=True).copy()
    one_coords = one_tls_cells[["centroid.0", "centroid.1"]].values.astype(float)
    one_type_codes = one_tls_cells["celltype_code"].values.astype(int)

    one_tree = cKDTree(one_coords)
    one_neighbors = one_tree.query_ball_point(one_coords, r=CN_RADIUS)

    one_counts = np.zeros(
        (one_tls_cells.shape[0], len(cn_cell_type_order)),
        dtype=np.int32,
    )
    one_neighbor_total = np.zeros(one_tls_cells.shape[0], dtype=np.int32)

    for row_i, neighbor_idx in enumerate(one_neighbors):
        if len(neighbor_idx) <= 1:
            continue

        neighbor_idx = np.asarray(neighbor_idx, dtype=np.int32)
        neighbor_idx = neighbor_idx[neighbor_idx != row_i]

        if neighbor_idx.shape[0] == 0:
            continue

        this_count = np.bincount(
            one_type_codes[neighbor_idx],
            minlength=len(cn_cell_type_order),
        )
        one_counts[row_i, :] = this_count
        one_neighbor_total[row_i] = int(this_count.sum())

    one_has_neighbor = one_neighbor_total > 0

    if one_has_neighbor.sum() < one_tls_cells.shape[0]:
        cn_cells_no_neighbor.append(
            one_tls_cells.loc[
                ~one_has_neighbor,
                ["cellid", "tls_id", "celltype_l2", "tls_subtype_rank02"],
            ].copy()
        )

    if one_has_neighbor.sum() == 0:
        continue

    one_prop = one_counts[one_has_neighbor, :] / one_neighbor_total[one_has_neighbor, None]
    one_prop_df = pd.DataFrame(
        one_prop,
        columns=cn_cell_type_order,
    )
    one_prop_df.insert(0, "cellid", one_tls_cells.loc[one_has_neighbor, "cellid"].values)
    one_prop_df.insert(1, "tls_id", one_tls_id)
    one_prop_df.insert(2, "slide", one_tls_cells.loc[one_has_neighbor, "slide"].values)
    one_prop_df.insert(3, "fov", one_tls_cells.loc[one_has_neighbor, "fov"].values)
    one_prop_df.insert(
        4,
        "celltype_l2",
        one_tls_cells.loc[one_has_neighbor, "celltype_l2"].values,
    )
    one_prop_df.insert(
        5,
        "tls_subtype_rank02",
        one_tls_cells.loc[one_has_neighbor, "tls_subtype_rank02"].values,
    )
    one_prop_df.insert(6, "n_neighbors", one_neighbor_total[one_has_neighbor])

    cn_neighborhood_records.append(one_prop_df)

cn_neighbor_prop = pd.concat(cn_neighborhood_records, ignore_index=True)
del cn_neighborhood_records

if len(cn_cells_no_neighbor) > 0:
    cn_cells_no_neighbor = pd.concat(cn_cells_no_neighbor, ignore_index=True)
else:
    cn_cells_no_neighbor = pd.DataFrame(
        columns=["cellid", "tls_id", "celltype_l2", "tls_subtype_rank02"]
    )

cn_neighbor_prop_path = os.path.join(cn_save_dir, "neighbor_prop.csv.gz")
cn_neighbor_prop.to_csv(
    cn_neighbor_prop_path,
    index=False,
    encoding="utf-8-sig",
    compression="gzip",
)
cn_cells_no_neighbor.to_csv(
    os.path.join(cn_save_dir, "cells_no_neighbor.csv"),
    index=False,
    encoding="utf-8-sig",
)

print("Saved neighbor composition:", cn_neighbor_prop_path)
print("Cells with neighbors:", cn_neighbor_prop.shape)
print("Cells without neighbors:", cn_cells_no_neighbor.shape)


# %%
# 10.3 Fit the fixed nine-CN model

cn_x = cn_neighbor_prop[cn_cell_type_order].astype(float).values
cn_model = MiniBatchKMeans(
    n_clusters=CN_FINAL_K,
    random_state=CN_RANDOM_STATE,
    n_init=20,
    batch_size=CN_BATCH_SIZE,
)
cn_final_raw = cn_model.fit_predict(cn_x).astype(np.int16)

cn_final_assignment = cn_neighbor_prop[[
    "cellid", "tls_id", "slide", "fov", "celltype_l2",
    "tls_subtype_rank02", "n_neighbors",
]].copy()
cn_final_assignment["CN"] = np.array(["CN_{}".format(x + 1) for x in cn_final_raw])
cn_final_assignment.to_csv(
    os.path.join(cn_save_dir, "cell_CN.csv.gz"),
    index=False,
    encoding="utf-8-sig",
    compression="gzip",
)


# %%
# 10.4 Summarize CN composition and TLS-level abundance after fitting

cn_final_count = (
    cn_final_assignment["CN"]
    .value_counts()
    .sort_index()
    .reset_index()
)
cn_final_count.columns = ["CN", "n_cells"]

cn_final_mean = (
    pd.concat(
        [
            cn_final_assignment[["CN"]].reset_index(drop=True),
            cn_neighbor_prop[cn_cell_type_order].reset_index(drop=True),
        ],
        axis=1,
    )
    .groupby("CN")[cn_cell_type_order]
    .mean()
    .reindex(cn_final_count["CN"].tolist())
)

cn_final_mean.to_csv(
    os.path.join(cn_save_dir, "cn_mean.csv"),
    encoding="utf-8-sig",
)
cn_final_count.to_csv(
    os.path.join(cn_save_dir, "cn_count.csv"),
    index=False,
    encoding="utf-8-sig",
)

cn_tls_count = (
    cn_final_assignment
    .groupby(["tls_id", "tls_subtype_rank02", "CN"])
    .size()
    .reset_index(name="n_cells")
)
cn_tls_wide = (
    cn_tls_count
    .pivot_table(
        index=["tls_id", "tls_subtype_rank02"],
        columns="CN",
        values="n_cells",
        aggfunc="sum",
        fill_value=0,
    )
    .reset_index()
)
cn_final_cn_cols = cn_final_count["CN"].tolist()
for one_cn_name in cn_final_cn_cols:
    if one_cn_name not in cn_tls_wide.columns:
        cn_tls_wide[one_cn_name] = 0

cn_tls_total = cn_tls_wide[cn_final_cn_cols].sum(axis=1)
cn_tls_prop = cn_tls_wide.copy()
cn_tls_prop[cn_final_cn_cols] = cn_tls_prop[cn_final_cn_cols].div(
    cn_tls_total.replace(0, np.nan),
    axis=0,
)
cn_tls_prop[cn_final_cn_cols] = cn_tls_prop[cn_final_cn_cols].fillna(0)

cn_tls_prop.to_csv(
    os.path.join(cn_save_dir, "tls_CN_prop.csv"),
    index=False,
    encoding="utf-8-sig",
)

cn_tls_prop_long = cn_tls_prop.melt(
    id_vars=["tls_id", "tls_subtype_rank02"],
    value_vars=cn_final_cn_cols,
    var_name="CN",
    value_name="CN_prop",
)
cn_tls_prop_long.to_csv(
    os.path.join(cn_save_dir, "tls_CN_prop_long.csv"),
    index=False,
    encoding="utf-8-sig",
)

cn_final_stat_rows = []
for one_cn_name in cn_final_cn_cols:
    b_values = cn_tls_prop.loc[
        cn_tls_prop["tls_subtype_rank02"] == "B_rich",
        one_cn_name,
    ].astype(float)
    t_values = cn_tls_prop.loc[
        cn_tls_prop["tls_subtype_rank02"] == "T_Chol_Neu_rich",
        one_cn_name,
    ].astype(float)

    try:
        one_p = mannwhitneyu(
            b_values,
            t_values,
            alternative="two-sided",
        ).pvalue
    except Exception:
        one_p = np.nan

    cn_final_stat_rows.append({
        "CN": one_cn_name,
        "mean_prop_B_rich": b_values.mean(),
        "mean_prop_T_Chol_Neu_rich": t_values.mean(),
        "delta_TCN_minus_Brich": t_values.mean() - b_values.mean(),
        "p_value": one_p,
    })

cn_final_stats = pd.DataFrame(cn_final_stat_rows)

cn_final_p_values = cn_final_stats["p_value"].astype(float).values
cn_final_p_adj = np.full(cn_final_p_values.shape, np.nan)
cn_final_valid_p = np.isfinite(cn_final_p_values)
if cn_final_valid_p.sum() > 0:
    cn_final_p_work = cn_final_p_values[cn_final_valid_p]
    cn_final_order = np.argsort(cn_final_p_work)
    cn_final_ranked = cn_final_p_work[cn_final_order]
    cn_final_adj_ranked = cn_final_ranked * len(cn_final_ranked) / (
        np.arange(len(cn_final_ranked)) + 1
    )
    cn_final_adj_ranked = np.minimum.accumulate(cn_final_adj_ranked[::-1])[::-1]
    cn_final_adj_ranked = np.clip(cn_final_adj_ranked, 0, 1)
    cn_final_adj_work = np.empty_like(cn_final_adj_ranked)
    cn_final_adj_work[cn_final_order] = cn_final_adj_ranked
    cn_final_p_adj[cn_final_valid_p] = cn_final_adj_work

cn_final_stats["p_adj_BH"] = cn_final_p_adj
cn_final_stats["p_star"] = "ns"
cn_final_stats.loc[cn_final_stats["p_value"] < 0.05, "p_star"] = "*"
cn_final_stats.loc[cn_final_stats["p_value"] < 0.01, "p_star"] = "**"
cn_final_stats.loc[cn_final_stats["p_value"] < 0.001, "p_star"] = "***"
cn_final_stats.to_csv(
    os.path.join(cn_save_dir, "cn_stats.csv"),
    index=False,
    encoding="utf-8-sig",
)

cn_run_log = pd.DataFrame({
    "item": [
        "analysis", "neighbor_radius_pixels", "final_K", "K_mode",
        "n_init", "batch_size", "inertia",
        "n_cells_with_neighbors", "n_cells_without_neighbors", "n_TLS", "random_state",
    ],
    "value": [
        "fixed_K_cell_neighborhood", CN_RADIUS, CN_FINAL_K, "fixed_archived_model",
        cn_model.n_init, CN_BATCH_SIZE, cn_model.inertia_,
        int(cn_neighbor_prop.shape[0]), int(cn_cells_no_neighbor.shape[0]),
        int(cn_neighbor_prop["tls_id"].nunique()), CN_RANDOM_STATE,
    ],
})
cn_run_log.to_csv(
    os.path.join(cn_save_dir, "run_log.csv"),
    index=False,
    encoding="utf-8-sig",
)


# %%
# 10.5 CN heatmaps

cn_plot_order = sorted(cn_final_mean.index, key=lambda name: int(name.split("_")[1]))
cn_final_mean_plot = cn_final_mean.loc[cn_plot_order, cn_cell_type_order] * 100
cn_y_labels = cn_plot_order

plt.figure(figsize=(11, max(4.5, len(cn_plot_order) * 0.48)))
sns.heatmap(
    cn_final_mean_plot,
    cmap="Reds",
    linewidths=0.25,
    linecolor="white",
    annot=True,
    fmt=".1f",
    cbar_kws={"label": "Mean neighbor composition (%)"},
)
plt.yticks(
    ticks=np.arange(len(cn_y_labels)) + 0.5,
    labels=cn_y_labels,
    rotation=0,
)
plt.title("Final CN mean neighbor composition | K={}".format(CN_FINAL_K))
plt.xlabel("Neighbor cell type")
plt.ylabel("Cell neighborhood")
plt.tight_layout()
plt.savefig(
    os.path.join(cn_save_dir, "cn_mean_heatmap.png"),
    dpi=300,
    bbox_inches="tight",
)
plt.savefig(
    os.path.join(cn_save_dir, "cn_mean_heatmap.pdf"),
    bbox_inches="tight",
)
plt.close()

cn_background = cn_neighbor_prop[cn_cell_type_order].mean(axis=0)
cn_celltype_oe = cn_final_mean.loc[cn_plot_order, cn_cell_type_order].div(
    cn_background.replace(0, np.nan),
    axis=1,
)
cn_celltype_oe = cn_celltype_oe.replace([np.inf, -np.inf], np.nan)
cn_celltype_log2oe = np.log2(cn_celltype_oe + 1e-9)

cn_celltype_oe.to_csv(
    os.path.join(cn_save_dir, "cn_celltype_OE.csv"),
    encoding="utf-8-sig",
)
cn_celltype_log2oe.to_csv(
    os.path.join(cn_save_dir, "cn_celltype_log2OE.csv"),
    encoding="utf-8-sig",
)

plt.figure(figsize=(11, max(4.5, len(cn_plot_order) * 0.48)))
sns.heatmap(
    cn_celltype_log2oe,
    cmap="RdBu_r",
    center=0,
    vmin=-2,
    vmax=2,
    linewidths=0.25,
    linecolor="white",
    annot=True,
    fmt=".1f",
    cbar_kws={"label": "log2 observed / expected"},
)
plt.yticks(
    ticks=np.arange(len(cn_y_labels)) + 0.5,
    labels=cn_y_labels,
    rotation=0,
)
plt.title("Cell-type enrichment within each CN | K={}".format(CN_FINAL_K))
plt.xlabel("Neighbor cell type")
plt.ylabel("Cell neighborhood")
plt.tight_layout()
plt.savefig(
    os.path.join(cn_save_dir, "cn_celltype_enrichment_heatmap.png"),
    dpi=300,
    bbox_inches="tight",
)
plt.savefig(
    os.path.join(cn_save_dir, "cn_celltype_enrichment_heatmap.pdf"),
    bbox_inches="tight",
)
plt.close()

cn_group_mean = (
    cn_tls_prop
    .groupby("tls_subtype_rank02")[cn_final_cn_cols]
    .mean()
    .reindex(["B_rich", "T_Chol_Neu_rich"])
)
cn_expected = cn_tls_prop[cn_final_cn_cols].mean(axis=0)
cn_roe = cn_group_mean.div(cn_expected.replace(0, np.nan), axis=1)
cn_roe = cn_roe[cn_plot_order]
cn_roe.to_csv(
    os.path.join(cn_save_dir, "cn_ROE.csv"),
    encoding="utf-8-sig",
)

plt.figure(figsize=(max(6, len(cn_plot_order) * 0.62), 3.2))
sns.heatmap(
    cn_roe,
    cmap="Reds",
    annot=True,
    fmt=".2f",
    linewidths=0.25,
    linecolor="white",
    cbar_kws={"label": "ROE from TLS-level mean proportion"},
)
plt.title("Final CN enrichment by TLS subtype | K={}".format(CN_FINAL_K))
plt.xlabel("Cell neighborhood")
plt.ylabel("TLS subtype")
plt.tight_layout()
plt.savefig(
    os.path.join(cn_save_dir, "cn_ROE_heatmap.png"),
    dpi=300,
    bbox_inches="tight",
)
plt.savefig(
    os.path.join(cn_save_dir, "cn_ROE_heatmap.pdf"),
    bbox_inches="tight",
)
plt.close()


# %%
# 10.6 CN boxplots by TLS subtype
#
# Unit: TLS. Each point is one TLS and y-axis is the fraction of cells assigned
# to a given CN within that TLS.

cn_box_data = cn_tls_prop_long.copy()
cn_box_data["CN"] = pd.Categorical(
    cn_box_data["CN"],
    categories=cn_plot_order,
    ordered=True,
)
cn_box_data = cn_box_data.sort_values(["CN", "tls_subtype_rank02"]).copy()

cn_box_colors = {
    "B_rich": "#e8556a",
    "T_Chol_Neu_rich": "#4a7baa",
}
cn_box_subtype_order = ["B_rich", "T_Chol_Neu_rich"]
cn_box_subtype_labels = ["B_rich", "T_Chol\nNeu_rich"]
cn_rng_box = np.random.default_rng(CN_RANDOM_STATE)

cn_n_cols = 4
cn_n_rows = int(np.ceil(len(cn_plot_order) / cn_n_cols))
fig, axes = plt.subplots(
    cn_n_rows,
    cn_n_cols,
    figsize=(13.2, max(7, cn_n_rows * 2.5)),
    sharex=False,
)
axes = np.array(axes).reshape(-1)

cn_final_stats_index = cn_final_stats.set_index("CN")

for ax_i, one_cn_name in enumerate(cn_plot_order):
    ax = axes[ax_i]
    cn_values_for_plot = []
    for one_subtype in cn_box_subtype_order:
        cn_values_for_plot.append(
            cn_tls_prop.loc[
                cn_tls_prop["tls_subtype_rank02"] == one_subtype,
                one_cn_name,
            ].astype(float).values
        )

    bp = ax.boxplot(
        cn_values_for_plot,
        positions=[0, 1],
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "black", "linewidth": 1.2},
        whiskerprops={"color": "#333333", "linewidth": 0.9},
        capprops={"color": "#333333", "linewidth": 0.9},
    )

    for box_i, one_box in enumerate(bp["boxes"]):
        one_box.set_facecolor(cn_box_colors[cn_box_subtype_order[box_i]])
        one_box.set_edgecolor("#333333")
        one_box.set_alpha(0.55)

    for subtype_i, one_subtype in enumerate(cn_box_subtype_order):
        y_values = cn_values_for_plot[subtype_i]
        x_values = subtype_i + cn_rng_box.normal(0, 0.055, size=len(y_values))
        ax.scatter(
            x_values,
            y_values,
            s=7,
            alpha=0.20,
            linewidths=0,
            color=cn_box_colors[one_subtype],
        )

    one_stat = cn_final_stats_index.loc[one_cn_name]
    q_text = f"{one_stat['p_adj_BH']:.2g}" if pd.notna(one_stat["p_adj_BH"]) else "NA"

    ax.set_title(
        "{} {}\nBH q={}".format(one_cn_name, one_stat["p_star"], q_text),
        fontsize=9,
    )
    ax.set_xticks([0, 1])
    ax.set_xticklabels(cn_box_subtype_labels, fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(axis="y", color="#dddddd", linewidth=0.55, alpha=0.75)
    ax.set_axisbelow(True)

    if ax_i % cn_n_cols == 0:
        ax.set_ylabel("CN fraction per TLS", fontsize=9)

for empty_i in range(len(cn_plot_order), len(axes)):
    axes[empty_i].axis("off")

fig.suptitle(
    "TLS-level CN abundance by rank02 TLS subtype | K={}".format(CN_FINAL_K),
    fontsize=13,
    y=0.995,
)
plt.tight_layout(rect=[0, 0, 1, 0.975])
plt.savefig(
    os.path.join(cn_save_dir, "cn_boxplots.png"),
    dpi=300,
    bbox_inches="tight",
)
plt.savefig(
    os.path.join(cn_save_dir, "cn_boxplots.pdf"),
    bbox_inches="tight",
)
plt.close()

print("Saved CN analysis:", cn_save_dir)
print("Fixed CN K:", CN_FINAL_K)
