"""Prepare K9 CN tables for the 667-TLS cohort."""

from pathlib import Path

import argparse
import json
import pandas as pd

parser = argparse.ArgumentParser(description="Prepare CN tables for the final 667 manuscript TLSs.")
parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parent / "data")
parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "result/Figure2_CN")
args = parser.parse_args()
data_dir = args.data_dir
table_dir = args.output_dir / "01-tables"
table_dir.mkdir(parents=True, exist_ok=True)


Path(table_dir).mkdir(parents=True, exist_ok=True)


# ---------- read final TLS layer labels ----------
tls_layer = pd.read_csv(
    data_dir / "06-ROE/CODEX_celltype/tls_delaunay_layer_cells_final.csv",
    usecols=["cellid", "tls_id", "tls_subtype", "celltype_l2", "delaunay_layer"],
    dtype={
        "cellid": str,
        "tls_id": str,
        "tls_subtype": str,
        "celltype_l2": str,
        "delaunay_layer": str,
    },
)


# ---------- read K9 CN assignment ----------
cn = pd.read_csv(
    data_dir / "02-CODEX_Cohort2/CN/cn_cell_assignment_k9_final.csv",
    usecols=["cellid", "CN"],
    dtype={"cellid": str, "CN": str},
)

cn["CN"] = cn["CN"].str.replace("CN_", "CN", regex=False)


# ---------- merge final TLS cells with CN labels ----------
final_tls = pd.read_csv(
    data_dir / "02-CODEX_Cohort2/TLS/tls_subtype_assignment_final.csv",
    dtype={"tls_id": str},
)
if final_tls["tls_id"].duplicated().any() or len(final_tls) != 667:
    raise ValueError("Expected 667 unique final manuscript TLSs.")
if final_tls["tls_subtype"].value_counts().to_dict() != {"BF-TLS": 325, "TNC-TLS": 342}:
    raise ValueError("Final TLS subtype counts differ from the manuscript cohort.")
if tls_layer["cellid"].duplicated().any() or cn["cellid"].duplicated().any():
    raise ValueError("Duplicate cellid in a CN or layer input.")
if set(tls_layer["tls_id"]) != set(final_tls["tls_id"]):
    raise ValueError("Layer input TLS IDs differ from the final 667 TLSs.")
expected_subtype = tls_layer["tls_id"].map(final_tls.set_index("tls_id")["tls_subtype"])
if not expected_subtype.eq(tls_layer["tls_subtype"]).all():
    raise ValueError("TLS subtype differs between the layer and assignment inputs.")
actual_sizes = tls_layer.groupby("tls_id").size()
if not actual_sizes.eq(final_tls.set_index("tls_id")["n_cells"].reindex(actual_sizes.index)).all():
    raise ValueError("Per-TLS cell counts differ from the final assignment input.")
total = tls_layer.merge(cn, on="cellid", how="left", validate="one_to_one")
if not total["CN"].isin([f"CN{x}" for x in range(1, 10)]).all():
    raise ValueError("Missing or invalid CN labels among final TLS cells.")
if len(total) != 1414248:
    raise ValueError("Expected 1,414,248 final TLS cells.")
total["Boundary_Class"] = total["delaunay_layer"]
total.loc[total["Boundary_Class"] == "Periphery", "Boundary_Class"] = "Edge"

total = total[
    [
        "cellid",
        "tls_id",
        "tls_subtype",
        "celltype_l2",
        "delaunay_layer",
        "Boundary_Class",
        "CN",
    ]
].copy()

total.to_csv(
    table_dir / 'cn_tls_layer_cell_table_final.csv',
    index=False,
    encoding="utf-8-sig",
)


# ---------- CN proportion per TLS ----------
cn_order = [f"CN{x}" for x in range(1, 10)]

cn_prop = (
    total[total["CN"].isin(cn_order)]
    .groupby(["tls_id", "tls_subtype", "CN"])
    .size()
    .reset_index(name="n_cells")
)

tls_total = (
    total[total["CN"].isin(cn_order)]
    .groupby(["tls_id", "tls_subtype"])
    .size()
    .reset_index(name="n_tls_cells")
)

all_tls = tls_total[["tls_id", "tls_subtype"]].drop_duplicates()
all_cn = pd.DataFrame({"CN": cn_order})
all_tls["_tmp"] = 1
all_cn["_tmp"] = 1

cn_prop = (
    all_tls.merge(all_cn, on="_tmp")
    .drop(columns="_tmp")
    .merge(cn_prop, on=["tls_id", "tls_subtype", "CN"], how="left")
    .merge(tls_total, on=["tls_id", "tls_subtype"], how="left")
)

cn_prop["n_cells"] = cn_prop["n_cells"].fillna(0).astype(int)
cn_prop["percentage"] = cn_prop["n_cells"] / cn_prop["n_tls_cells"] * 100

cn_prop.to_csv(
    table_dir / 'cn_tls_proportion_long_final.csv',
    index=False,
    encoding="utf-8-sig",
)

cn_wide = cn_prop.pivot_table(
    index=["tls_id", "tls_subtype"],
    columns="CN",
    values="percentage",
    fill_value=0,
).reset_index()

cn_wide.to_csv(
    table_dir / 'cn_tls_proportion_wide_final.csv',
    index=False,
    encoding="utf-8-sig",
)


# ---------- CN proportion per TLS layer ----------
layer_order = ["Periphery", "Transition", "Core"]

layer_count = (
    total[total["CN"].isin(cn_order) & total["delaunay_layer"].isin(layer_order)]
    .groupby(["tls_id", "tls_subtype", "delaunay_layer", "CN"])
    .size()
    .reset_index(name="n_cells")
)

all_layer = pd.DataFrame({"delaunay_layer": layer_order})
all_layer["_tmp"] = 1

all_tls_layer = all_tls.merge(all_layer, on="_tmp").drop(columns="_tmp")
all_tls_layer["_tmp"] = 1

all_cn_final = all_cn.drop(columns="_tmp")
all_cn_final["_tmp"] = 1

all_tls_layer_cn = all_tls_layer.merge(all_cn_final, on="_tmp").drop(columns="_tmp")

cn_layer_prop = all_tls_layer_cn.merge(
    layer_count,
    on=["tls_id", "tls_subtype", "delaunay_layer", "CN"],
    how="left",
)

cn_layer_prop["n_cells"] = cn_layer_prop["n_cells"].fillna(0).astype(int)
cn_layer_prop["n_layer_cells"] = cn_layer_prop.groupby(
    ["tls_id", "tls_subtype", "delaunay_layer"]
)["n_cells"].transform("sum")

cn_layer_prop["percentage"] = (
    cn_layer_prop["n_cells"] / cn_layer_prop["n_layer_cells"] * 100
)
cn_layer_prop.loc[cn_layer_prop["n_layer_cells"] == 0, "percentage"] = pd.NA

cn_layer_prop.to_csv(
    table_dir / 'CN_layer_proportion_by_TLS_final.csv',
    index=False,
    encoding="utf-8-sig",
)


# ---------- CN celltype count ----------
cn_celltype_count = (
    total[total["CN"].isin(cn_order)]
    .groupby(["CN", "celltype_l2"])
    .size()
    .reset_index(name="n_cells")
)

cn_celltype_count.to_csv(
    table_dir / 'cn_celltype_count_final.csv',
    index=False,
    encoding="utf-8-sig",
)


print("Final cells:", total.shape[0])
print("Final TLS:", total["tls_id"].nunique())
print("CN labels:")
print(total["CN"].value_counts(dropna=False).sort_index())
print("TLS subtype:")
print(total["tls_subtype"].value_counts(dropna=False))
print("Saved CN final input tables.")

summary = {
    "historical_assignment_cells": len(cn),
    "final_TLS": int(total["tls_id"].nunique()),
    "final_cells": len(total),
    "final_TLS_subtypes": final_tls["tls_subtype"].value_counts().to_dict(),
    "cells_without_layer": int(total["delaunay_layer"].isna().sum()),
    "CN_source": "Frozen historical K9 assignment, joined by cellid",
}
(args.output_dir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
