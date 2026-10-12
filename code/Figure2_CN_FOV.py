"""CN maps using segmentation masks and cell assignments."""
import argparse
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--segmentation-root', type=Path, required=True,
                    help='Contains <slide>/deepcell_output/<fov>_whole_cell.tiff')
parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parent / 'data')
parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parent / 'result/Figure2_CN/FOV_667')
parser.add_argument('--slide')
parser.add_argument('--fov')
parser.add_argument('--historical-scope', action='store_true',
                    help='Include all 696 historical regions, matching the original FOV script')
args = parser.parse_args()
data_dir = args.data_dir
table_dir = args.output_dir / '01-tables'
figure_dir = args.output_dir / '02-figures'

import os
from pathlib import Path
from tempfile import TemporaryDirectory

os.environ["MPLBACKEND"] = "Agg"

import matplotlib

matplotlib.use("Agg", force=True)

import numpy as np
import pandas as pd
import tifffile

import ark.settings as settings
from ark.utils.plot_utils import cohort_cluster_plot


# ---------- read final CN assignment ----------
cn = pd.read_csv(
    data_dir / "02-CODEX_Cohort2/CN/cn_cell_assignment_k9_final.csv",
    usecols=["cellid", "slide", "fov", "CN"],
    dtype={"cellid": str, "slide": str, "fov": str, "CN": str},
)

if not args.historical_scope:
    final_cells = pd.read_csv(data_dir / '06-ROE/CODEX_CN/cn_tls_layer_cell_table_final.csv', usecols=['cellid'])
    cn = cn[cn['cellid'].isin(final_cells['cellid'])].copy()
if args.slide:
    cn = cn[cn['slide'] == args.slide].copy()
if args.fov:
    cn = cn[cn['fov'] == args.fov].copy()
if cn.empty:
    raise ValueError('No cells for the requested scope, slide and FOV.')

cn["CN"] = cn["CN"].str.replace("CN_", "CN", regex=False)
cn = cn[cn["CN"].isin([f"CN{x}" for x in range(1, 10)])].copy()
cn["label"] = cn["cellid"].str.rsplit("__", n=1).str[-1].astype(int)
cn = cn.rename(columns={"CN": "cell_meta_cluster"})


# ---------- CN colors from alluvial plot ----------
cn_colors = pd.DataFrame(
    {
        "cell_meta_cluster": [f"CN{x}" for x in range(1, 10)],
        "color": [
            "#ef8a8c",
            "#2c7fb8",
            "#359a35",
            "#e31a1c",
            "#b2d58a",
            "#a6cee3",
            "#f6bd6b",
            "#f47c13",
            "#6a3d9a",
        ],
    }
)

Path(table_dir).mkdir(
    parents=True,
    exist_ok=True,
)
cn_colors.to_csv(
    table_dir / 'CN_color_map_final.csv',
    index=False,
)


# ---------- FOV summary ----------
summary = (
    cn.groupby(["slide", "fov", "cell_meta_cluster"])
    .size()
    .reset_index(name="n_cells")
    .pivot_table(
        index=["slide", "fov"],
        columns="cell_meta_cluster",
        values="n_cells",
        fill_value=0,
    )
    .reset_index()
)
summary.to_csv(
    table_dir / 'CN_FOV_plot_summary_final.csv',
    index=False,
)


# ---------- plot each slide ----------
Path(figure_dir).mkdir(
    parents=True,
    exist_ok=True,
)

for slide in sorted(cn["slide"].unique()):
    seg_dir = Path(args.segmentation_root) / slide / "deepcell_output"
    if not seg_dir.exists():
        print("Skip missing segmentation:", slide)
        continue

    fovs = sorted(cn.loc[cn["slide"] == slide, "fov"].unique())
    fovs = [x for x in fovs if (seg_dir / f"{x}_whole_cell.tiff").exists()]
    if len(fovs) == 0:
        print("Skip slide without matching FOV masks:", slide)
        continue

    one_slide = cn[(cn["slide"] == slide) & (cn["fov"].isin(fovs))].copy()
    save_dir = Path(
        figure_dir
    ) / slide
    save_dir.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory(prefix=f"cn_fov_{slide}_") as temp_dir:
        temp_dir = Path(temp_dir)

        for fov in fovs:
            raw_seg = seg_dir / f"{fov}_whole_cell.tiff"
            out_seg = temp_dir / f"{fov}_whole_cell.tiff"
            mask = np.squeeze(tifffile.imread(raw_seg))
            tifffile.imwrite(out_seg, mask)

        cohort_cluster_plot(
            fovs=fovs,
            seg_dir=temp_dir,
            save_dir=save_dir,
            cell_data=one_slide[["fov", "label", "cell_meta_cluster"]],
            erode=True,
            fov_col=settings.FOV_ID,
            label_col=settings.CELL_LABEL,
            cluster_col=settings.CELL_TYPE,
            seg_suffix="_whole_cell.tiff",
            cmap=cn_colors,
            display_fig=False,
            fig_file_type="png",
            figsize=(10, 10),
            dpi=300,
        )

    print("Saved CN FOV plots:", slide, len(fovs))

print("Done.")
