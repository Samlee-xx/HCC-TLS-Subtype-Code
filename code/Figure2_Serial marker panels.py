###### Figure 2M: TLS1 across four serial sections ######

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tifffile


###### inputs and display settings ######

repo = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description="Serial TLS fluorescence panels")
parser.add_argument("--data-dir", type=Path, default=repo / "data/10-serial")
parser.add_argument("--output-dir", type=Path, default=repo / "result/Figure2_serial")
args = parser.parse_args()
args.output_dir.mkdir(parents=True, exist_ok=True)

panels = pd.read_csv(args.data_dir / "images/Figure2M_panels.csv").sort_values("display_order")
predictions = pd.read_csv(args.output_dir / "serial_TLS_predictions_TLS.csv")
panels = panels.merge(
    predictions[["TLS", "section_number", "fov", "pred_subtype_logistic"]],
    on=["TLS", "section_number", "fov"], how="left", validate="one_to_one",
)
if len(panels) != 4 or panels["pred_subtype_logistic"].isna().any():
    raise ValueError("Expected four TLS1 panels with matching subtype predictions.")
with open(args.data_dir / "images/display_parameters.json", encoding="utf-8") as handle:
    display = json.load(handle)

plt.rcParams.update({"font.family": "Arial", "font.size": 10, "pdf.fonttype": 42})
figure = plt.figure(figsize=(13.4, 5.0), facecolor="white")
grid = figure.add_gridspec(3, 8, height_ratios=[2.1, 1.0, 0.7], hspace=0.04, wspace=0.035)
figure.subplots_adjust(left=0.13, right=0.995, top=0.85, bottom=0.07)


###### fixed channel scaling and additive fluorescence overlays ######

for column, panel in enumerate(panels.itertuples(index=False)):
    image_dir = args.data_dir / "images" / panel.image_directory
    with open(image_dir / "region.json", encoding="utf-8") as handle:
        region = json.load(handle)
    if region["fov"] != panel.fov:
        raise ValueError("Image region and panel metadata disagree.")

    channels = {}
    for marker, settings in display["channels"].items():
        raw = tifffile.imread(image_dir / f"{marker}.tif").astype(np.float32)
        if raw.shape != (region["height"], region["width"]):
            raise ValueError(f"Unexpected image shape: {image_dir / marker}")
        intensity = np.clip((raw - settings["min"]) / (settings["max"] - settings["min"]), 0, 1)
        channels[marker] = intensity[..., None] * np.array(settings["rgb"])
    composite = np.clip(sum(channels.values()), 0, 1)

    main = figure.add_subplot(grid[0, 2 * column:2 * column + 2])
    main.imshow(composite, interpolation="nearest")
    main.set_axis_off()
    position = main.get_position(original=True)
    figure.text((position.x0 + position.x1) / 2, position.y1 + 0.015,
                f"Section {panel.section_number}", ha="center", fontsize=12)
    bar_pixels = display["scale_bar_um"] / region["pixel_width_um"]
    x_end = region["width"] * 0.96
    y_bar = region["height"] * 0.94
    main.plot([x_end - bar_pixels, x_end], [y_bar, y_bar], color="white", linewidth=2)
    if column == 0:
        main.text(-0.035, 0.5, "TLS 1", transform=main.transAxes, rotation=90,
                  ha="right", va="center", fontsize=13)

    for offset, marker in enumerate(["CD21", "Ki67"]):
        axis = figure.add_subplot(grid[1, 2 * column + offset])
        axis.imshow(channels[marker], interpolation="nearest")
        axis.set_axis_off()
        axis.text(0.025, 0.04, marker, color=display["channels"][marker]["rgb"],
                  transform=axis.transAxes, ha="left", va="bottom", fontsize=10)

    caption = figure.add_subplot(grid[2, 2 * column:2 * column + 2])
    caption.set_axis_off()
    caption.text(0.5, 0.70, panel.pred_subtype_logistic, ha="center", va="center")
    caption.text(0.5, 0.25, panel.figure_maturation_stage, ha="center", va="center")
    if column == 0:
        caption.text(-0.5, 0.70, "TLS subtype", ha="left", va="center")
        caption.text(-0.5, 0.25, "Maturation\nstage", ha="left", va="center")

legend = figure.add_axes([0.39, 0.94, 0.28, 0.05], facecolor="black")
legend.set_xticks([])
legend.set_yticks([])
for position, marker in zip([0.125, 0.375, 0.625, 0.875], ["CD3e", "CD21", "Ki67", "CD20"]):
    settings = display["channels"][marker]
    legend.text(position, 0.5, settings["label"], color=settings["rgb"], ha="center", va="center", fontsize=12)
figure.text(0.015, 0.96, "M", fontsize=18, weight="bold")
figure.text(0.99, 0.025, f"Scale bars: {display['scale_bar_um']} µm", ha="right", fontsize=9)

for extension in ["pdf", "png"]:
    figure.savefig(args.output_dir / f"Figure2M_serial_TLS1_markers.{extension}", dpi=300)
plt.close(figure)
panels.to_csv(args.output_dir / "Figure2M_panel_metadata.csv", index=False)
