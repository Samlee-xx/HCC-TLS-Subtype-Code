# HCC-TLS Subtype Code

This repository provides analysis scripts and selected processed data for identifying, characterizing, and subtyping tertiary lymphoid structures (TLSs) in hepatocellular carcinoma (HCC). The analyses span spatial transcriptomic, CODEX, Visium HD, CosMx, and immunotherapy-related datasets.

## Analysis scope

The figure-oriented scripts cover:

- TLS detection, quality control, cellular composition analysis, spatial clustering, and subtype classification;
- CODEX-based DBSCAN identification of TLS regions and assessment of candidate regions;
- spatial organization analyses, including Delaunay-based TLS layering, cellular neighborhoods, cell–cell proximity, Ripley's cross-K, and spatial co-occurrence;
- spatial transcriptomic TLS signature analysis, feature importance, and cell-type deconvolution using cell2location, with input preparation and analysis of externally generated CIBERSORTx results;
- Visium HD TLS profiling and integration of CODEX and Visium HD cellular composition and functional-program profiles;
- tumor-proximity axis analysis, gene-expression correlations, and functional enrichment;
- CosMx cell annotation, TLS identification, chemokine and receptor expression, and local spatial relationships;
- prognostic and immunotherapy-related analyses, including TLS signature scoring, survival analysis, logistic regression, and neutrophil composition analysis;
- serial-section TLS subtype prediction and marker-panel visualization;
- TLSProfiler feature selection, model fitting and evaluation, and analyses of HookNet-TLS outputs.

## Repository structure

```text
HCC-TLS-Subtype-Code/
├── code/       # R and Python analysis scripts and workflow notes
├── data/       # Selected processed inputs, annotations, model files, and reference outputs
└── README.md
```

Scripts are named by the associated figure and analysis. The `data/` directory contains:

| Directory | Contents |
| --- | --- |
| `01-ST_TLS` | Spatial transcriptomic TLS objects, metadata, and deconvolution-related tables |
| `02-CODEX_Cohort2` | Cohort 2 cell annotations, TLS inputs, and cellular neighborhood data |
| `03-CODEX_Cohort3` | Cohort 3 cell annotations and TLS inputs |
| `04-Visium_HD` | Visium HD TLS profiles, cell annotations, and spatial-axis inputs |
| `05-cross_K` | Inputs for Ripley's cross-K analysis |
| `06-ROE` | Cell-type and neighborhood inputs for observed-to-expected enrichment analysis |
| `07-HookNet_TLS` | Processed HookNet-TLS analysis inputs |
| `08-CosMx` | Processed CosMx TLS data |
| `09-TLSProfiler` | Model definitions, processed features, calibration/assessment inputs, and workflow metadata |
| `10-serial` | Serial-section TLS data, model files, marker panels, and reference outputs |
| `11-CODEX_Cohort4` | Cohort 4 TLS and clinical-analysis inputs |

## Downloading the data

Large data files are tracked with [Git Large File Storage (Git LFS)](https://git-lfs.com/). Install Git LFS before cloning, then run:

```bash
git lfs install
git clone https://github.com/Samlee-xx/HCC-TLS-Subtype-Code.git
cd HCC-TLS-Subtype-Code
git lfs pull
```

Use Git with LFS to obtain the complete data files. Without LFS, some files may be downloaded as small pointer files rather than their actual contents.

## Running the analyses

The scripts are individual, figure-oriented analyses rather than a single end-to-end pipeline. They use both the processed inputs included here and external inputs referenced in the scripts. Raw sequencing data, complete source images, and all upstream preprocessing steps are not bundled with this repository.

Run commands from the repository root and check the input/output configuration at the beginning of each script. Some scripts retain study-specific absolute paths or defaults that assume `data/` is beside the script. Supply explicit paths where command-line arguments are supported, or update the path variables for your local environment.

For example, the following commands specify input and output directories for the corresponding scripts:

```bash
python code/Figure1_DBSCAN_TLS.py --data-dir data --output-dir result/Figure1_DBSCAN
Rscript code/Figure7_02_TLSProfiler.R data/09-TLSProfiler result/Figure7_TLSProfiler/02-final_model
```

Install the R/Python packages imported by the analysis you intend to run. Workflow metadata, model parameters, and software-version information are provided within relevant data subdirectories where available. CIBERSORTx processing and HookNet-TLS inference additionally require their respective external tools and inputs.

Saved reference outputs support comparison with newly generated results; they should not be interpreted as independent reruns. Final subtype annotations are provided for downstream analyses, and their presence alone does not reproduce all upstream annotation decisions. Model-selection, calibration, and assessment steps should be interpreted according to the accompanying scripts and metadata.

## Contact

For code-related questions, please contact the first author, Shuai Li, at [21110880009@m.fudan.edu.cn](mailto:21110880009@m.fudan.edu.cn).
