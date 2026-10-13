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
├── data/       # Processed inputs, annotations, model files, and reference outputs
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

## Contact

For code-related questions, please contact the first author, Shuai Li, at [21110880009@m.fudan.edu.cn](mailto:21110880009@m.fudan.edu.cn).
