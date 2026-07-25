# TCGA-BRCA Data Preparation Contracts

This folder contains placeholder scripts for preparing real TCGA-BRCA example
data, prerendered figures, and package-ready objects for OmicsKit.

## Rules

- Do not download files unless explicitly requested.
- Do not create synthetic biological results.
- Do not write fake placeholder data objects.
- Use `here::here()` or project-relative paths.
- Store raw/path-based inputs under `inst/extdata/brca/`.
- Store package-ready R objects under `data/`.
- Store heavy processing scripts in `data-raw/brca/`.
- Store prerendered figures under `vignettes/figures/` and
  `vignettes/figures_PA/`.

## Expected External Input Locations

- `inst/extdata/brca/raw_xena/`: raw TCGA-BRCA Xena/GDC input files.
- `inst/extdata/brca/gsea_outputs/`: real GSEA or pathway pipeline outputs.
- `inst/extdata/brca/omics_layers/`: path-based files such as `.bw`, `.bed`,
  `.tsv`, and `.rds` for omics-layer examples.
- `inst/extdata/brca/intermediate/`: project-relative intermediate files
  created by data-raw scripts.

## Script Order

1. `00_download_xena_brca.R`: declare/download raw TCGA-BRCA source files.
2. `01_prepare_metadata_brca.R`: prepare BRCA sample and clinical metadata.
3. `02_prepare_transcriptomics_dea_brca.R`: prepare RNA DEA inputs/results.
4. `03_prepare_proteomics_dea_brca.R`: prepare proteomics DEA inputs/results.
5. `04_prepare_pathway_inputs_brca.R`: prepare ranked genes and GMT inputs.
6. `05_import_gsea_pipeline_outputs_brca.R`: import real pathway outputs.
7. `06_prepare_pathway_objects_brca.R`: create pathway package objects.
8. `07_prepare_omics_layer_objects_brca.R`: prepare omics-layer objects/files.
9. `08_prepare_modeling_objects_brca.R`: prepare clinical-only modeling data.
10. `09_render_brca_figures.R`: render real prerendered vignette figures.
11. `99_save_brca_package_data.R`: save approved package data objects.

Each script must document inputs, outputs, and expected file locations before it
is implemented.
