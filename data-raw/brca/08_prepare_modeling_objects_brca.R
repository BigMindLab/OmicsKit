## 08_prepare_modeling_objects_brca.R
##
## Purpose:
##   Placeholder contract for clinical-only TCGA-BRCA modeling objects.
##
## Inputs:
##   - Real TCGA-BRCA clinical data from prepared metadata objects or
##     inst/extdata/brca/intermediate/.
##
## Outputs:
##   - Package-ready clinical-only modeling object(s) for data/.
##   - Optional intermediate model input tables in inst/extdata/brca/intermediate/.
##
## Expected file locations:
##   - inst/extdata/brca/intermediate/
##   - data/
##
## TODO:
##   - Define clinical outcome, survival time, event, and predictor columns.
##   - Use clinical data only.
##   - Do not include omics_score, PAM50, RPPA clusters, methylation clusters,
##     CN clusters, or multiomics cluster variables.
##   - Add clear stop() messages for missing clinical inputs.

stop("TODO: prepare real TCGA-BRCA clinical-only modeling objects.")
