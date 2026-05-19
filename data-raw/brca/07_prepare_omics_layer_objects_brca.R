## 07_prepare_omics_layer_objects_brca.R
##
## Purpose:
##   Create a lightweight provenance object for the TCGA-BRCA multi-omics
##   layer examples used in OmicsKit vignettes.
##
## Scope:
##   This script does NOT regenerate heavy genome-track or circos inputs.
##   It documents the real pre-rendered figures generated externally from
##   TCGA-BRCA data for sample TCGA-A1-A0SH-01.
##
## Inputs expected in the repository:
##   - vignettes/figures/BRCA1_multiomics_track_8.jpg
##   - vignettes/figures/TCGA_A1_A0SH_genomewide_circos.jpg
##   - vignettes/figures/TCGA_A1_A0SH_chr17_circos.jpg
##
## Output:
##   - data/brca_omics_layer_tracks.rda

options(stringsAsFactors = FALSE)

project_file <- function(...) {
  if (requireNamespace("here", quietly = TRUE)) {
    here::here(...)
  } else {
    file.path(...)
  }
}

require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(
      "Required file is missing for BRCA omics-layer provenance object: ",
      label, "\nExpected location: ",
      normalizePath(path, mustWork = FALSE),
      call. = FALSE
    )
  }
}

figures_dir <- project_file("vignettes", "figures")

figure_files <- list(
  genome_track_brca1 = file.path(figures_dir, "BRCA1_multiomics_track_8.jpg"),
  circos_genomewide = file.path(figures_dir, "TCGA_A1_A0SH_genomewide_circos.jpg"),
  circos_chr17 = file.path(figures_dir, "TCGA_A1_A0SH_chr17_circos.jpg")
)

for (nm in names(figure_files)) {
  require_file(figure_files[[nm]], nm)
}

brca_omics_layer_tracks <- list(
  dataset = "TCGA-BRCA",
  sample_id = "TCGA-A1-A0SH-01",
  genome_build = "hg38",
  selected_region = list(
    label = "BRCA1 +/- 500 kb",
    chromosome = "chr17",
    approximate_region_hg38 = "chr17:42.9Mb-43.6Mb",
    focal_gene = "BRCA1"
  ),
  reason_sample_selected = paste(
    "TCGA-A1-A0SH-01 was selected because it was available across the",
    "multi-omics layers used for the example: RNA-seq, CNV, DNA methylation",
    "450k, and somatic mutations. The sample also carries a BRCA1 nonsense",
    "mutation annotated as BRCA1_p.Q934* in the pre-rendered genome-track",
    "example."
  ),
  omics_layers = list(
    RNA = list(
      source = "UCSC Xena TCGA.BRCA.sampleMap_HiSeqV2",
      representation = "RNA-seq expression signal"
    ),
    CNV = list(
      source = "UCSC Xena TCGA.BRCA.sampleMap_SNP6_nocnv_genomicSegment",
      representation = "copy-number segment signal"
    ),
    methylation = list(
      source = "UCSC Xena TCGA.BRCA.sampleMap_HumanMethylation450",
      representation = "DNA methylation beta-value signal"
    ),
    mutations = list(
      source = "MC3 BRCA somatic mutation callset",
      representation = "functional mutation annotations in the BRCA1 region"
    )
  ),
  processing_summary = list(
    original_build = "hg19 where applicable",
    target_build = "hg38",
    liftover = "hg19ToHg38.over.chain",
    note = paste(
      "Heavy processing steps such as liftOver, bigWig/BED generation,",
      "genome-wide binning, and circos plotting were performed outside the",
      "vignette. The package stores only a compact provenance object and",
      "pre-rendered real figures."
    )
  ),
  figures = list(
    genome_track_brca1 = list(
      file = "vignettes/figures/BRCA1_multiomics_track_8.jpg",
      description = paste(
        "Multi-omics genome-track view of the BRCA1 region for sample",
        "TCGA-A1-A0SH-01, including gene annotations, RNA-seq, CNV,",
        "methylation, and mutation layers."
      ),
      intended_function = "nice_GenomeTrack"
    ),
    circos_genomewide = list(
      file = "vignettes/figures/TCGA_A1_A0SH_genomewide_circos.jpg",
      description = paste(
        "Genome-wide circos visualization for sample TCGA-A1-A0SH-01",
        "showing RNA, CNV, methylation, and mutation burden tracks."
      ),
      intended_function = "nice_circos"
    ),
    circos_chr17 = list(
      file = "vignettes/figures/TCGA_A1_A0SH_chr17_circos.jpg",
      description = paste(
        "Chromosome 17 circos visualization for sample TCGA-A1-A0SH-01,",
        "focused on the genomic context relevant to BRCA1."
      ),
      intended_function = "nice_circos"
    )
  ),
  vignette_usage = list(
    recommended_chunk_strategy = paste(
      "Show calls to nice_GenomeTrack() and nice_circos() with eval = FALSE,",
      "then display the corresponding pre-rendered figures with",
      "knitr::include_graphics()."
    ),
    rationale = paste(
      "This avoids forcing users or continuous-integration jobs to download",
      "large TCGA/Xena source files or regenerate heavy genome-track objects."
    )
  ),
  created = as.character(Sys.time())
)

if (!requireNamespace("usethis", quietly = TRUE)) {
  stop(
    "Package \"usethis\" is required to save brca_omics_layer_tracks.",
    call. = FALSE
  )
}

usethis::use_data(
  brca_omics_layer_tracks,
  compress = "xz",
  overwrite = TRUE
)

message("Saved brca_omics_layer_tracks to data/brca_omics_layer_tracks.rda")
message("Documented pre-rendered BRCA omics-layer figures:")
message("  - ", figure_files$genome_track_brca1)
message("  - ", figure_files$circos_genomewide)
message("  - ", figure_files$circos_chr17)
