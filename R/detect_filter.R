##########################
# Function detect_filter #
##########################

#' Find detectable genes across comparisons.
#'
#' This function identifies genes with measurable expression levels across samples.
#' Detectable genes must meet two conditions: the baseMean and their mean normalized counts in the phenotypes of interest must be greater than a set threshold.
#' It returns a list of detectable genes and the comparisons in which they can be found.
#'
#' @param norm.counts Data frame of the normalized counts with Ensembl IDs as rows and Sample IDs as columns.
#' @param df.BvsA Data frame comparing the first condition to the baseline.
#' @param df.CvsA Data frame comparing the second condition to the baseline (optional).
#' @param df.DvsA Data frame comparing the third condition to the baseline (optional).
#' @param samples.baseline Vector of Sample IDs or indexes corresponding to the baseline.
#' @param samples.condition1 Vector of Sample IDs or indexes corresponding to the first condition.
#' @param samples.condition2 Vector of Sample IDs or indexes corresponding to the second condition (optional).
#' @param samples.condition3 Vector of Sample IDs or indexes corresponding to the third condition (optional).
#' @param cutoffs Vector containing threshold values for baseMean, mean normalized counts and Log2 Fold Change; respectively. Default: c(50, 50, 0).
#' @export


detect_filter <- function(norm.counts, df.BvsA, df.CvsA = NULL, df.DvsA = NULL, cutoffs = c(50, 50, 0),
			  samples.baseline, samples.condition1, samples.condition2 = NULL, samples.condition3 = NULL)

{
  # Ensure cutoffs length is correct
  if (length(cutoffs) != 3) {
    stop("Cutoffs vector must contain three values: baseMean, mean normalized counts, and Log2 Fold Change thresholds.")
  }
  
  # Create an empty vector to store genes temporary
  genes_vector <- c()
  
  # Obtain mean normalized counts per phenotype
  norm.counts$Mean.Baseline <- rowMeans(norm.counts[, samples.baseline])
  norm.counts$Mean.Condition1 <- rowMeans(norm.counts[, samples.condition1])
  
  # Filter data frames by baseMean
  df.BvsA.det <- df.BvsA[df.BvsA$baseMean > cutoffs[1], ]
  
  if (!is.null(df.CvsA)) {
    norm.counts$Mean.Condition2 <- rowMeans(norm.counts[, samples.condition2])
    df.CvsA.det <- df.CvsA[df.CvsA$baseMean > cutoffs[1], ]
  }
  
  if (!is.null(df.DvsA)) {
    norm.counts$Mean.Condition3 <- rowMeans(norm.counts[, samples.condition3])
    df.DvsA.det <- df.DvsA[df.DvsA$baseMean > cutoffs[1], ]
  }
  
  # Get detectable genes
  for (i in 1:length(df.BvsA.det[, "ensembl"])) {
    if (df.BvsA.det[i, "log2FoldChange"] > cutoffs[3]) {
      if (norm.counts[df.BvsA.det[i, "ensembl"], "Mean.Condition1"] > cutoffs[2]) {
        genes_vector <- c(genes_vector, df.BvsA.det[i, "ensembl"])
      }
    } else if (df.BvsA.det[i, "log2FoldChange"] < -cutoffs[3]) {
      if (norm.counts[df.BvsA.det[i, "ensembl"], "Mean.Baseline"] > cutoffs[2]) {
        genes_vector <- c(genes_vector, df.BvsA.det[i, "ensembl"])
      }
    }
  }
  
  if (!is.null(df.CvsA)) {
    for (i in 1:length(df.CvsA.det[, "ensembl"])) {
      if (df.CvsA.det[i, "log2FoldChange"] > cutoffs[3]) {
        if (norm.counts[df.CvsA.det[i, "ensembl"], "Mean.Condition2"] > cutoffs[2]) {
          genes_vector <- c(genes_vector, df.CvsA.det[i, "ensembl"])
        }
      } else if (df.CvsA.det[i, "log2FoldChange"] < -cutoffs[3]) {
        if (norm.counts[df.CvsA.det[i, "ensembl"], "Mean.Baseline"] > cutoffs[2]) {
          genes_vector <- c(genes_vector, df.CvsA.det[i, "ensembl"])
        }
      }
    }
  }
  
  if (!is.null(df.DvsA)) {
    for (i in 1:length(df.DvsA.det[, "ensembl"])) {
      if (df.DvsA.det[i, "log2FoldChange"] > cutoffs[3]) {
        if (norm.counts[df.DvsA.det[i, "ensembl"], "Mean.Condition3"] > cutoffs[2]) {
          genes_vector <- c(genes_vector, df.DvsA.det[i, "ensembl"])
        }
      } else if (df.DvsA.det[i, "log2FoldChange"] < -cutoffs[3]) {
        if (norm.counts[df.DvsA.det[i, "ensembl"], "Mean.Baseline"] > cutoffs[2]) {
          genes_vector <- c(genes_vector, df.DvsA.det[i, "ensembl"])
        }
      }
    }
  }
  
  # Remove duplicates
  genes_vector <- unique(genes_vector)
  
  # Subset data frames to only include detectable genes
  df.BvsA.det <- df.BvsA.det[df.BvsA.det$ensembl %in% genes_vector, ]
  detected_genes <- list(Comparison1 = df.BvsA.det, DetectGenes = genes_vector)
  
  if (!is.null(df.CvsA)) {
    df.CvsA.det <- df.CvsA.det[df.CvsA.det$ensembl %in% genes_vector, ]
    detected_genes$Comparison2 <- df.CvsA.det
  }
  
  if (!is.null(df.DvsA)) {
    df.DvsA.det <- df.DvsA.det[df.DvsA.det$ensembl %in% genes_vector, ]
    detected_genes$Comparison3 <- df.DvsA.det
  }
  
  return(detected_genes)
}
