# Find detectable genes across comparisons.

This function identifies genes with measurable expression levels across
samples. Detectable genes must meet two conditions: the baseMean and
their mean normalized counts in the phenotypes of interest must be
greater than a set threshold. It returns a list of detectable genes and
the comparisons in which they can be found.

## Usage

``` r
detect_filter(
  norm.counts,
  df.BvsA,
  df.CvsA = NULL,
  df.DvsA = NULL,
  cutoffs = c(50, 50, 0),
  samples.baseline,
  samples.condition1,
  samples.condition2 = NULL,
  samples.condition3 = NULL
)
```

## Arguments

- norm.counts:

  Data frame of the normalized counts with Ensembl IDs as rows and
  Sample IDs as columns.

- df.BvsA:

  Data frame comparing the first condition to the baseline.

- df.CvsA:

  Data frame comparing the second condition to the baseline (optional).

- df.DvsA:

  Data frame comparing the third condition to the baseline (optional).

- cutoffs:

  Vector containing threshold values for baseMean, mean normalized
  counts and Log2 Fold Change; respectively. Default: c(50, 50, 0).

- samples.baseline:

  Vector of Sample IDs or indexes corresponding to the baseline.

- samples.condition1:

  Vector of Sample IDs or indexes corresponding to the first condition.

- samples.condition2:

  Vector of Sample IDs or indexes corresponding to the second condition
  (optional).

- samples.condition3:

  Vector of Sample IDs or indexes corresponding to the third condition
  (optional).

## Value

A named list. Always contains:

- `$Comparison1`: Data frame of detectable genes from `df.BvsA`.

- `$DetectGenes`: Character vector of unique detectable gene IDs across
  all comparisons.

If `df.CvsA` is provided, also contains `$Comparison2`. If `df.DvsA` is
provided, also contains `$Comparison3`.

## See also

[`nice_VSB()`](https://danielgarbozo.github.io/OmicsKit/reference/nice_VSB.md)
to plot expression of detected genes;
[norm_counts](https://danielgarbozo.github.io/OmicsKit/reference/norm_counts.md)
for an example normalized counts matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
data(norm_counts)
data(deseq2_results)
data(sampledata)

# detect_filter requires an "ensembl" column in the results data frame
res <- deseq2_results
colnames(res)[colnames(res) == "gene_id"] <- "ensembl"
rownames(res) <- res$ensembl

# Get sample IDs per group
samples_normal <- sampledata$patient_id[sampledata$sample_type == "normal"]
samples_tumor  <- sampledata$patient_id[sampledata$sample_type == "tumor"]

detected <- detect_filter(
  norm.counts        = as.data.frame(norm_counts),
  df.BvsA            = res,
  samples.baseline   = samples_normal,
  samples.condition1 = samples_tumor,
  cutoffs            = c(50, 50, 0)
)

# Number of detectable genes
length(detected$DetectGenes)

# Subset results
head(detected$Comparison1)
} # }
```
