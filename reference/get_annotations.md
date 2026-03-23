# Get annotations from Ensembl.

This function annotates a column of transcripts or gene IDs (ENSEMBL)
with information of the Biomart. If transcript IDs are provided, they
are also annotated with information of the genes to which they belong.

## Usage

``` r
get_annotations(
  ensembl_ids,
  mode = "genes",
  filename = "gene_annotations",
  version = "Current",
  format = "csv"
)
```

## Arguments

- ensembl_ids:

  The column of transcripts to be used as input.

- mode:

  To specify the IDs provided, between "transcripts" or "genes". Default
  = genes.

- filename:

  The name of the output file, which is table. Default =
  gene_annotations.

- version:

  This function can use the versions 102, 103, and 112 of Ensembl.
  Default = "Current".

- format:

  The output is saved in .csv or .xlsx formats. Default = csv.

## Value

A data frame with one row per input ID and the following columns:
`geneID`, `symbol`, `biotype`, `chromosome`, `gene_start`, `gene_end`,
`gene_length`, `description`. For `mode = "transcripts"`, an additional
`transcriptID` column is included. The data frame is also saved to disk
as a `.csv` or `.xlsx` file (see `filename` and `format`).

## Details

The Gene information added include:

- Gene ENSEMBL ID, HGNC Symbol, Description, Biotype and Chromosome.

- Gene start, end and length

## Note

Requires an active internet connection to query the Ensembl BioMart.
`gene_length` is computed as `gene_end - gene_start + 1` (genomic
length). For TPM calculation with
[`tpm()`](https://danielgarbozo.github.io/OmicsKit/reference/tpm.md),
this is an approximation, use transcript-level lengths for higher
accuracy.

## See also

[`add_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/add_annotations.md)
to join annotations to a counts matrix;
[`tpm()`](https://danielgarbozo.github.io/OmicsKit/reference/tpm.md)
which requires gene lengths from this function.

## Examples

``` r
if (FALSE) { # \dontrun{
# Annotate genes from Normalized counts (requires internet connection)
data(norm_counts)

# Requires a reference table with a "geneID" column.
# Use get_annotations() to generate it:
annotations <- get_annotations(
  ensembl_ids = rownames(norm_counts),
  mode        = "genes"
)

head(annotations)

# Use with add_annotations()
norm_counts_annot <- add_annotations(
  object    = norm_counts,
  reference = annotations,
  variables = c("symbol", "biotype")
)
} # }
```
