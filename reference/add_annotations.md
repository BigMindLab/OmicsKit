# A function to add annotations to a table of gene counts.

A function to add annotations to a table of gene counts.

## Usage

``` r
add_annotations(object, reference, variables = NULL, data_frame = FALSE)
```

## Arguments

- object:

  A table of gene counts (rows: genes, columns: samples, rownames are
  ENSEMBL gene IDs).

- reference:

  A reference table with the annotations including a column named
  "geneID".

- variables:

  Character vector of columns in `reference` to add. If NULL (default),
  all columns except geneID are used.

- data_frame:

  Logical; if TRUE, coerce `object` to a data.frame first. Default:
  FALSE.

## Value

The input `object` as a data frame with additional columns from
`reference` joined by Ensembl gene ID. A `geneID` column is added
containing the row names of the original object.

## See also

[`get_annotations()`](https://danielgarbozo.github.io/OmicsKit/reference/get_annotations.md)
to generate the `reference` table;
[norm_counts](https://danielgarbozo.github.io/OmicsKit/reference/norm_counts.md)
for an example input matrix.

## Examples

``` r
if (FALSE) { # \dontrun{
data(norm_counts)

# Requires a reference table with a "geneID" column.
# Use get_annotations() to generate it:
annotations <- get_annotations(
  ensembl_ids = rownames(norm_counts),
  mode        = "genes"
)

# Add gene symbol and biotype columns to the counts matrix
norm_counts_annot <- add_annotations(
  object    = norm_counts,
  reference = annotations,
  variables = c("symbol", "biotype")
)

# Inspect result
head(norm_counts_annot[, c("geneID", "symbol", "biotype")])

# Add all annotation columns (variables = NULL uses everything)
norm_counts_full <- add_annotations(
  object    = norm_counts,
  reference = annotations
)
} # }
```
