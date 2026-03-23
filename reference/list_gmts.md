# Read GMT files from a directory into a named gene set list

Scans a directory for `.gmt` files, parses them, and returns a single
named list where each element is a character vector of gene symbols for
one gene set. The output is ready to be passed directly to
[`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md).

## Usage

``` r
list_gmts(dir)
```

## Arguments

- dir:

  Character. Path to the directory containing one or more `.gmt` files.
  The function searches the directory non-recursively.

## Value

A named list where each element is a character vector of gene symbols
for one gene set. Names correspond to gene set names as defined in
column 1 of the GMT files. If the same gene set name appears in multiple
files, the last occurrence overwrites the earlier one.

## Details

**GMT format:** each row contains the gene set name in column 1, an
optional description in column 2, and gene symbols from column 3 onward.
Empty fields are automatically removed. This is the standard format used
by MSigDB (Molecular Signatures Database) and other gene set annotation
resources.

## See also

[`geneset_similarity()`](https://danielgarbozo.github.io/OmicsKit/reference/geneset_similarity.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Read all GMT files from a directory
geneset_list <- list_gmts("path/to/gmt_files/")

# Inspect output
length(geneset_list)              # number of gene sets
names(geneset_list)[1:5]          # first five gene set names
geneset_list[["KEGG_APOPTOSIS"]]  # genes in a specific set

# Pass directly to geneset_similarity
jac <- geneset_similarity(geneset_list, results_df, fdr_th = 0.05)
} # }
```
