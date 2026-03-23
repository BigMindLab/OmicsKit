# Power analysis for RNA-seq differential expression with optional plotting

Power analysis for RNA-seq differential expression with optional
plotting

## Usage

``` r
power_analysis(
  effect_size = 1,
  dispersion = 0.1,
  n_genes = 20000,
  prop_de = 0.05,
  alpha = 0.05,
  power_target = 0.8,
  max_n = 20,
  plot = TRUE
)
```

## Arguments

- effect_size:

  Numeric. Log2 fold-change to detect.

- dispersion:

  Numeric. Biological variance (σ²).

- n_genes:

  Integer. Total number of genes tested.

- prop_de:

  Numeric. Proportion of truly DE genes (0–1).

- alpha:

  Numeric. Desired FDR (type I error rate).

- power_target:

  Numeric. Desired statistical power (1 – β).

- max_n:

  Integer. Maximum sample size per group to explore.

- plot:

  Logical. If TRUE, draws the power curve; if FALSE, skips plotting.

## Value

A named list. If `plot = TRUE`, contains three elements:

- `$min_sample_size`: Integer. Minimum sample size per group to reach
  `power_target`.

- `$power_table`: A data frame with columns `SampleSize` and `Power`.

- `$plot`: A ggplot2 object of the power curve.

If `plot = FALSE`, returns only `$min_sample_size` and `$power_table`.

## Examples

``` r
# Basic power analysis with default parameters
result <- power_analysis(
  effect_size  = 1,
  dispersion   = 0.1,
  n_genes      = 20000,
  prop_de      = 0.05,
  alpha        = 0.05,
  power_target = 0.8,
  max_n        = 20
)

# Minimum sample size to reach 80% power
result$min_sample_size
#> [1] 7

# Full power table
head(result$power_table)
#>   SampleSize      Power
#> 1          2 0.06234486
#> 2          3 0.20477735
#> 3          4 0.41078547
#> 4          5 0.61880322
#> 5          6 0.78217643
#> 6          7 0.88846759

# Higher effect size requires fewer samples
power_analysis(effect_size = 2, dispersion = 0.1, plot = FALSE)$min_sample_size
#> [1] 2

# See plot
#power_analysis$plot
```
