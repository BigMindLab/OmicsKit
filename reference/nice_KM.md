# Function to generate Kaplan Meier survival plots for a given binary gene variable

This function fits Kaplan Meier survival curves stratified by status
(e.g., 'No' vs 'Yes') for a specified gene column in a data frame. It
allows:

- Automatic placement of the p-value at the midpoint of the time range
  (if coord = NULL).

- Automatic handling of cases where only one category ('No' or 'Yes') is
  present, returning an 'empty' placeholder plot with a warning message.

## Usage

``` r
nice_KM(
  data,
  gene,
  time_var,
  event_var,
  coord = NULL,
  title_prefix = "Mut ",
  colors = c("#1F8FFF", "#ED4D4D"),
  conf_int = TRUE,
  risk_table = FALSE,
  legend_pos = NULL,
  title_size = 24,
  axis_text = c(x = 18, y = 18),
  axis_title_size = c(x = 20, y = 20),
  legend_text = c(title = 18, labels = 16),
  pval_size = 6,
  returnData = FALSE
)
```

## Arguments

- data:

  Data.frame. Must contain columns specified by `gene`, `time_var`, and
  `event_var`.

- gene:

  String. Name of the column in `data` indicating mutation status
  (values 'No'/'Yes').

- time_var:

  String. Name of the column for survival time.

- event_var:

  String. Name of the column for event indicator (0/1).

- coord:

  Numeric or NULL. X-coordinate at which to display the log-rank
  p-value. If NULL, placed at midpoint of time range.

- title_prefix:

  String. Prefix for the legend title. Default `"Mut "`.

- colors:

  Character vector of length 2. Two colors for the plot lines. Default
  `c("#1F8FFF", "#ED4D4D")`.

- conf_int:

  Logical. Whether to show confidence intervals. Default `TRUE`.

- risk_table:

  Logical. Whether to include risk-table below the plot. Default
  `FALSE`.

- legend_pos:

  Numeric vector of length 2 or NULL. X/Y inside-plot coordinates for
  legend placement. Default `NULL`.

- title_size:

  Numeric. Font size for the plot title.

- axis_text:

  Named numeric vector. Sizes for axis tick labels.

- axis_title_size:

  Named numeric vector. Sizes for axis titles.

- legend_text:

  Named numeric vector. Sizes for legend title and labels.

- pval_size:

  Numeric. Font size for the p-value annotation.

- returnData:

  Logical. If `TRUE`, returns a list with `km_fit` and `plot`; if
  `FALSE`, returns only the `ggplot` object. Default `FALSE`.

## Value

A `ggplot` object (or a list with `km_fit` and `plot` if
`returnData = TRUE`).

A ggplot2 object if `returnData = FALSE` (default). If
`returnData = TRUE`, a named list with two elements:

- `$km_fit`: The
  [`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)
  object.

- `$plot`: The ggplot2 survival curve.

## References

Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from
incomplete observations. *Journal of the American Statistical
Association*, 53(282), 457–481.
[doi:10.1080/01621459.1958.10501452](https://doi.org/10.1080/01621459.1958.10501452)
