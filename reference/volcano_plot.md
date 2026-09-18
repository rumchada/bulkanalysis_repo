# Volcano Plot Visualization

Uses `ggplot2`, `ggrepel`, and `tidyverse` to construct a volcano plot
from differential expression results.

## Usage

``` r
volcano_plot(
  diffexp_df,
  log2fc_thresh = 0,
  p.val_adj_thresh = 0.05,
  lower_xlim = -5,
  upper_xlim = 5,
  step = 1
)
```

## Arguments

- diffexp_df:

  A tibble or data frame containing differential expression results.
  Must contain, at minimum, `p.val_adj` and `log2FC` columns.

- log2fc_thresh:

  Numeric threshold for the absolute log2 fold-change used to identify
  differentially expressed genes for display.

- p.val_adj_thresh:

  Numeric threshold for the adjusted p-value used to identify
  statistically significant differentially expressed genes.

- lower_xlim:

  Numeric value specifying the lower limit of the x-axis. Should be less
  than 0. For example, `-5`.

- upper_xlim:

  Numeric value specifying the upper limit of the x-axis. Should be
  greater than 0. For example, `5`.

- step:

  Numeric value specifying the interval between x-axis ticks.

## Value

A list containing:

- `plot`: A ggplot object containing the volcano plot.

- `data`: A data frame containing the differential expression results
  used to construct the plot.
