# Plot Module Pseudobulk Expression

Computes per-metadata pseudobulk module expression and plots a bubble
heatmap.

## Usage

``` r
plot_module_pseudobulk_expression(
  module_expression_summary,
  cell_metadata,
  mtd_order,
  pseudobulk_function = mean,
  scale = TRUE,
  cap_value = 2,
  axis_text_size = 10,
  legend_text_size = 10,
  point_range = c(0.1, 1),
  continuous_colours = NULL
)
```

## Arguments

- module_expression_summary:

  A module-by-cell matrix/data frame, or list of module vectors.

- cell_metadata:

  Metadata vector used to define pseudobulk groups.

- mtd_order:

  Desired display order for metadata groups.

- pseudobulk_function:

  Function used to aggregate expression per group.

- scale:

  Logical indicating whether pseudobulk values are z-scored per module.

- cap_value:

  Maximum absolute value used for color limits.

- axis_text_size:

  Text size for axes.

- legend_text_size:

  Text size for legends.

- point_range:

  Range for point-size scaling.

- continuous_colours:

  Optional color palette for expression values.

## Value

A ggplot object.
