# Plot TF Bubble Plot

Creates a bubble plot of top TFs per module.

## Usage

``` r
plot_module_tfs_bubbleplot(
  tf_stats,
  n_top = 10,
  cap_value_intersection = NULL,
  cap_value_hub = NULL,
  font_size = 10,
  point_size_range = c(3, 10)
)
```

## Arguments

- tf_stats:

  A TF statistics data frame.

- n_top:

  Number of top TFs to keep per module.

- cap_value_intersection:

  Optional cap for n_genes.

- cap_value_hub:

  Optional cap for n_hub_genes.

- font_size:

  Base font size for plot labels and legends.

- point_size_range:

  Numeric vector of length 2 defining point size range.

## Value

A ggplot object.
