# Plot UMAP with Discrete or Continuous Coloring

This function plots a UMAP embedding with points colored by either
discrete or continuous cell information. The function provides options
for sorting the cells based on the coloring variable, adjusting point
size and alpha, and customizing the color scheme for both discrete and
continuous coloring.

## Usage

``` r
plot_umap(
  umap_embedding,
  cell_info,
  mtd_name,
  cell_sort_order = c("lowest", "highest", "default"),
  scale_values = FALSE,
  cell_size = 0.3,
  cell_alpha = 0.8,
  legend_text_size = 10,
  legend_key_size = 0.5,
  axis_text_size = 10,
  show_labels = TRUE,
  label_size = 10,
  colourbar_width = 50,
  discrete_colors = NULL,
  continuous_colors = NULL
)
```

## Arguments

- umap_embedding:

  A matrix or data frame containing the UMAP coordinates of the cells.
  It should have two columns corresponding to UMAP1 and UMAP2.

- cell_info:

  A vector containing the information to color the cells by. It can be
  either a factor/character vector for discrete coloring or a numeric
  vector for continuous coloring. The length of this vector should match
  the number of rows in `umap_embedding`.

- mtd_name:

  A string specifying the name of the metadata variable to be used in
  the legend.

- cell_sort_order:

  A character vector specifying the order of cells based on the
  `cell_info` variable. For discrete coloring, it should contain the
  unique values of `cell_info` in the desired order. For continuous
  coloring, it can be either "lowest", "highest" or "default" to sort
  cells by the `cell_info` values in ascending, descending or no
  particular order, respectively.

- scale_values:

  A logical indicating whether to scale the `cell_info` values between 0
  and 1 for continuous coloring. Defaults to FALSE.

- cell_size:

  A numeric value specifying the size of the points in the plot.
  Defaults to 0.3.

- cell_alpha:

  A numeric value between 0 and 1 specifying the transparency of the
  points in the plot. Defaults to 0.8.

- legend_text_size:

  A numeric value specifying the size of the text in the legend.
  Defaults to 10.

- legend_key_size:

  A numeric value specifying the size of the legend keys for discrete
  coloring. Defaults to 0.5.

- axis_text_size:

  A numeric value specifying the size of the axis text in the plot.
  Defaults to 10.

- show_labels:

  A logical indicating whether to show labels for discrete groups.
  Defaults to TRUE.

- label_size:

  A numeric value specifying the size of the labels for discrete groups.
  Defaults to 10.

- colourbar_width:

  A numeric value specifying the width of the colorbar for continuous
  coloring. Defaults to 50 points. If set to 0, the colorbar will be
  hidden.

- discrete_colors:

  A list of color vectors for discrete coloring, where the names of the
  list correspond to the number of unique values in `cell_info`. If
  NULL, the default ggplot2 colors will be used.

- continuous_colors:

  A vector of colors for continuous coloring. If NULL, the default
  ggplot2 colors will be used.

## Value

A ggplot object representing the UMAP plot with the specified coloring
and customization options.
