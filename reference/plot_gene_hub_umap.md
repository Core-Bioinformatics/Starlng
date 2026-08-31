# Plot Gene UMAP with Hub Highlighting

Plots genes in UMAP space with module-colored points and weighted edges.
If provided, the hub genes are highlighted with larger points and
labels.

## Usage

``` r
plot_gene_hub_umap(
  umap_df,
  filtered_gene_adj,
  edge_weight_range = c(0.1, 1),
  edge_alpha = 0.5,
  edge_colour = "#bbbbbb",
  point_size = 2,
  point_alpha = 0.85,
  module_colours = NULL,
  node_text_size = 5,
  legend_text_size = 4,
  hub_point_scale = 2,
  hub_stroke = 0.7,
  label_box_alpha = 0.85,
  label_box_stroke = 0.15
)
```

## Arguments

- umap_df:

  A gene embedding matrix or data frame.

- filtered_gene_adj:

  Output from `get_filtered_gene_adjacency`.

- edge_weight_range:

  Numeric range for edge linewidth scaling.

- edge_alpha:

  Edge transparency.

- edge_colour:

  Default edge colour.

- point_size:

  Base point size.

- point_alpha:

  Point transparency.

- module_colours:

  Optional named vector of module colours.

- node_text_size:

  Text size for hub labels.

- legend_text_size:

  Text size for legends.

- hub_point_scale:

  Size multiplier for hub genes.

- hub_stroke:

  Stroke width for hub points.

- label_box_alpha:

  Background transparency for hub labels.

- label_box_stroke:

  Border width for hub label boxes.

## Value

A ggplot object.
