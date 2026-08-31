# Plot TF-Gene Network with ggraph

Plots a TF-gene igraph object. The TFs and the module gene hubs are
highlighted using different shapes and sized of the points. The edges
between TFs and genes are coloured differently from the edges between
TFs.

## Usage

``` r
plot_module_tfs_ggraph(
  module_tf_g,
  tf_colour = "black",
  module_colours = NULL,
  tf_shape = 17,
  hub_shape = 23,
  gene_shape = 16,
  node_size = 4,
  edge_width = 0.5,
  edge_colour = "#bbbbbb",
  label_size = 3,
  axis_text_size = 8,
  exclude_non_hub_genes = FALSE
)
```

## Arguments

- module_tf_g:

  An igraph object produced by get_tf_gene_network.

- tf_colour:

  Colour used for TF nodes and TF-TF edges.

- module_colours:

  Named vector of module colours.

- tf_shape:

  Shape used for TF nodes.

- hub_shape:

  Shape used for hub-gene nodes.

- gene_shape:

  Shape used for non-hub gene nodes.

- node_size:

  Base node-size scaling factor.

- edge_width:

  Base edge width.

- edge_colour:

  Colour used for TF-gene edges.

- label_size:

  Text size of node labels.

- axis_text_size:

  Base text size for title and legends.

- exclude_non_hub_genes:

  Logical indicating whether non-hub genes are included in the plotting.

## Value

A ggplot object, or invisible(NULL) when plotting dependencies are
unavailable.
