# Plot Module Transition Graph

Plots module connectivity as a tree-like graph.

## Usage

``` r
plot_module_transitions(
  module_adj_matrix,
  closest_module,
  start_module = NULL,
  edge_size = 0.5,
  edge_alpha = 0.5,
  edge_colour = "black",
  node_size = 3,
  node_colours = "white",
  node_label_size = 3,
  node_label_vjust = -1
)
```

## Arguments

- module_adj_matrix:

  A module adjacency matrix.

- closest_module:

  Named vector mapping modules to closest nodes.

- start_module:

  Optional module name used as graph root.

- edge_size:

  Edge width.

- edge_alpha:

  Edge transparency.

- edge_colour:

  Edge colour.

- node_size:

  Node size.

- node_colours:

  Node fill colour.

- node_label_size:

  Node-label text size.

- node_label_vjust:

  Vertical adjustment for node labels.

## Value

A ggplot object.
