# Plot a Trajectory Graph

Plots trajectory edges and optionally overlays selected node degrees.
This plot is meant to be used as an additional layer in an UMAP plot.

## Usage

``` r
plot_trajectory_graph(
  trajectory_object,
  edge_size = 0.5,
  edge_alpha = 1,
  edge_colour = "black",
  plot_nodes = c(1, 2, 3),
  node_size = 3,
  node_colours = stats::setNames(c("red", "blue", "green"), c(1, 2, 3)),
  node_label_size = 3,
  node_label_vjust = -1
)
```

## Arguments

- trajectory_object:

  A trajectory object produced by `get_trajectory_object`.

- edge_size:

  Line width for trajectory edges.

- edge_alpha:

  Transparency of trajectory edges.

- edge_colour:

  Colour of trajectory edges.

- plot_nodes:

  Which node-degree groups to display: 1 for leaves, 2 for intermediate
  nodes, 3 for branching points. Setting to 0 will hide all nodes.

- node_size:

  Point size for trajectory nodes.

- node_colours:

  Named colours for node-degree groups.

- node_label_size:

  Text size of node labels.

- node_label_vjust:

  Vertical adjustment for node labels.

## Value

A ggplot object of the trajectory.
