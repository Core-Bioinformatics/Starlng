# Build a Trajectory Helper Object

Extracts graph nodes, edges and projected cell-to-node assignments from
a Monocle3 object for downstream plotting and analysis.

## Usage

``` r
get_trajectory_object(monocle_object, reduction_name = "UMAP")
```

## Arguments

- monocle_object:

  A Monocle3 object.

- reduction_name:

  Name of the trajectory reduction to use.

## Value

A list with graph, node positions, edge table and closest vertices.
