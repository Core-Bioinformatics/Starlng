# Find Closest Trajectory Node to Module

Finds the closest trajectory node to a module using the UMAP distance
between the module centroid and the node location.

## Usage

``` r
get_closest_node_to_module(
  trajectory_object,
  cell_umap,
  module_expr,
  expression_threshold = 0,
  expression_percentile = 0,
  scale = TRUE
)
```

## Arguments

- trajectory_object:

  A trajectory object containing node positions.

- cell_umap:

  A cell embedding matrix or data frame with two dimensions.

- module_expr:

  A numeric vector with the aggregate expression of the module. If a
  logical vector is provided, it will be interpreted as the mask
  indicating the module population. Providing a list will return a named
  vector of closest nodes for each module.

- expression_threshold:

  Expression threshold used to define module-active cells. Will not be
  used if `module_expr` is already logical. Defaults to 0.

- expression_percentile:

  Optional percentile threshold for module-active cells. Will be applied
  after `expression_threshold` if provided. Defaults to 0.

- scale:

  Logical indicating whether expression values should be scaled between
  0 and 1. Will not be used if `module_expr` is already logical.
  Defaults to TRUE.

## Value

The closest node name, or a named vector of node names for list input.
