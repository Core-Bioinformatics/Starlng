# Build Module Transition Adjacency

Infers module connectivity by mapping modules to trajectory nodes and
pruning triangles based on module similarity.

## Usage

``` r
get_module_transitions(
  trajectory_object,
  closest_module,
  start_node = NULL,
  similarity_values = NULL
)
```

## Arguments

- trajectory_object:

  A trajectory object containing graph and node data.

- closest_module:

  Named vector mapping modules to closest nodes.

- start_node:

  Optional node used to initialize graph traversal. If NULL, a node with
  degree 1 will be selected as the starting point.

- similarity_values:

  Optional matrix used to break triangles. The values will be used to
  calculate the distance between the triangle edges.

## Value

A symmetric module adjacency matrix.
