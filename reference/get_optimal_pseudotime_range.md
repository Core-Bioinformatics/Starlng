# Compute an Optimal Pseudotime Range

Chooses a root node based on branch-length balance over the trajectory
diameter and returns the associated pseudotime ordering.

## Usage

``` r
get_optimal_pseudotime_range(monocle_object)
```

## Arguments

- monocle_object:

  A Monocle3 object with a learned trajectory graph.

## Value

A list with selected start node and pseudotime values.
