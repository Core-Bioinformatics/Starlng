# Community detection wrapper

This function is a simplified version of the `clustering_pipeline`. It
starts with an already-generated adjacency matrix and applied the leiden
algorithm to identify the clusters. The details of how the algorithm is
run are the same as in the `clustering_pipeline` function.

## Usage

``` r
community_detection_master(
  adj_object,
  resolutions,
  number_iterations,
  seeds,
  graph_type = "snn",
  memory_log_file = NULL
)
```

## Arguments

- adj_object:

  The adjacency matrix to be used. The matrix should be square and
  should have defined rownames.

- resolutions:

  Check the `clustering_pipeline` documentation.

- number_iterations:

  Check the `clustering_pipeline` documentation.

- seeds:

  Check the `clustering_pipeline` documentation.

- graph_type:

  Check the `clustering_pipeline` documentation.

- memory_log_file:

  Check the `clustering_pipeline` documentation.

## Value

A list of lists containing the clusters found for each combination of
number of neighbours and quality functions.
