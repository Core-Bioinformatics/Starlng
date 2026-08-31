# Learn the graph of the Monocle object

This function follows the logic of the `learn_graph` function from the
monocle3 package. The main difference consists in providing more
flexibility to the user, such as the number of nodes per log10 cells and
the metadata used as partition.

## Usage

``` r
custom_learn_graph(
  mon_obj,
  nodes_per_log10_cells = 30,
  learn_graph_controls = list(eps = 1e-05, maxiter = 100),
  use_partition = FALSE,
  use_closed_loops = FALSE,
  verbose = FALSE,
  metadata_for_nodes = NULL
)
```

## Arguments

- mon_obj:

  The Monocle object.

- nodes_per_log10_cells:

  The number of trajectory nodes created per log10 cells. Defaults to
  30.

- learn_graph_controls:

  A list of control parameters, as defined in the `learn_graph` function
  from Monocle. Defaults to list(eps = 1e-5, maxiter = 100).

- use_partition:

  A logical value indicating if the partition should be used for
  learning the graph. Defaults to FALSE.

- use_closed_loops:

  A logical value indicating if circular paths can be formed in the
  trajectory graph. Defaults to FALSE.

- verbose:

  Parameter that is passed to the `learn_graph` function and prints the
  progress.

- metadata_for_nodes:

  The metadata column that should be used to update the partition of the
  Monocle object. This is used when `use_partition` is set to TRUE.
  Defaults to NULL, meaning the partition is not updated.

## Value

The Monocle object with the learned graph.
