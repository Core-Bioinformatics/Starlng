# Order the cells by pseudotime

This function follows the logic of the `order_cells` function from the
monocle3 package. The main difference consists in allowing the user to
define both the start and the end points of the trajectory. In this
case, the monocle object is subsetted.

## Usage

``` r
custom_pseudotime_ordering(
  monocle_object,
  start_cells = NULL,
  end_cells = NULL
)
```

## Arguments

- monocle_object:

  A monocle object.

- start_cells:

  A list of cells names that define the start of the ordering. If NULL,
  the end_cells must be provided. Defaults to NULL.

- end_cells:

  A list of cells names that define the end of the ordering. If NULL,
  the start_cells must be provided. If end_cells is not NULL and
  start_cells is not NULL, subsetting is performed. If end_cells is not
  NULL and start_cells is NULL, the ordering will be reversed. Defaults
  to NULL.

## Value

A named vector with the pseudotime values. In the case of subsetting,
the cells that are not part of the inferred subtrajectory will be
assigned the NA value.
