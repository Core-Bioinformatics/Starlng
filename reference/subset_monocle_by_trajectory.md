# Subset the Monocle object by a trajectory

This function subsets the Monocle object based on a subgraph of the
inferrred trajectory. The subgraph is defined using start and end
points, defined by cells.

## Usage

``` r
subset_monocle_by_trajectory(mon_obj, start_cells = NULL, end_cells = NULL)
```

## Arguments

- mon_obj:

  A monocle object.

- start_cells:

  A list of cells names that define the start of the selected
  trajectory. If NULL, no subsetting is performed. Defaults to NULL.

- end_cells:

  A list of cells names that define the end of the selected trajectory.
  If NULL, no subsetting is performed. Defaults to NULL.

## Value

The Monocle object subsetted by the trajectory.
