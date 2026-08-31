# Calculate Average UMAP Distance

Calculates an average distance for a selected set of cells in UMAP
space, either from the centroid or as the mean pairwise distance.

## Usage

``` r
calculate_umap_average_distance(
  umap_df,
  selected_cells = NULL,
  centroid = TRUE
)
```

## Arguments

- umap_df:

  A matrix or data frame containing UMAP coordinates.

- selected_cells:

  Optional vector of cell indices or names to include.

- centroid:

  Logical indicating whether to use the centroid-based distance
  calculation. If FALSE, pairwise mean distance is returned.

## Value

A numeric vector of distances, or a single numeric value when
`centroid = FALSE`.
