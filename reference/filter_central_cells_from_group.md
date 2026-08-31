# Find central points of a group

Given a population of cells and a UMAP embedding, this function attempts
to identify a group of cells that are located in the centre of this
group. For this, the function calculates the geometric median of the
UMAP and then calculates the nearest points to this point.

## Usage

``` r
filter_central_cells_from_group(cell_group, umap_embedding, n_points = 5)
```

## Arguments

- cell_group:

  A vector with the names of the cells to be analysed.

- umap_embedding:

  A matrix with the UMAP embedding of the cells. The matrix should have
  defined the rownames as the cell names.

- n_points:

  The number of points to be selected as the central points.

## Value

A vector with the names of the cells that are considered central.
