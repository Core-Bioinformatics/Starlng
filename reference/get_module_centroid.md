# Compute Module Centroid in Embedding Space

Computes the geometric median centroid of module-active cells in UMAP
space.

## Usage

``` r
get_module_centroid(
  module_expr,
  cell_umap,
  expression_threshold = 0,
  expression_percentile = 0,
  scale = TRUE
)
```

## Arguments

- module_expr:

  A numeric/logical vector, or a list of such vectors.

- cell_umap:

  A cell embedding matrix or data frame with at least two dimensions.

- expression_threshold:

  Expression threshold used to define active cells.

- expression_percentile:

  Optional percentile threshold for active cells.

- scale:

  Logical indicating whether expression values should be scaled between
  0 and 1.

## Value

A numeric centroid vector, or a centroid matrix for list input.
