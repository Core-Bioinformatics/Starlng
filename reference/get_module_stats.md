# Compute Module Statistics

Builds per-cell module statistics using pseudotime and UMAP distances.

## Usage

``` r
get_module_stats(module_summ, module_mask, psd_value, umap_df, centroid = TRUE)
```

## Arguments

- module_summ:

  A named list of module summary vectors.

- module_mask:

  A logical matrix of module masks.

- psd_value:

  A numeric pseudotime vector.

- umap_df:

  A matrix or data frame with UMAP coordinates.

- centroid:

  Logical indicating whether UMAP distance should be computed from the
  centroid.

## Value

A data frame with module statistics, or NULL if no cells are selected.
