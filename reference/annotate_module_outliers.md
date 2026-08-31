# Annotate Module Outliers

Adds outlier labels to module summaries using module masks and
pseudotime information.

## Usage

``` r
annotate_module_outliers(
  modules_stats_summary,
  module_mask,
  psd_value,
  umap_dist_threshold = NULL
)
```

## Arguments

- modules_stats_summary:

  A data frame of module summary statistics.

- module_mask:

  A logical matrix of module masks.

- psd_value:

  A numeric pseudotime vector.

- umap_dist_threshold:

  Optional threshold for median UMAP distance.

## Value

The module summary data frame with an is_outlier column.
