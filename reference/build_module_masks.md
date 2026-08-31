# Build Module Masks

Builds boolean masks for module-active cells from module summaries and a
pseudotime mask.

## Usage

``` r
build_module_masks(
  module_summ,
  psd_mask,
  scale_threshold = 0,
  top_cells_percent = 100
)
```

## Arguments

- module_summ:

  A named list of module summary vectors.

- psd_mask:

  A logical vector indicating cells to keep.

- scale_threshold:

  Numeric threshold used to select active cells.

- top_cells_percent:

  Percentage of top cells kept per module.

## Value

A logical matrix of module masks.
