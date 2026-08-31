# Compute Module Pairwise Tables

Computes pairwise overlap, union and Spearman statistics between the
populations described by the module masks and summaries.

## Usage

``` r
compute_module_pairwise_tables(module_mask, module_summ)
```

## Arguments

- module_mask:

  A logical matrix of module masks.

- module_summ:

  A named list of module summary vectors.

## Value

A list with pairwise overlap and correlation matrices.
