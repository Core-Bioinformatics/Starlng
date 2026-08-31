# Summarize Module Statistics

Summarizes per-cell module statistics into per-module summaries.

## Usage

``` r
summarise_module_stats(modules_stats, gene_modules)
```

## Arguments

- modules_stats:

  A data frame returned by get_module_stats.

- gene_modules:

  A named list mapping modules to genes.

## Value

A data frame of module-level summary statistics.
