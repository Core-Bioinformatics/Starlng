# Add Hub-Gene Counts to TF Statistics

Adds, for each TF, the number of associated genes that are hub genes.

## Usage

``` r
add_tf_hub_stats(existing_stats, hub_genes)
```

## Arguments

- existing_stats:

  A TF statistics data frame.

- hub_genes:

  A data frame containing hub genes and their module labels. The name of
  the columns should be "gene" and "module".

## Value

A TF statistics data frame with an added n_hub_genes column when
possible.
