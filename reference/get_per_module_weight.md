# Compute Per-Module Gene Weights

Computes normalized intra-module connectivity weights for genes. For
each gene, the sum of weights of the edges incident to it are summed.
These weights are normalised by the maximum weight of the module.

## Usage

``` r
get_per_module_weight(gene_adj_matrix, gene_modules)
```

## Arguments

- gene_adj_matrix:

  A gene adjacency matrix.

- gene_modules:

  A named vector or list defining gene-module assignments.

## Value

A named numeric vector of module-normalized weights per gene.
