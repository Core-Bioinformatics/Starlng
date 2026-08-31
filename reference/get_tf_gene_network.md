# Build TF-Gene Network Graph

Builds an igraph object representing TF-TF and TF-gene relationships.
The graph is built based on the associations done from the TF enrichment
analysis. TF-TF relationships are derived from the module ordering: TFs
within the same module are linked in descending n_genes order, and
inter-module links connect the top-ranked TF from each adjacent module.

## Usage

``` r
get_tf_gene_network(
  tf_stats,
  module_adjacency = NULL,
  top_n_factors = 3,
  hub_genes = NULL
)
```

## Arguments

- tf_stats:

  A TF statistics data frame including module, tf and genes.

- module_adjacency:

  Optional module adjacency matrix.

- top_n_factors:

  Number of top TFs to retain per module.

- hub_genes:

  Optional data frame of hub genes used for annotation.

## Value

An igraph object, or invisible(NULL) when no valid TFs are available.
