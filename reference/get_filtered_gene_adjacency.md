# Filter Gene Adjacency by Module Structure

Filters gene-gene edges based on module connectivity i.e. remove genes
between genes from different modules if those modules are not neighbours
in the trajectory graph. The function provides the option to sample
non-hub genes, as well as select a top percentage of edges based on
their weight.

## Usage

``` r
get_filtered_gene_adjacency(
  gene_modules,
  gene_adjacency,
  module_adjacency,
  hub_genes = NULL,
  closest_node_for_module = NULL,
  trajectory_object = NULL,
  percentage_non_hub_nodes = 1,
  percentage_edges = 0.2
)
```

## Arguments

- gene_modules:

  A list mapping modules to gene vectors.

- gene_adjacency:

  A gene adjacency matrix.

- module_adjacency:

  A module adjacency matrix, or NULL to infer one.

- hub_genes:

  Optional data frame containing hub genes.

- closest_node_for_module:

  Optional named vector mapping modules to trajectory nodes.

- trajectory_object:

  Optional trajectory object used when module_adjacency is NULL.

- percentage_non_hub_nodes:

  Fraction of non-hub nodes kept per module.

- percentage_edges:

  Fraction of strongest edges retained.

## Value

A list with node metadata and filtered edge table.
