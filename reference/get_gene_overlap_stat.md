# Compute Gene Hub Statistics

Computes precision, recall and F1 overlap statistics between the
population induced by each gene and the overall population described by
the module it belongs to. If the gene weight is provided, the function
also calculates the combined score (defined as normalised weight \* F1
score) to determine hub genes.

## Usage

``` r
get_gene_overlap_stat(
  expr_matrix,
  gene_modules,
  module_expr,
  gene_expression_threshold = 0,
  gene_expression_percentile = 0,
  module_expression_threshold = 0,
  module_expression_percentile = 0,
  scale = TRUE,
  total_weight = NULL
)
```

## Arguments

- expr_matrix:

  A gene-by-cell expression matrix.

- gene_modules:

  A character vector of genes, or module-specific gene list.

- module_expr:

  A numeric/logical vector, or a list of such vectors.

- gene_expression_threshold:

  Threshold for defining expressed genes.

- gene_expression_percentile:

  Optional percentile threshold for gene expression.

- module_expression_threshold:

  Threshold for module activity.

- module_expression_percentile:

  Optional percentile threshold for module activity.

- scale:

  Logical indicating whether values should be min-max scaled.

- total_weight:

  Optional vector of additional gene weights.

## Value

A data frame with overlap metrics per gene.
