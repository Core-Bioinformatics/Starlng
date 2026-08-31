# Select cells by gene expression

This function filters from the dataset a group of cells which are
expressing a set of genes above a given threshold. For this, the voting
scheme will be applied.

## Usage

``` r
select_cells_by_gene_expr(
  expression_matrix,
  genes,
  thresh_percentile = 0.25,
  thresh_value = 0,
  n_coexpressed_thresh = length(genes)
)
```

## Arguments

- expression_matrix:

  The gene by cell expression matrix to be used.

- genes:

  A vector with the names of the genes to be used for the filtering.

- thresh_percentile:

  The percentile to be used as threshold for the expression values to
  select the cells associated with each gene. If set to 0, the threshold
  will be the value defined in `thresh_value`. Defaults to 0.25.

- thresh_value:

  The value to be used as threshold for the expression values to select
  the cells associated with each gene. Defaults to 0.

- n_coexpressed_thresh:

  The number of genes that should be coexpressed in a cell for it to be
  considered as selected. Defaults to the number of provided genes.

## Value

A list of cells that are selected using the voting scheme.
