# Gene Expression voting scheme

This function performs a voting scheme to identify cells that are
expressing one or multiple genes above a given threshold. The function
allows the user to introduce a factor of relaxation in the voting scheme
by setting the minimum number of genes that should be coexpressed in a
cell for it to be considered. For the obtained set of cells, the
function summarizes the expression values of the genes using a
user-defined function.

## Usage

``` r
voting_scheme(
  expression_matrix,
  genes = rownames(expression_matrix),
  thresh_percentile = 0.25,
  thresh_value = 0,
  n_coexpressed_thresh = length(genes),
  summary_function = NULL
)
```

## Arguments

- expression_matrix:

  The gene by cell expression matrix to be used. The matrix should have
  defined rownames and colnames.

- genes:

  A vector with the genes to be used for the voting scheme.

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

- summary_function:

  A function that takes as argument a numeric vector and summarises it
  into a single value. If set to NULL, the function will return a factor
  with the selected cells. Defaults to NULL.

## Value

A list with the information of the cells selected by the voting scheme.
If a summary function is provided, the list will contain the summary
values for each cell.
