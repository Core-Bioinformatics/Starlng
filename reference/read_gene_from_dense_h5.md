# Read expression matrix from a dense HDF5 file

This function read specific rows, associated with some genes, from a
dense matrix stored in a HDF5 file.

## Usage

``` r
read_gene_from_dense_h5(
  gene_names,
  matrix_h5_path,
  index_genes = NULL,
  check_intersect = TRUE,
  add_rownames = FALSE
)
```

## Arguments

- gene_names:

  List of names of genes to be read.

- matrix_h5_path:

  Path to the HDF5 file where the expression matrix is stored. The file
  must contain the following fields:

  - genes: the names of the genes in the matrix.

  - cells: the names of the cells in the matrix.

  - expression_matrix: the expression matrix

- index_genes:

  Named vector where the names are the genes and the values the index of
  the row associated with a gene. If NULL, the function will calculate
  automatically this index. Defaults to NULL.

- check_intersect:

  Logical indicating if the function should check if all the provided
  genes are present in the matrix. Defaults to TRUE.

- add_rownames:

  Logical indicating if the function should add the cell names as column
  names for the extracted matrix. Defaults to FALSE.

## Value

The expression matrix having the requested genes as rows.
