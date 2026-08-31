# Write expression matrix to a sparse HDF5 file

This function writes a expression matrix to a HDF5 file in a sparse CSC
format. The created file will contain three fields:

- genes: the names of the genes in the matrix.

- cells: the names of the cells in the matrix.

- expression_matrix: the expression matrix stored in CSC format.

## Usage

``` r
write_gene_matrix_sparse_h5(
  expression_matrix,
  file_path,
  compression_level = 7
)
```

## Arguments

- expression_matrix:

  The expression matrix to be written. The matrix should have defined
  rownames and colnames. If they are not defined, they will be generated
  as `gene_{index}` and `cell_{index}` respectively.

- file_path:

  The path to the HDF5 file where the matrix will be stored.

- compression_level:

  The compression level to be used. Defaults to 7.

## Value

The function does not return anything. It will create the HDF5 file with
the matrix stored in it.
