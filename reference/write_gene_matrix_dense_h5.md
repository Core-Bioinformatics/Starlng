# Write expression matrix to a dense HDF5 file

This function writes a expression matrix to a HDF5 file in a dense
format. The created file will contain three fields:

- genes: the names of the genes in the matrix.

- cells: the names of the cells in the matrix.

- expression_matrix: the expression matrix

## Usage

``` r
write_gene_matrix_dense_h5(
  expression_matrix,
  file_path,
  compression_level = 7,
  chunk_size = 100,
  all_at_once = TRUE,
  garbage_thresh = 2
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

- chunk_size:

  The size of the chunks to be used. Defaults to 100.

- all_at_once:

  Logical indicating if the writing process will be done by providing
  the entire matrix or by writing chunks of the matrix. It is
  recommended to write the entire matrix at once, as no noticeable
  speedup is observed in the other case. Defaults to TRUE.

- garbage_thresh:

  The threshold, measured GB, to trigger a garbage collection while
  writing the matrix by chunks.

## Value

The function does not return anything. It will create the HDF5 file with
the matrix stored in it.
