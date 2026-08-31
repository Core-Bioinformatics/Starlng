# Get PCA reduction from an expression matrix

This function calculates the PCA embedding of an expression matrix.

## Usage

``` r
pca_reduction(expr_matrix, npcs = 30, approx = FALSE, ...)
```

## Arguments

- expr_matrix:

  The gene by cell expression matrix to be used. It is expected that the
  matrix is already normalized and scaled.

- npcs:

  The number of principal components to be used. Defaults to 30.

- approx:

  Logical indicating if deterministic PCA (with `prcomp`) or approximate
  PCA (with `irlba`) should be used. Defaults to FALSE.

- ...:

  Aditional parameters passed to the PCA function.

## Value

A gene by npcs matrix containing the PCA embedding.

## Note

The PCA is applied in this context to the genes, not to the cells.
