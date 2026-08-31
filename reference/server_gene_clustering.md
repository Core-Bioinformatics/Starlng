# Server - Gene clustering

Creates the backend interface for the Gene Clustering panel inside the
Starlng Shiny application.

## Usage

``` r
server_gene_clustering(id, filtered_genes)
```

## Arguments

- id:

  The id of the shiny module, used to access the UI elements.

- filtered_genes:

  A reactive expression that contains the filtered genes that will be
  used for the clustering.

## Note

This function is a shiny module function and should be used in the
context of the app created using the `starlng_write_app` function.
