# Server - Prepare Session

Creates the backend interface that initializes the variables and the
environment used inside the Starlng Shiny application.

## Usage

``` r
prepare_session(
  parent_session,
  reactive_dim,
  height_ratio = 0.7,
  enrichment_organism = "hsapiens"
)
```

## Arguments

- parent_session:

  An object linking to the Shiny session of the parent container.

- reactive_dim:

  A reactive variable that returns the dimensions of the browser window.
  This reactive variable is used to update the dimensions of the plots
  inside the Shiny app.

- height_ratio:

  Variable indicating the ratio between the height of the plots and the
  height of the browser window.

- enrichment_organism:

  The organism to be used for the gene enrichment analysis. Defaults to
  "hsapiens" (human). For the list of supported organisms, please
  consult the `gprofiler2` package documentation.

## Note

This function is a shiny module function and should be used in the
context of the app created using the `starlng_write_app` function.
