# Server - Tabs update

Creates the backend interface that dynamically shows the gene module
tabs inside the Starlng Shiny application after running the gene
clustering.

## Usage

``` r
update_tabs(session)
```

## Arguments

- session:

  The shiny session object used to update the UI elements.

## Note

This function is a shiny module function and should be used in the
context of the app created using the `starlng_write_app` function.
