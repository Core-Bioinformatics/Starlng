# Plot Top Enriched Terms

Builds a dot plot with the top enriched terms per module.

## Usage

``` r
plot_enrichment_top_terms(
  enrichment_result,
  top_n = 2,
  colour_column = NULL,
  point_size_range = c(3, 10),
  font_size = 10,
  words_per_line = 5
)
```

## Arguments

- enrichment_result:

  A data frame produced by enrichment analysis.

- top_n:

  Number of top terms to keep per module. The ranking is determined by
  p-value. Defaults to 2.

- colour_column:

  Optional column name used to map point size. We suggest using
  "intersection_size".

- point_size_range:

  Numeric vector of length 2 defining point size range.

- font_size:

  Base font size used in the plot theme.

- words_per_line:

  Integer indicating how many words should be displayed per line for
  every term name. Defaults to 5.

## Value

A ggplot object, or NULL if the input is empty.
