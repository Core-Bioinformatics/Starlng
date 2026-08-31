# Generate Cell Heatmap

Generates a heatmap of gene-family expression across cells, optionally
ordered by pseudotime and annotated by metadata.

## Usage

``` r
generate_cell_heatmap(
  expression_matrix,
  gene_family_list,
  metadata_df,
  metadata_name = NULL,
  legend_name = "",
  apply_scale = TRUE,
  k_smooth = 30,
  cap = 20,
  axis_text_size = 10,
  discrete_colour_list = NULL,
  continuous_colors = NULL,
  show_pseudotime_legend = FALSE,
  show_metadata_legend = TRUE,
  show_modules_legend = FALSE,
  show_expression_legend = TRUE
)
```

## Arguments

- expression_matrix:

  A gene-by-cell expression matrix.

- gene_family_list:

  A named list of genes grouped by gene family.

- metadata_df:

  A data frame with cell metadata; rownames must match cell names in
  `expression_matrix`.

- metadata_name:

  Optional metadata column used for the top annotation.

- legend_name:

  Title used for the heatmap legend.

- apply_scale:

  Logical indicating whether to z-score each gene across cells before
  plotting.

- k_smooth:

  Integer kernel size used to smooth expression along cells.

- cap:

  Numeric cap applied to scaled or raw expression values.

- axis_text_size:

  Font size used for row labels and annotations.

- discrete_colour_list:

  Optional named list of discrete colour vectors.

- continuous_colors:

  Optional vector of colors used for the heatmap body.

- show_pseudotime_legend:

  Logical indicating whether to show the pseudotime legend. Defaults to
  FALSE.

- show_metadata_legend:

  Logical indicating whether to show the metadata legend. Defaults to
  TRUE.

- show_modules_legend:

  Logical indicating whether to show the gene family legend. Defaults to
  FALSE.

- show_expression_legend:

  Logical indicating whether to show the expression legend. Defaults to
  TRUE.

## Value

A `ComplexHeatmap` heatmap object.
