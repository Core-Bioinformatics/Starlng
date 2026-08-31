# Create Starlng Shiny app from a normalized gene expression matrix

Runs the entire Starlng pipeline on a normalized gene expression matrix
and generates a folder with all necessary files to run a Shiny app.

## Usage

``` r
starlng_write_app_default(
  folder_path,
  expression_matrix,
  metadata_df = NULL,
  pca_embedding = NULL,
  umap_embedding = NULL,
  app_title_name = "",
  learn_graph_parameters = list(nodes_per_log10_cells = 30, learn_graph_controls =
    list(eps = 1e-05, maxiter = 100)),
  gene_filtering_function = function(info_gene_df) {
     rownames(info_gene_df %>%
    dplyr::filter(.data$morans_I > 0.1, .data$q_value < 0.05))
 },
  clustering_parameters = list(n_neighbours = seq(from = 5, to = 50, by = 5), graph_type
    = "snn", prune_value = -1, resolutions = list(RBConfigurationVertexPartition =
    seq(from = 0.1, to = 2, by = 0.1), RBERVertexPartition = NULL,
    ModularityVertexPartition = NULL), number_iterations = 5, number_repetitions = 100),
  ecc_threshold = 0.9,
  freq_threshold = 30,
  enrichment_organism = "hsapiens",
  skip_tf_identification = FALSE,
  save_entire_monocle = TRUE,
  discrete_colours = list(),
  continuous_colours = list(),
  max_n_colors = 40,
  verbose = FALSE,
  compression_level = 7,
  chunk_size = 100,
  nthreads = 1
)
```

## Arguments

- folder_path:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- expression_matrix:

  A normalized gene by cell expression matrix. The matrix should have
  the rownames and the colnames defined.

- metadata_df:

  A data frame with cell metadata specific to the experiment. If NULL,
  the function will create a one-group metadata with the name
  "one_level".

- pca_embedding:

  A matrix with the PCA embedding of the cells. If NULL, the function
  will calculate the PCA embedding while creating the monocle3 object.

- umap_embedding:

  A matrix with the UMAP embedding of the cells. If NULL, the function
  will calculate the UMAP embedding while creating the monocle3.

- app_title_name:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- learn_graph_parameters:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- gene_filtering_function:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- clustering_parameters:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- ecc_threshold:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- freq_threshold:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- enrichment_organism:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- skip_tf_identification:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- save_entire_monocle:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- discrete_colours:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- continuous_colours:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- max_n_colors:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- verbose:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- compression_level:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- chunk_size:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

- nthreads:

  Check
  [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  documentation.

## Note

For more details about the other parameters and the content of the shiny
folder, consult the documentation of the
[`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
function.
