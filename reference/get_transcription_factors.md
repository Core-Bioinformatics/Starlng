# Get Enriched Transcription Factors

Performs transcription-factor enrichment against the TRANSFAC database
using `gprofiler2`, groups the transcription factors by their name and
identifies the associated genes.

## Usage

``` r
get_transcription_factors(
  gene_list,
  organism = "hsapiens",
  p_value_threshold = 0.05,
  correction_method = "fdr",
  bg_genes = NULL,
  ...
)
```

## Arguments

- gene_list:

  A character vector of genes, or a named list of gene vectors.

- organism:

  Organism code used by gprofiler2.

- p_value_threshold:

  Numeric p-value threshold used for filtering.

- correction_method:

  Multiple-testing correction method. Defaults to 'fdr'.

- bg_genes:

  Optional background gene universe. The domain scope of the analysis is
  automatically adjusted based on whether this is provided or not.

- ...:

  Additional arguments passed to gprofiler2::gost.

## Value

A list with enrichment table and TF-associated genes, or NULL if no
significant enrichment is found.
