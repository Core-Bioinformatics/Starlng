# Summarize Transcription Factor Statistics

Converts transcription-factor enrichment output to a tabular summary.
The summary includes, for each TF, the number of associated genes and
optionally the list of associated genes as a comma-separated string.

## Usage

``` r
get_tf_stats(
  transcription_factors,
  module_name = NULL,
  include_intersection_set = FALSE
)
```

## Arguments

- transcription_factors:

  Output from get_transcription_factors.

- module_name:

  Optional module name used when processing one module.

- include_intersection_set:

  Logical indicating whether to include associated genes as a
  comma-separated string.

## Value

A data frame with TF statistics, or NULL if no transcription factors are
available.
