# Filter Enrichment Results

Filters enrichment results by p-value and optional intersection-size
bounds. It assumes that the results are from a run of
[`gprofiler2::gost`](https://rdrr.io/pkg/gprofiler2/man/gost.html).

## Usage

``` r
filter_enrichment_results(
  enrichment_result,
  p_value_threshold = 0.05,
  intersection_size_upper = NULL,
  intersection_size_lower = NULL
)
```

## Arguments

- enrichment_result:

  A data frame or list containing enrichment results. If a list is
  provided, each element should correspond to a module. The elements
  will be automatically combined into a one data frame.

- p_value_threshold:

  Numeric threshold for adjusted p-values. Defaults to 0.05.

- intersection_size_upper:

  Optional upper bound for intersection size.

- intersection_size_lower:

  Optional lower bound for intersection size.

## Value

A filtered enrichment data frame, or NULL if no rows pass filters.
