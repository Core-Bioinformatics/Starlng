# Get the by-cluster consistency of configurations

This function calculates the Element-Centric Consistency (ECC) of the
partitions that are grouped by the number of clusters.

## Usage

``` r
get_clusters_consistency(by_cluster_list, order_logic = "avg_agreement")
```

## Arguments

- by_cluster_list:

  A nested list of partitions that are grouped by the number of
  clusters.

- order_logic:

  Determines how to order the partitions associated with the same number
  of clusters. Can be either by frequency (most frequent first) or by
  the agreement (the partition which is the most similar with the
  others, first). Defaults to "agreement".

## Value

A nested list with the partitions ordered by agreement / frequency and
the overall ECC value.
