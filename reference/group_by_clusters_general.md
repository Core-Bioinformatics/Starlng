# Group the partitions by the number of clusters

This function groups the partitions by the number of clusters. The
grouping can be done at each level of the nested list by keeping the
hierarchy of parameters.

## Usage

``` r
group_by_clusters_general(clustering_list, start_level = 1)
```

## Arguments

- clustering_list:

  A nested list of partitions. Each level is associated with values of a
  specific parameters. The last level contains the partitions.

- start_level:

  The level at which the grouping should start. The default is 1.

## Value

A nested list of partitions grouped by the number of clusters.
