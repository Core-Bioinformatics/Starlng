# Parallel NN indexing

This function provides a parallel approach of finding the nearest
neighbours of each points on a given embedding. The parallelisation is
done by splitting the data into chunks and running the
[`RANN::nn2`](https://jefferislab.github.io/RANN/reference/nn2.html)
function on each chunk in parallel.

## Usage

``` r
parallel_nn2_idx(embedding, k)
```

## Arguments

- embedding:

  The embedding matrix where each row represents a point and each column
  a dimension.

- k:

  The number of nearest neighbours to be found.

## Value

A matrix where each row contains the indices of the k-nearest neighbours
of the corresponding point.
