# Recommend a pseudotime ordering

This function provides a recommendation of a pseudotime ordering based
on the metadata available in the monocle object. The recommendation is
done by selecting an end node from the diameter of the trajectory graph.
The function identifies the metadata group whose centre is the closest
to the starting point.

## Usage

``` r
get_pseudotime_recommendation(monocle_object)
```

## Arguments

- monocle_object:

  A Monocle3 object.

## Value

A list that contains the recommended metadata column, the subgroup and
the pseudotime values.
