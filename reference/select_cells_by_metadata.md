# Select cells by metadata

This function filters from the dataset a group of cells which is defined
as a combination of filters determined by multiple metadata.

## Usage

``` r
select_cells_by_metadata(metadata_df, metadata_combinations)
```

## Arguments

- metadata_df:

  A dataframe with the metadata information. The rows of the dataframe
  should be the cells and the columns the metadata information.

- metadata_combinations:

  A list with the metadata columns and the unique groups of the metadata
  that should be used to filter the cells. The names of the list should
  be in the columns of the `metadata_df`.

## Value

The list of cells that are defined by the intersection of provided
groups.
