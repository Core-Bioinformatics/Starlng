# Plot Module Trends over Pseudotime

Fits smoothed module-expression trajectories over pseudotime and plots
module trends. If the number of points is large (\> 7000 cells), the
function uses a GAM-based smoothing approach for efficiency, otherwise
it uses LOESS smoothing.

## Usage

``` r
plot_module_trends_over_pseudotime(
  expression_list,
  pseudotime,
  module_colours = NULL,
  linewidth = 1.5,
  axis_text_size = 10,
  legend_text_size = 10,
  label_size = 10,
  show_labels = TRUE,
  smooth_threshold_max = 0.2
)
```

## Arguments

- expression_list:

  A list, matrix, or data frame with module expression per cell.

- pseudotime:

  Numeric pseudotime vector.

- module_colours:

  Optional named vector of module colours.

- linewidth:

  Line width for smoothed trajectories.

- axis_text_size:

  Text size for axis labels and ticks.

- legend_text_size:

  Text size for legend labels and title.

- label_size:

  Text size for module labels.

- show_labels:

  Logical indicating whether module labels are shown. The labels will be
  located at the maximum point of the smoothed trajectory.

- smooth_threshold_max:

  Threshold used to discard modules where the maximum value of the
  smoothed trajectory is below this value.

## Value

A ggplot object.
