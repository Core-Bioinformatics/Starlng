#' @importFrom rlang .data
#' @importFrom dplyr %>% 
#' 
gene_expr_palette <- c(
    "grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
    "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"
)

plot_umap_discrete <- function(umap_df,

                               cell_sort_order = NULL,
                               colors_list = NULL,
                               show_labels = TRUE) {

    unique_vals <- unique(umap_df$cell_info)
    n_unique <- length(unique_vals)

    if (!is.null(cell_sort_order) && length(cell_sort_order) != n_unique) {
        cell_sort_order <- NULL
    }

    if (!is.null(cell_sort_order) && !all(cell_sort_order %in% unique_vals)) {
        cell_sort_order <- NULL
    }

    if (!is.null(cell_sort_order)) {
        umap_df$cell_info <- factor(umap_df$cell_info, levels = cell_sort_order)
        umap_df <- umap_df %>% dplyr::arrange(.data$cell_info)
    }

    gplot_obj <- ggplot2::ggplot(umap_df, ggplot2::aes(
        x = .data$UMAP1, y = .data$UMAP2, colour = .data$cell_info)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(colour = "Cluster")
    
    if (is.null(colors_list)) {
        return(gplot_obj)
    }

    used_colors <- colors_list[[as.character(n_unique)]]

    return(gplot_obj + ggplot2::scale_color_manual(values = used_colors))
}

plot_umap_continuous <- function(umap_df,

                                 cell_sort_order = c("lowest", "highest", "default"),
                                 color_scheme = NULL) {

    cell_sort_order <- cell_sort_order[1]
    if (cell_sort_order == "lowest") {
        umap_df <- umap_df %>% dplyr::arrange(dplyr::desc(.data$cell_info))
    } else if (cell_sort_order == "highest") {
        umap_df <- umap_df %>% dplyr::arrange(.data$cell_info)
    }
    gplot_obj <- ggplot2::ggplot(umap_df, ggplot2::aes(
        x = .data$UMAP1, y = .data$UMAP2, colour = .data$cell_info)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(colour = "Cluster")
    
    if (is.null(color_scheme)) {
        return(gplot_obj)
    }

    return(gplot_obj + ggplot2::scale_color_gradientn(colours = color_scheme))
}

plot_umap <- function(umap_embedding,
                      cell_info,
                      cell_sort_order = c("lowest", "highest", "default"),
                      discrete_colors = NULL,
                      continuous_colors = NULL) {
    df <- data.frame(umap_embedding)
    colnames(df) <- c("UMAP1", "UMAP2")
    df$cell_info <- cell_info
    
      if (inherits(df$cell_info, c("factor", "character"))) {
        return(plot_umap_discrete(
            umap_df = df,
            cell_sort_order = cell_sort_order,
            colors_list = discrete_colors)
        )
    }

    return(plot_umap_continuous(
        umap_df = df,
        cell_sort_order = cell_sort_order,
        color_scheme = continuous_colors
    ))
}
