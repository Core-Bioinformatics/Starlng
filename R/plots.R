#' @importFrom dplyr %>% .data
#' 
gene_expr_palette <- c(
    "grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
    "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"
)

plot_umap_discrete <- function(umap_df,
                               cell_sort_order = NULL,
                               cell_size = 0.3,
                               cell_alpha = 0.8,
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
        ggplot2::geom_point(size = cell_size, alpha = cell_alpha) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "bottom",
            aspect.ratio = 1,
            legend.text = ggplot2::element_text(size = 10),
        ) +
        ggplot2::labs(colour = "Cluster")
    
    if (is.null(colors_list)) {
        return(gplot_obj)
    }

    used_colors <- colors_list[[as.character(n_unique)]]

    return(gplot_obj + ggplot2::scale_color_manual(values = used_colors))
}

sort_umap_discrete <- function(umap_df,
                               cell_sort_order = NULL) {

    unique_vals <- as.character(unique(umap_df$cell_info))
    n_unique <- length(unique_vals)

    if (!is.null(cell_sort_order) && length(cell_sort_order) != n_unique) {
        cell_sort_order <- NULL
    }

    if (!is.null(cell_sort_order) && !all(cell_sort_order %in% unique_vals)) {
        cell_sort_order <- NULL
    }

    if (is.null(cell_sort_order)) {
        return(umap_df)
    }

    umap_df$cell_info <- factor(umap_df$cell_info, levels = cell_sort_order)
    umap_df <- umap_df %>% dplyr::arrange(.data$cell_info)

    return(umap_df)
}

sort_umap_continuous <- function(umap_df,
                                 cell_sort_order = NULL) {
    cell_sort_order <- cell_sort_order[1]
    if (cell_sort_order == "lowest") {
        return(umap_df %>% dplyr::arrange(dplyr::desc(.data$cell_info)))
    } 
    
    if (cell_sort_order == "highest") {
        return(umap_df %>% dplyr::arrange(.data$cell_info))
    }

    return(umap_df)
}

plot_umap_continuous <- function(umap_df,
                                 cell_sort_order = c("lowest", "highest", "default"),
                                 cell_size = 0.3,
                                 cell_alpha = 0.8,
                                 color_scheme = NULL) {

    cell_sort_order <- cell_sort_order[1]
    if (cell_sort_order == "lowest") {
        umap_df <- umap_df %>% dplyr::arrange(dplyr::desc(.data$cell_info))
    } else if (cell_sort_order == "highest") {
        umap_df <- umap_df %>% dplyr::arrange(.data$cell_info)
    }
    gplot_obj <- ggplot2::ggplot(umap_df, ggplot2::aes(
        x = .data$UMAP1, y = .data$UMAP2, colour = .data$cell_info)) +
        ggplot2::geom_point(size = cell_size, alpha = cell_alpha) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "bottom",
            aspect.ratio = 1
        ) +
        ggplot2::labs(colour = "Cluster")
    
    if (is.null(color_scheme)) {
        return(gplot_obj)
    }

    return(gplot_obj + ggplot2::scale_color_gradientn(colours = color_scheme))
}

# TODO add the option to plot text above the discrete groups
#' Plot UMAP with Discrete or Continuous Coloring
#' 
#' @description This function plots a UMAP embedding with points colored by either
#' discrete or continuous cell information. The function provides options for
#' sorting the cells based on the coloring variable, adjusting point size and
#' alpha, and customizing the color scheme for both discrete and continuous
#' coloring.
#' 
#' @param umap_embedding A matrix or data frame containing the UMAP coordinates
#' of the cells. It should have two columns corresponding to UMAP1 and UMAP2.
#' @param cell_info A vector containing the information to color the cells by. It
#' can be either a factor/character vector for discrete coloring or a numeric vector
#' for continuous coloring. The length of this vector should match the number of rows in
#' `umap_embedding`.
#' @param mtd_name A string specifying the name of the metadata variable to be
#' used in the legend.
#' @param cell_sort_order A character vector specifying the order of cells based on
#' the `cell_info` variable. For discrete coloring, it should contain the unique values
#' of `cell_info` in the desired order. For continuous coloring, it can be either
#' "lowest", "highest" or "default" to sort cells by the `cell_info` values in
#' ascending, descending or no particular order, respectively.
#' @param scale_values A logical indicating whether to scale the `cell_info`
#' values between 0 and 1 for continuous coloring. Defaults to FALSE.
#' @param cell_size A numeric value specifying the size of the points in the plot.
#' Defaults to 0.3.
#' @param cell_alpha A numeric value between 0 and 1 specifying the transparency of
#' the points in the plot. Defaults to 0.8.
#' @param legend_text_size A numeric value specifying the size of the text in
#' the legend. Defaults to 10.
#' @param legend_key_size A numeric value specifying the size of the legend keys
#' for discrete coloring. Defaults to 0.5.
#' @param axis_text_size A numeric value specifying the size of the axis text in
#' the plot. Defaults to 10.
#' @param show_labels A logical indicating whether to show labels for discrete
#' groups. Defaults to TRUE.
#' @param label_size A numeric value specifying the size of the labels for
#' discrete groups. Defaults to 10.
#' @param colourbar_width A numeric value specifying the width of the colorbar
#' for continuous coloring. Defaults to 50 points. If set to 0, the colorbar
#' will be hidden.
#' @param discrete_colors A list of color vectors for discrete coloring, where the
#' names of the list correspond to the number of unique values in `cell_info`.
#' If NULL, the default ggplot2 colors will be used.
#' @param continuous_colors A vector of colors for continuous coloring. If NULL, the
#' default ggplot2 colors will be used.
#' 
#' @return A ggplot object representing the UMAP plot with the specified
#' coloring and customization options.
#' @export
plot_umap <- function(umap_embedding,
                      cell_info,
                      mtd_name,
                      cell_sort_order = c("lowest", "highest", "default"),
                      scale_values = FALSE,
                      cell_size = 0.3,
                      cell_alpha = 0.8,
                      legend_text_size = 10,
                      legend_key_size = 0.5,
                      axis_text_size = 10,
                      show_labels = TRUE,
                      label_size = 10,
                      colourbar_width = 50,
                      discrete_colors = NULL,
                      continuous_colors = NULL) {
    df <- data.frame(umap_embedding)
    colnames(df) <- c("UMAP1", "UMAP2")
    df$cell_info <- cell_info
    is_discrete <- inherits(df$cell_info, c("factor", "character"))

    if (is_discrete) {
        n_unique <- length(unique(df$cell_info))
        if (n_unique > 200) {
            df$cell_info <- as.numeric(as.factor(df$cell_info))
            is_discrete <- FALSE
        } else if (isFALSE(as.character(n_unique) %in% names(discrete_colors))) {
            discrete_colors[[as.character(n_unique)]] <- qualpalr::qualpal(n = n_unique)$hex
        }
    }

    if (is_discrete) {
        df <- sort_umap_discrete(df, cell_sort_order)
    } else {
        df <- sort_umap_continuous(df, cell_sort_order)
    }

    if (!is_discrete && scale_values) {
        df$cell_info <- scale_min_max(df$cell_info)
        mtd_name <- paste(mtd_name, "(scaled)")
    }

    # TODO add fast_rendering option with scattermost and take care of faceting

    gplot_obj <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$UMAP1, y = .data$UMAP2, colour = .data$cell_info)) +
        ggplot2::geom_point(size = cell_size, alpha = cell_alpha) +
        # scattermore::geom_scattermore(pointsize = cell_size * 2, alpha = cell_alpha) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "bottom",
            # legend text above the legend symbols

            aspect.ratio = 1,
            legend.text = ggplot2::element_text(size = legend_text_size),
            legend.title = ggplot2::element_text(size = legend_text_size),
            axis.title = ggplot2::element_text(size = axis_text_size),
            axis.text = ggplot2::element_text(size = axis_text_size),
            plot.title = ggplot2::element_text(size = axis_text_size * 1.5)
        ) +
        ggplot2::labs(colour = mtd_name)

    if (is_discrete) {
        gplot_obj <- gplot_obj +
            ggplot2::guides(
                colour = ggplot2::guide_legend(
                    title.position = "top",
                    override.aes = list(size = legend_key_size),
                    
                    # FIXME there is a problem with the unique values displayed on legend
                    # if legend text size is too big, the text won't fit anymore
                    # it should create more rows
                    # We should automatically adjust the number of rows based on the width and text size
                    ncol = 6
                )
            )
        if (is.null(discrete_colors)) {
            return(gplot_obj)
        }

        unique_vals <- unique(df$cell_info)

        n_unique <- length(unique_vals)
        used_colors <- discrete_colors[[as.character(n_unique)]]

        if (show_labels) {
            summary_median <- df %>% 
                dplyr::group_by(.data$cell_info) %>%
                dplyr::summarise(gmed = list(Gmedian::Gmedian(cbind(.data$UMAP1, .data$UMAP2))[1, ]), .groups = "drop") %>%
                tidyr::unnest_wider(.data$gmed, names_sep = " ") 
            colnames(summary_median)[2:3] <- c("UMAP1", "UMAP2")

            gplot_obj <- gplot_obj +
                ggplot2::geom_text(
                    data = summary_median,
                    ggplot2::aes(label = .data$cell_info, x = .data$UMAP1, y = .data$UMAP2),
                    size = label_size,
                    colour = "black",
                    vjust = -0.5
                )
        }

        return(gplot_obj + ggplot2::scale_color_manual(values = used_colors))
    }

    if (is.null(continuous_colors) || colourbar_width == 0) {
        return(gplot_obj)
    }

    return(gplot_obj +
        ggplot2::scale_color_gradientn(colours = continuous_colors) +
        ggplot2::guides(
            colour = ggplot2::guide_colourbar(
                title.position = "top",
                direction = "horizontal",
                barwidth = grid::unit(colourbar_width, "points")
            )
        )
    )
}

plot_umap_gene <- function(umap_embedding,
                           gene_matrix,
                           gene_name,
                           thresh_percentile,
                           thresh_value,
                           top_perc_value,
                           n_coexpressed_threshold,
                           summarise_expr,
                           scale_values = FALSE,
                           trajectory_object = NULL,
                           trajectory_width = 0.75,
                           cell_sort_order = c("lowest", "highest", "default"),
                           cell_size = 0.3,
                           cell_alpha = 0.8,
                           legend_text_size = 10,
                           axis_text_size = 10,
                           colourbar_width = 50,
                           continuous_colors = NULL) {
    summary_function <- switch(
        summarise_expr,
        "binary" = NULL,
        "average expressed" = function(x) { mean(x[x > 0]) },
        "average" = mean,
        "#genes" = function(x) { sum(x > 0) },
        NULL
    )

    if (summarise_expr != "binary") {
        legend_name <- "Gene expression"
        cell_ordering <- cell_sort_order
    } else {
        legend_name <- "Cells above\nthreshold"
        cell_ordering <- c("not selected", "selected")
    }

    gene_info <- voting_scheme(
        expression_matrix = gene_matrix,
        genes = gene_name,
        thresh_percentile = thresh_percentile / 100,
        thresh_value = thresh_value,
        n_coexpressed_thresh = n_coexpressed_threshold,
        summary_function = summary_function
    )

    if (summarise_expr != "binary") {
        qval <- stats::quantile(gene_info[gene_info > 0], 1 - top_perc_value)
        gene_info[gene_info < qval] <- 0
    }

    return(
        plot_umap(
            umap_embedding = umap_embedding,
            cell_info = gene_info,
            mtd_name = legend_name,
            cell_sort_order = cell_ordering,
            scale_values = scale_values,
            cell_size = cell_size,
            cell_alpha = cell_alpha,
            legend_text_size = legend_text_size,
            axis_text_size = axis_text_size,
            discrete_colors = list(
                "1" = "red",
                "2" = c("gray", "red")
            ),
            colourbar_width = colourbar_width,
            continuous_colors = continuous_colors
        ) +
        plot_trajectory_graph(
            trajectory_object = trajectory_object,
            edge_size = trajectory_width,
            plot_nodes = 0
        )$layers
    )
}

# assumes the shiny context with the environment env
plot_umap_gene_shiny <- function(shiny_env,
                                 gene_name,
                                 thresh_percentile,
                                 thresh_value,
                                 top_perc_value,
                                 n_coexpressed_threshold,
                                 summarise_expr,
                                 scale_values = FALSE,
                                 trajectory_width = 0.75,
                                 cell_sort_order = c("lowest", "highest", "default"),
                                 cell_size = 0.3,
                                 cell_alpha = 0.8,
                                 legend_text_size = 10,
                                 axis_text_size = 10,
                                 colourbar_width = 50,
                                 continuous_colors = NULL) {
    gene_matrix <- read_gene_from_dense_h5(
        gene_names = gene_name,
        matrix_h5_path = file.path("objects", "expression.h5"),
        index_genes = shiny_env$genes[gene_name],
        check_intersect = FALSE
    )
    # gene_matrix <- read_gene_from_sparse_h5(
    #     gene_names = gene_name,
    #     index_genes = shiny_env$genes[gene_name],
    #     sparse_mat = shiny_env$expr_mat
    # )

    plot_umap_gene(
        umap_embedding = shiny_env$umap_df,
        gene_matrix = gene_matrix,
        gene_name = gene_name,
        thresh_percentile = thresh_percentile,
        thresh_value = thresh_value,
        top_perc_value = top_perc_value,
        scale_values = scale_values,
        n_coexpressed_threshold = n_coexpressed_threshold,
        summarise_expr = summarise_expr,
        trajectory_object = shiny_env$trajectory_object,
        trajectory_width = trajectory_width,
        cell_sort_order = cell_sort_order,
        cell_size = cell_size,
        cell_alpha = cell_alpha,
        legend_text_size = legend_text_size,
        axis_text_size = axis_text_size,
        colourbar_width = colourbar_width,
        continuous_colors = continuous_colors
    )
}

plot_umap_gene_modules_shiny <- function(shiny_env,
                                         gene_modules,
                                         summarise_expr,
                                         scale_values = TRUE,
                                         n_columns = 3,
                                         trajectory_width = 0.75,
                                         cell_sort_order = c("lowest", "highest", "default"),
                                         cell_size = 0.3,
                                         cell_alpha = 0.8,
                                         legend_text_size = 10,
                                         axis_text_size = 10,
                                         colourbar_height = 50,
                                         continuous_colors = NULL) {
    is_discrete_module <- summarise_expr == "binary"
    if (inherits(colourbar_height, "unit")) {
        bar_height <- colourbar_height
    } else {
        bar_height <- grid::unit(as.numeric(colourbar_height), "points")
    }

    plot_list <- lapply(names(gene_modules), function(gene_module) {
        plot_umap_gene_shiny(
            shiny_env = shiny_env,
            gene_name = gene_modules[[gene_module]],
            thresh_percentile = 0,
            thresh_value = 0,
            n_coexpressed_threshold = 0,
            scale_values = scale_values,
            summarise_expr = summarise_expr,
            trajectory_width = trajectory_width,
            cell_sort_order = cell_sort_order,
            cell_size = cell_size,
            cell_alpha = cell_alpha,
            legend_text_size = legend_text_size,
            axis_text_size = axis_text_size,
            colourbar_width = 0,
            continuous_colors = NULL
        ) + ggplot2::ggtitle(paste("module", gene_module))
    })

    names(plot_list) <- names(gene_modules)

    max_value <- max(sapply(plot_list, function(plot) {
        max(plot$data$cell_info)
    }))

    min_value <- min(sapply(plot_list, function(plot) {
        min(plot$data$cell_info)
    }))

    legend_plot <- plot_list[[1]] +
        ggplot2::theme(legend.position = "right")
    if (is_discrete_module) {
        legend_plot <- legend_plot +
            ggplot2::guides(
                colour = ggplot2::guide_legend(
                    title.position = "top"
                )
            )
    } else {
        legend_plot <- legend_plot +
            ggplot2::guides(
                colour = ggplot2::guide_colourbar(
                    title.position = "top",
                    barheight = bar_height
                )
            ) +
            ggplot2::scale_colour_gradientn(
                name = paste0("Gene\n", summarise_expr, "\n", ifelse(scale_values, "(scaled)", "")),
                colours = continuous_colors,
                limits = c(min_value, max_value)
            )
    }

    legend_gtable <- ggplot2::ggplotGrob(legend_plot)
    legend_idx <- grep("^guide-box", vapply(legend_gtable$grobs, function(g) g$name, character(1)))

    for (i in seq_along(plot_list)) {
        plot_list[[i]] <- plot_list[[i]] + ggplot2::theme(legend.position = "none")
    }

    patchwork_grid <- patchwork::wrap_plots(plotlist = plot_list, ncol = n_columns)

    if (!is_discrete_module) {
        patchwork_grid <- patchwork_grid &
            ggplot2::scale_colour_gradientn(
                name = paste0("Gene\n", summarise_expr, "\n", ifelse(scale_values, "(scaled)", "")),
                colours = continuous_colors,
                limits = c(min_value, max_value)
            )
    }

    if (length(legend_idx) == 0) {
        return(patchwork_grid)
    }

    (
        patchwork_grid |
            patchwork::wrap_elements(full = legend_gtable$grobs[[legend_idx[1]]])
    ) + patchwork::plot_layout(widths = c(1, 0.14))
}

#' Plot UMAP for Multiple Gene Modules
#' 
#' @description This function plots UMAP embeddings for multiple gene modules,
#' with options for discrete or continuous coloring based on the summarization
#' of gene expression. The function allows for customization of the plot
#' appearance, including point size, alpha, legend text size, axis text size,
#' and color schemes. It also provides the option to include trajectory graphs
#' if a trajectory object is provided.
#' 
#' @param module_summaries A data frame or a list of data frames containing the
#' summarized gene expression information for each module. If a data frame is
#' provided, it should have columns corresponding to each module. If a list is
#' provided, each element should be a data frame for a specific module.
#' @param umap_embedding A matrix or data frame containing the UMAP coordinates
#' of the cells. It should have two columns corresponding to UMAP1 and UMAP2.
#' @param trajectory_object An optional trajectory object that can be used to
#' overlay trajectory graphs on the UMAP plots. The trajectory object should be
#' compatible with the `plot_trajectory_graph` function.
#' @param legend_detail A string to be appended to the legend title for continuous
#' coloring, providing additional context about the summarization method used.
#' @param combine_legend A logical indicating whether to combine the legends of
#' all module plots into a single legend. Defaults to TRUE.
#' @param scale_values A logical indicating whether to scale the summarized gene
#' expression values between 0 and 1 for continuous coloring. Defaults to TRUE.
#' @param n_columns An integer specifying the number of columns to use when 
#' arranging the module plots. Defaults to 3.
#' @param trajectory_width A numeric value specifying the width of the
#' trajectory edges when overlaying the trajectory graph on the UMAP plots.
#' Defaults to 0.75.
#' @param cell_sort_order A character vector specifying the order of cells based on
#' the summarized gene expression values. It can be either "lowest", "highest" or
#' "default" to sort cells by the summarized values in ascending, descending or
#' no particular order, respectively.
#' @param cell_size A numeric value specifying the size of the points in the
#' UMAP plots. Defaults to 0.3.
#' @param cell_alpha A numeric value specifying the transparency of the points in the
#' UMAP plots. Defaults to 0.8.
#' @param legend_text_size A numeric value specifying the size of the text in the
#' legends. Defaults to 10.
#' @param axis_text_size A numeric value specifying the size of the axis text in the
#' UMAP plots. Defaults to 10.
#' @param legend_key_size A numeric value specifying the size of the legend keys
#' for discrete coloring. Defaults to 1.
#' @param colourbar_height A numeric value specifying the height of the colorbar
#' for continuous coloring. Defaults to 50 points. If set to 0, the colorbar will be hidden.
#' @param continuous_colors A vector of colors for continuous coloring. If NULL, the default ggplot2 colors will be used.
#' @param show_labels A logical indicating whether to show labels for discrete groups. Defaults to FALSE.
#' 
#' @return A patchwork object containing the UMAP plots for each gene module,
#' arranged according to the specified number of columns, with a combined legend
#' if `combine_legend` is TRUE.
#' @export
plot_umap_gene_modules <- function(
    module_summaries,
    umap_embedding,
    trajectory_object = NULL,
    legend_detail = "",
    combine_legend = TRUE,
    scale_values = TRUE,
    n_columns = 3,
    trajectory_width = 0.75,
    cell_sort_order = c("lowest", "highest", "default"),
    cell_size = 0.3,
    cell_alpha = 0.8,
    legend_text_size = 10,
    axis_text_size = 10,
    legend_key_size = 1,
    colourbar_height = 50,
    continuous_colors = NULL,
    show_labels = FALSE
) {
    is_list <- inherits(module_summaries, "list")
    if (!is_list) {
        is_discrete_module <- module_summaries[1, 1] %in% c("not selected", "selected")
    } else {
        is_discrete_module <- module_summaries[[1]][1] %in% c("not selected", "selected")
    }
    if (inherits(colourbar_height, "unit")) {
        bar_height <- colourbar_height
    } else {
        bar_height <- grid::unit(as.numeric(colourbar_height), "points")
    }

    if (is_discrete_module) {
        legend_name <- "Cells above\nthreshold"
        cell_ordering <- c("not selected", "selected")
    } else {
        legend_name <- "Gene expression"
        cell_ordering <- cell_sort_order
    }

    if (is_list) {
        n_columns <- min(n_columns, length(module_summaries))
        module_names <- names(module_summaries)
    } else {
        n_columns <- min(n_columns, ncol(module_summaries))
        module_names <- colnames(module_summaries)
    }

    plot_list <- lapply(module_names, function(gene_module) {
        gc()
        if (is_list) {
            cell_info <- module_summaries[[gene_module]]
        } else {
            cell_info <- module_summaries[, gene_module]
        }
        gplot_obj <- plot_umap(
            umap_embedding = umap_embedding,
            cell_info = cell_info,
            mtd_name = legend_name,
            cell_sort_order = cell_ordering,
            scale_values = scale_values,
            cell_size = cell_size,
            cell_alpha = cell_alpha,
            legend_text_size = legend_text_size,
            axis_text_size = axis_text_size,
            discrete_colors = list(
                "1" = "red",
                "2" = c("gray", "red")
            ),
            colourbar_width = 0,
            continuous_colors = NULL,
            show_labels = show_labels,
            legend_key_size = legend_key_size
        )

        if (!is.null(trajectory_object)) {
            gplot_obj <- gplot_obj +
                plot_trajectory_graph(
                    trajectory_object = trajectory_object,
                    edge_size = trajectory_width,
                    plot_nodes = 0
                )$layers
        }

        gplot_obj + ggplot2::ggtitle(paste("module", gene_module))
    })

    names(plot_list) <- module_names

    legend_plot <- plot_list[[1]] +
        ggplot2::theme(legend.position = "right")
    if (is_discrete_module) {
        legend_plot <- legend_plot +
            ggplot2::guides(
                colour = ggplot2::guide_legend(
                    override.aes = list(size = legend_key_size),
                    title.position = "top",
                    nrow = 2
                )
            )
    } else {
        max_value <- max(sapply(plot_list, function(plot) {
            max(plot$data$cell_info)
        }))

        min_value <- min(sapply(plot_list, function(plot) {
            min(plot$data$cell_info)
        }))
        legend_plot <- legend_plot +
            ggplot2::guides(
                colour = ggplot2::guide_colourbar(
                    title.position = "top",
                    barheight = bar_height
                )
            )
        if (!is.null(continuous_colors)) {
            legend_plot <- legend_plot +
                ggplot2::scale_colour_gradientn(
                    name = paste0("Gene", legend_detail, ifelse(scale_values, "\n(scaled)", "")),
                    colours = continuous_colors,
                    limits = c(min_value, max_value)
                )
        }
    }

    legend_gtable <- ggplot2::ggplotGrob(legend_plot)
    legend_idx <- grep("^guide-box", vapply(legend_gtable$grobs, function(g) g$name, character(1)))

    if (combine_legend) {
        for (i in seq_along(plot_list)) {
            plot_list[[i]] <- plot_list[[i]] + ggplot2::theme(legend.position = "none")
        }
    }

    patchwork_grid <- patchwork::wrap_plots(plotlist = plot_list, ncol = n_columns)
    if (!is_discrete_module) {
        patchwork_grid <- patchwork_grid &
            ggplot2::guides(
                colour = ggplot2::guide_colourbar(
                    title.position = "top",
                    direction = "horizontal"
                )
            )
        if (!is.null(continuous_colors)) {
            patchwork_grid <- patchwork_grid &
                ggplot2::scale_colour_gradientn(
                    name = paste0("Gene", legend_detail, ifelse(scale_values, "\n(scaled)", "")),
                    colours = continuous_colors,
                    limits = c(min_value, max_value)
                )
        }
    }

    if (!combine_legend) {
        return(patchwork_grid)
    }

    if (length(legend_idx) == 0) {
        return(patchwork_grid)
    }

    (
        patchwork_grid |
            patchwork::wrap_elements(full = legend_gtable$grobs[[legend_idx[1]]])
    ) + patchwork::plot_layout(widths = c(1, 0.14))
}

plot_umap_gene_modules_shiny_2 <- function(module_summaries,
                                           shiny_env,
                                           legend_detail = "",
                                           combine_legend = TRUE,
                                           scale_values = TRUE,
                                           n_columns = 3,
                                           trajectory_width = 0.75,
                                           cell_sort_order = c("highest", "lowest", "default"),
                                           cell_size = 0.3,
                                           cell_alpha = 0.8,
                                           legend_text_size = 10,
                                           axis_text_size = 10,
                                           legend_key_size = 1,
                                           colourbar_height = 50,
                                           continuous_colors = NULL) {
    print(paste(Sys.time(), "Generating plots for gene modules..."))
    plot_obj <- plot_umap_gene_modules(
        module_summaries = module_summaries,
        umap_embedding = shiny_env$umap_df,
        trajectory_object = shiny_env$trajectory_object,
        legend_detail = legend_detail,
        combine_legend = combine_legend,
        scale_values = scale_values,
        n_columns = n_columns,
        trajectory_width = trajectory_width,
        cell_sort_order = cell_sort_order,
        cell_size = cell_size,
        cell_alpha = cell_alpha,
        legend_text_size = legend_text_size,
        axis_text_size = axis_text_size,
        legend_key_size = legend_key_size,
        colourbar_height = colourbar_height,
        continuous_colors = continuous_colors
    )
    print(paste(Sys.time(), "Finalised plots for gene modules..."))
    plot_obj
}



#' Generate Cell Heatmap
#'
#' @description Generates a heatmap of gene-family expression across cells,
#' optionally ordered by pseudotime and annotated by metadata.
#'
#' @param expression_matrix A gene-by-cell expression matrix.
#' @param gene_family_list A named list of genes grouped by gene family.
#' @param metadata_df A data frame with cell metadata; rownames must match cell
#' names in `expression_matrix`.
#' @param metadata_name Optional metadata column used for the top annotation.
#' @param legend_name Title used for the heatmap legend.
#' @param apply_scale Logical indicating whether to z-score each gene across
#' cells before plotting.
#' @param k_smooth Integer kernel size used to smooth expression along cells.
#' @param cap Numeric cap applied to scaled or raw expression values.
#' @param axis_text_size Font size used for row labels and annotations.
#' @param discrete_colour_list Optional named list of discrete colour vectors.
#' @param continuous_colors Optional vector of colors used for the heatmap body.
#'
#' @return A `ComplexHeatmap` heatmap object.
#' @export
generate_cell_heatmap <- function(expression_matrix,
                                  gene_family_list,
                                  metadata_df,
                                  metadata_name = NULL,
                                  legend_name = "",
                                  apply_scale = TRUE,
                                  k_smooth = 30,
                                  cap = 20,
                                  axis_text_size = 10,
                                  discrete_colour_list = NULL,
                                  continuous_colors = NULL) {
    if (is.null(metadata_name) || !(metadata_name %in% colnames(metadata_df))) {
        metadata_name <- colnames(metadata_df)[1]
    }
    metadata_df <- metadata_df[colnames(expression_matrix), ]
    mtd_na_mask <- is.na(metadata_df[[metadata_name]])
    expression_matrix <- expression_matrix[, !mtd_na_mask, drop = FALSE]
    metadata_df <- metadata_df[!mtd_na_mask, , drop = FALSE]
    expression_matrix <- expression_matrix[unlist(gene_family_list), , drop = FALSE]

    if (inherits(expression_matrix, "dgCMatrix")) {
        expression_matrix <- as.matrix(expression_matrix)
    }

    if (apply_scale) {
        expression_matrix <- t(scale(t(expression_matrix)))
        gc()
        expression_matrix[expression_matrix < -cap] <- -cap
        expression_matrix[expression_matrix > cap] <- cap
    } else {
        expression_matrix[expression_matrix > cap] <- cap
    }



    if (is.null(continuous_colors)) {
        continuous_colors <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
    }

    ncells <- ncol(expression_matrix)

    # reorder using pseudotime
    present_psd <- FALSE
    if ("pseudotime" %in% colnames(metadata_df)) {
        present_psd <- TRUE
        metadata_df <- metadata_df %>% 
            dplyr::filter(!is.na(.data$pseudotime)) %>%
            dplyr::arrange(.data$pseudotime)
        psd <- metadata_df$pseudotime
        ordered_cells <- rownames(metadata_df)
        expression_matrix <- expression_matrix[ , ordered_cells, drop = FALSE]
    }

    if (k_smooth > 0) {
        weights <- rep(1 / k_smooth, k_smooth)

        matrix_convolved <- matrix(0, nrow = nrow(expression_matrix), ncol = ncol(expression_matrix) + k_smooth - 1)
        for (i in seq_len(nrow(expression_matrix))) {
            matrix_convolved[i, ] <- stats::convolve(expression_matrix[i, ], weights, type = "open")
        }
        matrix_convolved <- matrix_convolved[, -seq_len(k_smooth / 2), drop = FALSE]
        matrix_convolved <- matrix_convolved[, seq_len(ncol(expression_matrix)), drop = FALSE]
        rownames(matrix_convolved) <- rownames(expression_matrix)
        colnames(matrix_convolved) <- colnames(expression_matrix)
        expression_matrix <- matrix_convolved
        rm(matrix_convolved)
        gc()
    }

    gene_family_vector <- rep(NA, nrow(expression_matrix))
    index <- 1
    for (gene_fam in names(gene_family_list)) {
        for (gene in gene_family_list[[gene_fam]]) {
            gene_family_vector[index] <- gene_fam
            index <- index + 1
        }
    }

    nfams <- length(gene_family_list)
    if (is.null(discrete_colour_list) || !(as.character(nfams) %in% names(discrete_colour_list))) {
        if (nfams == 1) {
            cols <- c("#033d03")
        } else {
            cols <- qualpalr::qualpal(n = nfams)$hex
        }
    } else {
        cols <- discrete_colour_list[[as.character(nfams)]]
    }
    names(cols) <- names(gene_family_list)

    left_ha <- ComplexHeatmap::rowAnnotation(
        modules = gene_family_vector,
        col = list(
            modules = cols
        ),
        show_legend = TRUE,
        show_annotation_name = FALSE
    )

    mtd_values <- metadata_df[[metadata_name]]
    unique_mtd <- unique(mtd_values)
    n_unique_mtd <- length(unique_mtd)
    if (is.null(discrete_colour_list) || !(as.character(n_unique_mtd) %in% names(discrete_colour_list))) {
        if (n_unique_mtd == 1) {
            mtd_cols <- c("#033d03")
        } else {
            mtd_cols <- qualpalr::qualpal(n = n_unique_mtd)$hex
        }
    } else {
        mtd_cols <- discrete_colour_list[[as.character(n_unique_mtd)]]
    }
    names(mtd_cols) <- unique_mtd

    arg_list <- list()
    arg_list[[metadata_name]] <- mtd_values
    if (present_psd) {
        psd_colors <- circlize::colorRamp2(c(0, max(psd)), c("black", "#288deb"))
        col_list <- list(pseudotime = psd_colors)
        col_list[[metadata_name]] <- mtd_cols

        arg_list$"pseudotime" <- psd
    } else {
        col_list <- list()
        col_list[[metadata_name]] <- mtd_cols
    }
    arg_list$"col" <- col_list
    arg_list$"show_legend" <- TRUE
    arg_list$"show_annotation_name" <- TRUE
    arg_list$"annotation_name_side" <- "left"

    bottom_ha <- do.call(ComplexHeatmap::HeatmapAnnotation, arg_list)

    remapped <- as.character(
        sapply(seq_along(names(gene_family_list)), function(i) { sprintf("%.2d", i) })
    )
    names(remapped) <- names(gene_family_list)
    args <- c(list(as.character(gene_family_vector)), remapped)
    remapped <- do.call(dplyr::recode, args)

    ComplexHeatmap::Heatmap(
        expression_matrix,
        name = paste0(ifelse(apply_scale, "Scaled ", " "), legend_name),
        row_order = seq_len(nrow(expression_matrix)),
        column_order = seq_len(ncol(expression_matrix)),
        row_split = remapped,
        row_title = names(gene_family_list),
        row_gap = grid::unit(0, "mm"),
        row_names_gp = grid::gpar(fontsize = axis_text_size),
        show_column_names = FALSE,
        show_row_names = FALSE,
        col = continuous_colors,
        bottom_annotation = bottom_ha,
        left_annotation = left_ha,
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(
            direction = "vertical",
            legend_width = grid::unit(5, "cm")
        )
    )
}
