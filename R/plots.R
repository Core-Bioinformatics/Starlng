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

    unique_vals <- unique(umap_df$cell_info)
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
plot_umap <- function(umap_embedding,
                      cell_info,
                      mtd_name,
                      cell_sort_order = c("lowest", "highest", "default"),
                      scale_values = FALSE,
                      cell_size = 0.3,
                      cell_alpha = 0.8,
                      legend_text_size = 10,
                      axis_text_size = 10,
                      colourbar_width = 50,
                      discrete_colors = NULL,
                      continuous_colors = NULL) {
    df <- data.frame(umap_embedding)
    colnames(df) <- c("UMAP1", "UMAP2")
    df$cell_info <- cell_info
    is_discrete <- inherits(df$cell_info, c("factor", "character"))

    if (is_discrete) {
        df <- sort_umap_discrete(df, cell_sort_order)
    } else {
        df <- sort_umap_continuous(df, cell_sort_order)
    }

    if (!is_discrete && scale_values) {
        min_val <- min(df$cell_info)
        max_val <- max(df$cell_info)

        if (max_val > min_val) {
            df$cell_info <- (df$cell_info - min_val) / (max_val - min_val)
            mtd_name <- paste(mtd_name, "(scaled)")
        }
    }

    gplot_obj <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$UMAP1, y = .data$UMAP2, colour = .data$cell_info)) +
        ggplot2::geom_point(size = cell_size, alpha = cell_alpha) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            legend.position = "bottom",
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
                    override.aes = list(size = legend_text_size * 0.6),
                    # FIXME there is a problem with the unique values displayed on legend
                    # if legend text size is too big, the text won't fit anymore
                    # it should create more rows
                    # We should automatically adjust the number of rows based on the width and text size
                    nrow = 3
                )
            )
        if (is.null(discrete_colors)) {
            return(gplot_obj)
        }

        n_unique <- length(unique(df$cell_info))
        used_colors <- discrete_colors[[as.character(n_unique)]]


        return(gplot_obj + ggplot2::scale_color_manual(values = used_colors))
    }

    if (is.null(continuous_colors)) {
        return(gplot_obj)
    }

    return(gplot_obj +
        ggplot2::scale_color_gradientn(colours = continuous_colors) +
        ggplot2::guides(
            colour = ggplot2::guide_colourbar(
                title.position = "top",
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
                           n_coexpressed_threshold,
                           summarise_expr,
                           scale_values = FALSE,
                           trajectory_gplot = NULL,
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
        legend_name <- "Cells above threshold"
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

    gplot_obj <- plot_umap(
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
    )

    if (is.null(trajectory_gplot)) {
        return(gplot_obj)
    }

    traj_layer <- trajectory_gplot$layers[[1]]
    traj_layer$aes_params$size <- trajectory_width 
    gplot_obj$layers <- c(gplot_obj$layers, traj_layer)
    return(gplot_obj)
}

# assumes the shiny context with the environment env
plot_umap_gene_shiny <- function(shiny_env,
                                 gene_name,
                                 thresh_percentile,
                                 thresh_value,
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
        gene_name = gene_name,
        matrix_h5_path = file.path("objects", "expression.h5"),
        index_genes = shiny_env$genes[gene_name],
        check_intersect = FALSE
    )

    plot_umap_gene(
        umap_embedding = shiny_env$umap_df,
        gene_matrix = gene_matrix,
        gene_name = gene_name,
        thresh_percentile = thresh_percentile,
        thresh_value = thresh_value,
        scale_values = scale_values,
        n_coexpressed_threshold = n_coexpressed_threshold,
        summarise_expr = summarise_expr,
        trajectory_gplot = shiny_env$trajectory_gplot,
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

    for (i in seq_along(plot_list)) {
        plot_list[[i]] <- plot_list[[i]] +
            ggplot2::theme(
                legend.position = "right"
            ) +
            ggplot2::guides(
                colour = ggplot2::guide_colourbar(
                    title.position = "top",
                    barheight = grid::unit(colourbar_height, "points")
                )
            )
    }

    patchwork::wrap_plots(plotlist = plot_list, ncol = n_columns) +
        patchwork::plot_layout(guides = "collect") &
        ggplot2::scale_colour_gradientn(
            name = paste0("Gene\n", summarise_expr, "\n", ifelse(scale_values, "(scaled)", "")),
            colours = continuous_colors,
            limits = c(min_value, max_value)
        )
}


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
    if (apply_scale) {
        expression_matrix <- t(scale(t(expression_matrix)))
        expression_matrix[expression_matrix < -cap] <- -cap
        expression_matrix[expression_matrix > cap] <- cap
    } else {
        expression_matrix[expression_matrix > cap] <- cap
    }

    if (is.null(metadata_name) || !(metadata_name %in% colnames(metadata_df))) {
        metadata_name <- colnames(metadata_df)[1]
    }
    metadata_df <- metadata_df[colnames(expression_matrix), ]

    if (is.null(continuous_colors)) {
        continuous_colors <- c("white", "red")
    }

    # reorder using pseudotime
    present_psd <- FALSE
    if ("pseudotime" %in% colnames(metadata_df)) {
        present_psd <- TRUE
        metadata_df <- metadata_df %>% 
            dplyr::filter(!is.na(.data$pseudotime)) %>%
            dplyr::arrange(.data$pseudotime)
        psd <- metadata_df$pseudotime
        ordered_cells <- rownames(metadata_df)
        expression_matrix <- expression_matrix[ , ordered_cells]
    }

    if (k_smooth > 0) {
        weights <- rep(1 / k_smooth, k_smooth)

        matrix_convolved <- matrix(0, nrow = nrow(expression_matrix), ncol = ncol(expression_matrix) + k_smooth - 1)
        for (i in seq_len(nrow(expression_matrix))) {
            matrix_convolved[i, ] <- convolve(expression_matrix[i, ], weights, type = "open")
        }
        matrix_convolved <- matrix_convolved[, -seq_len(k_smooth / 2), drop = FALSE]
        matrix_convolved <- matrix_convolved[, seq_len(ncol(expression_matrix)), drop = FALSE]
        rownames(matrix_convolved) <- rownames(expression_matrix)
        colnames(matrix_convolved) <- colnames(expression_matrix)
        expression_matrix <- matrix_convolved
        rm(matrix_convolved)
        gc()
    }


    # sample_names <- rep(NA, ncol(expression_matrix))
    # sample_names[seq_len(ncells1)] <- get("dts_names", envir = env)[1]
    # sample_names[(ncells1 + 1):ncol(expression_matrix)] <- get("dts_names", envir = env)[2]

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
            cols <- qualpalr::qualpal(n = nfams, colorspace = "pretty_dark")$hex
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
            mtd_cols <- qualpalr::qualpal(n = n_unique_mtd, colorspace = "pretty_dark")$hex
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
        use_raster = TRUE,
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
