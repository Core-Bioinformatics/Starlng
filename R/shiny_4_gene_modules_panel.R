#' @importFrom dplyr %>% .data
NULL

###### UI ######
ui_module_table <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::splitLayout(
            shiny::sliderInput(
                inputId = ns("top_cells_percent"),
                label = "Top cells (%)",
                min = 0, max = 100, value = 100, step = 0.1
            ),
            shiny::sliderInput(
                inputId = ns("scale_threshold"),
                label = "Scale threshold",
                min = 0, max = 1, value = 0.5, step = 0.01
            ),
            shinyWidgets::radioGroupButtons(
                inputId = ns("summarise_expr"),
                label = "How to summarise?",
                choices = c("average", "average expressed", "#genes")
            ),
            cellWidths = c("25%", "25%", "50%")
        ),
        DT::DTOutput(ns("module_stats"), height = "auto", width = "90%")
    )
}

ui_module_pseudotime_expression <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::plotOutput(ns("expression_pseudotime_plot"), height = "auto"),
        shiny::plotOutput(ns("pseudotime_violin_plot"), height = "auto")
    )
}

ui_grid_umaps <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::splitLayout(
            cellWidths = c("40px", "40px", "140px"),
            gear_umaps(ns, "settings", "highest", TRUE, TRUE),
            gear_download(ns, "module_umap", "modules"),
            shiny::numericInput(
                inputId = ns("n_columns"),
                label = "Number of columns",
                min = 1, max = 20, value = 0, step = 1
            )
        ),
        shiny::plotOutput(ns("module_umap_plot"), height = "auto")
    )
}

ui_module_heatmap <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tagList(
        shiny::h2("Module relationship heatmap", id = "gene_modules_heatmap_panel"),
        shiny::splitLayout(
            cellWidths = "40px",
            shinyWidgets::dropdownButton(
                inputId = ns("heatmap_settings"),
                shiny::tagList(
                    shiny::numericInput(
                        inputId = ns("colour_clipping"),
                        label = "Colour clipping value",
                        min = 0, value = 500
                    ),
                    shiny::sliderInput(
                        inputId = ns("text_size"),
                        label = "Text size",
                        min = 0.00, max = 30.00, value = 15, step = 0.05
                    ),
                    shiny::sliderInput(
                        inputId = ns("axis_size"),
                        label = "Axis labels size",
                        min = 0.00, max = 30.00, value = 15, step = 0.05
                    ),
                    shiny::sliderInput(
                        inputId = ns("legend_size"),
                        label = "Legend size",
                        min = 0.00, max = 30.00, value = 15, step = 0.05
                    ),
                    shiny::sliderInput(
                        inputId = ns("dend_height"),
                        label = "Dendogram height (mm)",
                        min = 1, max = 500, value = 70, step = 1
                    ),
                    shiny::selectInput(
                        inputId = ns("colour_scheme"),
                        label = "Colour scheme",
                        choices = c(""),
                        selected = NULL
                    )
                ),
                label = "",
                icon = shiny::icon("cog"),
                status = "success",
                size = "sm",
                circle = TRUE
            ),
            gear_download(ns, "heatmap", "heatmap")
        ),
        shiny::splitLayout(
            shinyWidgets::radioGroupButtons(
                inputId = ns("heatmap_info"),
                label = "Displayed Information",
                choices = c("Expressed cells", "Expressed cells (JSI)", "Spearman correlation"),
                selected = "Expressed cells"
            )
        ),
        shiny::plotOutput(ns("modules_heatmap"), height = "auto")
    )
}

ui_gene_hub_scores <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Gene Hub Scores", id = "gene_hub_scores_panel"),
        DT::DTOutput(ns("gene_hub_scores_table"), height = "auto", width = "90%"),
        # download as csv
        shiny::downloadButton(ns("download_gene_hub_scores"), "Download gene hub scores")
    )
}

ui_gene_umap_hubs <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Gene UMAP", id = "gene_umap_panel"),
        shiny::splitLayout(
            cellWidths = c("40px", "40px", "440px"),
            shinyWidgets::dropdownButton(
                inputId = ns("gear_umap_hubs"),
                shiny::splitLayout(
                    shiny::tagList(
                        shiny::sliderInput(
                            inputId = ns("text_size"),
                            label = "Text size",
                            min = 5.00, max = 30.00, value = 6.00, step = 0.5
                        ),
                        shiny::sliderInput(
                            inputId = ns("axis_size"),
                            label = "Axis labels size",
                            min = 5.00, max = 30.00, value = 10.00, step = 0.5
                        ),
                        shiny::sliderInput(
                            inputId = ns("legend_size"),
                            label = "Legend size",
                            min = 5.00, max = 30.00, value = 10.00, step = 0.5
                        ),
                        shiny::sliderInput(
                            inputId = ns("pt_alpha"),
                            label = "Colour alpha",
                            min = 0.00, max = 1.00, value = 0.8, step = 0.05
                        )
                    ),
                    shiny::tagList(
                        shiny::sliderInput(
                            inputId = ns("pt_size"),
                            label = "Point size",
                            min = 0.05, max = 20.00, value = 2.00, step = 0.1
                        ),
                        shiny::sliderInput(
                            inputId = ns("edge_range_width"),
                            label = "Edge width range",
                            min = 0.00, max = 5.00, value = c(0.05, 0.50), step = 0.05
                        ),
                        shiny::sliderInput(
                            inputId = ns("top_edges"),
                            label = "Top edges (%)",
                            min = 0, max = 100, value = 10, step = 0.1
                        ),
                        shiny::sliderInput(
                            inputId = ns("percent_nodes"),
                            label = "Percent nodes (Non-hubs)",
                            min = 0, max = 100, value = 100, step = 0.1
                        )
                    )
                ),
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "gene_umap_hubs", "gene_umap_hubs"),
            shiny::splitLayout(
                shiny::selectizeInput(
                    inputId = ns("select_modules"),
                    label = "Select modules",
                    choices = NULL,
                    selected = NULL,
                    multiple = TRUE,
                    options = list(
                        placeholder = "Select modules",
                        plugins = list("remove_button", "drag_drop"),
                        delimiter = ",",
                        create = TRUE
                    )
                ),
                shiny::numericInput(
                    inputId = ns("n_hub_genes"),
                    label = "Number of hub genes",
                    min = 1, max = 20, value = 2, step = 1
                )
            )
        ),
        shiny::plotOutput(ns("gene_umap_hubs_plot"), height = "auto"),
        shiny::plotOutput(ns("graph_module_transition"), height = "auto")
    )
}

#' UI - Gene Module UMAP
#'
#' @description Creates the UI interface for the Gene Module UMAP panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
ui_module_umap <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Visualisation of module expression", id = "gene_modules_panel"),
        shiny::splitLayout(
            ui_module_table(ns("gene_modules_table")),
            shiny::tagList(
                ui_module_pseudotime_expression(ns("gene_modules_pseudotime")),
                shiny::splitLayout(
                    shiny::selectizeInput(
                        inputId = ns("select_modules"),
                        label = "Select modules",
                        choices = NULL,
                        selected = NULL,
                        multiple = TRUE,
                        options = list(
                            placeholder = "Select modules",
                            plugins = list("remove_button", "drag_drop"),
                            delimiter = ",",
                            create = TRUE
                    )),
                    shiny::actionButton(
                        inputId = ns("order_modules"),
                        label = "Order modules",
                        icon = shiny::icon("sort-numeric-up")
                    ),
                    shiny::actionButton(
                        inputId = ns("reset_modules"),
                        label = "Reset modules!",
                        icon = shiny::icon("undo")
                    ),
                    shiny::uiOutput(ns("clip"))
                )
            ),
            cellWidths = c("50%", "45%")
        ),
        ui_grid_umaps(ns("gene_modules_umap")),
        ui_module_heatmap(ns("gene_modules_heatmap")),
        ui_gene_hub_scores(ns("gene_hub_scores")),
        ui_gene_umap_hubs(ns("gene_umap_hubs"))
    )
}

###### UTILITIES - SERVER #######
util_load_module_summaries <- function(shiny_env, gene_modules, summ_expr, clustering_type) {
    summary_function <- switch(
        summ_expr,
        "binary" = NULL,
        "average expressed" = function(x) { mean(x[x > 0]) },
        "average" = mean,
        "#genes" = function(x) { sum(x > 0) },
        NULL
    )
    n_modules <- length(gene_modules)

    condition_loading <- clustering_type == "preloaded"

    if (condition_loading) {
        print(paste(Sys.time(), "Loading precomputed module summaries from HDF5 file."))
        has_loaded <- FALSE
        try({
            module_names <- as.character(rhdf5::h5read(
                file = file.path("objects", "module_summaries.h5"),
                name = paste0(n_modules, "/modules")
            ))
            summ_list <- rhdf5::h5read(
                file = file.path("objects", "module_summaries.h5"),
                name = paste0(n_modules, "/expression_summaries"),
                index = list(NULL, NULL)
            )
            colnames(summ_list) <- module_names
            rownames(summ_list) <- shiny_env$cells
            summ_list <- split(summ_list, col(summ_list))
            print(paste(Sys.time(), "Finished loading precomputed module summaries."))
            has_loaded <- TRUE
        }, silent = TRUE)
        if (has_loaded) {
            return(summ_list)
        }
        warning("Tried to load precomputed module summaries from HDF5 file but they were not found. Recomputing...")
    } 
    summ_list <- lapply(names(gene_modules), function(gene_module) {
        gc()
        print(paste(Sys.time(), "Processing module:", gene_module))
        gene_matrix <- read_gene_from_dense_h5(
            gene_names = gene_modules[[gene_module]],
            matrix_h5_path = file.path("objects", "expression.h5"),
            index_genes = shiny_env$genes[gene_modules[[gene_module]]],
            check_intersect = FALSE
        )

        return(voting_scheme(
            expression_matrix = gene_matrix,
            genes = gene_modules[[gene_module]],
            thresh_percentile = 0,
            thresh_value = 0,
            n_coexpressed_thresh = 1,
            summary_function = summary_function
        ))
    })
    print(paste(Sys.time(), "Finished processing modules."))
    names(summ_list) <- names(gene_modules)

    return(summ_list)
}

util_load_gene_hubs <- function(shiny_env, gene_modules, module_mask, psd_mask, summ_expr, scale_threshold, module_ordering = NULL, condition_loading) {
    n_modules <- nrow(module_mask)
    has_loaded <- FALSE
    try({
        score_df <- rhdf5::h5read(
            file = file.path("objects", "module_summaries.h5"),
            name = paste0(n_modules, "/gene_hub_stats")
        )
        rownames(score_df) <- as.character(rhdf5::h5read(
            file = file.path("objects", "module_summaries.h5"),
            name = "genes"
        ))
        for (i in seq_len(5)) {
            score_df[,i] <- as.numeric(score_df[,i])
        }
        score_df$module <- as.character(score_df$module)
        has_loaded <- TRUE
    }, silent = TRUE)

    if (!has_loaded) {
        warning("Tried to load precomputed gene hub scores from HDF5 file but they were not found. Recomputing...")
        condition_loading <- FALSE
        module_ordering <- rownames(module_mask)
        total_weight <- get_per_module_weight(shiny_env$chosen_graph(), gene_modules)
    } else {
        module_ordering <- stats::setNames(score_df$module, rownames(score_df))
        total_weight <- stats::setNames(score_df$total_weight, rownames(score_df))
    }

    if (!condition_loading) {
        score_df <- do.call(rbind, lapply(rownames(module_mask), function(module_name) {
            get_gene_overlap_stat(
                expr_matrix = read_gene_from_dense_h5(
                    gene_names = gene_modules[[module_name]],
                    matrix_h5_path = file.path("objects", "expression.h5"),
                    index_genes = shiny_env$genes[gene_modules[[module_name]]],
                    check_intersect = FALSE
                )[, psd_mask, drop = FALSE],
                gene_modules = gene_modules[[module_name]],
                module_expr = module_mask[module_name, ],
                gene_expression_percentile = 0.5,
                scale = FALSE,
                total_weight = total_weight
            )
        }))
        score_df$module <- module_ordering[rownames(score_df)]
    }

    score_df$module <- factor(score_df$module, levels = as.character(seq_len(n_modules)))
    return(score_df %>% dplyr::arrange(.data$module, dplyr::desc(.data$combined_score)))
}

util_load_closest_module <- function(shiny_env, module_mask, condition_loading) {
    if (condition_loading) {
        has_loaded <- FALSE

        try({
            closest_node_per_module <- stats::setNames(
                as.character(rhdf5::h5read(
                    file = file.path("objects", "module_summaries.h5"),
                    name = paste0(nrow(module_mask), "/closest_nodes_to_module")
                )),
                as.character(rhdf5::h5read(
                    file = file.path("objects", "module_summaries.h5"),
                    name = paste0(nrow(module_mask), "/modules")
                ))
            )
            has_loaded <- TRUE
        }, silent = TRUE)
        if (has_loaded) {
            return(closest_node_per_module)
        }
        warning("Tried to load precomputed closest nodes to modules from HDF5 file but they were not found. Recomputing...")
    }
    closest_node_per_module <- get_closest_node_to_module(
        trajectory_object = shiny_env$trajectory_object,
        cell_umap = shiny_env$umap_df,
        module_expr = lapply(rownames(module_mask), function(module_name) {
                module_mask[module_name, ]
            }) %>% stats::setNames(rownames(module_mask)),
    )
    return(closest_node_per_module)
}

util_load_pairwise_tables <- function(module_mask, module_summ, condition_loading) {
    if (condition_loading) {
        file_exists <- FALSE
        try({
            matrix_list <- list(
                module_intersect_cells = rhdf5::h5read(
                    file = file.path("objects", "module_summaries.h5"),
                    name = paste0(nrow(module_mask), "/module_intersect_cells")
                ),
                module_union_cells = rhdf5::h5read(
                    file = file.path("objects", "module_summaries.h5"),
                    name = paste0(nrow(module_mask), "/module_union_cells")
                ),
                module_spearman_matrix = rhdf5::h5read(
                    file = file.path("objects", "module_summaries.h5"),
                    name = paste0(nrow(module_mask), "/module_spearman_matrix")
                )
            )
            file_exists <- TRUE
        }, silent = TRUE)
        if (file_exists) {
            return(matrix_list)
        }
        warning("Tried to load precomputed pairwise tables from HDF5 file but they were not found. Recomputing...")
    }
    pairwise_tables <- compute_module_pairwise_tables(
        module_mask = module_mask,
        module_summ = module_summ
    )
    return(pairwise_tables)
}

util_load_other_table_objects <- function(object_name, object_params, processing_function, condition_loading) {
    if (condition_loading) {
        module_object <- NULL
        try({
            module_object <- rhdf5::h5read(
                file = file.path("objects", "module_summaries.h5"),
                name = object_name
            )
        }, silent = TRUE)
        if (!is.null(module_object)) {
            return(module_object)
        }
        warning(paste("Tried to load", object_name, "from precomputed file but it was not found. Recomputing..."))
    }

    module_object <- do.call(
        processing_function,
        object_params
    )

    return(module_object)
}


###### SERVER ######
server_module_table_prepare <- function(id, parent_session) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            thresholds_reactive <- shiny::reactive({
                list(
                    top_cells_percent = input$top_cells_percent,
                    scale_threshold = input$scale_threshold
                )
            }) %>% shiny::debounce(3000)

            # module summaries
            shiny::observe({
                gene_modules <- env$chosen_modules()
                summ_expr <- input$summarise_expr
                clustering_type <- env$clustering_type()

                shiny::isolate({
                    shiny::req(gene_modules, summ_expr, length(gene_modules) > 0)
                    summ_list <- util_load_module_summaries(env, gene_modules, summ_expr, clustering_type)

                    env$modules_summaries(summ_list)

                    scaled_summ <- lapply(summ_list, scale_min_max)
                    names(scaled_summ) <- names(summ_list)

                    env$modules_summaries_scaled(scaled_summ)
                })
            })

            shiny::observe({
                module_summ <- env$modules_summaries_scaled()
                thresh_vals <- thresholds_reactive()
                psd_value <- env$psd_value()
                clustering_type <- env$clustering_type()
                summ_expr <- input$summarise_expr
                gene_modules <- env$chosen_modules()

                shiny::isolate({
                    shiny::req(!is.null(module_summ), psd_value, clustering_type, summ_expr, gene_modules)

                    psd_mask <- !is.na(psd_value)
                    condition_loading <- clustering_type == "preloaded" &
                        rhdf5::h5read(file.path("objects", "module_summaries.h5"), "summary_method") == summ_expr &
                        rhdf5::h5read(file.path("objects", "module_summaries.h5"), "scale_threshold") == thresh_vals$scale_threshold &
                        sum(psd_mask) == length(module_summ[[1]])
                    module_mask <- build_module_masks(
                        module_summ = module_summ,
                        psd_mask = psd_mask,
                        scale_threshold = thresh_vals$scale_threshold,
                        top_cells_percent = thresh_vals$top_cells_percent
                    )
                    env$modules_mask(module_mask)
  
                    pairwise_tables <- util_load_pairwise_tables(module_mask, module_summ, condition_loading)

                    env$modules_table_cells(pairwise_tables$module_intersect_cells)
                    env$modules_union_cells(pairwise_tables$module_union_cells)
                    env$modules_table_spearman(pairwise_tables$module_spearman_matrix)

                    modules_stats <- util_load_other_table_objects(
                        object_name = paste0(length(module_summ), "/modules_stats"),
                        object_params = list(
                            module_summ = module_summ,
                            module_mask = module_mask,
                            psd_value = psd_value,
                            umap_df = env$umap_df
                        ),
                        processing_function = get_module_stats,
                        condition_loading = condition_loading
                    )
                    
                    modules_stats_summary <- util_load_other_table_objects(
                        object_name = paste0(length(module_summ), "/modules_stats_summary"),
                        object_params = list(
                            modules_stats = modules_stats,
                            gene_modules = gene_modules
                        ),
                        processing_function = summarise_module_stats,
                        condition_loading = condition_loading
                    )

                    modules_stats_summary <- util_load_other_table_objects(
                        object_name = paste0(length(module_summ), "/modules_stats_summary"),
                        object_params = list(
                            modules_stats_summary = modules_stats_summary,
                            module_mask = module_mask,
                            psd_value = psd_value,
                            umap_dist_threshold = 0.85 * env$umap_median_distance

                        ),
                        processing_function = annotate_module_outliers,
                        condition_loading = condition_loading
                    )

                    env$modules_stats(modules_stats)
                    env$modules_stats_summary(modules_stats_summary)

                    shiny::req(!is.null(modules_stats), !is.null(modules_stats_summary), nrow(modules_stats_summary) > 0)

                    env$gene_hub_scores(util_load_gene_hubs(
                        shiny_env = env,
                        gene_modules = gene_modules,
                        module_mask = module_mask,
                        psd_mask = psd_mask,
                        summ_expr = summ_expr,
                        scale_threshold = thresh_vals$scale_threshold,
                        module_ordering = modules_stats_summary$module,
                        condition_loading = condition_loading
                    ))
                    env$closest_node_per_module(util_load_closest_module(
                        shiny_env = env,
                        module_mask = module_mask,
                        condition_loading = condition_loading
                    ))
                })
            })
        }
    )
}

server_module_table <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::observe({
                modules_stats <- env$modules_stats_summary()
                shiny::isolate({
                    shiny::req(modules_stats, nrow(modules_stats) > 0, cancelOutput = TRUE)
                    output$module_stats <- DT::renderDT({
                        DT::datatable(
                            modules_stats,
                            options = list(
                                pageLength = nrow(modules_stats),
                                columnDefs = list(
                                    list(
                                        className = "dt-right",
                                        targets = "_all"
                                    )
                                )
                            ),
                            colnames = c("module", "n_genes", "n cells", "avg<br>summary", "median<br>psd", "iqr<br>psd", "median<br>UMAP dist", "outlier?"),
                            escape = FALSE,
                            rownames = FALSE
                        ) %>% DT::formatStyle(
                            "is_outlier",
                            target = "row",
                            backgroundColor = DT::styleEqual(
                                levels = c("yes", "redundant", "no"),
                                values = c("transparent", "rgba(252, 186, 3, 0.2)", "rgba(0, 255, 0, 0.2)")
                            )
                        )
                    })
                })
            }) %>% shiny::bindEvent(env$modules_stats_summary())
        }
    )
}

server_module_pseudotime_expression <- function(id, module_ordering) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::observe({
                psd_df <- env$modules_stats()
                wdim <- env$window_dim()
                md_order <- module_ordering()
                summary_stat_df <- env$modules_stats_summary()

                shiny::isolate({
                    shiny::req(psd_df, md_order, cancelOutput = TRUE)
                    plt_height <- wdim[2] * 0.5
                    plt_width <- wdim[1] / 0.45
                    summary_stat_df$module <- as.character(summary_stat_df$module)
                    psd_df$module <- as.character(psd_df$module)
                    rownames(summary_stat_df) <- summary_stat_df$module
                    unique_groups <- sort(rownames(summary_stat_df))
                    summary_stat_df <- summary_stat_df[md_order, ]
                    shiny::req(all(md_order %in% unique_groups))

                    psd_df <- psd_df %>% dplyr::filter(.data$module %in% md_order)

                    mdl_iqrs <- summary_stat_df %>%
                        dplyr::filter(.data$module %in% md_order) %>%
                        dplyr::pull(.data$iqr_pseudotime)
                    names(mdl_iqrs) <- md_order
                    low_iqrs <- names(mdl_iqrs[mdl_iqrs < 0.1])
                    big_iqrs <- setdiff(md_order, low_iqrs)
                    if (length(low_iqrs) > 0) {
                        plot_df_small <- psd_df %>%
                            dplyr::filter(.data$module %in% low_iqrs) %>%
                            dplyr::group_by(.data$module) %>%
                            dplyr::filter(.data$psd_value >= stats::quantile(.data$psd_value, 0.25) & .data$psd_value <= stats::quantile(.data$psd_value, 0.75))
                    } else {
                        plot_df_small <- NULL
                    }

                    if (length(big_iqrs) > 0) {
                        plot_df_big <- psd_df %>% dplyr::filter(.data$module %in% big_iqrs)
                    } else {
                        plot_df_big <- NULL
                    }

                    output$expression_pseudotime_plot <- shiny::renderPlot(
                        height = plt_height,
                        {
                            discrete_colours <- env$color_options$discrete[[as.character(length(unique_groups))]]
                            if (is.null(discrete_colours)) {
                                discrete_colours <- qualpalr::qualpal(length(unique_groups))$hex
                            }
                            names(discrete_colours) <- unique_groups
                            gplot_obj <- ggplot2::ggplot()
                            if (!is.null(plot_df_small)) {
                                gplot_obj <- gplot_obj +
                                    ggplot2::geom_line(
                                        data = plot_df_small,
                                        ggplot2::aes(x = .data$psd_value, y = .data$avg_summary, color = .data$module, group = .data$module),
                                        linewidth = 1.0
                                    )
                            }
                            if (!is.null(plot_df_big)) {
                                gplot_obj <- gplot_obj +
                                    ggplot2::geom_smooth(
                                        data = plot_df_big,
                                        ggplot2::aes(x = .data$psd_value, y = .data$avg_summary, color = .data$module, group = .data$module),
                                        method = "loess",
                                        se = FALSE,
                                        linewidth = 1.0
                                    )
                            }
                            gplot_obj +
                                ggplot2::scale_color_manual(values = discrete_colours) +
                                ggplot2::theme_classic() +
                                ggplot2::labs(x = "Pseudotime", y = "Summary") +
                                ggplot2::theme(
                                    legend.position = "bottom",
                                    axis.text = ggplot2::element_text(size = 16),
                                    legend.text = ggplot2::element_text(size = 16),
                                    legend.title = ggplot2::element_blank(),
                                    axis.title = ggplot2::element_text(size = 16)
                                ) +
                                ggplot2::ylim(0, 1.1)
                    })
                })
            })

            shiny::observe({
                psd_df <- env$modules_stats()
                wdim <- env$window_dim()
                mod_ord <- module_ordering()

                shiny::isolate({
                    shiny::req(wdim, cancelOutput = TRUE)
                    output$pseudotime_violin_plot <- shiny::renderPlot(
                        height = wdim[2] * 0.5,
                        {
                            shiny::req(psd_df, mod_ord, cancelOutput = TRUE)
                            nunique <- sort(rownames(env$modules_stats_summary()))
                            psd_df <- psd_df %>% dplyr::filter(.data$module %in% mod_ord)
                            psd_df$module <- factor(psd_df$module, levels = mod_ord)
                            discrete_colours <- env$color_options$discrete[[as.character(length(nunique))]]
                            if (is.null(discrete_colours)) {
                                discrete_colours <- qualpalr::qualpal(length(nunique))$hex
                            }
                            names(discrete_colours) <- nunique
                            ggplot2::ggplot(psd_df, ggplot2::aes(x = .data$module, y = .data$psd_value, fill = .data$module)) +
                                ggplot2::geom_violin(scale = "width") +
                                ggplot2::scale_fill_manual(values = discrete_colours) +
                                ggplot2::theme_classic() +
                                ggplot2::labs(x = "Module", y = "Pseudotime") +
                                ggplot2::theme(
                                    legend.position = "none",
                                    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 16),
                                    axis.text.y = ggplot2::element_text(size = 16),
                                    axis.title.x = ggplot2::element_text(size = 16),
                                    axis.title.y = ggplot2::element_text(size = 16)
                                )
                        }
                    )

                })

            })
        }
    )
}

server_grid_umaps <- function(id, module_ordering) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            default_col_scheme <- "white_red"
            if (!(default_col_scheme %in% names(env$color_options$continuous))) {
                default_col_scheme <- names(env$color_options$continuous)[1]
            }

            shiny::updateSelectInput(
                session,
                inputId = "settings_colour_scheme",
                choices = names(env$color_options$continuous),
                selected = default_col_scheme
            )

            shiny::observe({
                shiny::req(env$chosen_modules())
                print(paste(Sys.time(), "Update number of columns"))
                shiny::updateNumericInput(
                    session,
                    inputId = "n_columns",
                    value = min(6, length(env$chosen_modules()))
                )
            })

            composite_plot <- shiny::reactive({
                print(paste(Sys.time(), "Creating composite plot for gene modules."))
                win_dim <- env$window_dim()
                # module_summ <- env$modules_summaries()
                module_summ <- env$modules_summaries_scaled()
                patch_ncol <- input$n_columns
                summarise_expr <- input$summarise_expr
                module_mask <- env$modules_mask()
                selected_modules <- module_ordering()

                req_gear_umap(session, "settings")

                shiny::isolate({
                    shiny::req(patch_ncol > 0, !is.null(module_summ), !is.null(win_dim))
                    shiny::req(selected_modules, all(selected_modules %in% names(module_summ)))
                    gc()
                    plot_rows <- ceiling(length(env$chosen_modules()) / patch_ncol)
                    plt_height <- win_dim[2] / 1.05
                    plt_width <- win_dim[1] / 1.05
                    panel_width <- plt_width / patch_ncol
                    height_from_panels <- panel_width * plot_rows * 1.1
                    plt_height <- min(plt_height, height_from_panels)
                    colourbar_height <- plt_height * 0.75

                    for (i in seq_along(module_summ)) {
                        if (!(names(module_summ)[i] %in% selected_modules)) {
                            next
                        }
                        module_summ[[i]][!module_mask[i, ]] <- 0
                    }

                    pl_res <- plot_umap_gene_modules_shiny_2(
                        module_summaries = module_summ[selected_modules],
                        shiny_env = env,
                        legend_detail = summarise_expr,
                        n_columns = patch_ncol,
                        # scale_values = input$settings_scale,
                        scale_values = FALSE,
                        trajectory_width = input$settings_trajectory_width,
                        cell_sort_order = input$settings_pt_order,
                        cell_size = input$settings_pt_size,
                        cell_alpha = input$settings_pt_alpha,
                        legend_text_size = input$settings_legend_size,
                        axis_text_size = input$settings_axis_size,
                        colourbar_height = colourbar_height,
                        continuous_colors = env$color_options$continuous[[input$settings_colour_scheme]]
                    )
                    pl_res
                })
            })

            output$module_umap_plot <- shiny::renderPlot(
                height = function() {
                    shiny::req(env$window_dim(), input$n_columns > 0, env$chosen_modules())
                    plot_cols <- input$n_columns
                    plot_rows <- ceiling(length(env$chosen_modules()) / plot_cols)
                    plt_height <- env$window_dim()[2] / 1.05
                    plt_width <- env$window_dim()[1] / 1.05
                    panel_width <- plt_width / plot_cols
                    height_from_panels <- panel_width * plot_rows * 1.1
                    min(plt_height, height_from_panels)
                },
                width = function() {
                    shiny::req(env$window_dim())
                    env$window_dim()[1] / 1.05
                },
                {
                    cmp_plot <- composite_plot()
                    shiny::req(cmp_plot)
                    print(paste(Sys.time(), "Rendering module UMAP plot."))
                    print(cmp_plot)
                }
            )

            output$download_module_umap <- shiny::downloadHandler(
                filename = function() {
                    req_gear_download(session, "module_umap")
                    paste(input$filename_module_umap, tolower(input$filetype_module_umap), sep = ".")
                },
                content = function(file) {
                    req_gear_download(session, "module_umap")
                    cmp_plot <- composite_plot()
                    shiny::req(cmp_plot)
                    device_fun <- save_filetypes[[input$filetype_module_umap]]
                    shiny::req(!is.null(device_fun))
                    device_fun(file, width = input$width_module_umap, height = input$height_module_umap)
                    on.exit(grDevices::dev.off(), add = TRUE)
                    print(cmp_plot)
                }
            )
        }
    )
}

server_module_heatmap <- function(id, module_ordering) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            default_col_scheme <- "blue_red"
            if (!(default_col_scheme %in% names(env$color_options$continuous))) {
                default_col_scheme <- names(env$color_options$continuous)[1]
            }
            shiny::updateSelectInput(
                session,
                inputId = "colour_scheme",
                choices = names(env$color_options$continuous),
                selected = default_col_scheme
            )

            hclust_dends_ncell <- shiny::reactive({
                mod_ord <- module_ordering()
                shiny::req(mod_ord)
                dist_matrix <- env$modules_table_cells()
                union_tab <- env$modules_union_cells()
                shiny::req(dist_matrix, union_tab, all(mod_ord %in% rownames(dist_matrix)), identical(rownames(dist_matrix), rownames(union_tab), length(mod_ord) > 1))
                shiny::req(nrow(dist_matrix) > 1, nrow(union_tab) > 1)
                for (i in seq_len(nrow(dist_matrix))) {
                    union_tab[i, i] <- dist_matrix[i, i]
                }

                dist_matrix <- 1 - dist_matrix[mod_ord, mod_ord, drop = FALSE] / union_tab[mod_ord, mod_ord, drop = FALSE]

                shiny::req(nrow(dist_matrix) > 1)

                for (i in mod_ord) {
                    for (j in mod_ord) {
                        if (union_tab[i, j] == 0) {
                            dist_matrix[i, j] <- 1
                        }
                    }
                }
                return(stats::hclust(stats::as.dist(dist_matrix), method = "average"))
            })

            hclust_dends_spearman <- shiny::reactive({
                mod_ord <- module_ordering()
                shiny::req(mod_ord)
                dist_matrix <- env$modules_table_spearman()
                shiny::req(dist_matrix, all(mod_ord %in% rownames(dist_matrix)))
                shiny::req(nrow(dist_matrix) > 1, length(mod_ord) > 1)

                dist_matrix <- 1 - dist_matrix[mod_ord, mod_ord, drop = FALSE]
                dist_matrix[is.na(dist_matrix)] <- 100
                shiny::req(nrow(dist_matrix) > 1)
                return(stats::hclust(stats::as.dist(dist_matrix), method = "average"))
            })

            heatmap_obj <- shiny::reactive({
                mod_ord <- module_ordering()
                htmp_type <- input$heatmap_info
                is_spearman <- htmp_type == "Spearman correlation"
                union_tab <- env$modules_union_cells()

                if (!is_spearman) {
                    htmp_matrix <- env$modules_table_cells()
                    hclust_dend <- hclust_dends_ncell()
                } else {
                    htmp_matrix <- env$modules_table_spearman()
                    hclust_dend <- hclust_dends_spearman()
                }

                col_scheme <- input$colour_scheme
                axis_size <- input$axis_size
                text_size <- input$text_size
                legend_size <- input$legend_size
                clipping_val <- input$colour_clipping
                dend_height <- input$dend_height

                shiny::isolate({
                    shiny::req(htmp_matrix, !is.null(htmp_matrix), nrow(htmp_matrix) > 1, hclust_dend, union_tab)
                    shiny::req(mod_ord, all(mod_ord %in% rownames(htmp_matrix)))
                    htmp_matrix <- htmp_matrix[mod_ord, mod_ord, drop = FALSE]

                    col_scheme <- env$color_options$continuous[[col_scheme]]
                    if (!is_spearman) {
                        htmp_matrix[htmp_matrix > clipping_val] <- clipping_val
                    } else {
                        col_scheme <- circlize::colorRamp2(seq(from = -1, to = 1, length.out = length(col_scheme)), col_scheme) 
                    }
                    gc()

                    if (htmp_type == "Expressed cells (JSI)") {
                        union_tab <- union_tab[mod_ord, mod_ord, drop = FALSE]
                        for (i in seq_len(nrow(htmp_matrix))) {
                            union_tab[i, i] <- htmp_matrix[i, i]
                        }
                        htmp_matrix <- round(htmp_matrix / union_tab, 2)
                        for (i in seq_len(nrow(htmp_matrix))) {
                            htmp_matrix[i, i] <- NA
                        }
                    }

                    if (min(htmp_matrix, na.rm = TRUE) == max(htmp_matrix, na.rm = TRUE)) {
                        col_scheme <- "gray"
                    }

                    ComplexHeatmap::Heatmap(
                        matrix = htmp_matrix,
                        name = htmp_type,
                        col = col_scheme,
                        cluster_rows = hclust_dend,
                        row_dend_width = grid::unit(dend_height, "mm"),
                        cluster_columns = hclust_dend,
                        show_row_dend = TRUE,
                        show_column_dend = FALSE,
                        row_names_side = "left",
                        column_names_side = "top",
                        row_names_gp = grid::gpar(fontsize = axis_size),
                        column_names_gp = grid::gpar(fontsize = axis_size),
                        column_title_gp = grid::gpar(fontsize = axis_size),
                        row_title_gp = grid::gpar(fontsize = axis_size),
                        heatmap_legend_param = list(
                            title_gp = grid::gpar(fontsize = legend_size),
                            labels_gp = grid::gpar(fontsize = legend_size),
                            legend_height = grid::unit(legend_size * 2, "mm")
                        ),
                        width = ncol(htmp_matrix) * (text_size / 2 + 2),
                        height = nrow(htmp_matrix) * (text_size / 2 + 2),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                            if (text_size > 0.1 && !is.na(htmp_matrix[i, j])) {
                                grid::grid.text(htmp_matrix[i, j], x = x, y = y, gp = grid::gpar(fontsize = text_size * 0.8))
                            }
                        },
                        na_col = "lightgrey"
                    )
                })
            })

            shiny::observe({
                htmp_obj <- heatmap_obj()
                req_gear_download(session, "heatmap")

                output$modules_heatmap <- shiny::renderPlot(
                    height = env$window_dim()[2] / 1.05,
                    width = env$window_dim()[1] / 1.05,
                    {
                        shiny::req(htmp_obj)
                        ComplexHeatmap::draw(htmp_obj, heatmap_legend_side = "right", annotation_legend_side = "right")
                    }
                )

                output$download_heatmap <- shiny::downloadHandler(
                    filename = function() {
                        paste(input$filename_heatmap, tolower(input$filetype_heatmap), sep = ".")
                    },
                    content = function(file) {
                        shiny::req(htmp_obj)
                        device_fun <- save_filetypes[[input$filetype_heatmap]]
                        shiny::req(!is.null(device_fun))
                        device_fun(file, width = input$width_heatmap, height = input$height_heatmap)
                        on.exit(grDevices::dev.off(), add = TRUE)
                        ComplexHeatmap::draw(htmp_obj, heatmap_legend_side = "right", annotation_legend_side = "right")
                    }
                        
                )

            })
            
        }
    )
}

server_gene_hub_scores <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::observe({
                score_df <- env$gene_hub_scores()
                shiny::isolate({
                    shiny::req(score_df, cancelOutput = TRUE)
                    score_df <- as.data.frame(score_df)
                    shiny::req(nrow(score_df) > 0, ncol(score_df) > 0, cancelOutput = TRUE)

                    output$gene_hub_scores_table <- DT::renderDT({
                        DT::datatable(
                            score_df,
                            filter = "top",
                            rownames = TRUE 
                        )
                    })

                    output$download_gene_hub_scores <- shiny::downloadHandler(
                        filename = function() {
                            "gene_hub_scores.csv"
                        },
                        content = function(file) {
                            score_df <- env$gene_hub_scores()
                            shiny::req(score_df, nrow(score_df) > 0)
                            utils::write.table(score_df, file, sep = ",", quote = FALSE, row.names = TRUE)
                        }
                    )
                })
            }) %>% shiny::bindEvent(env$gene_hub_scores())
        }
    )
}

server_gene_umap_hubs <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            module_ordering <- shiny::reactive({
                input$select_modules
            }) %>% shiny::debounce(3000)
            module_adjacency <- shiny::reactive({
                closest_node <- env$closest_node_per_module()
                used_modules <- module_ordering()

                shiny::isolate({
                    shiny::req(closest_node, used_modules, all(used_modules %in% names(closest_node)), cancelOutput = TRUE)
                    closest_node <- closest_node[used_modules]

                    return(
                        get_module_transitions(
                            trajectory_object = env$trajectory_object,
                            closest_module = closest_node
                        )
                    )
                })
            })

            shiny::observe({
                mod_adj <- module_adjacency()
                closest_node <- env$closest_node_per_module()
                text_size <- input$text_size
                point_size <- input$pt_size
                proposed_order <- module_ordering()
                shiny::isolate({
                    shiny::req(mod_adj, closest_node, proposed_order)
                    output$graph_module_transition <- shiny::renderPlot(
                        height = env$window_dim()[2] / 1.05,
                        width = env$window_dim()[1] / 1.05,
                        {
                            plot_module_transitions(
                                module_adj_matrix = mod_adj,
                                closest_module = closest_node,
                                node_positions = env$trajectory_object$node_positions,
                                node_label_size = text_size,
                                node_size = point_size,
                                start_module = proposed_order[1]
                            )
                        }
                    )
                })
            })

            hub_genes <- shiny::reactive({
                used_modules <- module_ordering()
                n_genes <- input$n_hub_genes
                score_df <- env$gene_hub_scores()

                shiny::isolate({
                    shiny::req(used_modules, n_genes, score_df)
                    score_df$gene <- rownames(score_df)
                    score_df <- score_df %>%
                        dplyr::filter(.data$module %in% used_modules) %>%
                        dplyr::group_by(.data$module) %>%
                        dplyr::slice_max(order_by = .data$combined_score, n = n_genes) %>%
                        dplyr::ungroup()
                    return(score_df)
                })
            })

            filtered_module_adjacency <- shiny::reactive({
                module_adj_matrix <- module_adjacency()
                gene_adj_matrix <- env$chosen_graph()
                current_hubs <- hub_genes()
                perc_nodes <- input$percent_nodes
                perc_edges <- input$top_edges
                chosen_modules <- env$chosen_modules()

                shiny::isolate({
                    shiny::req(module_adj_matrix, gene_adj_matrix, current_hubs, perc_nodes, perc_edges, chosen_modules)
                    used_modules <- rownames(module_adj_matrix)
                    return(get_filtered_gene_adjacency(
                        gene_modules = chosen_modules[used_modules],
                        module_adjacency = module_adj_matrix,
                        gene_adjacency = gene_adj_matrix,
                        hub_genes = current_hubs,
                        percentage_non_hub_nodes = perc_nodes / 100,
                        percentage_edges = perc_edges / 100
                    ))
                })
            })

            shiny::observe({
                fltr_adj <- filtered_module_adjacency()
                edge_range <- input$edge_range_width
                edge_alpha <- input$pt_alpha
                point_size <- input$pt_size
                legend_size <- input$legend_size
                text_size <- input$text_size
                umap_df <- env$chosen_umap()
                shiny::isolate({
                    shiny::req(fltr_adj, edge_range, umap_df)
                    gplot_obj <- plot_gene_hub_umap(
                        umap_df = umap_df,
                        filtered_gene_adj = fltr_adj,
                        edge_weight_range = edge_range,
                        edge_alpha = edge_alpha,
                        point_size = point_size,
                        legend_text_size = legend_size,
                        node_text_size = text_size
                    )
                    output$gene_umap_hubs_plot <- shiny::renderPlot(
                        height = env$window_dim()[2] / 1.05,
                        width = env$window_dim()[1] / 1.05,
                        {
                            shiny::req(gplot_obj)
                            print(gplot_obj)
                        }
                    )

                    output$download_gene_umap_hubs <- shiny::downloadHandler(
                        filename = function() {
                            paste(input$filename_gene_umap_hubs, tolower(input$filetype_gene_umap_hubs), sep = ".")
                        },
                        content = function(file) {
                            ggplot2::ggsave(
                                filename = file,
                                plot = gplot_obj,
                                width = input$width_gene_umap_hubs,
                                height = input$height_gene_umap_hubs,
                                units = "in",
                                dpi = 300
                            )
                        }
                    )
                })
            })
        }
    )
}
#' Server - Gene Module UMAP
#'
#' @description Creates the backend interface for the Gene Module UMAP panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
server_module_umap <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            server_module_table_prepare("gene_modules_table", session)
            server_module_table("gene_modules_table")

            shiny::observe({
                modules_stats <- env$modules_stats_summary()
                shiny::isolate({
                    shiny::req(modules_stats, nrow(modules_stats) > 0)
                    shiny::updateSelectizeInput(
                        session,
                        inputId = "select_modules",
                        choices = modules_stats$module,
                        selected = modules_stats$module[modules_stats$is_outlier == "no"],
                        server = TRUE
                    )

                    shiny::updateSelectizeInput(
                        session,
                        inputId = "gene_umap_hubs-select_modules",
                        choices = modules_stats$module,
                        selected = modules_stats$module[modules_stats$is_outlier == "no"],
                        server = TRUE
                    )
                })
            }) %>% shiny::bindEvent(env$modules_stats_summary())

            output$clip <- shiny::renderUI({
                rclipboard::rclipButton(
                    inputId = "clipbtn",
                    label = "Copy to clipboard",
                    clipText = paste(input$select_modules, collapse = ","),
                    icon = shiny::icon("clipboard"),
                    class = "btn btn-primary",
                    options = list(delay = list(show = 800, hide = 100), trigger = "hover")
                )
            })

            shiny::observe({
                modules_stats <- env$modules_stats_summary()
                shiny::updateSelectizeInput(
                    session,
                    inputId = "select_modules",
                    selected = modules_stats$module[modules_stats$is_outlier == "no"]
                )
            }) %>% shiny::bindEvent(input$reset_modules)

            shiny::observe({
                modules_stats <- env$modules_stats_summary()
                slct_modules <- input$select_modules
                shiny::req(slct_modules)
                shiny::updateSelectizeInput(
                    session,
                    inputId = "select_modules",
                    selected = intersect(modules_stats$module, slct_modules)
                )
            }) %>% shiny::bindEvent(input$order_modules)

            module_ordering <- shiny::reactive({
                input$select_modules
            }) %>% shiny::debounce(3000)

            shiny::observeEvent(shiny::req(module_ordering()), env$trigger4(TRUE))

            server_module_pseudotime_expression("gene_modules_pseudotime", module_ordering)
            server_grid_umaps("gene_modules_umap", module_ordering)
            server_module_heatmap("gene_modules_heatmap", module_ordering)
            server_gene_hub_scores("gene_hub_scores")
            server_gene_umap_hubs("gene_umap_hubs")
        }
    )
}
