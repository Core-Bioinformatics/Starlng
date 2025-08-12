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
                min = 0, max = 100, value = 75, step = 0.1
            ),
            shiny::sliderInput(
                inputId = ns("scale_threshold"),
                label = "Scale threshold",
                min = 0, max = 1, value = 0.4, step = 0.01
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
                shiny::selectizeInput(
                    inputId = ns("select_modules"),
                    label = "Select modules",
                    choices = NULL,
                    selected = NULL,
                    multiple = TRUE,
                    options = list(
                        placeholder = "Select modules",
                        plugins = list("remove_button", "drag_drop"),
                        delimiter = ","
                ))
            ),
            cellWidths = c("50%", "45%")
        ),
        ui_grid_umaps(ns("gene_modules_umap"))
    )
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

            # read and save the matrices
            shiny::observe({
                gene_modules <- env$chosen_modules()
                summ_expr <- input$summarise_expr

                shiny::isolate({
                    shiny::req(gene_modules, summ_expr, length(gene_modules) > 0)
                    summary_function <- switch(
                        summ_expr,
                        "binary" = NULL,
                        "average expressed" = function(x) { mean(x[x > 0]) },
                        "average" = mean,
                        "#genes" = function(x) { sum(x > 0) },
                        NULL
                    )

                    summ_list <- lapply(names(gene_modules), function(gene_module) {
                        gc()
                        print(paste(Sys.time(), "Processing module:", gene_module))
                        gene_matrix <- read_gene_from_dense_h5(
                            gene_names = gene_modules[[gene_module]],
                            matrix_h5_path = file.path("objects", "expression.h5"),
                            index_genes = env$genes[gene_modules[[gene_module]]],
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

                    env$modules_summaries(summ_list)

                    scaled_summ <- lapply(summ_list, scale_min_max)
                    names(scaled_summ) <- names(summ_list)
                    env$modules_summaries_scaled(scaled_summ)
                })
            })

            # 
            shiny::observe({
                module_summ <- env$modules_summaries_scaled()
                thresh_vals <- thresholds_reactive()
                psd_value <- env$psd_value()

                shiny::isolate({
                    perc_cells <- 1 - thresh_vals$top_cells_percent / 100
                    scale_threshold <- thresh_vals$scale_threshold
                    shiny::req(perc_cells >= 0, perc_cells <= 1, scale_threshold >= 0, scale_threshold <= 1, !is.null(module_summ), psd_value)

                    mask_nrow <- length(module_summ)
                    mask_ncol <- length(psd_value)
                    module_mask <- matrix(rep(FALSE, mask_nrow * mask_ncol), nrow = mask_nrow, ncol = mask_ncol)
                    module_cells_table <- matrix(0, nrow = mask_nrow, ncol = mask_nrow)

                    rownames(module_mask) <- names(module_summ)
                    rownames(module_cells_table) <- names(module_summ)
                    colnames(module_cells_table) <- names(module_summ)


                    for (i in seq_along(module_summ)) {
                        module_mask[i, ] <- module_summ[[i]] > scale_threshold
                        value_thresh <- quantile(module_summ[[i]][module_mask[i, ]], perc_cells)
                        mx_val <- max(module_summ[[i]])
                        if (value_thresh == mx_val) {
                            value_thresh <- 0.99 * mx_val
                        }
                        module_mask[i, ] <- module_summ[[i]] > value_thresh
                    }

                    for (i in seq(from = 1, to = nrow(module_cells_table) - 1)) {
                        for (j in seq(from = i + 1, to = ncol(module_cells_table))) {
                            ncommon_cells <- sum(module_mask[i, ] & module_mask[j, ])
                            module_cells_table[i, j] <- ncommon_cells
                            module_cells_table[j, i] <- ncommon_cells
                        }
                    }
                    for (i in seq_len(nrow(module_cells_table))) {
                        unique_cells <- module_mask[i, ] 
                        for (j in seq(from = 1, to = ncol(module_cells_table))) {
                            if (i == j) {
                                next
                            }
                            unique_cells <- unique_cells & !module_mask[j, ]
                        }
                        module_cells_table[i, i] <- sum(unique_cells)
                    }

                    env$modules_mask(module_mask)

                    modules_stats <- NULL
                    modules_stats_summary <- NULL
                    
                    for (i in seq_along(module_summ)) {
                        temp_df <- data.frame(
                            avg_expression = module_summ[[i]][module_mask[i, ]],
                            psd_value = psd_value[module_mask[i, ]],
                            umap_distance = calculate_umap_average_distance(
                                umap_df = env$umap_df,
                                selected_cells = which(module_mask[i, ])
                            )
                        )
                        temp_df$module <- names(module_summ)[i]

                        temp_df_summary <- data.frame(
                            module = names(module_summ)[i],
                            n_cells = sum(module_mask[i, ]),
                            avg_expression = round(mean(temp_df$avg_expression, na.rm = TRUE), 3),
                            median_pseudotime = round(median(temp_df$psd_value, na.rm = TRUE), 3),
                            iqr_pseudotime = round(IQR(temp_df$psd_value, na.rm = TRUE), 3),
                            median_umap_distance = round(median(temp_df$umap_distance, na.rm = TRUE), 3)
                        )

                        if (is.null(modules_stats)) {
                            modules_stats <- temp_df
                            modules_stats_summary <- temp_df_summary
                        } else {
                            modules_stats <- rbind(modules_stats, temp_df)
                            modules_stats_summary <- rbind(modules_stats_summary, temp_df_summary)
                        }
                    }
                    rownames(modules_stats_summary) <- modules_stats_summary$module

                    modules_stats_summary$is_outlier <- mad_z_score(modules_stats_summary$iqr_pseudotime) > 3.5

                    i <- 1
                    while (i <= nrow(modules_stats_summary) - 1) {
                        diff_first <- abs(modules_stats_summary$median_pseudotime[i] - modules_stats_summary$median_pseudotime[i + 1])
                        if (diff_first < 0.1) {
                            j <- i + 1
                            while (j <= nrow(modules_stats_summary)) {
                                diff_next <- abs(modules_stats_summary$median_pseudotime[j] - modules_stats_summary$median_pseudotime[i])
                                if (diff_next > 0.1) {
                                    j <- j - 1
                                    break
                                }
                                j <- j + 1
                            }
                            if (j > nrow(modules_stats_summary)) {
                                j <- nrow(modules_stats_summary)
                            }

                            z_scores <- mad_z_score(modules_stats_summary$iqr_pseudotime[i:j]) > 3.5
                            modules_stats_summary$is_outlier[i:j] <- modules_stats_summary$is_outlier[i:j] | z_scores

                            i <- j + 1
                        } else {
                            i <- i + 1
                        }
                    }
                    modules_stats_summary <- modules_stats_summary %>%
                        dplyr::arrange(.data$median_pseudotime, .data$iqr_pseudotime, .data$median_umap_distance)

                    env$modules_stats(modules_stats)
                    env$modules_stats_summary(modules_stats_summary)
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
                                autoWidth = TRUE,
                                columnDefs = list(
                                    list(
                                        className = "dt-right",
                                        targets = "_all"
                                    )
                                )
                            ),
                            rownames = FALSE
                        ) %>% DT::formatStyle(
                            "is_outlier",
                            target = "row",
                            backgroundColor = DT::styleEqual(
                                c(TRUE, FALSE),
                                c("rgba(255, 0, 0, 0.2)", "transparent")
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

                shiny::isolate({
                    shiny::req(psd_df, md_order, cancelOutput = TRUE)
                    plt_height <- wdim[2] * 0.5
                    plt_width <- wdim[1] / 0.45
                    psd_df <- psd_df %>% dplyr::filter(.data$module %in% md_order)
                    summary_stat_df <- env$modules_stats_summary()
                    unique_groups <- sort(rownames(summary_stat_df))
                    summary_stat_df <- summary_stat_df[md_order, ]

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
                            dplyr::filter(.data$psd_value >= quantile(.data$psd_value, 0.25) & .data$psd_value <= quantile(.data$psd_value, 0.75))
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
                                        ggplot2::aes(x = .data$psd_value, y = .data$avg_expression, color = .data$module, group = .data$module),
                                        linewidth = 1.0
                                    )
                            }
                            if (!is.null(plot_df_big)) {
                                gplot_obj <- gplot_obj +
                                    ggplot2::geom_smooth(
                                        data = plot_df_big,
                                        ggplot2::aes(x = .data$psd_value, y = .data$avg_expression, color = .data$module, group = .data$module),
                                        method = "loess",
                                        se = FALSE,
                                        linewidth = 1.0
                                    )
                            }
                            gplot_obj +
                                ggplot2::scale_color_manual(values = discrete_colours) +
                                ggplot2::theme_classic() +
                                ggplot2::labs(x = "Pseudotime", y = "Expression") +
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
                trigger_val <- env$trigger4()

                req_gear_umap(session, "settings")

                shiny::isolate({
                    shiny::req(trigger_val)
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

            shiny::observe({
                cmp_plot <- composite_plot()
                wdim <- env$window_dim()

                shiny::isolate({
                    shiny::req(cmp_plot, !is.null(cmp_plot), !is.null(wdim))
                    print(paste(Sys.time(), "Rendering module UMAP plot."))
                    plot_cols <- input$n_columns
                    plot_rows <- ceiling(length(env$chosen_modules()) / plot_cols)
                    plt_height <- env$window_dim()[2] / 1.05
                    plt_width <- env$window_dim()[1] / 1.05
                    panel_width <- plt_width / plot_cols
                    height_from_panels <- panel_width * plot_rows * 1.1
                    plt_height <- min(plt_height, height_from_panels)
                    env$trigger4(FALSE)

                    output$module_umap_plot <- shiny::renderPlot(
                        height = plt_height,
                        width = plt_width,
                        cmp_plot
                    )
                })
            })

            shiny::observe({
                shiny::req(composite_plot(), input$filename_module_umap, input$filetype_module_umap, input$width_module_umap, input$height_module_umap)

                shiny::isolate({
                    output$download_module_umap <- shiny::downloadHandler(
                        filename = function() {
                            paste(input$filename_module_umap, tolower(input$filetype_module_umap), sep = ".")
                        },
                        content = function(file) {
                            ggplot2::ggsave(
                                filename = file,
                                plot = composite_plot(),
                                width = input$width_module_umap,
                                height = input$height_module_umap
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
                shiny::isolate({
                    modules_stats <- env$modules_stats_summary()
                    shiny::req(modules_stats, nrow(modules_stats) > 0)
                    shiny::updateSelectizeInput(
                        session,
                        inputId = "select_modules",
                        choices = modules_stats$module,
                        selected = modules_stats$module[!modules_stats$is_outlier],
                        server = TRUE
                    )
                    print("Changed trigger to TRUE")
                })
            }) %>% shiny::bindEvent(env$modules_stats_summary())

            module_ordering <- shiny::reactive({
                input$select_modules
            }) %>% shiny::debounce(4000)

            shiny::observeEvent(shiny::req(module_ordering()), env$trigger4(TRUE))

            server_module_pseudotime_expression("gene_modules_pseudotime", module_ordering)
            server_grid_umaps("gene_modules_umap", module_ordering)
        }
    )
}
