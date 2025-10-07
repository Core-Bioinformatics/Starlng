#' @importFrom dplyr %>% .data
NULL

###### UI ######
ui_pseudotime_table <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        DT::dataTableOutput(ns("mtd_psd_table"))
    )
}

ui_pseudotime_mtd_order <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::plotOutput(ns("psd_boxplot"), height = "auto")
    )
}

ui_pseudobulk_heatmap <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h3("Pseudobulk heatmap"),
        shiny::splitLayout(
            cellWidths = c("40px", "40px", "300px"),
            shinyWidgets::dropdownButton(
                inputId = ns("heatmap_settings"),
                shiny::tagList(
                    shinyWidgets::radioGroupButtons(
                        inputId = ns("summarise_expr"),
                        label = "How to summarise?",
                        choices = c("average", "average expressed", "#genes")
                    ),
                    shinyWidgets::prettySwitch(
                        inputId = ns("scale_values"),
                        label = "Scale values",
                        status = "success",
                        value = TRUE,
                        fill = TRUE
                    ),
                    shiny::numericInput(
                        inputId = ns("cap_value"),
                        label = "Clipping value",
                        min = 0, max = 20, value = 2, step = 0.1
                    ),
                    shiny::sliderInput(
                        inputId = ns("text_size"),
                        label = "Text size",
                        min = 0.00, max = 30.00, value = 12, step = 0.05
                    ),
                    shiny::sliderInput(
                        inputId = ns("point_size"),
                        label = "Point size",
                        min = 0.10, max = 30, value = c(0.5, 10), step = 0.01
                    ),
                    shiny::sliderInput(
                        inputId = ns("axis_size"),
                        label = "Axis labels size",
                        min = 0.00, max = 30.00, value = 12, step = 0.05
                    ),
                    shiny::sliderInput(
                        inputId = ns("legend_size"),
                        label = "Legend size",
                        min = 0.00, max = 30.00, value = 12, step = 0.05
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
            gear_download(ns, "heatmap", "heatmap"),
            shinyWidgets::radioGroupButtons(
                inputId = ns("plot_type"),
                label = "Select plot type",
                choices = c("heatmap", "dotplot"),
                selected = "dotplot"
            )
        ),
        shiny::plotOutput(ns("pseudobulk_heatmap"), height = "auto")
    )
}

ui_by_cell_heatmap <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h3("Cell-resolution heatmap"),
        shiny::fluidRow(
            shiny::column(2,
                shiny::splitLayout(
                    cellWidths = "40px",
                    shinyWidgets::dropdownButton(
                        inputId = ns("heatmap_settings"),
                        shiny::tagList(
                            shinyWidgets::prettySwitch(
                                inputId = ns("scale_values"),
                                label = "Scale values",
                                status = "success",
                                value = TRUE,
                                fill = TRUE
                            ),
                            shiny::numericInput(
                                inputId = ns("k_smooth"),
                                label = "Smoothing degree",
                                min = 0, max = 100, value = 30, step = 1
                            ),
                            shiny::numericInput(
                                inputId = ns("cap_value"),
                                label = "Value cap",
                                min = 0, max = 20, value = 2, step = 0.1
                            ),
                            shiny::sliderInput(
                                inputId = ns("axis_size"),
                                label = "Axis labels size",
                                min = 0.00, max = 30.00, value = 12, step = 0.05
                            ),
                            shiny::sliderInput(
                                inputId = ns("legend_size"),
                                label = "Legend size",
                                min = 0.00, max = 30.00, value = 12, step = 0.05
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
                )
            ),
            shiny::p("Make sure the pseudobulk heatmap is generated first before pressing the button!"),
            shiny::p("Warning! Creating this heatmap might require a lot of memory and time. Please be patient. Consider using a small number of genes."),

            shiny::column(2,
                shiny::actionButton(
                    inputId = ns("generate_heatmap"),
                    label = "Generate heatmap",
                    class = "btn-danger"
                )
            )
        ),
        shiny::plotOutput(ns("by_cell_heatmap"), height = "600px", width = "1400px")#, height = "auto")
    )
}

#' UI - Gene Module Heatmap
#'
#' @description Creates the UI interface for the Gene Module Heatmap panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
ui_module_metadata_heatmap <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Order metadata by pseudotime", id = "pseudotime_panel"),
        ui_pseudotime_table(ns("psd_table")),
        shiny::selectInput(
            inputId = ns("metadata"),
            label = "Select metadata",
            choices = c(""),
            selected = NULL
        ),
        ui_pseudotime_mtd_order(ns("mtd_ordering")),
        shiny::h2("Gene modules heatmaps", id = "heatmaps_panel"),
        ui_pseudobulk_heatmap(ns("pseudobulk")),
        # shiny::actionButton(
        #     inputId = ns("size_env"),
        #     label = "Check env size",
        #     class = "btn-info"
        # ),
        shiny::fluidRow(
            shiny::column(5,
                shiny::selectizeInput(
                    inputId = ns("module_order"),
                    label = "Predefined order of modules:",
                    width = "100%",
                    choices = c(""),
                    selected = NULL,
                    multiple = TRUE,
                    options = list(
                        plugins = list("remove_button", "drag_drop"),
                        delimiter = ",",
                        create = TRUE
                    )
                )
            ),
            shiny::column(1,
                shiny::actionButton(
                    inputId = ns("order_modules"),
                    label = "Order modules!",
                    icon = shiny::icon("sort-numeric-up")
                )
            ),
            shiny::column(1,
                shiny::actionButton(
                    inputId = ns("reset_modules"),
                    label = "Reset modules!",
                    icon = shiny::icon("undo")
                )
            ),
            shiny::column(1,
                shiny::uiOutput(ns("clip"))
            )
        ),
        ui_by_cell_heatmap(ns("by_cell"))
    )
}

###### SERVER ######
server_pseudotime_table <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::observe({
                psd_val <- env$psd_value()
                psd_ordering <- list()

                shiny::isolate({
                    shiny::req(psd_val)
                    na_mask <- is.na(psd_val)
                    are_na <- any(na_mask)
                    na_mask <- !na_mask
                    if (are_na) {
                        psd_val <- psd_val[na_mask]
                    }
                    min_psd <- min(psd_val, na.rm = TRUE)
                    max_psd <- max(psd_val, na.rm = TRUE)

                    for (dscr_mtd in names(env$discrete_mtd)) {
                        if (are_na) {
                            mtd_grps <- split(psd_val, env$mtd_df[na_mask, dscr_mtd])
                            # next
                        } else {
                            mtd_grps <- split(psd_val, env$mtd_df[, dscr_mtd])
                        }

                        if (length(mtd_grps) == 0) {
                            next
                        }
                
                        psd_ordering[[dscr_mtd]] <- list()
                        psd_ordering[[dscr_mtd]]$stats <- lapply(mtd_grps, function(x) {
                            if (length(x) == 0) {
                                return(NULL)
                            }
                            qs_vec <- stats::quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
                            return(c(qs_vec, mean(x), qs_vec[4] - qs_vec[2]))
                        })
                        # names(psd_ordering[[dscr_mtd]]$stats) <- names(mtd_grps)
                        psd_ordering[[dscr_mtd]]$stats <- as.data.frame(do.call(rbind, psd_ordering[[dscr_mtd]]$stats))
                        colnames(psd_ordering[[dscr_mtd]]$stats) <- c("q05", "Q1", "median", "Q3", "q95", "mean", "IQR")

                        psd_ordering[[dscr_mtd]]$stats <- psd_ordering[[dscr_mtd]]$stats %>% dplyr::arrange(.data$median)

                        psd_ordering[[dscr_mtd]]$coverage <- calculate_pseudotime_iqr_coverage(psd_ordering[[dscr_mtd]]$stats) / (max_psd - min_psd) * 100

                        psd_ordering[[dscr_mtd]]$iqr_stats <- stats::fivenum(psd_ordering[[dscr_mtd]]$stats[, "IQR"])
                    }

                    env$psd_ordering(psd_ordering)
                })
            }) %>% shiny::bindEvent(env$psd_value())

            shiny::observe({
                psd_order <- env$psd_ordering()

                shiny::isolate({
                    shiny::req(length(psd_order) > 0)
                    group_stats_df <- lapply(psd_order, function(x) {
                                c(x$coverage, x$iqr_stats)
                            })
                    names(group_stats_df) <- names(psd_order)
                    group_stats_df <- as.data.frame(do.call(rbind, group_stats_df))
                    colnames(group_stats_df) <- c("Coverage", "min IQR", "Q1 IQR", "median IQR", "Q3 IQR", "max IQR")
                    q_coverage <- max(group_stats_df$Coverage) * 0.85
                    group_stats_df <- rbind(
                        group_stats_df %>% dplyr::filter(.data$Coverage > q_coverage) %>% dplyr::arrange(.data$"median IQR"),
                        group_stats_df %>% dplyr::filter(.data$Coverage <= q_coverage) %>% dplyr::arrange(dplyr::desc(.data$Coverage))
                    )

                    output$mtd_psd_table <- DT::renderDataTable(
                        {
                            DT::datatable(
                                group_stats_df,
                                rownames = TRUE
                            )
                        },
                        server = FALSE
                    )

                    best_option <- rownames(group_stats_df)[1]
                    shiny::updateSelectInput(
                        session,
                        inputId = "select_mtd",
                        choices = rownames(group_stats_df),
                        selected = best_option
                    )
                })
            }) %>% shiny::bindEvent(env$psd_ordering())
        }
    )
}

server_pseudotime_mtd_order <- function(id, select_mtd) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::observe({
                group_mtd <- select_mtd()
                psd_order <- env$psd_ordering()
                wdim <- env$window_dim()

                shiny::isolate({
                    shiny::req(wdim)
                    shiny::req(length(psd_order) > 0)
                    shiny::req(group_mtd %in% names(psd_order))
                    shiny::req(psd_order[[group_mtd]]$stats)
                    psd_order <- psd_order[[group_mtd]]$stats

                    plt_height <- floor(env$height_ratio * wdim[2])
                    plt_width <- wdim[1] * 0.9
                    

                    output$psd_boxplot <- shiny::renderPlot(
                        height = plt_height,
                        width = plt_width,
                        {
                        shiny::req(psd_order)

                        psd_order$x <- as.character(rownames(psd_order))
                        psd_order$x <- factor(psd_order$x, levels = psd_order$x)

                        ggplot2::ggplot(psd_order, ggplot2::aes(x = .data$x)) +
                            ggplot2::geom_boxplot(
                                ggplot2::aes(ymin = .data$q05, lower = .data$Q1, middle = .data$median,
                                    upper = .data$Q3, ymax = .data$q95),
                                stat = "identity"
                            ) +
                            ggplot2::labs(x = group_mtd, y = "Pseudotime") +
                            ggplot2::theme_bw() +
                            ggplot2::theme(
                                axis.text = ggplot2::element_text(size = 12),
                                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
                            )
                    })
            
                })

            })
        }
    )
}

server_pseudobulk_heatmap <- function(id, metadata_name, module_ordering) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::updateSelectInput(
                session,
                inputId = "colour_scheme",
                choices = names(env$color_options$continuous),
                selected = names(env$color_options$continuous)[length(env$color_options$continuous)]
            )

            psdbulk_matrix <- shiny::reactive({
                current_modules <- env$chosen_modules()
                mtd_name <- metadata_name()
                summarise_expr <- input$summarise_expr
                psd_ordering <- env$psd_ordering()
                # module_ordering <- module_ordering()
                shiny::isolate({
                    shiny::req(mtd_name, summarise_expr)
                    shiny::req(psd_ordering, current_modules)
                    summary_function <- switch(
                        summarise_expr,
                        "binary" = NULL,
                        "average expressed" = function(x) { mean(x[x > 0]) },
                        "average" = mean,
                        "#genes" = function(x) { sum(x > 0) },
                        NULL
                    )

                    pseudobulk_matrix <- sapply(
                        current_modules,
                        function(module) {
                            cell_info <- voting_scheme(
                                read_gene_from_dense_h5(
                                    gene_names = module,
                                    matrix_h5_path = file.path("objects", "expression.h5"),
                                    index_genes = env$genes[module],
                                    check_intersect = FALSE
                                ),
                                genes = module,
                                thresh_percentile = 0,
                                thresh_value = 0,
                                n_coexpressed_thresh = 0,
                                summary_function = summary_function
                            )

                            pseudobulk_summary(
                                mtd_value = env$mtd_df[, mtd_name],
                                cell_info = cell_info
                            )
                        }
                    )
                    pseudobulk_matrix <- matrix(pseudobulk_matrix, ncol = length(current_modules))

                    colnames(pseudobulk_matrix) <- paste0("module ", names(current_modules))
                    rownames(pseudobulk_matrix) <- env$discrete_mtd[[mtd_name]]

                    mtd_ordering <- rownames(psd_ordering[[mtd_name]]$stats)

                    return(pseudobulk_matrix[mtd_ordering, , drop = FALSE])
                })
            })

            percentages_mat <- shiny::reactive({
                ptype <- input$plot_type
                current_modules <- module_ordering()
                mtd_name <- metadata_name()
                summarise_expr <- input$summarise_expr
                psd_ordering <- env$psd_ordering()
                shiny::isolate({
                    shiny::req(ptype == "dotplot")
                    shiny::req(mtd_name, summarise_expr)
                    shiny::req(mtd_name %in% colnames(env$mtd_df))
                    shiny::req(psd_ordering, current_modules)

                    mod_mask <- env$modules_mask()

                    index_mtd <- split(seq_len(nrow(env$mtd_df)), env$mtd_df[, mtd_name])
                    mtd_ordering <- rownames(psd_ordering[[mtd_name]]$stats)
                    # percentages <- matrix(0, nrow = length(index_mtd), ncol = length(current_modules))
                    percentages <- matrix(0, nrow = length(mtd_ordering), ncol = length(current_modules))
                    colnames(percentages) <- paste0("module ", current_modules)

                    rownames(percentages) <- mtd_ordering

                    for (i in mtd_ordering) { #names(index_mtd)) {
                        ncells <- length(index_mtd[[i]])
                        for (j in current_modules) {
                            ncom <- sum(mod_mask[j, index_mtd[[i]]])
                            percentages[i, paste0("module ", j)] <- ncom / ncells
                        }
                    }

                    return(percentages)
                })
            })

            htmp_obj <- shiny::reactive({
                ptype <- input$plot_type
                pseudobulk_matrix <- psdbulk_matrix()
                module_ord <- module_ordering()
                req_gear_heatmap(session)

                shiny::isolate({
                    shiny::req(ptype == "heatmap")
                    shiny::req(!is.null(pseudobulk_matrix), module_ord)
                    shiny::req(all(paste("module", module_ord) %in% colnames(pseudobulk_matrix)))
                    module_ord <- paste0("module ", module_ord)
                    pseudobulk_matrix <- pseudobulk_matrix[, module_ord, drop = FALSE]
                    gc()

                    if (input$scale_values && ncol(pseudobulk_matrix) > 1) {
                        pseudobulk_matrix <- t(scale(t(pseudobulk_matrix)))
                        pseudobulk_matrix[pseudobulk_matrix > input$cap_value] <- input$cap_value
                        pseudobulk_matrix[pseudobulk_matrix < -input$cap_value] <- -input$cap_value
                    } else {
                        pseudobulk_matrix[pseudobulk_matrix > input$cap_value] <- input$cap_value
                    }

                    legend_name <- paste0(
                        input$summarise_expr,
                        ifelse(startsWith(input$summarise_expr, "average"), "\nexpression", ""),
                        ifelse(input$scale_values, "\n(z-scaled)", "")
                    )

                    ComplexHeatmap::Heatmap(
                        pseudobulk_matrix,
                        row_order = seq_len(nrow(pseudobulk_matrix)),
                        column_order = seq_len(ncol(pseudobulk_matrix)),
                        name = legend_name,
                        row_names_side = "left",
                        col = env$color_options$continuous[[input$colour_scheme]],
                        cell_fun = function(j, i, x, y, width, height, fill) {
                            if (input$text_size > 0.1) {
                                grid::grid.text(
                                    sprintf("%.2f", pseudobulk_matrix[i, j]),
                                    x, y, just = "center",
                                    gp = grid::gpar(fontsize = input$text_size)
                                )
                            }
                        },
                        column_names_gp = grid::gpar(fontsize = input$axis_size),
                        row_names_gp = grid::gpar(fontsize = input$axis_size),
                        heatmap_legend_param = list(
                            title_gp = grid::gpar(fontsize = input$legend_size),
                            labels_gp = grid::gpar(fontsize = input$legend_size)
                        )
                    )
                })
            })

            dotplot <- shiny::reactive({
                ptype <- input$plot_type
                pseudobulk_matrix <- psdbulk_matrix()
                perc <- percentages_mat()
                module_ord <- module_ordering()
                req_gear_heatmap(session)

                shiny::isolate({
                    shiny::req(ptype == "dotplot")
                    shiny::req(!is.null(pseudobulk_matrix), module_ord)
                    shiny::req(!is.null(perc))
                    shiny::req(all(paste("module", module_ord) %in% colnames(pseudobulk_matrix)))
                    shiny::req(all(paste("module", module_ord) %in% colnames(perc)))
                    module_ord <- paste0("module ", module_ord)
                    pseudobulk_matrix <- pseudobulk_matrix[, module_ord, drop = FALSE]
                    gc()

                    if (input$scale_values && ncol(pseudobulk_matrix) > 1) {
                        pseudobulk_matrix <- t(scale(t(pseudobulk_matrix)))
                        pseudobulk_matrix[pseudobulk_matrix > input$cap_value] <- input$cap_value
                        pseudobulk_matrix[pseudobulk_matrix < -input$cap_value] <- -input$cap_value
                    } else {
                        pseudobulk_matrix[pseudobulk_matrix > input$cap_value] <- input$cap_value
                    }

                    perc <- perc[rownames(pseudobulk_matrix), colnames(pseudobulk_matrix), drop = FALSE]

                    colnames(pseudobulk_matrix) <- gsub("module ", "", colnames(pseudobulk_matrix))

                    df <- reshape2::melt(pseudobulk_matrix)
                    colnames(df) <- c("metadata", "module", "expression")
                    df$perc <- reshape2::melt(perc)$value
                    df$metadata <- factor(df$metadata, levels = rev(rownames(pseudobulk_matrix)))

                    df$module <- factor(df$module, levels = colnames(pseudobulk_matrix))

                    ggplot2::ggplot(df, ggplot2::aes(x = .data$module, y = .data$metadata, size = .data$perc, colour = .data$expression)) +
                        ggplot2::geom_point() +
                        ggplot2::scale_colour_gradientn(
                            colours = env$color_options$continuous[[input$colour_scheme]]
                        ) +
                        ggplot2::scale_size(range = input$point_size) +
                        ggplot2::theme_classic() +
                        ggplot2::theme(
                            axis.text = ggplot2::element_text(size = input$axis_size),
                            axis.title = ggplot2::element_text(size = input$axis_size),
                            legend.text = ggplot2::element_text(size = input$legend_size),
                            legend.title = ggplot2::element_text(size = input$legend_size)
                        )
                })

            })

            shiny::observe({
                plot_type <- input$plot_type
                if (plot_type == "heatmap") {
                    plot_htmp <- htmp_obj()
                } else {
                    plot_htmp <- dotplot()
                }
                wdim <- env$window_dim()

                shiny::isolate({
                    shiny::req(!is.null(plot_htmp), wdim, cancelOutput = TRUE)
                    psd_mat <- psdbulk_matrix()
                    plt_height <- min(
                        wdim[2] * 0.9,
                        200 + 35 * nrow(psd_mat)
                    )
                    plt_width <- min(
                        wdim[1] * 0.9,
                        300 + 70 * ncol(psd_mat)
                    )

                    output$pseudobulk_heatmap <- shiny::renderPlot(
                        height = plt_height,
                        width = plt_width,
                        # ifelse(plot_type == "heatmap", ComplexHeatmap::draw(plot_htmp), plot_htmp)
                        plot_htmp
                    )
                })
            })

            shiny::observe({
                plot_type <- input$plot_type
                if (plot_type == "heatmap") {
                    plot_htmp <- htmp_obj()
                } else {
                    plot_htmp <- dotplot()
                }
                shiny::req(plot_htmp)
                req_gear_download(session, "heatmap")

                output$download_heatmap <- shiny::downloadHandler(
                    filename = function() {
                        paste(input$filename_heatmap, tolower(input$filetype_heatmap), sep = ".")
                    },
                    content = function(file) {
                        save_filetypes[[input$filetype_heatmap]](file, width = input$width_heatmap, height = input$height_heatmap)
                        if (plot_type == "heatmap") {
                            ComplexHeatmap::draw(plot_htmp)
                        } else {
                            print(plot_htmp)
                        }
                        grDevices::dev.off()
                    }
                )
            })

        }
    )
}

server_by_cell_heatmap <- function(id, metadata_name, module_ordering) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::updateSelectInput(
                session,
                inputId = "colour_scheme",
                choices = names(env$color_options$continuous),
                selected = names(env$color_options$continuous)[length(env$color_options$continuous)]
            )

            htmp_obj <- shiny::reactive({
                module_order <- module_ordering()
                mtd_name <- metadata_name()
                req_gear_heatmap(session, FALSE)

                shiny::isolate({
                    shiny::req(!is.null(module_order), mtd_name)
                    shinyjs::disable("generate_heatmap")
                    sub_module_list <- env$chosen_modules()[module_order]
                    gene_matrix <- do.call(rbind,
                        lapply(
                            sub_module_list,
                            function(module) {
                                gene_matrix <- read_gene_from_dense_h5(
                                    gene_names = module,
                                    matrix_h5_path = file.path("objects", "expression.h5"),
                                    index_genes = env$genes[module],
                                    check_intersect = FALSE
                                )

                                if (length(env$genes[module]) == 1) {
                                    gene_matrix <- matrix(gene_matrix, nrow = 1)
                                    rownames(gene_matrix) <- module
                                }
                                gene_matrix <- gene_matrix[sort_genes_by_metadata(
                                    expression_matrix = gene_matrix,
                                    metadata_info = env$mtd_df$pseudotime,
                                    summary_function = function(x) { mean(x, na.rm = TRUE)},
                                    thresh_percentile = 0,
                                    thresh_value = 0
                                ), , drop = FALSE]

                                return(gene_matrix)
                            }
                        )
                    )
                    colnames(gene_matrix) <- rownames(env$mtd_df)
                    print("calling generate cell heatmap")

                    htmp_obj <- generate_cell_heatmap(
                        expression_matrix = gene_matrix,
                        gene_family_list = sub_module_list,
                        metadata_df = env$mtd_df,
                        metadata_name = mtd_name,
                        apply_scale = input$scale_values,
                        axis_text_size = input$axis_size,
                        legend_name = "\nexpression",
                        k_smooth = input$k_smooth,
                        cap = input$cap_value,
                        discrete_colour_list = env$color_options$discrete,
                        continuous_colors = env$color_options$continuous[[input$colour_scheme]]
                        # heatmap_width = grid::unit(5, "in"),
                        # heatmap_height = grid::unit(5, "in")
                    )
                    gc()
                    return(htmp_obj)
                })
            }) %>% shiny::bindEvent(input$generate_heatmap)

            shiny::observe({
                htmp_obj_curr <- htmp_obj()

                shiny::isolate({
                    shiny::req(!is.null(htmp_obj_curr), cancelOutput = TRUE)
                    output$by_cell_heatmap <- shiny::renderPlot(
                        {
                            ComplexHeatmap::draw(
                                htmp_obj_curr,
                                merge_legend = TRUE
                            )
                            shinyjs::enable("generate_heatmap")
                            gc()
                        }
                    )
                })
            })

            shiny::observe({
                h_obj <- htmp_obj()
                shiny::req(!is.null(h_obj))
                req_gear_download(session, "heatmap")

                output$download_heatmap <- shiny::downloadHandler(
                    filename = function() {
                        paste(input$filename_heatmap, tolower(input$filetype_heatmap), sep = ".")
                    },
                    content = function(file) {
                        shiny::req(!is.null(h_obj), cancelOutput = TRUE)
                        save_filetypes[[input$filetype_heatmap]](file, width = input$width_heatmap, height = input$height_heatmap)
                        ComplexHeatmap::draw(h_obj, merge_legend = TRUE)
                        grDevices::dev.off()
                        gc()
                    }
                )
            })

        }
    )
}

#' Server - Gene Module Heatmap
#'
#' @description Creates the backend interface for the Gene Module Heatmap panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
server_module_metadata_heatmap <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::updateSelectInput(
                session,
                inputId = "metadata",
                choices = names(env$discrete_mtd),
                selected = names(env$discrete_mtd)[1]
            )

            shiny::observe({
                psd_order <- env$psd_ordering()

                shiny::isolate({
                    shiny::req(length(psd_order) > 0)
                    group_stats_df <- lapply(psd_order, function(x) {
                                c(x$coverage, x$iqr_stats)
                            })
                    names(group_stats_df) <- names(psd_order)
                    group_stats_df <- as.data.frame(do.call(rbind, group_stats_df))
                    colnames(group_stats_df) <- c("Coverage", "min IQR", "Q1 IQR", "median IQR", "Q3 IQR", "max IQR")
                    q_coverage <- max(group_stats_df$Coverage) * 0.85
                    group_stats_df <- rbind(
                        group_stats_df %>% dplyr::filter(.data$Coverage > q_coverage) %>% dplyr::arrange(.data$"median IQR"),
                        group_stats_df %>% dplyr::filter(.data$Coverage <= q_coverage) %>% dplyr::arrange(dplyr::desc(.data$Coverage))
                    )
                    best_option <- rownames(group_stats_df)[1]
                    shiny::updateSelectInput(
                        session,
                        inputId = "metadata",
                        choices = rownames(group_stats_df),
                        selected = best_option
                    )
                })
            }) %>% shiny::bindEvent(env$psd_ordering())

            shiny::observe({
                mod_summ <- env$modules_stats_summary()
                shiny::isolate({
                    shiny::req(mod_summ, cancelOutput = TRUE)

                    shiny::updateSelectizeInput(
                        session,
                        inputId = "module_order",
                        choices = rownames(mod_summ),
                        selected = rownames(mod_summ)[mod_summ$is_outlier == "no"],
                        server = TRUE
                    )
                })
            })

            output$clip <- shiny::renderUI({
                rclipboard::rclipButton(
                    inputId = "clipbtn",
                    label = "Copy to clipboard",
                    clipText = paste(input$module_order, collapse = ","),
                    icon = shiny::icon("clipboard"),
                    class = "btn btn-primary",
                    options = list(delay = list(show = 800, hide = 100), trigger = "hover")
                )
            })

            shiny::observe({
                modules_stats <- env$modules_stats_summary()
                shiny::updateSelectizeInput(
                    session,
                    inputId = "module_order",
                    selected = modules_stats$module[modules_stats$is_outlier == "no"]
                )
            }) %>% shiny::bindEvent(input$reset_modules)

            shiny::observe({
                modules_stats <- env$modules_stats_summary()
                slct_modules <- input$module_order
                shiny::req(slct_modules)
                shiny::updateSelectizeInput(
                    session,
                    inputId = "module_order",
                    selected = intersect(modules_stats$module, slct_modules)
                )

            }) %>% shiny::bindEvent(input$order_modules)

            module_ordering <- shiny::reactive({
                input$module_order
            }) %>% shiny::debounce(3000)

            shiny::observe({
                mod_summ <- env$modules_stats_summary()
                shiny::updateSelectizeInput(
                    session,
                    inputId = "module_order",
                    choices = rownames(mod_summ),
                    selected = rownames(mod_summ)[mod_summ$is_outlier == "no"]
                )
            }) %>% shiny::bindEvent(input$reset_order)

            metadata_name <- shiny::reactive(input$metadata)

            server_pseudotime_table("psd_table")
            server_pseudotime_mtd_order("mtd_ordering", metadata_name)
            server_pseudobulk_heatmap("pseudobulk", metadata_name, module_ordering)
            server_by_cell_heatmap("by_cell", metadata_name, module_ordering)

            shiny::observe({
                print(paste("Size of umap_df", format(utils::object.size(env$umap_df), units = "Mb")))
                print(paste("Size of mtd_df", format(utils::object.size(env$mtd_df), units = "Mb")))
                print(paste("Size of discrete_mtd", format(utils::object.size(env$discrete_mtd), units = "Mb")))
                print(paste("Size of mon_obj", format(utils::object.size(env$mon_obj), units = "Mb")))
                print(paste("Size of trajectory_gplot", format(utils::object.size(env$trajectory_gplot), units = "Mb")))
                print(paste("Size of moran_df", format(utils::object.size(env$moran_df), units = "Mb")))
                print(paste("Size of modules_summaries", format(utils::object.size(env$modules_summaries()), units = "Mb")))
                print(paste("Size of modules_mask", format(utils::object.size(env$modules_mask()), units = "Mb")))
                print(paste("Size of modules_stats", format(utils::object.size(env$modules_stats()), units = "Mb")))
                print(paste("Size of psd_value", format(utils::object.size(env$psd_value()), units = "Mb")))
            }) %>% shiny::bindEvent(input$size_env)
        }
    )
}
