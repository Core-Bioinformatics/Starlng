#' @importFrom dplyr %>% .data
NULL

###### UI ######


ui_cell_select_filter <- function(id) {
    ns <- shiny::NS(id)

    shiny::tabsetPanel(
        id = ns("tabset"),
        shiny::tabPanel(
            title = "By metadata",
            value = "metadata",
            shiny::wellPanel(
                shiny::splitLayout(
                    shiny::selectInput(
                        inputId = ns("metadata"),
                        choices = NULL,
                        label = "Select metadata"
                    ),
                    shiny::verticalLayout(
                        shiny::tags$b("Select groups"),
                        shinyWidgets::pickerInput(
                            inputId = ns("metadata_groups"),
                            choices = NULL,
                            options = list(
                                `actions-box` = TRUE,
                                title = "Select/deselect groups",
                                size = 10,
                                width = "90%",
                                `selected-text-format` = "count > 3"
                            ),
                            multiple = TRUE
                        )
                    )
                )
            )
        ),
        shiny::tabPanel(
            title = "By gene expression",
            value = "gene",
            shiny::wellPanel(
                shiny::selectizeInput(
                    inputId = ns("gene"),
                    label = "Select gene",
                    choices = c(""),
                    selected = NULL,
                    width = "100%",
                    multiple = TRUE
                ),
                shiny::splitLayout(
                    shiny::numericInput(
                        inputId = ns("expr_threshold"),
                        label = "Expression threshold",
                        min = 0, max = 10, value = 0, step = 0.01,
                        width = "95%"
                    ),
                    shiny::sliderInput(
                        inputId = ns("expr_threshold_perc"),
                        label = "Expression %threshold",
                        min = 0, max = 100, value = 50, step = 0.5,
                        width = "45%"
                    )
                ),
                shiny::splitLayout(
                    shiny::numericInput(
                        inputId = ns("relaxation"),
                        label = "#genes not expressed",
                        min = 0, max = 10, value = 0, step = 1,
                        width = "95%"
                    ),
                    shinyWidgets::pickerInput(
                        inputId = ns("gene_modules"),
                        label = "Filter by modules",
                        choices = c(),
                        multiple = TRUE,
                        options = list(
                            `actions-box` = TRUE,
                            title = "Select/deselect groups",
                            size = 10,
                            width = "40%",
                            `selected-text-format` = "count > 3"
                        )
                    )
                )
            )
        )
    )
}

ui_pseudotime_select_cells_panel <- function(id) {
    ns <- shiny::NS(id)
    actual_id <- strsplit(id, "-")[[1]]
    actual_id <- actual_id[length(actual_id)]

    shiny::tagList(
        shiny::h3(paste0("Select the ", actual_id, " cells")),
        shiny::sliderInput(
            inputId = ns("n_filters"),
            label = "Number of filters",
            min = 1, max = 5, value = 1, step = 1
        ),
        shiny::uiOutput(ns("filters")),
        # ui_cell_select_filter(ns("filter")),
        shiny::splitLayout(
            gear_umaps(ns, "settings", "highest", FALSE),
            shiny::sliderInput(
                inputId = ns("dist_threshold"),
                label = "Percentage outlier distance",
                min = 0, max = 100, value = ifelse(actual_id == "start", 75, 100), step = 0.5
            ),
        ),
        shiny::plotOutput(ns("cells_highlight"), height = "auto")
    )
}

#' UI - Pseudotime
#'
#' @description Creates the UI interface for the Pseudotime panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
ui_pseudotime_select_cells <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Pseudotime inference"),
        shiny::splitLayout(
            cellWidths = c("50%", "50%"),
            ui_pseudotime_select_cells_panel(ns("start")),
            ui_pseudotime_select_cells_panel(ns("end"))
        ),
        shiny::fluidRow(
            shiny::column(width = 7, offset = 4,
                shiny::actionButton(
                    ns("psd_button"),
                    "Calculate pseudotime!",
                    class = "btn-danger"
                ),
                shiny::actionButton(
                    ns("psd_reset"),
                    "Reset pseudotime",
                    class = "btn-default"
                ),
                gear_umaps(ns, "settings", "highest"),
                shiny::plotOutput(ns("pseudotime_plot"), height = "auto")
            )
        )
    )
}

###### SERVER ######
server_cell_select_filter <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            cells_list <- list("gene" = c(), "metadata" = c())
            
            cells_list[["metadata"]] <- shiny::reactive({
                shiny::req(input$metadata %in% names(env$discrete_mtd))
                shiny::req(!is.null(input$metadata_groups))
                shiny::req(all(input$metadata_groups %in% env$discrete_mtd[[input$metadata]]))

                ret_list <- list(input$metadata_groups)
                names(ret_list) <- input$metadata
                return(select_cells_by_metadata(
                    env$mtd_df,
                    ret_list
                ))
            })

            cells_list[["gene"]] <- shiny::reactive({
                shiny::req(input$gene, input$expr_threshold, input$expr_threshold_perc, input$relaxation)
                input$gene_modules

                shiny::isolate({
                    if (!is.null(input$gene_modules) && length(input$gene_modules) > 0 && length(input$gene_modules) < length(env$chosen_modules())) {
                        gene_list <- do.call(c, env$chosen_modules()[input$gene_modules])
                        n_coexpressed_thresh <- 0
                    } else {
                        shiny::req(length(input$gene) > 0)
                        gene_list <- input$gene
                        n_coexpressed_thresh <- length(gene_list) - input$relaxation
                    }

                    return(select_cells_by_gene_expr(
                        expression_matrix = read_gene_from_dense_h5(
                            gene_names = gene_list,
                            matrix_h5_path = file.path("objects", "expression.h5"),
                            index_genes = env$genes[gene_list],
                            check_intersect = FALSE,
                            add_rownames = TRUE
                        ),
                        genes = gene_list,
                        thresh_percentile = input$expr_threshold_perc / 100,
                        thresh_value = input$expr_threshold,
                        n_coexpressed_thresh = n_coexpressed_thresh
                    ))
                })
            })

            shiny::reactive({
                shiny::req(input$tabset, cells_list[[input$tabset]]())

                shiny::isolate({
                    return(cells_list[[input$tabset]]())
                })
            })
        }
    )
}

server_pseudotime_select_cells_panel <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            highlight_colour <- ifelse(id == "start", "red", "blue")
            shiny::observe({
                shiny::req(env$pseudotime_changes() > 0, id == "start")

                shiny::updateSliderInput(
                    session,
                    inputId = "dist_threshold",
                    value = 1
                )
            })

            shiny::observe({
                ns <- session$ns
                nfilters <- input$n_filters

                shiny::isolate({
                    output$filters <- shiny::renderUI({
                        spsComps::onNextInput({
                            for (i in seq_len(nfilters)) {
                                recommended_option <- names(env$discrete_mtd)[1]
                                if (i == 1 && id == "start") {
                                    recommended_option <- env$recommended_psd$recommended_mtd_name
                                    env$pseudotime_changes(1)
                                }
                                shiny::updateSelectInput(
                                    session,
                                    inputId = paste0("filter_", i, "-metadata"),
                                    choices = names(env$discrete_mtd),
                                    selected = recommended_option
                                )

                                shiny::updateSelectizeInput(
                                    session,
                                    inputId = paste0("filter_", i, "-gene"),
                                    choices = names(env$genes),
                                    selected = names(env$genes)[1],
                                    server = TRUE,
                                    options = list(
                                        placeholder = "Select gene(s)",
                                        maxOptions = 7,
                                        create = TRUE,
                                        persist = TRUE
                                    )
                                )
                            }
                        })

                        do.call(
                            shiny::tagList,
                            lapply(seq_len(nfilters), function(i) {
                                ui_cell_select_filter(ns(paste0("filter_", i)))
                            })
                        )
                    })
                })

            }) %>% shiny::bindEvent(input$n_filters)

            shiny::observe({
                for (i in seq_len(input$n_filters)) {
                    mtd <- input[[paste0("filter_", i, "-metadata")]]
                    shiny::req(!is.null(mtd))
                    shiny::req(mtd %in% names(env$discrete_mtd))
                    shinyWidgets::updatePickerInput(
                        session,
                        inputId = paste0("filter_", i, "-metadata_groups"),
                        choices = env$discrete_mtd[[mtd]],
                        selected = env$discrete_mtd[[mtd]]
                    )
                }

                shiny::isolate({
                    psd_change <- env$pseudotime_changes()
                    if (id == "start" && psd_change > 0) {
                        env$pseudotime_changes(0)
                        shinyWidgets::updatePickerInput(
                            session,
                            inputId = "filter_1-metadata_groups",
                            selected = env$recommended_psd$recommended_mtd_group
                        )
                    }
                })
            })

            shiny::observe({
                chosen_modules <- env$chosen_modules()
                shiny::req(chosen_modules)
                for (i in seq_len(input$n_filters)) {
                    shinyWidgets::updatePickerInput(
                        session,
                        inputId = paste0("filter_", i, "-gene_modules"),
                        choices = names(chosen_modules),
                        selected = names(chosen_modules)
                    )
                }
            })

            per_filter_cells <- shiny::reactive({
                lapply(
                    seq_len(input$n_filters), function(i) {
                        server_cell_select_filter(paste0("filter_", i))
                    }
                )
            })

            filtered_cells <- shiny::reactive({
                shiny::req(per_filter_cells())
                nfilters <- length(per_filter_cells())
                shiny::req(nfilters > 0)

                for (i in seq_len(nfilters)) {
                    shiny::req(per_filter_cells()[[i]]())
                }

                is_first <- TRUE
                for (i in seq_len(nfilters)) {
                    current_cells <- per_filter_cells()[[i]]()
                    if (is_first) {
                        final_cells <- current_cells
                        is_first <- FALSE
                    } else {
                        final_cells <- intersect(final_cells, current_cells)
                    }
                }

                return(final_cells)
            })


            actual_cells <- shiny::reactive({
                shiny::req(filtered_cells(), input$dist_threshold)

                remove_outlier_cells(
                    filtered_cells(),
                    env$umap_df,
                    input$dist_threshold / 100
                )
            })

            shiny::observe({
                shiny::req(actual_cells(), env$window_dim())

                # shiny::isolate({
                    plt_height <- min(
                        floor(env$height_ratio * env$window_dim()[2]),
                        env$window_dim()[1] / 2.1
                    )

                    output$cells_highlight <- shiny::renderPlot(
                        height = plt_height,
                        width = plt_height,
                        {
                            mask <- rep("not selected", nrow(env$mtd_df))
                            names(mask) <- rownames(env$mtd_df)
                            mask[actual_cells()] <- "selected"
                            mask <- factor(mask)
                            gplot_obj <- plot_umap(
                                umap_embedding = env$umap_df,
                                cell_info = mask,
                                mtd_name = "Selected cells",
                                cell_sort_order = c("not selected", "selected"),
                                cell_size = input$settings_pt_size,
                                cell_alpha = input$settings_pt_alpha,
                                legend_text_size = input$settings_legend_size,
                                axis_text_size = input$settings_axis_size,
                                discrete_colors = list(
                                    "1" = highlight_colour,
                                    "2" = c("gray", highlight_colour)
                                )
                            )

                            traj_layer <- env$trajectory_gplot$layers[[1]]
                            traj_layer$aes_params$size <- input$settings_trajectory_width 
                            gplot_obj$layers <- c(gplot_obj$layers, traj_layer)
                            gplot_obj
                        }

                    )
                # })

            })
            
            return(actual_cells)
        }
    )
}

#' Server - Pseudotime
#'
#' @description Creates the backend interface for the Pseudotime panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
server_pseudotime_select_cells <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::updateSelectInput(
                session,
                inputId = "settings_colour_scheme",
                choices = names(env$color_options$continuous),
                selected = names(env$color_options$continuous)[1]
            )


            start_cells <- server_pseudotime_select_cells_panel("start")
            end_cells <- server_pseudotime_select_cells_panel("end")

            shiny::observe({
                shiny::req(start_cells(), end_cells())

                shiny::isolate({
                    shinyjs::disable("psd_button")
                    current_start <- start_cells()
                    current_end <- end_cells()

                    if (length(current_end) %in% c(0, nrow(env$mtd_df))) {
                        current_end <- NULL
                    } else{
                        if (length(current_start) %in% c(0, nrow(env$mtd_df))) {
                            current_start <- NULL
                        } else {
                            current_end <- setdiff(current_end, current_start)
                        }

                        if (length(current_end) == 0) {
                            current_end <- NULL
                        }
                    }

                    psd_value <- custom_pseudotime_ordering(
                        env$mon_obj,
                        start_cells = current_start,
                        end_cells = current_end
                    )

                    temp_mtd_df <- env$mtd_df
                    temp_mtd_df$pseudotime <- psd_value
                    assign("mtd_df", temp_mtd_df, envir = env)

                    env$psd_value(psd_value)
                })
            }) %>% shiny::bindEvent(input$psd_button)

            shiny::observe({
                env$psd_value(env$recommended_psd$recommended_pseudotime)
                temp_mtd_df <- env$mtd_df
                temp_mtd_df$pseudotime <- env$recommended_psd$recommended_pseudotime
                assign("mtd_df", temp_mtd_df, envir = env)

            }) %>% shiny::bindEvent(input$psd_reset)

            shiny::observe({
                pseudotime_value <- env$psd_value()
                wdim <- env$window_dim()

                shiny::isolate({
                    shiny::req(pseudotime_value, wdim)
                    plt_height <- min(
                        floor(env$height_ratio * wdim[2]),
                        wdim[1]
                    )

                    output$pseudotime_plot <- shiny::renderPlot(
                        width = plt_height,
                        height = plt_height,
                        {
                            gplot_obj <- plot_umap(
                                umap_embedding = env$umap_df,
                                cell_info = pseudotime_value,
                                mtd_name = "pseudotime",
                                cell_sort_order = input$settings_pt_order,
                                cell_size = input$settings_pt_size,
                                cell_alpha = input$settings_pt_alpha,
                                legend_text_size = input$settings_legend_size,
                                axis_text_size = input$settings_axis_size,
                                discrete_colors = env$color_options$discrete,
                                colourbar_width = plt_height,
                                continuous_colors = env$color_options$continuous[[input$settings_colour_scheme]]
                            )

                            traj_layer <- env$trajectory_gplot$layers[[1]]
                            traj_layer$aes_params$size <- input$settings_trajectory_width 
                            gplot_obj$layers <- c(gplot_obj$layers, traj_layer)
                            gplot_obj
                        }
                    )

                    shinyjs::enable("psd_button")
                })
            })
        }
    )
}
