#' @importFrom dplyr %>% .data
NULL

###### UI ######
ui_pseudobulk_heatmap <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h3("Pseudobulk heatmap"),
        shiny::splitLayout(
            cellWidths = "40px",
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
                        value = FALSE,
                        fill = TRUE
                    ),
                    shiny::sliderInput(
                        inputId = ns("text_size"),
                        label = "Text size",
                        min = 0.00, max = 30.00, value = 10, step = 0.05
                    ),
                    shiny::sliderInput(
                        inputId = ns("axis_size"),
                        label = "Axis labels size",
                        min = 0.00, max = 30.00, value = 10, step = 0.05
                    ),
                    shiny::sliderInput(
                        inputId = ns("legend_size"),
                        label = "Legend size",
                        min = 0.00, max = 30.00, value = 10, step = 0.05
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
                                value = FALSE,
                                fill = TRUE
                            ),
                            shiny::numericInput(
                                inputId = ns("k_smooth"),
                                label = "Smoothing degree",
                                min = 0, max = 100, value = 0, step = 1
                            ),
                            shiny::numericInput(
                                inputId = ns("cap_value"),
                                label = "Value cap",
                                min = 0, max = 20, value = 10, step = 0.1
                            ),
                            shiny::sliderInput(
                                inputId = ns("axis_size"),
                                label = "Axis labels size",
                                min = 0.00, max = 30.00, value = 10, step = 0.05
                            ),
                            shiny::sliderInput(
                                inputId = ns("legend_size"),
                                label = "Legend size",
                                min = 0.00, max = 30.00, value = 10, step = 0.05
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
            shiny::column(2,
                shiny::selectInput(
                    inputId = ns("module_order"),
                    label = "Predefined order of modules:",
                    choices = c(""),
                    selected = NULL,
                    multiple = TRUE
                )
            ),
            shiny::column(1,
                shiny::actionButton(
                    inputId = ns("reset_order"),
                    label = "Reset order!",
                    icon = shiny::icon("undo")
                )
            ),
            shiny::column(2,
                shiny::actionButton(
                    inputId = ns("generate_heatmap"),
                    label = "Generate heatmap",
                    class = "btn-danger"
                )
            )
        ),
        shiny::plotOutput(ns("by_cell_heatmap"), height = "auto")
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
        shiny::h2("Gene modules heatmaps", id = "heatmaps_panel"),
        shiny::selectInput(
            inputId = ns("metadata"),
            label = "Select metadata",
            choices = c(""),
            selected = NULL
        ),
        ui_pseudobulk_heatmap(ns("pseudobulk")),
        ui_by_cell_heatmap(ns("by_cell"))
    )
}

###### SERVER ######
server_pseudobulk_heatmap <- function(id, metadata_name) {
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
                shiny::req(current_modules, mtd_name, summarise_expr)
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
                                gene_name = module,
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

                return(pseudobulk_matrix)
            })

            htmp_obj <- shiny::reactive({
                pseudobulk_matrix <- psdbulk_matrix()
                req_gear_heatmap(session)

                shiny::isolate({
                    shiny::req(!is.null(pseudobulk_matrix))

                    if (input$scale_values) {
                        pseudobulk_matrix <- t(scale(t(pseudobulk_matrix)))
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

            shiny::observe({
                shiny::req(htmp_obj(), env$window_dim())

                shiny::isolate({
                    plt_height <- min(
                        env$window_dim()[2],
                        200 + 35 * nrow(htmp_obj())
                    )
                    plt_width <- min(
                        env$window_dim()[1],
                        300 + 70 * ncol(htmp_obj())
                    )

                    output$pseudobulk_heatmap <- shiny::renderPlot(
                        height = plt_height,
                        width = plt_width,
                        ComplexHeatmap::draw(htmp_obj())
                    )
                })
            })

            shiny::observe({
                shiny::req(htmp_obj())
                req_gear_download(session, "heatmap")

                output$download_heatmap <- shiny::downloadHandler(
                    filename = function() {
                        paste(input$filename_heatmap, tolower(input$filetype_heatmap), sep = ".")
                    },
                    content = function(file) {
                        save_filetypes[[input$filetype_heatmap]](file, width = input$width_heatmap, height = input$height_heatmap)
                        ComplexHeatmap::draw(htmp_obj())
                        dev.off()
                    }
                )
            })

        }
    )
}

server_by_cell_heatmap <- function(id, metadata_name) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shiny::updateSelectInput(
                session,
                inputId = "colour_scheme",
                choices = names(env$color_options$continuous),
                selected = names(env$color_options$continuous)[length(env$color_options$continuous)]
            )

            shiny::observe({
                current_modules <- env$chosen_modules()
                shiny::req(current_modules)
                shiny::isolate({
                    shiny::updateSelectInput(
                        session,
                        inputId = "module_order",
                        choices = names(current_modules),
                        selected = names(current_modules)[1]
                    )
                })
            })
            htmp_obj <- shiny::reactive({
                shiny::req(!is.null(input$module_order), metadata_name())
                req_gear_heatmap(session, FALSE)

                # for (module_name in names(env$chosen_modules())) {
                #     # print(paste(Sys.time(), module_name))
                #     module_genes <- env$chosen_modules()[[module_name]]
                #     gene_matrix <- read_gene_from_dense_h5(
                #         gene_name = module_genes,
                #         matrix_h5_path = file.path("objects", "expression.h5"),
                #         index_genes = env$genes[module_genes],
                #         check_intersect = FALSE
                #     )

                #     colnames(gene_matrix) <- rownames(env$mtd_df)

                #     selected_cells <- select_cells_by_gene_expr(
                #         expression_matrix = gene_matrix,
                #         genes = module_genes,
                #         thresh_percentile = 0,
                #         thresh_value = 0,
                #         n_coexpressed_thresh = 0
                #     )

                #     # print(length(selected_cells))
                #     # print(fivenum(env$mtd_df[selected_cells, "pseudotime"]))


                # }
                shiny::isolate({
                    print(paste(Sys.time(), "rendering by cell"))
                    module_order <- input$module_order
                    sub_module_list <- env$chosen_modules()[module_order]
                    gene_matrix <- do.call(rbind,
                        lapply(
                            sub_module_list,
                            function(module) {
                                gene_matrix <- read_gene_from_dense_h5(
                                    gene_name = module,
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

                    htmp_obj <- generate_cell_heatmap(
                        expression_matrix = gene_matrix,
                        gene_family_list = sub_module_list,
                        metadata_df = env$mtd_df,
                        metadata_name = metadata_name(),
                        apply_scale = input$scale_values,
                        axis_text_size = input$axis_size,
                        legend_name = "\nexpression",
                        k_smooth = input$k_smooth,
                        cap = input$cap_value,
                        discrete_colour_list = env$color_options$discrete,
                        continuous_colors = env$color_options$continuous[[input$colour_scheme]]
                    )
                    shinyjs::enable("generate_heatmap")
                    return(htmp_obj)
                })
            }) %>% shiny::bindEvent(input$generate_heatmap)

            shiny::observe({
                shiny::req(!is.null(htmp_obj()), env$window_dim(), cancelOutput = TRUE)

                shiny::isolate({
                    plt_height <- min(
                        env$window_dim()[2] / 1.1,
                        200 + 25 * nrow(htmp_obj())
                    )
                    plt_width <- min(
                        env$window_dim()[1] / 1.1,
                        300 + 70 * ncol(htmp_obj())
                    )
                    output$by_cell_heatmap <- shiny::renderPlot(
                        height = plt_height,
                        width = plt_width,
                        {
                            shiny::req(!is.null(htmp_obj()), cancelOutput = TRUE)
                            ComplexHeatmap::draw(
                                htmp_obj(),
                                merge_legend = TRUE
                            )
                        }
                    )
                })
            })

            shiny::observe({
                shiny::req(!is.null(htmp_obj()))
                req_gear_download(session, "heatmap")

                output$download_heatmap <- shiny::downloadHandler(
                    filename = function() {
                        paste(input$filename_heatmap, tolower(input$filetype_heatmap), sep = ".")
                    },
                    content = function(file) {
                        save_filetypes[[input$filetype_heatmap]](file, width = input$width_heatmap, height = input$height_heatmap)
                        ComplexHeatmap::draw(htmp_obj(), merge_legend = TRUE)
                        dev.off()
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

            metadata_name <- shiny::reactive(input$metadata)

            server_pseudobulk_heatmap("pseudobulk", metadata_name)
            server_by_cell_heatmap("by_cell", metadata_name)
        }
    )
}
