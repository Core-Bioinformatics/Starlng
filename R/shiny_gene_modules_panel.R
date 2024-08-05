#' @importFrom dplyr %>% .data
NULL

###### UI ######
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
            cellWidths = "40px",
            gear_umaps(ns, "settings", "highest", TRUE, TRUE),
            gear_download(ns, "module_umap", "modules")
        ),
        shiny::fluidRow(
            shiny::column(1,
                shiny::numericInput(
                    inputId = ns("n_columns"),
                    label = "Number of columns",
                    min = 1, max = 20, value = 0, step = 1
                )
            ),
            shiny::column(4,
                shinyWidgets::radioGroupButtons(
                    inputId = ns("summarise_expr"),
                    label = "How to summarise?",
                    choices = c("average", "average expressed", "#genes")
                )
            )
        ),
        shiny::plotOutput(ns("module_umap_plot"), height = "auto")
    )
}

###### SERVER ######
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
            default_col_scheme <- "white_red"
            if (!(default_col_scheme %in% names(env$color_options$continuous))) {
                default_col_scheme <- names(env$color_options$continuous)[1]
            }

            shiny::observe({
                shiny::req(env$chosen_modules())
                shiny::updateNumericInput(
                    session,
                    inputId = "n_columns",
                    value = min(6, length(env$chosen_modules()))
                )
            })

            shiny::updateSelectInput(
                session,
                inputId = "settings_colour_scheme",
                choices = names(env$color_options$continuous),
                selected = default_col_scheme
            )

            composite_plot <- shiny::reactive({
                shiny::req(env$window_dim(), env$chosen_modules(), input$n_columns > 0, input$summarise_expr)
                req_gear_umap(session, "settings")

                shiny::isolate({
                    colourbar_height <- min(env$window_dim()) / 2

                    plot_umap_gene_modules_shiny(
                        shiny_env = env,
                        gene_modules = env$chosen_modules(),
                        summarise_expr = input$summarise_expr,
                        n_columns = input$n_columns,
                        scale_values = input$settings_scale,
                        trajectory_width = input$settings_trajectory_width,
                        cell_sort_order = input$settings_pt_order,
                        cell_size = input$settings_pt_size,
                        cell_alpha = input$settings_pt_alpha,
                        legend_text_size = input$settings_legend_size,
                        axis_text_size = input$settings_axis_size,
                        colourbar_height = colourbar_height,
                        continuous_colors = env$color_options$continuous[[input$settings_colour_scheme]]
                    )
                })
            })

            shiny::observe({
                shiny::req(composite_plot(), env$window_dim())

                shiny::isolate({
                    plt_height <- env$window_dim()[2] / 1.05
                    plt_width <- env$window_dim()[1] / 1.05

                    output$module_umap_plot <- shiny::renderPlot(
                        height = plt_height,
                        width = plt_width,
                        composite_plot()
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
