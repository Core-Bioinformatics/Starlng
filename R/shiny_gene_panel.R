#' @importFrom dplyr %>% .data

###### UI ######
ui_gene_umap_panel <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::splitLayout(
            shiny::selectizeInput(
                inputId = ns("gene"),
                label = "Gene name(s)",
                choices = NULL,
                multiple = TRUE,
                width = "100%"
            )
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
            shinyWidgets::radioGroupButtons(
                inputId = ns("summarise_expr"),
                label = "How to summarise?",
                choices = c("average", "average expressed", "#genes", "binary")
            )
        ),
        shiny::splitLayout(
            cellWidths = "40px",
            gear_umaps(ns, "settings", "highest"),
            gear_download(ns, "umap", "gene")
        ),
        shiny::plotOutput(ns("umap_plot"), height = "auto")
    )
}

ui_gene_umap <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Visualisation of gene expression", id = "gene_panel"),
        shiny::splitLayout(
            cellWidths = c("50%", "50%"),
            ui_gene_umap_panel(ns("left")),
            ui_gene_umap_panel(ns("right"))
        )
    )
}

###### SERVER ######
server_gene_umap_panel <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            umap_ggplot <- shiny::reactive({
                shiny::req(input$gene, input$expr_threshold_perc, input$expr_threshold, input$relaxation, input$summarise_expr)
                req_gear_umap(session, "settings")

                colourbar_width <- min(
                    env$window_dim()[1] / 2.5,
                    env$window_dim()[2] * env$height_ratio
                )

                shiny::isolate({
                    plot_umap_gene_shiny(
                        shiny_env = env,
                        gene_name = input$gene,
                        thresh_percentile = input$expr_threshold_perc,
                        thresh_value = input$expr_threshold,
                        n_coexpressed_threshold = length(input$gene) - input$relaxation,
                        scale_values = input$settings_scale,
                        summarise_expr = input$summarise_expr,
                        trajectory_width = input$settings_trajectory_width,
                        cell_sort_order = input$settings_pt_order,
                        cell_size = input$settings_pt_size,
                        cell_alpha = input$settings_pt_alpha,
                        legend_text_size = input$settings_legend_size,
                        axis_text_size = input$settings_axis_size,
                        colourbar_width = colourbar_width,
                        continuous_colors = env$color_options$continuous[[input$settings_colour_scheme]]
                    )
                })
            })

            shiny::observe({
                shiny::req(umap_ggplot(), env$window_dim())

                plt_height <- min(
                    floor(env$height_ratio * env$window_dim()[2]),
                    env$window_dim()[1] / 2.1
                )
                output$umap_plot <- shiny::renderPlot(
                    height = plt_height,
                    width = plt_height,
                    umap_ggplot()
                )
            })

            shiny::observe({
                shiny::req(umap_ggplot())
                req_gear_download(session, "umap")

                shiny::isolate({
                    output$download_umap <- shiny::downloadHandler(
                        filename = function() {
                            paste(input$filename_umap, tolower(input$filetype_umap), sep = ".")
                        },
                        content = function(file) {
                            ggplot2::ggsave(
                                filename = file,
                                plot = umap_ggplot(),
                                width = input$width_umap,
                                height = input$height_umap
                            )
                        }
                    )
                })
            })

        }
    )
}

server_gene_umap <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            for (panel in c("left", "right")) {
                shiny::updateSelectizeInput(
                    session,
                    inputId = paste0(panel, "-gene"),
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

                default_col_scheme <- "white_red"
                if (!(default_col_scheme %in% names(env$color_options$continuous))) {
                    default_col_scheme <- names(env$color_options$continuous)[1]
                }

                shiny::updateSelectInput(
                    session,
                    inputId = paste0(panel, "-settings_colour_scheme"),
                    choices = names(env$color_options$continuous),
                    selected = default_col_scheme
                )
            }
            server_gene_umap_panel("left")
            server_gene_umap_panel("right")
        }
    )
}
