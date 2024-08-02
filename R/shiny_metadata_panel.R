#' @importFrom dplyr %>% .data

###### UI ######
ui_metadata_umap_panel <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::splitLayout(
            shiny::selectInput(
                inputId = ns("metadata"),
                label = "Select metadata",
                choices = c(""),
                selected = NULL
            )
        ),
        shiny::splitLayout(
            cellWidths = "40px",
            gear_umaps(ns, "settings", "highest"),
            gear_download(ns, "umap", "metadata")
        ),
        shiny::plotOutput(ns("umap_plot"), height = "auto")
    )
}

ui_metadata_umap <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Metadata visualisation", id = "metadata_panel"),
        shiny::splitLayout(
            cellWidths = c("50%", "50%"),
            ui_metadata_umap_panel(ns("left")),
            ui_metadata_umap_panel(ns("right"))
        )
    )
}

###### SERVER ######
server_metadata_umap_panel <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {

            umap_ggplot <- shiny::reactive({
                shiny::req(input$metadata)
                req_gear_umap(session, "settings")

                colourbar_width <- min(
                    env$window_dim()[1] / 2.2,
                    env$window_dim()[2] * env$height_ratio
                ) / 1.2

                shiny::isolate({
                    gplot_obj <- plot_umap(
                        umap_embedding = env$umap_df,
                        cell_info = env$mtd_df[, input$metadata],
                        mtd_name = input$metadata,
                        cell_sort_order = input$settings_pt_order,
                        scale_values = input$settings_scale,
                        cell_size = input$settings_pt_size,
                        cell_alpha = input$settings_pt_alpha,
                        legend_text_size = input$settings_legend_size,
                        axis_text_size = input$settings_axis_size,
                        discrete_colors = env$color_options$discrete,
                        colourbar_width = colourbar_width,
                        continuous_colors = env$color_options$continuous[[input$settings_colour_scheme]]
                    )

                    traj_layer <- env$trajectory_gplot$layers[[1]]
                    traj_layer$aes_params$size <- input$settings_trajectory_width 
                    gplot_obj$layers <- c(gplot_obj$layers, traj_layer)
                    gplot_obj
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
                            gplot_obj <- umap_ggplot()
                            is_discrete <- inherits(env$mtd_df[[input$metadata]], c("factor", "character"))
                            if (!is_discrete) {
                                gplot_obj <- gplot_obj +
                                    ggplot2::guides(colour = ggplot2::guide_colourbar(
                                        barwidth = grid::unit(input$width_umap / 1.5, "inches")
                                    ))
                            }
                            ggplot2::ggsave(
                                filename = file,
                                plot = gplot_obj,
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

server_metadata_umap <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            for (panel in c("left", "right")) {
                shiny::updateSelectInput(
                    session,
                    inputId = paste0(panel, "-settings_colour_scheme"),
                    choices = names(env$color_options$continuous),
                    selected = names(env$color_options$continuous)[1]
                )
            }
            shiny::observe({
                for (panel in c("left", "right")) {
                    shiny::updateSelectInput(
                        session,
                        inputId = paste0(panel, "-metadata"),
                        choices = colnames(env$mtd_df),
                        selected = colnames(env$mtd_df)[1]
                    )
                }
            }) %>% shiny::bindEvent(env$pseudotime_changes())
            server_metadata_umap_panel("left")
            server_metadata_umap_panel("right")
        }
    )
}
