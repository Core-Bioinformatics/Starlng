#' @importFrom dplyr %>% .data
NULL

###### UI ######
#' UI - Gene Info Table
#'
#' @description Creates the UI interface for the Gene Info Tabel panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
ui_gene_info_table <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Gene information table", id = "gene_info_panel"),
        shinyWidgets::dropdownButton(
            label = "",
            icon = shiny::icon("cog"),
            status = "success",
            size = "sm",
            shiny::sliderInput(
                inputId = ns("moran_q"),
                label = "Upper threshold for q_value",
                min = 0.00,
                max = 1,
                value = 0.05,
                step = 0.001
            ),
            shiny::sliderInput(
                inputId = ns("moran_i"),
                label = "Lower threshold for Moran I",
                min = -1,
                max = 1,
                value = 0.1,
                step = 0.01
            ),
            shiny::numericInput(
                inputId = ns("average"),
                label = "Lower threshold for average expression",
                min = 0,
                max = 20,
                value = 0,
                step = 0.1
            ),
            shiny::numericInput(
                inputId = ns("average_nonzero"),
                label = "Lower threshold for average expression (non-zero)",
                min = 0,
                max = 20,
                value = 0.5,
                step = 0.1
            ),
            shiny::sliderInput(
                inputId = ns("perc_cells"),
                label = "Upper threshold for % cells expressing",
                min = 0,
                max = 100,
                value = 75,
                step = 0.5
            )
        ),
        shiny::p("High values of I indicate that nearby cells have very similar values of the gene's expression."),
        shiny::p("Low q_value are a good indicator of good positive autocorrelation"),
        DT::dataTableOutput(ns("gene_table")),
        shiny::downloadButton(
            outputId = ns("download_genes"),
            label = "Download genes as CSV"
        )
    )
}

###### SERVER######
#' Server - Gene Info Table
#'
#' @description Creates the backend interface for the Gene Info Tabel panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
server_gene_info_table <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            filtered_moran_df <- shiny::reactive({
                shiny::req(input$moran_q, input$moran_i, input$average, input$average_nonzero, input$perc_cells)
                env$moran_df %>%
                    dplyr::filter(.data$q_value < input$moran_q, .data$morans_I > input$moran_i, .data$average_expression > input$average, .data$average_expression_nonzero > input$average_nonzero, .data$percent_expressed_cells < input$perc_cells / 100)
            })

            shiny::observe({
                shiny::req(filtered_moran_df())

                shiny::isolate({
                    output$gene_table <- DT::renderDataTable(filtered_moran_df())
                })
            })

            return(shiny::reactive(rownames(filtered_moran_df())))
        }
    )
}
