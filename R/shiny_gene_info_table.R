#' @importFrom dplyr %>% .data

###### UI ######
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
server_gene_info_table <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            filtered_moran_df <- shiny::reactive({
                shiny::req(input$moran_q, input$moran_i)
                env$moran_df %>%
                    dplyr::filter(.data$q_value < input$moran_q, .data$morans_I > input$moran_i)
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
