#' @importFrom dplyr %>% .data

###### UI ######
ui_module_enrichment <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Module enrichment", id = "enrichment_panel"),
        shiny::fluidRow(
            shiny::column(2, shinyWidgets::pickerInput(
                inputId = ns("module_enrichment"),
                label = "Select gene modules",
                choices = c(""),
                multiple = TRUE,
                options = list(
                    `actions-box` = TRUE,
                    title = "Select gene modules",
                    size = 10,
                    width = "90%",
                    `selected-text-format` = "count > 3"
                ),
            )),
            shiny::column(2, shinyWidgets::pickerInput(
                inputId = ns("sources"),
                label = "Select data sources",
                choices = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP"),
                selected = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF", "MIRNA"),
                multiple = TRUE,
                options = list(
                    title = "Select data sources",
                    `actions-box` = TRUE,
                    size = 10,
                    width = "90%",
                    `selected-text-format` = "count > 3"
                )
            )),
            shiny::column(1, shiny::actionButton(
                inputId = ns("run_enrichment"),
                label = "Perform enrichment analysis!",
                class = "btn-danger"
            ))
        ),
        shiny::textOutput(ns("enrichment_message")),
        plotly::plotlyOutput(ns("enrichment_plot"), height = "auto"),
        DT::dataTableOutput(ns("enrichment_table")),
        shiny::downloadButton(
            outputId = ns("download_enrichment"),
            label = "Download enrichment results as CSV"
        )
    )
}

###### SERVER ######
server_module_enrichment <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            shinyjs::hide("run_enrichment")
            shinyjs::hide("download_enrichment")
            shiny::observe({
                shiny::req(env$chosen_modules())
                shinyWidgets::updatePickerInput(
                    session,
                    inputId = "module_enrichment",
                    choices = names(env$chosen_modules()),
                    selected = names(env$chosen_modules())[1]
                )
            })

            shiny::observe({
                shiny::req(!is.null(input$module_enrichment), length(input$module_enrichment) > 0)
                shinyjs::show("run_enrichment")
            })

            gprof_result <- shiny::reactive({
                shiny::req(input$module_enrichment)
                shinyjs::disable("run_enrichment")
                shinyjs::hide("download_enrichment")

                selected_genes <- do.call(c, env$chosen_modules()[input$module_enrichment])
                bg_genes <- names(env$genes)

                gprf_res <- gprofiler2::gost(
                    query = selected_genes,
                    sources = input$enrichment_sources,
                    organism = env$organism,
                    evcodes = TRUE,
                    domain_scope = "custom",
                    custom_bg = bg_genes
                )

                if (!is.null(gprf_res$result)) {
                    gprf_res$result$parents <- sapply(gprf_res$result$parents, toString)
                    shinyjs::show("download_enrichment")
                    output$enrichment_message <- shiny::renderText("")
                } else {
                    output$enrichment_message <- shiny::renderText("No enrichment term found for this module(s).")
                }

                shinyjs::enable("run_enrichment")

                return(gprf_res)
            }) %>% shiny::bindEvent(input$run_enrichment)


            shiny::observe({
                shiny::req(!is.null(gprof_result()$result), cancelOutput = TRUE)

                shiny::isolate({
                    output$enrichment_table <- DT::renderDT({
                        shiny::req(!is.null(gprof_result()))
                        df <- gprof_result()$result
                        df <- df[, 3:(ncol(df) - 2)]
                        return(df)
                    })

                    output$enrichment_plot <- plotly::renderPlotly({
                        shiny::req(!is.null(gprof_result()))
                        gprofiler2::gostplot(gprof_result())
                    })

                    output$download_enrichment <- shiny::downloadHandler(
                        filename = function() {
                            "enrichment_results.csv"
                        },
                        content = function(file) {
                            utils::write.csv(gprof_result()$result, file)
                        }
                    )
                })
            })
        }
    )
}
