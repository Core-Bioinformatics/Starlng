#' @importFrom dplyr %>% .data
NULL

###### UI ######
ui_module_enrichment <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Module enrichment", id = "enrichment_panel"),
        shiny::fluidRow(
            shiny::column(1, shinyWidgets::dropdownButton(
                inputId = ns("enrichment_settings"),
                shiny::tagList(
                    shiny::numericInput(
                        inputId = ns("p_value_threshold"),
                        label = "P-value threshold",
                        min = 0, max = 1, value = 0.05, step = 0.001
                    ),
                    shiny::selectInput(
                        inputId = ns("correction_method"),
                        label = "Correction method",
                        choices = c("gSCS", "fdr", "bonferroni"),
                        selected = "fdr"
                    ),
                    shiny::selectInput(
                        inputId = ns("background_type"),
                        label = "Background genes",
                        choices = c("Annotated" = "annotated", "All expressed genes" = "all_expressed", "All clustered genes" = "all_clustered"),
                        selected = "annotated"
                    ),
                    shinyWidgets::prettySwitch(
                        inputId = ns("combine_modules"),
                        label = "Combine selected modules",
                        status = "success",
                        value = FALSE,
                        fill = TRUE
                    ),
                    shinyWidgets::prettySwitch(
                        inputId = ns("order_genes"),
                        label = "Order genes by hub score",
                        status = "success",
                        value = FALSE,
                        fill = TRUE
                    )
                ),
                label = "",
                icon = shiny::icon("cog"),
                status = "success",
                size = "sm",
                circle = TRUE
            )),
            shiny::column(2, shiny::selectizeInput(
                inputId = ns("module_enrichment"),
                label = "Select gene modules",
                choices = NULL,
                selected = NULL,
                multiple = TRUE,
                options = list(
                    title = "Select gene modules",
                    placeholder = "Select gene modules",
                    plugins = list("remove_button", "drag_drop"),
                    delimiter = ",",
                    create = TRUE
                ),
            )),
            shiny::column(2, shiny::numericInput(
                inputId = ns("n_top_genes"),
                label = "Select number of top hub genes (-1 for all)",
                value = -1,
                min = -1,
                step = 1
            )),
            shiny::column(2, shinyWidgets::pickerInput(
                inputId = ns("enrichment_sources"),
                label = "Select data sources",
                choices = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "TF", "MIRNA", "CORUM", "HP", "HPA", "WP"),
                selected = c("GO:BP"),
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
        shiny::selectizeInput(
            inputId = ns("select_module_for_plot"),
            label = "Select module for enrichment plot",
            choices = c(""),
            multiple = FALSE,
            options = list(
                title = "Select module for enrichment plot",
                size = 10,
                width = "90%"
            )
        ),
        plotly::plotlyOutput(ns("enrichment_plot")),
        DT::dataTableOutput(ns("enrichment_table")),
        shiny::downloadButton(
            outputId = ns("download_enrichment"),
            label = "Download enrichment results as CSV"
        ),
        shiny::h3("Enrichment DotPlot"),
        shiny::splitLayout(
            cellWidths = c("40px", "40px", "300px", "300px"),
            shinyWidgets::dropdownButton(
                inputId = ns("dotplot_settings"),
                shiny::tagList(
                    shiny::sliderInput(
                        inputId = ns("dotplot_point_size_range"),
                        label = "Dotplot point size range",
                        min = 1,
                        max = 15,
                        value = c(3, 10),
                        step = 0.5
                    ),
                    shiny::numericInput(
                        inputId = ns("dotplot_font_size"),
                        label = "Dotplot font size",
                        value = 10,
                        min = 6,
                        max = 30,
                        step = 1
                    )
                ),
                label = "",
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "enrichment_dotplot", "enrichment_dotplot"),
            shiny::numericInput(
                inputId = ns("n_top_enriched_terms"),
                label = "# top enriched terms per module",
                value = 10,
                min = 1,
                step = 1
            ),

            shiny::selectInput(
                inputId = ns("select_enrichment_source"),
                label = "Select data source for enrichment plot",
                choices = "",
                selected = NULL,
                multiple = FALSE
            )
        ),
        shiny::plotOutput(ns("enrichment_dotplot"), height = "auto")
    )
}

ui_transcription_factor_activity <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::h3("Activity Table"),
        DT::dataTableOutput(ns("tf_activity_table"), width = "100%"),
        shiny::downloadButton(
            outputId = ns("download_tf_activity_table"),
            label = "Download TF activity table as CSV"
        ),
        shiny::h3("Activity Heatmap"),
        shiny::splitLayout(
            cellWidths = c("40px", "40px", "300px", "300px", "300px"),
            shinyWidgets::dropdownButton(
                inputId = ns("tf_activity_bubbleplot_settings"),
                shiny::tagList(
                    shiny::numericInput(
                        inputId = ns("intersection_size_cap"),
                        label = "# genes cap for bubbleplot",
                        min = 1,
                        value = 100,
                        step = 1
                    ),
                    shiny::numericInput(
                        inputId = ns("hub_size_cap"),
                        label = "# hub genes cap for bubbleplot",
                        min = 1,
                        value = 20,
                        step = 1
                    ),
                    shiny::sliderInput(
                        inputId = ns("tf_bubbleplot_point_range"),
                        label = "Bubbleplot point size range",
                        min = 0,
                        max = 15,
                        value = c(1, 3),
                        step = 0.1
                    ),
                    shiny::numericInput(
                        inputId = ns("tf_activity_bubbleplot_font_size"),
                        label = "Bubbleplot font size",
                        value = 12,
                        min = 6,
                        max = 30,
                        step = 1
                    )
                ),
                label = "",
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "tf_activity_bubbleplot", "tf_activity_bubbleplot"),
            shiny::numericInput(
                inputId = ns("n_top_tfs_bubbleplot"),
                label = "# top TFs to show in bubbleplot",
                value = 20,
                min = 1,
                step = 1
            ),
            shiny::selectizeInput(
                inputId = ns("select_tf_activity_module"),
                label = "Select module for TF activity bubbleplot",
                choices = "",
                selected = NULL,
                multiple = TRUE,
                options = list(
                    title = "Select module for TF activity bubbleplot",
                    size = 10,
                    width = "90%",
                    plugins = list("remove_button", "drag_drop"),
                    delimiter = ",",
                    create = TRUE
                )
            ),
            shiny::numericInput(
                inputId = ns("n_hub_genes_tf_network"),
                label = "# hub genes to show in TF network",
                value = 2,
                min = 1,
                step = 1
            )
        ),
        shiny::plotOutput(ns("tf_activity_bubbleplot"), height = "auto"),
        shiny::h3("TF Network Visualization"),
        shiny::splitLayout(
            cellWidths = c("40px", "40px", "300px"),
            shinyWidgets::dropdownButton(
                inputId = ns("tf_network_settings"),
                shiny::tagList(
                    shiny::sliderInput(
                        inputId = ns("node_size_tf_network"),
                        label = "TF network node size",
                        min = 0.1,
                        max = 10,
                        value = 4,
                        step = 0.1
                    ),
                    shiny::sliderInput(
                        inputId = ns("edge_width_tf_network"),
                        label = "TF network edge width",
                        min = 0.01,
                        max = 5,
                        value = 0.1,
                        step = 0.05
                    ),
                    shiny::numericInput(
                        inputId = ns("tf_network_label_font_size"),
                        label = "TF network label font size",
                        value = 4,
                        min = 1,
                        max = 20,
                        step = 0.5
                    ),
                    shiny::numericInput(
                        inputId = ns("tf_network_axis_font_size"),
                        label = "TF network axis font size",
                        value = 15,
                        min = 1,
                        max = 20,
                        step = 0.5
                    ),
                    shinyWidgets::prettySwitch(
                        inputId = ns("exclude_non_hub_genes"),
                        label = "Exclude non-hub genes from network plot",
                        status = "success",
                        value = FALSE,
                        fill = TRUE
                    )
                ),
                label = "",
                circle = TRUE,
                status = "success",
                size = "sm",
                icon = shiny::icon("cog")
            ),
            gear_download(ns, "tf_network_plot", "tf_network_plot"),
            shiny::numericInput(
                inputId = ns("n_tfs_per_module"),
                label = "# TFs per module for network plot",
                value = 3,
                min = 1,
                step = 1
            )
        ),
        shiny::plotOutput(ns("tf_network_plot"), height = "auto")
    )
}

#' UI - Functional Assessment
#'
#' @description Creates the UI interface for the Functional Assessment panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
ui_functional_assessment <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::h2("Transcription factor activity", id = ns("tf_activity_panel")),
        ui_transcription_factor_activity(ns("tf_activity")),
        shiny::h2("Enrichment analysis", id = ns("enrichment_panel")),
        ui_module_enrichment(ns("module_enrichment"))
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
                module_summary <- env$modules_stats_summary()
                ordered_modules <- names(env$chosen_modules())
                if (!is.null(module_summary) && nrow(module_summary) > 0) {
                    ordered_modules <- rownames(module_summary)
                    ordered_modules <- ordered_modules[ordered_modules %in% names(env$chosen_modules())]
                    ordered_modules <- c(ordered_modules, setdiff(names(env$chosen_modules()), ordered_modules))
                }
                shiny::updateSelectizeInput(
                    session,
                    inputId = "module_enrichment",
                    choices = ordered_modules,
                    selected = ordered_modules[1],
                    server = TRUE
                )
            })

            shiny::observe({
                shiny::req(!is.null(input$module_enrichment), length(input$module_enrichment) > 0)
                shinyjs::show("run_enrichment")
            })

            gene_sets <- shiny::reactive({
                selected_modules <- input$module_enrichment
                n_top_genes <- input$n_top_genes
                gene_scores <- env$gene_hub_scores()
                combine_modules <- input$combine_modules
                module_summary <- env$modules_stats_summary()

                shiny::req(selected_modules, gene_scores)
                selected_modules <- as.character(selected_modules)
                module_list <- lapply(selected_modules, function(module_name) {
                    df <- gene_scores %>%
                        dplyr::filter(.data$module == module_name) %>%
                        dplyr::arrange(dplyr::desc(.data$combined_score)) %>%
                        rownames()
                    return(df)
                })
                names(module_list) <- selected_modules
                if (n_top_genes > 1) {
                    for (i in seq_along(module_list)) {
                        module_list[[i]] <- head(module_list[[i]], n_top_genes)
                    }
                }

                if (combine_modules) {
                    module_list <- list(unique(unlist(module_list, use.names = FALSE)))
                    names(module_list) <- paste(selected_modules, collapse = "+")
                }

                return(module_list)
            }) %>% shiny::bindEvent(input$run_enrichment)

            gprof_result <- shiny::reactive({
                input_gene_sets <- gene_sets()
                p_val_thresh <- input$p_value_threshold
                correction_method <- input$correction_method
                enrichment_sources <- input$enrichment_sources
                background_type <- input$background_type
                order_genes <- input$order_genes

                shiny::isolate({
                    shiny::req(input_gene_sets, p_val_thresh, correction_method, enrichment_sources, background_type)
                    shinyjs::disable("run_enrichment")
                    shinyjs::hide("download_enrichment")

                    bg_genes <- switch(
                        background_type,
                        "all_expressed" = names(env$genes),
                        "all_clustered" = unique(unlist(env$chosen_modules(), use.names = FALSE)),
                        NULL
                    )

                    annotation_domain <- if (is.null(bg_genes)) {
                        "annotated"
                    } else {
                        "custom_annotated"
                    }

                    # NOTE: in case it's slow, try https://bioconductor.org/packages//release/bioc/html/clusterProfiler.html
                    gprf_res <- lapply(names(input_gene_sets), function(module_name) {
                        query_genes <- input_gene_sets[[module_name]]
                        current_res <- gprofiler2::gost(
                            query = query_genes,
                            ordered_query = order_genes,
                            sources = enrichment_sources,
                            organism = env$organism,
                            evcodes = TRUE,
                            domain_scope = annotation_domain,
                            custom_bg = bg_genes,
                            correction_method = correction_method,
                            significant = FALSE
                        )

                        if (is.null(current_res$result)) {
                            return(NULL)
                        }

                        current_res$result <- current_res$result[current_res$result$p_value < p_val_thresh, , drop = FALSE]
                        if (nrow(current_res$result) == 0) {
                            return(NULL)
                         }
                        current_res$result$parents <- sapply(current_res$result$parents, toString)
                        current_res$result$module <- module_name
                        return(current_res)
                    })
                    names(gprf_res) <- names(input_gene_sets)

                    if (any(sapply(gprf_res, function(res) !is.null(res$result) && nrow(res$result) > 0))) {
                        shinyjs::show("download_enrichment")
                        output$enrichment_message <- shiny::renderText("")
                    } else {
                        output$enrichment_message <- shiny::renderText("No enrichment term found for current settings.")
                    }
                    shinyjs::enable("run_enrichment")

                    return(gprf_res)
                })
            }) %>% shiny::bindEvent(gene_sets())

            shiny::observe({
                gprf_res <- gprof_result()
                shiny::req(gprf_res)
                null_results_mask <- sapply(gprf_res, function(res) is.null(res$result) || nrow(res$result) == 0)
                gprf_res <- gprf_res[!null_results_mask]
                all_sources <- unique(unlist(lapply(gprf_res, function(res) {
                    return(unique(res$result$source))
                }), use.names = FALSE))

                shiny::updateSelectInput(
                    session,
                    inputId = "select_enrichment_source",
                    choices = all_sources,
                    selected = all_sources[1]
                )
                
                shiny::updateSelectizeInput(
                    session,
                    inputId = "select_module_for_plot",
                    choices = names(gprf_res),
                    selected = names(gprf_res)[1]
                )

                combined_df <- do.call(rbind, lapply(gprf_res, function(res) {
                    if (!is.null(res$result)) {
                        return(res$result)
                    } else {
                        return(NULL)
                    }
                }))
                output$enrichment_table <- DT::renderDT({
                    shiny::req(!is.null(combined_df), nrow(combined_df) > 0, cancelOutput = TRUE)
                    df <- combined_df
                    df$evidence_codes <- NULL
                    df$intersection <- NULL
                    df$source <- factor(df$source)
                    return(DT::datatable(df, filter = "top"))
                })

                output$download_enrichment <- shiny::downloadHandler(
                    filename = function() {
                        "enrichment_results.csv"
                    },
                    content = function(file) {
                        utils::write.csv(combined_df, file, row.names = FALSE)
                    }
                )
            }) %>% shiny::bindEvent(gprof_result())

            shiny::observe({
                wdim <- env$window_dim()
                gprf_res <- gprof_result()
                selected_module <- input$select_module_for_plot

                shiny::isolate({
                    shiny::req(gprf_res, selected_module %in% names(gprf_res), wdim, cancelOutput = TRUE)
                    gprf_res <- gprf_res[[selected_module]]
                    shiny::req(!is.null(gprf_res$result), nrow(gprf_res$result) > 0, cancelOutput = TRUE)

                    output$enrichment_plot <- plotly::renderPlotly(
                        {
                            shiny::req(!is.null(gprf_res))
                            gprofiler2::gostplot(gprf_res)
                        }
                    )
                })
            }) %>% shiny::bindEvent(input$select_module_for_plot)

            shiny::observe({
                gprf_res <- gprof_result()
                selected_source <- input$select_enrichment_source
                n_top_terms <- input$n_top_enriched_terms
                dotplot_point_size_range <- input$dotplot_point_size_range
                dotplot_font_size <- input$dotplot_font_size

                shiny::req(gprf_res, selected_source, n_top_terms, dotplot_point_size_range, dotplot_font_size, cancelOutput = TRUE)

                combined_df <- do.call(rbind, lapply(gprf_res, function(res) {
                    if (!is.null(res$result)) {
                        return(res$result)
                    }
                    return(NULL)
                }))

                shiny::req(!is.null(combined_df), nrow(combined_df) > 0, cancelOutput = TRUE)
                order_modules <- names(gprf_res)
                combined_df$module <- factor(combined_df$module, levels = order_modules)
                combined_df <- combined_df[combined_df$source == selected_source, , drop = FALSE]
                shiny::req(nrow(combined_df) > 0, cancelOutput = TRUE)

                dotplot_obj <- plot_enrichment_top_terms(
                    enrichment_result = combined_df,
                    top_n = n_top_terms,
                    colour_column = if ("intersection_size" %in% names(combined_df)) "intersection_size" else NULL,
                    point_size_range = dotplot_point_size_range,
                    font_size = dotplot_font_size
                )

                output$enrichment_dotplot <- shiny::renderPlot(
                    height = env$window_dim()[2] * 0.6,
                    width = env$window_dim()[1] * 0.9,
                {
                    dotplot_obj
                })

                output$download_enrichment_dotplot <- shiny::downloadHandler(
                    filename = function() {
                        req_gear_download(session, "enrichment_dotplot")
                        paste(input$filename_enrichment_dotplot, tolower(input$filetype_enrichment_dotplot), sep = ".")
                    },
                    content = function(file) {
                        req_gear_download(session, "enrichment_dotplot")
                        device_fun <- save_filetypes[[input$filetype_enrichment_dotplot]]
                        device_fun(file, width = input$width_enrichment_dotplot, height = input$height_enrichment_dotplot)
                        on.exit(grDevices::dev.off(), add = TRUE)
                        print(dotplot_obj)
                    }
                )
            })
        }
    )
}

server_transcription_factor_activity <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            tf_stats <- shiny::reactive({
                module_list <- env$chosen_modules()
                shiny::req(module_list, length(module_list) > 0, cancelOutput = TRUE)
                nmodules <- length(module_list)
                tfs <- rhdf5::h5read(
                    file.path("objects", "module_summaries.h5"),
                    paste0(nmodules, "/tfs")
                )
                tfs$module <- as.character(tfs$module)
                tfs$tf <- as.character(tfs$tf)
                tfs$n_genes <- as.integer(tfs$n_genes)
                tfs$genes <- as.character(tfs$genes)

                gene_hubs <- gene_hub_scores()
                shiny::req(gene_hubs, cancelOutput = TRUE)

                return(add_tf_hub_stats(tfs, gene_hubs))
            })

            gene_hub_scores <- shiny::reactive({
                gene_hub_scores <- env$gene_hub_scores()
                shiny::req(gene_hub_scores, cancelOutput = TRUE)
                n_top_hub_genes <- input$n_hub_genes_tf_network
                shiny::req(n_top_hub_genes, cancelOutput = TRUE)
                gene_hub_scores$gene <- rownames(gene_hub_scores)

                gene_hub_scores <- gene_hub_scores %>%
                    dplyr::group_by(.data$module) %>%
                    dplyr::arrange(dplyr::desc(.data$combined_score)) %>%
                    dplyr::slice_head(n = n_top_hub_genes) %>%
                    dplyr::ungroup()
                return(gene_hub_scores)
            })

            shiny::observe({
                tfs <- tf_stats()
                shiny::req(tfs, nrow(tfs) > 0, cancelOutput = TRUE)
                module_ordering <- env$modules_stats_summary()
                shiny::req(module_ordering, cancelOutput = TRUE)
                module_ordering <- intersect(rownames(module_ordering), unique(tfs$module))
                tfs$module <- factor(tfs$module, levels = module_ordering)

                shiny::updateSelectizeInput(
                    session,
                    inputId = "select_tf_activity_module",
                    choices = module_ordering,
                    selected = module_ordering,
                    server = TRUE
                )

                output$tf_activity_table <- DT::renderDT({
                    DT::datatable(tfs, filter = "top")
                })

                output$download_tf_activity_table <- shiny::downloadHandler(
                    filename = function() {
                        "tf_activity_table.csv"
                    },
                    content = function(file) {
                        utils::write.csv(tf_stats(), file, row.names = FALSE)
                    }
                )
            }) %>% shiny::bindEvent(tf_stats())

            shiny::observe({
                tfs <- tf_stats()
                shiny::req(tfs, nrow(tfs) > 0, cancelOutput = TRUE)
                selected_modules <- input$select_tf_activity_module
                n_top_tfs <- input$n_top_tfs_bubbleplot
                font_size <- input$tf_activity_bubbleplot_font_size
                point_size_range <- input$tf_bubbleplot_point_range
                
                cap_value_intersection <- input$intersection_size_cap
                cap_value_hub <- input$hub_size_cap

                shiny::req(selected_modules, n_top_tfs, font_size, point_size_range, cap_value_intersection, cap_value_hub, cancelOutput = TRUE)
                tfs <- tfs %>% dplyr::filter(.data$module %in% selected_modules)
                tfs$module <- factor(tfs$module, levels = selected_modules)

                bubbleplot <- plot_module_tfs_bubbleplot(
                    tf_stats = tfs, 
                    n_top = n_top_tfs,
                    cap_value_intersection = cap_value_intersection,
                    cap_value_hub = cap_value_hub,
                    font_size = font_size,
                    point_size_range = point_size_range
                )

                output$tf_activity_bubbleplot <- shiny::renderPlot(
                    height = env$window_dim()[2] * 0.6,
                    width = env$window_dim()[1] * 0.9,
                    {
                    bubbleplot
                })

                output$download_tf_activity_bubbleplot <- shiny::downloadHandler(
                    filename = function() {
                        req_gear_download(session, "tf_activity_bubbleplot")
                        paste(input$filename_tf_activity_bubbleplot, tolower(input$filetype_tf_activity_bubbleplot), sep = ".")
                    },
                    content = function(file) {
                        req_gear_download(session, "tf_activity_bubbleplot")
                        device_fun <- save_filetypes[[input$filetype_tf_activity_bubbleplot]]
                        device_fun(file, width = input$width_tf_activity_bubbleplot, height = input$height_tf_activity_bubbleplot)
                        on.exit(grDevices::dev.off(), add = TRUE)
                        print(bubbleplot)
                    }
                )
            })

            tf_sub_module_adjacency <- shiny::reactive({
                tfs <- tf_stats()
                closest_node <- env$closest_node_per_module()
                trajectory_obj <- env$trajectory_object
                shiny::req(tfs, closest_node, trajectory_obj, nrow(tfs) > 0, cancelOutput = TRUE)

                closest_node <- closest_node[unique(tfs$module)]
                get_module_transitions(
                    trajectory_obj,
                    closest_node
                )
            })

            tf_gene_graph <- shiny::reactive({
                tfs <- tf_stats()
                gene_hubs <- gene_hub_scores()
                tf_adj <- tf_sub_module_adjacency()

                n_tfs_per_module <- input$n_tfs_per_module

                shiny::req(tfs, gene_hubs, tf_adj, n_tfs_per_module, cancelOutput = TRUE)

                return(get_tf_gene_network(
                    tf_stats = tfs,
                    module_adjacency = tf_adj,
                    top_n_factors = n_tfs_per_module,
                    hub_genes = gene_hubs
                ))
            })

            shiny::observe({
                tf_gene_g <- tf_gene_graph()
                edge_width <- input$edge_width_tf_network
                node_size <- input$node_size_tf_network
                font_size <- input$tf_network_label_font_size
                axis_text_size <- input$tf_network_axis_font_size * 0.8
                exclude_non_hub_genes <- input$exclude_non_hub_genes

                shiny::req(tf_gene_g, edge_width, node_size, font_size, axis_text_size, !is.null(exclude_non_hub_genes), cancelOutput = TRUE)

                graph_plot <- plot_module_tfs_ggraph(
                    module_tf_g = tf_gene_g,
                    node_size = node_size,
                    edge_width = edge_width,
                    label_size = font_size,
                    axis_text_size = axis_text_size,
                    exclude_non_hub_genes = exclude_non_hub_genes
                )


                output$tf_network_plot <- shiny::renderPlot(
                    height = env$window_dim()[2] * 0.6,
                    width = env$window_dim()[1] * 0.9,
                    {
                        print(graph_plot)
                    })

                output$download_tf_network_plot <- shiny::downloadHandler(
                    filename = function() {
                        req_gear_download(session, "tf_network_plot")
                        paste(input$filename_tf_network_plot, tolower(input$filetype_tf_network_plot), sep = ".")
                    },
                    content = function(file) {
                        req_gear_download(session, "tf_network_plot")
                        device_fun <- save_filetypes[[input$filetype_tf_network_plot]]
                        device_fun(file, width = input$width_tf_network_plot, height = input$height_tf_network_plot)
                        on.exit(grDevices::dev.off(), add = TRUE)
                        print(graph_plot)
                    }
                )
            })
        }
    )
}

#' Server - Functional Assessment
#'
#' @description Creates the backend interface for the Functional Assessment panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
server_functional_assessment <- function(id) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            server_module_enrichment("module_enrichment")
            server_transcription_factor_activity("tf_activity")
        }
    )
}
