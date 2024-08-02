#' @importFrom dplyr %>% .data

###### UI ######
ui_gene_clustering <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::h2("Gene clustering", id = "gene_clustering_panel"),
        shiny::tabsetPanel(
            id = ns("tabset"),
            shiny::tabPanel(
                title = "Preloaded",
                value = "preloaded",
                shiny::wellPanel(DT::dataTableOutput(ns("preloaded_gene_clusters")))
            ),
            shiny::tabPanel(
                title = "Real time gene clustering",
                value = "realtime",
                shiny::wellPanel(
                    shiny::splitLayout(
                        shiny::numericInput(
                            inputId = ns("n_seeds"),
                            label = "# different seeds",
                            min = 1, max = 1000, value = 50, step = 1
                        ),
                        shiny::numericInput(
                            inputId = ns("res_start"),
                            label = "Resolution - start",
                            min = 0.01, max = 20, value = 0.1, step = 0.01
                        ),
                        shiny::numericInput(
                            inputId = ns("res_stop"),
                            label = "Resolution - stop",
                            min = 0.01, max = 20, value = 1, step = 0.01
                        ),
                        shiny::numericInput(
                            inputId = ns("res_step"),
                            label = "Resolution - step",
                            min = 0.01, max = 20, value = 0.1, step = 0.01
                        )
                    ),
                    shiny::splitLayout(
                        shiny::numericInput(
                            inputId = ns("n_neighbours"),
                            label = "# neighbours",
                            min = 1, max = 1000, value = 25, step = 1
                        ),
                        shiny::numericInput(
                            inputId = ns("n_iterations"),
                            label = "# clustering iters",
                            min = 1, max = 100, value = 5, step = 1
                        ),
                        shinyWidgets::radioGroupButtons(
                            inputId = ns("qfunc"),
                            label = "Quality function",
                            choices = c("RBConfiguration", "RBER", "CPM"),
                            selected = "RBConfiguration"
                        ),
                        shinyWidgets::radioGroupButtons(
                            inputId = ns("graph_type"),
                            label = "Graph type",
                            choices = c("SNN", "NN"),
                            selected = "SNN"
                        )
                    ),
                    shiny::textOutput(ns("clustering_info")),
                    shiny::actionButton(
                        inputId = ns("run_clustering"),
                        label = "Cluster genes!",
                        class = "btn-danger"
                    ),
                    shiny::p("Filter the stable modules using the ECC and frequency thresholds"),
                    shiny::splitLayout(
                        shiny::sliderInput(
                            inputId = ns("ecc_threshold"),
                            label = "ECC threshold",
                            min = 0.00, max = 1.00, value = 0.9, step = 0.01
                        ),
                        shiny::sliderInput(
                            inputId = ns("freq_threshold"),
                            label = "#appearances threshold",
                            min = 0, max = 100, value = 30, step = 1
                        )
                    ),
                    DT::dataTableOutput(ns("realtime_gene_clusters"))
                )
            ),
            shiny::tabPanel(
                title = "Custom gene modules",
                value = "custom",

                shiny::wellPanel(
                    shiny::p("How many gene modules do you need?"),
                    shiny::splitLayout(
                        shiny::sliderInput(ns("n_families"), "Number of gene modules", value = 1, min = 1, max = 20, step = 1),
                        shiny::textInput(ns("family_name"), "Name of the gene clustering")
                    ),
                    shiny::uiOutput(ns("gene_families")),
                    shiny::actionButton(
                        inputId = ns("create_family"),
                        label = "Create gene modules!",
                        class = "btn-danger"
                    ),
                    DT::dataTableOutput(ns("custom_gene_clusters"))
                ),
            )
        ),
        shiny::fluidRow(
            shiny::column(3,
                shiny::downloadButton(
                    outputId = ns("download_gene_clusters"),
                    label = "Download clusters as CSV"
                )
            ),
            shiny::column(
                2,
                shiny::selectInput(
                    inputId = ns("gene_clusters_options"),
                    choices = NULL,
                    label = "Select a gene module"
                )
            ),
            shiny::column(
                2,
                shiny::actionButton(
                    inputId = ns("select_module"),
                    label = "Fix the clustering!",
                    class = "btn-danger"
                )
            ),
            shiny::column(3)
        )

    )
}

###### SERVER ######
server_gene_clustering <- function(id, filtered_genes) {
    shiny::moduleServer(
        id,
        function(input, output, session) {
            clust_df <- list()
            shinyjs::hide("ecc_threshold")
            shinyjs::hide("freq_threshold")

            if (is.null(env$preloaded_stable_modules)) {
                shiny::updateTabsetPanel(
                    session,
                    inputId = "tabset",
                    selected = "custom"
                )

                shiny::hideTab(
                    inputId = "tabset",
                    target = "preloaded"
                )
            } else {
                clust_df[["preloaded"]] <- env$preloaded_stable_modules
            }

            shiny::observe({
                shiny::req(filtered_genes())
                shiny::isolate({
                    output$clustering_info <- shiny::renderText({
                        paste0(
                            "Performing clustering on ", length(filtered_genes()), " genes. Modify the filters in the Gene information table section to increase / decrease the number of genes."
                        )
                    })
                })
            })

            realtime_clust_result <- shiny::reactive({
                shinyjs::disable("run_clustering")
                shinyjs::hide("ecc_threshold")
                shinyjs::hide("freq_threshold")
                print(Sys.time())
                cl_result <- group_by_clusters_general(clustering_pipeline(
                    # TODO change this after you reach a conclusion about the best option for the input matrix
                    embedding = get_feature_loading(
                        read_gene_from_dense_h5(
                            gene_name = filtered_genes(),
                            matrix_h5_path = file.path("objects", "expression.h5"),
                            index_genes = env$genes[filtered_genes()],
                            check_intersect = FALSE
                        )
                    ),
                    n_neighbours = input$n_neighbours,
                    graph_type = tolower(input$graph_type),
                    resolutions = seq(from = input$res_start, to = input$res_stop, by = input$res_step),
                    quality_functions = paste0(input$qfunc, "VertexPartition"),
                    number_interations = input$n_iterations,
                    number_repetitions = input$n_seeds
                ), 3)

                cl_result <- get_clusters_consistency(cl_result)[[1]][[1]]
                print(Sys.time())
                shinyjs::enable("run_clustering")
                shinyjs::show("ecc_threshold")
                shinyjs::show("freq_threshold")
                return(cl_result)

            }) %>% shiny::bindEvent(input$run_clustering)

            clust_df[["realtime"]] <- shiny::reactive({
                shiny::req(realtime_clust_result(), input$ecc_threshold, input$freq_threshold)

                shiny::isolate({
                    cl_result <- ClustAssess::choose_stable_clusters(
                        realtime_clust_result()$k,
                        ecc_threshold = input$ecc_threshold,
                        freq_threshold = input$freq_threshold
                    )

                    stb_df <- data.frame(genes = filtered_genes())
                    for (k in names(cl_result)) {
                        stb_df[[paste0("stable_modules_", k)]] <- factor(cl_result[[k]]$partitions[[1]]$mb)
                    }
                    rownames(stb_df) <- filtered_genes()
                    stb_df <- stb_df[,2:ncol(stb_df)]

                    return(stb_df)
                })
            })

            shiny::observe({
                shiny::req(input$n_families)
                shiny::isolate({
                    ns <- shiny::NS(id)
                    output$gene_families <- shiny::renderUI({
                        do.call(
                            shiny::tagList,
                            lapply(seq_len(input$n_families), function(i) {
                                shiny::splitLayout(
                                    cellWidths = c("15%", "85%"),
                                    shiny::textInput(
                                        inputId = ns(paste0("family_", i)),
                                        label = "Gene module name",
                                        value = paste0("module_", i)
                                    ),
                                    shiny::selectizeInput(
                                        inputId = ns(paste0("family_genes_", i)),
                                        label = "Select genes",
                                        choices = NULL,
                                        width = "95%",
                                        multiple = TRUE
                                    )
                                )
                            })
                        )
                    })

                    if (input$family_name == "") {
                        shiny::updateTextInput(
                            session,
                            inputId = "family_name",
                            value = paste0("custom_", input$n_families, "_modules")
                        )
                    }
                    
                    for (i in seq_len(input$n_families)) {
                        shiny::updateSelectizeInput(
                            session,
                            inputId = paste0("family_genes_", i),
                            choices = names(get("genes", envir = env)),
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
            })

            family_genes <- shiny::reactive({
                shiny::req(input$n_families > 0)
                family_names <- rep(NA, input$n_families)
                for (i in seq_len(input$n_families)) {
                    shiny::req(
                        !is.null(input[[paste0("family_", i)]]),
                        input[[paste0("family_", i)]] != ""
                    )
                    family_names[i] <- input[[paste0("family_", i)]]
                }

                family_list <- lapply(seq_len(input$n_families), function(i) {
                    input[[paste0("family_genes_", i)]]
                })
                names(family_list) <- family_names
                family_list
            }) #%>% shiny::debounce(1000)

            shiny::observe({
                shiny::req(length(family_genes()) > 0)
                shinyjs::disable("create_family")

                shiny::isolate({
                    all_genes <- c()
                    for (i in seq_len(input$n_families)) {
                        module_gene <- family_genes()[[i]]
                        shiny::req(!is.null(module_gene))
                        shiny::req(length(module_gene) > 0)
                        module_gene <- setdiff(module_gene, all_genes)
                        shiny::updateSelectizeInput(
                            session,
                            inputId = paste0("family_genes_", i),
                            selected = module_gene
                        )
                        shiny::req(length(module_gene) > 0)
                        all_genes <- c(all_genes, module_gene)
                    }
                    shinyjs::enable("create_family")
                })
            })

            clust_df[["custom"]] <- shiny::reactive({
                shiny::req(length(family_genes()) > 0, input$family_name)

                shiny::isolate({
                    all_genes <- do.call(c, family_genes())
                    module_name <- rep(NA, length(all_genes))
                    index <- 1
                    for (i in seq_along(family_genes())) {
                        ngenes <- length(family_genes()[[i]])
                        index_update <- seq(from = index, by = 1, length.out = ngenes)
                        module_name[index_update] <- names(family_genes())[i]
                        index <- index + ngenes
                    }

                    stb_df <- data.frame(genes = all_genes)
                    stb_df[[input$family_name]] <- factor(module_name)
                    rownames(stb_df) <- all_genes
                    stb_df <- stb_df[, 2:ncol(stb_df), drop = FALSE]

                    return(stb_df)
                })
            }) %>% shiny::bindEvent(input$create_family)

            shiny::observe({
                shiny::req(
                    input$tabset,
                    input$tabset %in% names(clust_df),
                    clust_df[[input$tabset]]()
                )

                shiny::isolate({
                    output[[paste0(input$tabset, "_gene_clusters")]] <- DT::renderDataTable(
                        DT::datatable(
                            clust_df[[input$tabset]](),
                            filter = "top"
                        )
                    )

                    module_options <- colnames(clust_df[[input$tabset]]())
                    shiny::updateSelectInput(
                        session,
                        inputId = "gene_clusters_options",
                        choices = module_options,
                        selected = module_options[1]
                    )

                    output$download_gene_clusters <- shiny::downloadHandler(
                        filename = function() {
                            paste0(input$tabset, "_gene_clusters.csv")
                        },
                        content = function(file) {
                            write.csv(clust_df[[input$tabset]](), file, row.names = TRUE)
                        }
                    )
                })
            })

            shiny::observe({
                shinyjs::disable("gene_clusters_options")
                shinyjs::disable("select_module")

                shiny::req(input$gene_clusters_options)
                shiny::req(input$tabset %in% names(clust_df))
                shiny::req(all(input$gene_clusters_options %in% colnames(clust_df[[input$tabset]]())))

                shinyjs::enable("select_module")
                shinyjs::enable("gene_clusters_options")
            })

            shiny::observe({
                shiny::req(input$gene_clusters_options, clust_df[[input$tabset]]())

                shiny::isolate({
                    module <- clust_df[[input$tabset]]()[, input$gene_clusters_options]
                    names(module) <- rownames(clust_df[[input$tabset]]())
                    env$chosen_modules(split(names(module), module))
                })
            }) %>% shiny::bindEvent(input$select_module)
        }
    )
}
