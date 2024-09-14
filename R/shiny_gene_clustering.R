#' @importFrom dplyr %>% .data
NULL

###### UI ######

#' UI - Gene clustering
#'
#' @description Creates the UI interface for the Gene Clustering panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
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
        ),
        # TODO imrpove the graphical params
        shiny::splitLayout(
            cellWidths = "40px",
            shinyWidgets::dropdownButton(
                inputId = ns("heatmap_settings"),
                shiny::tagList(
                    shiny::sliderInput(
                        inputId = ns("scale_threshold"),
                        label = "Scale threshold",
                        min = 0.00, max = 1.00, value = 0.75, step = 0.005
                    ),
                    shiny::numericInput(
                        inputId = ns("colour_clipping"),
                        label = "Colour clipping value",
                        min = 0, value = 1000
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
        shiny::splitLayout(
            shiny::selectInput(
                inputId = ns("gene_clusters_first"),
                choices = NULL,
                label = "First group - select a gene clustering"
            ),
            shiny::selectInput(
                inputId = ns("gene_clusters_second"),
                choices = NULL,
                label = "Second group - select a gene clustering"
            ),
            shinyWidgets::radioGroupButtons(
                inputId = ns("heatmap_info"),
                label = "Base of the intersection",
                choices = c("Gene set", "Expressed cells"),
                selected = "Gene set"
            )
        ),
        shiny::plotOutput(ns("modules_heatmap"), height = "auto"),
        shiny::div(class = "empty-space")
    )
}

###### SERVER ######
#' Server - Gene clustering
#'
#' @description Creates the backend interface for the Gene Clustering panel
#' inside the Starlng Shiny application.
#'
#' @param id The id of the shiny module, used to access the UI elements.
#' @param filtered_genes A reactive expression that contains the filtered genes
#' that will be used for the clustering.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
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
            } 
            
            clust_df[["preloaded"]] <- shiny::reactive({
                shiny::req(env$preloaded_stable_modules)
                return(env$preloaded_stable_modules)
            })

            shiny::observe({
                shiny::req(filtered_genes())
                shiny::isolate({
                    output$clustering_info <- shiny::renderText({
                        paste0(
                            "Performing clustering on ", length(filtered_genes()), " genes. Modify the filters in the Gene Table section to increase / decrease the number of genes."
                        )
                    })
                })
            })

            realtime_clust_result <- shiny::reactive({
                shinyjs::disable("run_clustering")
                shinyjs::hide("ecc_threshold")
                shinyjs::hide("freq_threshold")

                print(Sys.time())
                res_list <- list(seq(from = input$res_start, to = input$res_stop, by = input$res_step))
                names(res_list) <- paste0(input$qfunc, "VertexPartition")
                cl_result <- group_by_clusters_general(clustering_pipeline(
                    # TODO change this after you reach a conclusion about the best option for the input matrix
                    embedding = get_feature_loading(
                        expr_matrix = read_gene_from_dense_h5(
                            gene_names = filtered_genes(),
                            matrix_h5_path = file.path("objects", "expression.h5"),
                            index_genes = env$genes[filtered_genes()],
                            check_intersect = FALSE
                        ),
                        approx = TRUE
                    ),
                    n_neighbours = input$n_neighbours,
                    graph_type = tolower(input$graph_type),
                    resolutions = res_list,
                    number_iterations = input$n_iterations,
                    number_repetitions = input$n_seeds,
                    merge_identical_partitions = FALSE
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
                    rownames(stb_df) <- filtered_genes()

                    if (is.null(cl_result)) {
                        return(NULL)
                    }

                    for (k in names(cl_result)) {
                        stb_df[[paste0("stable_modules_", k)]] <- factor(cl_result[[k]]$partitions[[1]]$mb)
                    }
                    stb_df <- stb_df[, 2:ncol(stb_df), drop = FALSE]

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
                    input$tabset %in% names(clust_df)
                )
                clust_df[[input$tabset]]()

                shiny::isolate({
                    module_options <- c("-")
                    if (!is.null(clust_df[[input$tabset]]())) {
                        module_options <- c("-", colnames(clust_df[[input$tabset]]()))
                    }
                    
                    shiny::updateSelectInput(
                        session,
                        inputId = "gene_clusters_options",
                        choices = module_options,
                        selected = module_options[1]
                    )
                })
            })

            available_modules <- shiny::reactive({
                tmp_available_modules <- list()
                for (tab_type in c("preloaded", "realtime", "custom")) {
                    if (tab_type == "realtime" && input$run_clustering == 0) {
                        next
                    }

                    if (tab_type == "custom" && input$create_family == 0) {
                        next
                    }
                    appr_clust_df <- clust_df[[tab_type]]()
                    if (is.null(appr_clust_df)) {
                        next
                    }

                    tmp_available_modules[[tab_type]] <- paste0(tab_type, "---", colnames(appr_clust_df))
                }
                tmp_available_modules <- unlist(tmp_available_modules)
                names(tmp_available_modules) <- NULL
                return(tmp_available_modules)
            })


            shiny::observe({
                tmp_available_modules <- available_modules()
                names(tmp_available_modules) <- NULL
                shiny::req(tmp_available_modules)
                shiny::updateSelectInput(
                    session,
                    inputId = "gene_clusters_first",
                    choices = tmp_available_modules,
                    selected = tmp_available_modules[1]
                )

                shiny::updateSelectInput(
                    session,
                    inputId = "gene_clusters_second",
                    choices = tmp_available_modules,
                    selected = tmp_available_modules[1]
                )
            })

            heatmap_information <- shiny::reactive({
                first_clusters <- input$gene_clusters_first
                second_clusters <- input$gene_clusters_second
                htmp_info_type <- input$heatmap_info

                shiny::req(first_clusters, second_clusters)
                scale_threshold <- input$scale_threshold
                colour_clipping <- input$colour_clipping
                # available_modules_cell_associations()
                shiny::isolate({
                    tab_first <- strsplit(first_clusters, "---")[[1]][1]
                    first_clusters <- strsplit(first_clusters, "---")[[1]][2]
                    first_clusters <- clust_df[[tab_first]]()[ , first_clusters]
                    unique_first <- unique(first_clusters)
                    if (all(grepl("^[0-9]+$", unique_first))) {
                        unique_first <- sort(as.numeric(unique_first))
                    }

                    tab_second <- strsplit(second_clusters, "---")[[1]][1]
                    second_clusters <- strsplit(second_clusters, "---")[[1]][2]
                    second_clusters <- clust_df[[tab_second]]()[ , second_clusters]
                    unique_second <- unique(second_clusters)
                    if (all(grepl("^[0-9]+$", unique_second))) {
                        unique_second <- sort(as.numeric(unique_second))
                    }

                    if (htmp_info_type == "Gene set") {
                        gene_contingency <- table(first_clusters, second_clusters)
                    } else {
                        gene_set_split_first <- split(rownames(clust_df[[tab_first]]()), first_clusters)
                        cell_list_first <- lapply(gene_set_split_first, function(gene_set) {
                            voting_result <- scale_min_max(voting_scheme(
                                expression_matrix = scale_min_max(read_gene_from_dense_h5(
                                    gene_names = gene_set,
                                    matrix_h5_path = file.path("objects", "expression.h5"),
                                    index_genes = env$genes[gene_set],
                                    check_intersect = FALSE
                                )),
                                genes = gene_set,
                                thresh_percentile = 0,
                                thresh_value = 0,
                                n_coexpressed_thresh = 1,
                                summary_function = mean
                            ))

                            index <- which(voting_result > scale_threshold)
                            return(names(voting_result)[index])
                        })

                        if (input$gene_clusters_first == input$gene_clusters_second) {
                            cell_list_second <- cell_list_first
                        } else {
                            gene_set_split_second <- split(rownames(clust_df[[tab_second]]()), second_clusters)
                            cell_list_second <- lapply(gene_set_split_second, function(gene_set) {
                                voting_result <- scale_min_max(voting_scheme(
                                    expression_matrix = scale_min_max(read_gene_from_dense_h5(
                                        gene_names = gene_set,
                                        matrix_h5_path = file.path("objects", "expression.h5"),
                                        index_genes = env$genes[gene_set],
                                        check_intersect = FALSE
                                    )),
                                    genes = gene_set,
                                    thresh_percentile = 0,
                                    thresh_value = 0,
                                    n_coexpressed_thresh = 1,
                                    summary_function = mean
                                ))

                                index <- which(voting_result > scale_threshold)
                                return(names(voting_result)[index])
                            })
                        }

                        gene_contingency <- matrix(0, nrow = length(cell_list_first), ncol = length(cell_list_second))
                        for (i in seq_along(cell_list_first)) {
                            for (j in seq_along(cell_list_second)) {
                                gene_contingency[i, j] <- length(intersect(cell_list_first[[i]], cell_list_second[[j]]))
                            }
                        }
                        rownames(gene_contingency) <- unique_first
                        colnames(gene_contingency) <- unique_second

                    }

                    clipped_gene_contingency <- gene_contingency
                    colour_clipping <- min(max(gene_contingency), colour_clipping)
                    clipped_gene_contingency[clipped_gene_contingency > colour_clipping] <- colour_clipping

                    return(ComplexHeatmap::Heatmap(
                        clipped_gene_contingency,
                        row_order = unique_first,
                        column_order = unique_second,
                        row_names_side = "left",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                            grid::grid.text(
                                label = gene_contingency[i, j],
                                x = x,
                                y = y,
                                just = "center",
                                gp = grid::gpar(fontsize = 10)
                            )
                        }
                    ))
                })
            })

            shiny::observe({
                shiny::req(heatmap_information())
                window_dims <- env$window_dim()
                shiny::isolate({
                    plt_height <- min(
                        window_dims[2] * 0.8,
                        250 + 40 * nrow(heatmap_information())
                    )

                    plt_width <- min(
                        window_dims[1] * 0.8,
                        250 + 40 * ncol(heatmap_information())
                    )
                    output$modules_heatmap <- shiny::renderPlot(
                        height = plt_height,
                        width = plt_width,
                        {
                            heatmap_information()
                        }
                    )
                })
            })

            shiny::observe({
                shinyjs::disable("select_module")

                shiny::req(input$gene_clusters_options)
                shiny::req(input$tabset %in% names(clust_df))
                module_options <- input$gene_clusters_options
                shiny::req(!is.null(clust_df[[input$tabset]]()))
                clust_df[[input$tabset]]()

                shiny::isolate({
                    shiny::req(length(colnames(clust_df[[input$tabset]]())) > 0)
                    shiny::req(!is.null(clust_df[[input$tabset]]()))
                    if (module_options != "-") {
                        shinyjs::enable("select_module")
                    } else {
                        module_options <- colnames(clust_df[[input$tabset]]())
                    }
                    shiny::req(all(module_options %in% colnames(clust_df[[input$tabset]]())))

                    df <- clust_df[[input$tabset]]()[ , module_options, drop = FALSE]
                    df <- cbind(env$moran_df[rownames(df), ], df)
                    output[[paste0(input$tabset, "_gene_clusters")]] <- DT::renderDataTable({
                        shiny::req(!is.null(clust_df[[input$tabset]]()))
                        DT::datatable(df, filter = "top")
                    })

                    output$download_gene_clusters <- shiny::downloadHandler(
                        filename = function() {
                            paste0(input$tabset, "_gene_clusters.csv")
                        },
                        content = function(file) {
                            utils::write.csv(df, file, row.names = TRUE)
                        }
                    )
                })
            })

            shiny::observe({
                shiny::req(input$gene_clusters_options, clust_df[[input$tabset]]())

                shiny::isolate({
                    shinyjs::disable("select_module")
                    module <- clust_df[[input$tabset]]()[, input$gene_clusters_options]
                    names(module) <- rownames(clust_df[[input$tabset]]())
                    env$chosen_modules(split(names(module), module))
                })
            }) %>% shiny::bindEvent(input$select_module)

            
        }
    )
}
