#' @importFrom dplyr %>% .data

generate_discrete_colours <- function(metadata_df,
                                      existing_colours_list = NULL,
                                      max_n_colors = 40,
                                      qualpal_pallete = "pretty_dark") {
    if (is.null(existing_colours_list)) {
        existing_colours_list <- list()
    }

    for (mtd_names in colnames(metadata_df)) {
        if (!(inherits(metadata_df[[mtd_names]], c("factor", "character")))) {
            next
        }

        n_colors <- length(unique(metadata_df[[mtd_names]]))
        max_n_colors <- max(n_colors, max_n_colors)
    }

    for (i in seq_len(max_n_colors)) {
        col_name <- as.character(i)
        if (col_name %in% names(existing_colours_list)) {
            next
        }

        if (i == 1) {
            existing_colours_list[[col_name]] <- "#025147"
            next
        }

        existing_colours_list[[col_name]] <- qualpalr::qualpal(i, qualpal_pallete)$hex
    }

    return(existing_colours_list)
}

generate_continuous_colours <- function(existing_colours_list = NULL) {
    predefined_colours <- list(
        "viridis" = viridis::viridis(100),
        "plasma" = viridis::plasma(100),
        "white_red" = c("grey85", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84",
            "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"),
        "blue_red" = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
            "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"),
        "blue_red2" = grDevices::colorRampPalette(c("blue", "white", "red"))(100)
    )

    if (is.null(existing_colours_list)) {
        existing_colours_list <- list()
    }

    for (col_name in names(predefined_colours)) {
        if (col_name %in% names(existing_colours_list)) {
            next
        }

        existing_colours_list[[col_name]] <- predefined_colours[[col_name]]
    }

    return(existing_colours_list)
}

starlng_write_app_file <- function(file_path, title_name = "", height_ratio = 0.7) {
    title_name <- paste0("Starlng ShinyApp - ", title_name)
    content <- paste0("library(Starlng)

ui <- shiny::fluidPage(
    shinyjs::useShinyjs(),
    ui_global_setttings(),
    shiny::navbarPage(
        title = \"", title_name, "\",
        windowTitle = \"", title_name, "\",
        id = \"nav\",
        shiny::tabPanel(
            title = \"Metadata\",
            value = \"metadata\",
            ui_metadata_umap(\"metadata_umap\")
        ),
        shiny::tabPanel(
            title = \"Gene Table\",
            value = \"gene_info_table\",
            ui_gene_info_table(\"gene_info_table\")
        ),
        shiny::tabPanel(
            title = \"Gene UMAP\",
            value = \"gene_umap\",
            ui_gene_umap(\"gene_umap\")
        ),
        shiny::tabPanel(
            title = \"Gene Clustering\",
            value = \"gene_clustering\",
            ui_gene_clustering(\"gene_clustering\")
        ),
        shiny::tabPanel(
            title = \"Module UMAP\",
            value = \"module_umap\",
            ui_module_umap(\"module_umap\")
        ),
        shiny::tabPanel(
            title = \"Module Enrichment\",
            value = \"module_enrichment\",
            ui_module_enrichment(\"module_enrichment\")
        ),
        shiny::tabPanel(
            title = \"Pseudotime Ordering\",
            value = \"pseudotime\",
            ui_pseudotime_select_cells(\"pseudotime_select_cells\")
        ),
        shiny::tabPanel(
            title = \"Gene Heatmaps\",
            value = \"gene_heatmaps\",
            ui_module_metadata_heatmap(\"module_metadata_heatmap\")
        ),
        position = \"fixed-top\",
        inverse = TRUE
    )
)

server <- function(input, output, session) {
    prepare_session(shiny::reactive(input$dimension), ", height_ratio, ")
    update_gears_width(session)
    update_tabs(session)
    server_metadata_umap(\"metadata_umap\")
    filtered_genes <- server_gene_info_table(\"gene_info_table\")
    server_gene_umap(\"gene_umap\")
    server_gene_clustering(\"gene_clustering\", filtered_genes)
    server_module_umap(\"module_umap\")
    server_module_enrichment(\"module_enrichment\")
    server_pseudotime_select_cells(\"pseudotime_select_cells\")
    server_module_metadata_heatmap(\"module_metadata_heatmap\")
}

shiny::shinyApp(ui = ui, server = server)

    ")
    write(content, file = file_path)
}

#' @export
starlng_write_app_monocle <- function(folder_path,
                              monocle_object,
                              app_title_name = "",
                              learn_graph_parameters = list(),
                              gene_filtering_function = function(info_gene_df) {
                                  rownames(info_gene_df %>% dplyr::filter(.data$morans_I > 0.1))
                              },
                              clustering_parameters = list(
                                  "n_neighbours" = seq(from = 5, to = 50, by = 5),
                                  "graph_type" = "snn",
                                  "prune_value" = -1,
                                  "resolutions" = seq(from = 0.1, to = 1, by = 0.1),
                                  "quality_functions" = c("RBConfigurationVertexPartition"),
                                  "number_iterations" = 5,
                                  "number_repetitions" = 100
                              ),
                              ecc_threshold = 0.9,
                              freq_threshold = 30,
                              save_entire_monocle = TRUE,
                              discrete_colours = list(),
                              continuous_colours = list(),
                              max_n_colors = 40,
                              verbose = FALSE,
                              compression_level = 7,
                              chunk_size = 100,
                              nthreads = 1) {
    if (!dir.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)
    }
    if (verbose) print("Writing app.R...")
    starlng_write_app_file(file.path(folder_path, "app.R"), app_title_name)
    folder_path <- file.path(folder_path, "objects")
    if (!dir.exists(folder_path)) {
        dir.create(folder_path, showWarnings = FALSE)
    }

    # entire monocle object
    monocle_object <- do.call(
        custom_learn_graph,
        c(list("mon_obj" = monocle_object), learn_graph_parameters)
    )

    mon_path <- file.path(folder_path, "monocle_object.qs")
    if (save_entire_monocle) {
        if (verbose) print("Writing the monocle object...")
        qs::qsave(monocle_object, file = mon_path, nthreads = nthreads)
    }

    # autocorrelation object
    if (verbose) print("Writing the gene info table...")
    info_gene_path <- file.path(folder_path, "genes_info.csv")
    autocorr_result <- monocle3::graph_test(
        monocle_object,
        neighbor_graph = "principal_graph",
        cores = nthreads,
        verbose = verbose
    )
    rownames(autocorr_result) <- autocorr_result$gene_short_name
    autocorr_result <- autocorr_result[, c("morans_test_statistic", "morans_I", "q_value")] %>%
        dplyr::arrange(dplyr::desc(.data$morans_I))
    write.csv(autocorr_result, file = info_gene_path)

    # trajectory ggplot object
    if (verbose) print("Writing the trajectory ggplot object...")
    trajectory_path <- file.path(folder_path, "trajectory_ggplot.qs")
    trajectory_ggplot <- monocle3::plot_cells(
        monocle_object,
        cell_size = 0,
        label_roots = FALSE,
        label_leaves = FALSE,
        label_branch_point = FALSE,
        label_cell_groups = FALSE,
        label_principal_points = FALSE
    )

    trajectory_ggplot$layers <- trajectory_ggplot$layers[length(trajectory_ggplot$layers)]
    for (other_names in names(trajectory_ggplot)) {
        if (other_names == "layers") {
            next
        }

        trajectory_ggplot[[other_names]] <- NULL
    }
    qs::qsave(trajectory_ggplot, file = trajectory_path)

    # metadata object
    if (verbose) print("Writing the metadata object...")
    metadata_path <- file.path(folder_path, "metadata.qs")
    mtd_df <- as.data.frame(monocle_object@colData)
    mtd_df <- cbind(mtd_df, monocle_object@int_colData$reducedDims$UMAP)
    qs::qsave(mtd_df, file = metadata_path)
    
    chosen_slot <- "normalized_data"
    if (!(chosen_slot %in% names(monocle_object@assays@data))) {
        chosen_slot <- "counts"
    }
    expr_matrix <- monocle_object@assays@data[[chosen_slot]]
    # diet monocle object
    if (verbose) print("Writing the diet monocle object...")
    mon_path <- file.path(folder_path, "digest_monocle_object.qs")
    monocle_object <- diet_monocle_object(monocle_object)
    gc()
    qs::qsave(monocle_object, file = mon_path, nthreads = nthreads)

    # gene clustering
    skip_clustering <- length(clustering_parameters) == 0
    if (!skip_clustering) {
        if (verbose) print("Applying clustering and stability assessment...")
        assess_path <- file.path(folder_path, "full_stability_assessment.qs")
        clusters_path <- file.path(folder_path, "stable_modules.csv")

        clustering_parameters[["merge_identical_partitions"]] <- TRUE
        chosen_genes <- gene_filtering_function(autocorr_result)
        RhpcBLASctl::blas_set_num_threads(nthreads)
        clustering_parameters[["embedding"]] <- get_feature_loading(
            expr_matrix = expr_matrix[chosen_genes, ],
            30
        )
        RhpcBLASctl::blas_set_num_threads(1)

        created_cluster <- FALSE
        if (nthreads > 1 && foreach::getDoParWorkers() == 1) {
            created_cluster <- TRUE
            par_cluster <- parallel::makeCluster(nthreads)

            doParallel::registerDoParallel(par_cluster)
        }

        clust_results <- do.call(clustering_pipeline, clustering_parameters)
        if (created_cluster) {
            parallel::stopCluster(par_cluster)
            foreach::registerDoSEQ()
        }
        if (verbose) print("Writing clustering results...")

        qs::qsave(clust_results, file = assess_path, nthreads = nthreads)
        best_config <- select_best_configuration(clust_results)
        clust_results <- clust_results[[best_config[[1]]]][[best_config[[2]]]]

        stb_clust <- ClustAssess::choose_stable_clusters(
            clust_results$k,
            ecc_threshold = ecc_threshold,
            freq_threshold = freq_threshold
        )
        stb_df <- data.frame(genes = chosen_genes)
        for (k in names(stb_clust)) {
            stb_df[[paste0("stable_modules_", k)]] <- factor(stb_clust[[k]]$partitions[[1]]$mb)
        }

        write.csv(stb_df, file = clusters_path, row.names = FALSE, quote = FALSE)
    }

    # expression matrix
    if (verbose) print("Writing the expression matrix...")
    expr_path <- file.path(folder_path, "expression.h5")
    write_gene_matrix_dense_h5(
        expr_matrix,
        expr_path,
        compression_level = compression_level,
        chunk_size = chunk_size,
        all_at_once = TRUE
    )

    # colours object
    if (verbose) print("Writing the colours object...")
    colour_path <- file.path(folder_path, "colours.qs")
    colour_list <- list(
        "discrete" = generate_discrete_colours(
            mtd_df,
            discrete_colours,
            max_n_colors
        ),
        "continuous" = generate_continuous_colours(continuous_colours)
    )
    qs::qsave(colour_list, file = colour_path)
}


#' @export
starlng_write_app_default <- function(folder_path,
                              expression_matrix,
                              metadata_df = NULL,
                              pca_embedding = NULL,
                              umap_embedding = NULL,
                              app_title_name = "",
                              learn_graph_parameters = list(),
                              gene_filtering_function = function(info_gene_df) {
                                  rownames(info_gene_df %>% dplyr::filter(.data$morans_I > 0.1))
                              },
                              clustering_parameters = list(
                                  "n_neighbours" = seq(from = 5, to = 50, by = 5),
                                  "graph_type" = "snn",
                                  "prune_value" = -1,
                                  "resolutions" = seq(from = 0.1, to = 1, by = 0.1),
                                  "quality_functions" = c("RBConfigurationVertexPartition"),
                                  "number_iterations" = 5,
                                  "number_repetitions" = 100
                              ),
                              ecc_threshold = 0.9,
                              freq_threshold = 30,
                              save_entire_monocle = TRUE,
                              discrete_colours = list(),
                              continuous_colours = list(),
                              max_n_colors = 40,
                              verbose = FALSE,
                              compression_level = 7,
                              chunk_size = 100,
                              nthreads = 1) {
    monocle_object <- ClustAssess::create_monocle_default(
        normalized_expression_matrix = expression_matrix,
        metadata_df = metadata_df,
        pca_embedding = pca_embedding,
        umap_embedding = umap_embedding
    )

    starlng_write_app_monocle(
        folder_path = folder_path,
        monocle_object = monocle_object,
        app_title_name = app_title_name,
        learn_graph_parameters = learn_graph_parameters,
        gene_filtering_function = gene_filtering_function,
        clustering_parameters = clustering_parameters,
        ecc_threshold = ecc_threshold,
        freq_threshold = freq_threshold,
        save_entire_monocle = save_entire_monocle,
        discrete_colours = discrete_colours,
        continuous_colours = continuous_colours,
        max_n_colors = max_n_colors,
        verbose = verbose,
        compression_level = compression_level,
        chunk_size = chunk_size,
        nthreads = nthreads
    )
}

#' @export
starlng_write_app_clustassess <- function(folder_path,
                                  expression_matrix,
                                  clustassess_object,
                                  metadata_df,
                                  stable_feature_type,
                                  stable_feature_set_size,
                                  stable_clustering_method,
                                  stable_n_clusters = NULL,
                                  use_all_genes = TRUE,
                                  app_title_name = "",
                                  learn_graph_parameters = list(),
                                  gene_filtering_function = function(info_gene_df) {
                                      rownames(info_gene_df %>% dplyr::filter(.data$morans_I > 0.1))
                                  },
                                  clustering_parameters = list(
                                      "n_neighbours" = seq(from = 5, to = 50, by = 5),
                                      "graph_type" = "snn",
                                      "prune_value" = -1,
                                      "resolutions" = seq(from = 0.1, to = 1, by = 0.1),
                                      "quality_functions" = c("RBConfigurationVertexPartition"),
                                      "number_iterations" = 5,
                                      "number_repetitions" = 100
                                  ),
                                  ecc_threshold = 0.9,
                                  freq_threshold = 30,
                                  save_entire_monocle = TRUE,
                                  discrete_colours = list(),
                                  continuous_colours = list(),
                                  max_n_colors = 40,
                                  verbose = FALSE,
                                  compression_level = 7,
                                  chunk_size = 100,
                                  nthreads = 1) {
    monocle_object <- ClustAssess::create_monocle_from_clustassess(
        normalized_expression_matrix = expression_matrix,
        metadata_df = metadata_df,
        clustassess_object = clustassess_object,
        stable_feature_type = stable_feature_type,
        stable_feature_set_size = stable_feature_set_size,
        stable_clustering_method = stable_clustering_method,
        stable_n_clusters = stable_n_clusters,
        use_all_genes = use_all_genes
    )

    starlng_write_app_monocle(
        folder_path = folder_path,
        monocle_object = monocle_object,
        app_title_name = app_title_name,
        learn_graph_parameters = learn_graph_parameters,
        gene_filtering_function = gene_filtering_function,
        clustering_parameters = clustering_parameters,
        ecc_threshold = ecc_threshold,
        freq_threshold = freq_threshold,
        save_entire_monocle = save_entire_monocle,
        discrete_colours = discrete_colours,
        continuous_colours = continuous_colours,
        max_n_colors = max_n_colors,
        verbose = verbose,
        compression_level = compression_level,
        chunk_size = chunk_size,
        nthreads = nthreads
    )
}

#' @export
starlng_write_app_clustassess_app <- function(folder_path,
                                      ca_app_folder,
                                      stable_feature_type,
                                      stable_feature_set_size,
                                      stable_clustering_method,
                                      stable_n_clusters = NULL,
                                      use_all_genes = TRUE,
                                      app_title_name = "",
                                      learn_graph_parameters = list(),
                                      gene_filtering_function = function(info_gene_df) {
                                          rownames(info_gene_df %>% dplyr::filter(.data$morans_I > 0.1))
                                      },
                                      clustering_parameters = list(
                                          "n_neighbours" = seq(from = 5, to = 50, by = 5),
                                          "graph_type" = "snn",
                                          "prune_value" = -1,
                                          "resolutions" = seq(from = 0.1, to = 1, by = 0.1),
                                          "quality_functions" = c("RBConfigurationVertexPartition"),
                                          "number_iterations" = 5,
                                          "number_repetitions" = 100
                                      ),
                                      ecc_threshold = 0.9,
                                      freq_threshold = 30,
                                      save_entire_monocle = TRUE,
                                      discrete_colours = list(),
                                      continuous_colours = list(),
                                      max_n_colors = 40,
                                      verbose = FALSE,
                                      compression_level = 7,
                                      chunk_size = 100,
                                      nthreads = 1) {
    monocle_object <- ClustAssess::create_monocle_from_clustassess_app(
        app_folder = ca_app_folder,
        stable_feature_type = stable_feature_type,
        stable_feature_set_size = stable_feature_set_size,
        stable_clustering_method = stable_clustering_method,
        stable_n_clusters = stable_n_clusters,
        use_all_genes = use_all_genes
    )

    starlng_write_app_monocle(
        folder_path = folder_path,
        monocle_object = monocle_object,
        app_title_name = app_title_name,
        learn_graph_parameters = learn_graph_parameters,
        gene_filtering_function = gene_filtering_function,
        clustering_parameters = clustering_parameters,
        ecc_threshold = ecc_threshold,
        freq_threshold = freq_threshold,
        save_entire_monocle = save_entire_monocle,
        discrete_colours = discrete_colours,
        continuous_colours = continuous_colours,
        max_n_colors = max_n_colors,
        verbose = verbose,
        compression_level = compression_level,
        chunk_size = chunk_size,
        nthreads = nthreads
    )
}
