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

#' Create Starlng Shiny app from Monocle object
#' 
#' @description Runs the entire Starlng pipeline on a Monocle object and
#' generates a folder with all necessary files to run a Shiny app.
#' 
#' @param folder_path The path of the folder where the app will be saved.
#' @param monocle_object The monocle3 object that will be used as a starting
#' point.
#' @param app_title_name The title of the shiny app.
#' @param learn_graph_parameters A list of parameters that will be passed to
#' the `custom_learn_graph` function. Possible options:
#' - `nodes_per_log10_cells` - the number of nodes per log10 cells. It is used
#' to influence the level of detail in the inferred trajectory.
#' - `learn_graph_controls` - a list of parameters that will be passed to the
#' `learn_graph` function of the monocle3 package. These parameters are usually
#' used to influence the convergence of the algorithm or the pruning of the
#' resulting trajectory.
#' @param gene_filtering_function A function that will be used to filter the
#' genes that will be used in the clustering analysis. The function should take
#' only one parameter, which is a data frame. The data frame will contain the
#' following columns:
#' - `morans_I` - the Moran's I value of the gene.
#' - `q_value` - the adjusted p-value of the gene being spatially autocorrelated.
#' - `average_expression` - the average expression of the gene.
#' - `n_expressed_cells` - the number of cells that express the gene.
#' - `average_expression_nonzero` - the average expression of the gene in cells
#' that express it.
#' - `percent_expressed_cells` - the percentage of cells that express the gene.
#' The function should return a list of gene names that will be considered to be
#' the result of the filtering.
#' @param clustering_parameters A list of parameters that will be passed to the
#' `clustering_pipeline` function. Consult the documentation of the function for
#' further details.
#' @param ecc_threshold The threshold applied on the Element-Centric Consistency
#' (ECC) score to determine if a number of cluster is stable or not. By default,
#' it is set to 0.9.
#' @param freq_threshold The threshold applied on the number of times a number
#' of clusters is obtained after running the `clustering_pipeline` function. By
#' default, it is set to 30.
#' @param save_entire_monocle If TRUE, saves the monocle object in the app
#' folder. This object could be used to perform additional changes to the app
#' or to continue the downstream analysis.
#' @param discrete_colours A list of colours that will be used to colour the
#' groups from the discrete metadata. The list should consist in an
#' association between the number of groups and a colour palette. For the
#' missing number of groups, the `qualpalr` package will be used. Defaults
#' to an empty list.
#' @param continuous_colours A list of colours that will be used to colour the
#' continuous metadata and the gene expression. The list should consist in an
#' association between the name of the colour palette and a list of colours.
#' To this list the following palettes will be automatically added: viridis,
#' plasam, white-red and blue-red. Defaults to an empty list.
#' @param max_n_colors The maximum number of colours used for discrete groups.
#' Defaults to 40.
#' @param verbose If TRUE, prints intermediate progress messages of the
#' pipeline. Defaults to FALSE.
#' @param compression_level The compression level used to save the expression
#' matrix in the HDF5 format. Defaults to 7.
#' @param chunk_size The chunk size used to store the expression matrix in the
#' HDF5 format. It is recommended to use values between 10 and 1000, as they
#' might improve the speed of reading the rows from the file. Defaults to 100.
#' @param nthreads The number of threads used to perform the calculations. This
#' parameter will be automatically passed to the `graph_test` function that
#' will calculate the Moran's I test. For the clustering pipeline, if no
#' parallel backend is already registered, the function will generate
#' a PSOCK cluster with the number of threads specified. Defaults to 1.
#'
#' @return The function does not return anything, but it creates a folder which
#' contains all necessary files to run the Starlng Shiny app.
#' - `app.R` - the main file of the Shiny app.
#' - `objects` - a folder that contains the following files:
#'  - `monocle_object.qs` - the monocle object that was used to generate the app.
#'  - `recommended_pseudotime.qs` - the recommended pseudotime values.
#'  - `genes_info.csv` - a table with the gene information.
#'  - `trajectory_ggplot.qs` - the ggplot object of the trajectory.
#'  - `metadata.qs` - the metadata object.
#'  - `diet_monocle_object.qs` - the monocle object with the diet information.
#'  - `expression.h5` - the expression matrix in the HDF5 format.
#'  - `colours.qs` - the colours object.
#'
#' @export
starlng_write_app_monocle <- function(folder_path,
                              monocle_object,
                              app_title_name = "",
                              learn_graph_parameters = list(
                                  nodes_per_log10_cells = 30,
                                  learn_graph_controls = list(
                                    eps = 1e-5,
                                    maxiter = 100
                                  )
                              ),
                              gene_filtering_function = function(info_gene_df) {
                                  rownames(info_gene_df %>% dplyr::filter(.data$morans_I > 0.1, .data$q_value < 0.05))
                              },
                              clustering_parameters = list(
                                  "n_neighbours" = seq(from = 5, to = 50, by = 5),
                                  "graph_type" = "snn",
                                  "prune_value" = -1,
                                  "resolutions" = list(
                                      "RBConfigurationVertexPartition" = seq(from = 0.1, to = 2, by = 0.1),
                                      "RBERVertexPartition" = NULL,
                                      "ModularityVertexPartition" = NULL
                                  ),
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

    # recommended pseudotime
    if (verbose) print("Writing the recommended pseudotime...")
    psd_path <- file.path(folder_path, "recommended_pseudotime.qs")
    psd_values <- get_pseudotime_recommendation(monocle_object)
    qs::qsave(psd_values, file = psd_path)

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
    monocle_object <- monocle_object[rownames(autocorr_result), ]

    chosen_slot <- "normalized_data"
    if (!(chosen_slot %in% names(monocle_object@assays@data))) {
        chosen_slot <- "counts"
    }
    expr_matrix <- monocle_object@assays@data[[chosen_slot]]

    sum_function <- base::rowSums
    if (inherits(expr_matrix, "dgCMatrix")) {
        sum_function <- Matrix::rowSums
    }
    autocorr_result$average_expression <- sum_function(expr_matrix)
    autocorr_result$n_expressed_cells <- sum_function(expr_matrix > 0)
    autocorr_result$average_expression_nonzero <- autocorr_result$average_expression / autocorr_result$n_expressed_cells
    autocorr_result$average_expression <- autocorr_result$average_expression / ncol(expr_matrix)
    autocorr_result$percent_expressed_cells <- autocorr_result$n_expressed_cells / ncol(expr_matrix)
    utils::write.csv(autocorr_result, file = info_gene_path)

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
    mtd_df$pseudotime <- psd_values$recommended_pseudotime[rownames(mtd_df)]
    mtd_df <- cbind(mtd_df, monocle_object@int_colData$reducedDims$UMAP)
    qs::qsave(mtd_df, file = metadata_path)
    
    # diet monocle object
    if (verbose) print("Writing the diet monocle object...")
    mon_path <- file.path(folder_path, "diet_monocle_object.qs")
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
            npcs = 30,
            approx = FALSE
        )
        RhpcBLASctl::blas_set_num_threads(1)

        created_cluster <- FALSE
        if (nthreads > 1 && foreach::getDoParWorkers() == 1) {
            created_cluster <- TRUE
            par_cluster <- parallel::makePSOCKcluster(nthreads)

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
        message("Best configuration: n_neigh - ", best_config[[1]], ", quality function - ", best_config[[2]])
        qs::qsave(best_config, file = file.path(folder_path, "best_configuration.qs"))
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

        utils::write.csv(stb_df, file = clusters_path, row.names = FALSE, quote = FALSE)
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

#' Create Starlng Shiny app from a normalized gene expression matrix
#'
#' @description Runs the entire Starlng pipeline on a normalized gene expression
#' matrix and generates a folder with all necessary files to run a Shiny app.
#'
#' @note For more details about the other parameters and the content of the
#' shiny folder, consult the documentation of the `starlng_write_app_monocle`
#' function.
#'
#' @param expression_matrix A normalized gene by cell expression matrix. The
#' matrix should have the rownames and the colnames defined.
#' @param metadata_df A data frame with cell metadata specific to the experiment.
#' If NULL, the function will create a one-group metadata with the name
#' "one_level".
#' @param pca_embedding A matrix with the PCA embedding of the cells. If NULL,
#' the function will calculate the PCA embedding while creating the monocle3
#' object.
#' @param umap_embedding A matrix with the UMAP embedding of the cells. If NULL,
#' the function will calculate the UMAP embedding while creating the monocle3.
#' @param folder_path Check `starlng_write_app_monocle` documentation.
#' @param app_title_name Check `starlng_write_app_monocle` documentation.
#' @param learn_graph_parameters Check `starlng_write_app_monocle` documentation.
#' @param gene_filtering_function Check `starlng_write_app_monocle` documentation.
#' @param clustering_parameters Check `starlng_write_app_monocle` documentation.
#' @param ecc_threshold Check `starlng_write_app_monocle` documentation.
#' @param freq_threshold Check `starlng_write_app_monocle` documentation.
#' @param save_entire_monocle Check `starlng_write_app_monocle` documentation.
#' @param discrete_colours Check `starlng_write_app_monocle` documentation.
#' @param continuous_colours Check `starlng_write_app_monocle` documentation.
#' @param max_n_colors Check `starlng_write_app_monocle` documentation.
#' @param verbose Check `starlng_write_app_monocle` documentation.
#' @param compression_level Check `starlng_write_app_monocle` documentation.
#' @param chunk_size Check `starlng_write_app_monocle` documentation.
#' @param nthreads Check `starlng_write_app_monocle` documentation.
#'
#' @export
starlng_write_app_default <- function(folder_path,
                              expression_matrix,
                              metadata_df = NULL,
                              pca_embedding = NULL,
                              umap_embedding = NULL,
                              app_title_name = "",
                              learn_graph_parameters = list(
                                  nodes_per_log10_cells = 30,
                                  learn_graph_controls = list(
                                    eps = 1e-5,
                                    maxiter = 100
                                  )
                              ),
                              gene_filtering_function = function(info_gene_df) {
                                  rownames(info_gene_df %>% dplyr::filter(.data$morans_I > 0.1, .data$q_value < 0.05))
                              },
                              clustering_parameters = list(
                                  "n_neighbours" = seq(from = 5, to = 50, by = 5),
                                  "graph_type" = "snn",
                                  "prune_value" = -1,
                                  "resolutions" = list(
                                      "RBConfigurationVertexPartition" = seq(from = 0.1, to = 2, by = 0.1),
                                      "RBERVertexPartition" = NULL,
                                      "ModularityVertexPartition" = NULL
                                  ),
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

#' Create Starlng Shiny app from a ClustAssess object
#'
#' @description Runs the entire Starlng pipeline on a normalized gene expression
#' matrix togheter with the output of the ClustAssess automatic pipeline 
#' and generates a folder with all necessary files to run a Shiny app.
#'
#' @note For more details about the other parameters and the content of the
#' shiny folder, consult the documentation of the `starlng_write_app_monocle`
#' function.
#'
#' @param expression_matrix A normalized gene by cell expression matrix. The
#' matrix should have the rownames and the colnames defined.
#' @param metadata_df A data frame with cell metadata specific to the experiment.
#' If NULL, the function will create a one-group metadata with the name
#' "one_level".
#' @param clustassess_object The output of the automatic pipeline of the
#' ClustAssess package.
#' @param stable_feature_type The feature type that is chosen from the
#' ClustAssess analysis.
#' @param stable_feature_set_size The size of the feature set that is chosen
#' from the ClustAssess analysis.
#' @param stable_clustering_method The clustering method that is chosen from the
#' ClustAssess analysis.
#' @param stable_n_clusters The number of clusters that is chosen from the
#' ClustAssess analysis. If NULL, the function will use all the clusters that
#' are obtained. Defaults to NULL.
#' @param use_all_genes If TRUE, uses all the genes in the expression matrix.
#' If FALSE, uses only the stable feature set as determined by ClustAssess.
#' Defaults to TRUE.
#' @param folder_path Check `starlng_write_app_monocle` documentation.
#' @param app_title_name Check `starlng_write_app_monocle` documentation.
#' @param learn_graph_parameters Check `starlng_write_app_monocle` documentation.
#' @param gene_filtering_function Check `starlng_write_app_monocle` documentation.
#' @param clustering_parameters Check `starlng_write_app_monocle` documentation.
#' @param ecc_threshold Check `starlng_write_app_monocle` documentation.
#' @param freq_threshold Check `starlng_write_app_monocle` documentation.
#' @param save_entire_monocle Check `starlng_write_app_monocle` documentation.
#' @param discrete_colours Check `starlng_write_app_monocle` documentation.
#' @param continuous_colours Check `starlng_write_app_monocle` documentation.
#' @param max_n_colors Check `starlng_write_app_monocle` documentation.
#' @param verbose Check `starlng_write_app_monocle` documentation.
#' @param compression_level Check `starlng_write_app_monocle` documentation.
#' @param chunk_size Check `starlng_write_app_monocle` documentation.
#' @param nthreads Check `starlng_write_app_monocle` documentation.
#'
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
                                  learn_graph_parameters = list(
                                        nodes_per_log10_cells = 30,
                                        learn_graph_controls = list(
                                            eps = 1e-5,
                                            maxiter = 100
                                        )
                                    ),
                                  gene_filtering_function = function(info_gene_df) {
                                      rownames(info_gene_df %>% dplyr::filter(.data$morans_I > 0.1, .data$q_value < 0.05))
                                  },
                                  clustering_parameters = list(
                                      "n_neighbours" = seq(from = 5, to = 50, by = 5),
                                      "graph_type" = "snn",
                                      "prune_value" = -1,
                                      "resolutions" = list(
                                          "RBConfigurationVertexPartition" = seq(from = 0.1, to = 2, by = 0.1),
                                          "RBERVertexPartition" = NULL,
                                          "ModularityVertexPartition" = NULL
                                      ),
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

#' Create Starlng Shiny app from a ClustAssess object
#'
#' @description Runs the entire Starlng pipeline on a normalized gene expression
#' matrix togheter with the content of the ClustAssess Shiny app
#' and generates a folder with all necessary files to run a Shiny app.
#'
#' @note For more details about the other parameters and the content of the
#' shiny folder, consult the documentation of the `starlng_write_app_monocle`
#' function.
#'
#' @param ca_app_folder The path of the folder containing the files associated
#' with the ClustAssess Shiny app.
#' @param stable_feature_type The feature type that is chosen from the
#' ClustAssess analysis.
#' @param stable_feature_set_size The size of the feature set that is chosen
#' from the ClustAssess analysis.
#' @param stable_clustering_method The clustering method that is chosen from the
#' ClustAssess analysis.
#' @param stable_n_clusters The number of clusters that is chosen from the
#' ClustAssess analysis. If NULL, the function will use all the clusters that
#' are obtained. Defaults to NULL.
#' @param use_all_genes If TRUE, uses all the genes in the expression matrix.
#' If FALSE, uses only the stable feature set as determined by ClustAssess.
#' Defaults to TRUE.
#' @param folder_path Check `starlng_write_app_monocle` documentation.
#' @param app_title_name Check `starlng_write_app_monocle` documentation.`
#' @param learn_graph_parameters Check `starlng_write_app_monocle` documentation.`
#' @param gene_filtering_function Check `starlng_write_app_monocle` documentation.`
#' @param clustering_parameters Check `starlng_write_app_monocle` documentation.`
#' @param ecc_threshold Check `starlng_write_app_monocle` documentation.`
#' @param freq_threshold Check `starlng_write_app_monocle` documentation.`
#' @param save_entire_monocle Check `starlng_write_app_monocle` documentation.`
#' @param discrete_colours Check `starlng_write_app_monocle` documentation.`
#' @param continuous_colours Check `starlng_write_app_monocle` documentation.`
#' @param max_n_colors Check `starlng_write_app_monocle` documentation.`
#' @param verbose Check `starlng_write_app_monocle` documentation.
#' @param compression_level Check `starlng_write_app_monocle` documentation.`
#' @param chunk_size Check `starlng_write_app_monocle` documentation.
#' @param nthreads Check `starlng_write_app_monocle` documentation.
#'
#' @export
starlng_write_app_clustassess_app <- function(folder_path,
                                      ca_app_folder,
                                      stable_feature_type,
                                      stable_feature_set_size,
                                      stable_clustering_method,
                                      stable_n_clusters = NULL,
                                      use_all_genes = TRUE,
                                      app_title_name = "",
                                       learn_graph_parameters = list(
                                        nodes_per_log10_cells = 30,
                                        learn_graph_controls = list(
                                            eps = 1e-5,
                                            maxiter = 100
                                        )
                                    ),
                                      gene_filtering_function = function(info_gene_df) {
                                          rownames(info_gene_df %>% dplyr::filter(.data$morans_I > 0.1, .data$q_value < 0.05))
                                      },
                                      clustering_parameters = list(
                                          "n_neighbours" = seq(from = 5, to = 50, by = 5),
                                          "graph_type" = "snn",
                                          "prune_value" = -1,
                                          "resolutions" = list(
                                              "RBConfigurationVertexPartition" = seq(from = 0.1, to = 2, by = 0.1),
                                              "RBERVertexPartition" = NULL,
                                              "ModularityVertexPartition" = NULL
                                          ),
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
