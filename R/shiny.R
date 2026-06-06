#' @importFrom dplyr %>% .data

generate_discrete_colours <- function(metadata_df,
                                      existing_colours_list = NULL,
                                      max_n_colors = 40,
                                      qualpal_pallete = list(h = c(0, 360), s = c(0.2, 0.5), l = c(0.6, 0.85))) {
    if (is.null(existing_colours_list)) {
        existing_colours_list <- list()
    }

    for (mtd_names in colnames(metadata_df)) {
        if (!(inherits(metadata_df[[mtd_names]], c("factor", "character")))) {
            next
        }

        n_colors <- length(unique(metadata_df[[mtd_names]]))
        if (n_colors > 200) {
            warning("The number of colors for the metadata column ", mtd_names, " is higher than 200. Consider providing a custom color palette for this column or reducing the number of groups.")
            next
        }
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

starlng_write_app_file <- function(file_path, title_name = "", height_ratio = 0.7, enrichment_organism = "hsapiens") {
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
            title = \"Pseudotime Ordering\",
            value = \"pseudotime\",
            ui_pseudotime_select_cells(\"pseudotime_select_cells\")
        ),
        shiny::tabPanel(
            title = \"Gene Heatmaps\",
            value = \"gene_heatmaps\",
            ui_module_metadata_heatmap(\"module_metadata_heatmap\")
        ),
        shiny::tabPanel(
            title = \"Functional Assessment\",
            value = \"functional_assessment\",
            ui_functional_assessment(\"functional_assessment\")
        ),
        position = \"fixed-top\",
        inverse = TRUE
    )
)

server <- function(input, output, session) {
    prepare_session(session, shiny::reactive(input$dimension), ", height_ratio, ", '", enrichment_organism, "')
    update_gears_width(session)
    update_tabs(session)
    server_metadata_umap(\"metadata_umap\")
    # filtered_genes <- server_gene_info_table(\"gene_info_table\")
    filtered_genes <- server_gene_umap(\"gene_umap\")
    server_gene_clustering(\"gene_clustering\", filtered_genes)
    server_module_umap(\"module_umap\")
    server_pseudotime_select_cells(\"pseudotime_select_cells\")
    server_module_metadata_heatmap(\"module_metadata_heatmap\")
    server_functional_assessment(\"functional_assessment\")
}

shiny::shinyApp(ui = ui, server = server)

    ")
    write(content, file = file_path)
}

starlng_write_module_summaries <- function(stb_clust,
                                           folder_path,
                                           expr_matrix,
                                           chosen_genes,
                                           gene_adj_matrix,
                                           trajectory_object,
                                           cell_umap_df,
                                           scale_threshold,
                                           enrichment_organism = "hsapiens",
                                           verbose = FALSE) {
    clusters_path <- file.path(folder_path, "module_summaries.h5")
    if (file.exists(clusters_path)) {
        warning("Overwriting existing file. Waiting for 5 seconds...")
        Sys.sleep(5.5)
        file.remove(clusters_path)
    }
    rhdf5::h5createFile(clusters_path)
    rhdf5::h5write(names(stb_clust), clusters_path, "all_modules")
    rhdf5::h5write(chosen_genes, clusters_path, "genes")
    rhdf5::h5write(scale_threshold, clusters_path, "scale_threshold")
    rhdf5::h5write("average", clusters_path, "summary_method")

    for (k in names(stb_clust)) {
        rhdf5::h5createGroup(clusters_path, k)

        # partitions 
        module_partition <- stb_clust[[k]]$partitions[[1]]$mb
        names(module_partition) <- chosen_genes
        module_partition <- module_partition
        rhdf5::h5write(module_partition, clusters_path, paste0(k, "/clustering"))

        # module summary
        verbose_print(paste("Calculating module summaries for", k, "modules..."), verbose)
        module_names <- as.character(seq_len(as.integer(k)))
        module_summaries_list <- do.call(cbind, lapply(module_names, function(module) {
            module_genes <- names(module_partition)[module_partition == module]
            module_expr <- expr_matrix[module_genes, , drop = FALSE]
            module_summary <- voting_scheme(
                expression_matrix = module_expr,
                genes = module_genes,
                thresh_percentile = 0,
                thresh_value = 0,
                n_coexpressed_thresh = 1,
                summary_function = mean
            )
            return(module_summary)
        }))
        colnames(module_summaries_list) <- module_names
        rhdf5::h5createDataset(
            file = clusters_path,
            dataset = paste0(k, "/expression_summaries"),
            dims = dim(module_summaries_list),
            storage.mode = "double"
        )

        rhdf5::h5write(module_summaries_list, file = clusters_path, name = paste0(k, "/expression_summaries"))

      
        # gene hubs
        genes_per_module <- split(chosen_genes, module_partition)
        genes_per_module <- genes_per_module[module_names]

        ## convert to list
        verbose_print(paste("Calculating gene hub statistics for", k, "modules..."), verbose)
        module_summaries_list <- split(module_summaries_list, col(module_summaries_list))
        module_summaries_list <- module_summaries_list[module_names]
        module_weights <- get_per_module_weight(gene_adj_matrix, genes_per_module)
        rank_df <- get_gene_overlap_stat(
            expr_matrix = expr_matrix,
            gene_modules = genes_per_module,
            module_expr = module_summaries_list,
            module_expression_threshold = scale_threshold,
            gene_expression_threshold = 0,
            gene_expression_percentile = 0.5,
            scale = TRUE,
            total_weight = module_weights
        )
        rank_df$module <- module_partition[rownames(rank_df)]
        rank_df <- rank_df[chosen_genes, ]

        rhdf5::h5write(rank_df, file = clusters_path, name = paste0(k, "/gene_hub_stats"))

        # module adjacency
        verbose_print(paste("Infer the module adjacency for", k, "modules..."), verbose)
        closest_module <- get_closest_node_to_module(
            trajectory_object = trajectory_object,
            cell_umap = cell_umap_df,
            module_expr = module_summaries_list,
            expression_threshold = scale_threshold,
            scale = TRUE)
        closest_module <- closest_module[module_names]
        rhdf5::h5write(closest_module, file = clusters_path, name = paste0(k, "/closest_nodes_to_module"))

        module_adjacency <- get_module_transitions(trajectory_object, closest_module, start_node = NULL)
        rhdf5::h5write(module_adjacency, file = clusters_path, name = paste0(k, "/module_adjacency"))

        # tfs
        verbose_print(paste("Performing identification of transcription factors for", k, "modules..."), verbose)
        tfs <- get_transcription_factors(
            gene_list = genes_per_module,
            organism = enrichment_organism,
            p_value_threshold = 0.05,
            correction_method = "fdr",
            bg_genes = rownames(expr_matrix)
        )
        tf_stats <- get_tf_stats(tfs, include_intersection_set = TRUE)

        rhdf5::h5write(tf_stats, file = clusters_path, name = paste0(k, "/tfs"))

        rhdf5::h5write(module_names, file = clusters_path, name = paste0(k, "/modules"))
    }
    on.exit(rhdf5::h5closeAll())
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
#' @param gene_umap_parameters A list of parameters that will be passed to the
#' `uwot::umap` function to generate the gene UMAP. Consult the documentation
#' of the function to obtain the list of parameters.
#' @param ecc_threshold The threshold applied on the Element-Centric Consistency
#' (ECC) score to determine if a number of cluster is stable or not. By default,
#' it is set to 0.9.
#' @param freq_threshold The threshold applied on the number of times a number
#' of clusters is obtained after running the `clustering_pipeline` function. By
#' default, it is set to 30.
#' @param scale_threshold The threshold applied on the scaled summarised
#' expression to determine the cell population associated to each gene module.
#' By default, it is set to 0.5.
#' @param enrichment_organism The organism used for the enrichment analysis.
#' @param save_entire_monocle If TRUE, saves the monocle object in the app
#' folder. This object could be used to perform additional changes to the app
#' or to continue the downstream analysis. If memory is a concern, set this
#' parameter to FALSE. Defaults to TRUE.
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
#'  - `monocle_object.qs2` - the monocle object that was used to generate the app.
#'  - `recommended_pseudotime.qs2` - the recommended pseudotime values.
#'  - `genes_info.csv` - a table with the gene information.
#'  - `trajectory_ggplot.qs2` - the ggplot object of the trajectory.
#'  - `metadata.qs2` - the metadata object.
#'  - `diet_monocle_object.qs2` - the monocle object with the diet information.
#'  - `expression.h5` - the expression matrix in the HDF5 format.
#'  - `colours.qs2` - the colours object.
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
                              gene_umap_parameters = list(
                                "min_dist" = 0.3,
                                "n_neighbors" = 30,
                                "metric" = "cosine",
                                "n_components" = 2
                              ),
                              ecc_threshold = 0.9,
                              freq_threshold = 30,
                              scale_threshold = 0.5,
                              enrichment_organism = "hsapiens",
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
    verbose_print("Creating app.R file...", verbose)
    starlng_write_app_file(file.path(folder_path, "app.R"), app_title_name, enrichment_organism = enrichment_organism)
    folder_path <- file.path(folder_path, "objects")
    if (!dir.exists(folder_path)) {
        dir.create(folder_path, showWarnings = FALSE)
    }

    for (mtd_column in colnames(monocle_object@colData)) {
        if (!(inherits(monocle_object@colData[[mtd_column]], c("factor", "character")))) {
            next
        }
        mtd_value <- monocle_object@colData[[mtd_column]]
        mtd_value <- as.character(mtd_value)
        mtd_value[is.na(mtd_value)] <- "N/A"
        mtd_value <- factor(mtd_value)

        if (nlevels(mtd_value) == nrow(monocle_object@colData)) {
            monocle_object@colData[[mtd_column]] <- NULL
            next
        }
        if (nlevels(mtd_value) > 200) {
            mtd_value <- as.numeric(mtd_value)
        }

        monocle_object@colData[[mtd_column]] <- mtd_value
    }

    # entire monocle object
    monocle_object <- do.call(
        custom_learn_graph,
        c(list("mon_obj" = monocle_object), learn_graph_parameters)
    )

    mon_path <- file.path(folder_path, "monocle_object.qs2")
    if (save_entire_monocle) {
        verbose_print("Writing the monocle object...", verbose)
        qs2::qs_save(monocle_object, file = mon_path, nthreads = nthreads)
    }

    # recommended pseudotime
    verbose_print("Writing the recommended pseudotime...", verbose)
    psd_path <- file.path(folder_path, "recommended_pseudotime.qs2")
    # TODO improve the pseudotime recomm
    psd_values <- get_pseudotime_recommendation(monocle_object)
    qs2::qs_save(psd_values, file = psd_path)
    gc()

    # autocorrelation object
    verbose_print("Writing the gene info table...", verbose)
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
    verbose_print("Extracting the trajectory information...", verbose)
    trajectory_path <- file.path(folder_path, "trajectory_object.qs2")
    trajectory_object <- get_trajectory_object(monocle_object, "UMAP")
    qs2::qs_save(trajectory_object, file = trajectory_path)

    # metadata object
    verbose_print("Writing the metadata object...", verbose)
    metadata_path <- file.path(folder_path, "metadata.qs2")
    mtd_df <- as.data.frame(monocle_object@colData)
    mtd_df$pseudotime <- psd_values$recommended_pseudotime[rownames(mtd_df)]
    mtd_df <- cbind(mtd_df, monocle_object@int_colData$reducedDims$UMAP)
    qs2::qs_save(mtd_df, file = metadata_path)
    
    # diet monocle object
    verbose_print("Writing the diet monocle object...", verbose)
    mon_path <- file.path(folder_path, "diet_monocle_object.qs2")
    monocle_object <- diet_monocle_object(monocle_object)
    gc()
    qs2::qs_save(monocle_object, file = mon_path, nthreads = nthreads)

    # gene clustering
    skip_clustering <- length(clustering_parameters) == 0
    if (!skip_clustering) {
        verbose_print("Applying clustering and stability assessment...", verbose)
        assess_path <- file.path(folder_path, "full_stability_assessment.qs2")

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
        verbose_print("Writing clustering results...", verbose)

        qs2::qs_save(clust_results, file = assess_path, nthreads = nthreads)
        best_config <- select_best_configuration(clust_results$clusters_list)
        message("Best configuration: n_neigh - ", best_config[[1]], ", quality function - ", best_config[[2]])
        qs2::qs_save(best_config, file = file.path(folder_path, "best_configuration.qs2"))
        selected_embedding <- clust_results$embedding_list
        selected_embedding$adj_matrix <- selected_embedding$adj_matrix[[as.character(best_config[[1]])]]
        gene_umap_parameters$X <- selected_embedding$embedding
        set.seed(42)
        selected_embedding$umap <- as.data.frame(do.call(
            uwot::umap,
            gene_umap_parameters
        ))
        colnames(selected_embedding$umap) <- c("UMAP_1", "UMAP_2")
        qs2::qs_save(selected_embedding, file = file.path(folder_path, "gene_embedding.qs2"), nthreads = nthreads)

        clust_results <- clust_results$clusters_list[[best_config[[1]]]][[best_config[[2]]]]

        stb_clust <- ClustAssess::choose_stable_clusters(
            clust_results$k,
            ecc_threshold = ecc_threshold,
            freq_threshold = freq_threshold
        )

        if (!is.null(stb_clust)) {
            starlng_write_module_summaries(
                stb_clust = stb_clust,
                folder_path = folder_path,
                expr_matrix = expr_matrix,
                chosen_genes = chosen_genes,
                gene_adj_matrix = selected_embedding$adj_matrix,
                trajectory_object = trajectory_object,
                cell_umap_df = mtd_df[, c(ncol(mtd_df) - 1, ncol(mtd_df))],
                scale_threshold = scale_threshold,
                enrichment_organism = enrichment_organism,
                verbose = verbose
            )
        }
    }

    # expression matrix
    verbose_print("Writing the expression matrix...", verbose)
    expr_path <- file.path(folder_path, "expression.h5")
    # write_gene_matrix_sparse_h5(
    write_gene_matrix_dense_h5(
        expr_matrix,
        expr_path,
        compression_level = compression_level,
        chunk_size = chunk_size,
        all_at_once = TRUE
    )

    # colours object
    verbose_print("Writing the colours object...", verbose)
    colour_path <- file.path(folder_path, "colours.qs2")
    colour_list <- list(
        "discrete" = generate_discrete_colours(
            mtd_df,
            discrete_colours,
            max_n_colors
        ),
        "continuous" = generate_continuous_colours(continuous_colours)
    )
    qs2::qs_save(colour_list, file = colour_path)
}

#' Create Starlng Shiny app from a normalized gene expression matrix
#'
#' @description Runs the entire Starlng pipeline on a normalized gene expression
#' matrix and generates a folder with all necessary files to run a Shiny app.
#'
#' @note For more details about the other parameters and the content of the
#' shiny folder, consult the documentation of the `starlng_write_app_monocle()`
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
#' @param folder_path Check `starlng_write_app_monocle()` documentation.
#' @param app_title_name Check `starlng_write_app_monocle()` documentation.
#' @param learn_graph_parameters Check `starlng_write_app_monocle()` documentation.
#' @param gene_filtering_function Check `starlng_write_app_monocle()` documentation.
#' @param clustering_parameters Check `starlng_write_app_monocle()` documentation.
#' @param ecc_threshold Check `starlng_write_app_monocle()` documentation.
#' @param freq_threshold Check `starlng_write_app_monocle()` documentation.
#' @param enrichment_organism Check `starlng_write_app_monocle()` documentation.
#' @param save_entire_monocle Check `starlng_write_app_monocle()` documentation.
#' @param discrete_colours Check `starlng_write_app_monocle()` documentation.
#' @param continuous_colours Check `starlng_write_app_monocle()` documentation.
#' @param max_n_colors Check `starlng_write_app_monocle()` documentation.
#' @param verbose Check `starlng_write_app_monocle()` documentation.
#' @param compression_level Check `starlng_write_app_monocle()` documentation.
#' @param chunk_size Check `starlng_write_app_monocle()` documentation.
#' @param nthreads Check `starlng_write_app_monocle()` documentation.
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
                              enrichment_organism = "hsapiens",
                              save_entire_monocle = TRUE,
                              discrete_colours = list(),
                              continuous_colours = list(),
                              max_n_colors = 40,
                              verbose = FALSE,
                              compression_level = 7,
                              chunk_size = 100,
                              nthreads = 1) {
    verbose_print("Creating the monocle object...", verbose)
    monocle_object <- ClustAssess::create_monocle_default(
        normalized_expression_matrix = expression_matrix,
        metadata_df = metadata_df,
        pca_embedding = pca_embedding,
        umap_embedding = umap_embedding
    )
    gc()

    starlng_write_app_monocle(
        folder_path = folder_path,
        monocle_object = monocle_object,
        app_title_name = app_title_name,
        learn_graph_parameters = learn_graph_parameters,
        gene_filtering_function = gene_filtering_function,
        clustering_parameters = clustering_parameters,
        ecc_threshold = ecc_threshold,
        freq_threshold = freq_threshold,
        enrichment_organism = enrichment_organism,
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
#' shiny folder, consult the documentation of the `starlng_write_app_monocle()`
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
#' @param folder_path Check `starlng_write_app_monocle()` documentation.
#' @param app_title_name Check `starlng_write_app_monocle()` documentation.
#' @param learn_graph_parameters Check `starlng_write_app_monocle()` documentation.
#' @param gene_filtering_function Check `starlng_write_app_monocle()` documentation.
#' @param clustering_parameters Check `starlng_write_app_monocle()` documentation.
#' @param ecc_threshold Check `starlng_write_app_monocle()` documentation.
#' @param freq_threshold Check `starlng_write_app_monocle()` documentation.
#' @param enrichment_organism Check `starlng_write_app_monocle()` documentation.`
#' @param save_entire_monocle Check `starlng_write_app_monocle()` documentation.
#' @param discrete_colours Check `starlng_write_app_monocle()` documentation.
#' @param continuous_colours Check `starlng_write_app_monocle()` documentation.
#' @param max_n_colors Check `starlng_write_app_monocle()` documentation.
#' @param verbose Check `starlng_write_app_monocle()` documentation.
#' @param compression_level Check `starlng_write_app_monocle()` documentation.
#' @param chunk_size Check `starlng_write_app_monocle()` documentation.
#' @param nthreads Check `starlng_write_app_monocle()` documentation.
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
                                  enrichment_organism = "hsapiens",
                                  save_entire_monocle = TRUE,
                                  discrete_colours = list(),
                                  continuous_colours = list(),
                                  max_n_colors = 40,
                                  verbose = FALSE,
                                  compression_level = 7,
                                  chunk_size = 100,
                                  nthreads = 1) {
    verbose_print("Creating the monocle object...", verbose)
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
    gc()

    starlng_write_app_monocle(
        folder_path = folder_path,
        monocle_object = monocle_object,
        app_title_name = app_title_name,
        learn_graph_parameters = learn_graph_parameters,
        gene_filtering_function = gene_filtering_function,
        clustering_parameters = clustering_parameters,
        ecc_threshold = ecc_threshold,
        freq_threshold = freq_threshold,
        enrichment_organism = enrichment_organism,
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
#' shiny folder, consult the documentation of the `starlng_write_app_monocle()`
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
#' @param folder_path Check `starlng_write_app_monocle()` documentation.
#' @param app_title_name Check `starlng_write_app_monocle()` documentation.`
#' @param learn_graph_parameters Check `starlng_write_app_monocle()` documentation.`
#' @param gene_filtering_function Check `starlng_write_app_monocle()` documentation.`
#' @param clustering_parameters Check `starlng_write_app_monocle()` documentation.`
#' @param ecc_threshold Check `starlng_write_app_monocle()` documentation.`
#' @param freq_threshold Check `starlng_write_app_monocle()` documentation.`
#' @param enrichment_organism Check `starlng_write_app_monocle()` documentation.`
#' @param save_entire_monocle Check `starlng_write_app_monocle()` documentation.`
#' @param discrete_colours Check `starlng_write_app_monocle()` documentation.`
#' @param continuous_colours Check `starlng_write_app_monocle()` documentation.`
#' @param max_n_colors Check `starlng_write_app_monocle()` documentation.`
#' @param verbose Check `starlng_write_app_monocle()` documentation.
#' @param compression_level Check `starlng_write_app_monocle()` documentation.`
#' @param chunk_size Check `starlng_write_app_monocle()` documentation.
#' @param nthreads Check `starlng_write_app_monocle()` documentation.
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
                                      enrichment_organism = "hsapiens",
                                      save_entire_monocle = TRUE,
                                      discrete_colours = list(),
                                      continuous_colours = list(),
                                      max_n_colors = 40,
                                      verbose = FALSE,
                                      compression_level = 7,
                                      chunk_size = 100,
                                      nthreads = 1) {
    verbose_print("Creating the monocle object...", verbose)
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
        enrichment_organism = enrichment_organism,
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
