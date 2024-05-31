library(shiny)
library(monocle3)
library(Seurat)
# library(ClustAssess)
library(foreach)

# seurat_clustassess_csv is a dataframe read from a csv that contains two columns:
# seurat - path to the Seurat object
# clustassess - path to the ClustAssess object
find_corresp_index <- function(target_id, seurat_clustassess_csv, prefix = "", suffix = "") {
    for (i in seq_len(nrow(seurat_clustassess_csv))) {
        if (!file.exists(seurat_clustassess_csv$seurat[i])) {
            next
        }

        if (!file.exists(seurat_clustassess_csv$clustassess[i])) {
            next
        }

        basename_seurat <- basename(seurat_clustassess_csv$seurat[i])

        if (basename_seurat == glue::glue("{prefix}{target_id}{suffix}.rds")) {
            return(i)
        }
    }

    return(-1)
}


# to concatenate them back, `do.call(file.path, as.list(dir_list))`
get_file_hierarchy <- function(file_path) {
    dir_list <- c()

    while (file_path != "/") {
        file_path <- dirname(file_path)

        dir_list <- c(basename(file_path), dir_list)
    }

    return(dir_list)
}

find_common_directory <- function(path_list, return_full = FALSE) {
    base_comp <- get_file_hierarchy(path_list[1])

    for (i in seq_along(path_list)) {
        if (i == 1) {
            next
        }

        current_path <- get_file_hierarchy(path_list[i])
        min_length <- min(length(base_comp), length(current_path))

        for (j in seq_len(min_length)) {
            if (base_comp[j] != current_path[j]) {
                j <- j - 1
                break
            }
        }

        base_comp <- base_comp[seq_len(j)]
    }

    if (return_full) {
        return(base_comp)
    }

    return(base_comp[length(base_comp)])
}


ppi <- 72
shiny_plot_k_n_partitions <- function(summary_df,
                                      plt_height,
                                      distance_factor = 0.6,
                                      text_size = 1,
                                      pt_size_range = c(0.1, 1.5),
                                      filtered_cl_methods = NULL,
                                      threshold_ecc = 0,
                                      threshold_occurences = 0,
                                      display_legend = FALSE,
                                      threshold_nparts = NULL) {
    available_pchs <- 21:24
    if (is.null(filtered_cl_methods)) {
        filtered_cl_methods <- c("RBConfigurationVertexPartition", "RBERVertexPartition", "CPMVertexPartition")
    }

    if (is.null(threshold_nparts)) {
        threshold_nparts <- max(summary_df$n_partitions) + 1
    }

    summary_df <- summary_df %>% dplyr::filter(
        .data$ecc > threshold_ecc &
            .data$cl_method %in% filtered_cl_methods &
            .data$n_partitions < threshold_nparts &
            .data$total_freq > threshold_occurences
    )

    if (nrow(summary_df) == 0) {
        return()
    }

    cl_method <- summary_df$cl_method
    unique_cl_method <- unique(cl_method)
    cl_method <- match(cl_method, unique_cl_method)
    pchs <- available_pchs[factor(cl_method)]
    n_methods <- length(unique_cl_method)
    cl_method <- cl_method - mean(seq_len(n_methods))
    mask <- (cl_method < 0)
    cl_method[mask] <- floor(cl_method[mask])
    cl_method[!mask] <- ceiling(cl_method[!mask])
    max_distance <- ifelse(n_methods > 1, 0.9 / (n_methods - n_methods %% 2), 0)

    color_values <- viridis::viridis(min(50, nrow(summary_df)))

    if (max(summary_df$ecc) - min(summary_df$ecc) < 1e-3) {
        color_info <- factor(rep(summary_df$ecc[1], nrow(summary_df)))
    } else {
        color_info <- cut(summary_df$ecc, breaks = min(50, nrow(summary_df)), dig.lab = 3)
    }

    color_info <- color_values[color_info]

    x_axis <- summary_df$k
    unique_x_axis <- unique(x_axis)
    x_axis <- match(x_axis, unique_x_axis)
    x_axis <- x_axis + cl_method * max_distance * distance_factor
    actual_min <- min(summary_df$total_freq)
    actual_max <- max(summary_df$total_freq)

    if (actual_max > actual_min) {
        pt_sizes <- (pt_size_range[2] - pt_size_range[1]) * (summary_df$total_freq - actual_min) / (actual_max - actual_min) + pt_size_range[1]
    } else {
        pt_sizes <- pt_size_range[2]
    }

    if (display_legend) {
        plt_height <- plt_height / ppi

        predicted_height <- graphics::strheight(paste(rep("TEXT", 4), collapse = "\n"), units = "inches", cex = text_size) + graphics::strheight("TE\nXT", units = "inches", cex = pt_size_range[2] / 1.5)
        layout.matrix <- matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, ncol = 3, byrow = TRUE)

        graphics::layout(
            mat = layout.matrix,
            heights = c(
                lcm((plt_height - predicted_height - 0.1) * 2.54),
                lcm(predicted_height * 2.54)
            )
        )
    }

    plot(
        x = x_axis,
        y = summary_df$n_partitions,
        pch = pchs,
        bg = color_info,
        xaxt = "n",
        xlab = "k",
        ylab = "# different partitions",
        cex = pt_sizes,
        cex.axis = text_size,
        cex.lab = text_size
    )
    abline(v = seq(from = 0.5, by = 1, length.out = length(unique_x_axis)), lty = "dashed", col = "grey")
    axis(side = 1, at = seq_along(unique_x_axis), labels = unique_x_axis, las = 2, cex.axis = text_size)

    if (display_legend) {
        # point shape
        plot(seq_len(n_methods) - 1, rep(1, n_methods), axes = FALSE, bty = "n", ylab = "", xlab = "", cex = pt_size_range[2], pch = available_pchs[seq_len(n_methods)], main = "clustering\nmethods", ylim = c(-1, 1.5), cex.main = text_size)
        text(y = -0.15, x = seq_len(n_methods) - 1, labels = unique_cl_method, cex = text_size, xpd = TRUE)


        # point size
        nfreqs <- sum(summary_df$total_freq)
        chosen_numbers <- seq(from = actual_min, to = actual_max, length.out = 5)
        if (actual_max > actual_min) {
            pt_sizes <- (pt_size_range[2] - pt_size_range[1]) * (chosen_numbers - actual_min) / (actual_max - actual_min) + pt_size_range[1]
        } else {
            pt_sizes <- rep(pt_size_range[2], 5)
        }
        plot(0:4, rep(1, 5), axes = FALSE, bty = "n", ylab = "", xlab = "", pch = 19, cex = pt_sizes, main = "k frequency", ylim = c(-1, 1.5), cex.main = text_size)
        text(y = -0.15, x = 0:4, labels = round(chosen_numbers / nfreqs, digits = 3), cex = text_size, xpd = TRUE)

        # point colour
        unique_values <- c(min(summary_df$ecc), max(summary_df$ecc))
        legend_image <- as.raster(matrix(color_values, nrow = 1))
        plot(c(0, 1), c(-1, 1), type = "n", axes = F, bty = "n", ylab = "", xlab = "", main = "partition\nfrequency", cex.main = text_size)
        text(y = -0.5, x = seq(0, 1, l = 5), labels = round(seq(from = unique_values[1], to = unique_values[2], length.out = 5), digits = 3), cex = text_size)
        rasterImage(legend_image, 0, 0, 1, 1)
    }
}


# adapted from the `cluster_genes` function of monocle3
# https://github.com/cole-trapnell-lab/monocle3/blob/b545460966874948eb11a57a225594a107f1694d/R/cluster_genes.R
my_cluster_genes <- function(mon_obj,
                             resolution,
                             subset_genes = NULL,
                             quality_function = c("RBConfigurationVertexPartition", "RBERVertexPartition", "CPMVertexPartition"),
                             num_iter = 10,
                             n_neighbours = 25, # TODO choose the most stable
                             graph_type = c("snn", "nn"), # TODO choose the most stable
                             graph_base_embedding = c("PCA", "UMAP"), # TODO choose the most stable, although it's most probable PCA
                             n_repetitions = 100,
                             seed_sequence = NULL,
                             umap_arguments = list(),
                             keep_first = FALSE,
                             verbose = FALSE) {
    graph_type <- graph_type[1]
    graph_base_embedding <- graph_base_embedding[1]
    quality_function <- quality_function[1]
    resolution <- sort(resolution)

    if (is.null(subset_genes)) {
        subset_genes <- rownames(mon_obj)
    }

    if (is.null(seed_sequence)) {
        seed_sequence <- seq(
            from = 1,
            by = 100,
            length.out = n_repetitions
        )
    } else {
        n_repetitions <- length(seed_sequence)
    }

    # get the number of availabel workers
    ncores <- foreach::getDoParWorkers()

    # create the loading matrix
    print(glue::glue("[{Sys.time()}] Calculate the gene loading matrix."))
    preprocess_mat <- mon_obj@reduce_dim_aux[["PCA"]][["model"]]$svd_v[subset_genes, ] %*% diag(mon_obj@reduce_dim_aux[["PCA"]][["model"]]$svd_sdev)

    # create the umap
    if (!("seed" %in% names(umap_arguments))) {
        umap_arguments$seed <- 42 # if the user doesn't provide a seed, we'll fix at 42 (Seurat's default)
    }

    print(glue::glue("[{Sys.time()}] Calculate the gene UMAP."))
    umap_res <- do.call(
        what = uwot::umap,
        args = c(
            list(X = as.matrix(preprocess_mat)),
            umap_arguments
        )
    )
    row.names(umap_res) <- row.names(preprocess_mat)
    colnames(umap_res) <- paste0("UMAP_", seq_len(ncol(umap_res)))

    needed_vars <- c("shared_adj", "resolution", "quality_function", "graph_type", "num_iter", "ncores")
    print(glue::glue("[{Sys.time()}] Calculate the (S)NN-graph."))
    # calculate the graph
    adj_matrix <- Seurat::FindNeighbors(
        object = switch(graph_base_embedding,
            "UMAP" = umap_res,
            "PCA" = preprocess_mat
        ),
        k.param = n_neighbours,
        verbose = FALSE,
        nn.method = "rann",
        compute.SNN = FALSE
    )$nn
    highest_prune_param <- 0

    if (graph_type == "snn") {
        highest_prune_param <- ClustAssess::get_highest_prune_param(
            adj_matrix,
            n_neighbours
        )
        adj_matrix <- highest_prune_param$adj_matrix
        rownames(adj_matrix) <- row.names(umap_res)
        colnames(adj_matrix) <- row.names(umap_res)
    }

    if (ncores > 1) {
        shared_adj <- SharedObject::share(adj_matrix)
    } else {
        shared_adj <- adj_matrix

        if (graph_type == "snn") {
            g <- igraph::graph_from_adjacency_matrix(
                adjmatrix = shared_adj,
                mode = "undirected",
                weighted = TRUE
            )
        } else {
            g <- igraph::graph_from_adjacency_matrix(
                adjmatrix = shared_adj,
                mode = "directed",
                weighted = FALSE
            )
        }
    }

    all_vars <- ls()

    print(glue::glue("[{Sys.time()}] Start gene clustering."))
    # apply the clustering
    clustering_results <- foreach::foreach(
        seed = seed_sequence,
        .inorder = FALSE,
        .noexport = all_vars[!(all_vars %in% needed_vars)]
    ) %dopar% {
        if (ncores > 1) {
            if (graph_type == "snn") {
                g <- igraph::graph_from_adjacency_matrix(
                    adjmatrix = shared_adj,
                    mode = "undirected",
                    weighted = TRUE
                )
            } else {
                g <- igraph::graph_from_adjacency_matrix(
                    adjmatrix = shared_adj,
                    mode = "directed",
                    weighted = FALSE
                )
            }
        }

        seed_result <- lapply(
            resolution,
            function(res) {
                leidenbase::leiden_find_partition(
                    igraph = g,
                    resolution_parameter = res,
                    partition_type = quality_function,
                    edge_weights = igraph::E(g)$weight,
                    num_iter = num_iter,
                    seed = seed
                )$membership
            }
        )

        names(seed_result) <- as.character(resolution)

        seed_result <- list(seed_result)
        names(seed_result) <- quality_function

        return(seed_result)
    }

    if (ncores > 1) {
        shared_adj <- SharedObject::unshare(adj_matrix)
    }

    print(glue::glue("[{Sys.time()}] Grouping the partitions by the number of clusters."))

    final_clustering_result <- list()
    final_summary_df <- NULL
    for (qfunc in names(clustering_results[[1]])) {
        final_clustering_result[[qfunc]] <- list()
        indices <- list()

        for (res in as.character(resolution)) {
            for (i in seq_len(n_repetitions)) {
                mb <- clustering_results[[i]][[qfunc]][[res]]
                k <- as.character(length(unique(mb)))

                if (!(k %in% names(final_clustering_result[[qfunc]]))) {
                    final_clustering_result[[qfunc]][[k]] <- list(mb)
                    indices[[k]] <- 2
                } else {
                    final_clustering_result[[qfunc]][[k]][[indices[[k]]]] <- mb
                    indices[[k]] <- indices[[k]] + 1
                }
            }
        }

        temp_summary_df <- data.frame(matrix(0, nrow = length(indices), ncol = 6))
        colnames(temp_summary_df) <- c("k", "cl_method", "n_partitions", "total_freq", "first_freq", "ecc")
        index <- 1

        for (k in names(final_clustering_result[[qfunc]])) {
            final_clustering_result[[qfunc]][[k]] <- ClustAssess::merge_partitions(final_clustering_result[[qfunc]][[k]])

            if (length(final_clustering_result[[qfunc]][[k]]) == 1) {
                temp_summary_df[index, ] <- c(k, qfunc, 1, final_clustering_result[[qfunc]][[k]][[1]]$freq, final_clustering_result[[qfunc]][[k]][[1]]$freq, 1)
                index <- index + 1

                # final_clustering_result[[qfunc]][[k]] <- final_clustering_result[[qfunc]][[k]][[1]]$mb
                next
            }

            ecs_sim_matrix <- ClustAssess::element_sim_matrix(lapply(final_clustering_result[[qfunc]][[k]], function(x) {
                x$mb
            }))

            max_freq <- final_clustering_result[[qfunc]][[k]][[1]]$freq

            for (i in seq(from = 2, to = nrow(ecs_sim_matrix), by = 1)) {
                if (final_clustering_result[[qfunc]][[k]][[i]]$freq != max_freq) {
                    i <- i - 1
                    break
                }
            }
            nmax <- i

            # calculate the weighted ECC of the partitions
            ecc_val <- 0
            total_weights <- 0
            for (i in seq_len(nrow(ecs_sim_matrix) - 1)) {
                w1 <- final_clustering_result[[qfunc]][[k]][[i]]$freq
                total_weights <- total_weights + w1
                ecc_val <- ecc_val + w1 * (w1 - 1) / 2 # ECS between identical partitions

                for (j in seq(from = i + 1, to = nrow(ecs_sim_matrix), by = 1)) {
                    w2 <- final_clustering_result[[qfunc]][[k]][[j]]$freq
                    ecc_val <- ecc_val + w1 * w2 * ecs_sim_matrix[i, j]
                }
            }

            i <- i + 1
            total_weights <- total_weights + final_clustering_result[[qfunc]][[k]][[i]]$freq
            ecc_val <- ecc_val / (total_weights * (total_weights - 1) / 2)

            # in case of ties on the highest frequency, put at the first position the partition that has the highest similarity with the others
            if (nmax > 1) {
                average_agreement <- (rowSums(ecs_sim_matrix, na.rm = TRUE) + colSums(ecs_sim_matrix, na.rm = TRUE) - 2) / (nrow(ecs_sim_matrix) - 1)
                highest_sim_index <- which.max(average_agreement[seq_len(nmax)])

                temp_mb <- final_clustering_result[[qfunc]][[k]][[1]]$mb
                final_clustering_result[[qfunc]][[k]][[1]]$mb <- final_clustering_result[[qfunc]][[k]][[highest_sim_index]]$mb
                final_clustering_result[[qfunc]][[k]][[highest_sim_index]]$mb <- temp_mb
            }


            sum_freqs <- sum(sapply(final_clustering_result[[qfunc]][[k]], function(x) {
                x$freq
            }))
            temp_summary_df[index, ] <- c(k, qfunc, length(final_clustering_result[[qfunc]][[k]]), sum_freqs, final_clustering_result[[qfunc]][[k]][[1]]$freq, ecc_val)
            index <- index + 1

            if (keep_first) {
                final_clustering_result[[qfunc]][[k]] <- final_clustering_result[[qfunc]][[k]][[1]]
            }
        }

        if (is.null(final_summary_df)) {
            final_summary_df <- temp_summary_df
        } else {
            final_summary_df <- rbind(final_summary_df, temp_summary_df)
        }
    }

    final_summary_df$k <- as.integer(final_summary_df$k)
    final_summary_df$n_partitions <- as.integer(final_summary_df$n_partitions)
    final_summary_df$total_freq <- as.integer(final_summary_df$total_freq)
    final_summary_df$first_freq <- as.integer(final_summary_df$first_freq)
    final_summary_df$ecc <- as.numeric(final_summary_df$ecc)

    print(glue::glue("[{Sys.time()}] Done."))
    return(
        list(
            partitions = final_clustering_result,
            summary = final_summary_df,
            gene_umap = umap_res
        )
    )
}


write_object <- function(mon_obj,
                         trajectory_id,
                         start_genes,
                         start_expression_thresh,
                         start_relax_ngenes,
                         end_genes,
                         end_expression_thresh,
                         end_relax_ngenes,
                         output_dir,
                         #  resolution,
                         #  quality_function = "RBConfigurationVertexPartition",
                         #  n_repetitions = 100,
                         #  n_neighbours = 25,
                         #  graph_type = "snn",
                         #  graph_base_embedding = "PCA",
                         #  umap_arguments = list(),
                         #  num_iter = 10,
                         gene_variance_threshold = 0,
                         use_closed_loops = FALSE,
                         use_partitions = FALSE,
                         learn_graph_controls = NULL,
                         nodes_per_log10_cells = 30,
                         chunk_size = 100,
                         compression_level = 6,
                         ncores = 1) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # select the cells expressing the starting genes
    if (length(start_genes) > 1) {
        color_info <- DelayedMatrixStats::colSums2(mon_obj@assays@data$counts[start_genes, ] > start_expression_thresh) >= (length(start_genes) - start_relax_ngenes)
    } else {
        color_info <- mon_obj@assays@data$counts[start_genes, ] > start_expression_thresh
    }

    umap_emb <- SingleCellExperiment::reducedDims(mon_obj)$UMAP[color_info, ]
    if (nrow(umap_emb) > 0) {
        gmedian <- Gmedian::Gmedian(umap_emb)
        overall_distance <- ((umap_emb[, 1] - gmedian[1, 1])^2 + (umap_emb[, 2] - gmedian[1, 2])^2)^0.5
        dist_stats <- fivenum(overall_distance)
        color_info[seq_len(ncol(mon_obj))[color_info][overall_distance > dist_stats[3]]] <- FALSE
    }

    start_cells <- colnames(mon_obj)[color_info]

    # select the cells expressing the end genes
    if (length(end_genes) > 1) {
        color_info <- DelayedMatrixStats::colSums2(mon_obj@assays@data$counts[end_genes, ] > end_expression_thresh) >= (length(end_genes) - end_relax_ngenes)
    } else {
        color_info <- mon_obj@assays@data$counts[end_genes, ] > end_expression_thresh
    }

    umap_emb <- SingleCellExperiment::reducedDims(mon_obj)$UMAP[color_info, ]
    if (nrow(umap_emb) > 0) {
        gmedian <- Gmedian::Gmedian(umap_emb)
        overall_distance <- ((umap_emb[, 1] - gmedian[1, 1])^2 + (umap_emb[, 2] - gmedian[1, 2])^2)^0.5
        dist_stats <- fivenum(overall_distance)
        color_info[seq_len(ncol(mon_obj))[color_info][overall_distance > dist_stats[3]]] <- FALSE
    }
    end_cells <- colnames(mon_obj)[color_info]


     # find the set of cells between start and stop on the trajectory graph
    nodes_start <- monocle3::principal_graph_aux(mon_obj)$UMAP$pr_graph_cell_proj_closest_vertex[, 1][start_cells]
    nodes_end <- monocle3::principal_graph_aux(mon_obj)$UMAP$pr_graph_cell_proj_closest_vertex[, 1][end_cells]

    cells_between_start <- monocle3::choose_graph_segments(
        mon_obj,
        starting_pr_node = nodes_start[1],
        ending_pr_nodes = nodes_end,
        clear_cds = FALSE
    )

    cells_between_end <- monocle3::choose_graph_segments(
        mon_obj,
        starting_pr_node = nodes_end[1],
        ending_pr_nodes = nodes_start,
        clear_cds = FALSE
    )

    all_vertices_between <- union(
        igraph::vertex_attr(cells_between_end@principal_graph$UMAP, "name"),
        igraph::vertex_attr(cells_between_start@principal_graph$UMAP, "name")
    )

    all_cells_between <- union(
        colnames(cells_between_start),
        colnames(cells_between_end)
    )

    # subset the Monocle3 object using only this subtrajectory
    sub_mon_obj <- mon_obj[, all_cells_between]
    for (mtd in colnames(mon_obj@colData)) {
        if (!(is.character(sub_mon_obj@colData[[mtd]]) || is.factor(sub_mon_obj@colData[[mtd]]))) {
            next
        }

        sub_mon_obj@colData[[mtd]] <- factor(as.character(sub_mon_obj@colData[[mtd]]))
    }

    # keep the associated trajectory graph
    monocle3::principal_graph(sub_mon_obj)$UMAP <- igraph::induced_subgraph(monocle3::principal_graph(sub_mon_obj)$UMAP, all_vertices_between)

    # refactor the clusters
    sub_mon_obj@clusters@listData$UMAP$clusters <- setNames(factor(as.numeric(monocle3::clusters(sub_mon_obj))), colnames(sub_mon_obj))
    levels(sub_mon_obj@clusters@listData$UMAP$clusters) <- seq_len(nlevels(sub_mon_obj@clusters@listData$UMAP$clusters))
    sub_mon_obj@clusters@listData$UMAP$partitions <- setNames(factor(as.numeric(monocle3::partitions(sub_mon_obj))), colnames(sub_mon_obj))
    levels(sub_mon_obj@clusters@listData$UMAP$partitions) <- seq_len(nlevels(sub_mon_obj@clusters@listData$UMAP$partitions))

    # learn the trajectory graph on the subset
    sub_mon_obj <- monocle3::learn_graph(
        sub_mon_obj,
        close_loop = use_closed_loops,
        use_partition = use_partitions,
        learn_graph_control = learn_graph_controls
    )

    # calculate pseudotime
    sub_mon_obj <- monocle3::order_cells(sub_mon_obj, root_cells = start_cells[start_cells %in% all_cells_between])

    saveRDS(sub_mon_obj, file.path(output_dir, glue::glue("monocle_object.rds")))
    write.csv(monocle3::graph_test(sub_mon_obj, "principal_graph", cores = ncores), file.path(output_dir, "graph_test.csv"))

    # create and save the expression matrix
    expr_file_name <- file.path(output_dir, "expression.h5")
    expression_matrix <- as.matrix(sub_mon_obj@assays@data$normalized_data)
    gene_vars <- matrixStats::rowVars(expression_matrix)

    ## filter by var threshold
    expression_matrix <- expression_matrix[gene_vars >= gene_variance_threshold, ]

    gene_avg_expression <- matrixStats::rowMeans2(expression_matrix)
    order_expression <- order(gene_avg_expression, decreasing = TRUE)
    ## order the genes based on their average expression
    expression_matrix <- expression_matrix[order_expression, ]
    gene_avg_expression <- gene_avg_expression[order_expression]

    print(glue::glue("[{Sys.time()}] Writing the gene matrix"))

    rhdf5::h5createFile(expr_file_name)

    rhdf5::h5createDataset(expr_file_name, "rank_matrix",
        level = compression_level,
        dims = c(nrow(expression_matrix), ncol(expression_matrix)),
        storage.mode = "integer",
        chunk = c(chunk_size, ncol(expression_matrix))
    )

    rhdf5::h5createDataset(expr_file_name, "expression_matrix",
        level = compression_level,
        dims = c(nrow(expression_matrix), ncol(expression_matrix)),
        maxdims = c(nrow(expression_matrix), ncol(expression_matrix)),
        storage.mode = "double",
        chunk = c(1, ncol(expression_matrix))
    )

    rhdf5::h5write(rownames(expression_matrix), expr_file_name, "genes")
    rhdf5::h5write(colnames(expression_matrix), expr_file_name, "cells")
    rhdf5::h5write(gene_avg_expression, expr_file_name, "average_expression")
    rhdf5::h5write(chunk_size, expr_file_name, "chunk_size")
    rhdf5::h5write(expression_matrix, expr_file_name, "expression_matrix")

    # create and save the expression rank matrix
    print(glue::glue("[{Sys.time()}] Writing the rank matrix"))
    nchunks <- ceiling(nrow(expression_matrix) / chunk_size)
    for (i in seq_len(nchunks)) {
        if (i == nchunks) {
            index_list <- seq(from = chunk_size * (i - 1) + 1, by = 1, to = nrow(expression_matrix))
        } else {
            index_list <- seq(from = chunk_size * (i - 1) + 1, by = 1, length.out = chunk_size)
        }

        rhdf5::h5write(
            t(apply(
                expression_matrix[index_list, ],
                1,
                function(x) {
                    rank(x, ties.method = "min")
                }
            )),
            expr_file_name,
            "rank_matrix",
            index = list(index_list, NULL)
        )
    }

    write_app_file(output_dir, trajectory_id)
    styler::style_file(file.path(output_dir, "app.R"), indent_by = 4)

    print(glue::glue("[{Sys.time()}] Done."))
}

write_object_old <- function(original_output_dir, id, trajectory_id, paths_metadata, ncores = 1) {
    samples_metadata <- read.csv(paths_metadata, comment.char = "#")

    output_dir <- file.path(original_output_dir, id)
    dir.create(output_dir, showWarnings = FALSE)
    output_dir <- file.path(output_dir, trajectory_id)
    dir.create(output_dir, showWarnings = FALSE)

    found_id <- FALSE

    print("Searching")
    for (i in seq_along(samples_metadata[, 1])) {
        sid <- samples_metadata$SampleID[i]
        tid <- samples_metadata$PseudotimeID[i]

        if (sid != id || trajectory_id != tid) {
            next
        }
        print(samples_metadata[i, ])

        found_id <- TRUE
        ftype <- samples_metadata$Ftype[i]
        fsize <- samples_metadata$Fsize[i]
        clmethod <- samples_metadata$Cl_method[i]
        k <- samples_metadata$k[i]
        seurat_path <- samples_metadata$Seurat_path[i]
        clustassess_path <- samples_metadata$ClustAssess_path[i]

        start_genes <- strsplit(samples_metadata$StartGenes[i], ";")[[1]]
        start_expression <- samples_metadata$ThresholdStart[i]
        start_relax <- samples_metadata$RelaxStart[i]
        end_genes <- strsplit(samples_metadata$EndGenes[i], ";")[[1]]
        end_expression <- samples_metadata$ThresholdEnd[i]
        end_relax <- samples_metadata$RelaxEnd[i]
    }

    if (!found_id) {
        stop(glue::glue("Id {id} not found"))
    }

    seurat_obj <- readRDS(seurat_path)
    clustassess_obj <- readRDS(clustassess_path)
    assay_name <- "SCT"

    mon_obj <- ClustAssess::create_monocle_object(
        normalized_expression_matrix = seurat_obj@assays[[assay_name]]@data,
        clustassess_object = clustassess_obj,
        metadata = seurat_obj@meta.data,
        stable_feature_type = ftype,
        stable_feature_set_size = fsize,
        stable_clustering_method = clmethod,
        stable_n_clusters = k,
        cell_names = colnames(seurat_obj),
        gene_names = rownames(seurat_obj),
        use_all_genes = TRUE
    )
    seurat_obj@reductions$umap@cell.embeddings <- mon_obj@int_colData$reducedDims$UMAP

    mon_obj <- monocle3::cluster_cells(mon_obj)
    mon_obj@clusters@listData$UMAP$clusters <- setNames(mon_obj@colData[[paste0("stable_", k[1], "_clusters")]], colnames(mon_obj))
    mon_obj@clusters@listData$UMAP$partitions <- setNames(mon_obj@colData[[paste0("stable_", k[1], "_clusters")]], colnames(mon_obj))
    mon_obj <- monocle3::learn_graph(mon_obj, close_loop = FALSE, use_partition = FALSE)

    if (length(start_genes) > 1) {
        color_info <- colSums2(mon_obj@assays@data$counts[start_genes, ] > start_expression) >= (length(start_genes) - start_relax)
    } else {
        color_info <- mon_obj@assays@data$counts[start_genes, ] > start_expression
    }
    umap_emb <- SingleCellExperiment::reducedDims(mon_obj)$UMAP[color_info, ]
    if (nrow(umap_emb) > 0) {
        gmedian <- Gmedian::Gmedian(umap_emb)
        overall_distance <- ((umap_emb[, 1] - gmedian[1, 1])^2 + (umap_emb[, 2] - gmedian[1, 2])^2)^0.5
        dist_stats <- fivenum(overall_distance)
        color_info[seq_len(ncol(mon_obj))[color_info][overall_distance > dist_stats[3]]] <- FALSE
    }

    start_cells <- colnames(mon_obj)[color_info]

    if (length(end_genes) > 1) {
        color_info <- colSums2(mon_obj@assays@data$counts[end_genes, ] > end_expression) >= (length(end_genes) - end_relax)
    } else {
        color_info <- mon_obj@assays@data$counts[end_genes, ] > end_expression
    }

    umap_emb <- SingleCellExperiment::reducedDims(mon_obj)$UMAP[color_info, ]
    if (nrow(umap_emb) > 0) {
        gmedian <- Gmedian::Gmedian(umap_emb)
        overall_distance <- ((umap_emb[, 1] - gmedian[1, 1])^2 + (umap_emb[, 2] - gmedian[1, 2])^2)^0.5
        dist_stats <- fivenum(overall_distance)
        color_info[seq_len(ncol(mon_obj))[color_info][overall_distance > dist_stats[3]]] <- FALSE
    }
    end_cells <- colnames(mon_obj)[color_info]

    print(start_cells)
    print(end_cells)

    nodes_start <- principal_graph_aux(mon_obj)$UMAP$pr_graph_cell_proj_closest_vertex[, 1][start_cells]
    nodes_end <- principal_graph_aux(mon_obj)$UMAP$pr_graph_cell_proj_closest_vertex[, 1][end_cells]

    cells_between_start <- choose_graph_segments(
        mon_obj,
        starting_pr_node = nodes_start[1],
        ending_pr_nodes = nodes_end,
        clear_cds = FALSE
    )

    cells_between_end <- choose_graph_segments(
        mon_obj,
        starting_pr_node = nodes_end[1],
        ending_pr_nodes = nodes_start,
        clear_cds = FALSE
    )

    all_vertices_between <- union(
        igraph::vertex_attr(cells_between_end@principal_graph$UMAP, "name"),
        igraph::vertex_attr(cells_between_start@principal_graph$UMAP, "name")
    )

    all_cells_between <- union(
        colnames(cells_between_start),
        colnames(cells_between_end)
    )

    sub_mon_obj <- mon_obj[, all_cells_between]
    for (mtd in colnames(mon_obj@colData)) {
        if (!(is.character(sub_mon_obj@colData[[mtd]]) || is.factor(sub_mon_obj@colData[[mtd]]))) {
            next
        }

        sub_mon_obj@colData[[mtd]] <- factor(as.character(sub_mon_obj@colData[[mtd]]))
    }

    principal_graph(sub_mon_obj)$UMAP <- igraph::induced_subgraph(principal_graph(sub_mon_obj)$UMAP, all_vertices_between)
    sub_mon_obj@clusters@listData$UMAP$clusters <- setNames(factor(as.numeric(monocle3::clusters(sub_mon_obj))), colnames(sub_mon_obj))
    levels(sub_mon_obj@clusters@listData$UMAP$clusters) <- seq_len(nlevels(sub_mon_obj@clusters@listData$UMAP$clusters))
    sub_mon_obj@clusters@listData$UMAP$partitions <- setNames(factor(as.numeric(monocle3::partitions(sub_mon_obj))), colnames(sub_mon_obj))
    levels(sub_mon_obj@clusters@listData$UMAP$partitions) <- seq_len(nlevels(sub_mon_obj@clusters@listData$UMAP$partitions))
    sub_mon_obj <- learn_graph(sub_mon_obj, close_loop = FALSE, use_partition = FALSE)
    sub_mon_obj <- order_cells(sub_mon_obj, root_cells = start_cells[start_cells %in% all_cells_between])

    saveRDS(sub_mon_obj, file.path(output_dir, glue::glue("monocle_object.rds")))
    write.csv(graph_test(sub_mon_obj, "principal_graph", cores = ncores), file.path(output_dir, "graph_test.csv"))

    write_app_file(output_dir, paste(id, trajectory_id, sep = "-"))
    styler::style_file(file.path(output_dir, "app.R"), indent_by = 4)
}


write_app_file <- function(output_dir, id) {
    file_content <- paste0("
        library(monocle3)
        library(shiny)
        library(shinyWidgets)
        library(plotly)
        library(shinyjs)
        library(ggplot2)
        library(ggplotify)
        library(dplyr)
        library(tidyr)
        library(foreach)

        package_env <- new.env(parent = emptyenv())
        ppi <- 72

        # adapted from ClustAssess
        shiny_plot_k_n_partitions <- function(summary_df,
                                            plt_height,
                                            distance_factor = 0.6,
                                            text_size = 1,
                                            pt_size_range = c(0.1, 1.5),
                                            filtered_cl_methods = NULL,
                                            threshold_ecc = 0,
                                            threshold_occurences = 0,
                                            plot_title = \"\",
                                            display_legend = FALSE,
                                            threshold_nparts = NULL) {
            available_pchs <- 21:24
            if (is.null(filtered_cl_methods)) {
                filtered_cl_methods <- c(\"RBConfigurationVertexPartition\", \"RBERVertexPartition\", \"CPMVertexPartition\")
            }

            if (is.null(threshold_nparts)) {
                threshold_nparts <- max(summary_df$n_partitions) + 1
            }

            summary_df <- summary_df %>% dplyr::filter(
                .data$ecc > threshold_ecc &
                    .data$cl_method %in% filtered_cl_methods &
                    .data$n_partitions < threshold_nparts &
                    .data$total_freq > threshold_occurences
            )

            if (nrow(summary_df) == 0) {
                return()
            }

            cl_method <- summary_df$cl_method
            unique_cl_method <- unique(cl_method)
            cl_method <- match(cl_method, unique_cl_method)
            pchs <- available_pchs[factor(cl_method)]
            n_methods <- length(unique_cl_method)
            cl_method <- cl_method - mean(seq_len(n_methods))
            mask <- (cl_method < 0)
            cl_method[mask] <- floor(cl_method[mask])
            cl_method[!mask] <- ceiling(cl_method[!mask])
            max_distance <- ifelse(n_methods > 1, 0.9 / (n_methods - n_methods %% 2), 0)

            color_values <- viridis::viridis(min(50, nrow(summary_df)))

            if (max(summary_df$ecc) - min(summary_df$ecc) < 1e-3) {
                color_info <- factor(rep(summary_df$ecc[1], nrow(summary_df)))
            } else {
                color_info <- cut(summary_df$ecc, breaks = min(50, nrow(summary_df)), dig.lab = 3)
            }

            color_info <- color_values[color_info]

            x_axis <- summary_df$k
            unique_x_axis <- unique(x_axis)
            x_axis <- match(x_axis, unique_x_axis)
            x_axis <- x_axis + cl_method * max_distance * distance_factor
            actual_min <- min(summary_df$total_freq)
            actual_max <- max(summary_df$total_freq)

            if (actual_max > actual_min) {
                pt_sizes <- (pt_size_range[2] - pt_size_range[1]) * (summary_df$total_freq - actual_min) / (actual_max - actual_min) + pt_size_range[1]
            } else {
                pt_sizes <- pt_size_range[2]
            }

            if (display_legend) {
                plt_height <- plt_height / ppi

                predicted_height <- graphics::strheight(paste(rep(\"TEXT\", 4), collapse = \"\n\"), units = \"inches\", cex = text_size) + graphics::strheight(\"TE\nXT\", units = \"inches\", cex = pt_size_range[2] / 1.5)
                layout.matrix <- matrix(c(1, 1, 1, 2, 3, 4), nrow = 2, ncol = 3, byrow = TRUE)

                graphics::layout(
                    mat = layout.matrix,
                    heights = c(
                        lcm((plt_height - predicted_height - 0.1) * 2.54),
                        lcm(predicted_height * 2.54)
                    )
                )
            }

            plot(
                x = x_axis,
                y = summary_df$n_partitions,
                pch = pchs,
                bg = color_info,
                main = plot_title,
                xaxt = \"n\",
                xlab = \"k\",
                ylab = \"# different partitions\",
                cex = pt_sizes,
                cex.main = text_size,
                cex.axis = text_size,
                cex.lab = text_size
            )
            abline(v = seq(from = 0.5, by = 1, length.out = length(unique_x_axis)), lty = \"dashed\", col = \"grey\")
            axis(side = 1, at = seq_along(unique_x_axis), labels = unique_x_axis, las = 2, cex.axis = text_size)

            if (display_legend) {
                # point shape
                plot(seq_len(n_methods) - 1, rep(1, n_methods), axes = FALSE, bty = \"n\", ylab = \"\", xlab = \"\", cex = pt_size_range[2], pch = available_pchs[seq_len(n_methods)], main = \"clustering\nmethods\", ylim = c(-1, 1.5), cex.main = text_size)
                text(y = -0.15, x = seq_len(n_methods) - 1, labels = unique_cl_method, cex = text_size, xpd = TRUE)

                # point size
                nfreqs <- sum(summary_df$total_freq)
                chosen_numbers <- seq(from = actual_min, to = actual_max, length.out = 5)
                if (actual_max > actual_min) {
                    pt_sizes <- (pt_size_range[2] - pt_size_range[1]) * (chosen_numbers - actual_min) / (actual_max - actual_min) + pt_size_range[1]
                } else {
                    pt_sizes <- rep(pt_size_range[2], 5)
                }
                plot(0:4, rep(1, 5), axes = FALSE, bty = \"n\", ylab = \"\", xlab = \"\", pch = 19, cex = pt_sizes, main = \"k frequency\", ylim = c(-1, 1.5), cex.main = text_size)
                text(y = -0.15, x = 0:4, labels = round(chosen_numbers / nfreqs, digits = 3), cex = text_size, xpd = TRUE)

                # point colour
                unique_values <- c(min(summary_df$ecc), max(summary_df$ecc))
                legend_image <- as.raster(matrix(color_values, nrow = 1))
                plot(c(0, 1), c(-1, 1), type = \"n\", axes = F, bty = \"n\", ylab = \"\", xlab = \"\", main = \"partition\nfrequency\", cex.main = text_size)
                text(y = -0.5, x = seq(0, 1, l = 5), labels = round(seq(from = unique_values[1], to = unique_values[2], length.out = 5), digits = 3), cex = text_size)
                rasterImage(legend_image, 0, 0, 1, 1)
            }
        }

        init_env <- function(dimension_reactive) {
            rlang::env_bind(
                .env = package_env,
                mon_obj = readRDS(\"monocle_object.rds\"),
                dimension = dimension_reactive,
                filtered_graph_test_res = shiny::reactiveVal(NULL),
                gene_module = shiny::reactiveVal(NULL)
            )

            discrete_metadata <- sapply(colnames(package_env$mon_obj@colData), function(x) {
                is.character(package_env$mon_obj@colData[[x]]) || is.factor(package_env$mon_obj@colData[[x]])
            })
            metadata_names <- colnames(package_env$mon_obj@colData)
            discrete_metadata <- names(discrete_metadata[discrete_metadata])
            gene_list <- rhdf5::h5read(\"expression.h5\", \"genes\")
            umap_df <- data.frame(SingleCellExperiment::reducedDims(package_env$mon_obj)$UMAP)

            graph_test_res <- read.table(\"graph_test.csv\", sep = \",\", header = TRUE, row.names = 1)[c(3, 4, 2, 6)][gene_list, ]
            graph_test_res$avg_expression <- round(rhdf5::h5read(\"expression.h5\", \"average_expression\"), digits = 5)
            graph_test_res$module <- NA

            rlang::env_bind(
                .env = package_env,
                discrete_metadata = discrete_metadata,
                metadata_names = metadata_names,
                gene_list = gene_list,
                umap_df = umap_df,
                height_ratio = 0.65,
                graph_test_res = graph_test_res,
            )

            rlang::env_bind(
                .env = package_env,
                plt_height = shiny::reactive(
                    floor(package_env$height_ratio * dimension_reactive()[2])
                ),
                plt_width = shiny::reactive(dimension_reactive()[1]),
                plt_width_half = shiny::reactive(dimension_reactive()[1] / 2.1)
            )
        }

        # adapted from the `cluster_genes` function of monocle3
        # https://github.com/cole-trapnell-lab/monocle3/blob/b545460966874948eb11a57a225594a107f1694d/R/cluster_genes.R
        my_cluster_genes <- function(mon_obj,
                                    resolution,
                                    subset_genes = NULL,
                                    quality_function = c(\"RBConfigurationVertexPartition\", \"RBERVertexPartition\", \"CPMVertexPartition\"),
                                    num_iter = 10,
                                    n_neighbours = 25, # TODO choose the most stable
                                    graph_type = c(\"snn\", \"nn\"), # TODO choose the most stable
                                    graph_base_embedding = c(\"PCA\", \"UMAP\"), # TODO choose the most stable, although it's most probable PCA
                                    n_repetitions = 100,
                                    seed_sequence = NULL,
                                    umap_arguments = list(),
                                    keep_first = FALSE,
                                    verbose = FALSE) {
            graph_type <- graph_type[1]
            graph_base_embedding <- graph_base_embedding[1]
            quality_function <- quality_function[1]
            resolution <- sort(resolution)

            if (is.null(subset_genes)) {
                subset_genes <- rownames(mon_obj)
            }

            if (is.null(seed_sequence)) {
                seed_sequence <- seq(
                    from = 1,
                    by = 100,
                    length.out = n_repetitions
                )
            } else {
                n_repetitions <- length(seed_sequence)
            }

            # get the number of availabel workers
            ncores <- foreach::getDoParWorkers()

            # create the loading matrix
            print(glue::glue(\"[{Sys.time()}] Calculate the gene loading matrix.\"))
            preprocess_mat <- mon_obj@reduce_dim_aux[[\"PCA\"]][[\"model\"]]$svd_v[subset_genes, ] %*% diag(mon_obj@reduce_dim_aux[[\"PCA\"]][[\"model\"]]$svd_sdev)

            # create the umap
            if (!(\"seed\" %in% names(umap_arguments))) {
                umap_arguments$seed <- 42 # if the user doesn't provide a seed, we'll fix at 42 (Seurat's default)
            }

            print(glue::glue(\"[{Sys.time()}] Calculate the gene UMAP.\"))
            umap_res <- do.call(
                what = uwot::umap,
                args = c(
                    list(X = as.matrix(preprocess_mat)),
                    umap_arguments
                )
            )
            row.names(umap_res) <- row.names(preprocess_mat)
            colnames(umap_res) <- paste0(\"UMAP_\", seq_len(ncol(umap_res)))

            needed_vars <- c(\"shared_adj\", \"resolution\", \"quality_function\", \"graph_type\", \"num_iter\", \"ncores\")
            print(glue::glue(\"[{Sys.time()}] Calculate the (S)NN-graph.\"))
            # calculate the graph
            adj_matrix <- Seurat::FindNeighbors(
                object = switch(graph_base_embedding,
                    \"UMAP\" = umap_res,
                    \"PCA\" = preprocess_mat
                ),
                k.param = n_neighbours,
                verbose = FALSE,
                nn.method = \"rann\",
                compute.SNN = FALSE
            )$nn
            highest_prune_param <- 0

            if (graph_type == \"snn\") {
                highest_prune_param <- ClustAssess::get_highest_prune_param(
                    adj_matrix,
                    n_neighbours
                )
                adj_matrix <- highest_prune_param$adj_matrix
                rownames(adj_matrix) <- row.names(umap_res)
                colnames(adj_matrix) <- row.names(umap_res)
            }

            if (ncores > 1) {
                shared_adj <- SharedObject::share(adj_matrix)
            } else {
                shared_adj <- adj_matrix

                if (graph_type == \"snn\") {
                    g <- igraph::graph_from_adjacency_matrix(
                        adjmatrix = shared_adj,
                        mode = \"undirected\",
                        weighted = TRUE
                    )
                } else {
                    g <- igraph::graph_from_adjacency_matrix(
                        adjmatrix = shared_adj,
                        mode = \"directed\",
                        weighted = FALSE
                    )
                }
            }

            all_vars <- ls()

            print(glue::glue(\"[{Sys.time()}] Start gene clustering.\"))
            # apply the clustering
            clustering_results <- foreach::foreach(
                seed = seed_sequence,
                .inorder = FALSE,
                .noexport = all_vars[!(all_vars %in% needed_vars)]
            ) %dopar% {
                if (ncores > 1) {
                    if (graph_type == \"snn\") {
                        g <- igraph::graph_from_adjacency_matrix(
                            adjmatrix = shared_adj,
                            mode = \"undirected\",
                            weighted = TRUE
                        )
                    } else {
                        g <- igraph::graph_from_adjacency_matrix(
                            adjmatrix = shared_adj,
                            mode = \"directed\",
                            weighted = FALSE
                        )
                    }
                }

                seed_result <- lapply(
                    resolution,
                    function(res) {
                        leidenbase::leiden_find_partition(
                            igraph = g,
                            resolution_parameter = res,
                            partition_type = quality_function,
                            edge_weights = igraph::E(g)$weight,
                            num_iter = num_iter,
                            seed = seed
                        )$membership
                    }
                )

                names(seed_result) <- as.character(resolution)

                seed_result <- list(seed_result)
                names(seed_result) <- quality_function

                return(seed_result)
            }

            if (ncores > 1) {
                shared_adj <- SharedObject::unshare(adj_matrix)
            }

            print(glue::glue(\"[{Sys.time()}] Grouping the partitions by the number of clusters.\"))

            final_clustering_result <- list()
            final_summary_df <- NULL
            for (qfunc in names(clustering_results[[1]])) {
                final_clustering_result[[qfunc]] <- list()
                indices <- list()

                for (res in as.character(resolution)) {
                    for (i in seq_len(n_repetitions)) {
                        mb <- clustering_results[[i]][[qfunc]][[res]]
                        k <- as.character(length(unique(mb)))

                        if (!(k %in% names(final_clustering_result[[qfunc]]))) {
                            final_clustering_result[[qfunc]][[k]] <- list(mb)
                            indices[[k]] <- 2
                        } else {
                            final_clustering_result[[qfunc]][[k]][[indices[[k]]]] <- mb
                            indices[[k]] <- indices[[k]] + 1
                        }
                    }
                }

                temp_summary_df <- data.frame(matrix(0, nrow = length(indices), ncol = 6))
                colnames(temp_summary_df) <- c(\"k\", \"cl_method\", \"n_partitions\", \"total_freq\", \"first_freq\", \"ecc\")
                index <- 1

                for (k in names(final_clustering_result[[qfunc]])) {
                    final_clustering_result[[qfunc]][[k]] <- ClustAssess::merge_partitions(final_clustering_result[[qfunc]][[k]])

                    if (length(final_clustering_result[[qfunc]][[k]]) == 1) {
                        temp_summary_df[index, ] <- c(k, qfunc, 1, final_clustering_result[[qfunc]][[k]][[1]]$freq, final_clustering_result[[qfunc]][[k]][[1]]$freq, 1)
                        index <- index + 1

                        next
                    }

                    ecs_sim_matrix <- ClustAssess::element_sim_matrix(lapply(final_clustering_result[[qfunc]][[k]], function(x) {
                        x$mb
                    }))

                    max_freq <- final_clustering_result[[qfunc]][[k]][[1]]$freq

                    for (i in seq(from = 2, to = nrow(ecs_sim_matrix), by = 1)) {
                        if (final_clustering_result[[qfunc]][[k]][[i]]$freq != max_freq) {
                            i <- i - 1
                            break
                        }
                    }
                    nmax <- i

                    # calculate the weighted ECC of the partitions
                    ecc_val <- 0
                    total_weights <- 0
                    for (i in seq_len(nrow(ecs_sim_matrix) - 1)) {
                        w1 <- final_clustering_result[[qfunc]][[k]][[i]]$freq
                        total_weights <- total_weights + w1
                        ecc_val <- ecc_val + w1 * (w1 - 1) / 2 # ECS between identical partitions

                        for (j in seq(from = i + 1, to = nrow(ecs_sim_matrix), by = 1)) {
                            w2 <- final_clustering_result[[qfunc]][[k]][[j]]$freq
                            ecc_val <- ecc_val + w1 * w2 * ecs_sim_matrix[i, j]
                        }
                    }

                    i <- i + 1
                    total_weights <- total_weights + final_clustering_result[[qfunc]][[k]][[i]]$freq
                    ecc_val <- ecc_val / (total_weights * (total_weights - 1) / 2)

                    # in case of ties on the highest frequency, put at the first position the partition that has the highest similarity with the others
                    if (nmax > 1) {
                        average_agreement <- (rowSums(ecs_sim_matrix, na.rm = TRUE) + colSums(ecs_sim_matrix, na.rm = TRUE) - 2) / (nrow(ecs_sim_matrix) - 1)
                        highest_sim_index <- which.max(average_agreement[seq_len(nmax)])

                        temp_mb <- final_clustering_result[[qfunc]][[k]][[1]]$mb
                        final_clustering_result[[qfunc]][[k]][[1]]$mb <- final_clustering_result[[qfunc]][[k]][[highest_sim_index]]$mb
                        final_clustering_result[[qfunc]][[k]][[highest_sim_index]]$mb <- temp_mb
                    }


                    sum_freqs <- sum(sapply(final_clustering_result[[qfunc]][[k]], function(x) {
                        x$freq
                    }))
                    temp_summary_df[index, ] <- c(k, qfunc, length(final_clustering_result[[qfunc]][[k]]), sum_freqs, final_clustering_result[[qfunc]][[k]][[1]]$freq, ecc_val)
                    index <- index + 1

                    if (keep_first) {
                        final_clustering_result[[qfunc]][[k]] <- final_clustering_result[[qfunc]][[k]][1]
                    }
                }

                if (is.null(final_summary_df)) {
                    final_summary_df <- temp_summary_df
                } else {
                    final_summary_df <- rbind(final_summary_df, temp_summary_df)
                }
            }

            final_summary_df$k <- as.integer(final_summary_df$k)
            final_summary_df$n_partitions <- as.integer(final_summary_df$n_partitions)
            final_summary_df$total_freq <- as.integer(final_summary_df$total_freq)
            final_summary_df$first_freq <- as.integer(final_summary_df$first_freq)
            final_summary_df$ecc <- as.numeric(final_summary_df$ecc)

            print(glue::glue(\"[{Sys.time()}] Done.\"))
            return(
                list(
                    partitions = final_clustering_result,
                    summary = final_summary_df,
                    gene_umap = data.frame(umap_res)
                )
            )
        }

        gear_trajectory <- function(id) {
            dropdownButton(
                sliderInput(
                    inputId = paste0(id, \"_cellsize\"),
                    label = \"Point size\",
                    min = 0.1,
                    max = 5,
                    step = 0.05,
                    value = 0.1
                ),
                sliderInput(
                    inputId = paste0(id, \"_alpha\"),
                    label = \"Color alpha\",
                    min = 0.1,
                    max = 1,
                    step = 0.05,
                    value = 0.5
                ),
                sliderInput(
                    inputId = paste0(id, \"_edgesize\"),
                    label = \"Size of edges\",
                    min = 0.1,
                    max = 5,
                    step = 0.05,
                    value = 0.75
                ),
                shiny::sliderInput(
                    inputId = paste0(id, \"_textsize\"),
                    label = \"Text size\",
                    min = 0.5,
                    max = 50,
                    value = 20,
                    step = 0.5
                ),
                circle = TRUE,
                status = \"success\",
                size = \"sm\",
                icon = shiny::icon(\"cog\")
            )
        }

        gear_umaps <- function(id, narrow_options = TRUE) {
            shinyWidgets::dropdownButton(
                shiny::tagList(
                    shiny::sliderInput(
                        inputId = paste0(id, \"_pt_size\"),
                        label = \"Point size\",
                        min = 0.05, max = 5.00, value = 0.10, step = 0.05
                    ),
                    shiny::sliderInput(
                        inputId = paste0(id, \"_textsize\"),
                        label = \"Text size\",
                        min = 0.5,
                        max = 50,
                        value = 20,
                        step = 0.5
                    ),
                    if (narrow_options) {
                        shinyWidgets::prettySwitch(
                            inputId = paste0(id, \"_narrow\"),
                            label = \"Narrow down the area for points of interest\",
                            value = TRUE,
                            status = \"success\",
                            fill = TRUE
                        )
                    }
                ),
                circle = TRUE,
                status = \"success\",
                size = \"sm\",
                icon = shiny::icon(\"cog\")
            )
        }

        gear_download <- function(ns, id, label = \"\") {
            shinyWidgets::dropdownButton(
                label = \"\",
                icon = shiny::icon(\"download\"),
                status = \"success\",
                size = \"sm\",
                shiny::em(\"Note: Use one of the following extensions: PDF, PNG, SVG.\"),
                shiny::textInput(ns(paste0(\"filename_\", id)), \"File name:\", width = \"80%\", value = label),
                shiny::numericInput(ns(paste0(\"width_\", id)), \"Width (in):\", 7, 3, 100, 0.1),
                shiny::numericInput(ns(paste0(\"height_\", id)), \"Height (in):\", 7, 3, 100, 0.1),
                shiny::selectInput(ns(paste0(\"filetype_\", id)), \"Filetype\", choices = c(\"PDF\", \"PNG\", \"SVG\"), selected = \"PDF\", width = \"80%\"),
                shiny::downloadButton(ns(paste0(\"download_\", id)), label = \"Download Plot\")
            )
        }

        ##### UI #####
        ui_umaps_trajectory_expr <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::h2(\"Trajectory graph and distribution of cells expressing particular genes\"),
                shiny::splitLayout(
                    shiny::wellPanel(
                        shiny::splitLayout(
                            cellWidths = c(\"40px\", \"350px\"),
                            gear_trajectory(ns(\"traj1\")),
                            shiny::selectInput(
                                inputId = ns(\"metadata1\"),
                                label = \"Select metadata\",
                                choices = c(\"\"),
                                selected = NULL
                            )
                        ),
                        shiny::br(),
                        shiny::br(),
                        shiny::br(),
                        shiny::br(),
                        shiny::br(),
                        shiny::plotOutput(ns(\"graph_metadata\"), height = \"auto\")
                    ),
                    shiny::wellPanel(
                        shiny::selectizeInput(
                            inputId = ns(\"genes2\"),
                            label = \"Genes\",
                            choices = NULL,
                            multiple = TRUE
                        ),
                        shiny::splitLayout(
                            shiny::numericInput(
                                inputId = ns(\"expr_threshold2\"),
                                label = \"Gene expression threshold\",
                                min = 0, max = 10, value = 0, step = 0.01,
                                width = \"95%\"
                            ),
                            shiny::numericInput(
                                inputId = ns(\"relaxation2\"),
                                label = \"#genes not expressed\",
                                min = 0, max = 10, value = 0, step = 1,
                                width = \"95%\"
                            )
                        ),
                        gear_umaps(ns(\"umap2\"), TRUE),
                        shiny::plotOutput(ns(\"start_cells\"), height = \"auto\")
                    )
                )
            )
        }

        ui_rugplot <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::h2(\"Rugplots\"),
                shiny::splitLayout(
                    shiny::selectizeInput(
                        inputId = ns(\"genes3\"),
                        label = \"Gene to plot\",
                        choices = NULL
                    ),
                    shiny::selectInput(
                        inputId = ns(\"metadata3\"),
                        label = \"Color cells by\",
                        choices = c(\"\")
                    )
                ),
                shinyWidgets::dropdownButton(
                    label = \"\",
                    icon = shiny::icon(\"cog\"),
                    status = \"success\",
                    size = \"sm\",
                    shiny::sliderInput(
                        inputId = ns(\"rugplot_text\"),
                        label = \"Text size\",
                        min = 0.05,
                        max = 50,
                        value = 20,
                        step = 0.05
                    ),
                    shiny::sliderInput(
                        inputId = ns(\"rugplot_point\"),
                        label = \"Point size\",
                        min = 0.01,
                        max = 5,
                        value = 0.75,
                        step = 0.01
                    )
                ),
                shiny::plotOutput(ns(\"rugplot\"), height = \"auto\")
            )
        }

        ui_morans_table <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::h2(\"Moran's I test on all genes\"),
                shinyWidgets::dropdownButton(
                    label = \"\",
                    icon = shiny::icon(\"cog\"),
                    status = \"success\",
                    size = \"sm\",
                    shiny::sliderInput(
                        inputId = ns(\"moran_q\"),
                        label = \"Upper threshold for q_value\",
                        min = 0.00,
                        max = 1,
                        value = 0.05,
                        step = 0.001
                    ),
                    shiny::sliderInput(
                        inputId = ns(\"moran_i\"),
                        label = \"Lower threshold for Moran I\",
                        min = -1,
                        max = 1,
                        value = 0.1,
                        step = 0.01
                    )
                ),
                shiny::p(\"High values of I indicate that nearby cells have very similar values of the gene's expression.\"),
                shiny::p(\"Low q_value are a good indicator of good positive autocorrelation\"),
                DT::dataTableOutput(ns(\"moran_results\")),
                shiny::downloadButton(
                    outputId = ns(\"download_genes\"),
                    label = \"Download genes as CSV\"
                )
            )
        }

        ui_comparison_markers_panel <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::selectInput(
                    inputId = ns(\"select_k_markers\"),
                    label = \"Select the number of clusters (k) or metadata\",
                    choices = \"\"
                ),
                shinyWidgets::pickerInput(
                    inputId = ns(\"select_clusters_markers\"),
                    label = \"Select the groups of cells\",
                    choices = \"\",
                    inline = FALSE,
                    options = list(
                        `actions-box` = TRUE,
                        title = \"Select/deselect subgroups\",
                        size = 10,
                        width = \"90%\",
                        `selected-text-format` = \"count > 3\"
                    ),
                    multiple = TRUE
                )
            )
        }

        ui_comparison_markers <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::h2(\"Identification of markers\"),
                shiny::splitLayout(
                    cellWidths = c(\"90%\", \"10%\"),
                    shiny::plotOutput(ns(\"avg_expression_violin\"), height = \"auto\"),
                    shiny::tableOutput(ns(\"avg_expression_table\"))
                ),
                shinyWidgets::dropdownButton(
                    shiny::tagList(
                        shiny::sliderInput(
                            inputId = ns(\"logfc\"),
                            label = \"logFC threshold\",
                            min = 0.00, max = 10.00, value = 0.50, step = 0.01
                        ),
                        shiny::sliderInput(
                            inputId = ns(\"avg_expr_thresh\"),
                            label = \"Average expression threshold\",
                            min = 0.00, max = 0.01, value = 0.00
                        ),
                        shiny::sliderInput(
                            inputId = ns(\"min_pct\"),
                            label = \"Minimum gene frequency\",
                            min = 0.01, max = 1.00, value = 0.10, step = 0.01
                        ),
                        shiny::sliderInput(
                            inputId = ns(\"pval\"),
                            label = \"Maximum adj-pval\",
                            min = 0.001, max = 1.00, value = 0.01, step = 0.001
                        ),
                        shinyWidgets::prettySwitch(
                            inputId = ns(\"norm_type\"),
                            label = \"Data is normalised\",
                            value = TRUE,
                            status = \"success\",
                            fill = TRUE
                        )
                    ),
                    circle = TRUE,
                    status = \"success\",
                    size = \"sm\",
                    icon = shiny::icon(\"cog\")
                ),
                shiny::htmlOutput(ns(\"marker_text\")),
                shiny::fluidRow(
                    shiny::column(
                        6,
                        ui_comparison_markers_panel(ns(\"group_left\"))
                    ),
                    shiny::column(
                        6,
                        ui_comparison_markers_panel(ns(\"group_right\"))
                    )
                ),
                shiny::actionButton(ns(\"markers_button\"), \"Find markers!\", class = \"btn-danger\"),
                DT::dataTableOutput(ns(\"markers_dt\")),
                shiny::downloadButton(ns(\"markers_download_button\"), \"Download markers!\", class = \"btn-info\")
            )
        }

        ui_gene_clustering <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::h2(\"Gene clustering\"),
                shiny::p(\"Using a similar approach to cell clustering, the PhenoGraph pipeline is applied on the space of gene expression.\"),
                shiny::p(\"To assess the stability of the number of gene clusters (k), the user should provide:\"),
                shiny::HTML(\"
                <ul>
                    <li>a sequence of resolution values (start, stop, step) - this sequence will affect the range of possible values for k</li>
                    <li>a quality function - changing the default one will require noticeable changes on the resolution interval</li>
                    <li>the number of iteration of the gene clustering - increasing this will help in obtaining more qualitative results, but will slow the assessment; values between 5 and 10 are recommended</li>
                    <li>the number of neighbours for building the graph - values between 25 and 50 can be considered reasonable; increasing the number of neighbours will slow the assessment</li>
                    <li>the graph type (NN - unweighted oriented, SNN - weighted unoriented) - SNN is recommended</li>
                    <li>the embedding the graph is built on: PCA or UMAP - PCA is recommended</li>
                    <li>the number of different seeds the stability will be tested on - increasing the number will slow the assessment</li>
                </ul>
                \"),
                shiny::HTML(\"<em>Note</em>: <p>The clustering is done only on the genes filtered by the Moran's test (see the previous section).</p>\"),
                shiny::splitLayout(
                    cellWidths = \"20%\",
                    shiny::numericInput(
                        inputId = ns(\"res_start\"),
                        label = \"Resolution - start\",
                        value = 0.1,
                        min = 0
                    ),
                    shiny::numericInput(
                        inputId = ns(\"res_stop\"),
                        label = \"Resolution - stop\",
                        value = 1,
                        min = 0
                    ),
                    shiny::numericInput(
                        inputId = ns(\"res_step\"),
                        label = \"Resolution - step\",
                        value = 0.1,
                        min = 0
                    )
                ),
                shiny::splitLayout(
                    cellWidths = \"20%\",
                    shinyWidgets::radioGroupButtons(
                        inputId = ns(\"qfunc\"),
                        label = \"Quality function\",
                        choices = c(\"RBConfiguration\", \"RBER\", \"CPM\")
                    ),
                    shiny::numericInput(
                        inputId = ns(\"num_iter\"),
                        label = \"Number of iterations\",
                        value = 5,
                        min = 0
                    )
                ),
                shiny::splitLayout(
                    cellWidths = \"20%\",
                    shiny::numericInput(
                        inputId = ns(\"n_neighbours\"),
                        label = \"# neighbours\",
                        value = 25,
                        min = 5
                    ),
                    shinyWidgets::radioGroupButtons(
                        inputId = ns(\"graph_type\"),
                        label = \"Graph type\",
                        choices = c(\"SNN\", \"NN\")
                    ),
                    shinyWidgets::radioGroupButtons(
                        inputId = ns(\"graph_embedding\"),
                        label = \"Graph embedding\",
                        choices = c(\"PCA\", \"UMAP\")
                    )
                ),
                shiny::numericInput(
                    inputId = ns(\"n_repetitions\"),
                    label = \"Number of different seeds\",
                    value = 30,
                    min = 5
                ),
                shiny::actionButton(
                    inputId = ns(\"cluster_genes_button\"),
                    label = \"Cluster genes!\",
                    style = \"font-size:20px;\",
                    class = \"btn-danger\"
                ),
                shinyWidgets::dropdownButton(
                    shiny::checkboxGroupInput(
                        inputId = ns(\"k_select_groups\"),
                        label = \"Select groups\",
                        width = \"100%\",
                        choices = \"\"
                    ),
                    shiny::sliderInput(
                        inputId = ns(\"k_point_range\"),
                        label = \"Size point range\",
                        min = 0.1, max = 5.00, value = c(1.35, 2.35)
                    ),
                    shiny::sliderInput(
                        inputId = ns(\"k_distance_factor\"),
                        label = \"Distance between groups\",
                        min = 0.01, max = 1.00, value = 0.7
                    ),
                    shiny::sliderInput(
                        inputId = ns(\"k_ecc_thresh\"),
                        label = \"Low ECC threshold\",
                        min = 0.00, max = 1.00, value = 0.00, step = 0.005
                    ),
                    shiny::sliderInput(
                        inputId = ns(\"k_occ_thresh\"),
                        label = \"Low frequency threshold\",
                        min = 0.00, max = 1.00, value = 0.00
                    ),
                    shiny::sliderInput(
                        inputId = ns(\"k_nparts_thresh\"),
                        label = \"High # partitions threshold\",
                        min = 0.00, max = 1.00, value = 0.00
                    ),
                    shiny::sliderInput(
                        inputId = ns(\"k_text_size\"),
                        label = \"Text size\",
                        min = 0.50, max = 10.00, value = 1.5
                    ),
                    circle = TRUE,
                    status = \"success\",
                    size = \"sm\",
                    icon = shiny::icon(\"cog\")
                ),
                shiny::plotOutput(ns(\"clustering_assessment\"), height = \"auto\"),
                shiny::h2(\"Gene UMAP\"),
                shiny::splitLayout(
                    cellWidths = \"20%\",
                    shiny::selectInput(
                        inputId = ns(\"select_k_umap\"),
                        label = \"Select the number of gene clusters\",
                        choices = c(\"\")
                    ),
                    selectizeInput(
                        inputId = ns(\"search_gene\"),
                        label = \"Locate gene\",
                        choices = c(\"\")
                    )
                ),
                plotly::plotlyOutput(ns(\"gene_umap\"), height = \"auto\")
            )
        }

        ui_enrichment <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::h2(\"Enrichment analysis of gene clusters\"),
                shiny::splitLayout(
                    cellWidths = \"20%\",
                    shiny::selectInput(
                        inputId = ns(\"select_k_enrichment\"),
                        label = \"Select the number of gene clusters\",
                        choices = c(\"\")
                    ),
                    shiny::selectInput(
                        inputId = ns(\"module_id\"),
                        label = \"Select the gene cluster to perform enrichment analysis on\",
                        choices = c(\"\")
                    ),
                ),
                shiny::checkboxGroupInput(
                    inputId = ns(\"gprofilerSources\"),
                    label = \"Select data sources\",
                    choices = c(\"GO:BP\", \"GO:MF\", \"GO:CC\", \"KEGG\", \"REAC\", \"TF\", \"MIRNA\", \"CORUM\", \"HP\", \"HPA\", \"WP\"),
                    selected = c(\"GO:BP\", \"GO:MF\", \"GO:CC\", \"KEGG\", \"REAC\", \"TF\", \"MIRNA\")
                ),
                shiny::actionButton(
                    inputId = ns(\"enrichment_button\"),
                    label = \"Perform enrichment analysis!\",
                    style = \"font-size:20px;\",
                    class = \"btn-danger\"
                ),
                plotly::plotlyOutput(
                    outputId = ns(\"gost_plot\"),
                    height = \"auto\"
                ),
                DT::DTOutput(outputId = ns(\"gost_table\")),
                shiny::downloadButton(
                    outputId = ns(\"download_gost\"),
                    label = \"Download enriched terms as CSV\",
                    class = \"btn-info\"
                )
            )
        }

        ui_heatmap <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::h2(\"Gene clusters - cell clusters heatmap\"),
                shiny::p(\"The color indicates the scaled average expression of the genes across different cell groups.\"),
                shiny::splitLayout(
                    cellWidths = \"20%\",
                    shiny::selectInput(
                        inputId = ns(\"select_k_heatmap\"),
                        label = \"Select the number of gene clusters\",
                        choices = c(\"\")
                    ),
                    shiny::selectInput(
                        inputId = ns(\"metadata_heatmap\"),
                        label = \"Select the cell groups\",
                        choices = c(\"\")
                    ),
                ),
                plotly::plotlyOutput(ns(\"heatmap\"), height = \"auto\"),
                shiny::downloadButton(
                    outputId = ns(\"download_correlation\"),
                    label = \"Download heatmap as CSV\",
                    class = \"btn-info\"
                )
            )
        }

        ui_distribution_expr_panel <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::wellPanel(
                    shinyWidgets::dropdownButton(
                        label = \"\",
                        icon = shiny::icon(\"cog\"),
                        status = \"success\",
                        size = \"sm\",
                        shiny::sliderInput(
                            inputId = ns(\"distr_point\"),
                            label = \"Point size\",
                            min = 0.1,
                            max = 5,
                            value = 1,
                            step = 0.05
                        ),
                        shiny::sliderInput(
                            inputId = ns(\"distr_text\"),
                            label = \"Text_size\",
                            min = 0.5,
                            max = 50,
                            value = 20,
                            step = 0.5
                        )
                    ),
                    shiny::selectizeInput(
                        inputId = ns(\"gene4\"),
                        label = \"Select gene\",
                        choices = c(\"\")
                    ),
                    shiny::plotOutput(ns(\"distr_gene\"), height = \"auto\")
                )
            )
        }

        ui_distribution_expr <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::h2(\"Distribution of expression of genes\"),
                shiny::splitLayout(
                    ui_distribution_expr_panel(ns(\"left\")),
                    ui_distribution_expr_panel(ns(\"right\"))
                )
            )
        }

        ui_module_expr <- function(id) {
            ns <- shiny::NS(id)

            shiny::tagList(
                shiny::h2(\"Distribution of expression of gene clusters\"),
                shiny::selectInput(
                    inputId = ns(\"select_k_module_expr\"),
                    label = \"Select the number of gene clusters\",
                    choices = c(\"\")
                ),
                shinyWidgets::dropdownButton(
                    label = \"\",
                    icon = shiny::icon(\"cog\"),
                    status = \"success\",
                    size = \"sm\",
                    shiny::sliderInput(
                        inputId = ns(\"module_point\"),
                        label = \"Point size\",
                        min = 0.1,
                        max = 5,
                        value = 1,
                        step = 0.05
                    ),
                    shiny::sliderInput(
                        inputId = ns(\"module_text\"),
                        label = \"Text size\",
                        min = 0.5,
                        max = 50,
                        value = 20,
                        step = 0.5
                    )
                ),
                gear_download(ns, \"grid\", \"module_expression\"),
                shiny::plotOutput(ns(\"distr_module\"), height = \"auto\")
            )
        }

        ui <- fluidPage(
            tags$style(
                HTML(
                    \"#tabset_id {
                                    position: fixed;
                                    width: 100%;
                                    background-color: white;
                                    top: 0;
                                    z-index: 100;
                                    font-size: 25px;
                                    # margin-top: 72px;
                                    }\",
                    \"#graph_clust-show_config {
                                    z-index: 900;
                                    position: fixed;
                                    top: 10px;
                                    font-size: 15px;
                                    right: 25px;
                                }\",
                    \".shiny-split-layout > div {
                                    overflow: visible;
                                }\",
                )
            ),
            tags$head(tags$script(HTML('
                        var dimension = [0, 0];
                        var resizeId;
                        $(document).on(\"shiny:connected\", function(e) {
                            dimension[0] = window.innerWidth - 20;
                            dimension[1] = window.innerHeight - 30;
                            Shiny.onInputChange(\"dimension\", dimension);
                        });

                        function transferWindowSize() {
                            console.log(dimension);
                            console.log(window.innerHeight);

                            let dif_width = Math.abs(window.innerWidth - 20 - dimension[0]);
                            let dif_height = Math.abs(window.innerHeight - 30 - dimension[1]);
                            console.log(dif_height);

                            if (dif_width >= 200 || dif_height >= 200) {
                            console.log(\"Changed\")
                            dimension[0] = window.innerWidth - 20;
                            dimension[1] = window.innerHeight - 30;
                            Shiny.onInputChange(\"dimension\", dimension);
                            }
                        }

                        $(window).resize(function() {
                            clearTimeout(resizeId);
                            resizeId = setTimeout(transferWindowSize, 500);
                        });
                    '))),
            shinyjs::useShinyjs(),
            h1(\"Pseudotime analysis using Monocle3 - ", id, "\"),
            ui_umaps_trajectory_expr(\"trajectory_expression\"),
            ui_comparison_markers(\"markers\"),
            ui_rugplot(\"rugplot\"),
            ui_morans_table(\"morans_table\"),
            ui_gene_clustering(\"gene_clustering\"),
            ui_enrichment(\"enrichment\"),
            ui_heatmap(\"heatmap\"),
            ui_distribution_expr(\"distribution_expr\"),
            ui_module_expr(\"module_expr\")
        )

        ##### SERVER #####
        server_umaps_trajectory_expr <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    shiny::updateSelectInput(
                        session = session,
                        inputId = \"metadata1\",
                        choices = c(colnames(package_env$mon_obj@colData), \"pseudotime\")
                    )

                    shiny::updateSelectizeInput(
                        session = session,
                        inputId = \"genes2\",
                        server = TRUE,
                        choices = package_env$gene_list,
                        selected = package_env$gene_list[1],
                        options = list(
                            maxOptions = 7,
                            create = TRUE,
                            persist = TRUE
                        )
                    )

                    graph_plot_data <- shiny::reactive({
                        req(input$metadata1 != \"\")
                        shiny::req(package_env$dimension()[1] > 0)
                        shiny::req(package_env$dimension()[2] > 0)
                        p <- plot_cells(package_env$mon_obj,
                            trajectory_graph_color = \"black\",
                            alpha = input$traj1_alpha,
                            cell_size = input$traj1_cellsize,
                            trajectory_graph_segment_size = input$traj1_edgesize,
                            color_cells_by = input$metadata1,
                            label_cell_groups = FALSE,
                            graph_label_size = 0
                        ) + ggplot2::theme(text = ggplot2::element_text(size = input$traj1_textsize))
                        p_built <- ggplot_build(p)
                        p_built$data[[1]]$stroke <- NA

                        p_built
                    })
                    output$graph_metadata <- renderPlot(
                        {
                            as.ggplot(ggplot_gtable(graph_plot_data()))
                        },
                        width = function() {
                            shiny::req(package_env$plt_width(), package_env$plt_height())
                            max(min(package_env$plt_width_half(), package_env$plt_height()), 400)
                        },
                        height = function() {
                            shiny::req(package_env$plt_width(), package_env$plt_height())
                            max(min(package_env$plt_width_half(), package_env$plt_height()), 400)
                        }
                    )

                    start_mask <- shiny::reactive({
                        req(input$genes2 != \"\")
                        req(input$expr_threshold2)
                        if (length(input$genes2) == 0) {
                            color_info <- rep(FALSE, nrow(package_env$umap_df))
                        } else {
                            if (length(input$genes2) > 1) {
                                color_info <- DelayedMatrixStats::colSums2(package_env$mon_obj@assays@data$counts[input$genes2, ] > input$expr_threshold2) >= (length(input$genes2) - input$relaxation2)
                            } else {
                                color_info <- package_env$mon_obj@assays@data$counts[input$genes2, ] > input$expr_threshold2
                            }
                        }

                        if (input$\"umap2_narrow\") {
                            umap_emb <- package_env$umap_df[color_info, ]
                            if (nrow(umap_emb) == 0) {
                                break
                            }
                            gmedian <- Gmedian::Gmedian(umap_emb)
                            overall_distance <- ((umap_emb[, 1] - gmedian[1, 1])^2 + (umap_emb[, 2] - gmedian[1, 2])^2)^0.5
                            dist_stats <- fivenum(overall_distance)
                            color_info[seq_len(nrow(package_env$umap_df))[color_info][overall_distance > dist_stats[3]]] <- FALSE
                        }

                        color_info
                    })

                    output$start_cells <- renderPlot(
                        {
                            shiny::req(start_mask())
                            package_env$umap_df$color_val <- start_mask()

                            ggplot2::ggplot(rbind(package_env$umap_df[!package_env$umap_df$color_val, ], package_env$umap_df[package_env$umap_df$color_val, ]), aes(x = UMAP_1, y = UMAP_2, color = color_val)) +
                                ggplot2::geom_point(size = input$\"umap2_pt_size\") +
                                ggplot2::theme_bw() +
                                ggplot2::scale_color_manual(values = c(\"lightgray\", \"red\")) +
                                ggplot2::theme(
                                    legend.position = \"None\",
                                    text = ggplot2::element_text(size = input$\"umap2_textsize\")
                                )
                        },
                        width = function() {
                            shiny::req(package_env$plt_width(), package_env$plt_height())
                            max(min(package_env$plt_width_half(), package_env$plt_height()), 400)
                        },
                        height = function() {
                            shiny::req(package_env$plt_width(), package_env$plt_height())
                            max(min(package_env$plt_width_half(), package_env$plt_height()), 400)
                        }
                    )
                }
            )
        }

        server_rugplot <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    shiny::updateSelectInput(
                        session = session,
                        inputId = \"metadata3\",
                        choices = c(package_env$metadata_names, \"pseudotime\")
                    )

                    shiny::updateSelectizeInput(
                        session = session,
                        inputId = \"genes3\",
                        server = TRUE,
                        choices = package_env$gene_list,
                        selected = package_env$gene_list[1],
                        options = list(
                            maxOptions = 7,
                            create = TRUE,
                            persist = TRUE
                        )
                    )

                    output$rugplot <- shiny::renderPlot(
                        {
                            shiny::req(input$genes3 != \"\", input$metadata3 != \"\")
                            plot_genes_in_pseudotime(
                                cds_subset = package_env$mon_obj[input$genes3],
                                min_expr = 0,
                                cell_size = input$rugplot_point,
                                color_cells_by = input$metadata3
                            ) + theme(text = element_text(size = input$rugplot_text))
                        },
                        height = function() {
                            shiny::req(package_env$plt_width(), package_env$plt_height())
                            max(min(package_env$plt_width_half(), package_env$plt_height()), 400)
                        }
                    )
                }
            )
        }

        server_comparison_markers_panels <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    # also add categorical metadata
                    available_choices <- package_env$discrete_metadata

                    shiny::updateSelectInput(
                        session = session,
                        inputId = \"select_k_markers\",
                        choices = available_choices,
                        selected = available_choices[1]
                    )

                    shiny::observe({
                        shiny::req(input$select_k_markers %in% available_choices)

                        if (is.factor(package_env$mon_obj@colData[[input$select_k_markers]])) {
                            available_subgroups <- levels(package_env$mon_obj@colData[[input$select_k_markers]])
                        } else {
                            available_subgroups <- unique(package_env$mon_obj@colData[[input$select_k_markers]])
                        }

                        shinyWidgets::updatePickerInput(
                            session = session,
                            inputId = \"select_clusters_markers\",
                            choices = available_subgroups,
                            selected = available_subgroups[1]
                        )
                    }) %>% shiny::bindEvent(input$select_k_markers)
                }
            )
        }

        server_morans_table <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    shiny::observe({
                        shiny::req(package_env$graph_test_res$avg_expression[1])
                        package_env$filtered_graph_test_res(package_env$graph_test_res %>% filter(.data$morans_I > input$moran_i & .data$q_value < input$moran_q))
                    })

                    output$moran_results <- DT::renderDataTable(package_env$filtered_graph_test_res())

                    output$download_genes <- shiny::downloadHandler(
                        filename = function() {
                            \"genes.csv\"
                        },
                        content = function(file) {
                            write.csv(package_env$filtered_graph_test_res(), file)
                        }
                    )
                }
            )
        }

        server_comparison_markers <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    output$avg_expression_violin <- shiny::renderPlot(
                        {
                            vioplot::vioplot(
                                x = rhdf5::h5read(\"expression.h5\", \"average_expression\"),
                                horizontal = TRUE,
                                xlab = \"Average expression\",
                                main = \"Average gene expression\",
                                ylab = \"\",
                                xaxt = \"n\"
                            )
                        },
                        height = function() {
                            shiny::req(package_env$plt_height())
                            package_env$plt_height() / 2
                        }
                    )

                    avg_stats <- fivenum(rhdf5::h5read(\"expression.h5\", \"average_expression\"))

                    output$avg_expression_table <- shiny::renderTable(
                        {
                            data.frame(
                                row.names = c(\"min\", \"Q1\", \"median\", \"Q3\", \"max\"),
                                b = avg_stats
                            )
                        },
                        colnames = FALSE,
                        rownames = TRUE
                    )

                    shiny::updateSliderInput(
                        session = session,
                        inputId = \"avg_expr_thresh\",
                        min = round(avg_stats[1], digits = 3),
                        max = round(avg_stats[5], digits = 3),
                        step = 0.01
                    )

                    shinyjs::hide(\"markers_download_button\")
                    shinyjs::hide(\"markers_dt\")

                    server_comparison_markers_panels(\"group_left\")
                    server_comparison_markers_panels(\"group_right\")

                    markers_val <- shiny::reactive({
                        shiny::req(
                            input$\"group_left-select_clusters_markers\",
                            input$\"group_right-select_clusters_markers\"
                        )

                        shinyjs::disable(\"markers_button\")
                        subgroup_left <- input$\"group_left-select_k_markers\"
                        subgroup_right <- input$\"group_right-select_k_markers\"

                        mb1 <- package_env$mon_obj@colData[[subgroup_left]]
                        mb2 <- package_env$mon_obj@colData[[subgroup_right]]

                        cells_index_left <- which(mb1 %in% input$\"group_left-select_clusters_markers\")
                        cells_index_right <- which(mb2 %in% input$\"group_right-select_clusters_markers\")

                        markers_result <- ClustAssess::calculate_markers_shiny(
                            cells1 = cells_index_left,
                            cells2 = cells_index_right,
                            norm_method = ifelse(input$norm_type, \"LogNormalize\", \"\"),
                            used_slot = \"data\",
                            min_pct_threshold = input$min_pct,
                            logfc_threshold = input$logfc,
                            average_expression_threshold = input$avg_expr_thresh
                        )
                        if (nrow(markers_result) > 1) {
                            markers_result <- markers_result %>% dplyr::filter(.data$p_val_adj <= input$pval)
                        }
                        shinyjs::show(\"markers_dt\")
                        shinyjs::show(\"markers_download_button\")
                        shinyjs::enable(\"markers_button\")

                        return(markers_result)
                    }) %>% shiny::bindEvent(input$markers_button)


                    shiny::observe(
                        output$markers_dt <- DT::renderDataTable(
                            {
                                shiny::req(markers_val())
                                markers_val()
                            },
                            rownames = FALSE
                        )
                    ) %>% shiny::bindEvent(markers_val())

                    output$markers_download_button <- shiny::downloadHandler(
                        filename = function() {
                            \"markers.csv\"
                        },
                        content = function(file) {
                            write.csv(markers_val(), file)
                        }
                    )


                    shiny::observe({
                        shiny::req(markers_val())
                        shinyjs::show(\"markers_download_button\")
                    }) %>% shiny::bindEvent(markers_val())
                }
            )
        }

        server_gene_clustering <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    shiny::observe({
                        shinyjs::disable(\"cluster_genes_button\")
                        package_env$gene_module(my_cluster_genes(
                            mon_obj = package_env$mon_obj,
                            resolution = seq(from = input$res_start, to = input$res_stop, by = input$res_step),
                            quality_function = paste0(input$qfunc, \"VertexPartition\"),
                            num_iter = input$num_iter,
                            n_neighbours = input$n_neighbours,
                            graph_type = stringr::str_to_lower(input$graph_type),
                            graph_base_embedding = input$graph_embedding,
                            n_repetitions = input$n_repetitions,
                            keep_first = TRUE,
                            umap_arguments = list(
                                min_dist = 0.3,
                                n_neighbors = 30,
                                metric = \"cosine\",
                                init = \"spca\"
                            ),
                            subset_genes = rownames(package_env$filtered_graph_test_res())
                        ))

                        shiny::updateSliderInput(
                            session = session,
                            inputId = \"k_occ_thresh\",
                            max = max(package_env$gene_module()$summary$total_freq),
                            min = min(package_env$gene_module()$summary$total_freq) - 1,
                            step = 1
                        )

                        shiny::updateSliderInput(
                            session = session,
                            inputId = \"k_nparts_thresh\",
                            max = max(package_env$gene_module()$summary$n_partitions) + 1,
                            min = min(package_env$gene_module()$summary$n_partitions),
                            value = max(package_env$gene_module()$summary$n_partitions),
                            step = 1
                        )

                        shiny::updateCheckboxGroupInput(
                            session = session,
                            inputId = \"k_select_groups\",
                            choices = names(package_env$gene_module()$partitions),
                            selected = names(package_env$gene_module()$partitions)
                        )

                        shiny::updateSelectInput(
                            session = session,
                            inputId = \"select_k_umap\",
                            choices = setdiff(stringr::str_sort(names(package_env$gene_module()$partitions[[1]]), numeric = TRUE), \"1\")
                        )

                        shiny::updateSelectizeInput(
                            session = session,
                            inputId = \"search_gene\",
                            server = TRUE,
                            choices = rownames(package_env$gene_module()$gene_umap),
                            selected = rownames(package_env$gene_module()$gene_umap)[1],
                            options = list(
                                maxOptions = 7,
                                create = TRUE,
                                persist = TRUE
                            )
                        )

                        shinyjs::enable(\"cluster_genes_button\")
                    }) %>% shiny::bindEvent(input$cluster_genes_button)

                    shiny::observe({
                        temp_df <- package_env$filtered_graph_test_res()
                        temp_df$module <- package_env$gene_module()$partitions[[1]][[as.character(input$select_k_umap)]][[1]]$mb
                        package_env$filtered_graph_test_res(temp_df)
                    }) %>% shiny::bindEvent(input$select_k_umap)

                    output$clustering_assessment <- shiny::renderPlot(
                        {
                            shiny::req(input$k_nparts_thresh > 0, input$k_select_groups, \"summary\" %in% names(package_env$gene_module()))
                            shiny_plot_k_n_partitions(
                                summary_df = package_env$gene_module()$summary,
                                distance_factor = input$k_distance_factor,
                                filtered_cl_methods = input$k_select_groups,
                                text_size = input$k_text_size,
                                pt_size_range = input$k_point_range,
                                threshold_ecc = input$k_ecc_thresh,
                                threshold_occurences = input$k_occ_thresh,
                                threshold_nparts = input$k_nparts_thresh,
                                plt_height = package_env$plt_height(),
                                plot_title = \"Stability assessment of the number of gene clusters\",
                                display_legend = TRUE
                            )
                        },
                        height = function() {
                            shiny::req(package_env$dimension()[2] > 0)
                            package_env$plt_height()
                        }
                    )

                    output$gene_umap <- renderPlotly({
                        shiny::req(
                            \"summary\" %in% names(package_env$gene_module()),
                            input$search_gene != \"\",
                            input$select_k_umap %in% names(package_env$gene_module()$partitions[[1]])
                        )

                        a <- list(
                            x = package_env$gene_module()$gene_umap[input$search_gene, 1],
                            y = package_env$gene_module()$gene_umap[input$search_gene, 2],
                            text = input$search_gene,
                            xref = \"x\",
                            yref = \"y\",
                            showarrow = TRUE,
                            arrowhead = 7,
                            ax = 20,
                            ay = -40
                        )

                        plot_ly(
                            data = cbind(
                                package_env$gene_module()$gene_umap,
                                id = rownames(package_env$gene_module()$gene_umap),
                                module = factor(package_env$gene_module()$partitions[[1]][[as.character(input$select_k_umap)]][[1]]$mb),
                                avg_expression = package_env$graph_test_res[rownames(package_env$gene_module()$gene_umap), \"avg_expression\"]
                            ),
                            x = ~UMAP_1,
                            y = ~UMAP_2,
                            color = ~module,
                            text = ~ paste(\"Gene name:\", id, \"<br>Avg expression:\", avg_expression),
                            height = package_env$plt_height()
                        ) %>%
                            plotly::add_markers() %>%
                            plotly::layout(annotations = a)
                    })
                }
            )
        }

        server_enrichment <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    shiny::observe({
                        shiny::req(\"summary\" %in% names(package_env$gene_module()))
                        shiny::updateSelectInput(
                            session = session,
                            inputId = \"select_k_enrichment\",
                            choices = setdiff(stringr::str_sort(names(package_env$gene_module()$partitions[[1]]), numeric = TRUE), \"1\")
                        )
                    })

                    shiny::observe({
                        shiny::req(input$select_k_enrichment != \"\")

                        shiny::updateSelectInput(
                            session = session,
                            inputId = \"module_id\",
                            choices = as.character(seq_len(as.numeric(input$select_k_enrichment)))
                        )
                    }) %>% shiny::bindEvent(input$select_k_enrichment)
                    gprofiler_result <- shiny::reactiveVal(NULL)

                    shiny::observe({
                        shiny::req(input$module_id != \"\", \"summary\" %in% names(package_env$gene_module()))
                        shinyjs::disable(\"enrichment_button\")

                        gprf_res <- gprofiler2::gost(
                            query = rownames(package_env$gene_module()$gene_umap)[package_env$gene_module()$partitions[[1]][[input$select_k_enrichment]][[1]]$mb == as.numeric(input$module_id)],
                            sources = input$gprofilerSources,
                            organism = \"hsapiens\",
                            evcodes = TRUE
                        )

                        gprf_res$result$parents <- sapply(gprf_res$result$parents, toString)
                        gprofiler_result(gprf_res)

                        output$gost_table <- DT::renderDT({
                            gprofiler_result()$result[, seq_len(ncol(gprofiler_result()$result) - 2)]
                        })

                        output$gost_plot <- plotly::renderPlotly(
                            gprofiler2::gostplot(gprofiler_result())
                        )


                        shinyjs::show(\"download_gost\")
                        shinyjs::enable(\"enrichment_button\")
                    }) %>% shiny::bindEvent(input$enrichment_button)

                    output$download_gost <- shiny::downloadHandler(
                        filename = function() {
                            \"enrichment_results.csv\"
                        },
                        content = function(file) {
                            write.csv(gprofiler_result()$result, file)
                        }
                    )
                }
            )
        }

        server_heatmap <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    shiny::observe({
                        shiny::req(\"summary\" %in% names(package_env$gene_module()))
                        shiny::updateSelectInput(
                            session = session,
                            inputId = \"select_k_heatmap\",
                            choices = setdiff(stringr::str_sort(names(package_env$gene_module()$partitions[[1]]), numeric = TRUE), \"1\")
                        )
                    })

                    shiny::updateSelectInput(
                        session = session,
                        inputId = \"metadata_heatmap\",
                        choices = package_env$discrete_metadata
                    )

                    agg_mat <- shiny::reactive({
                        shiny::req(input$metadata_heatmap != \"\", input$select_k_heatmap != \"\")

                        cell_group_df <- tibble::tibble(
                            cell = row.names(package_env$mon_obj@colData),
                            cell_group = package_env$mon_obj@colData[[input$metadata_heatmap]]
                        )

                        agg_mat_res <- aggregate_gene_expression(
                            cds = package_env$mon_obj,
                            gene_group_df = data.frame(
                                id = rownames(package_env$gene_module()$gene_umap),
                                module = package_env$gene_module()$partitions[[1]][[input$select_k_heatmap]][[1]]$mb
                            ),
                            cell_group_df = cell_group_df
                        )

                        row.names(agg_mat_res) <- stringr::str_c(\"Gene cluster \", row.names(agg_mat_res))
                        colnames(agg_mat_res) <- stringr::str_c(\"Cluster \", colnames(agg_mat_res))

                        agg_mat_res
                    })

                    output$heatmap <- renderPlotly({
                        plot_ly(
                            z = agg_mat(),
                            y = row.names(agg_mat()),
                            x = colnames(agg_mat()),
                            type = \"heatmap\",
                            height = package_env$plt_height()
                        )
                    })

                    output$download_correlation <- shiny::downloadHandler(
                        filename = function() {
                            \"module_clusters_heatmap.csv\"
                        },
                        content = function(file) {
                            write.csv(agg_mat(), file)
                        }
                    )
                }
            )
        }

        server_distribution_expr_panel <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    shiny::observe({
                        shiny::req(package_env$graph_test_res)
                        shiny::updateSelectizeInput(
                            session = session,
                            inputId = \"gene4\",
                            server = TRUE,
                            choices = rownames(package_env$graph_test_res),
                            selected = rownames(package_env$graph_test_res)[1],
                            options = list(
                                maxOptions = 7,
                                create = TRUE,
                                persist = TRUE
                            )
                        )
                    })

                    output$distr_gene <- shiny::renderPlot(
                        {
                            shiny::req(input$gene4 != \"\")
                            plot_cells(
                                cds = package_env$mon_obj,
                                genes = input$gene4,
                                show_trajectory_graph = FALSE,
                                label_cell_groups = FALSE,
                                cell_size = input$distr_point,
                                label_leaves = FALSE
                            ) + theme(text = element_text(size = input$distr_text))
                        },
                        width = function() {
                            shiny::req(package_env$plt_width(), package_env$plt_height())
                            max(min(package_env$plt_width_half(), package_env$plt_height()), 400)
                        },
                        height = function() {
                            shiny::req(package_env$plt_width(), package_env$plt_height())
                            max(min(package_env$plt_width_half(), package_env$plt_height()), 400)
                        }
                    )
                }
            )
        }

        server_distribution_expr <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    server_distribution_expr_panel(\"left\")
                    server_distribution_expr_panel(\"right\")
                }
            )
        }

        server_module_expr <- function(id) {
            shiny::moduleServer(
                id,
                function(input, output, session) {
                    shiny::observe({
                        shiny::req(\"summary\" %in% names(package_env$gene_module()))
                        shiny::updateSelectInput(
                            session = session,
                            inputId = \"select_k_module_expr\",
                            choices = setdiff(stringr::str_sort(names(package_env$gene_module()$partitions[[1]]), numeric = TRUE), \"1\")
                        )
                    })

                    output$distr_module <- shiny::renderPlot(
                        {
                            shiny::req(input$select_k_module_expr != \"\", \"summary\" %in% names(package_env$gene_module()))
                            plot_cells(
                                cds = package_env$mon_obj,
                                genes = data.frame(
                                    id = rownames(package_env$gene_module()$gene_umap),
                                    module = package_env$gene_module()$partitions[[1]][[input$select_k_module_expr]][[1]]$mb
                                ),
                                show_trajectory_graph = FALSE,
                                label_cell_groups = FALSE,
                                label_leaves = FALSE,
                                cell_size = input$module_point
                            ) + ggplot2::theme(text = ggplot2::element_text(size = input$module_text))
                        },
                        height = function() {
                            shiny::req(package_env$plt_width(), package_env$plt_height())
                            package_env$plt_height()
                        }
                    )

                    output$download_grid <- shiny::downloadHandler(
                        filename = function() {
                            paste(input$filename_grid, tolower(input$filetype_grid), sep = \".\")
                        },
                        content = function(file) {
                            shiny::req(input$select_k_module_expr != \"\", \"summary\" %in% names(package_env$gene_module()))
                            ggplot2::ggsave(
                                filename = file,
                                plot = plot_cells(
                                    cds = package_env$mon_obj,
                                    genes = data.frame(
                                        id = rownames(package_env$gene_module()$gene_umap),
                                        module = package_env$gene_module()$partitions[[1]][[input$select_k_module_expr]][[1]]$mb
                                    ),
                                    show_trajectory_graph = FALSE,
                                    label_cell_groups = FALSE,
                                    label_leaves = FALSE,
                                    cell_size = input$module_point
                                ) + theme(text = element_text(size = input$module_text)),
                                width = input$width_grid,
                                height = input$height_grid
                            )
                        }
                    )
                }
            )
        }

        server <- function(input, output, session) {
            shinyjs::hide(\"download_gost\")
            shinyjs::hide(\"calculate_pseudotime\")

            init_env(shiny::reactive(input$dimension))
            server_comparison_markers(\"markers\")
            server_umaps_trajectory_expr(\"trajectory_expression\")
            server_rugplot(\"rugplot\")
            server_morans_table(\"morans_table\")
            server_gene_clustering(\"gene_clustering\")
            server_enrichment(\"enrichment\")
            server_heatmap(\"heatmap\")
            server_distribution_expr(\"distribution_expr\")
            server_module_expr(\"module_expr\")
        }

        shinyApp(ui, server)
    ")
    write(file_content, file.path(output_dir, "app.R"))
}



# target_folder_first <- "/servers/sutherland-scratch/andi/projects/0_2208_Floris/output/pseudotime_shiny_apps/first_part_analysis"
# target_folder_second <- "/servers/sutherland-scratch/andi/projects/0_2208_Floris/output/pseudotime_shiny_apps/second_part_analysis"
# samples_metadata_path <- "/servers/sutherland-scratch/andi/projects/0_2208_Floris/metadata/pseudotime_genes.csv"
# project_folder <- "/servers/sutherland-scratch/andi/projects/0_2208_Floris"
# samples_metadata <- read.csv(samples_metadata_path, comment.char = "#")
# seurat_ca_metadata_path <- "/servers/sutherland-scratch/andi/projects/0_2208_Floris/metadata/seurat_paths_to_clustasess.csv"
# seurat_ca_metadata <- read.csv(seurat_ca_metadata_path, comment.char = "#")
# se_prefix <- "aggr_filtered_coding_highrp-"

# for (i in seq_along(samples_metadata[, 1])) {
#     id <- samples_metadata$SampleID[i]
#     found_corresp <- FALSE

#     for (j in seq_len(nrow(seurat_ca_metadata))) {
#         seurat_path <- seurat_ca_metadata$seurat[j]

#         if (basename(seurat_path) == glue::glue("{id}.rds") || basename(seurat_path) == glue::glue("{se_prefix}{id}.rds")) {
#             found_corresp <- TRUE
#             print(paste(id, seurat_path))
#             basedir_seurat <- basename(dirname(seurat_path))
#             if (basedir_seurat == "HV") {
#                 target_dir_first <- file.path(target_folder_first, id)
#                 target_dir_second <- file.path(target_folder_second, id)
#             } else {
#                 target_dir_first <- file.path(target_folder_first, basedir_seurat, id)
#                 target_dir_second <- file.path(target_folder_second, basedir_seurat, id)
#             }
#             break
#         }
#     }

#     if (!found_corresp) {
#         next
#     }
    
#     psd_id <- samples_metadata$PseudotimeID[i]
#     start_genes <- strsplit(samples_metadata$StartGenes[i], ";")[[1]]
#     end_genes <- strsplit(samples_metadata$EndGenes[i], ";")[[1]]

#     start_thresh <- samples_metadata$ThresholdStart[i]
#     end_thresh <- samples_metadata$ThresholdEnd[i]

#     start_relax <- samples_metadata$RelaxStart[i]
#     end_relax <- samples_metadata$RelaxEnd[i]


#     print(paste(id, psd_id))
#     write_object(
#         mon_obj = readRDS(file.path(target_dir_first, "monocle_object.rds")),
#         trajectory_id = paste(id, psd_id, sep = "-"),
#         start_genes = start_genes,
#         end_genes = end_genes,
#         start_expression_thresh = start_thresh,
#         end_expression_thresh = end_thresh,
#         start_relax_ngenes = start_relax,
#         end_relax_ngenes = end_relax,
#         ncores = 50,
#         output_dir = file.path(target_dir_second, psd_id)
#     )
# }


# mdn <- readRDS(file.path(project_folder, "output", "pseudotime_shiny_apps", "first_part_analysis", id, "monocle_object.rds"))

# write_object.data_frame(ids_df = read.csv("/servers/sutherland-scratch/andi/projects/0_2304_Floris_snrna/metadata/pseudotime_genes.csv"), seurat_ca_corresp = read.csv("/servers/sutherland-scratch/andi/projects/0_2304_Floris_snrna/metadata/seurat_paths_to_clustassess.csv"), target_app_dir_first_part = "/servers/sutherland-scratch/andi/projects/0_2304_Floris_snrna/output/pseudotime_apps/first_part_analysis", target_app_dir_second_part = "/servers/sutherland-scratch/andi/projects/0_2304_Floris_snrna/output/pseudotime_apps/second_part_analysis")
write_object.data_frame <- function(ids_df,
                                    seurat_ca_corresp,
                                    target_app_dir_second_part,
                                    target_app_dir_first_part,
                                    ncores = 1,
                                    learn_graph_controls = NULL,
                                    use_closed_loops = FALSE,
                                    use_partitions = FALSE,
                                    nodes_per_log10_cells = 30,
                                    prefix = "",
                                    suffix = "") {
    common_seurat_dir <- find_common_directory(seurat_ca_corresp$seurat)
    for (i in seq_along(ids_df[, 1])) {
        id <- ids_df$SampleID[i]
        associated_index <- find_corresp_index(id, seurat_ca_corresp, prefix, suffix)

        if (associated_index == -1) {
            next
        }

        target_folder_first <- target_app_dir_first_part
        target_folder_second <- target_app_dir_second_part
        seurat_path <- dirname(seurat_ca_corresp$seurat[associated_index])
        # to_add_folders <- c(id)
        to_add_folders <- c(id, ifelse(use_closed_loops, "loop", "no_loop"), ids_df$Ftype[i], ids_df$Fsize[i])

        while (basename(seurat_path) != common_seurat_dir) {
            to_add_folders <- c(basename(seurat_path), to_add_folders)
            seurat_path <- dirname(seurat_path)
        }

        to_add_folders <- as.list(to_add_folders)

        for (j in seq_along(to_add_folders)) {
            target_folder_first <- file.path(target_folder_first, to_add_folders[[j]])
            target_folder_second <- file.path(target_folder_second, to_add_folders[[j]])
        }

        dir.create(target_folder_second, showWarnings = FALSE, recursive = TRUE)

        psd_id <- ids_df$PseudotimeID[i]
        start_genes <- strsplit(ids_df$StartGenes[i], ";")[[1]]
        end_genes <- strsplit(ids_df$EndGenes[i], ";")[[1]]

        start_thresh <- ids_df$ThresholdStart[i]
        end_thresh <- ids_df$ThresholdEnd[i]

        start_relax <- ids_df$RelaxStart[i]
        end_relax <- ids_df$RelaxEnd[i]

        print(glue::glue("Creating app for {id} - {psd_id}"))

        write_object(
            mon_obj = readRDS(file.path(target_folder_first, "monocle_object.rds")),
            trajectory_id = paste(id, psd_id, sep = "-"),
            start_genes = start_genes,
            end_genes = end_genes,
            start_expression_thresh = start_thresh,
            end_expression_thresh = end_thresh,
            start_relax_ngenes = start_relax,
            end_relax_ngenes = end_relax,
            ncores = ncores,
            learn_graph_controls = learn_graph_controls,
            use_closed_loops = use_closed_loops,
            use_partitions = use_partitions,
            nodes_per_log10_cells = nodes_per_log10_cells,
            output_dir = file.path(target_folder_second, psd_id)
        )
    }
}
