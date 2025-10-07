#' @importFrom foreach %dopar% %do% %:%
NULL

#' Parallel NN indexing
#' 
#' @description This function provides a parallel approach of finding the
#' nearest neighbours of each points on a given embedding. The parallelisation
#' is done by splitting the data into chunks and running the `RANN::nn2`
#' function on each chunk in parallel.
#' 
#' @param embedding The embedding matrix where each row represents a point
#' and each column a dimension.
#' @param k The number of nearest neighbours to be found.
#' 
#' @return A matrix where each row contains the indices of the k-nearest
#' neighbours of the corresponding point.
#' @export
parallel_nn2_idx <- function(embedding, k) {
    ncores <- foreach::getDoParWorkers()

    if (ncores == 1) {
        return(RANN::nn2(embedding, k = k)$nn.idx)
    }

    chunk_size <- nrow(embedding) %/% ncores
    chunks <- split(seq_len(chunk_size * ncores), rep(seq_len(ncores), each = chunk_size))

    if (chunk_size * ncores != nrow(embedding)) {
        chunks[[ncores]] <- c(chunks[[ncores]], (chunk_size * ncores + 1):nrow(embedding))
    }

    shared_emb <- SharedObject::share(embedding)
    i <- NULL

    # HACK check if there is another way to not get the `already exported variable` warning
    nn_result <- suppressWarnings(foreach::foreach (
        i = seq_len(ncores),
        .inorder = TRUE,
        .combine = rbind,
        .export = c("shared_emb", "chunks", "k")
    ) %dopar% {
        RANN::nn2(
            data = shared_emb,
            query = shared_emb[chunks[[i]], , drop = FALSE],
            k = k,
            eps = 0
        )$nn.idx
    })

    shared_emb <- SharedObject::unshare(shared_emb)

    nn_result
}

build_adj_from_idx <- function(nn_idx, k, graph_type = "snn", prune_value = -1) {
    nn_adj_matrix <- ClustAssess::getNNmatrix(nn_idx, k, 0, prune_value)

    if (graph_type == "nn") {
        return(nn_adj_matrix$nn)
    }

    if (prune_value >= 0) {
        return(nn_adj_matrix$snn)
    }

    prune_result <- ClustAssess::get_highest_prune_param(
        nn_adj_matrix$nn,
        k
    )

    return(prune_result$adj_matrix)
}

build_ig_from_adj <- function(adj_matrix, graph_type = "snn") {
    if (graph_type == "snn") {
        return(igraph::graph_from_adjacency_matrix(
            adjmatrix = adj_matrix,
            mode = "undirected",
            weighted = TRUE
        ))
    }

    return(igraph::graph_from_adjacency_matrix(
        adjmatrix = adj_matrix,
        mode = "directed"
    ))
}

preproc_shared_adj_obj <- function(shared_adj_object, graph_type) {
    if (igraph::is_igraph(shared_adj_object)) {
        return(shared_adj_object)
    }

    # if it's parallel, the object is an adjacency matrix
    # that should be converted to an igraph
    build_ig_from_adj(shared_adj_object, graph_type)
}

get_workers_res_index <- function(resolution_list, nseeds) {
    first_df <- TRUE
    seed_index <- seq_len(nseeds)
    for (i_res in seq_along(resolution_list)) {
        res_index <- seq_along(resolution_list[[i_res]])
        current_df <- expand.grid(seed_index, res_index)
        colnames(current_df) <- c("seed", "res")
        current_df$quality_function <- i_res

        if (first_df) {
            workers_res_index <- current_df
            first_df <- FALSE
        } else {
            workers_res_index <- rbind(workers_res_index, current_df)
        }
    }

    workers_res_index <- workers_res_index[, c(3, 2, 1)]

    return(workers_res_index)
}

create_hierarchy_from_master <- function(master_result, workers_res_index, resolution_list) {
    hier_list <- list()

    for (qfunc in names(resolution_list)) {
        hier_list[[qfunc]] <- list()
        for (res in as.character(resolution_list[[qfunc]])) {
            hier_list[[qfunc]][[res]] <- list()
        }
    }

    for (i in seq_len(nrow(workers_res_index))) {
        qfunc <- names(resolution_list)[workers_res_index$quality_function[i]]
        res <- as.character(resolution_list[[qfunc]][workers_res_index$res[i]])
        seed <- workers_res_index$seed[i]

        hier_list[[qfunc]][[res]][[seed]] <- master_result[[i]]
    }

    return(hier_list)
}


community_detection_worker <- function(shared_adj_object,
                                       graph_type,
                                       resolution,
                                       quality_function,
                                       number_iterations,
                                       seed) {
    g <- preproc_shared_adj_obj(shared_adj_object, graph_type)

    mb <- leidenbase::leiden_find_partition(
        igraph = g,
        edge_weights = igraph::E(g)$weight,
        resolution_parameter = resolution,
        partition_type = quality_function,
        num_iter = number_iterations,
        seed = seed,
        verbose = FALSE
    )
    mb$membership
}

#' Community detection wrapper
#' 
#' @description This function is a simplified version of the `clustering_pipeline`.
#' It starts with an already-generated adjacency matrix and applied the leiden
#' algorithm to identify the clusters. The details of how the algorithm is run
#' are the same as in the `clustering_pipeline` function.
#' 
#' @param adj_object The adjacency matrix to be used. The matrix should be
#' square and should have defined rownames.
#' @param resolutions Check the `clustering_pipeline` documentation.
#' @param number_iterations Check the `clustering_pipeline` documentation.
#' @param seeds Check the `clustering_pipeline` documentation.
#' @param graph_type Check the `clustering_pipeline` documentation.
#' @param memory_log_file Check the `clustering_pipeline` documentation.
#' 
#' @return A list of lists containing the clusters found for each combination
#' of number of neighbours and quality functions.
#' @export
community_detection_master <- function(adj_object,
                                       resolutions,
                                       number_iterations,
                                       seeds,
                                       graph_type = "snn",
                                       memory_log_file = NULL) {
    ncores <- foreach::getDoParWorkers()
    if (ncores > 1) {
        shared_adj_object <- SharedObject::share(adj_object)
    } else {
        shared_adj_object <- build_ig_from_adj(adj_object, graph_type)
    }

    quality_functions <- names(resolutions)
    is_empty <- sapply(resolutions, is.null)
    quality_functions <- quality_functions[!is_empty]
    resolutions <- resolutions[!is_empty]

    index <- NULL
    workers_indices <- get_workers_res_index(resolutions, length(seeds))

    used_functions <- c("community_detection_worker")
    used_objects <- c("shared_adj_object", "graph_type", "seeds", "resolutions", "quality_functions")

    # NOTE I've used indices instead of actual values to make the transition easier
    # in case doRNG might need to be used; for now, the results seem to be consistent
    # TODO investigate the use of dofuture, whether it will help optimising the code
    clusters_list <- foreach::foreach(
        index = seq_len(nrow(workers_indices)),
        .inorder = TRUE,
        .export = c(used_objects, used_functions),
        .packages = c("Starlng", "SharedObject")
    ) %dopar% {
        qfunc <- quality_functions[workers_indices$quality_function[index]]
        res <- as.numeric(resolutions[[qfunc]][workers_indices$res[index]])
        seed <- seeds[workers_indices$seed[index]]

        worker_res <- community_detection_worker(
            shared_adj_object = shared_adj_object,
            graph_type = graph_type,
            resolution = res,
            quality_function = qfunc,
            number_iterations = number_iterations,
            seed = seed
        )

        worker_res
    }

    if (ncores > 1) {
        shared_adj_object <- SharedObject::unshare(shared_adj_object)
    }

    return(create_hierarchy_from_master(
        master_result = clusters_list,
        workers_res_index = workers_indices,
        resolution_list = resolutions
    ))
}

#' Community detection pipeline
#' 
#' @description This function applies the PhenoGraph pipeline given an
#' embedding matrix. The pipeline creates the exact nearest-neighbour graph
#' using the RANN package, then applies the Leiden algorithm to find the
#' clusters. The pipeline can be run multiple times while changing the random
#' seed, to ensure the robustness of the resulting clusters.
#' 
#' @param embedding The embedding matrix where each row represents a point
#' and each column a dimension.
#' @param n_neighbours A vector of integers representing the number of number
#' of neighbours to be used to build the adjacency graph.
#' @param graph_type The type of graph. Can be either "nn" - directed unweighted
#' graph or "snn" - undirected weighted graph, where the weights are calculated
#' using the JSI score.
#' @param prune_value The value used to prune the edges of the graph. If the
#' value is negative, the function will determine the highest pruning value
#' that will keep the graph connected. Defaults to -1.
#' @param resolutions A list of vectors representing the resolution parameters
#' used in the community detection algorithm. Each name of the list is
#' associated with a quality function. The values of the list are the resolution
#' parameters to be used. The available choices in terms of quality functions
#' are "RBConfigurationVertexPartition", "RBERVertexPartition" and
#' "ModularityVertexPartition".
#' @param number_iterations How many iterations the Leiden algorithm should
#' run. Defaults to 5.
#' @param seeds A vector of integers representing the random seeds to be used
#' in the community detection algorithm. If NULL, the function will generate
#' the seeds based on the number of repetitions. Defaults to NULL.
#' @param number_repetitions The number of repetitions the community detection
#' pipeline should run. Defaults to 100.
#' @param merge_identical_partitions Logical indicating if the function should
#' merge the partitions that are identical and group them by the number of
#' clusters. Defaults to TRUE. 
#' @param memory_log_file The path to the file where the memory usage of the
#' function should be logged. Defaults to NULL.
#' 
#' @return A list of lists containing the clusters found for each combination
#' of number of neighbours and quality functions. If `merge_identical_partitions`
#' is set to TRUE, the last level will contain the partitions grouped by the
#' number of clusters, alongside with their overall Element-Centric Consistency
#' (ECC) score.
#' @export
clustering_pipeline <- function(embedding,
                                n_neighbours = seq(from = 5, to = 50, by = 5),
                                graph_type = "snn",
                                prune_value = -1,
                                resolutions = list(
                                    "RBConfigurationVertexPartition" = seq(from = 0.1, to = 2, by = 0.1),
                                    "RBERVertexPartition" = NULL,
                                    "ModularityVertexPartition" = NULL
                                ),
                                number_iterations = 5,
                                seeds = NULL,
                                number_repetitions = 100,
                                merge_identical_partitions = TRUE,
                                memory_log_file = NULL) {
    # TODO add progress bar
    # TODO check if it's worth to add the option of memory logging
    # TODO check if the nested foreach don't add overhead and if there are other better ways to do it
    if (is.null(seeds)) {
        seeds <- seq(from = 1, by = 100, length.out = number_repetitions)
    }

    if (length(seeds) < number_repetitions) {
        seeds <- c(seeds, seq(from = seeds[length(seeds)] + 1, by = 100, length.out = number_repetitions - length(seeds)))
    }

    point_names <- rownames(embedding)
    n_neighbours <- sort(unique(n_neighbours), decreasing = TRUE)
    nn_idx <- parallel_nn2_idx(embedding, k = n_neighbours[1])

    clusters_list <- list()
    for (n_neigh in n_neighbours) {
        # TODO check if you can build the adj matrix based on the previous one by deelting some neighbours
        nn_adj_matrix <- build_adj_from_idx(nn_idx, n_neigh, graph_type, prune_value)
        gc()
        rownames(nn_adj_matrix) <- point_names
        colnames(nn_adj_matrix) <- point_names

        clusters_list[[as.character(n_neigh)]] <- community_detection_master(
            adj_object = nn_adj_matrix,
            resolutions = resolutions,
            number_iterations = number_iterations,
            seeds = seeds,
            graph_type = graph_type,
            memory_log_file = memory_log_file
        )
        # clear_psock_memory()
    }

    if (!merge_identical_partitions) {
        return(clusters_list)
    }

    clusters_list <- group_by_clusters_general(clusters_list, 3)
    return(get_clusters_consistency(clusters_list))
}
