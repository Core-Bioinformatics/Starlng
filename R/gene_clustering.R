#' @importFrom foreach %dopar% %do% %:%
#' @importFrom doFuture %dofuture%

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

    return(ClustAssess::get_highest_prune_param(
        nn_adj_matrix$nn,
        k
    )$adj_matrix)
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

community_detection_worker <- function(shared_adj_object,
                                       graph_type,
                                       resolution,
                                       quality_function,
                                       number_interations,
                                       seed) {
    g <- preproc_shared_adj_obj(shared_adj_object, graph_type)

    mb <- leidenbase::leiden_find_partition(
        igraph = g,
        edge_weights = igraph::E(g)$weight,
        resolution_parameter = resolution,
        partition_type = quality_function,
        num_iter = number_interations,
        seed = seed,
        verbose = FALSE
    )
    mb$membership
}

community_detection_master <- function(adj_object,
                                       resolutions,
                                       quality_functions,
                                       number_interations,
                                       seeds,
                                       graph_type = "snn",
                                       memory_log_file = NULL) {
    ncores <- foreach::getDoParWorkers()
    if (ncores > 1) {
        shared_adj_object <- SharedObject::share(adj_object)
    } else {
        shared_adj_object <- build_ig_from_adj(adj_object, graph_type)
    }

    used_functions <- c("community_detection_worker")
    used_objects <- c("shared_adj_object", "graph_type", "seeds", "resolutions", "quality_functions")

    # NOTE I've used indices instead of actual values to make the transition easier
    # in case doRNG might need to be used; for now, the results seem to be consistent
    # TODO investigate the use of dofuture, whether it will help optimising the code
    clusters_list <- suppressWarnings(foreach::foreach(
        i_qfunc = seq_along(quality_functions),
        .inorder = TRUE,
        .final = function(x) setNames(x, quality_functions)
    ) %:%
    foreach::foreach(
        i_res = seq_along(resolutions),
        .inorder = TRUE,
        .final = function(x) setNames(x, as.character(resolutions))
    ) %:%
    foreach::foreach(
        i_seed = seq_along(seeds),
        .export = c(used_objects, used_functions),
        .packages = c("Starlng", "SharedObject"),
        .inorder = FALSE
    ) %dopar% {
        worker_res <- community_detection_worker(
            shared_adj_object = shared_adj_object,
            graph_type = graph_type,
            resolution = resolutions[[i_res]],
            quality_function = quality_functions[[i_qfunc]],
            number_interations = number_interations,
            seed = seeds[[i_seed]]
        )

        worker_res
    })

    if (ncores > 1) {
        shared_adj_object <- SharedObject::unshare(shared_adj_object)
    }

    return(clusters_list)
}

clustering_pipeline <- function(embedding,
                                n_neighbours = seq(from = 5, to = 50, by = 5),
                                graph_type = "snn",
                                prune_value = -1,
                                resolutions = seq(from = 0.1, to = 1, by = 0.1),
                                quality_functions = c("RBConfigurationVertexPartition"),
                                number_interations = 5,
                                seeds = NULL,
                                number_repetitions = 30,
                                memory_log_file = NULL) {
    # TODO add progress bar
    # TODO check if it's worth to add the option of memory logging
    # TODO check if the nested foreach don't add overhead and if there are other better ways to do it
    if (is.null(seeds)) {
        seeds <- seq(from = 1, by = 100, length.out = number_repetitions)
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
            quality_functions = quality_functions,
            number_interations = number_interations,
            seeds = seeds,
            graph_type = graph_type,
            memory_log_file = memory_log_file
        )
        clear_psock_memory()
    }

    return(clusters_list)
}
