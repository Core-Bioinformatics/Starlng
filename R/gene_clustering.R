
parallel_nn2_idx <- function(embedding, k) {
    ncores <- foreach::getDoParWorkers()

    if (ncores == 1) {
        return(RANN::nn2(embedding, k = k)$nn.idx)
    }

    chunk_size <- nrow(embedding) %/% ncores
    chunks <- suppressWarnings(
        split(seq_len(nrow(embedding)), rep(seq_len(ncores), each = chunk_size))
    )

    shared_emb <- SharedObject::share(embedding)

    nn_result <- foreach(i = seq_len(ncores)) %dopar% {
        RANN::nn2(
            data = shared_emb,
            query = shared_emb[chunks[[i]], , drop = FALSE],
            k = k,
            eps = 0
        )$nn.idx
    }

    shared_emb <- SharedObject::unshare(shared_emb)

    return(do.call(rbind, nn_result))
}

process_shared_adj_object_clustering <- function(shared_adj_object, graph_type) {
    ncores <- foreach::getDoParWorkers()
    if (ncores == 1) {
        return(shared_adj_object)
    }

    # if it's parallel, the object is an adjacency matrix
    # that should be converted to an igraph
    if (graph_type == "snn") {
        return(igraph::graph_from_adjacency_matrix(
            adjmatrix = shared_adj_object,
            mode = "undirected",
            weighted = TRUE
        ))
    }

    return(igraph::graph_from_adjacency_matrix(
        adjmatrix = shared_adj_object,
        mode = "directed"
    ))
}

community_detection_worker <- function(shared_adj_object,
                                       graph_type,
                                       resolution,
                                       quality_function,
                                       number_interations,
                                       seed) {
    g <- process_shared_adj_object_clustering(shared_adj_object, graph_type)

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
                                       vars_to_exclude) {
    ncores <- foreach::getDoParWorkers()
    if (ncores > 1) {
        shared_adj_object <- SharedObject::share(adj_object)
    } else {
        if (graph_type == "snn") {
            shared_adj_object <- igraph::graph_from_adjacency_matrix(
                adjmatrix = adj_object,
                mode = "undirected",
                weighted = TRUE
            )
        } else {
            shared_adj_object <- igraph::graph_from_adjacency_matrix(
                adjmatrix = adj_object,
                mode = "directed"
            )
        }
    }

    clusters_list <- foreach::foreach(
        qfunc = quality_functions,
        .inorder = TRUE,
        .final = function(x) setNames(x, quality_functions)
    ) %:%
    foreach::foreach(
        res = resolutions,
        .inorder = TRUE,
        .final = function(x) setNames(x, as.character(resolutions))
    ) %:%
    foreach::foreach(
        seed = seeds,
        .inorder = FALSE
    ) %do% {
        community_detection_worker(
            shared_adj_object = shared_adj_object,
            graph_type = graph_type,
            resolution = res,
            quality_function = qfunc,
            number_interations = number_interations,
            seed = seed
        )
    }

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
                                number_repetitions = 30) {
    if (is.null(seeds)) {
        seeds <- seq(from = 1, by = 100, length.out = number_repetitions)
    }

    point_names <- rownames(embedding)
    n_neighbours <- sort(unique(n_neighbours), decreasing = TRUE)
    nn_idx <- parallel_nn2_idx(embedding, k = n_neighbours[1])
    # clear_psock_memory()

    clusters_list <- list()
    for (n_neigh in n_neighbours) {
        nn_adj_matrix <- ClustAssess::getNNmatrix(nn_idx, n_neigh, 0, prune_value)

        if (graph_type == "snn") {
            if (prune_value == -1) {
                nn_adj_matrix <- ClustAssess::get_highest_prune_param(
                    nn_adj_matrix$nn,
                    n_neigh
                )$adj_matrix
            } else {
                nn_adj_matrix <- nn_adj_matrix$snn
            }
        } else {
            nn_adj_matrix <- nn_adj_matrix$nn
        }
        gc()
        rownames(nn_adj_matrix) <- point_names
        colnames(nn_adj_matrix) <- point_names

        clusters_list[[as.character(n_neigh)]] <- community_detection_master(
            adj_object = nn_adj_matrix,
            resolutions = resolutions,
            quality_functions = quality_functions,
            number_interations = number_interations,
            seeds = seeds,
            graph_type = graph_type
        )
    }

    return(clusters_list)
}
