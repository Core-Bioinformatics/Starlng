update_mononcle_partition <- function(mon_obj, new_partition) {
    new_partition <- factor(as.numeric(factor(new_partition)))

    if (!("partitions" %in% names(mon_obj@clusters@listData$UMAP))) {
        mon_obj <- monocle3::cluster_cells(mon_obj)
    }

    mon_obj@clusters@listData$UMAP$partitions <- setNames(new_partition, colnames(mon_obj))
    mon_obj@clusters@listData$UMAP$clusters <- setNames(new_partition, colnames(mon_obj))

    return(mon_obj)
}

subset_monocle_by_trajectory <- function(mon_obj, start_cells = NULL, end_cells = NULL) {
    if (is.null(start_cells) || is.null(end_cells)) {
        return(mon_obj)
    }

    if (length(start_cells) == 0 || length(end_cells) == 0) {
        return(mon_obj)
    }

    if (!("UMAP" %in% names(monocle3::principal_graph_aux(mon_obj)))) {
        stop("UMAP coordinates not found in the Monocle object. Please run UMAP identification first.")
    }

    nodes_start <- monocle3::principal_graph_aux(mon_obj)$UMAP$pr_graph_cell_proj_closest_vertex[, 1][start_cells]
    nodes_end <- monocle3::principal_graph_aux(mon_obj)$UMAP$pr_graph_cell_proj_closest_vertex[, 1][end_cells]

    # NOTE for now, the best way to select the subtrajectory given a set of start and end cells is to find the cells
    # between the first starting cell and end cells and the cells between the first ending cells and start cells.
    # However, we do not know if the "first" end or start cell is actually at the beginning or end of the trajectory and
    # we actually lose some cells because of that.
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

    all_cells_between <- union(
        colnames(cells_between_start),
        colnames(cells_between_end)
    )

    all_vertices_between <- union(
        igraph::vertex_attr(cells_between_end@principal_graph$UMAP, "name"),
        igraph::vertex_attr(cells_between_start@principal_graph$UMAP, "name")
    )

    mon_obj <- mon_obj[, all_cells_between]
    gc()

    # update the principal graph
    monocle3::principal_graph(mon_obj)$UMAP <- igraph::induced_subgraph(
        monocle3::principal_graph(mon_obj)$UMAP,
        all_vertices_between
    )

    # refactor the metadata
    for (mtd in colnames(mon_obj@colData)) {
        if (!(is.character(mon_obj@colData[[mtd]]) || is.factor(mon_obj@colData[[mtd]]))) {
            next
        }

        mon_obj@colData[[mtd]] <- factor(as.character(mon_obj@colData[[mtd]]))
    }

    # refactor the special fields `clusters` and `partitions`
    mon_obj@clusters@listData$UMAP$clusters <- setNames(factor(as.numeric(monocle3::clusters(mon_obj))), colnames(mon_obj))
    levels(mon_obj@clusters@listData$UMAP$clusters) <- seq_len(nlevels(mon_obj@clusters@listData$UMAP$clusters))
    mon_obj@clusters@listData$UMAP$partitions <- setNames(factor(as.numeric(monocle3::partitions(mon_obj))), colnames(mon_obj))
    levels(mon_obj@clusters@listData$UMAP$partitions) <- seq_len(nlevels(mon_obj@clusters@listData$UMAP$partitions))

    return(mon_obj)
}

# the code follows the logic from the monocle3 repository https://github.dev/cole-trapnell-lab/monocle3/blob/master/R/learn_graph.R
# the only addition is the number of nodes per log 10 cells that can be schosen by the user (not accessible via monocle)
custom_learn_graph <- function(mon_obj,
                               nodes_per_log10_cells = 30,
                               learn_graph_controls = list(
                                    eps = 1e-5,
                                    maxiter = 100
                               ),
                               use_partition = FALSE,
                               use_closed_loops = FALSE,
                               verbose = FALSE,
                               metadata_for_nodes = NULL) {
    cal_ncenter <- function(num_cell_communities, ncells, nodes_per_log10_cells=30) {
        round(num_cell_communities * nodes_per_log10_cells * log10(ncells))
    }

    if (is.null(metadata_for_nodes)) {
        n_nodes <- cal_ncenter(1, ncol(mon_obj), nodes_per_log10_cells)
    } else {
        if (!(metadata_for_nodes %in% colnames(mon_obj@colData))) {
            stop("The metadata column for the nodes was not found in the Monocle object")
        }

        sum_coms <- summary(mon_obj@colData[[metadata_for_nodes]])
        n_nodes <- sum(sapply(sum_coms, function(x) cal_ncenter(1, x, nodes_per_log10_cells)))
    }

    if (is.null(learn_graph_controls)) {
        learn_graph_controls <- list()
    }

    learn_graph_controls$ncenter <- n_nodes

    if (!("clusters" %in% names(mon_obj@clusters@listData$UMAP))) {
        mon_obj <- monocle3::cluster_cells(mon_obj)
    }

    if (!is.null(metadata_for_nodes)) {
        mon_obj <- update_mononcle_partition(mon_obj, mon_obj@colData[[metadata_for_nodes]])
    }

    if (use_partition) {
        learn_graph_controls$ncenter <- NULL

        
    }

    return(monocle3::learn_graph(
        cds = mon_obj,
        close_loop = use_closed_loops,
        use_partition = use_partition,
        learn_graph_control = learn_graph_controls,
        verbose = verbose
    ))
}

custom_pseudotime_ordering <- function(monocle_object,
                                       start_cells = NULL,
                                       end_cells = NULL) {

    if (is.null(start_cells) && is.null(end_cells)) {
        stop("Both start and end cells cannot be NULL")
    }

    is_reverse <- FALSE
    if (is.null(start_cells)) {
        start_cells <- end_cells
        end_cells <- NULL
        is_reverse <- TRUE
    }

    if (!is.null(end_cells)) {
        temp_mon_obj <- subset_monocle_by_trajectory(monocle_object, start_cells, end_cells)
        temp_mon_obj <- monocle3::order_cells(temp_mon_obj, root_cells = start_cells)

        psd_values <- monocle3::pseudotime(temp_mon_obj)
    } else {
        monocle_object <- monocle3::order_cells(monocle_object, root_cells = start_cells)
        psd_values <- monocle3::pseudotime(monocle_object)
    }

    if (is_reverse) {
        psd_values <- max(psd_values) - psd_values
    }

    entire_psd <- rep(NA, ncol(monocle_object))
    names(entire_psd) <- colnames(monocle_object)
    entire_psd[names(psd_values)] <- psd_values

    return(entire_psd)
}

get_pseudotime_recommendation <- function(monocle_object) {
    discrete_groups <- list()
    for (mtd_name in colnames(monocle_object@colData)) {
        if (inherits(monocle_object@colData[[mtd_name]], c("factor", "character"))) {
            discrete_groups[[mtd_name]] <- split(colnames(monocle_object), monocle_object@colData[[mtd_name]])
        }
    }

    umap_df <- monocle_object@int_colData$reducedDims$UMAP
    current_iqr <- 0

    for (mtd_name in names(discrete_groups)) {
        for (mtd_group in names(discrete_groups[[mtd_name]])) {
            cell_group <- discrete_groups[[mtd_name]][[mtd_group]]
            if (length(cell_group) == 0) {
                next
            }

            filtered_cells <- filter_central_cells_from_group(
                cell_group = cell_group,
                umap_embedding = umap_df,
                n_points = 5
            )
            monocle_object <- monocle3::order_cells(monocle_object, root_cells = filtered_cells)

            fvn <- fivenum(monocle3::pseudotime(monocle_object))
            iqr_psd <- fvn[4] - fvn[2]

            if (iqr_psd > current_iqr) {
                current_iqr <- iqr_psd
                recommended_mtd_group <- mtd_group
                recommended_mtd_name <- mtd_name
                recommended_pseudotime <- monocle3::pseudotime(monocle_object)
            }
        }
    }

    return(list(
        recommended_mtd_name = recommended_mtd_name,
        recommended_mtd_group = recommended_mtd_group,
        recommended_pseudotime = recommended_pseudotime
    ))
}
