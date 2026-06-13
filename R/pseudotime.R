#' Update the partition of the Monocle object
#' 
#' @description This function updates the partition of the Monocle object,
#' which is used for tasks such as trajectory analysis.
#' 
#' @param mon_obj The Monocle object.
#' @param new_partition The new partition to be updated.
#' 
#' @return The Monocle object with the updated partition.
update_monocle_partition <- function(mon_obj, new_partition) {
    new_partition <- factor(as.numeric(factor(new_partition)))

    if (!("partitions" %in% names(mon_obj@clusters@listData$UMAP))) {
        mon_obj <- monocle3::cluster_cells(mon_obj)
    }

    mon_obj@clusters@listData$UMAP$partitions <- stats::setNames(new_partition, colnames(mon_obj))
    mon_obj@clusters@listData$UMAP$clusters <- stats::setNames(new_partition, colnames(mon_obj))

    return(mon_obj)
}

#' Subset the Monocle object by a trajectory
#' 
#' @description This function subsets the Monocle object based on a subgraph
#' of the inferrred trajectory. The subgraph is defined using start and end
#' points, defined by cells.
#' 
#' @param mon_obj A monocle object.
#' @param start_cells A list of cells names that define the start of the
#' selected trajectory. If NULL, no subsetting is performed. Defaults to NULL.
#' @param end_cells A list of cells names that define the end of the
#' selected trajectory. If NULL, no subsetting is performed. Defaults to NULL.
#' 
#' @return The Monocle object subsetted by the trajectory.
#' @export
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
    mon_obj@clusters@listData$UMAP$clusters <- stats::setNames(factor(as.numeric(monocle3::clusters(mon_obj))), colnames(mon_obj))
    levels(mon_obj@clusters@listData$UMAP$clusters) <- seq_len(nlevels(mon_obj@clusters@listData$UMAP$clusters))
    mon_obj@clusters@listData$UMAP$partitions <- stats::setNames(factor(as.numeric(monocle3::partitions(mon_obj))), colnames(mon_obj))
    levels(mon_obj@clusters@listData$UMAP$partitions) <- seq_len(nlevels(mon_obj@clusters@listData$UMAP$partitions))

    return(mon_obj)
}

#' Learn the graph of the Monocle object
#' 
#' @description This function follows the logic of the `learn_graph` function
#' from the monocle3 package. The main difference consists in providing more
#' flexibility to the user, such as the number of nodes per log10 cells and
#' the metadata used as partition.
#' 
#' @param mon_obj The Monocle object.
#' @param nodes_per_log10_cells The number of trajectory nodes created per
#' log10 cells. Defaults to 30.
#' @param learn_graph_controls A list of control parameters, as defined in the
#' `learn_graph` function from Monocle. Defaults to list(eps = 1e-5,
#' maxiter = 100).
#' @param use_partition A logical value indicating if the partition should be
#' used for learning the graph. Defaults to FALSE.
#' @param use_closed_loops A logical value indicating if circular paths can be
#' formed in the trajectory graph. Defaults to FALSE.
#' @param verbose Parameter that is passed to the `learn_graph` function and
#' prints the progress.
#' @param metadata_for_nodes The metadata column that should be used to
#' update the partition of the Monocle object. This is used when `use_partition`
#' is set to TRUE. Defaults to NULL, meaning the partition is not updated.
#' 
#' @return The Monocle object with the learned graph.
#' @export
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
        mon_obj <- update_monocle_partition(mon_obj, mon_obj@colData[[metadata_for_nodes]])
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

#' Build a Trajectory Helper Object
#'
#' @description Extracts graph nodes, edges and projected cell-to-node
#' assignments from a Monocle3 object for downstream plotting and analysis.
#'
#' @param monocle_object A Monocle3 object.
#' @param reduction_name Name of the trajectory reduction to use.
#'
#' @return A list with graph, node positions, edge table and closest vertices.
#' @export
get_trajectory_object <- function(monocle_object, reduction_name = "UMAP") {
    if (is.null(reduction_name) ||!(reduction_name %in% names(monocle_object@principal_graph))) {
        reduction_name <- names(monocle_object@principal_graph)[1]
        warning(paste("Reduction name", reduction_name, "not found in monocle object. Using", reduction_name, "instead."))
    }

    if (is.null(reduction_name)) {
        stop("No valid reduction name provided and no reductions found in monocle object.")
    }

    trajectory_object <- list()
    trajectory_object$graph <- monocle_object@principal_graph[[reduction_name]]

    node_positions <- t(monocle_object@principal_graph_aux[[reduction_name]]$dp_mst)
    node_positions <- as.data.frame(node_positions)
    node_positions$node <- rownames(node_positions)
    node_positions$degree <- igraph::degree(trajectory_object$graph, mode = "in")

    trajectory_object$node_positions <- node_positions

    edges_df <- igraph::as_data_frame(trajectory_object$graph)
    edges_df$weight <- NULL # NOTE I don't think the weights are needed in our context
    edges_df$x_from <- node_positions[edges_df$from, 1]
    edges_df$y_from <- node_positions[edges_df$from, 2]
    edges_df$x_to <- node_positions[edges_df$to, 1]
    edges_df$y_to <- node_positions[edges_df$to, 2]
    trajectory_object$edges_df <- edges_df

    trajectory_object$closest_vertex <- monocle_object@principal_graph_aux[[reduction_name]]$pr_graph_cell_proj_closest_vertex[,1]

    trajectory_object$n_nodes <- nrow(node_positions)
    trajectory_object$n_edges <- nrow(edges_df)

    return(trajectory_object)
}

#' Plot a Trajectory Graph
#'
#' @description Plots trajectory edges and optionally overlays selected node
#' degrees. This plot is meant to be used as an additional layer in an UMAP
#' plot.
#'
#' @param trajectory_object A trajectory object produced by
#' `get_trajectory_object`.
#' @param edge_size Line width for trajectory edges.
#' @param edge_alpha Transparency of trajectory edges.
#' @param edge_colour Colour of trajectory edges.
#' @param plot_nodes Which node-degree groups to display: 1 for leaves, 2 for
#' intermediate nodes, 3 for branching points. Setting to 0 will hide all nodes.
#' @param node_size Point size for trajectory nodes.
#' @param node_colours Named colours for node-degree groups.
#' @param node_label_size Text size of node labels.
#' @param node_label_vjust Vertical adjustment for node labels.
#'
#' @return A ggplot object of the trajectory.
#' @export
plot_trajectory_graph <- function(
    trajectory_object,
    edge_size = 0.5,
    edge_alpha = 1,
    edge_colour = "black",
    plot_nodes = c(1, 2, 3),
    node_size = 3,
    node_colours = stats::setNames(c("red", "blue", "green"), c(1, 2, 3)),
    node_label_size = 3,
    node_label_vjust = -1
) {
    if (is.null(trajectory_object)) {
        return(ggplot2::ggplot() + ggplot2::theme_void())
    }
    edges_df <- trajectory_object$edges_df

    edge_gplot <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = edges_df, ggplot2::aes(x = .data$x_from, y = .data$y_from, xend = .data$x_to, yend = .data$y_to), linewidth = edge_size, alpha = edge_alpha, color = edge_colour) +
        ggplot2::theme_void()

    if (is.null(plot_nodes) || (length(plot_nodes) == 1 && plot_nodes == 0)) {
        return(edge_gplot)
    }

    plot_nodes <- intersect(plot_nodes, 1:3)

    node_positions <- trajectory_object$node_positions
    node_positions$degree[node_positions$degree > 3] <- 3
    node_positions <- node_positions %>%
        dplyr::filter(.data$degree %in% plot_nodes) %>%
        dplyr::mutate(color = factor(.data$degree, levels = plot_nodes))
    dimension_columns <- colnames(node_positions)[seq_len(2)]

    return(
        edge_gplot +
        ggplot2::geom_point(data = node_positions, ggplot2::aes(x = .data[[dimension_columns[1]]], y = .data[[dimension_columns[2]]], color = .data$color), size = node_size) +
        ggplot2::scale_color_manual(values = node_colours) +
        ggplot2::geom_text(data = node_positions, ggplot2::aes(x = .data[[dimension_columns[1]]], y = .data[[dimension_columns[2]]], label = .data$node), size = node_label_size, vjust = node_label_vjust)
    )
}


#' Order the cells by pseudotime
#' 
#' @description This function follows the logic of the `order_cells` function
#' from the monocle3 package. The main difference consists in allowing the user
#' to define both the start and the end points of the trajectory. In this case,
#' the monocle object is subsetted.
#' 
#' @param monocle_object A monocle object.
#' @param start_cells A list of cells names that define the start of the
#' ordering. If NULL, the end_cells must be provided. Defaults to NULL.
#' @param end_cells A list of cells names that define the end of the ordering.
#' If NULL, the start_cells must be provided. If end_cells is not NULL and
#' start_cells is not NULL, subsetting is performed. If end_cells is not NULL
#' and start_cells is NULL, the ordering will be reversed. Defaults to NULL.
#' 
#' @return A named vector with the pseudotime values. In the case of subsetting,
#' the cells that are not part of the inferred subtrajectory will be assigned
#' the NA value.
#' @export
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

#' Compute an Optimal Pseudotime Range
#'
#' @description Chooses a root node based on branch-length balance over the
#' trajectory diameter and returns the associated pseudotime ordering.
#'
#' @param monocle_object A Monocle3 object with a learned trajectory graph.
#'
#' @return A list with selected start node and pseudotime values.
#' @keywords internal
get_optimal_pseudotime_range <- function(monocle_object) {
    trajectory_graph <- monocle_object@principal_graph$UMAP
    node_degrees <- igraph::degree(trajectory_graph, mode = "all")

    longest_path <- names(igraph::get_diameter(trajectory_graph, directed = FALSE, weights = igraph::E(trajectory_graph)$weight))
    node_branch_length <- rep(0, length(longest_path))
    names(node_branch_length) <- longest_path
    for (i in longest_path) {
        if (node_degrees[i] <= 2) {
            next
        }

        # create a subgraph without the main branch
        deleted_vertices <- setdiff(longest_path, i)
        subgraph <- igraph::delete_vertices(trajectory_graph, deleted_vertices)
        current_distances <- igraph::distances(subgraph, v = i)
        node_branch_length[i] <- max(current_distances[current_distances != Inf], na.rm = TRUE)
    }

    # NOTE weighted approach using distances on the UMAP projection
    # dp_mst <- monocle_object@principal_graph_aux$UMAP$dp_mst
    # dist_dp_mst <- as.matrix(1 / dist(t(dp_mst)))
    # edge_list <- igraph::get.edgelist(trajectory_graph)
    # weights <- apply(edge_list, 1, function(x) {
    #     dist_dp_mst[x[1], x[2]]
    # })
    # igraph::E(trajectory_graph)$weight <- weights

    degree_first_half <- sum(node_branch_length[longest_path[seq_len(length(longest_path) / 2)]])
    degree_second_half <- sum(node_branch_length[longest_path[seq(length(longest_path) / 2 + 1, length(longest_path))]])

    start_node <- longest_path[length(longest_path)]
    if (degree_first_half < degree_second_half) {
        start_node <- longest_path[1]
    }
    monocle_object <- monocle3::order_cells(monocle_object, root_pr_nodes = start_node)

    return(list(
        start_node = start_node,
        pseudotime = monocle3::pseudotime(monocle_object)
    ))
}

#' Recommend a pseudotime ordering
#' 
#' @description This function provides a recommendation of a pseudotime ordering
#' based on the metadata available in the monocle object. The recommendation is
#' done by selecting an end node from the diameter of the trajectory graph. The
#' function identifies the metadata group whose centre is the closest to the
#' starting point.
#'
#' @param monocle_object A Monocle3 object.
#' @return A list that contains the recommended metadata column, the subgroup
#' and the pseudotime values.
#' @export
get_pseudotime_recommendation <- function(monocle_object) {
    optimal_range <- get_optimal_pseudotime_range(monocle_object)
    closest_vertex <- monocle_object@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
    discrete_groups <- list()
    for (mtd_name in colnames(monocle_object@colData)) {
        if (inherits(monocle_object@colData[[mtd_name]], c("factor", "character"))) {
            discrete_groups[[mtd_name]] <- split(colnames(monocle_object), monocle_object@colData[[mtd_name]])
        }
    }

    umap_df <- monocle_object@int_colData$reducedDims$UMAP
    node_positions <- monocle_object@principal_graph_aux$UMAP$dp_mst
    best_criteria <- NULL

    for (mtd_name in names(discrete_groups)) {
        for (mtd_group in names(discrete_groups[[mtd_name]])) {
            cell_group <- discrete_groups[[mtd_name]][[mtd_group]]
            if (length(cell_group) == 0) {
                next
            }

            filtered_cells <- filter_central_cells_from_group(
                cell_group = cell_group,
                umap_embedding = umap_df,
                n_points = 19
            )

            # approach 1: count the most overlaps to the starting point
            # filtered_node_ids <- paste0("Y_", closest_vertex[filtered_cells, 1])
            # criteria_value <- sum(filtered_node_ids == optimal_range$start_node) / length(filtered_node_ids)

            # approach 2: calculate the average UMAP distance to the starting point
            start_node_position <- t(node_positions)[optimal_range$start_node, ]
            umap_dist <- umap_df[filtered_cells, , drop = FALSE]
            umap_dist[, 1] <- umap_dist[, 1] - start_node_position[1]
            umap_dist[, 2] <- umap_dist[, 2] - start_node_position[2]
            umap_dist <- sqrt(umap_dist[, 1]^2 + umap_dist[, 2]^2)
            criteria_value <- 1 / mean(umap_dist, na.rm = TRUE)

            # NOTE uncomment if you want to calculate the pseudotime again
            # monocle_object <- monocle3::order_cells(monocle_object, root_cells = filtered_cells)
            # criteria_value <- recommendation_criteria(monocle3::pseudotime(monocle_object))

            if (is.null(best_criteria) || criteria_value > best_criteria) {
                best_criteria <- criteria_value
                recommended_mtd_group <- mtd_group
                recommended_mtd_name <- mtd_name

                if (criteria_value == 1) {
                    break
                }
            }
        }
    }

    return(list(
        recommended_mtd_name = recommended_mtd_name,
        recommended_mtd_group = recommended_mtd_group,
        recommended_pseudotime = optimal_range$pseudotime
    ))
}


#' Order Metadata Groups by Pseudotime
#'
#' @description Orders categorical metadata levels by their median pseudotime
#' value.
#'
#' @param metadata_value A vector with categorical metadata values.
#' @param pseudotime A numeric vector with pseudotime values.
#'
#' @return A character vector of metadata levels ordered by median pseudotime.
#' @export
order_metadata_groups_by_pseudotime <- function(
    metadata_value,
    pseudotime
) {
    grouped_df <- data.frame(
        mtd = as.character(metadata_value),
        psd = pseudotime
    )
    grouped_df %>%
        dplyr::group_by(.data$mtd) %>%
        dplyr::summarise(median_psd = stats::median(.data$psd, na.rm = TRUE)) %>%
        dplyr::arrange(.data$median_psd) %>%
        dplyr::pull(.data$mtd)
}

