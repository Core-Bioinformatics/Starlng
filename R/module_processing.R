
#' Build Module Masks
#'
#' @description Builds boolean masks for module-active cells from module
#' summaries and a pseudotime mask.
#'
#' @param module_summ A named list of module summary vectors.
#' @param psd_mask A logical vector indicating cells to keep.
#' @param scale_threshold Numeric threshold used to select active cells.
#' @param top_cells_percent Percentage of top cells kept per module.
#'
#' @return A logical matrix of module masks.
#' @export
build_module_masks <- function(module_summ, psd_mask, scale_threshold = 0, top_cells_percent = 100) {
    perc_cells <- 1 - top_cells_percent / 100
    stopifnot(
        perc_cells >= 0,
        perc_cells <= 1,
        scale_threshold >= 0,
        scale_threshold <= 1,
        !is.null(module_summ),
        !is.null(psd_mask)
    )

    mask_nrow <- length(module_summ)
    mask_ncol <- length(psd_mask)
    module_mask <- matrix(FALSE, nrow = mask_nrow, ncol = mask_ncol)
    rownames(module_mask) <- names(module_summ)

    for (i in seq_along(module_summ)) {
        module_mask[i, ] <- (module_summ[[i]] > scale_threshold) & psd_mask
        if (sum(module_mask[i, ]) < 5) {
            next
        }

        value_thresh <- stats::quantile(module_summ[[i]][module_mask[i, ]], perc_cells)
        mx_val <- max(module_summ[[i]])
        if (value_thresh == mx_val) {
            value_thresh <- 0.99 * mx_val
        }
        module_mask[i, ] <- (module_summ[[i]] >= value_thresh) & psd_mask
    }

    module_mask
}

#' Compute Module Pairwise Tables
#'
#' @description Computes pairwise overlap, union and Spearman statistics between
#' the populations described by the module masks and summaries.
#'
#' @param module_mask A logical matrix of module masks.
#' @param module_summ A named list of module summary vectors.
#'
#' @return A list with pairwise overlap and correlation matrices.
#' @export
compute_module_pairwise_tables <- function(module_mask, module_summ) {
    module_names <- rownames(module_mask)
    mask_nrow <- nrow(module_mask)

    module_intersect_cells <- matrix(0, nrow = mask_nrow, ncol = mask_nrow)
    rownames(module_intersect_cells) <- module_names
    colnames(module_intersect_cells) <- module_names

    if (mask_nrow > 1) {
        for (i in seq(from = 1, to = mask_nrow - 1)) {
            for (j in seq(from = i + 1, to = mask_nrow)) {
                ncommon_cells <- sum(module_mask[i, ] & module_mask[j, ])
                module_intersect_cells[i, j] <- ncommon_cells
                module_intersect_cells[j, i] <- ncommon_cells
            }
        }
    }

    for (i in seq_len(mask_nrow)) {
        unique_cells <- module_mask[i, ]
        for (j in seq_len(mask_nrow)) {
            if (i == j) {
                next
            }
            unique_cells <- unique_cells & (!module_mask[j, ])
        }
        module_intersect_cells[i, i] <- sum(unique_cells)
    }

    module_spearman_matrix <- matrix(NA, nrow = mask_nrow, ncol = mask_nrow)
    module_union_cells <- matrix(0, nrow = mask_nrow, ncol = mask_nrow)
    rownames(module_spearman_matrix) <- module_names
    colnames(module_spearman_matrix) <- module_names
    rownames(module_union_cells) <- module_names
    colnames(module_union_cells) <- module_names

    for (i in seq_along(module_summ)) {
        for (j in seq_along(module_summ)) {
            if (i >= j) {
                next
            }

            mask1 <- module_mask[i, ]
            mask2 <- module_mask[j, ]
            united_mask <- mask1 | mask2
            expr1 <- module_summ[[i]][united_mask]
            expr2 <- module_summ[[j]][united_mask]

            module_union_cells[i, j] <- sum(united_mask)
            module_union_cells[j, i] <- sum(united_mask)

            if (length(expr1) < 5 || length(expr2) < 5) {
                module_spearman_matrix[i, j] <- NA
                module_spearman_matrix[j, i] <- NA
                next
            }

            sp_val <- suppressWarnings(
                stats::cor(expr1, expr2, method = "spearman")
            )
            sp_val <- round(sp_val, 2)
            module_spearman_matrix[i, j] <- sp_val
            module_spearman_matrix[j, i] <- sp_val
        }
    }

    list(
        module_intersect_cells = module_intersect_cells,
        module_union_cells = module_union_cells,
        module_spearman_matrix = module_spearman_matrix
    )
}

#' Compute Module Statistics
#'
#' @description Builds per-cell module statistics using pseudotime and UMAP
#' distances.
#'
#' @param module_summ A named list of module summary vectors.
#' @param module_mask A logical matrix of module masks.
#' @param psd_value A numeric pseudotime vector.
#' @param umap_df A matrix or data frame with UMAP coordinates.
#' @param centroid Logical indicating whether UMAP distance should be computed
#' from the centroid.
#'
#' @return A data frame with module statistics, or NULL if no cells are selected.
#' @export
get_module_stats <- function(module_summ, module_mask, psd_value, umap_df, centroid = TRUE) {
    modules_stats <- NULL

    for (i in seq_along(module_summ)) {
        if (sum(module_mask[i, ]) == 0) {
            next
        }
        temp_psd_val <- psd_value[module_mask[i, ]]
        if (all(is.na(temp_psd_val))) {
            next
        }

        temp_df <- data.frame(
            avg_summary = module_summ[[i]][module_mask[i, ]],
            psd_value = temp_psd_val,
            umap_distance = calculate_umap_average_distance(
                umap_df = umap_df,
                selected_cells = which(module_mask[i, ]),
                centroid = centroid
            )
        ) %>% dplyr::filter(!is.na(.data$psd_value))
        temp_df$module <- names(module_summ)[i]

        if (is.null(modules_stats)) {
            modules_stats <- temp_df
        } else {
            modules_stats <- rbind(modules_stats, temp_df)
        }
    }

    if (is.null(modules_stats)) {
        return(NULL)
    }

    modules_stats$module <- factor(modules_stats$module, levels = names(module_summ))
    return(modules_stats)
}

#' Summarize Module Statistics
#'
#' @description Summarizes per-cell module statistics into per-module summaries.
#'
#' @param modules_stats A data frame returned by get_module_stats.
#' @param gene_modules A named list mapping modules to genes.
#'
#' @return A data frame of module-level summary statistics.
#' @export
summarise_module_stats <- function(modules_stats, gene_modules) {
    modules_stats$module <- as.character(modules_stats$module)
    summary_stats <- modules_stats %>%
        dplyr::group_by(.data$module) %>%
        dplyr::summarise(
            n_cells = length(.data$avg_summary),
            avg_summary = round(mean(.data$avg_summary, na.rm = TRUE), 3),
            median_pseudotime = round(stats::median(.data$psd_value, na.rm = TRUE), 3),
            iqr_pseudotime = round(stats::IQR(.data$psd_value, na.rm = TRUE), 3),
            median_umap_distance = round(stats::median(.data$umap_distance, na.rm = TRUE), 3)
        )
    summary_stats$n_genes <- sapply(summary_stats$module, function(mod) length(gene_modules[[mod]]))
    current_cols <- colnames(summary_stats)
    new_cols <- c(current_cols[1], "n_genes", current_cols[2:(length(current_cols) - 1)])
    summary_stats <- as.data.frame(summary_stats)[, new_cols]
    rownames(summary_stats) <- summary_stats$module
    return(summary_stats)
}


#' Detect Module Outliers
#'
#' @description Classifies modules as outliers, redundant, or valid based on
#' pseudotime spread, UMAP distance and coverage heuristics. A module is
#' classified as outlier if its IQR and UMAP distance have a MAD z-score greater
#' than 3.5 in absolute value. Modules that fit in the 85% of the pseudotime
#' IQR or the average UMAP distance of the baseline (the entire dataset) are
#' considered outliers as well.
#' Redundancy is based on the percentage of new cells covered by the module
#' compared to the already covered cells (the modules are ordered based on their
#' median UMAP distance and IQR pseudotime, so the best modules are evaluated
#' first). A second round allows redundant modules to be labeled as non-redundant
#' if they provide a significant percentage of new cells compared to the already
#' covered population (by default, the threshold is the median F1 score of
#' the non-redundant modules).
#'
#' @param modules_stats A data frame of module summary statistics.
#' @param cell_masks A logical matrix of module-to-cell memberships.
#' @param psd_value A numeric pseudotime vector.
#' @param thresh_psd_good Threshold for accepting a module as non-redundant.
#' If NULL, it will be set to the 25% quantile of the IQR pseudotime of the
#' non-outlier modules.
#' @param thresh_psd_bad Threshold for flagging a module as an outlier. If NULL,
#' it will be set to the 85% quantile of the IQR pseudotime of the entire dataset.
#' @param umap_dist_threshold Optional threshold for median UMAP distance.
#'
#' @return A list with outlier labels and coverage evolution data.
#' @export
detect_outlier <- function(modules_stats, cell_masks, psd_value, thresh_psd_good = NULL, thresh_psd_bad = NULL, umap_dist_threshold = NULL) {
    if (is.null(thresh_psd_bad)) {
        thresh_psd_bad <- stats::IQR(psd_value, na.rm = TRUE) * 0.85
    }
    modules_stats <- modules_stats %>% dplyr::arrange(
        .data$median_umap_distance,
        .data$iqr_pseudotime
    )
    cell_masks <- cell_masks[modules_stats$module, ]

    mad_z_scores_psd <- stats::setNames(mad_z_score(modules_stats$iqr_pseudotime), modules_stats$module)
    mad_z_scores_umap <- stats::setNames(mad_z_score(modules_stats$median_umap_distance), modules_stats$module)

    outlier_output <- rep("no", nrow(modules_stats))
    names(outlier_output) <- modules_stats$module

    ncells_per_module <- apply(cell_masks, 1, sum)

    outlier_output[abs(mad_z_scores_psd) > 3.5 | abs(mad_z_scores_umap) > 3.5] <- "redundant"
    outlier_output[abs(mad_z_scores_psd) > 3.5 & abs(mad_z_scores_umap) > 3.5] <- "yes"
    outlier_output[modules_stats$iqr_pseudotime >= thresh_psd_bad] <- "yes"
    outlier_output[ncells_per_module <= 10] <- "yes"

    if (!is.null(umap_dist_threshold)) {
        outlier_output[modules_stats$median_umap_distance >= umap_dist_threshold] <- "yes"
    }

    if (is.null(thresh_psd_good)) {
        thresh_psd_good <- stats::quantile(modules_stats$iqr_pseudotime[outlier_output == "no"], 0.25, na.rm = TRUE)
    }

    covered_mask <- rep(FALSE, ncol(cell_masks))

    coverage_evolution_df <- NULL
    for (i in seq_len(nrow(modules_stats))) {
        module_name <- modules_stats$module[i]
        is_already_labeled <- outlier_output[module_name] != "no"
        if (is_already_labeled) {
            next
        }

        current_mask <- cell_masks[module_name, ]
        ncells <- sum(current_mask)
        if (ncells < 10) {
            outlier_output[module_name] <- "redundant"
            next
        }
        temp_mask <- covered_mask | current_mask
        nunique <- sum(temp_mask) - sum(covered_mask)
        iqr <- modules_stats$iqr_pseudotime[i]
        eligible <- TRUE

        if (nunique / ncells <= 0.5) {
            eligible <- iqr < thresh_psd_good
        }

        potential_coverage_added <- nunique / length(covered_mask) + 1e-10
        potential_unique_percentage <- nunique / ncells + 1e-10
        percentage_iqr <- 1 / (iqr + 1e-10)
        potential_f1_score <- 2 / (1 / potential_unique_percentage + 1 / potential_coverage_added) #+ 1 / percentage_iqr)

        if (!eligible) {
            outlier_output[module_name] <- "redundant"
        } else {
            new_coverage <- data.frame(
                added_module = module_name,
                coverage_added = potential_coverage_added,
                percentage_unique = potential_unique_percentage,
                percentage_iqr = percentage_iqr,
                f1_score = potential_f1_score,
                module_iqr = iqr,
                median_pseudotime = modules_stats$median_pseudotime[i],
                median_umap_distance = modules_stats$median_umap_distance[i]
            )
            if (is.null(coverage_evolution_df)) {
                coverage_evolution_df <- new_coverage
            } else {
                coverage_evolution_df <- rbind(coverage_evolution_df, new_coverage)
            }
            covered_mask <- temp_mask
        }
    }

    # provide a second chance for the redundant modules
    # although most of their population overlap with other better modules, there's
    # a chance that the specific difference is not covered at all by the other modules
    # allow the module if the percentage of new cells is above 25%
    threshold_f1_score <- NULL
    if (!is.null(coverage_evolution_df)) {
        threshold_f1_score <- stats::median(coverage_evolution_df$f1_score, na.rm = TRUE)
    }

    for (i in seq_len(nrow(modules_stats))) {
        module_name <- modules_stats$module[i]
        if (outlier_output[module_name] != "redundant") {
            next
        }

        current_mask <- cell_masks[module_name, ]
        ncells <- sum(current_mask)
        if (ncells < 10) {
            next
        }
        temp_mask <- covered_mask | current_mask
        nunique <- sum(temp_mask) - sum(covered_mask)
        iqr <- modules_stats$iqr_pseudotime[i]

        potential_percentage_unique <- nunique / ncells + 1e-10
        potential_coverage_added <- nunique / length(covered_mask) + 1e-10
        percentage_iqr <- 1 - 1 / (iqr + 1e-10)
        potential_f1_score <- 2 / (1 / potential_percentage_unique + 1 / potential_coverage_added)# + 1 / percentage_iqr)

        if (!is.null(threshold_f1_score) && potential_f1_score < threshold_f1_score) {
            next
        }
        
        outlier_output[module_name] <- "no"
        new_coverage <- data.frame(
            added_module = module_name,
            coverage_added = potential_coverage_added,
            percentage_unique = potential_percentage_unique,
            percentage_iqr = percentage_iqr,
            f1_score = potential_f1_score,
            module_iqr = iqr,
            median_pseudotime = modules_stats$median_pseudotime[i],
            median_umap_distance = modules_stats$median_umap_distance[i]
        )
        if (is.null(coverage_evolution_df)) {
            coverage_evolution_df <- new_coverage
        } else {
            coverage_evolution_df <- rbind(coverage_evolution_df, new_coverage)
        }
        covered_mask <- temp_mask
        threshold_f1_score <- stats::median(coverage_evolution_df$f1_score, na.rm = TRUE)
    }

    return(list(
        "outlier_output" = outlier_output,
        "coverage_evolution_df" = coverage_evolution_df
    ))
}

#' Annotate Module Outliers
#'
#' @description Adds outlier labels to module summaries using module masks and
#' pseudotime information.
#'
#' @param modules_stats_summary A data frame of module summary statistics.
#' @param module_mask A logical matrix of module masks.
#' @param psd_value A numeric pseudotime vector.
#' @param umap_dist_threshold Optional threshold for median UMAP distance.
#'
#' @return The module summary data frame with an is_outlier column.
#' @keywords internal
annotate_module_outliers <- function(modules_stats_summary, module_mask, psd_value, umap_dist_threshold = NULL) {
    if (is.null(modules_stats_summary) || nrow(modules_stats_summary) == 0) {
        return(modules_stats_summary)
    }

    outlier_result <- detect_outlier(
        modules_stats = modules_stats_summary,
        cell_masks = module_mask,
        psd_value = psd_value,
        thresh_psd_good = NULL, #psd_span / 10,
        thresh_psd_bad = NULL, #psd_span / 3
        umap_dist_threshold = umap_dist_threshold
    )

    modules_stats_summary$is_outlier <- outlier_result$outlier_output[as.character(modules_stats_summary$module)]

    modules_stats_summary %>%
        dplyr::arrange(.data$median_pseudotime, .data$iqr_pseudotime, .data$median_umap_distance)
}


#### module connectivity ####

#' Find Closest Trajectory Node to Module
#'
#' @description Finds the closest trajectory node to a module using the UMAP
#' distance between the module centroid and the node location.
#'
#' @param trajectory_object A trajectory object containing node positions.
#' @param cell_umap A cell embedding matrix or data frame with two dimensions.
#' @param module_expr A numeric vector with the aggregate expression of the module.
#' If a logical vector is provided, it will be interpreted as the mask indicating
#' the module population. Providing a list will return a named vector of closest
#' nodes for each module.
#' @param expression_threshold Expression threshold used to define module-active
#' cells. Will not be used if `module_expr` is already logical. Defaults to 0.
#' @param expression_percentile Optional percentile threshold for module-active
#' cells. Will be applied after `expression_threshold` if provided. Defaults to 0.
#' @param scale Logical indicating whether expression values should be scaled 
#' between 0 and 1. Will not be used if `module_expr` is already logical.
#' Defaults to TRUE.
#'
#' @return The closest node name, or a named vector of node names for list
#' input.
#' @export
get_closest_node_to_module <- function(
    trajectory_object,
    cell_umap,
    module_expr,
    expression_threshold = 0,
    expression_percentile = 0,
    scale = TRUE
) {
    if (inherits(module_expr, "list")) {
        return(sapply(module_expr, function(expr) {
            get_closest_node_to_module(trajectory_object, cell_umap, expr, expression_threshold, expression_percentile, scale)
        }))
    }
    centroid_pop <- get_module_centroid(module_expr, cell_umap, expression_threshold, expression_percentile, scale)
    node_distance <- sqrt((trajectory_object$node_positions[, 1] - centroid_pop[1])^2 + (trajectory_object$node_positions[, 2] - centroid_pop[2])^2)
    closest_node <- trajectory_object$node_positions$node[which.min(node_distance)]
    return(closest_node)
}

#' Compute Module Centroid in Embedding Space
#'
#' @description Computes the geometric median centroid of module-active cells in
#' UMAP space.
#'
#' @param module_expr A numeric/logical vector, or a list of such vectors.
#' @param cell_umap A cell embedding matrix or data frame with at least two
#' dimensions.
#' @param expression_threshold Expression threshold used to define active cells.
#' @param expression_percentile Optional percentile threshold for active cells.
#' @param scale Logical indicating whether expression values should be scaled 
#' between 0 and 1.
#'
#' @return A numeric centroid vector, or a centroid matrix for list input.
#' @export
get_module_centroid <- function(module_expr, cell_umap, expression_threshold = 0, expression_percentile = 0, scale = TRUE) {
    if (inherits(module_expr, "list")) {
        df <- do.call(rbind, lapply(module_expr, function(expr) {
            get_module_centroid(expr, cell_umap, expression_threshold, expression_percentile, scale)
        }))
        colnames(df) <- colnames(cell_umap)[1:2]
        rownames(df) <- names(module_expr)
        return(df)
    }

    if (!is.logical(module_expr)) {
        if (scale) {
            module_expr <- (module_expr - min(module_expr)) / (max(module_expr) - min(module_expr))
        }

        mask_expression <- module_expr > expression_threshold

        if (expression_percentile > 0) {
            expression_threshold <- stats::quantile(module_expr[mask_expression], probs = expression_percentile)
            module_expr <- module_expr > expression_threshold
        } else {
            module_expr <- mask_expression
        }
    }
    centroid_pop <- Gmedian::Gmedian(cell_umap[module_expr, 1:2, drop = FALSE])
    return(centroid_pop)
}

#' Build Module Transition Adjacency
#'
#' @description Infers module connectivity by mapping modules to trajectory
#' nodes and pruning triangles based on module similarity.
#'
#' @param trajectory_object A trajectory object containing graph and node data.
#' @param closest_module Named vector mapping modules to closest nodes.
#' @param start_node Optional node used to initialize graph traversal. If NULL,
#' a node with degree 1 will be selected as the starting point.
#' @param similarity_values Optional matrix used to break triangles. The values
#' will be used to calculate the distance between the triangle edges.
#'
#' @return A symmetric module adjacency matrix.
#' @export
get_module_transitions <- function(
    trajectory_object,
    closest_module,
    start_node = NULL,
    similarity_values = NULL
) {
    if (is.null(start_node)) {
        start_node <- trajectory_object$node_positions %>%
            dplyr::filter(.data$degree == 1) %>%
            dplyr::slice_head(n = 1) %>%
            dplyr::pull(.data$node)
    }

    node_adj_matrix <- igraph::as_adjacency_matrix(trajectory_object$graph)
    queue <- c(start_node)
    visited <- stats::setNames(rep(FALSE, nrow(node_adj_matrix)), rownames(node_adj_matrix))
    unique_nodes <- unique(closest_module)

    while(length(queue) > 0) {
        current_node <- queue[1]
        visited[current_node] <- TRUE
        queue <- queue[-1]
        neighbors <- names(igraph::neighbors(trajectory_object$graph, current_node))
        queue <- c(queue, neighbors[!visited[neighbors]])

        if (current_node %in% unique_nodes) {
            next
        }

        for (neighbor in neighbors) {
            node_adj_matrix[neighbor, ] <- node_adj_matrix[neighbor, ] | node_adj_matrix[current_node, ]
            node_adj_matrix[, neighbor] <- node_adj_matrix[, neighbor] | node_adj_matrix[, current_node]
            node_adj_matrix[neighbor, neighbor] <- 0
        }
        node_adj_matrix[current_node, ] <- 0
        node_adj_matrix[, current_node] <- 0
    }


    module_names <- names(closest_module)
    module_adj_matrix <- matrix(0, nrow = length(module_names), ncol = length(module_names), dimnames = list(module_names, module_names))
    for (i in seq_along(module_names)) {
        node_i <- closest_module[i]
        for (j in seq(from = i, to = length(module_names))) {
            node_j <- closest_module[j]
            if (i == j) {
                next
            }

            if (node_i == node_j) {
                module_adj_matrix[i, j] <- 1
                module_adj_matrix[j, i] <- 1
                next
            }

            module_adj_matrix[i, j] <- node_adj_matrix[node_i, node_j]
            module_adj_matrix[j, i] <- module_adj_matrix[i, j]
        }
    }

    if (is.null(similarity_values) || isFALSE(all(module_names %in% rownames(similarity_values)))) {
        return(module_adj_matrix)
    }

    # remove triangles by deleting edge of most distant modules
    while (TRUE) {
        g <- igraph::graph_from_adjacency_matrix(module_adj_matrix, mode = "undirected")
        triangles <- igraph::triangles(g)
        if (length(triangles) == 0) {
            break
        }
        for (i in seq(1, length(triangles), by = 3)) {
            tri_nodes <- triangles[i:(i+2)]
            tri_node_names <- names(igraph::V(g))[tri_nodes]
            tri_similarities <- similarity_values[tri_node_names, ]
            distance_sims <- as.matrix(stats::dist(tri_similarities))

            max_sim <- max(distance_sims, na.rm = TRUE)
            if (max_sim == 0) {
                max_sim_node1 <- tri_node_names[1]
                max_sim_node2 <- tri_node_names[2]
            } else {
                max_sim_indices <- which(distance_sims == max_sim, arr.ind = TRUE)[1, ]
                max_sim_node1 <- tri_node_names[max_sim_indices[1]]
                max_sim_node2 <- tri_node_names[max_sim_indices[2]]
            }

            module_adj_matrix[max_sim_node1, max_sim_node2] <- 0
            module_adj_matrix[max_sim_node2, max_sim_node1] <- 0
        }
    }

    return(module_adj_matrix)
}

#' Plot Module Transition Graph
#'
#' @description Plots module connectivity as a tree-like graph.
#'
#' @param module_adj_matrix A module adjacency matrix.
#' @param closest_module Named vector mapping modules to closest nodes.
#' @param start_module Optional module name used as graph root.
#' @param edge_size Edge width.
#' @param edge_alpha Edge transparency.
#' @param edge_colour Edge colour.
#' @param node_size Node size.
#' @param node_colours Node fill colour.
#' @param node_label_size Node-label text size.
#' @param node_label_vjust Vertical adjustment for node labels.
#'
#' @return A ggplot object.
#' @export
plot_module_transitions <- function(
    module_adj_matrix,
    closest_module,
    start_module = NULL,
    edge_size = 0.5,
    edge_alpha = 0.5,
    edge_colour = "black",
    node_size = 3,
    node_colours = "white",
    node_label_size = 3,
    node_label_vjust = -1
) {
    if (is.null(start_module)) {
        degree_modules <- rowSums(module_adj_matrix)
        start_module <- names(degree_modules)[which.min(degree_modules)][1]
    }
    g_modules <- igraph::graph_from_adjacency_matrix(module_adj_matrix, mode = "undirected")

    ggraph::ggraph(g_modules, layout = "tree", root = start_module) +
        ggraph::geom_edge_link(edge_width = edge_size, edge_alpha = edge_alpha, edge_colour = edge_colour) +
        ggraph::geom_node_point(size = node_size, fill = node_colours) +
        ggraph::geom_node_text(ggplot2::aes(label = igraph::V(g_modules)$name), size = node_label_size, vjust = node_label_vjust) +
        ggplot2::theme_void() +
        ggplot2::ggtitle("Module connectivity graph") +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = node_label_size * 1.5)
        )
}


####  hub genes (sorting the genes per module) #####

#' Compute Per-Module Gene Weights
#'
#' @description Computes normalized intra-module connectivity weights for
#' genes. For each gene, the sum of weights of the edges incident to it are
#' summed. These weights are normalised by the maximum weight of the module.
#'
#' @param gene_adj_matrix A gene adjacency matrix.
#' @param gene_modules A named vector or list defining gene-module assignments.
#'
#' @return A named numeric vector of module-normalized weights per gene.
#' @export
get_per_module_weight <- function(
    gene_adj_matrix,
    gene_modules
) {
    if (!inherits(gene_modules, "list")) {
        if (is.null(names(gene_modules))) {
            names(gene_modules) <- as.character(seq_along(gene_modules))
        }
        gene_modules <- split(names(gene_modules), gene_modules)
        return(get_per_module_weight(gene_adj_matrix, gene_modules))
    }

    gene_names <- rownames(gene_adj_matrix)
    module_weight <- stats::setNames(rep(0, length(gene_names)), gene_names)
    self_loops <- gene_adj_matrix[1,1] == 1

    for (module_genes in gene_modules) {
        if (length(module_genes) == 1) {
            module_weight[module_genes] <- 0
            next
        }
        current_module_weight <- Matrix::rowSums(gene_adj_matrix[module_genes, module_genes, drop = FALSE]) - self_loops
        max_weight <- max(current_module_weight)
        if (max_weight == 0) {
            module_weight[module_genes] <- 0
            next
        }
        module_weight[module_genes] <- current_module_weight / max_weight
    }
    
    return(module_weight)
}

#' Compute Gene Hub Statistics
#'
#' @description Computes precision, recall and F1 overlap statistics between
#' the population induced by each gene and the overall population described
#' by the module it belongs to. If the gene weight is provided, the function
#' also calculates the combined score (defined as normalised weight * F1 score)
#' to determine hub genes.
#'
#' @param expr_matrix A gene-by-cell expression matrix.
#' @param gene_modules A character vector of genes, or module-specific gene
#' list.
#' @param module_expr A numeric/logical vector, or a list of such vectors.
#' @param gene_expression_threshold Threshold for defining expressed genes.
#' @param gene_expression_percentile Optional percentile threshold for gene
#' expression.
#' @param module_expression_threshold Threshold for module activity.
#' @param module_expression_percentile Optional percentile threshold for module
#' activity.
#' @param scale Logical indicating whether values should be min-max scaled.
#' @param total_weight Optional vector of additional gene weights.
#'
#' @return A data frame with overlap metrics per gene.
#' @export
get_gene_overlap_stat <- function(
    expr_matrix,
    gene_modules,
    module_expr,
    gene_expression_threshold = 0,
    gene_expression_percentile = 0,
    module_expression_threshold = 0,
    module_expression_percentile = 0,
    scale = TRUE,
    total_weight = NULL
) {
    if (inherits(module_expr, "list")) {
        module_names <- names(module_expr)
        df <- do.call(rbind, lapply(module_names, function(module_name) {
            get_gene_overlap_stat(
                expr_matrix,
                gene_modules[[module_name]],
                module_expr[[module_name]],
                gene_expression_threshold,
                gene_expression_percentile,
                module_expression_threshold,
                module_expression_percentile,
                scale,
                total_weight
            )
        }))
        return(df)
    }

    if (!is.logical(module_expr)) {
        # if a mask is already provided, skip thresholding
        if (scale) {
            module_expr <- (module_expr - min(module_expr)) / (max(module_expr) - min(module_expr))
        }
        mask_module <- module_expr > module_expression_threshold

        if (module_expression_percentile > 0) {
            module_expression_threshold <- stats::quantile(module_expr[mask_module], probs = module_expression_percentile)
            module_expr <- module_expr > module_expression_threshold
        } else {
            module_expr <- mask_module
        }
    }

    gene_names <- gene_modules
    submat <- expr_matrix[gene_names, , drop = FALSE]

    if (!is.logical(submat)) {
        if (scale) {
            submat <- t(apply(submat, 1, function(x) (x - min(x)) / (max(x) - min(x))))
        }
        
        mask_gene <- submat > gene_expression_threshold

        if (gene_expression_percentile > 0) {
            thresh_value <- sapply(gene_names, function(gene) {
                gene_values <- submat[gene, module_expr]
                if (length(gene_values) == 0) {
                    return(Inf)
                }
                return(stats::quantile(gene_values, probs = gene_expression_percentile))
            })
            submat <- submat > thresh_value
        } else {
            submat <- mask_gene
        }
    }

    precision_population <- rowSums(submat[, module_expr, drop = FALSE]) / rowSums(submat) 
    recall_population <- rowSums(submat[, module_expr, drop = FALSE]) / sum(module_expr) 

    df <- data.frame(precision_population = precision_population, recall_population = recall_population)
    df$f1_score_population <- 2 * (df$precision_population * df$recall_population) / (df$precision_population + df$recall_population + 1e-10)

    if (!is.null(total_weight)) {
        if (!is.null(names(total_weight))) {
            total_weight <- total_weight[gene_names]
        } else {
            total_weight <- total_weight[seq_along(gene_names)]
        }
        
        df$total_weight <- total_weight
        df$combined_score <- df$f1_score_population * total_weight
    }
    return(df)
}


#### module plots ####


#' Filter Gene Adjacency by Module Structure
#'
#' @description Filters gene-gene edges based on module connectivity i.e. remove
#' genes between genes from different modules if those modules are not neighbours
#' in the trajectory graph. The function provides the option to sample non-hub genes,
#' as well as select a top percentage of edges based on their weight.
#'
#' @param gene_modules A list mapping modules to gene vectors.
#' @param gene_adjacency A gene adjacency matrix.
#' @param module_adjacency A module adjacency matrix, or NULL to infer one.
#' @param hub_genes Optional data frame containing hub genes.
#' @param closest_node_for_module Optional named vector mapping modules to
#' trajectory nodes.
#' @param trajectory_object Optional trajectory object used when
#' module_adjacency is NULL.
#' @param percentage_non_hub_nodes Fraction of non-hub nodes kept per module.
#' @param percentage_edges Fraction of strongest edges retained.
#'
#' @return A list with node metadata and filtered edge table.
#' @export
get_filtered_gene_adjacency <- function(
    gene_modules, # a list
    gene_adjacency,
    module_adjacency,
    hub_genes = NULL,
    closest_node_for_module = NULL,
    trajectory_object = NULL,
    percentage_non_hub_nodes = 1,
    percentage_edges = 0.2
) {
    used_genes <- unlist(gene_modules)
    used_modules <- names(gene_modules)
    gene_adjacency <- gene_adjacency[used_genes, used_genes, drop = FALSE]

    if (is.null(module_adjacency)) {
        if (is.null(closest_node_for_module) || is.null(trajectory_object)) {
            stop("If `module_adjacency` is set to NULL, both `closest_module` and `trajectory_object` should be provided.")
        }
        closest_node_for_module <- closest_node_for_module[used_modules]
        module_adjacency <- get_module_transitions(trajectory_object, closest_node_for_module)
    }

    for (i in seq_len(nrow(gene_adjacency))) {
        gene_adjacency[i, i] <- 0
    }
    for (i in used_modules) {
        genes_i <- gene_modules[[i]]
        for (j in used_modules) { 
            genes_j <- gene_modules[[j]]
            if (i == j) {
                next
            }

            if (module_adjacency[i, j] == 1) {
                next
            }

            gene_adjacency[genes_i, genes_j] <- 0
            gene_adjacency[genes_j, genes_i] <- 0
        }
    }

    if (percentage_non_hub_nodes < 1) {
        for (module in used_modules) {
            module_genes <- gene_modules[[module]]
            module_hub_genes <- intersect(hub_genes$gene, module_genes)
            module_non_hub_genes <- setdiff(module_genes, module_hub_genes)
            n_sample <- ceiling(length(module_non_hub_genes) * percentage_non_hub_nodes)
            module_genes <- c(
                module_hub_genes,
                sample(module_non_hub_genes, min(n_sample, length(module_non_hub_genes)), replace = FALSE)
            )
            if (length(module_genes) == 0) {
                gene_modules[[module]] <- NULL
                next
            }
            gene_modules[[module]] <- module_genes
        }
            
        used_genes <- unlist(gene_modules)
        used_modules <- names(gene_modules)
        gene_adjacency <- gene_adjacency[used_genes, used_genes]
    }

    nodes_df <- reshape2::melt(gene_modules)
    colnames(nodes_df) <- c("gene", "module")
    nodes_df$is_hub <- nodes_df$gene %in% hub_genes$gene
    rownames(nodes_df) <- nodes_df$gene

    if (percentage_edges < 1) {
        # get top intra
        for (module in used_modules) {
            module_genes <- gene_modules[[module]]
            sub_adj <- gene_adjacency[module_genes, module_genes, drop = FALSE]
            edge_values <- sub_adj@x
            edge_threshold <- stats::quantile(edge_values[edge_values > 0], probs = 1 - percentage_edges)
            gene_adjacency[module_genes, module_genes][gene_adjacency[module_genes, module_genes] < edge_threshold] <- 0
        }

        # get top inter
        for (i in seq_len(nrow(module_adjacency))) {
            module_i <- rownames(module_adjacency)[i]
            genes_i <- gene_modules[[module_i]]
            for (j in seq(from = i, to = nrow(module_adjacency))) {
                module_j <- rownames(module_adjacency)[j]
                genes_j <- gene_modules[[module_j]]
                if (i == j) {
                    next
                }
                if (module_adjacency[i, j] == 0) {
                    next
                }

                sub_adj <- gene_adjacency[genes_i, genes_j, drop = FALSE]
                edge_values <- sub_adj@x
                if (length(edge_values) == 0) {
                    next
                }
                edge_threshold <- stats::quantile(edge_values[edge_values > 0], probs = 1 - percentage_edges)
                gene_adjacency[genes_i, genes_j][gene_adjacency[genes_i, genes_j] < edge_threshold] <- 0
                gene_adjacency[genes_j, genes_i][gene_adjacency[genes_j, genes_i] < edge_threshold] <- 0
            }
        }
    }

    edges_df <- igraph::as_data_frame(igraph::graph_from_adjacency_matrix(gene_adjacency, mode = "undirected", weighted = TRUE))
    edges_df$module_from <- nodes_df[edges_df$from, "module"]
    edges_df$module_to <- nodes_df[edges_df$to, "module"]
    edges_df$module_assigned <- edges_df$module_from
    edges_df$module_assigned[edges_df$module_assigned != edges_df$module_to] <- "inter"

    return(list(
        nodes_df = nodes_df,
        edges_df = edges_df
    ))
}

#' Plot Gene UMAP with Hub Highlighting
#'
#' @description Plots genes in UMAP space with module-colored points and
#' weighted edges. If provided, the hub genes are highlighted with larger points
#' and labels.
#'
#' @param umap_df A gene embedding matrix or data frame.
#' @param filtered_gene_adj Output from `get_filtered_gene_adjacency`.
#' @param edge_weight_range Numeric range for edge linewidth scaling.
#' @param edge_alpha Edge transparency.
#' @param edge_colour Default edge colour.
#' @param point_size Base point size.
#' @param module_colours Optional named vector of module colours.
#' @param node_text_size Text size for hub labels.
#' @param legend_text_size Text size for legends.
#' @param hub_point_scale Size multiplier for hub genes.
#' @param hub_stroke Stroke width for hub points.
#'
#' @return A ggplot object.
#' @export
plot_gene_hub_umap <- function(
    umap_df,
    filtered_gene_adj,
    edge_weight_range = c(0.1, 1),
    edge_alpha = 0.5,
    edge_colour = "#bbbbbb",
    point_size = 2,
    module_colours = NULL,
    node_text_size = 5,
    legend_text_size = 4,
    hub_point_scale = 2,
    hub_stroke = 0.7
) {
    has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

    used_genes <- filtered_gene_adj$nodes_df$gene
    used_modules <- unique(filtered_gene_adj$nodes_df$module)

    umap_df <- umap_df[used_genes, , drop = FALSE]

    if (is.null(module_colours)) {
        n_modules <- length(used_modules)
        if (n_modules == 1) {
            module_colours <- "blue"
        } else {
            module_colours <- qualpalr::qualpal(n_modules, list(h = c(0, 360), s = c(0.2, 0.8), l = c(0.4, 0.65)))$hex
        }
        names(module_colours) <- as.character(used_modules)
    }
    module_colours[["inter"]] <- edge_colour
    
    umap_df <- as.data.frame(umap_df)
    umap_df$gene <- rownames(umap_df)
    umap_df$module <- filtered_gene_adj$nodes_df[umap_df$gene, "module"]
    umap_df$is_hub <- FALSE
    if ("is_hub" %in% colnames(filtered_gene_adj$nodes_df)) {
        umap_df$is_hub <- filtered_gene_adj$nodes_df[umap_df$gene, "is_hub"]
        umap_df$is_hub[is.na(umap_df$is_hub)] <- FALSE
    }
    hub_df <- umap_df[umap_df$is_hub, , drop = FALSE]

    umap_columns <- colnames(umap_df)[1:2]

    filtered_gene_adj$edges_df$x <- umap_df[filtered_gene_adj$edges_df$from, umap_columns[1]]
    filtered_gene_adj$edges_df$y <- umap_df[filtered_gene_adj$edges_df$from, umap_columns[2]]
    filtered_gene_adj$edges_df$xend <- umap_df[filtered_gene_adj$edges_df$to, umap_columns[1]]
    filtered_gene_adj$edges_df$yend <- umap_df[filtered_gene_adj$edges_df$to, umap_columns[2]]

    gplot_obj <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = filtered_gene_adj$edges_df, ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend, linewidth = .data$weight, color = .data$module_assigned), alpha = edge_alpha) +
        ggplot2::geom_point(data = umap_df[!umap_df$is_hub, , drop = FALSE], ggplot2::aes(x = .data[[umap_columns[1]]], y = .data[[umap_columns[2]]], color = .data$module), size = point_size, alpha = 0.85) +
        ggplot2::geom_point(data = hub_df, ggplot2::aes(x = .data[[umap_columns[1]]], y = .data[[umap_columns[2]]], color = .data$module), size = point_size * hub_point_scale, shape = 21, stroke = hub_stroke, fill = "white") +
        ggplot2::scale_color_manual(values = module_colours) +
        ggplot2::scale_linewidth_continuous(range = edge_weight_range) +
        ggplot2::theme_void() +
        ggplot2::theme(
            legend.text = ggplot2::element_text(size = legend_text_size),
            legend.title = ggplot2::element_text(size = legend_text_size * 1.2),
            plot.title = ggplot2::element_text(size = legend_text_size * 1.6)
        ) +
        ggplot2::ggtitle("UMAP of gene embedding with module colors")

    if (nrow(hub_df) > 0) {
        if (has_ggrepel) {
            gplot_obj <- gplot_obj +
                ggrepel::geom_label_repel(
                    data = hub_df,
                    ggplot2::aes(x = .data[[umap_columns[1]]], y = .data[[umap_columns[2]]], label = .data$gene, color = .data$module),
                    size = node_text_size,
                    fill = ggplot2::alpha("white", 0.85),
                    label.size = 0.15,
                    fontface = "bold",
                    max.overlaps = Inf,
                    box.padding = 0.3,
                    point.padding = 0.2,
                    min.segment.length = 0
                )
        } else {
            gplot_obj <- gplot_obj +
                ggplot2::geom_text(
                    data = hub_df,
                    ggplot2::aes(x = .data[[umap_columns[1]]], y = .data[[umap_columns[2]]], label = .data$gene, color = .data$module),
                    size = node_text_size,
                    fontface = "bold",
                    vjust = -0.8
                )
        }
    }

    return(gplot_obj)
}



#' Plot Module Trends over Pseudotime
#'
#' @description Fits smoothed module-expression trajectories over pseudotime
#' and plots module trends. If the number of points is large (> 7000 cells),
#' the function uses a GAM-based smoothing approach for efficiency, otherwise it
#' uses LOESS smoothing.
#'
#' @param expression_list A list, matrix, or data frame with module
#' expression per cell.
#' @param pseudotime Numeric pseudotime vector.
#' @param module_colours Optional named vector of module colours.
#' @param linewidth Line width for smoothed trajectories.
#' @param axis_text_size Text size for axis labels and ticks.
#' @param legend_text_size Text size for legend labels and title.
#' @param label_size Text size for module labels.
#' @param show_labels Logical indicating whether module labels are shown. The
#' labels will be located at the maximum point of the smoothed trajectory.
#'
#' @return A ggplot object.
#' @export
plot_module_trends_over_pseudotime <- function(expression_list,
    pseudotime,
    module_colours = NULL,
    linewidth = 1.5,
    axis_text_size = 10,
    legend_text_size = 10,
    label_size = 10,
    show_labels = TRUE
) {
    if (inherits(expression_list, "list")) {
        module_names <- names(expression_list)
        df <- do.call(rbind, lapply(module_names, function(module_name) {
            data.frame(
                pseudotime = pseudotime,
                expression = expression_list[[module_name]],
                module = module_name
            )
        }))
    } else if (inherits(expression_list, c("matrix", "data.frame"))) {
        module_names <- colnames(expression_list)
        df <- reshape2::melt(as.data.frame(expression_list), varnames = "module", value.name = "expression")

        df$pseudotime <- rep(pseudotime, each = ncol(expression_list))
    } else {
        stop("Unsupported input type for expression_list. Please provide a list of vectors or a matrix/data.frame.")
    }

    if (isFALSE(all(module_names %in% names(module_colours)))) {
        module_colours <- NULL
    }
    if (is.null(module_colours)) {
        module_colours <- stats::setNames(qualpalr::qualpal(length(module_names), list(h = c(0, 360), s = c(0.2, 0.8), l = c(0.4, 0.65)))$hex, module_names)
    }
    ncells <- length(pseudotime)

    # get the smooth approximation
    smooth_approx <- df %>%
        dplyr::group_by(.data$module) %>%
        dplyr::reframe(pseudotime = seq(min(.data$pseudotime), max(.data$pseudotime), length.out = ncells)) %>%
        dplyr::group_by(.data$module)
    if (ncells < 7000) {
        # use loess
        smooth_approx <- smooth_approx %>%
            dplyr::group_modify(~ {
                module_name <- .y$module
                module_data <- df[df$module == module_name, ]
                loess_fit <- stats::loess(expression ~ pseudotime, data = module_data)
                .x$expression <- stats::predict(loess_fit, newdata = .x)
                return(.x)
            })
    } else {
        # use gam
        smooth_approx <- smooth_approx %>%
            dplyr::group_modify(~ {
                module_name <- .y$module
                module_data <- df[df$module == module_name, ]
                gam_fit <- mgcv::gam(expression ~ s(pseudotime), data = module_data)
                .x$expression <- stats::predict(gam_fit, newdata = .x)
                return(.x)
            })
    }

    ggplot_obj <- ggplot2::ggplot(df, ggplot2::aes(x = .data$pseudotime, y = .data$expression, color = .data$module, group = .data$module)) +
        ggplot2::geom_line(data = smooth_approx, linewidth = linewidth) +
        ggplot2::scale_color_manual(values = module_colours) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Pseudotime", y = "Module Expression", color = "Module") +
        ggplot2::theme(
            legend.text = ggplot2::element_text(size = legend_text_size),
            legend.title = ggplot2::element_text(size = legend_text_size),
            axis.text = ggplot2::element_text(size = axis_text_size),
            axis.title = ggplot2::element_text(size = axis_text_size),
            plot.title = ggplot2::element_text(size = label_size)
        )
    
    if (!show_labels) {
        return(ggplot_obj)
    }

    text_pos <- smooth_approx %>%
        dplyr::group_by(.data$module) %>%
        dplyr::slice_max(order_by = .data$expression, n = 1) %>%
        dplyr::ungroup() %>%
        dplyr::select(.data$module, pseudotime_at_max = .data$pseudotime, max_expr = .data$expression) %>%
        dplyr::mutate(label = .data$module)
    
    min_expr <- min(smooth_approx$expression, na.rm = TRUE)
    max_expr <- max(smooth_approx$expression, na.rm = TRUE)

    return(
        ggplot_obj +
            ggplot2::geom_text(data = text_pos, ggplot2::aes(x = .data$pseudotime_at_max, y = .data$max_expr, label = .data$label), size = label_size, vjust = -0.5) +
            ggplot2::coord_cartesian(ylim = c(min_expr * 0.97, max_expr * 1.1))
    )
}

#' Plot Module Pseudobulk Expression
#'
#' @description Computes per-metadata pseudobulk module expression and plots a
#' bubble heatmap.
#'
#' @param module_expression_summary A module-by-cell matrix/data frame, or list
#' of module vectors.
#' @param cell_metadata Metadata vector used to define pseudobulk groups.
#' @param mtd_order Desired display order for metadata groups.
#' @param pseudobulk_function Function used to aggregate expression per group.
#' @param scale Logical indicating whether pseudobulk values are z-scored per
#' module.
#' @param cap_value Maximum absolute value used for color limits.
#' @param axis_text_size Text size for axes.
#' @param legend_text_size Text size for legends.
#' @param point_range Range for point-size scaling.
#' @param continuous_colours Optional color palette for expression values.
#'
#' @return A ggplot object.
#' @export
plot_module_pseudobulk_expression <- function(
    module_expression_summary,
    cell_metadata,
    mtd_order,
    pseudobulk_function = mean,
    scale = TRUE,
    cap_value = 2,
    axis_text_size = 10,
    legend_text_size = 10,
    point_range = c(0.1, 1),
    continuous_colours = NULL
) {
    is_na_mask <- !is.na(cell_metadata)
    mtd_order <- mtd_order[!is.na(mtd_order)]

    if (!inherits(module_expression_summary, "list")) {
        module_names <- colnames(module_expression_summary)
        module_expression_summary <- lapply(module_names, function(module_name) {
            module_expression_summary[, module_name]
        })
        names(module_expression_summary) <- module_names
        return(plot_module_pseudobulk_expression(module_expression_summary, cell_metadata, mtd_order, scale, cap_value, axis_text_size, legend_text_size, point_range))
    }

    module_names <- names(module_expression_summary)
    ncells <- sum(is_na_mask)
    df <- do.call(rbind, lapply(module_names, function(module_name) {
        temp_df <- data.frame(
            expression = module_expression_summary[[module_name]][is_na_mask],
            mtd = as.character(cell_metadata[is_na_mask])
        )

        temp_df %>%
            dplyr::group_by(.data$mtd) %>%
            dplyr::summarise(
                psdbk_expr = pseudobulk_function(.data$expression, na.rm = TRUE),
                perc_expressed = sum(.data$expression > 0, na.rm = TRUE) / ncells,
                module = module_name
            )
    }))
    df$mtd <- factor(df$mtd, levels = mtd_order)
    df$module <- factor(df$module, levels = module_names)

    if (is.null(continuous_colours)) {
        continuous_colours = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
            "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
    }

    colour_limits <- c(0, cap_value)
    if (scale) {
        # scale values to mimick t(scale(t(psd_matrix)))
        df <- df %>%
            dplyr::group_by(.data$module) %>%
            dplyr::mutate(psdbk_expr = scale(.data$psdbk_expr)) %>%
            dplyr::ungroup()
        colour_limits <- c(-cap_value, cap_value)
    }

    return(
        ggplot2::ggplot(df, ggplot2::aes(y = .data$mtd, x = .data$module, color = .data$psdbk_expr, size = .data$perc_expressed)) +
            ggplot2::geom_point() +
            ggplot2::theme_classic() +
            ggplot2::labs(x = "Module", y = "Metadata Group", color = "Pseudobulk\nExpression", size = "Percentage\nExpressed") +
            ggplot2::scale_color_gradientn(colors = continuous_colours, limits = colour_limits) +
            ggplot2::theme(
                legend.text = ggplot2::element_text(size = legend_text_size),
                legend.title = ggplot2::element_text(size = legend_text_size),
                axis.text = ggplot2::element_text(size = axis_text_size),
                axis.title = ggplot2::element_text(size = axis_text_size)
            )
    )
}
