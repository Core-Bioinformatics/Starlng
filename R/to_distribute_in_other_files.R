
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


#### get the pseudotime trajectory  ####
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

#### module connectivity ####
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


#### gene umap ####

# gene umap
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


#### general enrichment (dotplot of top enriched terms) ####
filter_enrichment_results <- function(
    enrichment_result,
    p_value_threshold = 0.05,
    intersection_size_upper = NULL,
    intersection_size_lower = NULL
) {
    if (is.null(enrichment_result)) {
        return(NULL)
    }

    if (inherits(enrichment_result, "list")) {
        if ("result" %in% names(enrichment_result)) {
            return(filter_enrichment_results(enrichment_result$result, p_value_threshold, intersection_size_upper, intersection_size_lower))
        }

        module_names <- names(enrichment_result)
        if (is.null(module_names)) {
            module_names <- as.character(seq_along(enrichment_result))
        }
        names(enrichment_result) <- module_names
        unified_enrichment <- do.call(rbind, lapply(module_names, function(module_name) {
            enr <- enrichment_result[[module_name]]
            if (is.null(enr)) {
                return(NULL)
            }

            if ("result" %in% names(enr)) {
                enr <- enr$result
            }
            enr$module <- module_name
            return(enr)
        }))
        return(filter_enrichment_results(unified_enrichment, p_value_threshold, intersection_size_upper, intersection_size_lower))
    }

    if (nrow(enrichment_result) == 0) {
        return(NULL)
    }

    filtered_result <- enrichment_result %>%
        dplyr::filter(.data$p_value < p_value_threshold) %>%
        dplyr::arrange(.data$p_value)

    if (!is.null(intersection_size_upper)) {
        filtered_result <- filtered_result %>% dplyr::filter(.data$intersection_size <= intersection_size_upper)
    }

    if (!is.null(intersection_size_lower)) {
        filtered_result <- filtered_result %>% dplyr::filter(.data$intersection_size >= intersection_size_lower)
    }

    if (nrow(filtered_result) == 0) {
        return(NULL)
    }

    return(filtered_result)
}

plot_enrichment_top_terms <- function(
    enrichment_result,
    top_n = 2,
    colour_column = NULL,
    point_size_range = c(3, 10),
    font_size = 10
) {
    if (is.null(enrichment_result) || nrow(enrichment_result) == 0) {
        return(NULL)
    }
    top_terms <- enrichment_result %>%
        dplyr::group_by(.data$module) %>%
        dplyr::arrange(.data$p_value) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::pull(.data$term_name)
    enrichment_result <- enrichment_result %>%
        dplyr::filter(.data$term_name %in% top_terms) %>%
        dplyr::arrange(.data$module, dplyr::desc(.data$p_value))

    if (is.null(colour_column)) {
        enrichment_result$colour <- 1
    } else {
        enrichment_result$colour <- enrichment_result[[colour_column]]
    }

    # order the term names to follow the ordering of the module
    # terms should have 5 words per line
    enrichment_result$term_name <- sapply(enrichment_result$term_name, function(term) {
        words <- unlist(strsplit(term, " "))
        if (length(words) <= 4) {
            return(term)
        }
        new_term <- ""
        for (i in seq(1, length(words), by = 4)) {
            new_term <- paste(new_term, paste(words[i:min(i+3, length(words))], collapse = " "), "\n")
        }
        return(new_term)
    })
    enrichment_result$term_name <- factor(enrichment_result$term_name, levels = unique(enrichment_result$term_name))



    gplot_obj <- ggplot2::ggplot(enrichment_result, ggplot2::aes(x = .data$module, y = .data$term_name, color = -log10(.data$p_value), size = .data$colour)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Module", y = "Term", title = "Top Enriched Terms") +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = font_size),
            axis.text.y = ggplot2::element_text(size = font_size),
            axis.title = ggplot2::element_text(size = font_size * 1.2),
            legend.title = ggplot2::element_text(size = font_size),
            legend.text = ggplot2::element_text(size = font_size),
            plot.title = ggplot2::element_text(size = font_size * 1.4)
        )

    if (is.null(colour_column)) {
        gplot_obj <- gplot_obj +
            ggplot2::scale_color_continuous(name = "-log10(p-value)") +
            ggplot2::scale_size_continuous(name = NULL, labels = NULL, breaks = NULL)
    } else {
        gplot_obj <- gplot_obj +
            ggplot2::scale_color_continuous(name = "-log10(p-value)") +
            ggplot2::scale_size_continuous(name = colour_column, range = point_size_range)
    }
    return(gplot_obj)
}

#### tf functions ####

clean_tf_names <- function(tf) {
    tf <- toupper(tf)
    tf_clean <- gsub("ALPHA", "A", tf)
    tf_clean <- gsub("BETA", "B", tf_clean)
    tf_clean <- gsub("GAMMA", "G", tf_clean)
    tf_clean <- gsub("DELTA", "D", tf_clean)
    tf_clean <- gsub("EPSILON", "E", tf_clean)
    tf_clean <- gsub("KAPPA", "K", tf_clean)
    tf_clean <- gsub("-", "", tf_clean)
    return(tf_clean)
}

get_transcription_factors <- function(
    gene_list,
    organism = "hsapiens",
    p_value_threshold = 0.05,
    correction_method = "fdr",
    bg_genes = NULL,
    ...
) {
    if (inherits(gene_list, "list")) {
        tf_list <- lapply(gene_list, function(genes) {
            get_transcription_factors(genes, organism, p_value_threshold, correction_method, bg_genes, ...)
        })
        for (module_name in names(tf_list)) {
            if (is.null(tf_list[[module_name]]) || is.null(tf_list[[module_name]]$tf_associated_genes) || length(tf_list[[module_name]]$tf_associated_genes) == 0) {
                tf_list[[module_name]] <- NULL
                next
            }
        }

        tf_names <- lapply(tf_list, function(x) {
            current_tfs <- names(x$tf_associated_genes)
            current_tfs <- sapply(current_tfs, clean_tf_names)
        })
        tf_mapping <- list()

        for (module_name in names(tf_names)) {
            current_tfs <- tf_names[[module_name]]
            for (i in seq_along(current_tfs)) {
                clean_name <- current_tfs[i]
                if (clean_name %in% names(tf_mapping)) {
                    names(current_tfs)[i] <- tf_mapping[[clean_name]]
                } else {
                    tf_mapping[[clean_name]] <- names(current_tfs)[i]
                }
            }
            names(tf_list[[module_name]]$tf_associated_genes) <- names(current_tfs)
        }
        return(tf_list)
    }

    result <- list()

    if (is.null(bg_genes)) {
        annotation_domain <- "annotated"
    } else {
        annotation_domain <- "custom_annotated"
    }

    tf_enrichment <- suppressMessages(gprofiler2::gost(
        query = gene_list,
        organism = organism,
        sources = c("TF"),
        significant = FALSE,
        domain_scope = annotation_domain,
        custom_bg = bg_genes,
        correction_method = correction_method,
        evcodes = TRUE,
        ...
    ))

    if (is.null(tf_enrichment)) {
        return(NULL)
    }
    tf_enrichment <- tf_enrichment$result %>% 
        dplyr::filter(.data$p_value < p_value_threshold) %>%
        dplyr::arrange(.data$p_value)
    result$data_frame <- tf_enrichment

    tf_enrichment$tf_name <- sapply(tf_enrichment$term_name, function(term_id) {
        strsplit(strsplit(term_id, ";")[[1]][1], "Factor: ")[[1]][2]
    })
    tf_associated_genes <- list()
    for (tf in unique(tf_enrichment$tf_name)) {
        associated_genes <- unlist(strsplit(tf_enrichment$intersection[tf_enrichment$tf_name == tf], ","))
        tf_associated_genes[[tf]] <- associated_genes
    }

    tf_clean_names <- sapply(names(tf_associated_genes), clean_tf_names)
    names(tf_clean_names) <- names(tf_associated_genes)

    for (tf in names(tf_clean_names)) {
        duplicates <- setdiff(names(tf_clean_names[tf_clean_names == tf_clean_names[tf]]), tf)
        duplicates <- intersect(duplicates, names(tf_associated_genes))

        if (length(duplicates) == 0) {
            next
        }

        for (duplicate in duplicates) {
            tf_associated_genes[[tf]] <- unique(c(tf_associated_genes[[tf]], tf_associated_genes[[duplicate]]))
            tf_associated_genes[[duplicate]] <- NULL
        }
    }
    for (tf in names(tf_associated_genes)) {
        tf_associated_genes[[tf]] <- unique(tf_associated_genes[[tf]])
    }
    result$tf_associated_genes <- tf_associated_genes

    return(result)
}

get_tf_stats <- function(transcription_factors, module_name = NULL, include_intersection_set = FALSE) {
    if (is.null(transcription_factors)) {
        return(NULL)
    }
    if ("tf_associated_genes" %in% names(transcription_factors)) {
        if( is.null(transcription_factors$tf_associated_genes) || length(transcription_factors$tf_associated_genes) == 0) {
            return(NULL)
        }
        current_tf_names <- names(transcription_factors$tf_associated_genes)
        tf_stats <- do.call(rbind, lapply(seq_along(transcription_factors$tf_associated_genes), function(i) {
            genes <- transcription_factors$tf_associated_genes[[i]]
            df <- data.frame(module = module_name, tf = current_tf_names[i], n_genes = length(genes))
            if (include_intersection_set) {
                df$genes <- paste(genes, collapse = ",")
            }
            return(df)
        }))
        return(tf_stats)
    }

    return(do.call(rbind, lapply(names(transcription_factors), function(module) {
        get_tf_stats(transcription_factors[[module]], module_name = module, include_intersection_set = include_intersection_set)
    })))
}

add_tf_hub_stats <- function(existing_stats, hub_genes) {
    if (is.null(hub_genes)) {
        return(existing_stats)
    }
    if (isFALSE("genes" %in% colnames(existing_stats))) {
        return(existing_stats)
    }
    hub_genes <- split(hub_genes$gene, hub_genes$module)
    n_hub_genes <- rep(0, nrow(existing_stats))
    for (i in seq_len(nrow(existing_stats))) {
        module_i <- existing_stats$module[i]
        tf_genes <- unlist(strsplit(existing_stats$genes[i], ","))
        n_hub_genes <- sum(tf_genes %in% hub_genes[[module_i]])
        existing_stats$n_hub_genes[i] <- n_hub_genes
    }
    # put the intersection column at the end
    existing_stats <- existing_stats[, c(setdiff(colnames(existing_stats), "genes"), "genes")]
    return(existing_stats)
}

get_tf_gene_network <- function(
    # transcription_factors,
    tf_stats,
    module_adjacency = NULL,
    top_n_factors = 3,
    hub_genes = NULL
) {
    if (is.null(tf_stats) || nrow(tf_stats) == 0) {
        warning("No transcription factors available for plotting.")
        return(invisible(NULL))
    }
    modules_from_tfs <- unique(tf_stats$module)
    modules_from_adj <- colnames(module_adjacency)
    common_modules <- intersect(modules_from_tfs, modules_from_adj)

    tf_stats <- tf_stats %>%
        dplyr::filter(.data$module %in% common_modules)
    module_adjacency <- module_adjacency[common_modules, common_modules, drop = FALSE]

    top_tfs <- tf_stats %>%
        dplyr::group_by(.data$module) %>%
        dplyr::arrange(dplyr::desc(.data$n_genes)) %>%
        dplyr::slice_head(n = top_n_factors) %>%
        dplyr::ungroup() %>%
        dplyr::pull(.data$tf) %>%
        unique()

    if (length(top_tfs) == 0) {
        warning("No top TFs passed filtering; skipping TF network plot.")
        return(invisible(NULL))
    }

    tf_adjacency <- matrix(0, nrow = length(top_tfs), ncol = length(top_tfs), dimnames = list(top_tfs, top_tfs))
    for (i in seq_along(top_tfs)) {
        tf_i <- top_tfs[i]
        modules_i <- tf_stats$module[tf_stats$tf == tf_i]
        for (j in seq(from = i, to = length(top_tfs))) {
            tf_j <- top_tfs[j]
            modules_j <- tf_stats$module[tf_stats$tf == tf_j]
            if (i == j) {
                next
            }

            are_linked <- any(match(modules_i, modules_j, nomatch = 0) > 0)
            if (!is.null(module_adjacency)) {
                are_linked <- are_linked || any(module_adjacency[modules_i, modules_j, drop = FALSE] == 1)
            }
            if (are_linked) {
                tf_adjacency[tf_i, tf_j] <- 1
                tf_adjacency[tf_j, tf_i] <- 1
            }
        }
    }

    # tf - gene adjacency
    # used_modules <- unique(tf_stats$module)
    tf_stats <- tf_stats %>%
        dplyr::filter(.data$tf %in% top_tfs)
    tf_gene_edges <- do.call(rbind, lapply(seq_len(nrow(tf_stats)), function(i) {
        module_i <- tf_stats$module[i]
        tf <- tf_stats$tf[i]
        genes <- unlist(strsplit(tf_stats$genes[i], ","))
        if (length(genes) == 0) {
            return(NULL)
        }
        return(data.frame(from = tf, to = genes, type = "tf-gene", module = module_i))
    }))
    
    hub_gene_names <- character(0)
    if (!is.null(hub_genes) && "gene" %in% colnames(hub_genes)) {
        hub_gene_names <- hub_genes$gene
    }
    tf_gene_edges$is_hub_gene <- tf_gene_edges$to %in% hub_gene_names

    gene_module_map <- tapply(tf_gene_edges$module, tf_gene_edges$to, function(x) as.character(x[1]))
    # create the graph
    g_tf <- igraph::graph_from_adjacency_matrix(tf_adjacency, mode = "undirected")
    igraph::V(g_tf)$type <- "tf"

    gene_nodes <- unique(tf_gene_edges$to)
    gene_nodes <- setdiff(gene_nodes, igraph::V(g_tf)$name)
    if (length(gene_nodes) > 0) {
        g_tf <- igraph::add_vertices(g_tf, nv = length(gene_nodes), name = gene_nodes, type = rep("gene", length(gene_nodes)))
    }

    g_tf <- igraph::add_edges(
        g_tf,
        t(as.matrix(tf_gene_edges[, c("from", "to")])),
        type = rep("tf-gene", nrow(tf_gene_edges)),
        is_hub_gene = tf_gene_edges$is_hub_gene,
        module = tf_gene_edges$module
    )

    if (igraph::vcount(g_tf) == 0) {
        warning("Graph has no vertices; skipping plot.")
        return(invisible(g_tf))
    }

    igraph::V(g_tf)$is_hub_gene <- igraph::V(g_tf)$type == "gene" & igraph::V(g_tf)$name %in% hub_gene_names
    igraph::V(g_tf)$module <- ifelse(igraph::V(g_tf)$type == "gene", gene_module_map[igraph::V(g_tf)$name], "")

    return(g_tf)
}

plot_module_tfs_ggraph <- function(
    module_tf_g,
    tf_colour = "black",
    module_colours = NULL,
    tf_shape = 17,
    hub_shape = 23,
    gene_shape = 16,
    node_size = 4,
    edge_width = 0.5,
    edge_colour = "#bbbbbb",
    label_size = 3,
    axis_text_size = 8,
    exclude_non_hub_genes = FALSE
) {
    gene_node_size <- node_size
    hub_node_size <- node_size * 1.7
    tf_node_size <- node_size * 2.3
    has_ggraph_or_ggplot <- requireNamespace("ggraph", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)
    if (!has_ggraph_or_ggplot) {
        warning("ggraph and ggplot2 packages are required for plotting the TF-gene network. Please install them to see the plot.")
        return(invisible(NULL))
    }
    has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
    if (is.null(module_colours)) {
        unique_modules <- setdiff(unique(igraph::V(module_tf_g)$module), "")
        module_colours <- stats::setNames(grDevices::hcl.colors(length(unique_modules), palette = "Dark 3"), unique_modules)
    }

    vertex_colour_key <- ifelse(igraph::V(module_tf_g)$type == "tf", "TF", as.character(igraph::V(module_tf_g)$module))
    vertex_colour_key[is.na(vertex_colour_key)] <- "Unknown"
    colour_values <- c("TF" = tf_colour, module_colours)
    if (any(vertex_colour_key == "Unknown")) {
        colour_values <- c(colour_values, "Unknown" = "grey70")
    }

    vertex_shapes <- ifelse(
        igraph::V(module_tf_g)$type == "tf",
        "TF",
        ifelse(igraph::V(module_tf_g)$is_hub_gene, "Hub gene", "Gene")
    )
    vertex_sizes <- ifelse(igraph::V(module_tf_g)$type == "tf", tf_node_size, ifelse(igraph::V(module_tf_g)$is_hub_gene, hub_node_size, gene_node_size))
    vertex_labels <- ifelse(igraph::V(module_tf_g)$type == "tf" | igraph::V(module_tf_g)$is_hub_gene, igraph::V(module_tf_g)$name, "")

    edge_type <- igraph::E(module_tf_g)$type
    edge_type[is.na(edge_type)] <- "tf-tf"

    igraph::V(module_tf_g)$vertex_colour <- vertex_colour_key
    igraph::V(module_tf_g)$vertex_shape <- vertex_shapes
    igraph::V(module_tf_g)$vertex_size <- vertex_sizes
    igraph::V(module_tf_g)$vertex_label <- vertex_labels
    igraph::E(module_tf_g)$edge_type <- edge_type
    igraph::E(module_tf_g)$edge_width <- ifelse(edge_type == "tf-tf", edge_width * 4, edge_width)

    if (exclude_non_hub_genes) {
        module_tf_g <- igraph::delete_vertices(module_tf_g, igraph::V(module_tf_g)[.data$type == "gene" & !.data$is_hub_gene])
    }

    layout_df <- ggraph::create_layout(module_tf_g, layout = "stress")
    layout_df$x <- layout_df$x * 1.5
    layout_df$y <- layout_df$y * 1.5

    gplot_obj <- ggraph::ggraph(layout_df) +
        ggraph::geom_edge_link(ggplot2::aes(colour = .data$edge_type, edge_width = .data$edge_width), alpha = 0.55) +
        ggraph::geom_node_point(ggplot2::aes(colour= .data$vertex_colour, size = .data$vertex_size, shape = .data$vertex_shape)) +
        ggraph::scale_edge_colour_manual(values = c("tf-gene" = edge_colour, "tf-tf" = tf_colour), name = "Edge type") +
        ggraph::scale_edge_width(range = c(edge_width, edge_width * 2), guide = "none") +
        ggplot2::scale_shape_manual(values = c("TF" = tf_shape, "Gene" = gene_shape, "Hub gene" = hub_shape), name = "Node shape") +
        ggplot2::scale_colour_manual(values = colour_values, name = "Module") +
        ggplot2::scale_size_identity(guide = "none") +
        ggplot2::coord_equal() +
        ggplot2::theme_void() +
        ggplot2::ggtitle("TF-Gene Regulatory Network") +
        ggplot2::theme(
            legend.text = ggplot2::element_text(size = axis_text_size),
            legend.title = ggplot2::element_text(size = axis_text_size * 1.2),
            plot.title = ggplot2::element_text(size = axis_text_size * 1.4)
        )

    if (has_ggrepel) {
        gplot_obj <- gplot_obj +
            ggrepel::geom_label_repel(
                data = layout_df |>
                    dplyr::filter(.data$vertex_label != ""),
                ggplot2::aes(x = .data$x, y = .data$y, label = .data$vertex_label, colour = .data$vertex_colour),
                size = label_size,
                fill = ggplot2::alpha("white", 0.75),
                label.size = 0.15,
                fontface = "bold",
                max.overlaps = Inf,
                box.padding = 0.35,
                point.padding = 0.25,
                min.segment.length = 0
            )
    } else {
        gplot_obj <- gplot_obj +
            ggraph::geom_node_text(ggplot2::aes(label = .data$vertex_label, colour = .data$vertex_colour), size = 3, repel = TRUE, check_overlap = TRUE, max.overlaps = Inf) 
    }

    return(gplot_obj)
}

plot_module_tfs_bubbleplot <- function(
    tf_stats,
    n_top = 10,
    cap_value_intersection = NULL,
    cap_value_hub = NULL,
    font_size = 10,
    point_size_range = c(3, 10)
) {
    top_tfs <- tf_stats %>%
        dplyr::group_by(.data$module) %>%
        dplyr::arrange(dplyr::desc(.data$n_genes)) %>%
        dplyr::slice_head(n = n_top) %>%
        dplyr::ungroup() %>%
        dplyr::pull(.data$tf) %>%
        unique()
    tf_stats <- tf_stats %>%
        dplyr::filter(.data$tf %in% top_tfs) %>%
        dplyr::arrange(.data$module)
    tf_stats$tf <- factor(tf_stats$tf, levels = unique(tf_stats$tf))
    if (!is.null(cap_value_intersection)) {
        tf_stats$n_genes <- pmin(tf_stats$n_genes, cap_value_intersection)
    }
    if (isFALSE("n_hub_genes" %in% colnames(tf_stats))) {
        tf_stats$n_hub_genes <- 0
    }
    if (!is.null(cap_value_hub)) {
        tf_stats$n_hub_genes <- pmin(tf_stats$n_hub_genes, cap_value_hub)
    }

    ggplot2::ggplot(tf_stats, ggplot2::aes(x = .data$module, y = .data$tf, size = .data$n_genes, fill = .data$n_hub_genes)) +
        ggplot2::geom_point(shape = 21, colour = "black") +
        ggplot2::scale_fill_gradient(name = "Number of\nhub genes", low = "white", high = "steelblue") +
        ggplot2::scale_size_continuous(name = "Number of\nassociated genes", range = point_size_range) +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Module", y = "Transcription Factor", title = "TF-Gene Associations by Module") +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = font_size),
            axis.text.y = ggplot2::element_text(size = font_size),
            legend.title = ggplot2::element_text(size = font_size),
            legend.text = ggplot2::element_text(size = font_size),
            plot.title = ggplot2::element_text(size = font_size * 1.4),
            axis.title = ggplot2::element_text(size = font_size * 1.2)
        )
}
