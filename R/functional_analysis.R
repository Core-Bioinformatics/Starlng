
#### ENRICHMENT ANALYSIS ####

#' Filter Enrichment Results
#'
#' @description Filters enrichment results by p-value and optional
#' intersection-size bounds. It assumes that the results are from a run of
#' `gprofiler2::gost`.
#'
#' @param enrichment_result A data frame or list containing enrichment results.
#' If a list is provided, each element should correspond to a module. The
#' elements will be automatically combined into a one data frame.
#' @param p_value_threshold Numeric threshold for adjusted p-values. Defaults to
#' 0.05.
#' @param intersection_size_upper Optional upper bound for intersection size.
#' @param intersection_size_lower Optional lower bound for intersection size.
#'
#' @return A filtered enrichment data frame, or NULL if no rows pass filters.
#' @export
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

#' Plot Top Enriched Terms
#'
#' @description Builds a dot plot with the top enriched terms per module.
#'
#' @param enrichment_result A data frame produced by enrichment analysis.
#' @param top_n Number of top terms to keep per module. The ranking is determined
#' by p-value. Defaults to 2.
#' @param colour_column Optional column name used to map point size. We suggest
#' using "intersection_size".
#' @param point_size_range Numeric vector of length 2 defining point size range.
#' @param font_size Base font size used in the plot theme.
#'
#' @return A ggplot object, or NULL if the input is empty.
#' @export
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

#### TF ANALYSIS ####

#' Normalize Transcription Factor Names
#'
#' @description Normalizes transcription-factor names by replacing common Greek
#' words to the alphabet lettersand removing dashes.
#'
#' @param tf A character vector with transcription-factor names.
#'
#' @return A character vector with normalized names.
#' @keywords internal
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

#' Get Enriched Transcription Factors
#'
#' @description Performs transcription-factor enrichment against the TRANSFAC
#' database using `gprofiler2`, groups the transcription factors by their
#' name and identifies the associated genes.
#'
#' @param gene_list A character vector of genes, or a named list of gene vectors.
#' @param organism Organism code used by gprofiler2.
#' @param p_value_threshold Numeric p-value threshold used for filtering.
#' @param correction_method Multiple-testing correction method. Defaults to
#' 'fdr'.
#' @param bg_genes Optional background gene universe. The domain scope of the 
#' analysis is automatically adjusted based on whether this is provided or not.
#' @param ... Additional arguments passed to gprofiler2::gost.
#'
#' @return A list with enrichment table and TF-associated genes, or NULL if no
#' significant enrichment is found.
#' @export
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

#' Summarize Transcription Factor Statistics
#'
#' @description Converts transcription-factor enrichment output to a tabular
#' summary. The summary includes, for each TF, the number of associated genes
#' and optionally the list of associated genes as a comma-separated string.
#'
#' @param transcription_factors Output from get_transcription_factors.
#' @param module_name Optional module name used when processing one module.
#' @param include_intersection_set Logical indicating whether to include
#' associated genes as a comma-separated string.
#'
#' @return A data frame with TF statistics, or NULL if no transcription factors
#' are available.
#' @export
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

#' Add Hub-Gene Counts to TF Statistics
#'
#' @description Adds, for each TF, the number of associated genes that are hub
#' genes.
#'
#' @param existing_stats A TF statistics data frame.
#' @param hub_genes A data frame containing hub genes and their module labels.
#' The name of the columns should be "gene" and "module".
#'
#' @return A TF statistics data frame with an added n_hub_genes column when
#' possible.
#' @export
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

#' Build TF-Gene Network Graph
#'
#' @description Builds an igraph object representing TF-TF and TF-gene
#' relationships. The graph is built based on the associations done from the
#' TF enrichment analysis. TF-TF relationships are derived from neighbourhood
#' of the modules along the graph trajectory.
#'
#' @param tf_stats A TF statistics data frame including module, tf and genes.
#' @param module_adjacency Optional module adjacency matrix.
#' @param top_n_factors Number of top TFs to retain per module.
#' @param hub_genes Optional data frame of hub genes used for annotation.
#'
#' @return An igraph object, or invisible(NULL) when no valid TFs are available.
#' @export
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
    common_modules <- unique(tf_stats$module)
    if (!is.null(module_adjacency)) {
        modules_from_adj <- colnames(module_adjacency)
        common_modules <- intersect(common_modules, modules_from_adj)
        module_adjacency <- module_adjacency[common_modules, common_modules, drop = FALSE]
    }

    tf_stats <- tf_stats %>%
        dplyr::filter(.data$module %in% common_modules)

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

#' Plot TF-Gene Network with ggraph
#'
#' @description Plots a TF-gene igraph object. The TFs and the module gene
#' hubs are highlighted using different shapes and sized of the points. The
#' edges between TFs and genes are coloured differently from the edges between TFs.
#'
#' @param module_tf_g An igraph object produced by get_tf_gene_network.
#' @param tf_colour Colour used for TF nodes and TF-TF edges.
#' @param module_colours Named vector of module colours.
#' @param tf_shape Shape used for TF nodes.
#' @param hub_shape Shape used for hub-gene nodes.
#' @param gene_shape Shape used for non-hub gene nodes.
#' @param node_size Base node-size scaling factor.
#' @param edge_width Base edge width.
#' @param edge_colour Colour used for TF-gene edges.
#' @param label_size Text size of node labels.
#' @param axis_text_size Base text size for title and legends.
#' @param exclude_non_hub_genes Logical indicating whether non-hub genes are
#' included in the plotting.
#'
#' @return A ggplot object, or invisible(NULL) when plotting dependencies are
#' unavailable.
#' @export
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

#' Plot TF Bubble Plot
#'
#' @description Creates a bubble plot of top TFs per module.
#'
#' @param tf_stats A TF statistics data frame.
#' @param n_top Number of top TFs to keep per module.
#' @param cap_value_intersection Optional cap for n_genes.
#' @param cap_value_hub Optional cap for n_hub_genes.
#' @param font_size Base font size for plot labels and legends.
#' @param point_size_range Numeric vector of length 2 defining point size range.
#'
#' @return A ggplot object.
#' @export
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
