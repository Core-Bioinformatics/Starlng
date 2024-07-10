#' @importFrom dplyr %>% 
#' @importFrom rlang .data
 
select_cells_by_gene_expr <- function(expression_matrix,
                                      genes,
                                      thresh_percentile = 0.25,
                                      thresh_value = 0,
                                      n_coexpressed_thresh = length(genes)) {
    genes <- intersect(genes, rownames(expression_matrix))

    if (length(genes) == 0) {
        stop("No genes found in the expression matrix")
    }

    expression_matrix <- expression_matrix[genes, ]
    n_coexpressed <- max(min(n_coexpressed, nrow(expression_matrix)), 1)

    if (thresh_percentile > 0) {
        thresh_value <- sapply(genes, function(gene) {
            quantile(sum(expression_matrix[gene, ]), probs = thresh_percentile)
        })
    } else {
        thresh_value <- rep(thresh_value[1], length(genes))
    }

    mask_expressed <- expression_matrix > thresh_value

    index_cells <- which(colSums(mask_expressed) >= n_coexpressed_thresh)

    if (is.null(colnames(expression_matrix))) {
        return(index_cells)
    }

    return(colnames(expression_matrix)[index_cells])
}

select_cells_by_metadata <- function(metadata_df, metadata_combinations) {
    mtd_names <- intersect(names(metadata_combinations), colnames(metadata_df))

    if (length(mtd_names) == 0) {
        stop("No metadata columns found in the metadata dataframe")
    }

    if (is.null(rownames(metadata_df))) {
        rownames(metadata_df) <- seq_len(nrow(metadata_df))
    }
    

    for (mtd_name in mtd_names) {
        if (!(all(metadata_combinations[[mtd_name]] %in% unique(metadata_df[[mtd_name]])))) {
            next
        }

        metadata_df <- metadata_df %>% dplyr::filter(.data[[mtd_name]] %in% metadata_combinations[[mtd_name]])
    }

    return(rownames(metadata_df))
}

remove_outlier_cells <- function(cell_names, umap_emb, percentile_threshold = 0.75) {
    if (is.null(rownames(umap_emb))) {
        rownames(umap_emb) <- seq_len(nrow(umap_emb))
    }

    if (length(intersect(cell_names, rownames(umap_emb))) == 0) {
        stop("The cells were not found in the UMAP embedding")
    }

    umap_emb <- umap_emb[cell_names, ]

    gmedian <- Gmedian::Gmedian(umap_emb)
    distance_cells_gmedian <- (
        (umap_emb[, 1] - gmedian[1, 1]) ^ 2 +
        (umap_emb[, 2] - gmedian[1, 2]) ^ 2
    ) ^ 0.5

    threshold_distance <- quantile(distance_cells_gmedian, probs = percentile_threshold)

    return(cell_names[distance_cells_gmedian <= threshold_distance])
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

    mon_obj <- mon_obj[, all_cells_between]
    gc()

    # update the principal graph
    monocle3::principal_graph(mon_obj)$UMAP <- igraph::induced_subgraph(
        monocle3::principal_graph(mon_obj)$UMAP,
        all_cells_between
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
                               learn_graph_controls = NULL,
                               use_partition = FALSE,
                               use_closed_loops = FALSE,
                               verbose = FALSE,
                               metadata_for_nodes = NULL) {
    cal_ncenter <- function(num_cell_communities, ncells, nodes_per_log10_cells=15) {
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

    if (use_partition) {
        learn_graph_controls$ncenter <- NULL

        if (!is.null(metadata_for_nodes)) {
            mon_obj@clusters@listData$UMAP$partitions <- setNames(factor(as.numeric(mon_obj@colData[[metadata_for_nodes]])), colnames(mon_obj))
            levels(mon_obj@clusters@listData$UMAP$partitions) <- seq_len(nlevels(mon_obj@clusters@listData$UMAP$partitions))
        }
    }

    return(monocle3::learn_graph(
        cds = mon_obj,
        close_loop = use_closed_loops,
        use_partition = use_partition,
        learn_graph_control = learn_graph_controls,
        verbose = verbose
    ))
}

update_mononcle_partition <- function(mon_obj, new_partition) {
    new_partition <- factor(as.numeric(factor(new_partition)))

    if (!("partitions" %in% names(mon_obj@clusters@listData$UMAP))) {
        mon_obj <- monocle3::cluster_cells(mon_obj)
    }

    mon_obj@clusters@listData$UMAP$partitions <- setNames(new_partition, colnames(mon_obj))
    mon_obj@clusters@listData$UMAP$clusters <- setNames(new_partition, colnames(mon_obj))

    return(mon_obj)
}



