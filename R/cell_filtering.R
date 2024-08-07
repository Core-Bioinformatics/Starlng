#' @importFrom dplyr %>% .data
voting_scheme <- function(expression_matrix,
                          genes,
                          thresh_percentile = 0.25,
                          thresh_value = 0,
                          n_coexpressed_thresh = length(genes),
                          summary_function = NULL) {
    if (is.null(rownames(expression_matrix))) {
        stop("Expression matrix has no row names")
    }

    genes <- intersect(genes, rownames(expression_matrix))

    if (length(genes) == 0) {
        stop("No genes found in the expression matrix")
    }

    if (is.null(colnames(expression_matrix))) {
        colnames(expression_matrix) <- seq_len(ncol(expression_matrix))
    }

    expression_matrix <- expression_matrix[genes, , drop = FALSE]
    n_coexpressed_thresh <- max(min(n_coexpressed_thresh, nrow(expression_matrix)), 1)

    if (thresh_percentile > 0) {
        thresh_value <- sapply(genes, function(gene) {
            non_zero_values <- expression_matrix[gene, ]
            non_zero_values <- non_zero_values[non_zero_values > 0]
            quantile(non_zero_values, probs = thresh_percentile)
        })
    } else {
        thresh_value <- rep(thresh_value[1], length(genes))
    }

    mask_expressed <- expression_matrix > thresh_value
    index_cells <- colSums(mask_expressed) >= n_coexpressed_thresh
    if (is.null(summary_function)) {
        cell_info <- rep("not selected", ncol(expression_matrix))
        cell_info[index_cells] <- "selected"
        names(cell_info) <- colnames(expression_matrix)

        return(factor(cell_info))
    }

    index_cells <- which(index_cells)
    cell_info <- rep(0, ncol(expression_matrix))
    names(cell_info) <- colnames(expression_matrix)
    if (length(index_cells) == 0) {
        return(cell_info)
    }

    cell_info[index_cells] <- sapply(index_cells, function(cell_index) {
        cell_values <- expression_matrix[, cell_index]
        return(summary_function(cell_values))
    })

    return(cell_info)
}

select_cells_by_gene_expr <- function(expression_matrix,
                                      genes,
                                      thresh_percentile = 0.25,
                                      thresh_value = 0,
                                      n_coexpressed_thresh = length(genes)) {
    # genes <- intersect(genes, rownames(expression_matrix))

    # if (length(genes) == 0) {
    #     stop("No genes found in the expression matrix")
    # }

    # expression_matrix <- expression_matrix[genes, ]
    # n_coexpressed <- max(min(n_coexpressed, nrow(expression_matrix)), 1)

    # if (thresh_percentile > 0) {
    #     thresh_value <- sapply(genes, function(gene) {
    #         quantile(sum(expression_matrix[gene, ]), probs = thresh_percentile)
    #     })
    # } else {
    #     thresh_value <- rep(thresh_value[1], length(genes))
    # }

    # mask_expressed <- expression_matrix > thresh_value

    # index_cells <- which(colSums(mask_expressed) >= n_coexpressed_thresh)

    # if (is.null(colnames(expression_matrix))) {
    #     return(index_cells)
    # }
    voting_scheme_results <- voting_scheme(expression_matrix, genes, thresh_percentile, thresh_value, n_coexpressed_thresh)
    index_cells <- which(voting_scheme_results == "selected")

    return(names(voting_scheme_results)[index_cells])
}

# metadata combinations is a list metadata column - metadata unique values
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

remove_outlier_cells <- function(cell_names, umap_emb, percentile_threshold = 0.75, gmedian_point = NULL) {
    if (is.null(rownames(umap_emb))) {
        rownames(umap_emb) <- seq_len(nrow(umap_emb))
    }

    if (length(intersect(cell_names, rownames(umap_emb))) == 0) {
        stop("The cells were not found in the UMAP embedding")
    }

    umap_emb <- umap_emb[cell_names, ]

    if (is.null(gmedian_point)) {
        gmedian_point <- Gmedian::Gmedian(umap_emb)
    }
    distance_cells_gmedian <- (
        (umap_emb[, 1] - gmedian_point[1, 1]) ^ 2 +
        (umap_emb[, 2] - gmedian_point[1, 2]) ^ 2
    ) ^ 0.5

    threshold_distance <- quantile(distance_cells_gmedian, probs = percentile_threshold)

    return(cell_names[distance_cells_gmedian <= threshold_distance])
}

filter_central_cells_from_group <- function(cell_group, umap_embedding, n_points = 5) {
    if (!(all(cell_group %in% rownames(umap_embedding)))) {
        stop("Some cells were not found in the UMAP embedding")
    }

    umap_embedding <- umap_embedding[cell_group, , drop = FALSE]
    if (nrow(umap_embedding) == 0) {
        stop("No cells found in the UMAP embedding based on the provided group.")
    }

    if (nrow(umap_embedding) <= n_points) {
        return(cell_group)
    }

    gmedian_point <- Gmedian::Gmedian(umap_embedding)
    distance_cells_gmedian <- (
        (umap_embedding[, 1] - gmedian_point[1, 1]) ^ 2 +
        (umap_embedding[, 2] - gmedian_point[1, 2]) ^ 2
    ) ^ 0.5

    start_threshold <- 0
    end_threshold <- 1
    prev_distance <- 0
    while (start_threshold < end_threshold) {
        threshold_mid <- round((start_threshold + end_threshold) / 2, 2)
        if (prev_distance == threshold_mid) {
            break
        }

        threshold_distance <- quantile(distance_cells_gmedian, probs = threshold_mid)
        index_cells <- which(distance_cells_gmedian <= threshold_distance)

        if (length(index_cells) <= n_points) {
            start_threshold <- threshold_mid
        } else {
            end_threshold <- threshold_mid
        }
        prev_distance <- threshold_mid
    }

    return(cell_group[index_cells])
}
