#' @importFrom dplyr %>% .data
NULL

#' Gene Expression voting scheme
#' 
#' @description This function performs a voting scheme to identify cells that
#' are expressing one or multiple genes above a given threshold. The function
#' allows the user to introduce a factor of relaxation in the voting scheme
#' by setting the minimum number of genes that should be coexpressed in a cell
#' for it to be considered. For the obtained set of cells, the function
#' summarizes the expression values of the genes using a user-defined function.
#' 
#' @param expression_matrix The gene by cell expression matrix to be used. The
#' matrix should have defined rownames and colnames.
#' @param genes A vector with the genes to be used for the voting scheme.
#' @param thresh_percentile The percentile to be used as threshold for the
#' expression values to select the cells associated with each gene. If set to 0,
#' the threshold will be the value defined in `thresh_value`. Defaults to 0.25.
#' @param thresh_value The value to be used as threshold for the expression
#' values to select the cells associated with each gene. Defaults to 0.
#' @param n_coexpressed_thresh The number of genes that should be coexpressed
#' in a cell for it to be considered as selected. Defaults to the number of
#' provided genes.
#' @param summary_function A function that takes as argument a numeric vector
#' and summarises it into a single value. If set to NULL, the function will
#' return a factor with the selected cells. Defaults to NULL.
#' 
#' @return A list with the information of the cells selected by the voting
#' scheme. If a summary function is provided, the list will contain the summary
#' values for each cell.
#' @export
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
            stats::quantile(non_zero_values, probs = thresh_percentile)
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

#' Select cells by gene expression
#' 
#' @description This function filters from the dataset a group of cells which
#' are expressing a set of genes above a given threshold. For this, the voting
#' scheme will be applied.
#' 
#' @param expression_matrix The gene by cell expression matrix to be used.
#' @param genes A vector with the names of the genes to be used for the
#' filtering. 
#' @param thresh_percentile The percentile to be used as threshold for the
#' expression values to select the cells associated with each gene. If set to 0,
#' the threshold will be the value defined in `thresh_value`. Defaults to 0.25.
#' @param thresh_value The value to be used as threshold for the expression
#' values to select the cells associated with each gene. Defaults to 0.
#' @param n_coexpressed_thresh The number of genes that should be coexpressed
#' in a cell for it to be considered as selected. Defaults to the number of
#' provided genes.
#' 
#' @return A list of cells that are selected using the voting scheme.
#' @export
select_cells_by_gene_expr <- function(expression_matrix,
                                      genes,
                                      thresh_percentile = 0.25,
                                      thresh_value = 0,
                                      n_coexpressed_thresh = length(genes)) {
    voting_scheme_results <- voting_scheme(expression_matrix, genes, thresh_percentile, thresh_value, n_coexpressed_thresh)
    index_cells <- which(voting_scheme_results == "selected")

    return(names(voting_scheme_results)[index_cells])
}

#' Select cells by metadata
#' 
#' @description This function filters from the dataset a group of cells which
#' is defined as a combination of filters determined by multiple metadata.
#' 
#' @param metadata_df A dataframe with the metadata information. The rows of
#' the dataframe should be the cells and the columns the metadata information.
#' @param metadata_combinations A list with the metadata columns and the unique
#' groups of the metadata that should be used to filter the cells. The names of
#' the list should be in the columns of the `metadata_df`.
#' 
#' @return The list of cells that are defined by the intersection of provided
#' groups.
#' @export
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

#' Remove outlier cells from a group
#' 
#' @description This function removes the outlier cells from a group based on
#' the distance patterns, as observed on an UMAP embedding. The function
#' considers a point to be outlier if the distance to the geometric median
#' of the group is above a given threshold determined by a quantile.
#' 
#' @param cell_names A vector with the names of the cells defining the group.
#' @param umap_emb A matrix with the UMAP embedding of the cells. The matrix
#' should have defined the rownames as the cell names.
#' @param percentile_threshold The percentile to be used as threshold for the
#' distance to the geometric median. This value will be used to calculate
#' the quantile. Defaults to 0.75.
#' @param gmedian_point The geometric median point of the group. If NULL, the
#' function will calculate this point. Defaults to NULL.
#' 
#' @return A vector with the names of the cells that are not considered
#' outliers.
#' @export
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

    threshold_distance <- stats::quantile(distance_cells_gmedian, probs = percentile_threshold)

    return(cell_names[distance_cells_gmedian <= threshold_distance])
}

#' Find central points of a group
#' 
#' @description Given a population of cells and a UMAP embedding, this function
#' attempts to identify a group of cells that are located in the centre of this
#' group. For this, the function calculates the geometric median of the UMAP
#' and then calculates the nearest points to this point.
#' 
#' @param cell_group A vector with the names of the cells to be analysed.
#' @param umap_embedding A matrix with the UMAP embedding of the cells. The
#' matrix should have defined the rownames as the cell names.
#' @param n_points The number of points to be selected as the central points.
#' 
#' @return A vector with the names of the cells that are considered central.
#' @export
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

        threshold_distance <- stats::quantile(distance_cells_gmedian, probs = threshold_mid)
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
