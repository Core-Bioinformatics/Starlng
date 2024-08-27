
# assumes that the expr matrix is gene x cell
# assumes that the genes have been already filtered i.e only the genes that are meant to be used are present
#' Get feature loading from an expression matrix
#' 
#' @description This function calculates the feature loading of an expression
#' matrix using PCA.
#' 
#' @param expr_matrix The gene by cell expression matrix to be used. It is
#' expected that the matrix is already normalized and scaled.
#' @param npcs The number of principal components to be used. Defaults to 30.
#' @param approx Logical indicating if deterministic PCA (with `prcomp`) or
#' approximate PCA (with `irlba`) should be used. Defaults to FALSE.
#' @param ... Aditional parameters passed to the PCA function.
#' 
#' @return A gene by npcs matrix containing the feature loading.
#' @export
get_feature_loading <- function(expr_matrix, npcs = 30, approx = FALSE, ...) {
    transpose_f <- base::t
    if (inherits(expr_matrix, "dgCMatrix")) {
        transpose_f <- Matrix::t
    }
    
    if (!approx) {
        return(stats::prcomp(
            transpose_f(expr_matrix),
            center = TRUE,
            scale = TRUE,
            rank. = npcs,
            retx = TRUE,
            ...
        )$rotation)
    }

    return(irlba::irlba(
        transpose_f(expr_matrix),
        nv = npcs,
        ...
    )$v)
}

#' Get PCA reduction from an expression matrix
#' 
#' @description This function calculates the PCA embedding of an expression
#' matrix.
#' 
#' @note The PCA is applied in this context to the genes, not to the cells.
#' 
#' @param expr_matrix The gene by cell expression matrix to be used. It is
#' expected that the matrix is already normalized and scaled.
#' @param npcs The number of principal components to be used. Defaults to 30.
#' @param approx Logical indicating if deterministic PCA (with `prcomp`) or
#' approximate PCA (with `irlba`) should be used. Defaults to FALSE.
#' @param ... Aditional parameters passed to the PCA function.
#' 
#' @return A gene by npcs matrix containing the PCA embedding.
#' @export
pca_reduction <- function(expr_matrix, npcs = 30, approx = FALSE, ...) {
    if (!approx) {
        return(stats::prcomp(
            expr_matrix,
            center = TRUE,
            scale = TRUE,
            rank. = npcs,
            ...
        )$x)
    }

    pca_res <- irlba::irlba(
        expr_matrix,
        nv = npcs,
        ...
    )$u

    return(pca_res$u %*% diag(pca_res$d))
}

pseudobulk_reduction <- function(expr_matrix, cell_clusters) {
    if (is.null(colnames(expr_matrix))) {
        colnames(expr_matrix) <- seq_len(ncol(expr_matrix))
    }

    barcode_grouping <- split(colnames(expr_matrix), cell_clusters)

    pseudobulk_matrix <- lapply(
        barcode_grouping,
        function(barcode_group) {
            rowMeans(expr_matrix[, barcode_group])
        }
    )

    do.call(cbind, pseudobulk_matrix)
}
