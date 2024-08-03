
# assumes that the expr matrix is gene x cell
# assumes that the genes have been already filtered i.e only the genes that are meant to be used are present
#' @export
get_feature_loading <- function(expr_matrix, npcs = 30) {
    transpose_f <- base::t
    if (inherits(expr_matrix, "dgCMatrix")) {
        transpose_f <- Matrix::t
    }
    prcomp(
        transpose_f(expr_matrix),
        center = TRUE,
        scale = TRUE,
        rank. = npcs,
        retx = TRUE
    )$rotation
}

pca_reduction <- function(expr_matrix, npcs = 30) {
    prcomp(
        expr_matrix,
        center = TRUE,
        scale = TRUE,
        rank. = npcs
    )$x
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
