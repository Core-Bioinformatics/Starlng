
read_gene_from_dense_h5 <- function(gene_names,
                                    matrix_h5_path,
                                    index_genes = NULL,
                                    check_intersect = TRUE,
                                    add_rownames = FALSE) {
    if (check_intersect) {
        gene_names <- intersect(gene_names, rhdf5::h5read(matrix_h5_path, "genes"))
    }

    if (length(gene_names) == 0) {
        stop("No genes found in the expression matrix")
    }

    if (is.null(index_genes)) {
        index_genes <- which(rhdf5::h5read(matrix_h5_path, "genes") %in% gene_names)
    }

    gene_matrix <- rhdf5::h5read(
        matrix_h5_path,
        "expression_matrix",
        index = list(index_genes, NULL)
    )
    rownames(gene_matrix) <- gene_names

    if (add_rownames) {
        colnames(gene_matrix) <- rhdf5::h5read(matrix_h5_path, "cells")
    }

    return(gene_matrix)
}

read_gene_from_sparse_h5 <- function() {

}
