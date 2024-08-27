
#' Read expression matrix from a dense HDF5 file
#' 
#' @description This function read specific rows, associated with some genes,
#' from a dense matrix stored in a HDF5 file.
#' 
#' @param gene_names List of names of genes to be read.
#' @param matrix_h5_path Path to the HDF5 file where the expression matrix
#' is stored. The file must contain the following fields:
#' - genes: the names of the genes in the matrix.
#' - cells: the names of the cells in the matrix.
#' - expression_matrix: the expression matrix
#' @param index_genes Named vector where the names are the genes and the
#' values the index of the row associated with a gene. If NULL, the function
#' will calculate automatically this index. Defaults to NULL.
#' @param check_intersect Logical indicating if the function should check if
#' all the provided genes are present in the matrix. Defaults to TRUE.
#' @param add_rownames Logical indicating if the function should add the cell
#' names as column names for the extracted matrix. Defaults to FALSE.
#' 
#' @return The expression matrix having the requested genes as rows.
#' @export
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
        genes <- rhdf5::h5read(matrix_h5_path, "genes")
        index_genes <- which(genes %in% gene_names)
        gene_names <- genes[index_genes]
    }

    gene_matrix <- rhdf5::h5read(
        matrix_h5_path,
        "expression_matrix",
        index = list(index_genes, NULL)
    )
    rownames(gene_matrix) <- gene_names

    if (add_rownames) {
        colnames(gene_matrix) <- as.character(rhdf5::h5read(matrix_h5_path, "cells"))
    }

    return(gene_matrix)
}

read_gene_from_sparse_h5 <- function() {

}

#' Write expression matrix to a dense HDF5 file
#' 
#' @description This function writes a expression matrix to a HDF5 file in a
#' dense format. The created file will contain three fields:
#' - genes: the names of the genes in the matrix.
#' - cells: the names of the cells in the matrix.
#' - expression_matrix: the expression matrix
#' 
#' @param expression_matrix The expression matrix to be written. The matrix
#' should have defined rownames and colnames. If they are not defined, they
#' will be generated as `gene_{index}` and `cell_{index}` respectively.
#' @param file_path The path to the HDF5 file where the matrix will be stored.
#' @param compression_level The compression level to be used. Defaults to 7.
#' @param chunk_size The size of the chunks to be used. Defaults to 100.
#' @param all_at_once Logical indicating if the writing process will be done
#' by providing the entire matrix or by writing chunks of the matrix. It is
#' recommended to write the entire matrix at once, as no noticeable speedup
#' is observed in the other case. Defaults to TRUE.
#' @param garbage_thresh The threshold, measured GB, to trigger a garbage
#' collection while writing the matrix by chunks.
#' 
#' @return The function does not return anything. It will create the HDF5 file
#' with the matrix stored in it.
#' @export
write_gene_matrix_dense_h5 <- function(expression_matrix,
                                       file_path,
                                       compression_level = 7,
                                       chunk_size = 100,
                                       all_at_once = TRUE,
                                       garbage_thresh = 2) {

    if (file.exists(file_path)) {
        warning("Overwriting existing file. Waiting for 5 seconds...")
        Sys.sleep(5.5)
        file.remove(file_path)
    }

    if (is.null(rownames(expression_matrix))) {
        rownames(expression_matrix) <- paste0("gene_", seq_len(nrow(expression_matrix)))
    }

    if (is.null(colnames(expression_matrix))) {
        colnames(expression_matrix) <- paste0("cell_", seq_len(ncol(expression_matrix)))
    }
    rhdf5::h5createFile(file_path)
    rhdf5::h5createDataset(
        file_path,
        "expression_matrix",
        dims = dim(expression_matrix),
        maxdims = dim(expression_matrix),
        storage.mode = "double",
        chunk = c(chunk_size, ncol(expression_matrix)),
        level = compression_level
    )

    rhdf5::h5write(rownames(expression_matrix), file_path, "genes")
    rhdf5::h5write(colnames(expression_matrix), file_path, "cells")

    if (all_at_once) {
        if (inherits(expression_matrix, "dgCMatrix")) {
            expression_matrix <- as.matrix(expression_matrix)
        }
        rhdf5::h5write(expression_matrix, file_path, "expression_matrix")
        return()
    }

    nchunks <- ceiling(nrow(expression_matrix) / chunk_size)
    ngenes <- nrow(expression_matrix)

    garbage_sum <- 0
    for (i in seq_len(nchunks)) {
        start <- (i - 1) * chunk_size + 1
        end <- min(i * chunk_size, ngenes)
        written_mat <- expression_matrix[start:end, ]
        if (garbage_thresh > 0) garbage_sum <- garbage_sum + utils::object.size(written_mat) / 1024^3
        if (inherits(written_mat, "dgCMatrix")) {
            written_mat <- as.matrix(written_mat)
            if (garbage_thresh > 0) garbage_sum <- garbage_sum + utils::object.size(written_mat) / 1024^3
        }
        if (garbage_thresh > 0 && garbage_sum > garbage_thresh) {
            gc()
            garbage_sum <- 0
        }

        rhdf5::h5write(
            obj = written_mat,
            file = file_path,
            name = "expression_matrix",
            index = list(start:end, NULL)
        )
    }
    rhdf5::h5closeAll()
}