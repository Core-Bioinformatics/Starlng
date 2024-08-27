
#' Sort genes by metadata information
#' 
#' @description This function sorts the genes in an expression matrix based on
#' a summary that is performed on a population of cells. For each gene,
#' the population is selected as the cells expressing that gene above a given
#' threshold. The summary is calculated on the metadata information of the
#' obtained population.
#' 
#' @param expression_matrix The gene by cell expression matrix to be used.
#' @param metadata_info The metadata information to be used for the sorting.
#' The metadata should be a numeric vector.
#' @param summary_function A function that takes as argument a numeric vector
#' and summarises it into a single value. Defaults to `mean`.
#' @param thresh_percentile The percentile to be used as threshold for the
#' expression values to select the cells associated with each gene. If set to 0,
#' the threshold will be the value defined in `thresh_value`. Defaults to 0.25.
#' @param thresh_value The value to be used as threshold for the expression
#' values to select the cells associated with each gene. Defaults to 0.
#' @param decreasing Logical indicating if the genes should be sorted in
#' decreasing or increasing order. Defaults to FALSE.
#' 
#' @return A vector with the indices of the genes sorted by the metadata
#' summaries.
#' @export
sort_genes_by_metadata <- function(expression_matrix,
                                   metadata_info,
                                   summary_function = mean,
                                   thresh_percentile = 0.25,
                                   thresh_value = 0,
                                   decreasing = FALSE) {
    gene_summary_values <- sapply(rownames(expression_matrix), function(gene) {
        gene_values <- expression_matrix[gene, , drop = FALSE]

        if (thresh_percentile > 0) {
            thresh_value <- stats::quantile(gene_values, probs = thresh_percentile)
        }
        index_cells <- which(gene_values > thresh_value)

        return(summary_function(metadata_info[index_cells]))
    })

    return(order(gene_summary_values, decreasing = decreasing))
}
