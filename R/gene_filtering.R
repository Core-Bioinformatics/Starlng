
sort_genes_by_metadata <- function(expression_matrix,
                                   metadata_info,
                                   summary_function = mean,
                                   thresh_percentile = 0.25,
                                   thresh_value = 0,
                                   decreasing = FALSE) {
    gene_summary_values <- sapply(rownames(expression_matrix), function(gene) {
        gene_values <- expression_matrix[gene, ]

        if (thresh_percentile > 0) {
            thresh_value <- quantile(gene_values, probs = thresh_percentile)
        }
        index_cells <- which(gene_values > thresh_value)

        return(summary_function(metadata_info[index_cells]))
    })

    return(order(gene_summary_values, decreasing = decreasing))
}