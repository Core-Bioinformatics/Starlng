library(dplyr)
library(gprofiler2)
library(Starlng)

setwd("paper/comparison_with_scenic")

stab_modules <- read.csv(file.path("starlng_app_low_nclust_high", "objects", "stable_modules.csv"))
stab_modules <- stab_modules[, c(1, ncol(stab_modules))]
associated_genes <- split(stab_modules$gene, stab_modules[, 2])
for (i in seq_along(associated_genes)) {
    associated_genes[[i]] <- list(genes = associated_genes[[i]])
}

expr_matrix_path <- file.path("starlng_app_low_nclust_high", "objects", "expression.h5")
expr_matrix <- rhdf5::h5read(expr_matrix_path, "expression_matrix")
rownames(expr_matrix) <- rhdf5::h5read(expr_matrix_path, "genes")
colnames(expr_matrix) <- rhdf5::h5read(expr_matrix_path, "cells")

for (i in seq_along(associated_genes)) {
    temp_enrich <- gprofiler2::gost(
        query = associated_genes[[i]]$genes,
        organism = "hsapiens",
        sources = c("GO", "TF", "KEGG", "REAC"),
        significant = TRUE,
        domain_scope = "custom",
        custom_bg = as.vector(rownames(expr_matrix)),
        evcodes = TRUE
    )
    if (!is.null(temp_enrich)) {
        temp_enrich <- temp_enrich$result
    }


    associated_genes[[i]]$enrichment <- temp_enrich
    associated_genes[[i]]$avg_expr <- voting_scheme(
        expression_matrix = expr_matrix[associated_genes[[i]]$genes, , drop = FALSE],
        genes = associated_genes[[i]]$genes,
        thresh_percentile = 0,
        thresh_value = 0,
        n_coexpressed_thresh = 1,
        summary_function = mean
    )
}
str(associated_genes)

qs2::qs_save(associated_genes, "starlng_high_results.qs2", nthreads = 5)
