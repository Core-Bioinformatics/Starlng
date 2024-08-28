devtools::load_all()
library(monocle3)
library(future)
library(doFuture)
library(foreach)

setwd("paper/stability_analysis")

dts_name <- "immune"
app_path <- paste0("/mnt/h/Bioinf_data/Pipelines/starlng_app_", dts_name)
objects_path <- file.path(app_path, "objects")

monocle_path <- paste0(dts_name, "_monocle.qs")
# follows the vignette: https://cole-trapnell-lab.github.io/monocle3/docs/differential/

if (file.exists(monocle_path) && file.size(monocle_path) > 0) {
    sifi_object <- qs::qread(monocle_path, nthreads = 30)
} else {
    expr_matrix <- rhdf5::h5read(file.path(objects_path, "expression.h5"), "expression_matrix")
    rownames(expr_matrix) <- rhdf5::h5read(file.path(objects_path, "expression.h5"), "genes")
    colnames(expr_matrix) <- rhdf5::h5read(file.path(objects_path, "expression.h5"), "cells")

    monocle_object <- ClustAssess::create_monocle_default(expr_matrix)
    monocle_object <- preprocess_cds(monocle_object, num_dim = 100)
    monocle_object <- reduce_dimension(monocle_object)
    monocle_object <- cluster_cells(monocle_object, resolution = 1e-5)

    pr_graph_test_res <- graph_test(monocle_object, neighbor_graph = "knn", cores = 8)
    pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
    monocle_object <- monocle_object[pr_deg_ids, ]

    qs::qsave(monocle_object, file = monocle_path, nthreads = 30)
}

set.seed(202408)
seeds <- sample(1:1000, 30)

cluster_outputs <- foreach(seed = seeds, .options.future = list(seed = TRUE)) %dofuture% {
    set.seed(seed)
    find_gene_modules(monocle_object, resolution = 1e-3, random_seed = seed)
}

cluster_outputs <- lapply(cluster_outputs, function(x) {
    as.numeric(x$module)
})

qs::qsave(cluster_outputs, file = paste0(dts_name, "_monocle_cluster_outputs.qs"), nthreads = 30)




