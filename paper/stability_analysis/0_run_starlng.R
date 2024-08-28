devtools::load_all()
library(future)
library(doFuture)
library(foreach)
library(dplyr)

setwd("paper/stability_analysis")

dts_name <- "immune"
app_path <- paste0("/mnt/h/Bioinf_data/Pipelines/starlng_app_", dts_name)
objects_path <- file.path(app_path, "objects")

starlng_path <- paste0(dts_name, "_starlng.qs")

if (file.exists(starlng_path) && file.size(starlng_path) > 0) {
    starlng_object <- qs::qread(starlng_path, nthreads = 30)
} else {
    expr_matrix <- rhdf5::h5read(file.path(objects_path, "expression.h5"), "expression_matrix")
    rownames(expr_matrix) <- rhdf5::h5read(file.path(objects_path, "expression.h5"), "genes")
    colnames(expr_matrix) <- rhdf5::h5read(file.path(objects_path, "expression.h5"), "cells")

    mon_obj <- ClustAssess::create_monocle_default(expr_matrix)
    mon_obj <- custom_learn_graph(mon_obj)

    autocorr_result <- monocle3::graph_test(
        mon_obj,
        neighbor_graph = "principal_graph",
        cores = 30,
        verbose = TRUE
    )
    rownames(autocorr_result) <- autocorr_result$gene_short_name
    autocorr_result <- autocorr_result[, c("morans_test_statistic", "morans_I", "q_value")] %>%
        arrange(desc(morans_I)) %>% filter(q_value < 0.05, morans_I > 0.1)
    chosen_genes <- rownames(autocorr_result)

    RhpcBLASctl::blas_set_num_threads(30)
    clustering_parameters <- list(
        "merge_identical_partitions" = TRUE,
        "embedding" = get_feature_loading(expr_matrix[chosen_genes, ], 30)
    )
    RhpcBLASctl::blas_set_num_threads(1)
    qs::qsave(clustering_parameters, file = starlng_path, nthreads = 30)
}

set.seed(202408)
seeds <- sample(1:1000, 30)

par_cluster <- parallel::makePSOCKcluster(30)
doParallel::registerDoParallel(par_cluster)

cluster_outputs <- foreach(seed = seeds) %do% {
    print(seed)
    set.seed(seed)
    clustering_parameters[["seeds"]] <- sample(1:1e5, 100)
    clust_results <- do.call(clustering_pipeline, clustering_parameters)
    best_config <- select_best_configuration(clust_results)
    clust_results <- clust_results[[best_config[1]]][[best_config[2]]]
    clust_results
}

qs::qsave(cluster_outputs, file = paste0(dts_name, "_starlng_cluster_outputs.qs"), nthreads = 30)
