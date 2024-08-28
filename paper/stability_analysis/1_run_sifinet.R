
# library(Starlng)
devtools::load_all()
library(SiFINeT)
library(future)
library(doFuture)
library(foreach)

setwd("paper/stability_analysis")

dts_name <- "immune"
app_path <- paste0("/mnt/h/Bioinf_data/Pipelines/starlng_app_", dts_name)
objects_path <- file.path(app_path, "objects")

sifi_path <- paste0(dts_name, "_sifi.qs")

if (file.exists(sifi_path) && file.size(sifi_path) > 0) {
    sifi_object <- qs::qread(sifi_path, nthreads = 30)
} else {
    expr_matrix <- rhdf5::h5read(file.path(objects_path, "expression.h5"), "expression_matrix")
    rownames(expr_matrix) <- rhdf5::h5read(file.path(objects_path, "expression.h5"), "genes")
    colnames(expr_matrix) <- rhdf5::h5read(file.path(objects_path, "expression.h5"), "cells")

    sifi_object <- create_SiFINeT_object(
        counts = expr_matrix,
        gene.name = NULL,
        meta.data = NULL,
        data.name = NULL,
        sparse = FALSE,
        rowfeature = TRUE
    )

    sifi_object <- quantile_thres(sifi_object)
    sifi_object <- feature_coexp(sifi_object)

    sifi_object <- create_network(
        so = sifi_object,
        alpha = 0.05,
        manual = FALSE,
        least_edge_prop = 0.01
    )

    sifi_object <- filter_lowexp(
        so = sifi_object,
        t1 = 10,
        t2 = 0.9,
        t3 = 0.9
    )

    qs::qsave(sifi_object, file = sifi_path, nthreads = 30)
}

set.seed(202408)
seeds <- sample(1:1000, 30)
print(seeds)

plan(multicore, workers = 5)

cluster_outputs <- foreach(seed = seeds) %dofuture% {
    set.seed(seed)
    sifi_object <- cal_connectivity(
        so = sifi_object,
        m = 10,
        niter = 100
    )

    sifi_object <- find_unique_feature(
        so = sifi_object,
        t1 = 5,
        t2 = 0.4,
        t3 = 0.3,
        t1p = 5,
        t2p = 0.7,
        t3p = 0.5,
        resolution = 1,
        min_set_size = 5
    )

    return(sifi_object@uni_cluster)
}

qs::qsave(cluster_outputs, file = paste0(dts_name, "_sifi_clusters.qs"), nthreads = 30)
plan(sequential)
