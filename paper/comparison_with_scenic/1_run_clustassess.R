library(Seurat)
library(ClustAssess)

setwd("paper/comparison_with_scenic")
so <- readRDS("immuneCellsSCTransformed.rds")
so <- SCTransform(so, variable.features.n = 3000, return.only.var.genes = FALSE, verbose = FALSE)

features <- dimnames(so@assays$RNA)[[1]]
var_features <- so@assays[["SCT"]]@var.features
n_abundant <- 3000
most_abundant_genes <- rownames(so@assays$SCT)[order(Matrix::rowSums(so@assays$SCT),
    decreasing = TRUE
)]

RhpcBLASctl::blas_set_num_threads(1)
ncores <- 10
my_cluster <- parallel::makeCluster(
    ncores,
    type = "PSOCK"
)

steps <- seq(from = 500, to = 3000, by = 500)
ma_hv_genes_intersection_sets <- sapply(steps, function(x) intersect(most_abundant_genes[1:x], var_features[1:x]))
ma_hv_genes_intersection <- Reduce(union, ma_hv_genes_intersection_sets)
ma_hv_steps <- sapply(ma_hv_genes_intersection_sets, length)

automm_output <- automatic_stability_assessment(
    expression_matrix = so@assays$SCT@scale.data,
    n_repetitions = 100,
    n_neigh_sequence = seq(from = 5, to = 50, by = 5),
    resolution_sequence = seq(from = 0.1, to = 2, by = 0.1),
    features_sets = list(
        "HV" = var_features,
        "MA" = most_abundant_genes[seq_len(3000)]
    ),
    steps = list(
        "HV" = steps,
        "MA" = steps
    ),
    n_top_configs = 2,
    umap_arguments = list(
        min_dist = 0.3,
        n_neighbors = 30,
        metric = "cosine"
    ),
    save_temp = FALSE,
    verbose = TRUE
)

saveRDS(automm_output, "clustassess_output.rds")
parallel::stopCluster(my_cluster)
foreach::registerDoSEQ()

write_shiny_app(
    object = so,
    assay_name = "SCT",
    clustassess_object = automm_output,
    project_folder = "clustassess_app",
    shiny_app_title = "immune dataset",
    qualpalr_colorspace = list(h = c(0, 360), s = c(0.2, 0.5), l = c(0.6, 0.85)),
    prompt_feature_choice = FALSE
)
