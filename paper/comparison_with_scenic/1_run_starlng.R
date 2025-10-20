library(ClustAssess)
library(Starlng)
library(dplyr)

if (basename(getwd()) != "comparison_with_scenic") {
    setwd("paper/comparison_with_scenic")
}
ca_folder <- "clustassess_app"
starlng_write_app_clustassess_app(
    folder_path = "starlng_app_low_nclust_low",
    app_title_name = "Immune cells - MASLD/MASH",
    ca_app_folder = ca_folder,
    stable_feature_type = "HV",
    stable_feature_set_size = 2000,
    stable_clustering_method = "SLM",
    stable_n_clusters = c(24),
    use_all_genes = TRUE,
    gene_filtering_function = function(df) {
        df %>% dplyr::filter(morans_I > 0.01, q_value < 0.05) %>% rownames
    },
    verbose = TRUE,
    nthreads = 10,
    learn_graph_parameters = list(
        nodes_per_log10_cells = 30,
        learn_graph_controls = list(
            eps = 1e-5,
            maxiter = 10,
            prune_graph = TRUE,
            minimal_branch_len = 5
        )
    ),
    clustering_parameters = list(
        "n_neighbours" = seq(from = 5, to = 50, by = 5),
        "graph_type" = "snn",
        "prune_value" = -1,
        "resolutions" = list(
            "RBConfigurationVertexPartition" = seq(from = 0.1, to = 5, by = 0.1)
        ),
        "number_iterations" = 5,
        "number_repetitions" = 100

    )
)
