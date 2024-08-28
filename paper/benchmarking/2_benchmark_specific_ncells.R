
devtools::load_all("/mnt/h/Bioinf_data/Pipelines/Starlng")
library(microbenchmark)
library(dplyr)
# library(Starlng)

dts_name <- commandArgs(trailingOnly = TRUE)[1]
dir.create(dts_name, showWarnings = FALSE)
ncells <- commandArgs(trailingOnly = TRUE)[2]
ncores <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
nthreads_learn <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
print(paste(dts_name, ncells, ncores, nthreads_learn))

data_path <- "/mnt/h/Bioinf_data/Pipelines/ClustAssess-paper/benchmarking/immune"
expression_matrix <- qs::qread(file.path(data_path, paste0("zinbSim_normalized_", ncells, ".qs")), nthreads = 30)

mcb_res <- microbenchmark::microbenchmark({
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    starlng_write_app_default(
        expression_matrix = expression_matrix,
        metadata_df = NULL,
        pca_embedding = NULL,
        umap_embedding = NULL,
        app_title = paste0(dts_name, " - ", ncells),
	folder_path = file.path(dts_name, ncells),
	gene_filtering_function = function(info_gene_df) {
		rownames(info_gene_df %>% dplyr::filter(.data$morans_I > 0, .data$q_value < 0.05))
	},
        nthreads = nthreads_learn
    )
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    gc()
}, times = 10)

qs::qsave(mcb_res, file = file.path(dts_name, paste0("mcb_res_", ncells, "_ncores_", ncores, "_for_graph_learn_", nthreads_learn, ".qs")), nthreads = 30)
