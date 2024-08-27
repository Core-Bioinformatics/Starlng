test_that("nn2 parallel produces same output - perfect chunks", {
    set.seed(1234)
    emb <- matrix(runif(1e3 * 30), nrow = 1e3, ncol = 30)
    rownames(emb) <- paste0("cell", seq_len(1e3))
    k <- 25

    default_nn <- parallel_nn2_idx(emb, k)

    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    parallel_nn <- parallel_nn2_idx(emb, k)
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()

    expect_identical(default_nn, parallel_nn)
})

test_that("nn2 parallel produces same output - imbalanced chunks", {
    set.seed(1234)
    emb <- matrix(runif(1001 * 30), nrow = 1001, ncol = 30)
    rownames(emb) <- paste0("cell", seq_len(1001))
    k <- 25

    default_nn <- parallel_nn2_idx(emb, k)

    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    parallel_nn <- parallel_nn2_idx(emb, k)
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()

    expect_identical(default_nn, parallel_nn)
})

test_that("clustering pipeline is robust to parallelization - perfect chunks", {
    set.seed(1234)
    emb <- matrix(runif(1e3 * 30), nrow = 1e3, ncol = 30)
    rownames(emb) <- paste0("cell", seq_len(1e3))
    
    one_core_cl_output <- clustering_pipeline(emb, n_neighbours = c(5, 50), resolutions = list("RBConfigurationVertexPartition" = c(0.1, 1)))

    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    parallel_cl_output <- clustering_pipeline(emb, n_neighbours = c(5, 50), resolutions = list("RBConfigurationVertexPartition" = c(0.1, 1)))
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()

    expect_identical(one_core_cl_output, parallel_cl_output)
})

test_that("clustering pipeline is robust to parallelization - imbalanced chunks", {
    set.seed(1234)
    emb <- matrix(runif(1001 * 30), nrow = 1001, ncol = 30)
    rownames(emb) <- paste0("cell", seq_len(1001))
    
    one_core_cl_output <- clustering_pipeline(emb, n_neighbours = c(5, 50), resolutions = list("RBConfigurationVertexPartition" = c(0.1, 1)))

    cl <- parallel::makeCluster(2)
    doParallel::registerDoParallel(cl)
    parallel_cl_output <- clustering_pipeline(emb, n_neighbours = c(5, 50), resolutions = list("RBConfigurationVertexPartition" = c(0.1, 1)))
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()

    expect_identical(one_core_cl_output, parallel_cl_output)
})

# TODO add tests for `group_by_clusters_general`
# TODO add tests for `assess_stability_by_clusters`