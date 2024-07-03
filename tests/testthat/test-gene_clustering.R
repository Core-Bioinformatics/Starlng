test_that("nn2 parallel produces same output", {
    set.seed(2024)
    emb <- matrix(runif(1e3 * 30), nrow = 1e3, ncol = 30)
    k <- 30

    default_nn <- parallel_nn2_idx(emb, k)

    doParallel::registerDoParallel(cores = 2)
    parallel_nn <- parallel_nn2_idx(emb, k)
    foreach::registerDoSEQ()

    expect_identical(default_nn, parallel_nn)
})

