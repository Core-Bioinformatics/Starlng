
clear_psock_memory <- function() {
    ncores <- foreach::getDoParWorkers()

    if (ncores == 1) {
        return()
    }

    foreach (i = seq_len(ncores)) %dopar% {
        rm(list = ls())
    }
}