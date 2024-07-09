
# TODO you will need to check if the process is actually fork or psock
# if it's fork maybe it's not that safe to remove all the content
clear_psock_memory <- function() {
    ncores <- foreach::getDoParWorkers()

    if (ncores == 1) {
        return()
    }

    foreach::foreach (i = seq_len(ncores)) %dopar% {
        rm(list = ls())
        gc()
    }
}