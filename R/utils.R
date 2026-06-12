
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

scale_min_max <- function(x) {
    min_val <- min(x)
    max_val <- max(x)

    if (min_val == max_val) {
        return(x)
    }

    return((x - min_val) / (max_val - min_val))
}

verbose_print <- function(message_str, verbose = TRUE) {
    if (verbose) {
        print(paste0("[", Sys.time(), "] ", message_str))
    }
    utils::flush.console()
}

mad_z_score <- function(x) {
    if (length(x) < 2) {
        return(rep(FALSE, length(x)))
    }

    med <- stats::median(x, na.rm = TRUE)
    mad_val <- stats::mad(x, center = med, constant = 1.4826, na.rm = TRUE)

    if (mad_val == 0) {
        return(rep(0, length(x)))
    }

    z_scores <- (x - med) / mad_val
    return(z_scores)
}

#' Diet Monocle object
#' 
#' @description This function simplifies the structure of the Monocle object.
#' The purpose is to create a lightweight object that can be used in the
#' context of the Shiny application. The resulting object will only contain
#' the needed information to perform the pseudotime ordering. Therefore,
#' this function removes information such as: the expression matrices, the
#' auxiliary data related to the dimensionality reduction and the clustering
#' and the metadata columns.
#' 
#' @param mon_obj Monocle object.
#' 
#' @return Diet Monocle object.
#' @export
diet_monocle_object <- function(mon_obj) {
    mask <- rep(TRUE, nrow(mon_obj) * ncol(mon_obj))
    # remove all count information
    mon_obj@assays@data@listData$counts[mask] <- 0
    # delete the other assays
    for (slot_names in names(mon_obj@assays@data@listData)) {
        if (slot_names == "counts") {
            next
        }
        mon_obj@assays@data@listData[[slot_names]] <- NULL
    }

    # keep only one metadata column
    mon_obj@colData[, 2:ncol(mon_obj@colData)] <- NULL

    # keep only the UMAP dim reduction
    for (dim_reduc_name in names(mon_obj@reduce_dim_aux@listData)) {
        if (dim_reduc_name == "UMAP") {
            next
        }
        mon_obj@reduce_dim_aux@listData[[dim_reduc_name]] <- NULL
    }

    for (dim_reduc_name in names(mon_obj@int_colData@listData$reducedDims@listData)) {
        if (dim_reduc_name == "UMAP") {
            next
        }
        mon_obj@int_colData@listData$reducedDims@listData[[dim_reduc_name]] <- NULL
    }

    for (dim_reduc_name in names(mon_obj@clusters)) {
        if (dim_reduc_name == "UMAP") {
            next
        }

        mon_obj@clusters[[dim_reduc_name]] <- NULL
    }

    mon_obj@clusters$UMAP$cluster_result <- NULL

    return(mon_obj)
}
