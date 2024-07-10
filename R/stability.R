
group_by_clusters_general <- function(clustering_list, start_level = 1) {
    if (start_level > 1) {
        if (is.null(names(clustering_list))) {
            stop(paste0("The list doesn't have as much levels as specified (", start_level, ")"))
        }

        for (sublist in names(clustering_list)) {
            clustering_list[[sublist]] <- group_by_clusters_general(clustering_list[[sublist]], start_level - 1)
        }

        return(clustering_list)
    }

    by_cluster_list <- list()
    if (is.null(names(clustering_list))) {
        for (mb in clustering_list) {
            k_value <- length(unique(mb))
            if (!k_value %in% names(by_cluster_list)) {
                by_cluster_list[[as.character(k_value)]] <- list(mb)
            } else {
                by_cluster_list[[as.character(k_value)]] <- c(by_cluster_list[[as.character(k_value)]], list(mb))
            }
        }
        
        return(list(k = by_cluster_list))
    }

    for (sublist in names(clustering_list)) {
        by_cluster_sublist <- group_by_clusters_general(clustering_list[[sublist]], start_level)$k

        for (k_value in names(by_cluster_sublist)) {
            if (!k_value %in% names(by_cluster_list)) {
                by_cluster_list[[k_value]] <- by_cluster_sublist[[k_value]]
            } else {
                by_cluster_list[[k_value]] <- c(by_cluster_list[[k_value]], by_cluster_sublist[[k_value]])
            }
        }
    }

    k_values <- names(by_cluster_list)
    k_values <- stringr::str_sort(k_values, numeric = TRUE)
    by_cluster_list <- by_cluster_list[k_values]

    return(list(k = by_cluster_list))
}

get_clusters_consistency <- function(by_cluster_list, order_logic = "agreement") {
    if (!("k" %in% names(by_cluster_list))) {
        for (sub_name in names(by_cluster_list)) {
            by_cluster_list[[sub_name]] <- get_clusters_consistency(by_cluster_list[[sub_name]], order_logic)
        }
        return(by_cluster_list)
    }

    k_values <- names(by_cluster_list$k)
    by_cluster_list$k <- foreach::foreach (
        k_list = by_cluster_list$k,
        .final = function(x) setNames(x, k_values)
    ) %dopar% { 
        ClustAssess::merge_partitions(k_list, order_logic = order_logic)
    }

    k_freqs <- sapply(by_cluster_list$k, function(k_list) {
        sum(sapply(k_list$partitions, function(part) part$freq))
    })

    total_freqs <- sum(k_freqs)
    weighted_ecc_by_k <- lapply(k_values, function(k_value) {
        by_cluster_list$k[[k_value]]$ecc * k_freqs[[k_value]]
    })
    weighted_ecc <- Reduce("+", weighted_ecc_by_k) / total_freqs

    by_cluster_list$overall_ecc <- weighted_ecc
    by_cluster_list$total_freqs <- total_freqs

    return(by_cluster_list)
}

get_all_configurations <- function(grouped_by_k_list,
                                   prefix = c(),
                                   sep_char = ";",
                                   first_summary_function = median,
                                   second_summary_function = function(x) { quantile(x, 0.75) - quantile(x, 0.25) }) {
    k_values <- names(grouped_by_k_list$k)

    sublist_names <- names(grouped_by_k_list)

    if ("overall_ecc" %in% sublist_names) {
        overall_ecc <- grouped_by_k_list$overall_ecc
        summary_results <- c(
            first_summary_function(overall_ecc),
            second_summary_function(overall_ecc)
        )
        names(summary_results) <- NULL
        return_list <- list(prefix = summary_results)
        names(return_list) <- paste(prefix, collapse = sep_char)
        return(return_list)
    }

    return(do.call(c, lapply(sublist_names, function(sublist_name) {
        get_all_configurations(
            grouped_by_k_list[[sublist_name]],
            c(prefix, sublist_name),
            sep_char,
            first_summary_function,
            second_summary_function
        )
    })))
}

select_best_configuration <- function(grouped_by_k_list,
                                      sep_char = ";",
                                      first_ranking_function = median,
                                      second_ranking_function = function(x) { quantile(x, 0.75) - quantile(x, 0.25) }) {
    all_configurations <- get_all_configurations(
        grouped_by_k_list,
        c(),
        sep_char = sep_char,
        first_summary_function = first_ranking_function,
        second_summary_function = second_ranking_function
    )

    first_ranking_values <- sapply(all_configurations, function(config) config[1])
    mean_first_value <- mean(first_ranking_values)
    eligible_indices <- which(first_ranking_values >= mean_first_value)

    all_configurations <- all_configurations[eligible_indices]
    second_ranking_values <- sapply(all_configurations, function(config) config[2])
    max_second_value <- min(second_ranking_values)

    for (config_name in names(all_configurations)) {
        if (second_ranking_values[[config_name]] == max_second_value) {
            return(strsplit(config_name, sep_char)[[1]])
        }
    }
}
