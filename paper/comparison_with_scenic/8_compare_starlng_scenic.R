library(Starlng)
library(ggplot2)
library(dplyr)
library(rhdf5)
library(monocle3)
library(qs2)
library(qualpalr)

if (basename(getwd()) != "comparison_with_scenic") {
    setwd("paper/comparison_with_scenic")
}

mon_obj <- qs_read(file.path("starlng_app_low_nclust_high", "objects", "monocle_object.qs2"), nthreads = 10)
psd <- qs_read(file.path("starlng_app_low_nclust_high", "objects", "recommended_pseudotime.qs2"), nthreads = 10)$recommended_pseudotime
mon_obj <- order_cells(mon_obj, root_cells = colnames(mon_obj)[1])
mon_obj@principal_graph_aux@listData$UMAP$pseudotime <- psd

umap_df <- reducedDim(mon_obj, "UMAP")

expr_matrix_path <- file.path("starlng_app_low_nclust_high", "objects", "expression.h5")
expr_matrix <- rhdf5::h5read(expr_matrix_path, "expression_matrix")
rownames(expr_matrix) <- rhdf5::h5read(expr_matrix_path, "genes")
colnames(expr_matrix) <- rhdf5::h5read(expr_matrix_path, "cells")

results <- list(
    "scenic" = qs_read("scenic_results.qs2", nthreads = 10),
    "starlng_low" = qs_read("starlng_low_results.qs2", nthreads = 10),
    "starlng_high" = qs_read("starlng_high_results.qs2", nthreads = 10)
)

scale_min_max <- function(x) {
    min_val <- min(x)
    max_val <- max(x)

    if (min_val == max_val) {
        return(x)
    }

    return((x - min_val) / (max_val - min_val))
}

within_group_distance <- function(umap_temp_df) {
    if (nrow(umap_temp_df) < 2) {
        return(0)
    }
    dist_matrix <- as.matrix(dist(umap_temp_df))
    return(mean(dist_matrix[lower.tri(dist_matrix)]))
}

mad_z_score <- function(x) {
    if (length(x) < 2) {
        return(rep(FALSE, length(x)))
    }

    med <- stats::median(x, na.rm = TRUE)
    mad_val <- stats::mad(x, center = med, constant = 1, na.rm = TRUE)

    if (mad_val == 0) {
        return(rep(0, length(x)))
    }

    z_scores <- (x - med) / mad_val
    return(z_scores)
}

calculate_umap_average_distance <- function(umap_df, selected_cells = NULL) {
    if (is.null(selected_cells)) {
        selected_cells <- seq_len(nrow(umap_df))
    }

    if (length(selected_cells) < 2) {
        return(0)
    }

    umap_df <- as.matrix(umap_df[selected_cells, ])
    # calculate centroid 
    centroid <- colMeans(umap_df, na.rm = TRUE)
    umap_df[, 1] <- umap_df[, 1] - centroid[1]
    umap_df[, 2] <- umap_df[, 2] - centroid[2]
    umap_df <- umap_df ^ 2
    umap_df <- rowSums(umap_df, na.rm = TRUE) ^ 0.5

    return(median(umap_df))
}

### PSD

ov_psd_stats <- list()

thresh_psd_good <- (quantile(psd, 0.95) - quantile(psd, 0.05)) / 10
thresh_psd_bad <- (quantile(psd, 0.95) - quantile(psd, 0.05)) / 3

for (ob_name in names(results)) {
    ov_psd_stats[[ob_name]] <- list()
    obj <- results[[ob_name]]

    psd_stats <- c()
    umap_stats <- c()
    for (i in seq_along(obj)) {
        obj[[i]]$avg_expr <- scale_min_max(obj[[i]]$avg_expr)
        cell_mask <- obj[[i]]$avg_expr > 0.5
        temp_st <- fivenum(psd[cell_mask])
        psd_stats <- c(psd_stats, temp_st[4] - temp_st[2])
        # umap_stats <- c(umap_stats, within_group_distance(umap_df[cell_mask, , drop = FALSE]))
        umap_stats <- c(umap_stats, calculate_umap_average_distance(umap_df[cell_mask, , drop = FALSE]))
    }
    psd_stats <- round(psd_stats, 3)
    umap_stats <- round(umap_stats, 3)

    is_outlier <- mad_z_score(psd_stats) > 3.5 & mad_z_score(umap_stats) > 3.5
    is_outlier[psd_stats >= thresh_psd_bad] <- TRUE
    ov_psd_stats[[ob_name]]$n_modules <- length(is_outlier)
    ov_psd_stats[[ob_name]]$n_outliers <- sum(is_outlier)
    module_order <- order(psd_stats, umap_stats)

    perc_cells <- c()
    cell_mask <- rep(FALSE, length(psd))
    selected_modules <- c()
    total_ncells <- c()
    for (i in module_order) {
        obj[[i]]$avg_expr <- scale_min_max(obj[[i]]$avg_expr)
        current_mask <- obj[[i]]$avg_expr > 0.5
        total_ncells <- c(total_ncells, sum(current_mask))
        if (is_outlier[i]) {
            next
        }
        temp_mask <- cell_mask | current_mask
        eligible <- TRUE
        ncells <- sum(current_mask)
        if (ncells < 10) {
            next
        }
        iqr <- psd_stats[i]
        if (length(perc_cells) > 0) {
            nunique <- sum(temp_mask) - sum(cell_mask)
            if (ncells < 100 && nunique / ncells < 0.6) {
                eligible <- FALSE
            }
            if (ncells >= 100 && nunique / ncells < 0.2) {
                eligible <- FALSE
            }
            
            if (iqr < thresh_psd_good) {
                eligible <- TRUE
            }
            if (!eligible) {
                next
            }
        }
        cell_mask <- temp_mask
        selected_modules <- c(selected_modules, i)
        perc_cells <- c(perc_cells, sum(cell_mask) / length(psd))
    }
    names(total_ncells) <- module_order
    names(perc_cells) <- selected_modules
    ov_psd_stats[[ob_name]]$total_ncells <- total_ncells
    ov_psd_stats[[ob_name]]$not_eligible <- length(is_outlier) - sum(is_outlier) - length(selected_modules)
    ov_psd_stats[[ob_name]]$selected_modules <- selected_modules
    ov_psd_stats[[ob_name]]$perc_cells <- perc_cells
    ov_psd_stats[[ob_name]]$psd_stats <- psd_stats
    ov_psd_stats[[ob_name]]$umap_stats <- umap_stats
}

df <- rbind(
    data.frame(
        method = names(ov_psd_stats),
        perc_cells = sapply(ov_psd_stats, function(x) x$n_outliers / x$n_modules),
        type_cells = "outliers"
    ),
    data.frame(
        method = names(ov_psd_stats),
        perc_cells = sapply(ov_psd_stats, function(x) x$not_eligible / x$n_modules),
        type_cells = "not_eligible"
    ),
    data.frame(
        method = names(ov_psd_stats),
        perc_cells = sapply(ov_psd_stats, function(x) length(x$selected_modules) / x$n_modules),
        type_cells = "selected"
    )
)

fig_path <- file.path("..", "figures")

pdf(file.path(fig_path, "comparison_usable_modules.pdf"), width = 7, height = 4)
ggplot(df, aes(y = method, x = perc_cells, fill = type_cells)) +
    geom_bar(stat = "identity") +
    labs(title = "", y = "", x = "Proportion of modules") +
    theme_classic() +
    theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.direction = "horizontal",
        legend.position = "bottom"
    ) +
    scale_fill_manual(values = qualpalr::qualpal(3)$hex, name = "Module types")
dev.off()


df <- NULL
for (i in seq_along(ov_psd_stats)) {
    temp_df <- data.frame(
        method = names(ov_psd_stats)[i],
        perc_cells = ov_psd_stats[[i]]$perc_cells,
        iqr = ov_psd_stats[[i]]$psd_stats[ov_psd_stats[[i]]$selected_modules]
    )
    if (is.null(df)) {
        df <- temp_df
    } else{
        df <- rbind(df, temp_df)
    }
}
pdf(file.path(fig_path, "comparison_cell_coverage.pdf"), width = 7, height = 4)
ggplot(df, aes(x = iqr, y = perc_cells, color = method)) +
    geom_point() +
    geom_line() +
    labs(title = "", x = "Pseudotime IQR", y = "% Cell Coverage") +
    scale_colour_manual(values = qualpalr::qualpal(3)$hex, name = "Method") +
    theme_classic() +
    theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)
    )
dev.off()

qs2::qs_save(ov_psd_stats, "overall_pseudotime_stats.qs2")


### ENRICH

library(UpSetR)

method_terms <- lapply(results, function(x) {
    unique(unlist(lapply(x, function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment$term_id)
    })))
})
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_all_unfiltered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale = 2))
dev.off()

method_terms <- lapply(names(results), function(x) {
    unique(unlist(lapply(results[[x]][ov_psd_stats[[x]]$selected_modules], function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment$term_id)
    })))
})
names(method_terms) <- names(results)
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_all_filtered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale = 2))
dev.off()

method_terms <- lapply(names(results), function(x) {
    unique(unlist(lapply(results[[x]][ov_psd_stats[[x]]$selected_modules], function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment %>% filter(source != "TF") %>% pull(term_id))
    })))
})
names(method_terms) <- names(results)
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_not_tf_filtered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale = 2))
dev.off()

method_terms <- lapply(names(results), function(x) {
    unique(unlist(lapply(results[[x]], function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment %>% filter(source != "TF") %>% pull(term_id))
    })))
})
names(method_terms) <- names(results)
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_not_tf_unfiltered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale =2))
dev.off()

# "BICC1" %in% unique(unlist(lapply(results[[1]][ov_psd_stats[[1]]$selected_modules], function(x) x$genes))) 


method_terms <- lapply(names(results), function(x) {
    unique(unlist(lapply(results[[x]][ov_psd_stats[[x]]$selected_modules], function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment %>% filter(source == "GO:BP") %>% pull(term_id))
    })))
})
names(method_terms) <- names(results)
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_gobp_filtered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale = 2))
dev.off()


method_terms <- lapply(names(results), function(x) {
    unique(unlist(lapply(results[[x]], function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment %>% filter(source == "GO:BP") %>% pull(term_id))
    })))
})
names(method_terms) <- names(results)
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_gobp_unfiltered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale = 2))
dev.off()

method_terms <- lapply(names(results), function(x) {
    unique(unlist(lapply(results[[x]][ov_psd_stats[[x]]$selected_modules], function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment %>% filter(source %in% c("KEGG", "REAC")) %>% pull(term_id))
    })))
})
names(method_terms) <- names(results)
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_kegg_reac_filtered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale = 2))
dev.off()

method_terms <- lapply(names(results), function(x) {
    unique(unlist(lapply(results[[x]][ov_psd_stats[[x]]$selected_modules], function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment %>% filter(source %in% c("REAC")) %>% pull(term_id))
    })))
})
names(method_terms) <- names(results)
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_reac_filtered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale = 2))
dev.off()

method_terms <- lapply(names(results), function(x) {
    unique(unlist(lapply(results[[x]][ov_psd_stats[[x]]selected_modules], function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment %>% filter(source %in% c("KEGG")) %>% pull(term_id))
    })))
})
names(method_terms) <- names(results)
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_kegg_filtered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale = 2))
dev.off()

method_terms <- lapply(names(results), function(x) {
    unique(unlist(lapply(results[[x]][ov_psd_stats[[x]]$selected_modules], function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment %>% filter(source %in% c("TF")) %>% pull(term_id))
    })))
})
names(method_terms) <- names(results)
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_tf_filtered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale = 2))
dev.off()


method_terms <- lapply(names(results), function(x) {
    unique(unlist(lapply(results[[x]], function(y) {
        if (is.null(y$enrichment)) {
            return(NULL)
        }
        return(y$enrichment %>% filter(source %in% c("KEGG", "REAC")) %>% pull(term_id))
    })))
})
names(method_terms) <- names(results)
str(method_terms)

# upset plot
pdf(file.path(fig_path, "comparison_upset_kegg_reac_unfiltered.pdf"), width = 5, height = 6)
print(upset(fromList(method_terms), sets = names(method_terms), keep.order = TRUE, text.scale = 2))
dev.off()


### ENRICHMENT PLOTS
ov_psd_stats <- qs_read("overall_pseudotime_stats.qs2")

defined_cells <- list()
for (i in names(results)) {
    defined_cells[[i]] <- list()
    for (m in seq_along(results[[i]])) {
        defined_cells[[i]][[m]] <- colnames(expr_matrix)[scale_min_max(results[[i]][[m]]$avg_expr) > 0.5]
    }
}

psd_ordering <- lapply(names(defined_cells), function(mt) {
    used_mods <- ov_psd_stats[[mt]]$selected_modules
    return(used_mods[order(sapply(used_mods, function(cell) {
        median(psd[defined_cells[[mt]][[cell]]], na.rm = TRUE)
    }))])
})
names(psd_ordering) <- names(defined_cells)

top_terms <- list()
id_to_name_mapping <- c()
for (mt in names(results)) {
    top_terms[[mt]] <- list()
    for (i in ov_psd_stats[[mt]]$selected_modules) {#seq_along(results[[mt]])) {
        if (is.null(results[[mt]][[i]]$enrichment)) {
            next
        }
        top_ter <- results[[mt]][[i]]$enrichment %>%
            filter(source == "GO:BP") %>%
            arrange(p_value) %>%
            head(5)
        i <- as.character(i)
        if (nrow(top_ter) == 0) {
            top_terms[[mt]][[i]] <- NULL
            next
        }
        top_terms[[mt]][[i]] <- top_ter
        top_terms[[mt]][[i]]$module <- i
        
        ids <- top_terms[[mt]][[i]] %>% pull(term_id)
        nms <- top_terms[[mt]][[i]] %>% pull(term_name)
        names(nms) <- ids 
        id_to_name_mapping <- c(id_to_name_mapping, nms)
    }
}
id_to_name_mapping <- id_to_name_mapping[unique(names(id_to_name_mapping))]


for (mt in names(top_terms)) {
    combined_df <- do.call(rbind, top_terms[[mt]]) %>% arrange(p_value)
    top_t <- 25
    i <- 1
    used_terms <- c()
    while (i <= nrow(combined_df) && top_t > 0) {
        current_term <- combined_df[i, "term_id"]
        if (!(current_term %in% used_terms)) {
            used_terms <- c(used_terms, current_term)
            top_t <- top_t - 1
        }

        i <- i + 1
    }
    combined_df$module <- factor(combined_df$module, levels = intersect(psd_ordering[[mt]], unique(combined_df$module)))


    combined_df[,c("module", "term_name")]
    pdf(file.path(fig_path, paste0("enrichment_", mt, ".pdf")), width = 12)
    print(ggplot(combined_df %>% filter(term_id %in% used_terms), aes(x = module, y = term_name, size = -log10(p_value), colour =term_size )) +
        geom_point() +
        theme_bw() +
        scale_colour_viridis_c())
    dev.off()


}
