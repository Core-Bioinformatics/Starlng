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

fig_path <- file.path('../figures')
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

ov_psd_stats <- qs_read("overall_pseudotime_stats.qs2")

scale_min_max <- function(x) {
    min_val <- min(x)
    max_val <- max(x)

    if (min_val == max_val) {
        return(x)
    }

    return((x - min_val) / (max_val - min_val))
}

get_gene_umap <- function(expr_matrix, result_list, selected_modules = NULL) {
    if (is.null(selected_modules)) {
        selected_modules <- names(result_list)
    }
    df <- NULL
    for (i in selected_modules) {
        temp_df <- data.frame(genes = result_list[[i]]$genes)
        temp_df$module <- i
        if (is.null(df)) {
            df <- temp_df
        } else {
            df <- rbind(df, temp_df)
        }
    }
    cols <- qualpal(length(result_list))$hex
    names(cols) <- names(result_list)
    df$module <- factor(df$module, selected_modules)
    genes <- intersect(df$genes, rownames(expr_matrix))
    expr_matrix <- expr_matrix[genes, , drop = FALSE]

    fload <- get_feature_loading(expr_matrix)
    umap_emb <- uwot::umap(fload, n_neighbors = 30, min_dist = 0.3, seed = 42)
    df$UMAP1 <- umap_emb[,1]
    df$UMAP2 <- umap_emb[,2]

    text_position <- df %>% group_by(module) %>% summarise(x_text = mean(UMAP1), y_text = mean(UMAP2))

    return(
        ggplot() +
            geom_point(data = df, mapping = aes(x = UMAP1, y = UMAP2, colour = module)) +
            scale_colour_manual(values = cols[levels(df$module)]) +
            geom_text(data = text_position, aes(x = x_text, y = y_text, label = module), size = 10) +
            theme_classic() +
            theme(
                legend.position = "none"
            )
    )
}

umap_mod <- get_gene_umap(expr_matrix, results$starlng_high, ov_psd_stats$starlng_high$selected_modules)
pdf(file.path(fig_path, "gene_umap_starlng_high_filtered.pdf"), height = 7, width = 7)
umap_mod
dev.off()

umap_mod <- get_gene_umap(expr_matrix, results$starlng_high)
pdf(file.path(fig_path, "gene_umap_starlng_high_unfiltered.pdf"), height = 7, width = 7)
umap_mod
dev.off()

umap_mod <- get_gene_umap(expr_matrix, results$starlng_low, ov_psd_stats$starlng_low$selected_modules)
pdf(file.path(fig_path, "gene_umap_starlng_low_filtered.pdf"), height = 7, width = 7)
umap_mod
dev.off()

umap_mod <- get_gene_umap(expr_matrix, results$starlng_low)
pdf(file.path(fig_path, "gene_umap_starlng_low_unfiltered.pdf"), height = 7, width = 7)
umap_mod
dev.off()

# PSD

defined_cells <- list()
for (i in names(results)) {
    defined_cells[[i]] <- list()
    for (m in names(results[[i]])) {
        defined_cells[[i]][[m]] <- colnames(expr_matrix)[scale_min_max(results[[i]][[m]]$avg_expr) > 0.5]
    }
}

psd <- qs_read(file.path("starlng_app_low_nclust_high", "objects", "recommended_pseudotime.qs2"), nthreads = 10)$recommended_pseudotime

psd_violin_plot <- function(psd_vec, module_cell_list, selected_modules = NULL) {
    if (is.null(selected_modules)) {
        selected_modules <- names(module_cell_list)
    }
    df <- NULL
    for (i in selected_modules) {
        temp_df <- data.frame(pseudotime = psd_vec[module_cell_list[[i]]])
        temp_df$module <- i
        if (is.null(df)) {
            df <- temp_df
        } else {
            df <- rbind(df, temp_df)
        }
    }
    cols <- qualpal(length(module_cell_list))$hex
    names(cols) <- names(module_cell_list)
    df$module <- factor(df$module, selected_modules)
    mds <- df %>% group_by(module) %>% summarise(med = median(pseudotime)) %>% arrange(med)
    df$module <- factor(df$module, levels = mds$module)

    return(
        ggplot(df, aes(x = module, y = pseudotime, fill = module)) +
            geom_violin(scale = "width", draw_quantiles = c(0.25, 0.5, 0.75), adjust = 2, size = 0.2) +
            theme_classic() +
            theme(
                legend.position = "none",
                axis.text = element_text(size = 18),
                axis.title = element_text(size = 18),
                axis.text.x = element_text(angle = 45, hjust = 1)
            ) +
            labs(x = "Module", y = "Pseudotime") +
            scale_fill_manual(values = cols[levels(df$module)])
    )
}

pdf(file.path(fig_path, "psd_violin_starlng_high_filtered.pdf"), height = 7, width = 10)
psd_violin_plot(psd, defined_cells$starlng_high, ov_psd_stats$starlng_high$selected_modules)
dev.off()

pdf(file.path(fig_path, "psd_violin_starlng_high_unfiltered.pdf"), height = 7, width = 10)
psd_violin_plot(psd, defined_cells$starlng_high)
dev.off()

pdf(file.path(fig_path, "psd_violin_starlng_low_filtered.pdf"), height = 7, width = 10)
psd_violin_plot(psd, defined_cells$starlng_low, ov_psd_stats$starlng_low$selected_modules)
dev.off()

pdf(file.path(fig_path, "psd_violin_starlng_low_unfiltered.pdf"), height = 7, width = 10)
psd_violin_plot(psd, defined_cells$starlng_low)
dev.off()

# stability plot
type_app <- "low"
stab <- qs_read(file.path(paste0("starlng_app_low_nclust_", type_app), "objects", "full_stability_assessment.qs2"))
df <- do.call(
    rbind,
    lapply(names(stab), function(x) {
        data.frame(
            ecc = stab[[x]][[1]]$overall_ecc,
            k = rep(x, length(stab[[x]][[1]]$overall_ecc))
        )
    })
)
df$k <- factor(df$k, levels = rev(names(stab)))

pdf(file.path(fig_path, paste0("stability_plot_neigh_", type_app, ".pdf")))
ggplot(df, aes(x = k, y = ecc)) +
    geom_boxplot() +
    theme_bw() +
    theme(
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 18)
    ) +
    labs(x = "number of nearest neighbours", y = "Overall ECC")
dev.off()

df <- do.call(
    rbind,
    lapply(names(stab$"40"[[1]]$k), function(x) {
        data.frame(
            mean_ecc = mean(stab$"40"[[1]]$k[[x]]$ecc),
            iqr_ecc = IQR(stab$"40"[[1]]$k[[x]]$ecc),
            freq = sum(sapply(stab$"40"[[1]]$k[[x]]$partitions, function(x) x$freq)),
            k = x
        )
    })
)
rownames(df) <- df$k

df$k <- factor(df$k, levels = names(stab$"40"[[1]]$k))
cap_val <- 200
df$freq[df$freq > cap_val] <- cap_val

pdf(file.path(fig_path, paste0("stability_plot_k_", type_app, ".pdf")), width = 18)
print(ggplot() +
    geom_point(data = df %>% filter(mean_ecc <= 0.9), aes(x = k, colour = iqr_ecc, y = freq), size = 1) +
    geom_point(data = df %>% filter(mean_ecc > 0.9), aes(x = k, colour = iqr_ecc, y = freq, size = mean_ecc)) +
    theme_bw() +
    scale_colour_viridis_c(direction = -1) +
    scale_size_continuous(range = c(4, 7)) +
    geom_hline(yintercept = 30, linetype = "dashed", colour = "red") +
    scale_x_discrete(
        limits = levels(df$k),
        labels = sapply(rownames(df), function(x) {
            if (df[x, "mean_ecc"] < 0.9 || df[x, "freq"] < 30) {
                return(paste0("<span style='color:red'>", x, "</span>"))
            }
            return(x)
        })
    ) +
    theme(
        axis.text = element_text(size = 18),
        axis.text.x = ggtext::element_markdown(),
        axis.title = element_text(size = 18)
    )
)
dev.off()


# umap with gene - induced clusters

umap_df <- data.frame(umap_df)
umap_df$clusters <- NA

defined_cells <- list()
for (i in names(results)) {
    defined_cells[[i]] <- list()
    for (m in names(results[[i]])) {
        defined_cells[[i]][[m]] <- colnames(expr_matrix)[scale_min_max(results[[i]][[m]]$avg_expr) > 0.5]
    }
}
defined_cells[["starlng_high"]][["5_31"]] <- intersect(defined_cells[["starlng_high"]][["5"]], defined_cells[["starlng_high"]][["31"]])
defined_cells[["starlng_high"]][["9"]] <- setdiff(defined_cells[["starlng_high"]][["9"]], defined_cells[["starlng_high"]][["5_31"]])

defined_cells[["starlng_high"]][["24_6_1_32"]] <- union(intersect(intersect(defined_cells[["starlng_high"]][["24"]], defined_cells[["starlng_high"]][["6"]]), defined_cells[["starlng_high"]][["1"]]), defined_cells[["starlng_high"]][["32"]])
defined_cells[["starlng_high"]][["10_17_13"]] <- setdiff(union(union(defined_cells[["starlng_high"]][["10"]], defined_cells[["starlng_high"]][["17"]]), defined_cells[["starlng_high"]][["13"]]), defined_cells[["starlng_high"]][["24_6_1_32"]])
defined_cells[["starlng_high"]][["8_14"]] <- union(defined_cells[["starlng_high"]][["8"]], defined_cells[["starlng_high"]][["14"]])


module_order <- c("29", "15", "25", "23", "2", "22", "9", "4", "5_31", "26",
"24_6_1_32", "10_17_13",
"11", "19", "20", "21", "18", "27", "8_14"
# "10", "32", "24", "17", "13", "33", "6", "1"
)
actual_module_order <- c()

for(i in module_order) {
    current_cells <- defined_cells$starlng_high[[i]]
    mask <- rownames(umap_df) %in% current_cells & is.na(umap_df$clusters)
    if (sum(mask) == 0) {
        next
    }
    actual_module_order <- c(actual_module_order, i)
    umap_df$clusters[mask] <- i
}

umap_df$clusters <- factor(umap_df$clusters, levels = actual_module_order)
text_loc <- umap_df %>% group_by(clusters) %>% summarise(x = mean(UMAP_1), y = mean(UMAP_2))

pdf(file.path(fig_path, "umap_gene_induced_clusters_starlng_high.pdf"), width = 8, height = 7)
ggplot() +
    geom_point(data = umap_df %>% filter(is.na(clusters)), aes(x = UMAP_1, y = UMAP_2), colour = "gray90", size = 0.75) +
    geom_point(data = umap_df %>% filter(!is.na(clusters)), aes(x = UMAP_1, y = UMAP_2, colour = clusters), size = 1.1) +
    geom_text(data = text_loc, aes(x = x, y = y, label = clusters), size = 6) +
    scale_colour_manual(values = qualpal(length(actual_module_order), bg = "grey90", cvd = c(protan = 0.5, deutan = 0.5, tritan = 0.3))$hex) +
    theme_classic() 
dev.off()

### Low filtering

umap_df <- data.frame(umap_df)
umap_df$clusters <- NA

defined_cells <- list()
for (i in names(results)) {
    defined_cells[[i]] <- list()
    for (m in names(results[[i]])) {
        defined_cells[[i]][[m]] <- colnames(expr_matrix)[scale_min_max(results[[i]][[m]]$avg_expr) > 0.5]
    }
}

defined_cells[["starlng_low"]][["6_14"]] <- union(defined_cells[["starlng_low"]][["6"]], defined_cells[["starlng_low"]][["14"]])
defined_cells[["starlng_low"]][["27"]] <- setdiff(defined_cells[["starlng_low"]][["27"]], defined_cells[["starlng_low"]][["6_14"]])
defined_cells[["starlng_low"]][["33_13"]] <- setdiff(
    union(defined_cells[["starlng_low"]][["33"]], defined_cells[["starlng_low"]][["13"]]),
    unique(unlist(defined_cells[["starlng_low"]][c("12", "7", "25", "24", "26", "7")]))
)

module_order <- c("31", "21", "28", "16", "20", "10",
"27", "6_14",
"30", "15", "22", "8", "5",
"12", "25", "24", "26", "7", "33_13"
)
actual_module_order <- c()

for(i in module_order) {
    current_cells <- defined_cells$starlng_low[[i]]
    mask <- rownames(umap_df) %in% current_cells & is.na(umap_df$clusters)
    if (sum(mask) == 0) {
        print(i)
        next
    }
    actual_module_order <- c(actual_module_order, i)
    umap_df$clusters[mask] <- i
}

umap_df$clusters <- factor(umap_df$clusters, levels = actual_module_order)
text_loc <- umap_df %>% group_by(clusters) %>% summarise(x = mean(UMAP_1), y = mean(UMAP_2))

pdf(file.path(fig_path, "umap_gene_induced_clusters_starlng_low.pdf"), width = 8, height = 7)
ggplot() +
    geom_point(data = umap_df %>% filter(is.na(clusters)), aes(x = UMAP_1, y = UMAP_2), colour = "gray90", size = 0.75) +
    geom_point(data = umap_df %>% filter(!is.na(clusters)), aes(x = UMAP_1, y = UMAP_2, colour = clusters), size = 1.1) +
    geom_text(data = text_loc, aes(x = x, y = y, label = clusters), size = 6) +
    scale_colour_manual(values = qualpal(length(actual_module_order), bg = "grey90", cvd = c(protan = 0.5, deutan = 0.5, tritan = 0.3))$hex) +
    theme_classic() 
dev.off()



