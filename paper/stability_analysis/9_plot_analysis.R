library(qs)
library(Starlng)
library(dplyr)
library(ggplot2)
library(UpSetR)

setwd(file.path(getwd(), "paper", "stability_analysis"))
starlng_clusters <- qread("immune_starlng_cluster_outputs.qs")

########### SUPP - Consistency of clusters picked by Starlng #######
starlng_cluster_list <- list(freqs = list(), parts = list())
for (i in seq_len(30)) {
    current_k <- names(starlng_clusters[[i]]$k)
    for (k in current_k) {
        if (k %in% names(starlng_cluster_list$freqs)) {
            starlng_cluster_list$freqs[[k]] <- starlng_cluster_list$freqs[[k]] + 1
            starlng_cluster_list$parts[[k]] <- c(starlng_cluster_list$parts[[k]], list(starlng_clusters[[i]]$k[[k]]$partitions[[1]]$mb))
        } else {
            starlng_cluster_list$freqs[[k]] <- 1
            starlng_cluster_list$parts[[k]] <- list(starlng_clusters[[i]]$k[[k]]$partitions[[1]]$mb)
        }
    }
}

for (k in names(starlng_cluster_list$parts)) {
    starlng_cluster_list$ecc[[k]] <- ClustAssess::element_consistency(starlng_cluster_list$parts[[k]])
}

# plot as boxplot the distribution of ECC by k, by showing above the boxplot the number in freqs just once
ecc_df <- data.frame(ECC = unlist(starlng_cluster_list$ecc), k = rep(names(starlng_cluster_list$ecc), sapply(starlng_cluster_list$ecc, length)))
ecc_df$k <- factor(ecc_df$k, levels = names(starlng_cluster_list$freqs))

supp_starlng_k_ecc <- ggplot(ecc_df, aes(x = k, y = ECC)) +
    geom_boxplot() +
    geom_text(data = data.frame(x = names(starlng_cluster_list$freqs), y = 1.05, label = unlist(starlng_cluster_list$freqs)), aes(x = x, y = y, label = label), vjust = -0.5) +
    theme_bw()
supp_starlng_k_ecc

########### MAIN - compare the stability of clusters - Starlng vs others ###########
monocle_output <- qread("immune_monocle_cluster_outputs.qs")
hotspot_top500 <- read.csv("immune_hotspot500_modules.csv", header = FALSE, row.names = 1)
hotspot_top1000 <- read.csv("immune_hotspot1000_modules.csv", header = FALSE, row.names = 1)
hotspot_top1500 <- read.csv("immune_hotspot1500_modules.csv", header = FALSE, row.names = 1)

methods_ecc <- list(
    "Starlng best" = starlng_cluster_list$ecc$"4",
    "Starlng overall" = Reduce("+", starlng_cluster_list$ecc) / length(starlng_cluster_list$ecc),
    "Starlng worst" = starlng_cluster_list$ecc$"17",
    "Monocle" = ClustAssess::element_consistency(monocle_output),
    "Hotspot top 500" = ClustAssess::element_consistency(lapply(seq_len(ncol(hotspot_top500)), function(i) hotspot_top500[,i])),
    "Hotspot top 1000" = ClustAssess::element_consistency(lapply(seq_len(ncol(hotspot_top1000)), function(i) hotspot_top1000[,i])),
    "Hotspot top 1500" = ClustAssess::element_consistency(lapply(seq_len(ncol(hotspot_top1500)), function(i) hotspot_top1500[,i]))
)

method_names <- names(methods_ecc)

methods_ecc_df <- data.frame(
    ECC = unlist(methods_ecc),
    method = do.call(c, sapply(method_names, function(x) rep(x, length(methods_ecc[[x]]))))
)
methods_ecc_df$method <- factor(methods_ecc_df$method, levels = method_names)

main_methods_ecc <- ggplot(methods_ecc_df, aes(x = method, y = ECC)) +
    geom_boxplot() +
    theme_bw()
main_methods_ecc

########### SUPP - Genes UpsetPlot #######
gene_list <- list(
    "Starlng" = rownames(qread("immune_starlng.qs")$embedding),
    "Monocle" = rownames(qread("immune_monocle.qs")),
    "Hotspot top 500" = rownames(hotspot_top500),
    "Hotspot top 1000" = rownames(hotspot_top1000),
    "Hotspot top 1500" = rownames(hotspot_top1500)
)

all_genes <- unique(unlist(gene_list))

gene_df <- data.frame(gene = all_genes)
for (method in names(gene_list)) {
    gene_df[[method]] <- as.numeric(all_genes %in% gene_list[[method]])
}

gene_upset <- upset(gene_df, nsets = length(gene_list), nintersects = 40, order.by = "freq")
gene_upset
dev.off()

########### SUPP - Enrichment UpsetPlot #######




########### MAIN - ECC per k ###########
starlng_app_clusters <- qs::qread("/mnt/d/Bioinf_data/Pipelines/starlng_app_immune/objects/full_stability_assessment.qs")
best_config <- qs::qread("/mnt/d/Bioinf_data/Pipelines/starlng_app_immune/objects/best_configuration.qs")
starlng_app_clusters <- starlng_app_clusters[[best_config[1]]][[best_config[2]]]
k_vals <- names(starlng_app_clusters$k)
cl_ecc_df <- data.frame(
    ECC = unlist(lapply(k_vals, function(k) starlng_app_clusters$k[[k]]$ecc)),
    k = rep(k_vals, each = length(starlng_app_clusters$k[[k_vals[1]]]$ecc))
)
cl_ecc_df$k <- factor(cl_ecc_df$k, levels = k_vals)

k_freqs <- sapply(k_vals, function(k) sum(sapply(starlng_app_clusters$k[[k]]$partitions, function(part) part$freq)))

ggplot(cl_ecc_df, aes(x = k, y = ECC)) +
    geom_boxplot() +
    geom_text(data = data.frame(x = k_vals, y = 1.05, label = k_freqs), aes(x = x, y = y, label = label), vjust = -0.5) +
    theme_bw()

########### MAIN + SUPP - Enriched terms ###########
starlng_modules <- read.csv("/mnt/d/Bioinf_data/Pipelines/starlng_app_immune/objects/stable_modules.csv", header = TRUE, row.names = 1)["stable_modules_5"]
starlng_modules <- split(rownames(starlng_modules), starlng_modules$stable_modules_5)

modules <- list(
    "Starlng" = starlng_modules,
    "Monocle" = split(gene_list[["Monocle"]], as.numeric(monocle_output[[1]])),
    "Hotspot" = split(gene_list[["Hotspot top 500"]], as.numeric(hotspot_top500[[1]]))
)

enrichments <- lapply(names(modules), function(method) {
    enrichments <- lapply(modules[[method]], function(module_genes) {
        gprofiler2::gost(module_genes, organism = "hsapiens", sources = "GO:BP")$result[, c("p_value", "term_size", "precision", "recall", "term_name")]
    })
    names(enrichments) <- names(modules[[method]])
    enrichments
})
names(enrichments) <- names(modules)

top_max <- 5
enrichments_top <- lapply(enrichments, function(method_enrichments) {
    lapply(method_enrichments, function(module_enrichments) {
        if (is.null(module_enrichments)) {
            return(NULL)
        }
        ntop <- min(top_max, nrow(module_enrichments))
        module_enrichments[order(module_enrichments$p_value)[seq_len(ntop)], ]
    })
})
names(enrichments_top) <- names(enrichments)

enrich_dotplots <- list()
for (method_name in names(enrichments_top)) {
    all_terms <- unique(unlist(lapply(enrichments_top[[method_name]], function(x) x$term_name)))

    enrich_df <- do.call(rbind, lapply(names(enrichments_top[[method_name]]), function(module_name) {
        module_enrichments <- enrichments_top[[method_name]][[module_name]]
        if (is.null(module_enrichments)) {
            return(NULL)
        }
        module_enrichments$module_name <- module_name
        return(module_enrichments)
    }))

    # enrich_df$term_name <- factor(enrich_df$term_name, levels = all_terms)
    enrich_df$module_name <- factor(enrich_df$module_name, levels = names(enrichments_top[[method_name]]))
    enrich_df$p_value <- -log10(enrich_df$p_value)

    enrich_dotplots[[method_name]] <- ggplot(enrich_df, aes(y = term_name, x = module_name, size = p_value, fill = recall)) +
        geom_point(shape = 21) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))

}


enrich_dotplots[[1]]



