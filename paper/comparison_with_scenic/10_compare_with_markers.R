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

mks <- read.csv("markers_cell_cluster_8.csv") %>% filter(avg_log2FC > 0)
rownames(mks) <- mks$gene
moran_vals <- read.csv(file.path("starlng_app_low_nclust_high", "objects", "genes_info.csv"))
rownames(moran_vals) <- moran_vals$X
stable_modules <- read.csv(file.path("starlng_app_low_nclust_high", "objects", "stable_modules.csv"))
moran_vals$clusters <- NA
moran_vals[stable_modules$genes, "clusters"] <- stable_modules$stable_modules_37

moran_vals <- moran_vals[union(mks$gene, moran_vals %>% filter(clusters %in% c(15, 25)) %>% rownames), ]
moran_vals$avg_log2FC <- NA
moran_vals[rownames(mks), "avg_log2FC"] <- mks$"avg_log2FC"

moran_vals$clusters[!is.na(moran_vals$clusters) & !(moran_vals$clusters %in% c(15, 25))] <- -1
moran_vals$clusters <- factor(moran_vals$clusters)




ggplot(moran_vals, aes(x = morans_I, y = avg_log2FC, colour = clusters)) +
# ggplot(moran_vals %>% filter(!(clusters %in% c(15, 25))), aes(x = morans_I, y = avg_log2FC, colour = clusters)) +
    geom_point() +
    theme_classic() +
    theme(
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16)
    )

# cumulative percentages of captured genes
logfc_val <- c()
perc <- c()
current_count <- 0
current_passed <- 0
moran_vals <- moran_vals %>% arrange(desc(avg_log2FC))
for (i in seq_len(nrow(moran_vals))) {
    current_logfc <- moran_vals[i, "avg_log2FC"]
    current_cl <- moran_vals[i, "clusters"]
    if (!is.na(current_cl) && current_cl %in% c(15, 25)) {
        current_passed <- current_passed + 1
    }
    current_count <- current_count + 1
    logfc_val <- c(logfc_val, current_logfc)
    perc <- c(perc, current_passed / current_count)
}
perc <- perc[!is.na(logfc_val)]
logfc_val <- logfc_val[!is.na(logfc_val)]



# logfc thresh bins
interv_val <- 0.25
logfc_intervs <- c()
intervs <- list(
    "15" = c(),
    "25" = c(),
    "-1" = c(),
    "NA" = c()
)
n <- 0

for (i in seq(from = 7, to = 0.5, by = -interv_val)) {
    logfc_intervs <- c(logfc_intervs, i)
    temp_moran <- moran_vals %>% filter(!is.na(avg_log2FC), avg_log2FC >= i, avg_log2FC < i + interv_val)
    n_total <- nrow(temp_moran)
    temp_moran <- split(temp_moran, temp_moran$clusters)
    print(i)
    intervs[["15"]] <- c(intervs[["15"]], ifelse(n > 0, intervs[["15"]][n], 0) + nrow(temp_moran[["15"]]))
    intervs[["25"]] <- c(intervs[["25"]], ifelse(n > 0, intervs[["25"]][n], 0) + nrow(temp_moran[["25"]]))
    intervs[["-1"]] <- c(intervs[["-1"]], ifelse(n > 0, intervs[["-1"]][n], 0) + nrow(temp_moran[["-1"]]))
    intervs[["NA"]] <- c(intervs[["NA"]], ifelse(n > 0, intervs[["NA"]][n], 0) + n_total - sum(sapply(temp_moran, nrow)))
    # n <- n + 1
}
names(intervs[["15"]]) <- logfc_intervs
names(intervs[["25"]]) <- logfc_intervs
names(intervs[["-1"]]) <- logfc_intervs
names(intervs[["NA"]]) <- logfc_intervs
intervs

interv_df <- data.frame(
    logfc = rep(logfc_intervs, 4),
    count = c(intervs[["15"]], intervs[["25"]], intervs[["-1"]], intervs[["NA"]]),
    cluster = factor(rep(c("15", "25", "-1", "NA"), each = length(logfc_intervs)))
)

# convert to percentages grouped by logfc
interv_df <- interv_df %>% group_by(logfc) %>%
    mutate(perc = count / sum(count)) %>%
    ungroup()

max_logfc <- min(interv_df$logfc)
# add text with count on the geom bars
plt <- ggplot() +
    geom_bar(
        data = interv_df,
        mapping = aes(x = logfc, y = perc, fill = cluster),
        stat = "identity", 
        position = "stack",
        alpha = 0.5) +
    geom_text(
        data = interv_df %>% filter(count > 0),
        mapping = aes(x = logfc, y = perc, label = count, group = cluster),
        position = position_stack(vjust = 0.25),
        size = 5
    ) +
    geom_line(
        data = data.frame(logfc = logfc_val, perc = perc),
        mapping = aes(x = logfc, y = perc),
        colour = "black", size = 1.5
    ) +
    scale_x_reverse() +
    xlab("log2FC") +
    theme_classic() +
    theme(
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = "bottom",
        legend.direction = "horizontal"
    ) 
plt  
    ggsave(
        filename = file.path(fig_path, "marker_moran_logfc_comparison.pdf"),
        plt,
        width = 8, height = 6
    ) 

# ggplot() + 
#     geom_line(
#         data = data.frame(logfc = logfc_val, perc = perc),
#         mapping = aes(x = logfc, y = perc)) +
#     geom_bar(
#         data = data.frame(
#             logfc = rep(logfc_intervs, 4),
#             count = c(intervs[["15"]], intervs[["25"]], intervs[["-1"]], intervs[["NA"]]),
#             cluster = factor(rep(c("15", "25", "-1", "NA"), each = length(logfc_intervs)))
#         ),
#         mapping = aes(x = logfc, y = count / max(count), fill = cluster),
#         stat = "identity", position = "stack", alpha = 0.5
#     ) +
#     theme_classic() +
#     xlab("log2FC threshold") +
#     ylab("Percentage of captured markers being used in clustering") +
#     theme(
#         axis.text = element_text(size = 16),
#         legend.text = element_text(size = 16),
#         axis.title = element_text(size = 16),
#         legend.title = element_text(size = 16)
#     )
