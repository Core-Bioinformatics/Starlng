library(ggplot2)
library(microbenchmark)
library(reshape2)

ncells <- format(seq(from = 1e4, to = 1e5, by = 1e4), scientific = FALSE, justify = "none", trim = TRUE)
ncores <- c(1, seq(from = 5, to = 30, by = 5))
dts_name <- "immune"

setwd("paper/benchmarking")

######### Time across ncells - Pipeline #############
first_df <- TRUE
for (ncell in ncells) {
    mcb_res <- qs::qread(file.path(dts_name, paste0("mcb_res_", ncell, ".qs")), nthreads = 30)
    mcb_res <- mcb_res$time / (1e9 * 60)
    for (i in 2:10) {
        mcb_res[i] <- mcb_res[i] - 5 / 60 # remove the 5 seconds spent on overwriting the h5 file
    }

    if (first_df) {
        mcb_df <- data.frame(ncells = rep(ncell, 10), time = mcb_res)
        first_df <- FALSE
    } else {
        mcb_df <- rbind(mcb_df, data.frame(ncells = rep(ncell, 10), time = mcb_res))
    }
}

mcb_df$ncells <- factor(mcb_df$ncells, levels = ncells)

time_ncells <- ggplot(mcb_df, aes(x = ncells, y = time)) +
    geom_boxplot() +
    geom_jitter(width = 0.15) +
    # draw a line through the median
    stat_summary(fun = median, geom = "line", aes(group = 1), color = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Number of cells", y = "Time (minutes)") +
    ggtitle("Time taken to write Starlng app for different number of cells") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("time_ncells.pdf", time_ncells, width = 8, height = 6)
print(time_ncells)

######## Memory across ncells - Pipeline ###########
first_df <- TRUE
for (ncell in ncells) {
    memory_usage <- read.table(file.path(dts_name, paste0(dts_name, "_", ncell, "_memory.txt")))
    baseline_val <- median(memory_usage$V1[1:3])
    memory_usage <- (memory_usage$V1[-(1:3)] - baseline_val) / (1024 ^ 2)

    if (first_df) {
        first_df <- FALSE
        memory_df <- data.frame(ncells = rep(ncell, length(memory_usage)), memory = memory_usage)
    } else {
        memory_df <- rbind(memory_df, data.frame(ncells = rep(ncell, length(memory_usage)), memory = memory_usage)
        )
    }
}
memory_df$ncells <- factor(memory_df$ncells, levels = ncells)

memory_ncells <- ggplot(memory_df, aes(x = ncells, y = memory)) +
    geom_boxplot() +
    # geom_jitter(width = 0.15) +
    # draw a line through the median
    stat_summary(fun = max, geom = "line", aes(group = 1), color = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Number of cells", y = "Memory usage (GB)") +
    ggtitle("Memory usage of Starlng app for different number of cells") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("memory_ncells.pdf", memory_ncells, width = 8, height = 6)
print(memory_ncells)

####### Time across ncores ########
first_df <- TRUE
for (ncore in ncores) {
    mcb_res <- qs::qread(file.path(dts_name, paste0("mcb_res_100000_ncores_", ncore, "_for_graph_learn_", ncore, ".qs")), nthreads = 30)
    mcb_res <- mcb_res$time / (1e9 * 60)
    for (i in 1:10) {
        mcb_res[i] <- mcb_res[i] - 5 / 60 # remove the 5 seconds spent on overwriting the h5 file
    }

    if (first_df) {
        mcb_df <- data.frame(ncores = rep(ncore, 10), time = mcb_res)
        first_df <- FALSE
    } else {
        mcb_df <- rbind(mcb_df, data.frame(ncores = rep(ncore, 10), time = mcb_res))
    }
}

mcb_df$ncores <- factor(mcb_df$ncores, levels = ncores)

time_ncores <- ggplot(mcb_df, aes(x = ncores, y = time)) +
    geom_boxplot() +
    geom_jitter(width = 0.15) +
    # draw a line through the median
    stat_summary(fun = median, geom = "line", aes(group = 1), color = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Number of cores", y = "Time (minutes)") +
    ggtitle("Time taken to write Starlng app for different number of cores") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("time_ncores.pdf", time_ncores, width = 8, height = 6)
time_ncores

####### Memory across ncores ########
first_df <- TRUE
for (ncore in ncores) {
    memory_usage <- read.table(file.path(dts_name, paste0(dts_name, "_100000_memory_" , ncore, "_cores_", ncore, "_for_graph_learn.txt")))
    baseline_val <- median(memory_usage$V1[1:3])
    memory_usage <- (memory_usage$V1[-(1:3)] - baseline_val) / (1024 ^ 2)

    if (first_df) {
        first_df <- FALSE
        memory_df <- data.frame(ncores = rep(ncore, length(memory_usage)), memory = memory_usage)
    } else {
        memory_df <- rbind(memory_df, data.frame(ncores = rep(ncore, length(memory_usage)), memory = memory_usage)
        )
    }
}

memory_df$ncores <- factor(memory_df$ncores, levels = ncores)

memory_ncores <- ggplot(memory_df, aes(x = ncores, y = memory)) +
    geom_boxplot() +
    # geom_jitter(width = 0.15) +
    # draw a line through the median
    stat_summary(fun = max, geom = "line", aes(group = 1), color = "red") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Number of cores", y = "Memory usage (GB)") +
    ggtitle("Memory usage of Starlng app for different number of cores") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("memory_ncores.pdf", memory_ncores, width = 8, height = 6)
memory_ncores

####### 30-30 vs 30-1 vs 1-1 - Time Memory ######
first_usage <- read.table(file.path("immune", "immune_100000_memory_30_cores_30_for_graph_learn.txt"))
baseline_val <- median(first_usage$V1[1:3])
first_usage <- (first_usage$V1[-(1:3)] - baseline_val) / (1024 ^ 2)
first_df <- data.frame(configuration = rep("30 overall, 30 graph test", length(first_usage)), memory = first_usage)

second_usage <- read.table(file.path("immune", "immune_100000_memory_30_cores_1_for_graph_learn.txt"))
baseline_val <- median(second_usage$V1[1:3])
second_usage <- (second_usage$V1[-(1:3)] - baseline_val) / (1024 ^ 2)
second_df <- data.frame(configuration = rep("30 overall, 1 graph test", length(second_usage)), memory = second_usage)

third_usage <- read.table(file.path("immune", "immune_100000_memory_1_cores_1_for_graph_learn.txt"))
baseline_val <- median(third_usage$V1[1:3])
third_usage <- (third_usage$V1[-(1:3)] - baseline_val) / (1024 ^ 2)
third_df <- data.frame(configuration = rep("1 overall, 1 graph test", length(third_usage)), memory = third_usage)

memory_compar_df <- rbind(first_df, second_df, third_df)
memory_compar_df$configuration <- factor(memory_compar_df$configuration, levels = c("30 overall, 30 graph test", "30 overall, 1 graph test", "1 overall, 1 graph test"))

memory_comparison <- ggplot(memory_compar_df, aes(x = configuration, y = memory, colour = configuration)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Configuration", y = "Memory usage (GB)") +
    ggtitle("Memory usage of Starlng app for different configurations") +
    theme(plot.title = element_text(hjust = 0.5))
memory_comparison


first_time <- qs::qread(file.path("immune", "mcb_res_100000_ncores_30_for_graph_learn_30.qs"))
first_time <- first_time$time / (1e9 * 60)
for (i in 1:10) {
    first_time[i] <- first_time[i] - 5 / 60 # remove the 5 seconds spent on overwriting the h5 file
}
first_df <- data.frame(configuration = rep("30 overall, 30 graph test", length(first_time)), time = first_time)

second_time <- qs::qread(file.path("immune", "mcb_res_100000_ncores_30_for_graph_learn_1.qs"))
second_time <- second_time$time / (1e9 * 60)
for (i in 1:10) {
    second_time[i] <- second_time[i] - 5 / 60 # remove the 5 seconds spent on overwriting the h5 file
}
second_df <- data.frame(configuration = rep("30 overall, 1 graph test", length(second_time)), time = second_time)

third_time <- qs::qread(file.path("immune", "mcb_res_100000_ncores_1_for_graph_learn_1.qs"))
third_time <- third_time$time / (1e9 * 60)
for (i in 1:10) {
    third_time[i] <- third_time[i] - 5 / 60 # remove the 5 seconds spent on overwriting the h5 file
}
third_df <- data.frame(configuration = rep("1 overall, 1 graph test", length(third_time)), time = third_time)

time_compar_df <- rbind(first_df, second_df, third_df)
time_compar_df$configuration <- factor(time_compar_df$configuration, levels = c("30 overall, 30 graph test", "30 overall, 1 graph test", "1 overall, 1 graph test"))

time_comparison <- ggplot(time_compar_df, aes(x = configuration, y = time, colour = configuration)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Configuration", y = "Time (minutes)") +
    ggtitle("Time taken to write Starlng app for different configurations") +
    theme(plot.title = element_text(hjust = 0.5))
time_comparison


######### APP - Runtime #########
first_df <- TRUE
for (ncell in ncells) {
    current_memory_usage <- read.table(file.path("immune_app", paste0("memory_", ncell, ".txt")))$V1 / (1024 ^ 2)
    baseline_memory <- median(current_memory_usage[1:3])
    current_memory_usage <- current_memory_usage[-(1:3)] - baseline_memory

    if (first_df) {
        first_df <- FALSE
        memory_df <- data.frame(ncells = rep(ncell, length(current_memory_usage)), memory = current_memory_usage)
    } else {
        memory_df <- rbind(memory_df, data.frame(ncells = rep(ncell, length(current_memory_usage)), memory = current_memory_usage))
    }
}
memory_df$ncells <- factor(memory_df$ncells, levels = ncells)

memory_app <- ggplot(memory_df, aes(x = ncells, y = memory)) +
    geom_boxplot() +
    theme_bw() +
    geom_smooth(aes(group = 1), method = "loess", se = TRUE, color = "black", fill = "grey", level = 0.99) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Number of cells", y = "Memory usage (GB)") +
    ggtitle("Memory usage of Starlng app for different number of cells") +
    theme(plot.title = element_text(hjust = 0.5)
    )
memory_app

######### APP - Memory #########
time_to_minutes <- function(x) {
    x <- as.character(x)
    x <- strsplit(x, ":")[[1]]
    x <- as.numeric(x)
    x <- x[1] * 60 + x[2] + x[3] / 60
    return(x)
}

first_df <- TRUE 
for (ncell in ncells) {
    current_time <- read.csv(file.path("immune_app", paste0("immune_runtime_app_", ncell, ".csv")), header = FALSE)
    current_time <- sapply(current_time$V2, function(x) { time_to_minutes(x) })

    if (first_df) {
        first_df <- FALSE
        time_df <- data.frame(ncells = rep(ncell, 1), clustering_runtime = current_time[3], heatmap_runtime = current_time[5], total_runtime = sum(current_time))
    } else {
        time_df <- rbind(time_df, data.frame(ncells = rep(ncell, 1), clustering_runtime = current_time[3], heatmap_runtime = current_time[5], total_runtime = sum(current_time)))
    }
}

# barplot on three categories
time_df$ncells <- factor(time_df$ncells, levels = ncells)
time_df_melt <- melt(time_df, id.vars = "ncells")
colnames(time_df_melt) <- c("ncells", "category", "runtime")
time_df_melt$category <- factor(time_df_melt$category, levels = c("clustering_runtime", "heatmap_runtime", "total_runtime"))

ggplot(time_df_melt, aes(x = ncells, y = runtime, fill = category)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Number of cells", y = "Runtime (minutes)", fill = "Executed component") +
    ggtitle("Runtime of Starlng app for different number of cells") +
    theme(plot.title = element_text(hjust = 0.5))
