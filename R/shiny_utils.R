save_filetypes <- list(
    "PDF" = grDevices::pdf,
    "PNG" = function(filename, width, height) {
        ragg::agg_png(filename, width, height, units = "in")
    },
    "SVG" = grDevices::svg
)

gear_umaps <- function(ns, id, default_order = "default", contains_continuous = TRUE, scale_default = FALSE) {
    shinyWidgets::dropdownButton(
        inputId = ns(paste0(id, "_gear_umap")),
        shiny::splitLayout(
            shiny::tagList(
                shiny::sliderInput(
                    inputId = ns(paste0(id, "_text_size")),
                    label = "Text size",
                    min = 5.00, max = 30.00, value = 10.00, step = 0.5
                ),
                shiny::sliderInput(
                    inputId = ns(paste0(id, "_axis_size")),
                    label = "Axis labels size",
                    min = 5.00, max = 30.00, value = 10.00, step = 0.5
                ),
                shiny::sliderInput(
                    inputId = ns(paste0(id, "_legend_size")),
                    label = "Legend size",
                    min = 5.00, max = 30.00, value = 10.00, step = 0.5
                ),
                shiny::sliderInput(
                    inputId = ns(paste0(id, "_pt_alpha")),
                    label = "Colour alpha",
                    min = 0.00, max = 1.00, value = 0.8, step = 0.05
                )
            ),
            shiny::tagList(
                shiny::sliderInput(
                    inputId = ns(paste0(id, "_pt_size")),
                    label = "Point size",
                    min = 0.05, max = 5.00, value = 0.60, step = 0.05
                ),
                if (contains_continuous) {
                    shinyWidgets::radioGroupButtons(
                        inputId = ns(paste0(id, "_pt_order")),
                        label = "Point ordering",
                        choices = c("default", "highest", "lowest"),
                        selected = default_order
                    )
                },
                shiny::sliderInput(
                    inputId = ns(paste0(id, "_trajectory_width")),
                    label = "Trajectory width",
                    min = 0.00, max = 5.00, value = 0.70, step = 0.05
                ),
                shiny::splitLayout(
                    cellWidths = "110px",
                    shinyWidgets::prettySwitch(
                        inputId = ns(paste0(id, "_labels")),
                        label = "Show labels",
                        status = "success",
                        width = "50%",
                        value = TRUE,
                        fill = TRUE
                    ),
                    shinyWidgets::prettySwitch(
                        inputId = ns(paste0(id, "_scale")),
                        label = "Scale values to 0-1",
                        status = "success",
                        width = "50%",
                        value = scale_default,
                        fill = TRUE
                    ),
                ),
                if (contains_continuous) {
                    shiny::selectInput(
                        inputId = ns(paste0(id, "_colour_scheme")),
                        label = "Continuous colour scheme",
                        choices = c(""),
                        selected = NULL
                    )
                }
            )
        ),
        circle = TRUE,
        status = "success",
        size = "sm",
        icon = shiny::icon("cog")
    )
}



gear_download <- function(ns, id, label = "") {
    shinyWidgets::dropdownButton(
        label = "",
        icon = shiny::icon("download"),
        status = "success",
        size = "sm",
        shiny::em("Note: Use one of the following extensions: PDF, PNG, SVG."),
        shiny::textInput(ns(paste0("filename_", id)), "File name:", width = "80%", value = label),
        shiny::numericInput(ns(paste0("width_", id)), "Width (in):", 7, 3, 100, 0.1),
        shiny::numericInput(ns(paste0("height_", id)), "Height (in):", 7, 3, 100, 0.1),
        shiny::selectInput(ns(paste0("filetype_", id)), "Filetype", choices = c("PDF", "PNG", "SVG"), selected = "PDF", width = "80%"),
        shiny::downloadButton(ns(paste0("download_", id)), label = "Download Plot")
    )
}

pseudobulk_summary <- function(cell_info, mtd_value, summary_function = mean) {
    split_index <- split(seq_along(cell_info), mtd_value)
    pseudobulk_values <- rep(0, length(split_index))

    for (i in seq_along(split_index)) {
        pseudobulk_values[i] <- summary_function(cell_info[split_index[[i]]])
    }
    names(pseudobulk_values) <- names(split_index)
    
    return(pseudobulk_values)
}

req_gear_umap <- function(session, id = "settings") {
    for (input_val in paste0(id, "_", c("text_size", "axis_size", "legend_size", "pt_alpha", "pt_size", "pt_order", "trajectory_width", "labels", "scale", "colour_scheme"))) {
        shiny::req(!is.null(session$input[[input_val]]))
    }
}

req_gear_heatmap <- function(session, pseudobulk = TRUE) {
    used_vars <- c("scale_values", "axis_size", "legend_size", "colour_scheme", "cap_value")
    if (pseudobulk) {
        used_vars <- c(used_vars, "summarise_expr", "text_size", "point_size")
    } else {
        used_vars <- c(used_vars, "k_smooth" )
    }

    for (input_val in used_vars) {
        shiny::req(!is.null(session$input[[input_val]]))
    }
}

req_gear_download <- function(session, id) {
    for (input_val in paste0(c("filename_", "width_", "height_", "filetype_"), id)) {
        shiny::req(!is.null(session$input[[input_val]]))
    }
}

read_qs_format <- function(fl_path) {
    if (file.exists(fl_path)) {
        return(qs2::qs_read(fl_path, nthreads = 1))
    }

    fl_path <- sub("\\.qs2$", ".qs", fl_path)
    if (file.exists(fl_path)) {
        stop("Using deprecated `qs` format. Consider re-saving the file in `qs2` format for faster loading.")
    }

    stop(paste0("File not found: ", fl_path, "2"))
}

gene_name_transformation <- function(gene_name) {
    gene_name <- toupper(gene_name)
    gene_name <- gsub(" ", "", gene_name)
    gene_name <- gsub("-", "", gene_name)
    gene_name <- gsub("_", "", gene_name)
    gene_name <- gsub("\\.", "", gene_name)

    return(gene_name)
}

calculate_pseudotime_iqr_coverage <- function(mtd_df) {
    # this will calculate the sum of the area covered by the intervals defined by the q1 and q3 values
    mtd_df <- mtd_df[order(mtd_df$Q1, mtd_df$Q3), ]
    prev_start <- mtd_df$Q1[1]
    prev_end <- mtd_df$Q3[1]
    coverage <- prev_end - prev_start
    
    if (nrow(mtd_df) == 1) {
        return(coverage)
    }

    for (i in seq(from = 2, to = nrow(mtd_df))) {
        curr_start <- mtd_df$Q1[i]
        curr_end <- mtd_df$Q3[i]

        if (curr_start > prev_end) {
            coverage <- coverage + (curr_end - curr_start)
            prev_end <- curr_end
            prev_start <- curr_start
            next
        }

        if (curr_end > prev_end) {
            coverage <- coverage + (curr_end - prev_end)
            prev_end <- curr_end
        }

        prev_start <- curr_start
    }

    return(coverage)
}

calculate_umap_average_distance <- function(umap_df, selected_cells = NULL, centroid = TRUE) {
    if (is.null(selected_cells)) {
        selected_cells <- seq_len(nrow(umap_df))
    }

    if (length(selected_cells) < 2) {
        return(0)
    }

    umap_df <- as.matrix(umap_df[selected_cells, ])
    if (centroid) {
        # calculate centroid 
        centroid <- colMeans(umap_df, na.rm = TRUE)
        umap_df[, 1] <- umap_df[, 1] - centroid[1]
        umap_df[, 2] <- umap_df[, 2] - centroid[2]
        umap_df <- umap_df ^ 2
        umap_df <- rowSums(umap_df, na.rm = TRUE) ^ 0.5
        return(umap_df)
    }
    # pairwise distance
    dist_matrix <- as.matrix(stats::dist(umap_df))
    umap_df <- rowMeans(dist_matrix, na.rm = TRUE)

    return(umap_df)
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

# TODO you should work a bit more on improving this function
# idea- instead of using hard thresholds, try to define them based on the disitribution of iqrs you notice
# like the good thresh would be Q1 of the iqrs, and thresh bad would be Q3 + 1.5*IQR or something like that
detect_outlier <- function(modules_stats, cell_masks, psd_value, thresh_psd_good = NULL, thresh_psd_bad = NULL, umap_dist_threshold = NULL) {
    if (is.null(thresh_psd_bad)) {
        thresh_psd_bad <- stats::IQR(psd_value, na.rm = TRUE) * 0.85
    }
    modules_stats <- modules_stats %>% dplyr::arrange(
        .data$median_umap_distance,
        .data$iqr_pseudotime
    )
    cell_masks <- cell_masks[modules_stats$module, ]

    mad_z_scores_psd <- stats::setNames(mad_z_score(modules_stats$iqr_pseudotime), modules_stats$module)
    mad_z_scores_umap <- stats::setNames(mad_z_score(modules_stats$median_umap_distance), modules_stats$module)

    outlier_output <- rep("no", nrow(modules_stats))
    names(outlier_output) <- modules_stats$module

    ncells_per_module <- apply(cell_masks, 1, sum)

    outlier_output[abs(mad_z_scores_psd) > 3.5 | abs(mad_z_scores_umap) > 3.5] <- "redundant"
    outlier_output[abs(mad_z_scores_psd) > 3.5 & abs(mad_z_scores_umap) > 3.5] <- "yes"
    outlier_output[modules_stats$iqr_pseudotime >= thresh_psd_bad] <- "yes"
    outlier_output[ncells_per_module <= 10] <- "yes"

    if (!is.null(umap_dist_threshold)) {
        outlier_output[modules_stats$median_umap_distance >= umap_dist_threshold] <- "yes"
    }

    if (is.null(thresh_psd_good)) {
        thresh_psd_good <- stats::quantile(modules_stats$iqr_pseudotime[outlier_output == "no"], 0.25, na.rm = TRUE)
    }

    covered_mask <- rep(FALSE, ncol(cell_masks))

    coverage_evolution_df <- NULL
    for (i in seq_len(nrow(modules_stats))) {
        module_name <- modules_stats$module[i]
        is_already_labeled <- outlier_output[module_name] != "no"
        if (is_already_labeled) {
            next
        }

        current_mask <- cell_masks[module_name, ]
        ncells <- sum(current_mask)
        if (ncells < 10) {
            outlier_output[module_name] <- "redundant"
            next
        }
        temp_mask <- covered_mask | current_mask
        nunique <- sum(temp_mask) - sum(covered_mask)
        iqr <- modules_stats$iqr_pseudotime[i]
        eligible <- TRUE

        if (nunique / ncells <= 0.5) {
            eligible <- iqr < thresh_psd_good
        }

        potential_coverage_added <- nunique / length(covered_mask) + 1e-10
        potential_unique_percentage <- nunique / ncells + 1e-10
        percentage_iqr <- 1 / (iqr + 1e-10)
        potential_f1_score <- 2 / (1 / potential_unique_percentage + 1 / potential_coverage_added) #+ 1 / percentage_iqr)

        if (!eligible) {
            outlier_output[module_name] <- "redundant"
        } else {
            new_coverage <- data.frame(
                added_module = module_name,
                coverage_added = potential_coverage_added,
                percentage_unique = potential_unique_percentage,
                percentage_iqr = percentage_iqr,
                f1_score = potential_f1_score,
                module_iqr = iqr,
                median_pseudotime = modules_stats$median_pseudotime[i],
                median_umap_distance = modules_stats$median_umap_distance[i]
            )
            if (is.null(coverage_evolution_df)) {
                coverage_evolution_df <- new_coverage
            } else {
                coverage_evolution_df <- rbind(coverage_evolution_df, new_coverage)
            }
            covered_mask <- temp_mask
        }
    }

    # provide a second chance for the redundant modules
    # although most of their population overlap with other better modules, there's
    # a chance that the specific difference is not covered at all by the other modules
    # allow the module if the percentage of new cells is above 25%
    threshold_f1_score <- NULL
    if (!is.null(coverage_evolution_df)) {
        threshold_f1_score <- stats::median(coverage_evolution_df$f1_score, na.rm = TRUE)
    }

    for (i in seq_len(nrow(modules_stats))) {
        module_name <- modules_stats$module[i]
        if (outlier_output[module_name] != "redundant") {
            next
        }

        current_mask <- cell_masks[module_name, ]
        ncells <- sum(current_mask)
        if (ncells < 10) {
            next
        }
        temp_mask <- covered_mask | current_mask
        nunique <- sum(temp_mask) - sum(covered_mask)
        iqr <- modules_stats$iqr_pseudotime[i]

        potential_percentage_unique <- nunique / ncells + 1e-10
        potential_coverage_added <- nunique / length(covered_mask) + 1e-10
        potential_uncoverage_solved <- nunique / sum(!covered_mask) + 1e-10
        percentage_iqr <- 1 - 1 / (iqr + 1e-10)
        potential_f1_score <- 2 / (1 / potential_percentage_unique + 1 / potential_coverage_added)# + 1 / percentage_iqr)

        if (!is.null(threshold_f1_score) && potential_f1_score < threshold_f1_score) {
            next
        }
        
        outlier_output[module_name] <- "no"
        new_coverage <- data.frame(
            added_module = module_name,
            coverage_added = potential_coverage_added,
            percentage_unique = potential_percentage_unique,
            percentage_iqr = percentage_iqr,
            f1_score = potential_f1_score,
            module_iqr = iqr,
            median_pseudotime = modules_stats$median_pseudotime[i],
            median_umap_distance = modules_stats$median_umap_distance[i]
        )
        if (is.null(coverage_evolution_df)) {
            coverage_evolution_df <- new_coverage
        } else {
            coverage_evolution_df <- rbind(coverage_evolution_df, new_coverage)
        }
        covered_mask <- temp_mask
        threshold_f1_score <- stats::median(coverage_evolution_df$f1_score, na.rm = TRUE)
    }

    return(list(
        "outlier_output" = outlier_output,
        "coverage_evolution_df" = coverage_evolution_df
    ))
}

# ---------

build_module_masks <- function(module_summ, psd_mask, scale_threshold = 0, top_cells_percent = 100) {
    perc_cells <- 1 - top_cells_percent / 100
    stopifnot(
        perc_cells >= 0,
        perc_cells <= 1,
        scale_threshold >= 0,
        scale_threshold <= 1,
        !is.null(module_summ),
        !is.null(psd_mask)
    )

    mask_nrow <- length(module_summ)
    mask_ncol <- length(psd_mask)
    module_mask <- matrix(FALSE, nrow = mask_nrow, ncol = mask_ncol)
    rownames(module_mask) <- names(module_summ)

    for (i in seq_along(module_summ)) {
        module_mask[i, ] <- (module_summ[[i]] > scale_threshold) & psd_mask
        if (sum(module_mask[i, ]) < 5) {
            next
        }

        value_thresh <- stats::quantile(module_summ[[i]][module_mask[i, ]], perc_cells)
        mx_val <- max(module_summ[[i]])
        if (value_thresh == mx_val) {
            value_thresh <- 0.99 * mx_val
        }
        module_mask[i, ] <- (module_summ[[i]] >= value_thresh) & psd_mask
    }

    module_mask
}

compute_module_pairwise_tables <- function(module_mask, module_summ) {
    module_names <- rownames(module_mask)
    mask_nrow <- nrow(module_mask)

    module_intersect_cells <- matrix(0, nrow = mask_nrow, ncol = mask_nrow)
    rownames(module_intersect_cells) <- module_names
    colnames(module_intersect_cells) <- module_names

    if (mask_nrow > 1) {
        for (i in seq(from = 1, to = mask_nrow - 1)) {
            for (j in seq(from = i + 1, to = mask_nrow)) {
                ncommon_cells <- sum(module_mask[i, ] & module_mask[j, ])
                module_intersect_cells[i, j] <- ncommon_cells
                module_intersect_cells[j, i] <- ncommon_cells
            }
        }
    }

    for (i in seq_len(mask_nrow)) {
        unique_cells <- module_mask[i, ]
        for (j in seq_len(mask_nrow)) {
            if (i == j) {
                next
            }
            unique_cells <- unique_cells & (!module_mask[j, ])
        }
        module_intersect_cells[i, i] <- sum(unique_cells)
    }

    module_spearman_matrix <- matrix(NA, nrow = mask_nrow, ncol = mask_nrow)
    module_union_cells <- matrix(0, nrow = mask_nrow, ncol = mask_nrow)
    rownames(module_spearman_matrix) <- module_names
    colnames(module_spearman_matrix) <- module_names
    rownames(module_union_cells) <- module_names
    colnames(module_union_cells) <- module_names

    for (i in seq_along(module_summ)) {
        for (j in seq_along(module_summ)) {
            if (i >= j) {
                next
            }

            mask1 <- module_mask[i, ]
            mask2 <- module_mask[j, ]
            united_mask <- mask1 | mask2
            expr1 <- module_summ[[i]][united_mask]
            expr2 <- module_summ[[j]][united_mask]

            module_union_cells[i, j] <- sum(united_mask)
            module_union_cells[j, i] <- sum(united_mask)

            if (length(expr1) < 5 || length(expr2) < 5) {
                module_spearman_matrix[i, j] <- NA
                module_spearman_matrix[j, i] <- NA
                next
            }

            sp_val <- suppressWarnings(
                stats::cor(expr1, expr2, method = "spearman")
            )
            sp_val <- round(sp_val, 2)
            module_spearman_matrix[i, j] <- sp_val
            module_spearman_matrix[j, i] <- sp_val
        }
    }

    list(
        module_intersect_cells = module_intersect_cells,
        module_union_cells = module_union_cells,
        module_spearman_matrix = module_spearman_matrix
    )
}

get_module_stats <- function(module_summ, module_mask, psd_value, umap_df, centroid = TRUE) {
    modules_stats <- NULL

    for (i in seq_along(module_summ)) {
        if (sum(module_mask[i, ]) == 0) {
            next
        }
        temp_psd_val <- psd_value[module_mask[i, ]]
        if (all(is.na(temp_psd_val))) {
            next
        }

        temp_df <- data.frame(
            avg_summary = module_summ[[i]][module_mask[i, ]],
            psd_value = temp_psd_val,
            umap_distance = calculate_umap_average_distance(
                umap_df = umap_df,
                selected_cells = which(module_mask[i, ]),
                centroid = centroid
            )
        ) %>% dplyr::filter(!is.na(.data$psd_value))
        temp_df$module <- names(module_summ)[i]

        if (is.null(modules_stats)) {
            modules_stats <- temp_df
        } else {
            modules_stats <- rbind(modules_stats, temp_df)
        }
    }

    if (is.null(modules_stats)) {
        return(NULL)
    }

    modules_stats$module <- factor(modules_stats$module, levels = names(module_summ))
    return(modules_stats)
}

summarise_module_stats <- function(modules_stats, gene_modules) {
    modules_stats$module <- as.character(modules_stats$module)
    summary_stats <- modules_stats %>%
        dplyr::group_by(.data$module) %>%
        dplyr::summarise(
            n_cells = length(.data$avg_summary),
            avg_summary = round(mean(.data$avg_summary, na.rm = TRUE), 3),
            median_pseudotime = round(stats::median(.data$psd_value, na.rm = TRUE), 3),
            iqr_pseudotime = round(stats::IQR(.data$psd_value, na.rm = TRUE), 3),
            median_umap_distance = round(stats::median(.data$umap_distance, na.rm = TRUE), 3)
        )
    summary_stats$n_genes <- sapply(summary_stats$module, function(mod) length(gene_modules[[mod]]))
    current_cols <- colnames(summary_stats)
    new_cols <- c(current_cols[1], "n_genes", current_cols[2:(length(current_cols) - 1)])
    summary_stats <- as.data.frame(summary_stats)[, new_cols]
    rownames(summary_stats) <- summary_stats$module
    return(summary_stats)
}

annotate_module_outliers <- function(modules_stats_summary, module_mask, psd_value, umap_dist_threshold = NULL) {
    if (is.null(modules_stats_summary) || nrow(modules_stats_summary) == 0) {
        return(modules_stats_summary)
    }

    psd_span <- stats::quantile(psd_value, 0.95, na.rm = TRUE) - stats::quantile(psd_value, 0.05, na.rm = TRUE)
    outlier_result <- detect_outlier(
        modules_stats = modules_stats_summary,
        cell_masks = module_mask,
        psd_value = psd_value,
        thresh_psd_good = NULL, #psd_span / 10,
        thresh_psd_bad = NULL, #psd_span / 3
        umap_dist_threshold = umap_dist_threshold
    )

    modules_stats_summary$is_outlier <- outlier_result$outlier_output[as.character(modules_stats_summary$module)]

    modules_stats_summary %>%
        dplyr::arrange(.data$median_pseudotime, .data$iqr_pseudotime, .data$median_umap_distance)
}
