save_filetypes <- list(
    "PDF" = pdf,
    "PNG" = function(filename, width, height) {
        ragg::agg_png(filename, width, height, units = "in")
    },
    "SVG" = svg
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
        used_vars <- c(used_vars, "summarise_expr", "text_size")
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
        warning("Using deprecated `qs` format.")
        return(qs::qread(fl_path, nthreads = 1))
    }

    stop(paste0("File not found: ", fl_path, "2"))
}
