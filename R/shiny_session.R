env <- new.env(parent = emptyenv())

###### UI ######
#' UI - Global Settings
#'
#' @description Creates the UI elements with the global settings for the UI
#' components inside Starlng Shiny application.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
ui_global_setttings <- function() {
    shiny::tagList(
        shiny::tags$style(
            shiny::HTML(
                ".custom-empty {
                    margin-bottom: 40px;
                }",
                ".shiny-tab-input {
                    font-size: 20px;
                }",
                ".shiny-split-layout > div {
                    overflow: visible;
                }",
                ".btn-danger {
                    font-size: 20px;
                }",
                ".empty-space {
                    margin-bottom: 200px;
                }"
            )
        ),
        shiny::tags$head(shiny::tags$script(shiny::HTML('
            var dimension = [0, 0];
            var resizeId;
            $(document).on("shiny:connected", function(e) {
                dimension[0] = window.innerWidth - 20;
                dimension[1] = window.innerHeight - 30;
                Shiny.onInputChange("dimension", dimension);
            });

            function transferWindowSize() {
                console.log(dimension);
                console.log(window.innerHeight);

                let dif_width = Math.abs(window.innerWidth - 20 - dimension[0]);
                let dif_height = Math.abs(window.innerHeight - 30 - dimension[1]);
                console.log(dif_height);

                if (dif_width >= 200 || dif_height >= 200) {
                console.log("Changed")
                dimension[0] = window.innerWidth - 20;
                dimension[1] = window.innerHeight - 30;
                Shiny.onInputChange("dimension", dimension);
                }
            }

            $(window).resize(function() {
                clearTimeout(resizeId);
                resizeId = setTimeout(transferWindowSize, 500);
            });
        '))),
        shiny::div(class = "custom-empty")
    )
}

###### SERVER ######
#' Server - Prepare Session
#'
#' @description Creates the backend interface that initializes the variables
#' and the environment used inside the Starlng Shiny application.
#'
#' @param reactive_dim A reactive variable that returns the dimensions of the
#' browser window. This reactive variable is used to update the dimensions of
#' the plots inside the Shiny app.
#' @param height_ratio Variable indicating the ratio between the height of the
#' plots and the height of the browser window.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
prepare_session <- function(reactive_dim, height_ratio = 0.7) {
    # NOTE can we put this as a setting for the user to set in the app?
    assign("height_ratio", height_ratio, envir = env)
    assign("window_dim", reactive_dim, envir = env)

    mtd_path <- file.path("objects", "metadata.qs")
    if (!file.exists(mtd_path)) {
        stop("`objects/metadata.qs` file not found")
    }
    mtd_df <- qs::qread(mtd_path, nthreads = 1)
    ncols <- ncol(mtd_df)
    assign("umap_df", mtd_df[, (ncols - 1):ncols], envir = env)
    assign("mtd_df", mtd_df[, seq_len(ncols - 2)], envir = env)

    discrete_mask <- sapply(colnames(mtd_df), function(x_col) {
        inherits(mtd_df[[x_col]], c("factor", "character")) 
    })
    discrete_mtd <- lapply(colnames(mtd_df)[discrete_mask], function(x_col) {
        stringr::str_sort(unique(mtd_df[[x_col]]), numeric = TRUE)
    })
    names(discrete_mtd) <- colnames(mtd_df)[discrete_mask]
    assign("discrete_mtd", discrete_mtd, envir = env)

    recommended_psd_path <- file.path("objects", "recommended_pseudotime.qs")
    if (!file.exists(recommended_psd_path)) {
        stop("`objects/recommended_pseudotime.qs` file not found")
    }
    assign("recommended_psd", qs::qread(recommended_psd_path, nthreads = 1), envir = env)

    mon_obj_path <- file.path("objects", "diet_monocle_object.qs")
    if (!file.exists(mon_obj_path)) {
        stop("`objects/diet_monocle_object.qs` file not found")
    }
    assign("mon_obj", qs::qread(mon_obj_path, nthreads = 1), envir = env)

    colors_path <- file.path("objects", "colours.qs")
    if (!file.exists(colors_path)) {
        stop("`objects/colours.qs` file not found")
    }
    assign("color_options", qs::qread(colors_path, nthreads = 1), envir = env)

    traj_path <- file.path("objects", "trajectory_ggplot.qs")
    if (!file.exists(traj_path)) {
        stop("`objects/trajectory_gplot.qs` file not found")
    }
    assign("trajectory_gplot", qs::qread(traj_path, nthreads = 1), envir = env)

    moran_path <- file.path("objects", "genes_info.csv")
    if (!file.exists(moran_path)) {
        stop("`objects/genes_info.csv` file not found")
    }
    assign("moran_df", read.csv(moran_path, row.names = 1), envir = env)

    stable_modules_path <- file.path("objects", "stable_modules.csv")
    if (file.exists(stable_modules_path)) {
        pre_stable_modules <- read.csv(stable_modules_path, row.names = 1)
        for (i in colnames(pre_stable_modules)) {
            pre_stable_modules[, i] <- factor(pre_stable_modules[, i])
        }
        assign("preloaded_stable_modules", pre_stable_modules, envir = env)
        # assign("preloaded_options", shiny::reactiveVal(colnames(pre_stable_modules)), envir = env)
    }
    assign("chosen_modules", shiny::reactiveVal(NULL), envir = env)

    expr_path <- file.path("objects", "expression.h5")
    assign("genes", rhdf5::h5read(expr_path, "genes"), envir = env)
    index <- seq_along(env$genes)
    names(index) <- env$genes
    assign("genes", index, envir = env)
    assign("cells", rhdf5::h5read(expr_path, "cells"), envir = env)

    assign("pseudotime_changes", shiny::reactiveVal(0), envir = env)
    assign("organism", "hsapiens", envir = env)
}

#' Server - Gear Width
#'
#' @description Creates the backend interface that dynamically updates the
#' width of the gear settings window inside the Starlng Shiny application.
#'
#' @param session The shiny session object used to update the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
update_gears_width <- function(session) {
    shiny::observe({
        win_dims <- env$window_dim()
        shiny::req(win_dims)

        shinyjs::runjs(paste0(
            "$('.dropdown-menu .shiny-split-layout').css('width', '", win_dims[1] / 2.2, "px');"
        ))
    })
}

#' Server - Tabs update
#'
#' @description Creates the backend interface that dynamically shows the gene
#' module tabs inside the Starlng Shiny application after running the gene
#' clustering.
#'
#' @param session The shiny session object used to update the UI elements.
#'
#' @note This function is a shiny module function and should be used
#' in the context of the app created using the `starlng_write_app` function.
#'
#' @export
update_tabs <- function(session) {
    shiny::observe({
        shinyjs::disable("module_enrichment")
        for (target_tab in c("module_enrichment", "module_umap", "gene_heatmaps")) {
            shiny::hideTab(
                inputId = "nav",
                target = target_tab,
                session = session
            )
        }
        shiny::req(env$chosen_modules())
        for (target_tab in c("module_enrichment", "module_umap", "gene_heatmaps")) {
            shiny::showTab(
                inputId = "nav",
                target = target_tab,
                session = session
            )
        }
        shiny::updateNavbarPage(
            session = session,
            inputId = "nav",
            selected = "module_umap"
        )
    })
}
