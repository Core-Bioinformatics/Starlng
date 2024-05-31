library(shiny)
library(monocle3)
library(Seurat)
library(ClustAssess)

# example of calling
# write_object.data_frame(ids_df = read.csv("/servers/sutherland-scratch/andi/projects/0_2304_Floris_snrna/metadata/sample_combinations.csv"), seurat_ca_corresp = read.csv("/servers/sutherland-scratch/andi/projects/0_2304_Floris_snrna/metadata/seurat_paths_to_clustassess.csv"), target_app_dir = "/servers/sutherland-scratch/andi/projects/0_2304_Floris_snrna/output/pseudotime_apps/first_part_analysis")
write_object.data_frame <- function(ids_df,
                                    seurat_ca_corresp,
                                    target_app_dir,
                                    prefix = "",
                                    suffix = "",
                                    use_closed_loops = FALSE,
                                    use_partitions = FALSE,
                                    learn_graph_controls = NULL,
                                    nodes_per_log10_cells = 30) {
    common_seurat_dir <- find_common_directory(seurat_ca_corresp$seurat)
    print(common_seurat_dir)
    for (i in seq_along(ids_df[, 1])) {
        id <- ids_df$Id[i]

        if (is.na(ids_df$k[i])) {
            next
        }

        associated_index <- find_corresp_index(id, seurat_ca_corresp, prefix, suffix)
        print(associated_index)

        if (associated_index == -1) {
            next
        }

        print(glue::glue("Creating app for {id}"))
        target_folder <- target_app_dir
        seurat_path <- dirname(seurat_ca_corresp$seurat[associated_index])
        to_add_folders <- c(id, ifelse(use_closed_loops, "loop", "no_loop"), ids_df$Ftype[i], ids_df$Fsize[i])

        while (basename(seurat_path) != common_seurat_dir) {
            to_add_folders <- c(basename(seurat_path), to_add_folders)
            seurat_path <- dirname(seurat_path)
        }

        to_add_folders <- as.list(to_add_folders)

        for (j in seq_along(to_add_folders)) {
            target_folder <- file.path(target_folder, to_add_folders[[j]])
        }

        dir.create(target_folder, showWarnings = FALSE, recursive = TRUE)

        seurat_path <- seurat_ca_corresp$seurat[associated_index]
        ca_path <- seurat_ca_corresp$clustassess[associated_index]
        write_object.Seurat(
            seurat_object = readRDS(seurat_path),
            clustassess_object = readRDS(ca_path),
            output_dir = target_folder,
            app_name = id,
            stable_config = list(
                ftype = ids_df$Ftype[i],
                fsize = ids_df$Fsize[i],
                clmethod = ids_df$Cl_method[i],
                k = ids_df$k[i]
            ),
            use_closed_loops = use_closed_loops,
            use_partitions = use_partitions,
            learn_graph_controls = learn_graph_controls,
            nodes_per_log10_cells = nodes_per_log10_cells 
        )
    }
}


write_object.Seurat <- function(seurat_object,
                                clustassess_object,
                                output_dir,
                                app_name,
                                assay_used = "SCT",
                                stable_config = list(
                                    ftype = "",
                                    fsize = "",
                                    clmethod = "",
                                    k = ""
                                ),
                                use_closed_loops = FALSE,
                                use_partitions = FALSE,
                                learn_graph_controls = NULL,
                                nodes_per_log10_cells = 30) {
    write_object.default(
        expression_matrix = seurat_object@assays[[assay_used]]@data,
        meta_data = seurat_object@meta.data,
        clustassess_object = clustassess_object,
        output_dir = output_dir,
        app_name = app_name,
        stable_config = stable_config,
        use_closed_loops = use_closed_loops,
        use_partitions = use_partitions,
        learn_graph_controls = learn_graph_controls,
        nodes_per_log10_cells = nodes_per_log10_cells
    )
}

# write_object.default <- function(output_dir, id, paths_metadata) {
write_object.default <- function(expression_matrix,
                                 meta_data,
                                 clustassess_object,
                                 output_dir,
                                 app_name,
                                 stable_config = list(
                                     ftype = "",
                                     fsize = "",
                                     clmethod = "",
                                     k = ""
                                 ),
                                 use_closed_loops = FALSE,
                                 use_partitions = FALSE,
                                 learn_graph_controls = NULL,
                                 nodes_per_log10_cells = 30) {
    mon_obj <- create_monocle_object(
        normalized_expression_matrix = expression_matrix,
        clustassess_object = clustassess_object,
        metadata = meta_data,
        stable_feature_type = stable_config$ftype,
        stable_feature_set_size = stable_config$fsize,
        stable_clustering_method = stable_config$clmethod,
        stable_n_clusters = stable_config$k,
        cell_names = colnames(expression_matrix),
        gene_names = rownames(expression_matrix),
        use_all_genes = TRUE
    )

    mon_obj <- cluster_cells(mon_obj)
    mon_obj@clusters@listData$UMAP$clusters <- setNames(mon_obj@colData[[paste0("stable_", stable_config$k[1], "_clusters")]], colnames(mon_obj))
    mon_obj@clusters@listData$UMAP$partitions <- setNames(mon_obj@colData[[paste0("stable_", stable_config$k[1], "_clusters")]], colnames(mon_obj))

    if (is.list(learn_graph_controls)) {
        learn_graph_controls[["ncenter"]] <- cal_ncenter(as.integer(stable_config$k), ncol(mon_obj), nodes_per_log10_cells = nodes_per_log10_cells) 
    }

    mon_obj <- learn_graph(
        cds = mon_obj,
        close_loop = use_closed_loops,
        use_partition = use_partitions,
        learn_graph_control = learn_graph_controls
    )

    saveRDS(mon_obj, file.path(output_dir, glue::glue("monocle_object.rds")))

    write_app_file(output_dir, app_name)
    styler::style_file(file.path(output_dir, "app.R"), indent_by = 4)
}

write_app_file <- function(output_dir, id) {
    file_content <- paste0("
    library(monocle3)
    library(shiny)
    library(shinyWidgets)
    library(shinyjs)
    library(ggplot2)
    library(ggplotify)
    library(dplyr)
    library(tidyr)


    gear_trajectory <- function(id) {
        dropdownButton(
            sliderInput(
                inputId = paste0(\"cellsize\", id),
                label = \"Point size\",
                min = 0.1,
                max = 5,
                step = 0.05,
                value = 0.1
            ),
            sliderInput(
                inputId = paste0(\"alpha\", id),
                label = \"Color alpha\",
                min = 0.1,
                max = 1,
                step = 0.05,
                value = 0.5
            ),
            sliderInput(
                inputId = paste0(\"edgesize\", id),
                label = \"Size of edges\",
                min = 0.1,
                max = 5,
                step = 0.05,
                value = 0.75
            ),
            circle = TRUE,
            status = \"success\",
            size = \"sm\",
            icon = shiny::icon(\"cog\")
        )
    }

    gear_umaps <- function(id, narrow_options = TRUE) {
        shinyWidgets::dropdownButton(
            shiny::tagList(
                shiny::sliderInput(
                    inputId = paste0(id, \"_pt_size\"),
                    label = \"Point size\",
                    min = 0.05, max = 5.00, value = 0.10, step = 0.05
                ),
                if (narrow_options) {
                    shinyWidgets::prettySwitch(
                        inputId = paste0(id, \"_narrow\"),
                        label = \"Narrow down the area for points of interest\",
                        value = TRUE,
                        status = \"success\",
                        fill = TRUE
                    )
                }
            ),
            circle = TRUE,
            status = \"success\",
            size = \"sm\",
            icon = shiny::icon(\"cog\")
        )
    }

    gear_download <- function(id, label = \"\") {
        shinyWidgets::dropdownButton(
            label = \"\",
            icon = shiny::icon(\"download\"),
            status = \"success\",
            size = \"sm\",
            shiny::em(\"Note: Use one of the following extensions: PDF, PNG, SVG.\"),
            shiny::textInput(ns(paste0(\"filename_\", id)), \"File name:\", width = \"80%\", value = label),
            shiny::numericInput(ns(paste0(\"width_\", id)), \"Width (in):\", 7, 3, 100, 0.1),
            shiny::numericInput(ns(paste0(\"height_\", id)), \"Height (in):\", 7, 3, 100, 0.1),
            shiny::selectInput(ns(paste0(\"filetype_\", id)), \"Filetype\", choices = c(\"PDF\", \"PNG\", \"SVG\"), selected = \"PDF\", width = \"80%\"),
            shiny::downloadButton(ns(paste0(\"download_\", id)), label = \"Download Plot\")
        )
    }


    ui <- fluidPage(
        tags$style(
            HTML(
                \"#tabset_id {
                        position: fixed;
                        width: 100%;
                        background-color: white;
                        top: 0;
                        z-index: 100;
                        font-size: 25px;
                        # margin-top: 72px;
                        }\",
                \"#graph_clust-show_config {
                        z-index: 900;
                        position: fixed;
                        top: 10px;
                        font-size: 15px;
                        right: 25px;
                    }\",
                \".shiny-split-layout > div {
                        overflow: visible;
                    }\",
            )
        ),
        tags$head(tags$script(HTML('
            var dimension = [0, 0];
            var resizeId;
            $(document).on(\"shiny:connected\", function(e) {
                dimension[0] = window.innerWidth - 20;
                dimension[1] = window.innerHeight - 30;
                Shiny.onInputChange(\"dimension\", dimension);
            });

            function transferWindowSize() {
                console.log(dimension);
                console.log(window.innerHeight);

                let dif_width = Math.abs(window.innerWidth - 20 - dimension[0]);
                let dif_height = Math.abs(window.innerHeight - 30 - dimension[1]);
                console.log(dif_height);

                if (dif_width >= 200 || dif_height >= 200) {
                console.log(\"Changed\")
                dimension[0] = window.innerWidth - 20;
                dimension[1] = window.innerHeight - 30;
                Shiny.onInputChange(\"dimension\", dimension);
                }
            }

            $(window).resize(function() {
                clearTimeout(resizeId);
                resizeId = setTimeout(transferWindowSize, 500);
            });
        '))),
        useShinyjs(),
        h1(\"Pseudotime using Monocle3 - Choosing the start and ending points - ", id, "\"),
        h3(\"Trajectory graph\"),
        splitLayout(
            cellWidths = c(\"40px\", \"350px\"),
            gear_trajectory(\"1\"),
            selectInput(
                inputId = \"metadata1\",
                label = \"Select metadata\",
                choices = c(\"\"),
                selected = NULL
            )
        ),
        plotOutput(\"graph_metadata\", height = \"auto\"),
        h3(\"Input / Output cells\"),
        splitLayout(
            wellPanel(
                selectizeInput(
                    inputId = \"genes2\",
                    label = \"Genes\",
                    choices = NULL,
                    multiple = TRUE
                ),
                splitLayout(
                    shiny::numericInput(
                        inputId = \"expr_threshold2\",
                        label = \"Gene expression threshold\",
                        min = 0, max = 10, value = 0, step = 0.01,
                        width = \"95%\"
                    ),
                    shiny::numericInput(
                        inputId = \"relaxation2\",
                        label = \"#genes not expressed\",
                        min = 0, max = 10, value = 0, step = 1,
                        width = \"95%\"
                    )
                ),
                gear_umaps(\"umap2\", TRUE),
                plotOutput(\"start_cells\", height = \"auto\"),
                shiny::actionButton(
                    inputId = \"fix_start_cells\",
                    label = \"Fix the starting cells!\",
                    style = \"font-size:20px;\",
                    class = \"btn-danger\"
                )
            ),
            wellPanel(
                selectizeInput(
                    inputId = \"genes3\",
                    label = \"Genes\",
                    choices = NULL,
                    multiple = TRUE
                ),
                splitLayout(
                    shiny::numericInput(
                        inputId = \"expr_threshold3\",
                        label = \"Gene expression threshold\",
                        min = 0, max = 10, value = 0, step = 0.01,
                        width = \"95%\"
                    ),
                    shiny::numericInput(
                        inputId = \"relaxation3\",
                        label = \"#genes not expressed\",
                        min = 0, max = 10, value = 0, step = 1,
                        width = \"95%\"
                    )
                ),
                gear_umaps(\"umap3\", TRUE),
                plotOutput(\"end_cells\", height = \"auto\"),
                shiny::actionButton(
                    inputId = \"fix_end_cells\",
                    label = \"Fix the end cells!\",
                    style = \"font-size:20px;\",
                    class = \"btn-danger\"
                )
            )
        ),
        h3(\"Pseudotime - all data\"),
        shiny::actionButton(
            inputId = \"calculate_pseudotime\",
            label = \"Calculate pseudotime!\",
            style = \"font-size:20px\",
            class = \"btn-danger\"
        ),
        splitLayout(
            verticalLayout(
                gear_trajectory(\"start\"),
                plotOutput(\"pseudotime_start\", height = \"auto\"),
            ),
            verticalLayout(
                gear_trajectory(\"end\"),
                plotOutput(\"pseudotime_end\", height = \"auto\")
            )
        ),
        verticalLayout(
            gear_trajectory(\"zoomed_in\"),
            plotOutput(\"pseudotime_zoomed_in\", height = \"auto\"),
        ),
        shiny::p(\"To recalculate the trajectory graph on the area between starting and end genes in a new shiny app, provide the following info to the bioinformatician (along with a name for the chosen trajectory).\"),
        shiny::verbatimTextOutput(\"decision\")
    )

    server <- function(input, output, session) {
        shinyjs::hide(\"calculate_pseudotime\")
        mon_obj <- readRDS(\"monocle_object.rds\")
        gene_list <- rownames(mon_obj)
        umap_df <- data.frame(SingleCellExperiment::reducedDims(mon_obj)$UMAP)
        height_ratio <- 0.65

        shiny::updateSelectInput(
            session = session,
            inputId = \"metadata1\",
            choices = colnames(mon_obj@colData)
        )

        shiny::updateSelectizeInput(
            session = session,
            inputId = \"genes2\",
            server = TRUE,
            choices = gene_list,
            selected = gene_list[1],
            options = list(
                maxOptions = 7,
                create = TRUE,
                persist = TRUE
            )
        )

        shiny::updateSelectizeInput(
            session = session,
            inputId = \"genes3\",
            server = TRUE,
            choices = gene_list,
            selected = gene_list[1],
            options = list(
                maxOptions = 7,
                create = TRUE,
                persist = TRUE
            )
        )

        plt_height <- shiny::reactive(
            floor(height_ratio * input$dimension[2])
        )

        plt_width <- shiny::reactive(
            input$dimension[1]
        )

        graph_plot_data <- shiny::reactive({
            req(input$metadata1 != \"\")
            shiny::req(plt_width() > 0)
            shiny::req(plt_height() > 0)
            p <- plot_cells(mon_obj,
                trajectory_graph_color = \"black\",
                alpha = input$alpha1,
                cell_size = input$cellsize1,
                trajectory_graph_segment_size = input$edgesize1,
                color_cells_by = input$metadata1,
                label_cell_groups = FALSE,
                # label_leaves=FALSE,
                # label_branch_points=FALSE,
                # cell_stroke = 0,
                graph_label_size = 0
            )
            p_built <- ggplot_build(p)
            p_built$data[[1]]$stroke <- NA

            p_built
        })
        output$graph_metadata <- renderPlot(
            {
                as.ggplot(ggplot_gtable(graph_plot_data()))
            },
            width = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width() / 2, plt_height()), 400)
            },
            height = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width() / 2, plt_height()), 400)
            }
        )

        start_mask <- shiny::reactive({
            req(input$genes2 != \"\")
            req(input$expr_threshold2)
            if (length(input$genes2) == 0) {
                color_info <- rep(FALSE, nrow(umap_df))
            } else {
                if (length(input$genes2) > 1) {
                    color_info <- colSums2(mon_obj@assays@data$counts[input$genes2, ] > input$expr_threshold2) >= (length(input$genes2) - input$relaxation2)
                } else {
                    color_info <- mon_obj@assays@data$counts[input$genes2, ] > input$expr_threshold2
                }
            }

            if (input$\"umap2_narrow\") {
                umap_emb <- umap_df[color_info, ]
                if (nrow(umap_emb) == 0) {
                    break
                }
                gmedian <- Gmedian::Gmedian(umap_emb)
                overall_distance <- ((umap_emb[, 1] - gmedian[1, 1])^2 + (umap_emb[, 2] - gmedian[1, 2])^2)^0.5
                dist_stats <- fivenum(overall_distance)
                color_info[seq_len(nrow(umap_df))[color_info][overall_distance > dist_stats[3]]] <- FALSE
            }


            color_info
        })


        end_mask <- shiny::reactive({
            req(input$genes3 != \"\")
            req(input$expr_threshold3)
            if (length(input$genes3) == 0) {
                color_info <- rep(FALSE, nrow(umap_df))
            } else {
                if (length(input$genes3) > 1) {
                    color_info <- colSums2(mon_obj@assays@data$counts[input$genes3, ] > input$expr_threshold3) >= (length(input$genes3) - input$relaxation3)
                } else {
                    color_info <- mon_obj@assays@data$counts[input$genes3, ] > input$expr_threshold3
                }
            }

            if (input$\"umap3_narrow\") {
                umap_emb <- umap_df[color_info, ]
                if (nrow(umap_emb) > 0) {
                    gmedian <- Gmedian::Gmedian(umap_emb)
                    overall_distance <- ((umap_emb[, 1] - gmedian[1, 1])^2 + (umap_emb[, 2] - gmedian[1, 2])^2)^0.5
                    dist_stats <- fivenum(overall_distance)
                    color_info[seq_len(nrow(umap_df))[color_info][overall_distance > dist_stats[3]]] <- FALSE
                }
            }


            color_info
        })

        output$start_cells <- renderPlot(
            {
                req(start_mask())
                umap_df$color_val <- start_mask()

                ggplot(rbind(umap_df[!umap_df$color_val, ], umap_df[umap_df$color_val, ]), aes(x = UMAP_1, y = UMAP_2, color = color_val)) +
                    geom_point(size = input$\"umap2_pt_size\") +
                    theme_bw() +
                    scale_color_manual(values = c(\"lightgray\", \"red\")) +
                    theme(legend.position = \"None\")
            },
            width = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width() / 2.1, plt_height()), 400)
            },
            height = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width() / 2.1, plt_height()), 400)
            }
        )

        output$end_cells <- renderPlot(
            {
                shiny::req(end_mask())
                umap_df$color_val <- end_mask()

                ggplot(rbind(umap_df[!umap_df$color_val, ], umap_df[umap_df$color_val, ]), aes(x = UMAP_1, y = UMAP_2, color = color_val)) +
                    geom_point(size = input$\"umap2_pt_size\") +
                    theme_bw() +
                    scale_color_manual(values = c(\"lightgray\", \"red\")) +
                    theme(legend.position = \"None\")
            },
            width = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width() / 2.1, plt_height()), 400)
            },
            height = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width() / 2.1, plt_height()), 400)
            }
        )


        start_cells <- shiny::reactive(
            colnames(mon_obj)[seq_len(nrow(umap_df))[start_mask()]]
        ) %>% shiny::bindEvent(input$\"fix_start_cells\")

        end_cells <- shiny::reactive(
            colnames(mon_obj)[seq_len(nrow(umap_df))[end_mask()]]
        ) %>% shiny::bindEvent(input$\"fix_end_cells\")

        shiny::observe({
            shiny::req(start_cells())
            shiny::req(end_cells())

            shinyjs::show(\"calculate_pseudotime\")
        })

        shiny::observe({
        }) %>% shiny::bindEvent(input$calculate_pseudotime)


        start_pseudotime <- shiny::reactiveVal(NULL)
        end_pseudotime <- shiny::reactiveVal(NULL)
        pseudotime_between_cells <- shiny::reactiveVal(NULL)

        shiny::observe({
            shiny::req(start_cells(), end_cells())
            shinyjs::disable(\"calculate_pseudotime\")

            mon_obj <- order_cells(mon_obj, root_cells = end_cells())
            end_pseudotime(pseudotime(mon_obj))

            mon_obj <- order_cells(mon_obj, root_cells = start_cells())
            start_pseudotime(pseudotime(mon_obj))

            nodes_start <- principal_graph_aux(mon_obj)$UMAP$pr_graph_cell_proj_closest_vertex[, 1][start_cells()]
            nodes_end <- principal_graph_aux(mon_obj)$UMAP$pr_graph_cell_proj_closest_vertex[, 1][end_cells()]

            cells_between_start <- choose_graph_segments(
                mon_obj,
                starting_pr_node = nodes_start[1],
                ending_pr_nodes = nodes_end,
                clear_cds = FALSE
            )

            cells_between_end <- choose_graph_segments(
                mon_obj,
                starting_pr_node = nodes_end[1],
                ending_pr_nodes = nodes_start,
                clear_cds = FALSE
            )

            all_vertices_between <- union(
                igraph::vertex_attr(cells_between_end@principal_graph$UMAP, \"name\"),
                igraph::vertex_attr(cells_between_start@principal_graph$UMAP, \"name\")
            )

            all_cells_between <- union(
                colnames(cells_between_start),
                colnames(cells_between_end)
            )


            psdt <- pseudotime(mon_obj)
            psdt[!(colnames(mon_obj) %in% all_cells_between)] <- NA
            pseudotime_between_cells(psdt)

            shinyjs::enable(\"calculate_pseudotime\")
        }) %>% shiny::bindEvent(input$calculate_pseudotime)

        output$pseudotime_start <- shiny::renderPlot(
            {
                shiny::req(start_pseudotime())
                mon_obj@principal_graph_aux$UMAP$pseudotime <- start_pseudotime()
                p <- plot_cells(mon_obj,
                    trajectory_graph_color = \"black\",
                    alpha = input$alphastart,
                    cell_size = input$cellsizestart,
                    trajectory_graph_segment_size = input$edgesizestart,
                    color_cells_by = \"pseudotime\",
                    graph_label_size = 0
                ) + scale_color_viridis_c(\"pseudotime\")
                p_built <- ggplot_build(p)
                p_built$data[[1]]$stroke <- NA
                as.ggplot(ggplot_gtable(p_built)) + ggtitle(\"Pseudotime calculated using starting cells\")
            },
            width = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width() / 2.1, plt_height()), 400)
            },
            height = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width() / 2.1, plt_height()), 400)
            }
        )

        output$pseudotime_end <- shiny::renderPlot(
            {
                shiny::req(end_pseudotime())
                mon_obj@principal_graph_aux$UMAP$pseudotime <- end_pseudotime()
                p <- plot_cells(mon_obj,
                    trajectory_graph_color = \"black\",
                    alpha = input$alphaend,
                    cell_size = input$cellsizeend,
                    trajectory_graph_segment_size = input$edgesizeend,
                    color_cells_by = \"pseudotime\",
                    graph_label_size = 0
                ) + scale_color_viridis_c(\"pseudotime\")
                p_built <- ggplot_build(p)
                p_built$data[[1]]$stroke <- NA
                as.ggplot(ggplot_gtable(p_built)) + ggtitle(\"Pseudotime calculated using end cells\")
            },
            width = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width() / 2.1, plt_height()), 400)
            },
            height = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width() / 2.1, plt_height()), 400)
            }
        )

        output$pseudotime_zoomed_in <- shiny::renderPlot(
            {
                shiny::req(pseudotime_between_cells())
                mon_obj@principal_graph_aux$UMAP$pseudotime <- pseudotime_between_cells()
                p <- plot_cells(mon_obj,
                    trajectory_graph_color = \"black\",
                    alpha = input$alphazoomed_in,
                    cell_size = input$cellsizezoomed_in,
                    trajectory_graph_segment_size = input$edgesizezoomed_in,
                    color_cells_by = \"pseudotime\",
                    graph_label_size = 0
                ) + scale_color_viridis_c(\"pseudotime\")
                p_built <- ggplot_build(p)
                p_built$data[[1]]$stroke <- NA
                as.ggplot(ggplot_gtable(p_built)) + ggtitle(\"Pseudotime from start to end cells\")
            },
            width = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width(), plt_height()), 400)
            },
            height = function() {
                shiny::req(plt_width(), plt_height())
                max(min(plt_width(), plt_height()), 400)
            }
        )

        output$decision <- renderPrint({
            paste(
                \"", id, "\",
                paste(input$genes2, collapse = ';'),
                input$expr_threshold2,
                input$relaxation2,
                paste(input$genes3, collapse = ';'),
                input$expr_threshold3,
                input$relaxation3,
                sep = \",\"
            )
        })
    }

    shinyApp(ui, server)
    ")
    write(file_content, file.path(output_dir, "app.R"))
}



# seurat_clustassess_csv is a dataframe read from a csv that contains two columns:
# seurat - path to the Seurat object
# clustassess - path to the ClustAssess object
find_corresp_index <- function(target_id, seurat_clustassess_csv, prefix = "", suffix = "") {
    for (i in seq_len(nrow(seurat_clustassess_csv))) {
        basename_seurat <- basename(seurat_clustassess_csv$seurat[i])

        if (basename_seurat == glue::glue("{prefix}{target_id}{suffix}.rds")) {
            return(i)
        }

        if (basename_seurat == glue::glue("{target_id}.rds")) {
            return(i)
        }
    }

    return(-1)
}


cal_ncenter <- function(num_cell_communities, ncells,
                        nodes_per_log10_cells=15) {
  round(num_cell_communities * nodes_per_log10_cells * log10(ncells))
}


# to concatenate them back, `do.call(file.path, as.list(dir_list))`
get_file_hierarchy <- function(file_path) {
    dir_list <- c()

    while (file_path != "/") {
        file_path <- dirname(file_path)

        dir_list <- c(basename(file_path), dir_list)
    }

    return(dir_list)
}

find_common_directory <- function(path_list, return_full = FALSE) {
    base_comp <- get_file_hierarchy(path_list[1])

    for (i in seq_along(path_list)) {
        if (i == 1) {
            next
        }

        current_path <- get_file_hierarchy(path_list[i])
        min_length <- min(length(base_comp), length(current_path))

        for (j in seq_len(min_length)) {
            if (base_comp[j] != current_path[j]) {
                j <- j - 1
                break
            }
        }

        base_comp <- base_comp[seq_len(j)]
    }

    if (return_full) {
        return(base_comp)
    }

    return(base_comp[length(base_comp)])
}
