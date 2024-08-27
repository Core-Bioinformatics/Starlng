# ## plotting trajectory
#   if (show_trajectory_graph) {

#     ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
#       as.data.frame() %>%
#       dplyr::select(prin_graph_dim_1 = {{x}}, prin_graph_dim_2 = {{y}}) %>%
#       dplyr::mutate(sample_name = rownames(.),
#                     sample_state = rownames(.))

#     dp_mst <- cds@principal_graph[[reduction_method]]

#     edge_df <- dp_mst %>%
#       igraph::as_data_frame() %>%
#       dplyr::select(source = "from", target = "to") %>%
#       dplyr::left_join(ica_space_df %>%
#                          dplyr::select(
#                            source="sample_name",
#                            source_prin_graph_dim_1="prin_graph_dim_1",
#                            source_prin_graph_dim_2="prin_graph_dim_2"),
#                        by = "source") %>%
#       dplyr::left_join(ica_space_df %>%
#                          dplyr::select(
#                            target="sample_name",
#                            target_prin_graph_dim_1="prin_graph_dim_1",
#                            target_prin_graph_dim_2="prin_graph_dim_2"),
#                        by = "target")
#   }


#   if (show_trajectory_graph){
#     g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
#                                      y="source_prin_graph_dim_2",
#                                      xend="target_prin_graph_dim_1",
#                                      yend="target_prin_graph_dim_2"),
#                           size=trajectory_graph_segment_size,
#                           color=I(trajectory_graph_color),
#                           linetype="solid",
#                           na.rm=TRUE,
#                           data=edge_df)


#     if (label_principal_points) {
#       mst_branch_nodes <- branch_nodes(cds, reduction_method)
#       mst_leaf_nodes <- leaf_nodes(cds, reduction_method)
#       mst_root_nodes <- root_nodes(cds, reduction_method)
#       pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
#       princ_point_df <- ica_space_df %>%
#         dplyr::slice(match(names(pps), sample_name))

#       g <- g +
#         geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
#                    shape = 21, stroke=I(trajectory_graph_segment_size),
#                    color="white",
#                    fill="black",
#                    size=I(graph_label_size * 1.5),
#                    na.rm=TRUE, princ_point_df) +
#         ggrepel::geom_text_repel(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
#                              label="sample_name"),
#                   size=I(graph_label_size * 1.5), color="Black", na.rm=TRUE,
#                   princ_point_df)
#     }
#     if (label_branch_points){
#       mst_branch_nodes <- branch_nodes(cds, reduction_method)
#       branch_point_df <- ica_space_df %>%
#         dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
#         dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))

#       g <- g +
#         geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
#                    shape = 21, stroke=I(trajectory_graph_segment_size),
#                    color="white",
#                    fill="black",
#                    size=I(graph_label_size * 1.5),
#                    na.rm=TRUE, branch_point_df) +
#         geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
#                              label="branch_point_idx"),
#                   size=I(graph_label_size), color="white", na.rm=TRUE,
#                   branch_point_df)
#     }

#     if (label_leaves){
#       mst_leaf_nodes <- leaf_nodes(cds, reduction_method)
#       leaf_df <- ica_space_df %>%
#         dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
#         dplyr::mutate(leaf_idx = seq_len(dplyr::n()))

#       g <- g +
#         geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
#                    shape = 21, stroke=I(trajectory_graph_segment_size),
#                    color="black",
#                    fill="lightgray",
#                    size=I(graph_label_size * 1.5),
#                    na.rm=TRUE,
#                    leaf_df) +
#         geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
#                              label="leaf_idx"),
#                   size=I(graph_label_size), color="black", na.rm=TRUE, leaf_df)
#     }

#     if (label_roots){
#       mst_root_nodes <- root_nodes(cds, reduction_method)
#       root_df <- ica_space_df %>%
#         dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
#         dplyr::mutate(root_idx = seq_len(dplyr::n()))

#       g <- g +
#         geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
#                    shape = 21, stroke=I(trajectory_graph_segment_size),
#                    color="black",
#                    fill="white",
#                    size=I(graph_label_size * 1.5),
#                    na.rm=TRUE,
#                    root_df) +
#         geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
#                              label="root_idx"),
#                   size=I(graph_label_size), color="black", na.rm=TRUE, root_df)
#     }

#   }

#' Diet Monocle object
#' 
#' @description This function simplifies the structure of the Monocle object.
#' The purpose is to create a lightweight object that can be used in the
#' context of the Shiny application. The resulting object will only contain
#' the needed information to perform the pseudotime ordering. THerefore,
#' this function removes information such as: the expression matrices, the
#' auxiliary data related to the dimensionality reduction and the clustering
#' and the metadata columns.
#' 
#' @param mon_obj Monocle object.
#' 
#' @return Diet Monocle object.
#' @export
diet_monocle_object <- function(mon_obj) {
    mask <- rep(TRUE, nrow(mon_obj) * ncol(mon_obj))
    # remove all count information
    mon_obj@assays@data@listData$counts[mask] <- 0
    # delete the other assays
    for (slot_names in names(mon_obj@assays@data@listData)) {
        if (slot_names == "counts") {
            next
        }
        mon_obj@assays@data@listData[[slot_names]] <- NULL
    }

    # keep only one metadata column
    mon_obj@colData[, 2:ncol(mon_obj@colData)] <- NULL

    # keep only the UMAP dim reduction
    for (dim_reduc_name in names(mon_obj@reduce_dim_aux@listData)) {
        if (dim_reduc_name == "UMAP") {
            next
        }
        mon_obj@reduce_dim_aux@listData[[dim_reduc_name]] <- NULL
    }

    for (dim_reduc_name in names(mon_obj@int_colData@listData$reducedDims@listData)) {
        if (dim_reduc_name == "UMAP") {
            next
        }
        mon_obj@int_colData@listData$reducedDims@listData[[dim_reduc_name]] <- NULL
    }

    for (dim_reduc_name in names(mon_obj@clusters)) {
        if (dim_reduc_name == "UMAP") {
            next
        }

        mon_obj@clusters[[dim_reduc_name]] <- NULL
    }

    mon_obj@clusters$UMAP$cluster_result <- NULL

    return(mon_obj)
}
