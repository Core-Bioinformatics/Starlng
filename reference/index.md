# Package index

## Shiny Write Utilities

- [`starlng_write_app_clustassess()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_clustassess.md)
  : Create Starlng Shiny app from a ClustAssess object
- [`starlng_write_app_clustassess_app()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_clustassess_app.md)
  : Create Starlng Shiny app from a ClustAssess object
- [`starlng_write_app_default()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_default.md)
  : Create Starlng Shiny app from a normalized gene expression matrix
- [`starlng_write_app_monocle()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_monocle.md)
  : Create Starlng Shiny app from Monocle object

## Shiny App Components

- [`ui_functional_assessment()`](https://core-bioinformatics.github.io/Starlng/reference/ui_functional_assessment.md)
  : UI - Functional Assessment
- [`ui_gene_clustering()`](https://core-bioinformatics.github.io/Starlng/reference/ui_gene_clustering.md)
  : UI - Gene clustering
- [`ui_gene_info_table()`](https://core-bioinformatics.github.io/Starlng/reference/ui_gene_info_table.md)
  : UI - Gene Info Table
- [`ui_gene_umap()`](https://core-bioinformatics.github.io/Starlng/reference/ui_gene_umap.md)
  : UI - Gene UMAP
- [`ui_global_setttings()`](https://core-bioinformatics.github.io/Starlng/reference/ui_global_setttings.md)
  : UI - Global Settings
- [`ui_metadata_umap()`](https://core-bioinformatics.github.io/Starlng/reference/ui_metadata_umap.md)
  : UI - Metadata UMAP
- [`ui_module_metadata_heatmap()`](https://core-bioinformatics.github.io/Starlng/reference/ui_module_metadata_heatmap.md)
  : UI - Gene Module Heatmap
- [`ui_module_umap()`](https://core-bioinformatics.github.io/Starlng/reference/ui_module_umap.md)
  : UI - Gene Module UMAP
- [`ui_pseudotime_select_cells()`](https://core-bioinformatics.github.io/Starlng/reference/ui_pseudotime_select_cells.md)
  : UI - Pseudotime
- [`server_functional_assessment()`](https://core-bioinformatics.github.io/Starlng/reference/server_functional_assessment.md)
  : Server - Functional Assessment
- [`server_gene_clustering()`](https://core-bioinformatics.github.io/Starlng/reference/server_gene_clustering.md)
  : Server - Gene clustering
- [`server_gene_info_table()`](https://core-bioinformatics.github.io/Starlng/reference/server_gene_info_table.md)
  : Server - Gene Info Table
- [`server_gene_umap()`](https://core-bioinformatics.github.io/Starlng/reference/server_gene_umap.md)
  : Server - Gene UMAP
- [`server_metadata_umap()`](https://core-bioinformatics.github.io/Starlng/reference/server_metadata_umap.md)
  : Server - Metadata UMAP
- [`server_module_metadata_heatmap()`](https://core-bioinformatics.github.io/Starlng/reference/server_module_metadata_heatmap.md)
  : Server - Gene Module Heatmap
- [`server_module_umap()`](https://core-bioinformatics.github.io/Starlng/reference/server_module_umap.md)
  : Server - Gene Module UMAP
- [`server_pseudotime_select_cells()`](https://core-bioinformatics.github.io/Starlng/reference/server_pseudotime_select_cells.md)
  : Server - Pseudotime
- [`prepare_session()`](https://core-bioinformatics.github.io/Starlng/reference/prepare_session.md)
  : Server - Prepare Session
- [`update_gears_width()`](https://core-bioinformatics.github.io/Starlng/reference/update_gears_width.md)
  : Server - Gear Width
- [`update_tabs()`](https://core-bioinformatics.github.io/Starlng/reference/update_tabs.md)
  : Server - Tabs update

## Gene/Cell Utilities

- [`filter_central_cells_from_group()`](https://core-bioinformatics.github.io/Starlng/reference/filter_central_cells_from_group.md)
  : Find central points of a group
- [`remove_outlier_cells()`](https://core-bioinformatics.github.io/Starlng/reference/remove_outlier_cells.md)
  : Remove outlier cells from a group
- [`select_cells_by_gene_expr()`](https://core-bioinformatics.github.io/Starlng/reference/select_cells_by_gene_expr.md)
  : Select cells by gene expression
- [`select_cells_by_metadata()`](https://core-bioinformatics.github.io/Starlng/reference/select_cells_by_metadata.md)
  : Select cells by metadata
- [`calculate_umap_average_distance()`](https://core-bioinformatics.github.io/Starlng/reference/calculate_umap_average_distance.md)
  : Calculate Average UMAP Distance
- [`voting_scheme()`](https://core-bioinformatics.github.io/Starlng/reference/voting_scheme.md)
  : Gene Expression voting scheme
- [`read_gene_from_dense_h5()`](https://core-bioinformatics.github.io/Starlng/reference/read_gene_from_dense_h5.md)
  : Read expression matrix from a dense HDF5 file
- [`sort_genes_by_metadata()`](https://core-bioinformatics.github.io/Starlng/reference/sort_genes_by_metadata.md)
  : Sort genes by metadata information
- [`write_gene_matrix_dense_h5()`](https://core-bioinformatics.github.io/Starlng/reference/write_gene_matrix_dense_h5.md)
  : Write expression matrix to a dense HDF5 file
- [`write_gene_matrix_sparse_h5()`](https://core-bioinformatics.github.io/Starlng/reference/write_gene_matrix_sparse_h5.md)
  : Write expression matrix to a sparse HDF5 file
- [`plot_umap()`](https://core-bioinformatics.github.io/Starlng/reference/plot_umap.md)
  : Plot UMAP with Discrete or Continuous Coloring

## Monocle3 Utilities

- [`diet_monocle_object()`](https://core-bioinformatics.github.io/Starlng/reference/diet_monocle_object.md)
  : Diet Monocle object
- [`subset_monocle_by_trajectory()`](https://core-bioinformatics.github.io/Starlng/reference/subset_monocle_by_trajectory.md)
  : Subset the Monocle object by a trajectory
- [`update_monocle_partition()`](https://core-bioinformatics.github.io/Starlng/reference/update_monocle_partition.md)
  : Update the partition of the Monocle object

## Pseudotime Utilities

- [`custom_learn_graph()`](https://core-bioinformatics.github.io/Starlng/reference/custom_learn_graph.md)
  : Learn the graph of the Monocle object
- [`get_trajectory_object()`](https://core-bioinformatics.github.io/Starlng/reference/get_trajectory_object.md)
  : Build a Trajectory Helper Object
- [`plot_trajectory_graph()`](https://core-bioinformatics.github.io/Starlng/reference/plot_trajectory_graph.md)
  : Plot a Trajectory Graph
- [`custom_pseudotime_ordering()`](https://core-bioinformatics.github.io/Starlng/reference/custom_pseudotime_ordering.md)
  : Order the cells by pseudotime
- [`get_pseudotime_recommendation()`](https://core-bioinformatics.github.io/Starlng/reference/get_pseudotime_recommendation.md)
  : Recommend a pseudotime ordering
- [`order_metadata_groups_by_pseudotime()`](https://core-bioinformatics.github.io/Starlng/reference/order_metadata_groups_by_pseudotime.md)
  : Order Metadata Groups by Pseudotime

## Stability Utilities

- [`clustering_pipeline()`](https://core-bioinformatics.github.io/Starlng/reference/clustering_pipeline.md)
  : Community detection pipeline
- [`community_detection_master()`](https://core-bioinformatics.github.io/Starlng/reference/community_detection_master.md)
  : Community detection wrapper
- [`get_clusters_consistency()`](https://core-bioinformatics.github.io/Starlng/reference/get_clusters_consistency.md)
  : Get the by-cluster consistency of configurations
- [`get_feature_loading()`](https://core-bioinformatics.github.io/Starlng/reference/get_feature_loading.md)
  : Get feature loading from an expression matrix
- [`group_by_clusters_general()`](https://core-bioinformatics.github.io/Starlng/reference/group_by_clusters_general.md)
  : Group the partitions by the number of clusters
- [`parallel_nn2_idx()`](https://core-bioinformatics.github.io/Starlng/reference/parallel_nn2_idx.md)
  : Parallel NN indexing
- [`pca_reduction()`](https://core-bioinformatics.github.io/Starlng/reference/pca_reduction.md)
  : Get PCA reduction from an expression matrix
- [`select_best_configuration()`](https://core-bioinformatics.github.io/Starlng/reference/select_best_configuration.md)
  : Select the most stable configurations of parameters

## Module Processing

- [`build_module_masks()`](https://core-bioinformatics.github.io/Starlng/reference/build_module_masks.md)
  : Build Module Masks
- [`compute_module_pairwise_tables()`](https://core-bioinformatics.github.io/Starlng/reference/compute_module_pairwise_tables.md)
  : Compute Module Pairwise Tables
- [`detect_outlier()`](https://core-bioinformatics.github.io/Starlng/reference/detect_outlier.md)
  : Detect Module Outliers
- [`get_closest_node_to_module()`](https://core-bioinformatics.github.io/Starlng/reference/get_closest_node_to_module.md)
  : Find Closest Trajectory Node to Module
- [`get_filtered_gene_adjacency()`](https://core-bioinformatics.github.io/Starlng/reference/get_filtered_gene_adjacency.md)
  : Filter Gene Adjacency by Module Structure
- [`get_module_centroid()`](https://core-bioinformatics.github.io/Starlng/reference/get_module_centroid.md)
  : Compute Module Centroid in Embedding Space
- [`get_module_stats()`](https://core-bioinformatics.github.io/Starlng/reference/get_module_stats.md)
  : Compute Module Statistics
- [`get_module_transitions()`](https://core-bioinformatics.github.io/Starlng/reference/get_module_transitions.md)
  : Build Module Transition Adjacency
- [`plot_module_transitions()`](https://core-bioinformatics.github.io/Starlng/reference/plot_module_transitions.md)
  : Plot Module Transition Graph
- [`summarise_module_stats()`](https://core-bioinformatics.github.io/Starlng/reference/summarise_module_stats.md)
  : Summarize Module Statistics
- [`get_per_module_weight()`](https://core-bioinformatics.github.io/Starlng/reference/get_per_module_weight.md)
  : Compute Per-Module Gene Weights
- [`get_gene_overlap_stat()`](https://core-bioinformatics.github.io/Starlng/reference/get_gene_overlap_stat.md)
  : Compute Gene Hub Statistics
- [`plot_umap_gene_modules()`](https://core-bioinformatics.github.io/Starlng/reference/plot_umap_gene_modules.md)
  : Plot UMAP for Multiple Gene Modules
- [`plot_module_pseudobulk_expression()`](https://core-bioinformatics.github.io/Starlng/reference/plot_module_pseudobulk_expression.md)
  : Plot Module Pseudobulk Expression
- [`generate_cell_heatmap()`](https://core-bioinformatics.github.io/Starlng/reference/generate_cell_heatmap.md)
  : Generate Cell Heatmap
- [`plot_gene_hub_umap()`](https://core-bioinformatics.github.io/Starlng/reference/plot_gene_hub_umap.md)
  : Plot Gene UMAP with Hub Highlighting
- [`plot_module_trends_over_pseudotime()`](https://core-bioinformatics.github.io/Starlng/reference/plot_module_trends_over_pseudotime.md)
  : Plot Module Trends over Pseudotime

## Enrichment Analysis

- [`filter_enrichment_results()`](https://core-bioinformatics.github.io/Starlng/reference/filter_enrichment_results.md)
  : Filter Enrichment Results
- [`plot_enrichment_top_terms()`](https://core-bioinformatics.github.io/Starlng/reference/plot_enrichment_top_terms.md)
  : Plot Top Enriched Terms

## TF Analysis

- [`get_transcription_factors()`](https://core-bioinformatics.github.io/Starlng/reference/get_transcription_factors.md)
  : Get Enriched Transcription Factors
- [`get_tf_stats()`](https://core-bioinformatics.github.io/Starlng/reference/get_tf_stats.md)
  : Summarize Transcription Factor Statistics
- [`add_tf_hub_stats()`](https://core-bioinformatics.github.io/Starlng/reference/add_tf_hub_stats.md)
  : Add Hub-Gene Counts to TF Statistics
- [`get_tf_gene_network()`](https://core-bioinformatics.github.io/Starlng/reference/get_tf_gene_network.md)
  : Build TF-Gene Network Graph
- [`plot_module_tfs_ggraph()`](https://core-bioinformatics.github.io/Starlng/reference/plot_module_tfs_ggraph.md)
  : Plot TF-Gene Network with ggraph
- [`plot_module_tfs_bubbleplot()`](https://core-bioinformatics.github.io/Starlng/reference/plot_module_tfs_bubbleplot.md)
  : Plot TF Bubble Plot
