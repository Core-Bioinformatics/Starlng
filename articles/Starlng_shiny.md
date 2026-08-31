# Generating the Starlng Shiny Application on single-cell data

In this vignette we will demonstrate how to apply the Starlng shiny
function to single-cell RNA-seq data. We will exemplify this on the
PBMC3k dataset.

``` r

library(Starlng)
library(Seurat)
#> Loading required package: SeuratObject
#> Warning: package 'SeuratObject' was built under R version 4.4.1
#> Loading required package: sp
#> Warning: package 'sp' was built under R version 4.4.1
#> 
#> Attaching package: 'SeuratObject'
#> The following objects are masked from 'package:base':
#> 
#>     intersect, t
library(SeuratData)
#> Warning: package 'SeuratData' was built under R version 4.4.1
#> ── Installed datasets ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── SeuratData v0.2.2.9002 ──
#> ✔ pbmc3k       3.1.4                                                                                                            ✔ pbmcMultiome 0.1.4
#> ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── Key ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> ✔ Dataset loaded successfully
#> ❯ Dataset built with a newer version of Seurat than installed
#> ❓ Unknown version of Seurat installed

print(packageVersion("Starlng"))
#> [1] '1.0.0'
```

We repeated the procedures of downloading, qc filtering, normalising and
dimensionality reduction as described in the `ClustAssess`
[vignette](https://github.com/Core-Bioinformatics/ClustAssess/blob/main/vignettes/stability-pipeline-description.Rmd).

``` r

InstallData("pbmc3k")
#> Warning: The following packages are already installed and will not be reinstalled: pbmc3k
data("pbmc3k")
pbmc3k <- UpdateSeuratObject(pbmc3k)
#> Validating object structure
#> Updating object slots
#> Ensuring keys are in the proper structure
#> Warning: Assay RNA changing from Assay to Assay
#> Ensuring keys are in the proper structure
#> Ensuring feature names don't have underscores or pipes
#> Updating slots in RNA
#> Validating object structure for Assay 'RNA'
#> Object representation is consistent with the most current Seurat version
pbmc3k <- PercentageFeatureSet(pbmc3k, pattern = "^MT-", col.name = "percent.mito")
pbmc3k <- PercentageFeatureSet(pbmc3k, pattern = "^RP[SL][[:digit:]]", col.name = "percent.rp")
all.index <- seq_len(nrow(pbmc3k))
MT.index <- grep(pattern = "^MT-", x = rownames(pbmc3k), value = FALSE)
RP.index <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(pbmc3k), value = FALSE)
pbmc3k <- pbmc3k[!((all.index %in% MT.index) | (all.index %in% RP.index)), ]
pbmc3k <- subset(pbmc3k, nFeature_RNA < 2000 & nCount_RNA < 2500 & percent.mito < 7 & percent.rp > 7)
pbmc3k <- NormalizeData(pbmc3k, verbose = FALSE)
pbmc3k <- FindVariableFeatures(pbmc3k, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
pbmc3k <- ScaleData(pbmc3k, features = rownames(pbmc3k), verbose = FALSE)
pbmc3k <- RunPCA(pbmc3k,
    npcs = 30,
    approx = FALSE,
    verbose = FALSE
)
pbmc3k <- RunUMAP(pbmc3k,
    reduction = "pca",
    dims = 1:30,
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    verbose = FALSE
)
#> Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
#> To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
#> This message will be shown once per session
pbmc3k
#> An object of class Seurat 
#> 13607 features across 2401 samples within 1 assay 
#> Active assay: RNA (13607 features, 3000 variable features)
#>  3 layers present: counts, data, scale.data
#>  2 dimensional reductions calculated: pca, umap
```

Applying `Starlng` is straightforward, using the
[`starlng_write_app_default()`](https://core-bioinformatics.github.io/Starlng/reference/starlng_write_app_default.md)
function. The function requires a normalised expression matrix, a
metadata dataframe, a PCA and UMAP embedding. The function allows
controlling the number of points used in inferring the pseudotime
trajectory (`learn_graph_parameters`), filtering the genes based on
autocorrelation (`gene_filtering_function`) and the parameters used in
the stability assessment of gene clustering (`clustering_parameters`).

``` r

starlng_write_app_default(
    folder_path = "starlng_app_pbmc3k",
    expression_matrix = GetAssayData(pbmc3k, assay = "RNA", slot = "data"),
    metadata_df = pbmc3k@meta.data,
    pca_embedding = Embeddings(pbmc3k, reduction = "pca"),
    umap_embedding = Embeddings(pbmc3k, reduction = "umap"),
    app_title_name = "Starlng PBMC3K Example"
)
```

The function will create a Starlng Shiny application in the specified
folder. The application will contain most of the functionalities
described in the main `Starlng` vignette. By default, for each stable
gene module, the population summaries (using a scale threshold of 0.5)
and analysis will be already precomputed, but the user can choose to
perform the analysis on the fly with other parameters.

You can run the application using the following command:

``` r

shiny::runApp("starlng_app_pbmc3k")
```

## Session Info

``` r

sessionInfo()
#> R version 4.4.0 (2024-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 22.04.5 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C           LC_TELEPHONE=C        
#> [11] LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: Europe/Bucharest
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] pbmcMultiome.SeuratData_0.1.4 pbmc3k.SeuratData_3.1.4       SeuratData_0.2.2.9002         Seurat_5.1.0                  SeuratObject_5.0.2            sp_2.1-4                      Starlng_1.0.0                
#> 
#> loaded via a namespace (and not attached):
#>   [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3          rlang_1.3.0            magrittr_2.0.3         RcppAnnoy_0.0.22       otel_0.2.0             spatstat.geom_3.2-9    matrixStats_1.3.0      ggridges_0.5.6        
#>  [11] compiler_4.4.0         png_0.1-8              vctrs_0.6.5            reshape2_1.4.4         stringr_1.5.1          crayon_1.5.3           pkgconfig_2.0.3        fastmap_1.2.0          utf8_1.2.4             promises_1.3.0        
#>  [21] purrr_1.0.2            xfun_0.56              jsonlite_2.0.0         goftest_1.2-3          later_1.3.2            spatstat.utils_3.1-2   irlba_2.3.5.1          parallel_4.4.0         cluster_2.1.6          R6_2.6.1              
#>  [31] ica_1.0-3              stringi_1.8.4          RColorBrewer_1.1-3     spatstat.data_3.0-4    reticulate_1.37.0      parallelly_1.37.1      lmtest_0.9-40          scattermore_1.2        Rcpp_1.0.13            iterators_1.0.14      
#>  [41] knitr_1.51             tensor_1.5             future.apply_1.11.2    zoo_1.8-12             sctransform_0.4.1      httpuv_1.6.15          Matrix_1.7-0           splines_4.4.0          igraph_2.2.1           tidyselect_1.2.1      
#>  [51] abind_1.4-5            spatstat.random_3.2-3  codetools_0.2-19       miniUI_0.1.2           spatstat.explore_3.2-7 listenv_0.9.1          lattice_0.22-5         tibble_3.2.1           plyr_1.8.9             shiny_1.8.1.1         
#>  [61] S7_0.2.1               ROCR_1.0-11            evaluate_1.0.5         Rtsne_0.17             future_1.33.2          fastDummies_1.7.3      survival_3.5-8         polyclip_1.10-6        fitdistrplus_1.1-11    pillar_1.9.0          
#>  [71] KernSmooth_2.23-22     foreach_1.5.2          plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0         ggplot2_4.0.3          scales_1.4.0           globals_0.16.3         xtable_1.8-4           glue_1.7.0            
#>  [81] lazyeval_0.2.2         tools_4.4.0            data.table_1.15.4      RSpectra_0.16-2        RANN_2.6.1             leiden_0.4.3.1         dotCall64_1.1-1        cowplot_1.1.3          grid_4.4.0             tidyr_1.3.1           
#>  [91] colorspace_2.1-1       nlme_3.1-163           patchwork_1.2.0        cli_3.6.6              rappdirs_0.3.3         spatstat.sparse_3.0-3  spam_2.10-0            fansi_1.0.6            viridisLite_0.4.2      dplyr_1.1.4           
#> [101] uwot_0.2.2             gtable_0.3.6           digest_0.6.35          progressr_0.14.0       ggrepel_0.9.5          htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.8.1      lifecycle_1.0.5        httr_1.4.7            
#> [111] mime_0.12              MASS_7.3-60
```
