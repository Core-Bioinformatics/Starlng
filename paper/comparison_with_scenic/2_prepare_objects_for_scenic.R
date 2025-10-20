library(Seurat)
library(SCENIC)
library(SingleCellExperiment)
library(SCopeLoomR)

add_cell_annotation <- function(loom, cellAnnotation) {
  cellAnnotation <- data.frame(cellAnnotation)
  if(any(c("nGene", "nUMI") %in% colnames(cellAnnotation)))
  {
    warning("Columns 'nGene' and 'nUMI' will not be added as annotations to the loom file.")
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nGene", drop=FALSE]
    cellAnnotation <- cellAnnotation[,colnames(cellAnnotation) != "nUMI", drop=FALSE]
  }
  
  if(ncol(cellAnnotation)<=0) stop("The cell annotation contains no columns")
  if(!all(get_cell_ids(loom) %in% rownames(cellAnnotation))) stop("Cell IDs are missing in the annotation")
  
  cellAnnotation <- cellAnnotation[get_cell_ids(loom),,drop=FALSE]
  # Add annotation
  for(cn in colnames(cellAnnotation))
  {
    add_col_attr(loom=loom, key=cn, value=cellAnnotation[,cn])
  }
  
  invisible(loom)
}

setwd("paper/comparison_with_scenic")
so <- readRDS("immuneCellsSCTransformed.rds")
so <- SCTransform(so, variable.features.n = 3000, return.only.var.genes = FALSE, verbose = FALSE)
so <- RunPCA(so, npcs = 50, verbose = FALSE)
so <- RunUMAP(so, reduction = "pca", dims = 1:50, verbose = FALSE)

cc_genes <- cc.genes.updated.2019
so <- CellCycleScoring(so, s.features = cc_genes$s.genes, g2m.features = cc_genes$g2m.genes)
qs2::qs_save(so, "immuneCellsSCTransformed.qs2")

expr_matrix <- so@assays$SCT@data
meta_data <- so@meta.data

loom <- build_loom(file.path("seurat_info.loom"), dgem = expr_matrix)
loom <- add_cell_annotation(loom, meta_data)
close_loom(loom)
