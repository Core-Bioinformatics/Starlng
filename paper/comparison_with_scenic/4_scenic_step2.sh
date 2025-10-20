#!/bin/bash

pyscenic ctx \
    expr_mat.adjacencies.tsv \
    hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname seurat_info.loom \
    --output reg.csv \
    --mask_dropouts \
    --num_workers 10
