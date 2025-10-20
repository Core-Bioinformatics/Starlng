#!/bin/bash

pyscenic grn \
    --num_workers 10 \
    -o expr_mat.adjacencies.tsv \
    seurat_info.loom \
    allTFs_hg38.txt
