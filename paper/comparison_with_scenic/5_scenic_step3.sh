#!/bin/bash

pyscenic aucell \
    seurat_info.loom \
    reg.csv \
    --output pyscenic_output.loom \
    --num_workers 10
