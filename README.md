# STARLNG: STability Analysis of coReguLated Nests of Genes

[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/Starlng)](https://github.com/r-hub/cranlogs.app)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/grand-total/Starlng)](https://github.com/r-hub/cranlogs.app)
[![](https://img.shields.io/github/languages/code-size/Core-Bioinformatics/Starlng.svg)](https://github.com/Core-Bioinformatics/Starlng)
[![CRAN status](https://www.r-pkg.org/badges/version/Starlng?color=orange)](https://CRAN.R-project.org/package=Starlng)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/Core-Bioinformatics/Starlng/main?style=flat&color=white)
<!-- [![](https://img.shields.io/github/r-package/v/Core-Bioinformatics/Starlng%2Frelease-1.1.0?label=devel%20version&color=green)](https://github.com/Core-Bioinformatics/Starlng/tree/release-1.1.0) -->

This repository contains the Starlng R package, which identifies stable clusters of coexpressed genes and describes their position along the pseudotime trajectory. The package builds on top of the Monocle3 [1] and ClustAssess [2] frameworks.

## Installation
<!-- Starlng can be install from CRAN -->

<!-- `install.packages("Starlng")` -->

or from Github using the `remotes` package:

`remotes::install_github("Core-Bioinformatics/Starlng")`.

The following packages are required for Starlng:

* circlize
* ClustAssess
* ComplexHeatmap
* dplyr
* DT
* foreach
* ggplot2
* Gmedian
* gprofiler2
* HDF5Array
* igraph
* leidenbase
* Matrix (>= 1.5.0)
* methods
* monocle3
* patchwork
* qs
* qs2
* qualpalr
* RANN
* rclipboard
* reshape2
* RhpcBLASctl
* rhdf5
* shiny
* shinyjs
* shinyWidgets
* spsComps
* tidyr
* stringr
* viridis

We suggest installing the following packages for optimal performance:
* doFuture
* doParallel
* irlba
* testthat (>= 3.0.0)
* parallel
* plotly
* SharedObject


## Citing Starlng
To be added.

## References
[1] J. Cao, M. Spielmann, X. Qiu, X. Huang, D. M. Ibrahim, A. J. Hill, F. Zhang, S. Mundlos, L. Christiansen, F. J. Steemers, C. Trapnell, and J. Shendure, “The single-cell transcriptional landscape of mammalian organogenesis,” Nature, vol. 566, p. 496–502, Feb. 2019.

[2] A. Shahsavari, A. Munteanu, and I. Mohorianu, “Clustassess: tools for assessing the robustness of single-cell clustering,” bioRxiv, 2022.
