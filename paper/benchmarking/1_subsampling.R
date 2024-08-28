library(zinbwave)
library(Seurat)
# library(MASS)
library(digest)
library(mclust)
library(parallel)

simulateW <- function(zinb, ncells = 100, nclust = 3, ratioSSW_SSB = 1, colIni = 1){
  par(mfrow = c(2,2))
  # zinbW
  xlim = c(min(zinb@W[,1]) - 1, max(zinb@W[,1]) + 1) 
  ylim = c(min(zinb@W[,2]) - 1, max(zinb@W[,2]) + 1)
  plot(zinb@W, col = colIni, xlim = xlim, ylim = ylim,
       main = 'zinb fitted W\ncolor = brain area')
  
  # mclustW
  mclustW = Mclust(zinb@W, G = nclust)
  plot(mclustW$data, main = 'multivar gauss fit\nmclust K = 3', xlab = 'W1', ylab = 'W2', 
       xlim = xlim, ylim = ylim, col = mclustW$classification)
  
  # multivar gaussian
  clust = sample(mclustW$classification, ncells, replace = TRUE)
  stopifnot(length(unique(clust)) == nclust)
  # a = b = 1
  simW1 = lapply(clust, function(i){
    MASS::mvrnorm(n = 1, mu = mclustW$parameters$mean[, i], 
            Sigma = mclustW$parameters$variance$sigma[,, i])
  })
  simW1 = do.call(rbind, simW1)
  plot(simW1, col = clust,
       main = paste0('multivar gauss sim\nncells=', ncells, ', scaleRatioSSW_SSB = 1'),
       xlab = 'W1', ylab = 'W2', xlim = xlim, ylim = ylim)
  
  simW2 = computeNewW(simW1, clust, ratioSSW_SSB)
  plot(simW2, col = clust,
       main = paste0('multivar gauss sim\nncells=', ncells, ', scaleRatioSSW_SSB =' , ratioSSW_SSB),
       xlab = 'W1', ylab = 'W2', xlim = xlim, ylim = ylim)
  par(mfrow = c(1, 1))
  
  return(list(simW = simW2, bio = clust))
}

computeNewW = function(W, labels, ratioSSW_SSB = 1){
  Vtot = apply(W, 2, var)
  W_bar = apply(W, 2, mean)
  N = nrow(W)
  
  cc = table(labels)
  ccNames = names(cc)
  nk = as.vector(cc)
  W_bar_c = sapply(as.numeric(ccNames), function(i){
    apply(W[labels == i, ], 2, mean)
  })
  colnames(W_bar_c) = ccNames
  SS_between = colSums(nk*t((W_bar_c - W_bar)^2))
  
  SS_within = rowSums(sapply(seq_len(N), function(i){
    (W[i, ] - W_bar_c[, as.character(labels[i])])^2
  }))

  b2 = ratioSSW_SSB
  a = sqrt( ( (N-1) * Vtot ) / ( SS_between + b2 * SS_within ) )
  W_start = sapply(seq_len(N), function(i){
    (1 - a) * W_bar +  W_bar_c[, as.character(labels[i])] * a * (1 - sqrt(b2)) +  a * sqrt(b2) * W[i,]
  })
  t(W_start)
}

simulateGamma <- function(zinb, ncells = 100, gammapiOffset = 0, colIni = 1,
                          colSim = 1){
  # gamma zinb
  gamma = data.frame(gammaMu = zinb@gamma_mu[1, ],
                     gammaPi = zinb@gamma_pi[1, ])
  # mclustW
  mclustGamma = Mclust(gamma, G = 1)
  
  # multivar gaussian
  simGamma = MASS::mvrnorm(n = ncells,
                     mu = mclustGamma$parameters$mean[,1] + c(0, gammapiOffset), 
                     Sigma = mclustGamma$parameters$variance$sigma[,, 1])
  
  par(mfrow = c(1,2))
  xlim = c(min(c(gamma[,1], simGamma[,1])) - .5,
           max(c(gamma[,1], simGamma[,1])) + .5) 
  ylim = c(min(c(gamma[,2], simGamma[,2])) - .5,
           max(c(gamma[,2], simGamma[,2])) + .5) 
  plot(gamma[,1], gamma[,2], col = colIni,
       xlim = xlim, ylim = ylim, xlab = 'gamma_mu', ylab = 'gamma_pi', 
       main = 'zinb fitted Gamma\ncolor = brain area')
  plot(simGamma, col = colSim,
       main = paste0('bivar gauss sim\nncells=', ncells, 
                     ', gammaPi offset = ', gammapiOffset),
       xlab = 'gamma_mu', ylab = 'gamma_pi', xlim = xlim, ylim = ylim)
  par(mfrow = c(1, 1))
  
  return(simGamma)
}


zinbSimWrapper <- function(core, colIni, ncells = 100, ngenes = 1000, nclust = 3, 
                           ratioSSW_SSB = 1, gammapiOffset = 0, B = 1, ncores = NULL,
                           fileName = 'zinbSim.rda', model_fileName = NULL){
  # sample ngenes 
  set.seed(9128)
  if (ngenes > nrow(core)) repl = T else repl = F
  core = core[sample(1:nrow(core), ngenes, repl = repl),]
  
  # fit zinb (if you already fitted zinb, it is cached)
  d = digest(core, "md5")
  if (is.null(model_fileName)) {
    file_basename <- strsplit(fileName, "\\.")[[1]][1]
    # tmp = paste0(tempdir(), '/', d)
    model_fileName = sprintf("%s_zinb.rda", file_basename)
  }
  if (!file.exists(model_fileName)){
    print(glue::glue('run ZINB and saving to {model_fileName}'))
    if (is.null(ncores)) ncores = max(1, detectCores() - 1)
    zinb <- zinbFit(
      core,
      K = 2,
      commondispersion = FALSE,
      epsilon = ngenes,
      verbose = TRUE,
      BPPARAM = BiocParallel::SnowParam(workers = ncores))
    save(zinb, file = model_fileName)
  }else{
    load(model_fileName)
  }
  
  # sim W
  print("Simulating W")
  w = simulateW(zinb, ncells, nclust, ratioSSW_SSB, colIni)
  simW = w$simW
  bio = w$bio
  
  # sim gamma
  print("Simulating Gama")
  simGamma = simulateGamma(zinb, ncells, gammapiOffset,
                           colIni = colIni, colSim = bio)
  
  # sim model
  print("Simulating model")
  simModel = zinbModel(W=simW, gamma_mu = matrix(simGamma[,1], nrow = 1),
                  gamma_pi = matrix(simGamma[,2], nrow = 1),
                  alpha_mu=zinb@alpha_mu, alpha_pi=zinb@alpha_pi,
                  beta_mu=zinb@beta_mu, beta_pi=zinb@beta_pi, zeta = zinb@zeta)
  
  # sim data
  if (B == 1){
    simData = zinbSim(simModel, seed = 1)
  } else{
    simData = lapply(seq_len(B), function(j){
      zinbSim(simModel, seed = j)
    })
  }

  print("Saving")
  
  save(bio, simModel, simData, file = fileName)
}

setwd("benchmarking")
dts_name <- "immune"
so <- readRDS(paste0("../", dts_name, "CellsSCTransformed.rds"))
# pc_genes <- read.csv("~/protein_coding_gene_names.txt", header = FALSE)[, 1]
# so <- subset(so, features = pc_genes)
so <- SCTransform(so, return.only.var.genes = FALSE, variable.features.n = 3000)
most_abundant_genes <- rownames(so@assays$RNA@counts)[order(rowSums(so@assays$RNA@counts), decreasing = TRUE)[seq_len(3000)]]
keep_genes <- so@assays$SCT@var.features
keep_genes <- union(keep_genes, most_abundant_genes[seq_along(keep_genes)])
so <- subset(so, features = keep_genes)

cnt_matrix <- so@assays$RNA@counts
cnt_matrix <- cnt_matrix[rowSums(cnt_matrix) > 0,]

dir.create(dts_name, showWarnings = FALSE)

for (ncells in seq(from = 10000, to = 100000, by = 10000)) {
    print(ncells)
    zinbSimWrapper(
      cnt_matrix,
      colIni = 1,
      ncells = ncells,
      ngenes = 3000,
      nclust = 3,
      ratioSSW_SSB = 1,
      gammapiOffset = 0,
      B = 1,
      ncores = 30,
      fileName = file.path(dts_name, paste0('zinbSim_', ncells, ".rda")),
      model_fileName = file.path(dts_name, "zinbSim_zinb.rda")
    )
}

