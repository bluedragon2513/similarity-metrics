options(repos = c(CRAN = "https://cran.r-project.org"))

# imports described in Seurat's DESCRIPTION
my_list <- list(
    "Seurat", "Rcpp", "dplyr", "ggplot2", "cowplot",
    "Matrix", "reticulate", "RANN", "RColorBrewer", "igraph", "irlba",
    "leiden", "lmtest", "MASS", "matrixStats", "pbapply", "plotly",
    "png", "progressr", "purrr", "RcppAnnoy", "RcppHNSW", "RSpectra",
    "Rtsne", "scattermore", "sctransform", "shiny", "spatstat.explore",
    "spatstat.geom", "tibble", "uwot", "fitdistrplus", "future",
    "future.apply", "generics", "ggrepel", "ggridges", "httr", "ica",
    "jsonlite", "lifecycle", "miniUI", "patchwork", "ROCR",
    "rlang", "tools", "utils", "stats", "ape", "rsvd", "testthat",
    "hdf5r", "Rfast2", "VGAM", "enrichR", "mixtools", 
    "data.table", "R.utils", "harmony", "cluster",
    "fastDummies", "SeuratObject", "graphics", "grDevices",
    "grid", "scales", "RcppEigen", "RcppProgress",
    "lazyData", "arrow", "BiocManager", "KernSmooth"
)

for (pkg in my_list) {
    library(pkg, character.only = TRUE)
}

bioc_list <- c(
    "S4Vectors", "SummarizedExperiment", "SingleCellExperiment",
    "BiocGenerics", "GenomicRanges", "GenomeInfoDb", "IRanges",
    "Biobase", "DelayedArray"
)

BiocManager::install(bioc_list, ask = FALSE, update = FALSE)

# import from Seurat's library 
source("Seurat/seurat_code/R/utilities.R")
source("Seurat/seurat_code/R/clustering.R")
source("Seurat/seurat_code/R/convenience.R")
source("Seurat/seurat_code/R/data.R")
source("Seurat/seurat_code/R/dimensional_reduction.R")
source("Seurat/seurat_code/R/mixscape.R")
source("Seurat/seurat_code/R/objects.R")
source("Seurat/seurat_code/R/RcppExports.R")
source("Seurat/seurat_code/R/reexports.R")
source("Seurat/seurat_code/R/integration.R")

# Preprocess Seurat Objects in list
#
# @param list: list containing Seurat Objects (list)
#
# @return: return list of preprocessed data (list)
preprocessing <- function(list) {
    for (i in 1:length(list)) {
        list[[i]] <- ScaleData(list[[i]])
    }
    return(list)
}

# Use Seurat's library and FindIntegrationAnchors to obtain score between two datasets
# and write them to a csv file. 
#
# @param path: path that represents path to rds file (String)
# @param dims: dimensions to reduce by (int)
# @param score: used for k.score (int)
# @param anchor: used for k.anchor (int)
run_seurat <- function(path, dims, score, anchor) {
    # Read Surat Objects from rds file
    directory <- path
    file_names <- list.files(path = directory, pattern = ".rds")
    data.list <- list()
    for (i in seq_along(file_names)) {
        seurat_object <- readRDS(paste0(directory, "/", file_names[i]))
        data.list <- c(data.list, seurat_object)
    }
    # list of MNNs obtained in Seurat
    mnn_list <- FindIntegrationAnchors(data.list,
                                        dims = 1:dims,
                                        k.score = score,
                                        k.anchor = anchor)
    
    # Use scanorama's scoring formula
    cells_1 <- ncol(data.list[[1]])
    cells_2 <- ncol(data.list[[2]])
    score1 <- mnn_list[[1]][[1]] / cells_1
    score2 <- mnn_list[[1]][[2]] / cells_2
    score <- max(score1, score2)
    if (score > 1) {
        score <- min(score1, score2)
    }
    if (score > 1) {
        score <- 1
    }
    write.csv(score, "/Users/carsonnannini/Research/similarity-metrics/data/pancreas/celltype_scores_seurat.csv")
}
