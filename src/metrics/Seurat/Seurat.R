library(Seurat)
library(glue)
library(SeuratObject)
library(future)
library(pbapply)
library(future.apply)
options(repos = c(CRAN = "https://cran.r-project.org"))

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
    "hdf5r", "S4Vectors", "SummarizedExperiment", "SingleCellExperiment",
    "BiocGenerics", "GenomicRanges", "GenomeInfoDb", 
    "IRanges", "Rfast2", "Biobase", "VGAM",
    "enrichR", "mixtools", "data.table",
    "R.utils", "DelayedArray", "harmony"
)

for (pkg in my_list) {
    library(pkg, character.only = TRUE)
}
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

preprocessing <- function(list) {
    for (i in 1:length(list)) {
        list[[i]] <- NormalizeData(list[[i]])
        list[[i]] <- ScaleData(list[[i]])
        list[[i]] <- FindVariableFeatures(list[[i]])
    }
    return(list)
}

run_seurat <- function(path, dims) {
    directory <- path
    file_names <- list.files(path = directory, pattern = ".rds")
    data.list <- list()
    for (i in seq_along(file_names)) {
        seurat_object <- readRDS(paste0(directory, "/", file_names[i]))
        data.list <- c(data.list, list(seurat_object))
    }
    anchors <- FindIntegrationAnchors(data.list,
                                      l2.norm = FALSE,
                                      reduction = "cca",
                                      normalization.method = NULL,
                                      dims = 1:dims,
                                      k.anchor = 20,
                                      k.filter = 800,
                                      k.score = 20)
    cells_1 <- ncol(data.list[[1]])
    cells_2 <- ncol(data.list[[2]])
    score1 <- anchors[[1]][[1]] / cells_1
    score2 <- anchors[[1]][[2]] / cells_2
    score <- max(score1, score2)
    if (score > 1) {
        score <- min(score1, score2)
    }
    if (score > 1) {
        score <- 1
    }
    # score <- max(anchors[[1]][[1]] / cells_1, anchors[[1]][[2]] / cells_2)
    # matrix <- IntegrateData(anchors, normalization.method = NULL)
    # data_frame <- as.data.frame(matrix)
    # for(i in 1:length(file_names)) {
    #     names(data_frame)[names(data_frame) == glue("V{i}")] <- glue("{file_names[i]}")
    #     rownames(data_frame)[i] <- glue("{file_names[i]}")
    # }
    write.csv(score, "/Users/carsonnannini/Research/similarity-metrics/data/pancreas/celltype_scores_seurat.csv")
}
