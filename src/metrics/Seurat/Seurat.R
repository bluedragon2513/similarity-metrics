library(Seurat)
library(future)
library(pbapply)
library(future.apply)
library(glue)
source("Seurat/seurat_code/utilities.R")
source("Seurat/seurat_code/clustering.R")
source("Seurat/seurat_code/convenience.R")
source("Seurat/seurat_code/data.R")
source("Seurat/seurat_code/dimensional_reduction.R")
source("Seurat/seurat_code/mixscape.R")
source("Seurat/seurat_code/objects.R")
source("Seurat/seurat_code/RcppExports.R")
source("Seurat/seurat_code/reexports.R")


preprocessing <- function(list) {
    for (i in 1:length(list)) {
        list[[i]] <- NormalizeData(list[[i]])
        list[[i]] <- ScaleData(list[[i]])
        list[[i]] <- FindVariableFeatures(list[[i]])
    }
    return(list)
}

run_seurat <- function(path) {
    directory <- path
    file_names <- list.files(path = directory, pattern = ".rds")
    data.list <- list()
    for (i in seq_along(file_names)) {
        seurat_object <- readRDS(paste0(directory, "/", file_names[i]))
        data.list <- c(data.list, list(seurat_object))
    }
    anchors <- FindIntegrationAnchors(object.list = data.list, anchor.features = 30)
    source("Seurat/seurat_code/integration.R")
    matrix <- IntegrateData(anchors)
    data_frame <- as.data.frame(matrix)
    for(i in 1:length(file_names)) {
        names(data_frame)[names(data_frame) == glue("V{i}")] <- glue("{file_names[i]}")
        rownames(data_frame)[i] <- glue("{file_names[i]}")
    }
    write.csv(data_frame, "/Users/carsonnannini/Research/similarity-metrics/data/pancreas/celltype_scores_seurat.csv")
}