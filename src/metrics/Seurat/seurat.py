from scib import preprocessing
import anndata as ad
import pickle
import scanpy as sc
import rpy2.robjects as robjects
import importlib.util
import sys
sys.path.append("src/library/adata_preprocessing.py")
from library.adata_preprocessing import scib_normalize
import pandas as pd
import numpy as np


def seurat(datasets, normalize=False, **kwargs):
    if normalize:
        adata = ad.AnnData(np.concatenate((datasets[0], datasets[1])))
        adata = scib_normalize(adata)
        datasets[0] = ad.AnnData(adata.X[:len(datasets[0])])
        datasets[1] = ad.AnnData(adata.X[len(datasets[0]):])
        datasets[0].layers['counts'] = datasets[0].X
        datasets[1].layers['counts'] = datasets[1].X
    else:
        adata1 = ad.AnnData(datasets[0])
        adata1.layers['counts'] = adata1.X
        datasets[0] = adata1
        adata2 = ad.AnnData(datasets[1])
        adata2.layers['counts'] = adata2.X
        datasets[1] = adata2

    print(datasets[0].n_obs, datasets[1].n_obs)
    print(datasets[0].n_obs, datasets[1].n_obs)
    adata_to_seurat(datasets[0], "dataset_1")
    adata_to_seurat(datasets[1], "dataset_2")
    if(datasets[0].n_obs <= 50 or datasets[1].n_obs <= 50):
        run_seurat("Seurat/data", min(datasets[0].n_obs, datasets[1].n_obs))
    else:
         run_seurat("Seurat/data", 50)
    score = pd.read_csv('data/pancreas/celltype_scores_seurat.csv')
    print(score.iloc[0].iloc[1])
    return score.iloc[0].iloc[1]


def read_file(path):
    adata = sc.read_h5ad(path)
    return adata


def rds_filter_batches(adata):
    adatas = {}
    for c in adata.obs['tech'].cat.categories:
        new_adata = ad.AnnData(
            X=adata.layers['counts'][adata.obs['tech'] == c, :],
            obs=adata.obs[adata.obs['tech'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas[c] = new_adata
    return adatas

def adata_to_seurat(adatas, name):
    preprocessing.save_seurat(
        adata = adatas,
        path = f"Seurat/data/{name}.rds",
        batch="sample_id",
    )

def run_seurat(path, dims):
    robjects.r("source('src/metrics/Seurat/Seurat.R')")
    function = robjects.r('run_seurat')
    r_string = robjects.StrVector([path])
    r_int = robjects.StrVector([dims])
    function(r_string, dims)


if __name__ == "__main__":
    adata = read_file("data/unpreprocessed/human_pancreas_norm_complexBatch.h5ad")
    dict = rds_filter_batches(adata)
    retval = seurat([dict['inDrop1'], dict['inDrop1']], normalize=False)