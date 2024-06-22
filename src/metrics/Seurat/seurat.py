from scib import preprocessing
import anndata as ad
import pickle
import scanpy as sc
import rpy2.robjects as robjects
import importlib.util
import sys
sys.path.append("src/library/adata_preprocessing.py")
# from src.library.adata_preprocessing import scib_normalize
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
    datasets[0] = check_dataset(datasets[0])
    datasets[1] = check_dataset(datasets[1])
    print(datasets[0].n_obs, datasets[1].n_obs)
    adata_to_seurat(datasets[0], "dataset_1")
    adata_to_seurat(datasets[1], "dataset_2")
    if(datasets[0].n_obs == 6 or datasets[1].n_obs == 6):
        run_seurat("Seurat/datasets", 5, 5)
    elif(datasets[0].n_obs <= 30 or datasets[1].n_obs <= 30):
        run_seurat("Seurat/datasets", min(datasets[0].n_obs, datasets[1].n_obs) - 1, min(datasets[0].n_obs, datasets[1].n_obs) - 1)
    # elif(datasets[0].n_obs <= 50 or datasets[1].n_obs <= 50):
    #     run_seurat("Seurat/datasets", min(datasets[0].n_obs, datasets[1].n_obs) - 1, 30)
    else:
        run_seurat("Seurat/datasets", 30, 30)
    score = pd.read_csv('data/pancreas/celltype_scores_seurat.csv')
    retval = float(score.iloc[0].iloc[1])
    print(retval)
    return retval

def check_dataset(dataset):
    if(dataset.n_obs < 6):
        dataset = addCells(dataset, 6 - dataset.n_obs)
    return dataset

def addCells(dataset, add):
    num_features = dataset.n_vars
    new_data = np.concatenate((dataset.X, np.random.randn(add, num_features)), axis=0)
    new_dataset = ad.AnnData(new_data)
    new_dataset.layers['counts'] = new_dataset.X
    return new_dataset


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
        path = f"Seurat/datasets/{name}.rds",
        batch="sample_id",
    )

def run_seurat(path, dims, score):
    robjects.r("source('src/metrics/Seurat/Seurat.R')")
    function = robjects.r('run_seurat')
    r_string = robjects.StrVector([path])
    r_dims = robjects.IntVector([dims])
    r_score = robjects.IntVector([score])
    function(r_string, r_dims, r_score)


if __name__ == "__main__":
    adata = read_file("data/pancreas/human_pancreas_norm_complexBatch.h5ad")
    dict = rds_filter_batches(adata)
    retval = seurat([dict['inDrop2'], dict['inDrop4']], normalize=False)