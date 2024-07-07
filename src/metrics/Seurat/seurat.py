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

# Obtain the score between two datasets using Seurat's library
#
# @param datasets: the two datasets being evaluated (list)
# @param normalize: use scib_normalization if true (boolean)
#
# @return: return score between datasets (float)
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
    get_min = min(datasets[0].n_obs, datasets[1].n_obs)
    if(datasets[0].n_obs == 3 or datasets[1].n_obs == 3):
        run_seurat("Seurat/datasets", 2, 2, 2)
    elif(datasets[0].n_obs < 6 or datasets[1].n_obs < 6):
        run_seurat("Seurat/datasets", get_min - 1, get_min - 1, get_min - 1)
    elif(datasets[0].n_obs <= 30 or datasets[1].n_obs <= 30):
        run_seurat("Seurat/datasets", get_min - 1, get_min - 1, 5)
    # elif(datasets[0].n_obs <= 50 or datasets[1].n_obs <= 50):
    #     run_seurat("Seurat/datasets", min(datasets[0].n_obs, datasets[1].n_obs) - 1, 30)
    else:
        run_seurat("Seurat/datasets", 30, 30, 5)
    score = pd.read_csv('data/pancreas/celltype_scores_seurat.csv')
    retval = float(score.iloc[0].iloc[1])
    print(retval)
    return retval

# If dataset has less than 3 cells, call addCells. Otherwise return origional dataset. 
#
# @param dataset: dataset being checked. (AnnData)
#
# @return: Return appropriate dataset. (AnnData)
def check_dataset(dataset):
    if(dataset.n_obs < 3):
        dataset = addCells(dataset, 3 - dataset.n_obs)
    return dataset

# Replicate exsisting cells 'add' times. Used when the number of cells is less than 3. 
#
# @param dataset: dataset to modify (AnnData)
# @param add: cells to add (int)
# @param add_value: add 'add_value' to replicated cell so dataset can work with Seurat (float)
#
# @ return: modified dataset
def addCells(dataset, add, add_value=0.1):
    num_cells = dataset.n_obs
    indices_to_duplicate = np.random.choice(num_cells, add, replace=True)
    duplicated_cells = dataset.X[indices_to_duplicate, :]
    duplicated_cells += add_value
    duplicated_cells[duplicated_cells < 0] = 0
    new_data = np.concatenate((dataset.X, duplicated_cells), axis=0)
    new_dataset = ad.AnnData(new_data)
    new_dataset.layers['counts'] = new_dataset.X
    return new_dataset

# Read file path for data being evaluated
#
# @param path: relative path of file (String)
#
# @return: adata 
def read_file(path):
    adata = sc.read_h5ad(path)
    return adata

# Filter and store data from file into a dictionary
#
# @param adata: adata being stored
#
# @return: object containing data, (dictionary)
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

# Convert adata into a Seurat Object and store in an rds file.
#
# @param adatas: adata being converted to a Seurat Object (dictionary)
# @param name: name of file storing Seurat Objects
def adata_to_seurat(adatas, name):
    preprocessing.save_seurat(
        adata = adatas,
        path = f"Seurat/datasets/{name}.rds",
        batch="sample_id",
    )

# Convert parameters into R objects and call 'run_seurat' in Seurat.R
#
# @param path: path to rds file containing Seurat Objects (String)
# @param dims: parameter for dims (int)
# @param score: parameter for k.score (int)
# @param anchor: parameter for k.anchor (int)
def run_seurat(path, dims, score, anchor):
    robjects.r("source('src/metrics/Seurat/Seurat.R')")
    function = robjects.r('run_seurat')
    r_string = robjects.StrVector([path])
    r_dims = robjects.IntVector([dims])
    r_score = robjects.IntVector([score])
    r_anchor = robjects.IntVector([anchor])
    function(r_string, r_dims, r_score, r_anchor)

# main
if __name__ == "__main__":
    adata = read_file("data/pancreas/human_pancreas_norm_complexBatch.h5ad")
    dict = rds_filter_batches(adata)
    retval = seurat([dict['celseq'], dict['fluidigmc1']], normalize=False)