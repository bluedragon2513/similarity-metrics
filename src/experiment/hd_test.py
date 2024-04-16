import numpy as np 
# from sklearn.decomposition import PCA
# from sklearn.datasets import load_iris
import anndata as ad
import scanpy as sc # how to read the data
from concurrent.futures import ProcessPoolExecutor
import pickle

from modules.metrics import MetricsData

def filter_batches(adata):
    adatas = []
    for c in adata.obs['tech'].cat.categories:
        new_adata = ad.AnnData(
            X=adata.X[adata.obs['tech'] == c, :],
            obs=adata.obs[adata.obs['tech'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas.append(new_adata)
    return adatas

def filter_celltypes(adata):
    adatas = []
    for c in adata.obs['celltype'].cat.categories:
        new_adata = ad.AnnData(
            X=adata.X[adata.obs['celltype'] == c, :],
            obs=adata.obs[adata.obs['celltype'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas.append(new_adata)
    return adatas

def scorer_celltype(anndatas, folder="src/", algorithm="kruskal", dataset_name="pancreas", save_file="celltype-scores"):
    # get all of the cell types
    # cell_types = np.unique([anndata.obs['celltype'].cat.categories for anndata in anndatas])

    celltype_scores = {} # compute celltype scores
    for i in range(len(anndatas)):
        anndata_i = anndatas[i].X
        batch_i = anndatas[i].obs['tech'].iloc[0]

        for j in range(i, len(anndatas)):
            anndata_j = anndatas[j].X
            batch_j = anndatas[j].obs['tech'].iloc[0]
            M = MetricsData(np.transpose(anndata_i), np.transpose(anndata_j))
            celltype_scores[(batch_i,batch_j)] = M.compute_stress_kruskal()
    
    print("I am scoring itX")
    with open(f"{folder}experiment/{algorithm}/{dataset_name}/{save_file}-kruskal.pkl", "wb") as f:
        pickle.dump(celltype_scores, f)

def min_count(adata1, adata2, all_cell_types):
    mins = []
    for ct1 in all_cell_types:
        arr = []
        mins.append(arr)
        len_ct1 = len(adata1[adata1.obs['celltype'] == ct1])
        for ct2 in all_cell_types:
            len_ct2 = len(adata2[adata2.obs['celltype'] == ct2])
            arr.append(min(len_ct1,len_ct2))
    return mins

    # helper function for mwjmsi
def max_count(adata1, adata2, all_cell_types):
    count = np.sum([max(len(adata1[adata1.obs['celltype'] == ct]),len(adata2[adata2.obs['celltype'] == ct])) for ct in all_cell_types])
    return count

def mwjmsi(adata1, adata2, all_cell_types):
    adatas1 = [adata1.X[adata1.obs['celltype'] == c] for c in all_cell_types]
    adatas2 = [adata2.X[adata2.obs['celltype'] == c] for c in all_cell_types]

    mins = min_count(adata1, adata2, all_cell_types)
    sum = 0
    for ct in range(len(all_cell_types)):
        # if there are no cells of either adata, then ignore it
        if adatas1[ct].shape[0] == 0 or adatas2[ct].shape[0] == 0:
            continue
        M = MetricsData(np.transpose(adatas1[ct]), np.transpose(adatas2[ct]))
        sum += mins[ct][ct] * M.compute_stress_kruskal()

    return sum/max_count(adata1, adata2, all_cell_types)

if __name__ == "__main__":
    # get data
        # 
    adata = sc.read_h5ad("data/preprocessed/human_pancreas_preprocessed.h5ad")
    adata = sc.pp.pca(data=adata, copy=True, n_comps=50)
    adatas = filter_batches(adata)
    scorer_celltype(adatas)

    dict = {}
    for i in range(len(adatas)):
        for j in range(i,len(adatas)):
            batch_pair = (adatas[i].obs['tech'].iloc[0], adatas[j].obs['tech'].iloc[0])
            dict[batch_pair] = mwjmsi(adatas[i], adatas[j], np.unique([adata.obs['celltype'].cat.categories for adata in adatas]))
    # store dictionary
    with open("src/experiment/data/batch-stress-scores.pkl", "wb") as f:
        pickle.dump(dict, f)

    
    
