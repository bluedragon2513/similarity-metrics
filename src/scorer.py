# imports
    # standard libraries
import numpy as np
import anndata as ad
import pickle
from tqdm import trange
    # user libraries
from scanorama import scanorama
from library import mwjmsi, amwjmsi, gjsi
from scorer import *
from adata_preprocessing import *
import pandas as pd
# from seurat import * // for anthony to run

# functions
def scorer_batch(
        anndatas, 
        data_folder="data/", 
        algorithm=scanorama, 
        dataset_name="pancreas", 
        save_file="batch-scores", 
        batch_normalize=None,
        processing="",
        **kwargs):
    batch_scores = {} # compute batch scores
    batch_scores_jaccard = {}
    for i in trange(len(anndatas)):
        for j in range(i, len(anndatas)):
            b1, b2 = anndatas[i], anndatas[j]
            score = algorithm([b1.X,b2.X], normalize=batch_normalize, **kwargs)
            b1name, b2name = b1.obs['tech'].iloc[0],b2.obs['tech'].iloc[0]
            batch_scores[(b1name, b2name)] = score
            batch_scores_jaccard[(b1name, b2name)] = gjsi(b1, b2, np.unique([anndata.obs['celltype'].cat.categories for anndata in anndatas]))
    
    # save file
    with open(f"{data_folder}/{dataset_name}/{algorithm.__name__}/{processing}/{save_file}.pkl", "wb") as f:
        pickle.dump(batch_scores, f)
    with open(f"{data_folder}/{dataset_name}/{algorithm.__name__}/{processing}/{save_file}-jaccard.pkl", "wb") as f:
        pickle.dump(batch_scores_jaccard, f)

def scorer_celltype(
        anndatas, 
        data_folder="data/", 
        algorithm=scanorama, 
        dataset_name="pancreas", 
        save_file="celltype-scores", 
        combiner=amwjmsi, 
        batch_normalize=None,
        celltype_normalize=None,
        processing="",
        **kwargs):
    # get all of the cell types
    cell_types = np.unique([anndata.obs['celltype'].cat.categories for anndata in anndatas]) 

    celltype_scores = {} # compute celltype scores
    for i in trange(len(anndatas)):
        anndata_i = anndatas[i]
        batch_i = anndatas[i].obs['tech'].iloc[0]

        for j in range(i, len(anndatas)):
            anndata_j = anndatas[j]
            batch_j = anndatas[j].obs['tech'].iloc[0]
            celltype_scores[(batch_i,batch_j)] = combiner(
                anndata_i, 
                anndata_j, 
                cell_types, 
                algorithm=algorithm, 
                batch_normalize=batch_normalize, 
                celltype_normalize=celltype_normalize,
                **kwargs)
    
    print("I am scoring itX")
    with open(f"{data_folder}/{dataset_name}/{algorithm.__name__}/{processing}/{save_file}-{combiner.__name__}.pkl", "wb") as f:
        pickle.dump(celltype_scores, f)

def score_seurat(path):
    # adata = read_file(path)
    # adata_to_seurat(adata, 'tech')
    # string = "Seurat/data/tech"
    # run_seurat(string)
    df = pd.read_csv('data/pancreas/seurat_scores.csv')
    for column in df.columns:
        df[column].fillna(0, inplace=True)
    df = df.iloc[:, 1:]
    df_dict = df.to_dict(orient='records')
    dict_list = [{k[:-4]:v for k,v in x.items()} for x in df_dict]
    processed_dict = {}
    for i in range(len(dict_list)):
        store_key = list(dict_list[i].keys())[i]
        for key in dict_list[i].keys():
            if dict_list[i][key] != 0:
                new_key = (store_key, key)
                processed_dict[new_key] = dict_list[i][key] / 10
    return processed_dict

def score_celltypes_seurat(path):
    # adata = read_file(path)
    # adata_to_seurat(adata, 'celltype')
    # string = "Seurat/data/celltype"
    # run_seurat(string)
    df = pd.read_csv('data/pancreas/celltype_scores_seurat.csv')
    for column in df.columns:
        df[column].fillna(0, inplace=True)
    df = df.iloc[:, 1:]
    df_dict = df.to_dict(orient='records')
    dict_list = [{k[:-4]:v for k,v in x.items()} for x in df_dict]
    processed_dict = {}
    for i in range(len(dict_list)):
        store_key = list(dict_list[i].keys())[i]
        for key in dict_list[i].keys():
            if dict_list[i][key] != 0:
                new_key = (store_key, key)
                processed_dict[new_key] = dict_list[i][key] / 10
    return processed_dict

if __name__ == "__main__":
    with open(f"{0}p.pkl", "wb") as f:
        pickle.dump("", f)