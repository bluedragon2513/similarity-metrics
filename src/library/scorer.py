# imports
    # standard libraries
import os
import numpy as np
import pandas as pd
from anndata import AnnData
from tqdm import trange
from typing import List, Callable
    # user libraries
from src.metrics.Scanorama.scanorama import scanorama
from src.library.combiner import mwjmsi, amwjmsi, gjsi
from src.library.scorer import *
from src.library.adata_preprocessing import *
from src.library.json_handler import *
# from seurat import * // for anthony to run


# functions
def score(
        adatas: List[AnnData],
        data_folder: str="data/", 
        algorithm: Callable=scanorama, 
        dataset_name: str="pancreas", 
        batch_save_file: str="batch-scores", 
        celltype_save_file: str="celltype-scores",
        combiner: Callable=amwjmsi, 
        batch_normalize: Callable=None,
        celltype_normalize: Callable=None,
        processing: str="",
        rerun=True,
        **kwargs
    ) -> None:
    folder = f"{data_folder}/{dataset_name}/{algorithm.__name__}/{processing}/"
    if rerun or not os.path.isfile(folder + "gjsi-scores.json"):
        scorer_gjsi(adatas, 
                    data_folder=data_folder, 
                    algorithm=algorithm, 
                    dataset_name=dataset_name, 
                    processing=processing)
    if rerun or not os.path.isfile(folder + "batch-scores.json"):
        scorer_batch(adatas, 
                    data_folder=data_folder, 
                    algorithm=algorithm, 
                    dataset_name=dataset_name, 
                    save_file=batch_save_file,
                    batch_normalize=batch_normalize,
                    processing=processing,
                    **kwargs)
    if rerun or not os.path.isfile(folder + f"celltype-scores-{combiner.__name__}.json"):
        scorer_celltype(adatas, 
                        data_folder=data_folder, 
                        algorithm=algorithm, 
                        dataset_name=dataset_name, 
                        save_file=celltype_save_file, 
                        combiner=combiner,
                        batch_normalize=batch_normalize,
                        celltype_normalize=celltype_normalize,
                        processing=processing,
                        **kwargs)

def scorer_gjsi(
        anndatas: List[AnnData], 
        data_folder: str="data/",
        algorithm: Callable=scanorama, 
        dataset_name: str="pancreas", 
        save_file: str="gjsi-scores",
        processing: str="") -> None:
    dict = {}
    cell_types = np.unique([anndata.obs['celltype'].cat.categories for anndata in anndatas])
    for i in trange(len(anndatas), desc="Scoring GJSI: "):
        for j in range(i, len(anndatas)):
            b1, b2 = anndatas[i], anndatas[j]
            b1name, b2name = b1.obs['tech'].iloc[0], b2.obs['tech'].iloc[0]
            score = gjsi(b1, b2, cell_types)

            dict[(b1name, b2name)] = score

    folder = f"{data_folder}/{dataset_name}/{algorithm.__name__}/{processing}/"
    json_writer(f"{folder}/{save_file}.json", dict)

def scorer_batch(
        anndatas: List[AnnData], 
        data_folder: str="data/",
        algorithm: Callable=scanorama, 
        dataset_name: str="pancreas", 
        save_file: str="batch-scores", 
        batch_normalize: Callable=None,
        processing: str="",
        **kwargs) -> None:
    batch_scores = {} # compute batch scores
    # batch_scores_jaccard = {}
    for i in trange(len(anndatas), desc="Scoring Batches: "):
        for j in range(i, len(anndatas)):
            b1, b2 = anndatas[i], anndatas[j]

            b1name, b2name = b1.obs['tech'].iloc[0], b2.obs['tech'].iloc[0]
            score = algorithm([b1.X,b2.X], normalize=batch_normalize, **kwargs)

            batch_scores[(b1name, b2name)] = score
    
    # save file
    folder = f"{data_folder}/{dataset_name}/{algorithm.__name__}/{processing}/"
    json_writer(f"{folder}/{save_file}.json", batch_scores)

def scorer_celltype(
        anndatas: List[AnnData],
        data_folder: str="data/", 
        algorithm: Callable=scanorama, 
        dataset_name: str="pancreas", 
        save_file: str="celltype-scores", 
        combiner: Callable=amwjmsi, 
        batch_normalize: Callable=None,
        celltype_normalize: Callable=None,
        processing: str="",
        **kwargs) -> None:
    # get all of the cell types
    cell_types = np.unique([anndata.obs['celltype'].cat.categories for anndata in anndatas]) 

    celltype_scores = {} # compute celltype scores
    for i in trange(len(anndatas), desc="Scoring Cell Types: "):
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
    
    folder = f"{data_folder}/{dataset_name}/{algorithm.__name__}/{processing}/"
    json_writer(f"{folder}/{save_file}-{combiner.__name__}.json", celltype_scores)

def score_seurat(path: str) -> None:
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

def score_celltypes_seurat(path: str) -> None:
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
