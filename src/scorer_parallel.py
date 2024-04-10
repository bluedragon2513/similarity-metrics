# imports
    # standard libraries
import numpy as np
import anndata as ad
import pickle
from concurrent.futures import ProcessPoolExecutor # parallel computation (kinda)
    # user libraries
from scanorama import find_alignments
from library import mwjmsi, amwjmsi

# constants
FOLDER = "data/pancreas/"

# functions
def score_batches(adata1,adata2):
    _, _, t = find_alignments([adata1.X,adata2.X])
    return 0 if (0,1) not in t else t[(0,1)]

def scorer_batch(anndatas, dataset_name="pancreas", folder="data/", save_file="batch-scores"):
    futures = {}
    with ProcessPoolExecutor() as executor:
        for i in range(len(anndatas)):
            for j in range(i, len(anndatas)):
                futures[(anndatas[i].obs['tech'].iloc[0],anndatas[j].obs['tech'].iloc[0])] = executor.submit(score_batches,anndatas[i],anndatas[j])
        executor.shutdown(wait=True)

    batch_scores = {k: v.result() for k, v in futures.items()} # compute batch scores

    # save file
    with open(f"{folder}{dataset_name}/{save_file}.pkl", "wb") as f:
        pickle.dump(batch_scores, f)

def scorer_celltype(anndatas, dataset_name="pancreas", folder="data/", save_file="celltype-scores", similarity_function=amwjmsi):
    # get all of the cell types
    cell_types = np.unique([anndata.obs['celltype'].cat.categories for anndata in anndatas]) 

    celltype_scores = {} # compute celltype scores

    with ProcessPoolExecutor() as executor:
        for i in range(len(anndatas)):
            anndata_i = anndatas[i]
            batch_i = anndatas[i].obs['tech'].iloc[0]

            for j in range(i, len(anndatas)):
                anndata_j = anndatas[j]
                batch_j = anndatas[j].obs['tech'].iloc[0]
                celltype_scores[(batch_i,batch_j)] = executor.submit(similarity_function, anndata_i, anndata_j, cell_types)
        
        executor.shutdown(wait=True)
    
    celltype_scores = {k: v.result() for k, v in celltype_scores.items()} # compute celltype_scores

    with open(f"{folder}{dataset_name}/{save_file}.pkl", "wb") as f:
        pickle.dump(celltype_scores, f)

if __name__ == "__main__":
    with open(f"{FOLDER}p.pkl", "wb") as f:
        pickle.dump("", f)