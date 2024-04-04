# imports
    # standard libraries
import numpy as np
import anndata as ad
import pickle
from concurrent.futures import ProcessPoolExecutor # parallel computation (kinda)
    # user libraries
from scanorama import find_alignments
from library import mwjmsi

# constants
FOLDER = "data/pancreas/"

# functions
def score_batches(adata1,adata2):
    _, _, t = find_alignments([adata1.X,adata2.X])
    return 0 if (0,1) not in t else t[(0,1)]

def scorer_batch(anndatas):

    futures = {}
    with ProcessPoolExecutor() as executor:
        for i in range(len(anndatas)):
            for j in range(i+1, len(anndatas)):
                futures[(anndatas[i].obs['tech'][0],anndatas[j].obs['tech'][0])] = executor.submit(score_batches,anndatas[i],anndatas[j])
        executor.shutdown(wait=True)

    batch_scores = {k: v.result() for k, v in futures.items()} # compute batch scores

    with open(f"{FOLDER}pancreas-batch-scores.pkl", "wb") as f:
        pickle.dump(batch_scores, f)

def scorer_celltype(anndatas):
    cell_types = [] # get all of the cell types
    for anndata in anndatas:
        cell_types.append(anndata.obs['celltype'].cat.categories)
    cell_types = np.unique(cell_types)

    celltype_scores = {} # compute celltype scores

    with ProcessPoolExecutor() as executor:
        for i in range(len(anndatas)):
            anndata_i = anndatas[i]
            batch_i = anndatas[i].obs['tech'][0]
            for j in range(i+1, len(anndatas)):
                anndata_j = anndatas[j]
                batch_j = anndatas[j].obs['tech'][0]
                celltype_scores[(batch_i,batch_j)] = executor.submit(mwjmsi, anndata_i, anndata_j, cell_types)
        executor.shutdown(wait=True)
    
    celltype_scores = {k: v.result() for k, v in celltype_scores.items()} # compute celltype_scores

    with open(f"{FOLDER}pancreas-celltype-scores.pkl", "wb") as f:
        pickle.dump(celltype_scores, f)

if __name__ == "__main__":
    with open(f"{FOLDER}p.pkl", "wb") as f:
        pickle.dump("", f)