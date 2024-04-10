# imports
    # standard libraries
import numpy as np
import anndata as ad
import pickle
    # user libraries
from scanorama import find_alignments
from library import mwjmsi, amwjmsi

# functions
def scorer_batch(anndatas, folder="data/", algorithm="scanorama", dataset_name="pancreas", save_file="batch-scores"):
    batch_scores = {} # compute batch scores
    for i in range(len(anndatas)):
        for j in range(i, len(anndatas)):
            _, _, t = find_alignments([anndatas[i].X,anndatas[j].X])
            batch_scores[(anndatas[i].obs['tech'].iloc[0],anndatas[j].obs['tech'].iloc[0])] = 0 if (0,1) not in t else t[(0,1)]
    
    # save file
    with open(f"{folder}{algorithm}/{dataset_name}/{save_file}.pkl", "wb") as f:
        pickle.dump(batch_scores, f)

def scorer_celltype(anndatas, folder="data/", algorithm="scanorama", dataset_name="pancreas", save_file="celltype-scores", similarity_function=amwjmsi):
    # get all of the cell types
    cell_types = np.unique([anndata.obs['celltype'].cat.categories for anndata in anndatas]) 

    celltype_scores = {} # compute celltype scores
    for i in range(len(anndatas)):
        anndata_i = anndatas[i]
        batch_i = anndatas[i].obs['tech'].iloc[0]

        for j in range(i, len(anndatas)):
            anndata_j = anndatas[j]
            batch_j = anndatas[j].obs['tech'].iloc[0]
            celltype_scores[(batch_i,batch_j)] = similarity_function(anndata_i, anndata_j, cell_types)
    
    print("I am scoring itX")
    with open(f"{folder}{algorithm}/{dataset_name}/{save_file}-{similarity_function.__name__}.pkl", "wb") as f:
        pickle.dump(celltype_scores, f)

if __name__ == "__main__":
    with open(f"{0}p.pkl", "wb") as f:
        pickle.dump("", f)