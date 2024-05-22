# imports
    # standard libraries
import numpy as np
import networkx as nx
from anndata import AnnData
from typing import List, Callable
    # user libraries
from src.metrics.Scanorama.scanorama import find_alignments, scanorama


# functions
    # helper function
def count(
        adata: AnnData, 
        ct: str) -> int:
    return len(adata[adata.obs['celltype'] == ct])

    # helper function for mwjmsi
def min_count(
        adata1: AnnData, 
        adata2: AnnData, 
        all_cell_types: List[str]) -> List[List[int]]:
    mins = []
    for ct1 in all_cell_types:
        len_ct1 = count(adata1, ct1)
        mins.append([min(len_ct1, count(adata2, ct2)) for ct2 in all_cell_types])
    return mins

    # helper function for mwjmsi
def max_count(
        adata1: AnnData, 
        adata2: AnnData, 
        all_cell_types: List[str]) -> int:
    return sum([max(count(adata1, ct), count(adata2, ct)) for ct in all_cell_types])

# Generalized Jaccard Similarity Index
def gjsi(
        adata1: AnnData, 
        adata2: AnnData, 
        all_cell_types: List[str]) -> float:
    mins = sum([min(count(adata1,ct),count(adata2,ct)) for ct in all_cell_types])
    maxs = max_count(adata1, adata2, all_cell_types)
    return mins/maxs

# Modified Weighted-Jaccard-Multiset Similarity Index
def mwjmsi(
        adata1: AnnData, 
        adata2: AnnData, 
        all_cell_types: List[str], 
        algorithm: Callable=scanorama, 
        batch_normalize: Callable=None, 
        celltype_normalize: Callable=None, 
        **kwargs) -> float:
    assert \
        (batch_normalize and not celltype_normalize) or \
        (not batch_normalize and celltype_normalize) or \
        (not batch_normalize and not celltype_normalize)
    
    if batch_normalize:
        adata = AnnData(np.concatenate((adata1.X, adata2.X)))
        adata = batch_normalize(adata)
        adata1 = AnnData(adata.X[:len(adata1.X)], obs=adata1.obs)
        adata2 = AnnData(adata.X[len(adata1.X):], obs=adata2.obs)

    adatas1 = [adata1.X[adata1.obs['celltype'] == c] for c in all_cell_types]
    adatas2 = [adata2.X[adata2.obs['celltype'] == c] for c in all_cell_types]

    mins = min_count(adata1, adata2, all_cell_types)
    sum = 0
    for ct in range(len(all_cell_types)):
        # if there are no cells of either adata, then ignore it
        if adatas1[ct].shape[0] == 0 or adatas2[ct].shape[0] == 0:
            continue
        datas = [adatas1[ct], adatas2[ct]]
        sum += mins[ct][ct] * algorithm(datas, normalize=celltype_normalize, **kwargs)

    return sum/max_count(adata1, adata2, all_cell_types)

# Adjusted Modified Weighted-Jaccard-Multiset Similarity Index
def amwjmsi(
        adata1: AnnData, 
        adata2: AnnData, 
        all_cell_types: List[str], 
        algorithm: Callable=scanorama, 
        batch_normalize: Callable=None, 
        celltype_normalize: Callable=None, 
        **kwargs) -> float:
    assert \
        (batch_normalize and not celltype_normalize) or \
        (not batch_normalize and celltype_normalize) or \
        (not batch_normalize and not celltype_normalize)
    
    if batch_normalize:
        adata = AnnData(np.concatenate((adata1.X, adata2.X)))
        adata = batch_normalize(adata)
        adata1 = AnnData(adata.X[:len(adata1.X)], obs=adata1.obs)
        adata2 = AnnData(adata.X[len(adata1.X):], obs=adata2.obs)
    
    adatas1 = [adata1.X[adata1.obs['celltype'] == c] for c in all_cell_types]
    adatas2 = [adata2.X[adata2.obs['celltype'] == c] for c in all_cell_types]

    mins = min_count(adata1, adata2, all_cell_types)
    sum = 0
    maxs = 0
    unmatched1 = [] # indices in adatasX that are not matched
    unmatched2 = [] # indices in adatasX that are not matched
    # compute matched indices
    for ct in range(len(all_cell_types)):
        if adatas1[ct].shape[0] == 0:
            if adatas2[ct].shape[0] != 0:
                unmatched2.append(ct)
            continue
        if adatas2[ct].shape[0] == 0:
            unmatched1.append(ct)
            continue
        datas = [adatas1[ct], adatas2[ct]]
        sum += mins[ct][ct] * algorithm(datas, normalize=celltype_normalize, **kwargs)
        maxs += max(len(adatas1[ct]), len(adatas2[ct]))

    # get maximum weight bipartite graphs
    edge_dict = {}
    edges = []
    for i in unmatched1:
        for j in unmatched2:
            weight = algorithm([adatas1[i], adatas2[j]], normalize=celltype_normalize, **kwargs)
            edges.append((i,j,weight))
            edge_dict[(i,j)] = weight
    # compute maximum weight edges
    G = nx.Graph()
    G.add_weighted_edges_from(edges)
    non_matching_matches = nx.max_weight_matching(G)
    matches = []
    # compute score
    for (i,j) in non_matching_matches:
        sum += min(len(adatas1[i]),len(adatas2[j])) * edge_dict[(i,j)]
        maxs += max(len(adatas1[i]), len(adatas2[j]))
        matches.append(i)
        matches.append(j)

    # get remaining unmatched types, add to the count
    remaining1 = [i for i in unmatched1 if i not in matches]
    remaining2 = [j for j in unmatched2 if j not in matches]
    for r in remaining1:
        maxs += len(adatas1[r])
    for r in remaining2:
        maxs += len(adatas2[r])

    # compute 
    return sum/maxs

# Modified Weighted-Jaccard-Multiset Similarity Index
def experiment_amwjmsi(
        adata1: AnnData, 
        adata2: AnnData, 
        all_cell_types: List[str], 
        algorithm: Callable=scanorama, 
        normalize: bool=False) -> float:
    adatas1 = [adata1.X[adata1.obs['celltype'] == c] for c in all_cell_types]
    adatas2 = [adata2.X[adata2.obs['celltype'] == c] for c in all_cell_types]
    sum = 0
    maxs = 0
    unmatched1 = [i for i in range(len(adatas1))] # indices in adatasX that are not matched
    unmatched2 = [j for j in range(len(adatas2))] # indices in adatasX that are not matched

    # get maximum weight bipartite graphs
    edge_dict = {}
    edges = []
    for i in unmatched1:
        if (len(adatas1[i]) == 0):
            continue
        for j in unmatched2:
            if (len(adatas2[j]) == 0):
                continue
            
            weight = algorithm([adatas1[i], adatas2[j]], normalize=normalize)
            edge_dict[(i,j)] = weight
            edges.append((i,j,weight))
    # compute maximum weight edges
    G = nx.Graph()
    G.add_weighted_edges_from(edges)
    non_matching_matches = nx.max_weight_matching(G)
    matches = []
    # compute score
    for (i,j) in non_matching_matches:
        if (i,j) not in edge_dict.keys():
            i, j = j, i # nx.max_weight_matching uses undirected edges
        sum += min(len(adatas1[i]), len(adatas2[j])) * edge_dict[(i,j)]
        maxs += max(len(adatas1[i]), len(adatas2[j]))
        matches.append(i)
        matches.append(j)

    # get remaining unmatched types, add to the count
    remaining1 = [i for i in unmatched1 if i not in matches]
    remaining2 = [j for j in unmatched2 if j not in matches]
    for r in remaining1:
        maxs += len(adatas1[r])
    for r in remaining2:
        maxs += len(adatas2[r])

    # compute 
    return sum/maxs