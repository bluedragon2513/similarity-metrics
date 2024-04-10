# imports
    # standard libraries
import numpy as np
import networkx as nx
    # user libraries
from scanorama import find_alignments

    # helper function for mwjmsi
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

# Modified Weighted-Jaccard-Multiset Similarity Index
def mwjmsi(adata1, adata2, all_cell_types):
    adatas1 = [adata1.X[adata1.obs['celltype'] == c] for c in all_cell_types]
    adatas2 = [adata2.X[adata2.obs['celltype'] == c] for c in all_cell_types]

    mins = min_count(adata1, adata2, all_cell_types)
    sum = 0
    for ct in range(len(all_cell_types)):
        # if there are no cells of either adata, then ignore it
        if adatas1[ct].shape[0] == 0 or adatas2[ct].shape[0] == 0:
            continue
        datas = [adatas1[ct], adatas2[ct]]
        _, _, t = find_alignments(datas, verbose=1, knn=20)
        sum += mins[ct][ct] * t[(0,1)]

    return sum/max_count(adata1, adata2, all_cell_types)

# Adjusted Modified Weighted-Jaccard-Multiset Similarity Index
def amwjmsi(adata1, adata2, all_cell_types):
    adatas1 = [adata1.X[adata1.obs['celltype'] == c] for c in all_cell_types]
    adatas2 = [adata2.X[adata2.obs['celltype'] == c] for c in all_cell_types]
    # print(adatas1); print(adatas2)

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
        _, _, t = find_alignments(datas, verbose=1, knn=20)
        # print(f"\t{mins[ct][ct]}")
        sum += mins[ct][ct] * t[(0,1)]
        maxs += max(len(adatas1[ct]), len(adatas2[ct]))

    # get maximum weight bipartite graphs
    edge_dict = {}
    edges = []
    for i in unmatched1:
        for j in unmatched2:
            _, _, t = find_alignments([adatas1[i], adatas2[j]], verbose=1, knn=20)
            weight = t[(0,1)]
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
    # print(f"sum: {sum} - maxs: {maxs} - score: {sum/maxs}")
    return sum/maxs

# Modified Weighted-Jaccard-Multiset Similarity Index
def experiment_amwjmsi(adata1, adata2, all_cell_types):
    adatas1 = [adata1.X[adata1.obs['celltype'] == c] for c in all_cell_types]
    adatas2 = [adata2.X[adata2.obs['celltype'] == c] for c in all_cell_types]
    # print(adatas1); print(adatas2)

    mins = min_count(adata1, adata2, all_cell_types)
    sum = 0
    maxs = 0
    unmatched1 = [i for i in range(len(adatas1))] # indices in adatasX that are not matched
    unmatched2 = [j for j in range(len(adatas2))] # indices in adatasX that are not matched
    # compute matched indices
    # for ct in range(len(all_cell_types)):
    #     if adatas1[ct].shape[0] == 0:
    #         if adatas2[ct].shape[0] != 0:
    #             unmatched2.append(ct)
    #         continue
    #     if adatas2[ct].shape[0] == 0:
    #         unmatched1.append(ct)
    #         continue
    #     datas = [adatas1[ct], adatas2[ct]]
    #     _, _, t = find_alignments(datas, verbose=1, knn=20)
    #     # print(f"\t{mins[ct][ct]}")
    #     sum += mins[ct][ct] * t[(0,1)]
    #     maxs += max(len(adatas1[ct]), len(adatas2[ct]))

    # get maximum weight bipartite graphs
    edge_dict = {}
    edges = []
    for i in unmatched1:
        if (len(adatas1[i]) == 0):
            continue
        for j in unmatched2:
            if (len(adatas2[j]) == 0):
                continue
            _, _, t = find_alignments([adatas1[i], adatas2[j]], verbose=1, knn=20)
            weight = t[(0,1)]
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
    # print(f"sum: {sum} - maxs: {maxs} - score: {sum/maxs}")
    return sum/maxs