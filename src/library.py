# imports
    # standard libraries
import numpy as np
    # user libraries
from scanorama import find_alignments

    # helper function for mwjmsi
def min_count(adata1, adata2, all_cell_types):
    len_a1 = len(adata1.layers['counts'])
    len_a2 = len(adata2.layers['counts'])

    mins = []
    for ct1 in all_cell_types:
        arr = []
        mins.append(arr)
        len_ct1 = len(adata1[adata1.obs['celltype'] == ct1])
        for ct2 in all_cell_types:
            len_ct2 = len(adata2[adata2.obs['celltype'] == ct2])
            arr.append(min(len_ct1/len_a1,len_ct2/len_a2))
    return mins

    # helper function for mwjmsi
def max_count(adata1, adata2, all_cell_types):
    count = np.sum([max(len(adata1[adata1.obs['celltype'] == ct]),len(adata2[adata2.obs['celltype'] == ct])) for ct in all_cell_types])
    return count

# def mwjmsi_v2(adata1, adata2, all_cell_types, similarity_matrix):
#     mins = min_count(adata1, adata2, all_cell_types)
#     sum = 0
#     for (min,sim) in zip(mins,similarity_matrix):
#         for (m,s) in zip(min,sim):
#             sum += m * s
#     return sum/max_count(adata1, adata2, all_cell_types)

# Modified Weighted-Jaccard-Multiset Similarity Index
def mwjmsi(adata1, adata2, all_cell_types):
    # print(adata1.layers['counts'])
    adatas1 = [adata1.layers['counts'][adata1.obs['celltype'] == c] for c in all_cell_types]
    adatas2 = [adata2.layers['counts'][adata2.obs['celltype'] == c] for c in all_cell_types]

    mins = min_count(adata1, adata2, all_cell_types)
    sum = 0
    for ct in range(len(all_cell_types)):
        for ct1 in range(len(all_cell_types)):
            if adatas1[ct].shape[0] == 0 or adatas2[ct1].shape[0] == 0 or ct == ct1:
                # print(f"no, not running {adatas1[ct].shape[0]} {adatas2[ct1].shape[0]} {ct} {ct1}")
                continue
            # print("yes, running")
            datas = [adatas1[ct],adatas2[ct1]]
            a, m, t = find_alignments(datas, verbose=1, knn=20)
            sum += mins[ct][ct1] * t[(0,1)]
            # print(t[(0,1)])
            # print(f"{a} {m} {len(adatas1[ct])} {len(adatas2[ct1])}")
            # print(a)
            # print()
    return sum/max_count(adata1, adata2, all_cell_types)