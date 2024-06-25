import anndata as ad
from anndata import AnnData
from typing import List
import numpy as np
from numpy.linalg import norm
import networkx as nx
from sklearn.decomposition import PCA
from src.library.adata_preprocessing import scib_normalize
from src.library.combiner import count, max_count
from src.metrics.Stress.metrics import MetricsData


def get_cell_types(adata1, adata2):
    return np.unique(adata1.obs['celltype'].copy().extend(adata2.obs['celltype']))

def get_matching_cell_types(adata1, adata2, all_cell_types=None):
    if not all_cell_types:
        all_cell_types = get_cell_types(adata1, adata2)
    matching_cell_types = []

    for ct in all_cell_types:
        if count(adata1, ct) != 0 and count(adata2, ct) != 0:
            matching_cell_types.append(ct)

    return matching_cell_types

def get_non_matching_cell_types(adata1, adata2, all_cell_types=None):
    if not all_cell_types:
        all_cell_types = get_cell_types(adata1, adata2)
    non_matching_cell_types = []

    for ct in all_cell_types:
        if count(adata1, ct) == 0 or count(adata2, ct) == 0:
            non_matching_cell_types.append(ct)

    return non_matching_cell_types

def not_in(set1, set2):
    """
        Returns items in set1 that are not in set2
    """
    diff = []
    for item in set1:
        if item not in set2:
            diff.append(item)
    return diff

def get_indices_of_celltypes(adata, selected_cell_types):
    return [(ct in selected_cell_types) for ct in adata.obs['celltype']]

def get_unmatched_cells(adata1, adata2, matching_cell_types=None):
    if not matching_cell_types:
        matching_cell_types = get_matching_cell_types(adata1, adata2)

    # getting the cell types from the batches
    cell_types_1 = np.unique(adata1.obs['celltype'])
    cell_types_2 = np.unique(adata2.obs['celltype'])

    # getting the cell types that are not in the matching cell types
    cell_types_1 = not_in(cell_types_1, matching_cell_types)
    cell_types_2 = not_in(cell_types_2, matching_cell_types)

    # getting the cells that are not in the matching cell types
    b1_cells = adata1.X[get_indices_of_celltypes(adata1, cell_types_1)]
    b2_cells = adata2.X[get_indices_of_celltypes(adata2, cell_types_2)]

    return b1_cells, b2_cells

def compute_cosine_similarity(cell, cells2):
    X, Y = np.array(cells2), np.array(cell)
    return np.dot(X,Y)/(norm(X, axis=1)*norm(Y))

def compute_cosine_similarities(cells1, cells2):
    matrix = []
    for cell in cells1:
        matrix.append(compute_cosine_similarity(cell, cells2))
    return matrix

# currently this is not working because we are assuming combiner.py never passes along observation data to the algorithms
# this was working off of the assumption that the algorithms would never require those observations
def jaccard_on_matching_and_bipartite_cosine_similarity(datasets, normalize=False, **kwargs):
    if normalize:
        adata = ad.AnnData(np.concatenate((datasets[0], datasets[1])))
        adata = scib_normalize(adata)
        datasets[0] = adata.X[:len(datasets[0])]
        datasets[1] = adata.X[len(datasets[0]):]

    adata1, adata2 = datasets[0], datasets[1]
    all_cell_types = get_cell_types(adata1, adata2)
    matching_cell_types = get_matching_cell_types(adata1, adata2, all_cell_types=all_cell_types)
    non_matching_cell_types = get_non_matching_cell_types(adata1, adata2, all_cell_types=all_cell_types)

    # compute the scores for the matching cell types
    numerator = sum([min(count(adata1,ct),count(adata2,ct)) for ct in all_cell_types])
    denominator = max_count(adata1, adata2, matching_cell_types)

    # compute the cosine similarity matrix
    unmatched_1, unmatched_2 = get_unmatched_cells(adata1, adata2, matching_cell_types)
    cosine_similarity_matrix = compute_cosine_similarities(unmatched_1, unmatched_2)
    
    # compute the bipartite using the cosine similarity matrix
    graph = nx.Graph()
    edges = []

    offset = len(cosine_similarity_matrix)
    for i in range(len(cosine_similarity_matrix)):
        for j in range(len(cosine_similarity_matrix[i])):
            edges.append((i, j+offset, cosine_similarity_matrix[i][j]))

    graph.add_weighted_edges_from(edges)
    maximum_weight_graph = nx.max_weight_matching(graph)

    # add maximum weights to numerator, add min cells to denominator
    for (i,j) in maximum_weight_graph:
        if j < offset:
            i, j = j, i
        j -= offset

        numerator += cosine_similarity_matrix[i,j]
    denominator += min(len(unmatched_1), len(unmatched_2))

    return numerator/denominator

def bipartite_cosine_similarity(datasets, normalize=False, **kwargs):
    if normalize:
        adata = ad.AnnData(np.concatenate((datasets[0], datasets[1])))
        adata = scib_normalize(adata)
        datasets[0] = adata.X[:len(datasets[0])]
        datasets[1] = adata.X[len(datasets[0]):]
    
    total_weight = 0
    min_cells = min(datasets[0].X.shape[0], datasets[1].X.shape[0])

    # compute cosine similarity matrix
    cosine_similarity_matrix = compute_cosine_similarities(datasets[0], datasets[1])

    # compute the bipartite using the cosine similarity matrix
    graph = nx.Graph()
    edges = []

    offset = len(cosine_similarity_matrix)
    for i in range(len(cosine_similarity_matrix)):
        for j in range(len(cosine_similarity_matrix[i])):
            edges.append((i, j+offset, cosine_similarity_matrix[i][j]))

    graph.add_weighted_edges_from(edges)
    maximum_weight_graph = nx.max_weight_matching(graph)

    # add maximum weights to total_weight
    for (i,j) in maximum_weight_graph:
        if j < offset:
            i, j = j, i
        j -= offset

        total_weight += cosine_similarity_matrix[i,j]

    return total_weight/min_cells

# def cosine_similarity(datasets, normalize=False, **kwargs):
#     if normalize:
#         adata = ad.AnnData(np.concatenate((datasets[0], datasets[1])))
#         adata = scib_normalize(adata)
#         datasets[0] = adata.X[:len(datasets[0])]
#         datasets[1] = adata.X[len(datasets[0]):]

    

def kruskal(datasets, normalize=False, **kwargs):
    if normalize:
        adata = ad.AnnData(np.concatenate((datasets[0], datasets[1])))
        adata = scib_normalize(adata)
        datasets[0] = adata.X[:len(datasets[0])]
        datasets[1] = adata.X[len(datasets[0]):]
    
    M = MetricsData(np.transpose(datasets[0]), np.transpose(datasets[1]))
    return 1-M.compute_stress_kruskal()

def get_pca_components(matrix, num_components=100):
    pca = PCA(n_components=num_components)
    pca.fit(matrix)
    return pca.components_

def kruskal_with_pca(datasets, normalize=False, **kwargs):
    if normalize:
        adata = ad.AnnData(np.concatenate((datasets[0], datasets[1])))
        adata = scib_normalize(adata)
        datasets[0] = adata.X[:len(datasets[0])]
        datasets[1] = adata.X[len(datasets[0]):]

    d1, d2 = np.transpose(datasets[0]), np.transpose(datasets[1])
    num_components = min(100, d1.shape[0], d2.shape[0])
    d1, d2 = get_pca_components(d1, num_components), get_pca_components(d2, num_components)
    
    M = MetricsData(d1, d2)
    return 1-M.compute_stress_kruskal()