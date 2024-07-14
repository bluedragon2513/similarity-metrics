# imports
    # standard libraries
from annoy import AnnoyIndex
import numpy as np
from scipy.sparse import vstack, issparse
import anndata as ad
import networkx as nx
    # user libraries
# from src.library.adata_preprocessing import scib_normalize

# Default parameters
KNN = 20

def lcc(datasets, normalize=None, **kwargs):
    datasets = list(map(lambda dataset: dataset.toarray() if issparse(dataset) else dataset, datasets))
    if normalize:
        adata = ad.AnnData(np.concatenate((datasets[0], datasets[1])))
        adata = normalize(adata)
        datasets[0] = adata.X[:len(datasets[0])]
        datasets[1] = adata.X[len(datasets[0]):]

    num_lcc = len(largest_connected_components(datasets[0], datasets[1], knn=kwargs.get("knn", KNN)))
    num_points = datasets[0].shape[0] + datasets[1].shape[0]
    return num_lcc/num_points

def largest_connected_components(d1, d2, knn=KNN):
    G = nx.Graph()
    for n in d1:
        G.add_node(tuple(n))
    for n in d2:
        G.add_node(tuple(n))

    # creating the edges based on cosine similarity
    for (i, j) in nn_approx(d1, d2, knn=knn, metric='angular'):
        G.add_edge(tuple(d1[i]), tuple(d2[j]))

    for (j, i) in nn_approx(d2, d1, knn=knn, metric='angular'):
        G.add_edge(tuple(d2[j]), tuple(d1[i]))

    # for node in d1:
    #     approx_nns = nn_approx(np.array([node]), d2, knn=knn, metric='angular')
    #     for nn in approx_nns:
    #         G.add_edge(tuple(node), tuple(d2[nn]))
    # for node in d2:
    #     approx_nns = nn_approx(np.array([node]), d1, knn=knn, metric='angular')
    #     for nn in approx_nns:
    #         G.add_edge(tuple(node), tuple(d1[nn]))

    return max(nx.connected_components(G), key=len)

# taken from Scanorama
def nn_approx(ds1, ds2, knn=KNN, metric='manhattan', n_trees=10):
    # Build index.
    a = AnnoyIndex(ds2.shape[1], metric=metric)
    for i in range(ds2.shape[0]):
        a.add_item(i, ds2[i, :])
    a.build(n_trees)

    # Search index.
    ind = []
    for i in range(ds1.shape[0]):
        ind.append(a.get_nns_by_vector(ds1[i, :], knn, search_k=-1))
    ind = np.array(ind)

    # Match.
    match = set()
    for a, b in zip(range(ds1.shape[0]), ind):
        for b_i in b:
            match.add((a, b_i))

    return match


if __name__ == "__main__":
    d1 = [
        [1,2,3,4,5],
        [2,3,4,5,6],
        [-1,-2,-3,-4,-5],
        [500,600,700,800,900],
        [3,3,3,3,3],
        [10000,20000,30000,40000,50000]
    ]
    d2 = [
        [400, 500, 300, 200, 100],
        [-50000,-40000,-30000,-20000,-10000]
    ]
    d1, d2 = np.array(d1), np.array(d2)

    print(len(largest_connected_components(d1, d2))/(d1.shape[0] + d2.shape[0]))