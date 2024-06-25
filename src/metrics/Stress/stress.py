import anndata as ad
import numpy as np
from sklearn.decomposition import PCA
from src.library.adata_preprocessing import scib_normalize
from src.metrics.Stress.metrics import MetricsData


def kruskal(datasets, normalize=False, **kwargs):
    if normalize:
        adata = ad.AnnData(np.concatenate((datasets[0], datasets[1])))
        adata = scib_normalize(adata)
        datasets[0] = adata.X[:len(datasets[0])]
        datasets[1] = adata.X[len(datasets[0]):]
    
    M = MetricsData(np.transpose(datasets[0]), np.transpose(datasets[1]))
    return 1-M.compute_stress_kruskal()

def kruskal_with_pca(datasets, normalize=False, **kwargs):
    if normalize:
        adata = ad.AnnData(np.concatenate((datasets[0], datasets[1])))
        adata = scib_normalize(adata)
        datasets[0] = adata.X[:len(datasets[0])]
        datasets[1] = adata.X[len(datasets[0]):]

    d1, d2 = np.transpose(datasets[0]), np.transpose(datasets[1])
    num_components = min(100, d1.shape[0], d2.shape[0])
    pca = PCA(n_components=num_components)
    pca.fit(d1)
    d1 = pca.components_
    pca.fit(d2)
    d2 = pca.components_
    
    M = MetricsData(d1, d2)
    return 1-M.compute_stress_kruskal()