import anndata as ad
import numpy as np
from adata_preprocessing import scib_normalize
from experiment.modules.metrics import MetricsData

def kruskal(datasets, normalize=False, **kwargs):
    if normalize:
        adata = ad.AnnData(np.concatenate((datasets[0], datasets[1])))
        adata = scib_normalize(adata)
        datasets[0] = adata.X[:len(datasets[0])]
        datasets[1] = adata.X[len(datasets[0]):]
    
    M = MetricsData(np.transpose(datasets[0]), np.transpose(datasets[1]))
    return 1-M.compute_stress_kruskal()