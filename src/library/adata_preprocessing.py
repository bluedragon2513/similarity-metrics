# imports
    # standard libaries
from anndata import AnnData
import scanpy as sc
from sklearn.preprocessing import normalize
from typing import List, Dict
    # user libraries


# functions
def rds_filter_batches(adata: AnnData) -> Dict[str, AnnData]:
    adatas = {}
    for c in adata.obs['tech'].cat.categories:
        new_adata = AnnData(
            X=adata.layers['counts'][adata.obs['tech'] == c, :],
            obs=adata.obs[adata.obs['tech'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas[c] = new_adata
    return adatas

def filter_batches(adata: AnnData) -> List[AnnData]:
    adatas = []
    for c in adata.obs['tech'].cat.categories:
        new_adata = AnnData(
            X=adata.X[adata.obs['tech'] == c, :],
            obs=adata.obs[adata.obs['tech'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas.append(new_adata)
    return adatas

def filter_celltypes(adata: AnnData) -> List[AnnData]:
    adatas = []
    for c in adata.obs['celltype'].cat.categories:
        new_adata = AnnData(
            X=adata.X[adata.obs['celltype'] == c, :],
            obs=adata.obs[adata.obs['celltype'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas.append(new_adata)
    return adatas

def filter_hvg(adata: AnnData) -> AnnData:
    if (adata.X.shape[0] <= 2000):
        return adata
    sc.pp.filter_genes(adata, min_cells=5)
    hvg_filter = sc.pp.highly_variable_genes(adata, n_top_genes=2000, n_bins=20, flavor="cell_ranger", inplace=False)['highly_variable']
    anndata = AnnData(
        X=adata.X[:,hvg_filter],
        obs=adata.obs.copy(),
    )
    anndata.layers['counts'] = adata.X[:,hvg_filter]
    return anndata

def scib_normalize(adata: AnnData) -> AnnData:
    adata = sc.pp.normalize_total(adata, target_sum=1e6, copy=True) # each cell has the same total count after normalization
    sc.pp.log1p(adata) # so that every value is applied log(1+X)
    return filter_hvg(adata)

def scan_normalize(adata: AnnData, axis=1) -> None:
    normalize(adata.X, axis=axis)

def raw_counts(adata: AnnData) -> AnnData:
    anndata = AnnData(
        X=adata.layers['counts'],
        obs=adata.obs.copy(),
    )
    anndata.layers['counts'] = adata.layers['counts']
    return anndata