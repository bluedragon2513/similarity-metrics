from scib import preprocessing
import anndata as ad
import pickle
import scanpy as sc
import rpy2.robjects as robjects


def read_file(path):
    adata = sc.read_h5ad(path)
    return adata


def rds_filter_batches(adata, filter_by):
    adatas = {}
    for c in adata.obs[filter_by].cat.categories:
        new_adata = ad.AnnData(
            X=adata.layers['counts'][adata.obs[filter_by] == c, :],
            obs=adata.obs[adata.obs[filter_by] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas[c] = new_adata
    return adatas

def adata_to_seurat(adata, filter_by):
    adatas = rds_filter_batches(adata, filter_by)
    for k,adata in adatas.items():
        preprocessing.save_seurat(
            adata = adata,
            path = f"Seurat/data/{filter_by}/{k}.rds",
            batch=k,
        )

def run_seurat(path):
    robjects.r("source('src/Seurat.R')")
    function = robjects.r('run_seurat')
    r_string = robjects.StrVector([path])
    function(r_string)


