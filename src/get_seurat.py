from scib import preprocessing
import anndata as ad
import pickle
import scanpy as sc
import rpy2.robjects as robjects


def rds_filter_batches(adata):
    adatas = {}
    for c in adata.obs['tech'].cat.categories:
        new_adata = ad.AnnData(
            X=adata.layers['counts'][adata.obs['tech'] == c, :],
            obs=adata.obs[adata.obs['tech'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas[c] = new_adata
    return adatas


# def filter_batches(adata):
#     adatas = []
#     for c in adata.obs['tech'].cat.categories:
#         new_adata = ad.AnnData(
#             X=adata.layers['counts'][adata.obs['tech'] == c, :],
#             obs=adata.obs[adata.obs['tech'] == c]
#         )
#         new_adata.layers['counts'] = new_adata.X
#         adatas.append(new_adata)
#     return adatas


adata = sc.read_h5ad("/Users/carsonnannini/Research/similarity-metrics/data/preprocessed/human_pancreas_preprocessed.h5ad")


adatas = rds_filter_batches(adata)
for k,adata in adatas.items():
    preprocessing.save_seurat(
        adata = adata,
        path = f"Seurat/data/{k}.rds",
        batch=k,
    )
