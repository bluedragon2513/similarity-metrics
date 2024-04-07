from scib import preprocessing
import anndata
import rpy2.robjects as ro
import scanpy
import pickle
import scanpy as sc


adata = sc.read_h5ad("scib/sim1_1_norm.h5ad")

preprocessing.save_seurat(
    adata=adata,
    path="seurat_object.rds",
    batch="sample_id",
    hvgs=["batchname", "batchname_all", "final_call_label"] 
)

