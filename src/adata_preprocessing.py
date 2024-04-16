# imports
    # standard libaries
import anndata as ad
import scanpy as sc
    # user libraries


# functions
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

def filter_batches(adata):
    adatas = []
    for c in adata.obs['tech'].cat.categories:
        new_adata = ad.AnnData(
            X=adata.X[adata.obs['tech'] == c, :],
            obs=adata.obs[adata.obs['tech'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas.append(new_adata)
    return adatas

def filter_celltypes(adata):
    adatas = []
    for c in adata.obs['celltype'].cat.categories:
        new_adata = ad.AnnData(
            X=adata.X[adata.obs['celltype'] == c, :],
            obs=adata.obs[adata.obs['celltype'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas.append(new_adata)
    return adatas

def filter_hvg(adata):
    hvg_filter = sc.pp.highly_variable_genes(adata, layer='counts', n_top_genes=2000, n_bins=20, flavor="cell_ranger", inplace=False)['highly_variable']
    anndata = ad.AnnData(
        X=adata.X[:,hvg_filter],
        obs=adata.obs.copy(),
    )
    anndata.layers['counts'] = adata.layers['counts'][:,hvg_filter]
    return anndata

def scib_normalize(adata):
    sc.pp.normalize_total(adata, target_sum=1e6) # each cell has the same total count after normalization
    sc.pp.log1p(adata) # so that every value is applied log(1+X)

# main method
if __name__ == "__main__":
    adata = sc.read_h5ad("data/preprocessed/human_pancreas_norm_complexBatch.h5ad")
    hvg_filter = sc.pp.highly_variable_genes(adata, layer='counts', n_top_genes=2000, n_bins=20, flavor="cell_ranger", inplace=False)['highly_variable']
    # adata.obs.reset_index(drop=True, inplace=True)  # Reset index to ensure alignment
    anndata = ad.AnnData(
        X=adata.X[:,hvg_filter],
        obs=adata.obs.copy(),
    )
    anndata.layers['counts'] = adata.layers['counts'][:,hvg_filter]
    sc.write("data/preprocessed/human_pancreas_preprocessed.h5ad", anndata)