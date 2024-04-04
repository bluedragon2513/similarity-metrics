# imports
    # standard libaries
import anndata as ad
    # user libraries


# functions
def filter_batches(adata):
    adatas = []
    for c in adata.obs['tech'].cat.categories:
        new_adata = ad.AnnData(
            X=adata.layers['counts'][adata.obs['tech'] == c, :],
            obs=adata.obs[adata.obs['tech'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas.append(new_adata)
    return adatas

def filter_celltypes(adata):
    adatas = []
    for c in adata.obs['celltype'].cat.categories:
        new_adata = ad.AnnData(
            X=adata.layers['counts'][adata.obs['celltype'] == c, :],
            obs=adata.obs[adata.obs['celltype'] == c]
        )
        new_adata.layers['counts'] = new_adata.X
        adatas.append(new_adata)
    return adatas