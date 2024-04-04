# imports
    # standard
import scanpy as sc
from concurrent.futures import ProcessPoolExecutor
    # user
from scorer_parallel import *
from adata_preprocessing import *
# functions

# main
if __name__ == "__main__":
    adata = sc.read_h5ad("src/data/human_pancreas_norm_complexBatch.h5ad")
    adatas = filter_batches(adata)

    # with ProcessPoolExecutor() as executor:
    #     # executor.submit(scorer_batch, adatas)
    #     executor.submit(scorer_celltype, adatas)
    #     executor.shutdown(wait=True)
    scorer_celltype(adatas)
    print("done")