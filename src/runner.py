# imports
    # standard
import scanpy as sc
from concurrent.futures import ProcessPoolExecutor
    # user
# from scorer_parallel import *
from scorer import *
from adata_preprocessing import *
from plotter import plotter
from library import experiment_amwjmsi
from scanorama import scanorama
from stress import kruskal

# constants
data = "data/preprocessed/human_pancreas_norm_complexBatch.h5ad"
algorithm = scanorama
dataset_name = "pancreas"
folder = "data/"
batch_save_file = "batch-scores"
celltype_save_file = "celltype-scores"
similarity_function = mwjmsi
dataset_normalize = None
batch_normalize = None
celltype_normalize = None

# functions
def run(data=data, 
        algorithm=algorithm, 
        dataset_name=dataset_name, 
        folder=folder, 
        batch_save_file=batch_save_file, 
        celltype_save_file=celltype_save_file,
        similarity_function=similarity_function,
        dataset_normalize=dataset_normalize,
        batch_normalize=batch_normalize,
        celltype_normalize=celltype_normalize):
    if dataset_normalize:
        adatas = filter_batches(dataset_normalize(raw_counts(sc.read_h5ad(data))))
    else:
        adatas = filter_batches(raw_counts(sc.read_h5ad(data)))
    

    scorer_batch(adatas, 
                 folder=folder, 
                 algorithm=algorithm, 
                 dataset_name=dataset_name, 
                 save_file=batch_save_file,
                 batch_normalize=batch_normalize)
    scorer_celltype(adatas, 
                    folder=folder, 
                    algorithm=algorithm, 
                    dataset_name=dataset_name, 
                    save_file=celltype_save_file, 
                    similarity_function=similarity_function,
                    batch_normalize=batch_normalize,
                    celltype_normalize=celltype_normalize)
    
    with open(f"{folder}{algorithm.__name__}/{dataset_name}/{batch_save_file}.pkl", "rb") as f:
        batch_scores = pickle.load(f)
    with open(f"{folder}{algorithm.__name__}/{dataset_name}/{batch_save_file}-jaccard.pkl", "rb") as f:
        batch_scores_jaccard = pickle.load(f)
    with open(f"{folder}{algorithm.__name__}/{dataset_name}/{celltype_save_file}-{similarity_function.__name__}.pkl", "rb") as f:
        celltype_scores = pickle.load(f)

    batch_scores = {k:v for k, v in batch_scores.items() if k[0] != k[1]}
    celltype_scores = {k:v for k, v in celltype_scores.items() if k[0] != k[1]}
    plotter(
        batch_scores, 
        celltype_scores, 
        f"{folder}{algorithm.__name__}/{dataset_name}/{similarity_function.__name__}", 
        title=f"{dataset_name.title()}: batch vs. celltype ({algorithm.__name__.title()})",
        annotations=batch_scores_jaccard
    )

    print("done")

def seurat_run():
    batch_scores = score_seurat('data/preprocessed/human_pancreas_preprocessed.h5ad')
    celltype_scores = score_celltypes_seurat('data/preprocessed/human_pancreas_preprocessed.h5ad')
    plotter(batch_scores, celltype_scores, "/Users/carsonnannini/Research/similarity-metrics/data/seurat", "Batch vs Celltype (Seurat)")

# main
if __name__ == "__main__":
    run(algorithm=scanorama, dataset_normalize=scib_normalize)
    run(algorithm=scanorama, dataset_normalize=scib_normalize, similarity_function=amwjmsi)
    # run(algorithm=scanorama, batch_normalize=scib_normalize, celltype_normalize=None)
    # run(algorithm=scanorama, batch_normalize=scib_normalize, celltype_normalize=None, similarity_function=amwjmsi)
    # run(algorithm=scanorama, batch_normalize=None, celltype_normalize=scib_normalize)
    # run(algorithm=scanorama, batch_normalize=None, celltype_normalize=scib_normalize, similarity_function=amwjmsi)
    # run(algorithm=scanorama, data="data/preprocessed/human_pancreas_norm_complexBatch.h5ad")
    # run(algorithm=scanorama, data="data/preprocessed/human_pancreas_norm_complexBatch.h5ad", similarity_function=amwjmsi)

    print("FINISHED!")

    # run(algorithm=scanorama, similarity_function=amwjmsi, batch_normalize=scib_normalize)
    # run(algorithm=kruskal, similarity_function=experiment_amwjmsi, normalize=False)
