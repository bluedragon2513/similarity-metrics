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
data_folder = "data"
batch_save_file = "batch-scores"
celltype_save_file = "celltype-scores"
combiner = mwjmsi
dataset_normalize = None
batch_normalize = None
celltype_normalize = None

# functions
def run(data=data, 
        algorithm=algorithm, 
        dataset_name=dataset_name, 
        data_folder=data_folder, 
        batch_save_file=batch_save_file, 
        celltype_save_file=celltype_save_file,
        combiner=combiner,
        dataset_normalize=dataset_normalize,
        batch_normalize=batch_normalize,
        celltype_normalize=celltype_normalize,
        **kwargs):
    
    # normalization is mutually exclusive exclusive
    assert(sum([True if n else False for n in [dataset_normalize, batch_normalize, celltype_normalize]]) <= 1)
    if dataset_normalize:
        processing = "dataset_normalize"
    elif batch_normalize:
        processing = "batch_normalize"
    elif celltype_normalize:
        processing = "celltype_normalize"
    else:
        processing = "none_normalize"
    
    print("Getting the data")
    if dataset_normalize:
        adatas = filter_batches(dataset_normalize(raw_counts(sc.read_h5ad(data))))
    else:
        adatas = filter_batches(raw_counts(sc.read_h5ad(data)))

    print("Scoring Batches")
    # scorer_batch(adatas, 
    #              data_folder=data_folder, 
    #              algorithm=algorithm, 
    #              dataset_name=dataset_name, 
    #              save_file=batch_save_file,
    #              batch_normalize=batch_normalize,
    #              processing=processing,
    #              **kwargs)
    print("Scoring Cell Types")
    scorer_celltype(adatas, 
                    data_folder=data_folder, 
                    algorithm=algorithm, 
                    dataset_name=dataset_name, 
                    save_file=celltype_save_file, 
                    combiner=combiner,
                    batch_normalize=batch_normalize,
                    celltype_normalize=celltype_normalize,
                    processing=processing,
                    **kwargs)
    
    print("Plotting scores")
    folder = f"{data_folder}/{dataset_name}/{algorithm.__name__}/{processing}/"
    with open(f"{folder}{batch_save_file}.pkl", "rb") as f:
        batch_scores = pickle.load(f)
    with open(f"{folder}{batch_save_file}-jaccard.pkl", "rb") as f:
        batch_scores_jaccard = pickle.load(f)
    with open(f"{folder}{celltype_save_file}-{combiner.__name__}.pkl", "rb") as f:
        celltype_scores = pickle.load(f)

    batch_scores = {k:v for k, v in batch_scores.items() if k[0] != k[1]}
    celltype_scores = {k:v for k, v in celltype_scores.items() if k[0] != k[1]}
    plotter(
        batch_scores, 
        celltype_scores, 
        f"{folder}/{combiner.__name__}", 
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
    # run(algorithm=scanorama, dataset_normalize=scib_normalize, verbose=False)
    # run(algorithm=scanorama, dataset_normalize=scib_normalize, combiner=amwjmsi, verbose=False)
    run(algorithm=scanorama, batch_normalize=scib_normalize, celltype_normalize=None, verbose=False)
    run(algorithm=scanorama, batch_normalize=scib_normalize, celltype_normalize=None, similarity_function=amwjmsi, verbose=False)
    run(algorithm=scanorama, batch_normalize=None, celltype_normalize=scib_normalize, verbose=False)
    run(algorithm=scanorama, batch_normalize=None, celltype_normalize=scib_normalize, similarity_function=amwjmsi, verbose=False)
    run(algorithm=scanorama, data="data/preprocessed/human_pancreas_norm_complexBatch.h5ad", verbose=False)
    run(algorithm=scanorama, data="data/preprocessed/human_pancreas_norm_complexBatch.h5ad", similarity_function=amwjmsi, verbose=False)

    run(algorithm=kruskal, dataset_normalize=scib_normalize, verbose=False)
    run(algorithm=kruskal, dataset_normalize=scib_normalize, combiner=amwjmsi, verbose=False)
    run(algorithm=kruskal, batch_normalize=scib_normalize, celltype_normalize=None, verbose=False)
    run(algorithm=kruskal, batch_normalize=scib_normalize, celltype_normalize=None, similarity_function=amwjmsi, verbose=False)
    run(algorithm=kruskal, batch_normalize=None, celltype_normalize=scib_normalize, verbose=False)
    run(algorithm=kruskal, batch_normalize=None, celltype_normalize=scib_normalize, similarity_function=amwjmsi, verbose=False)
    run(algorithm=kruskal, data="data/preprocessed/human_pancreas_norm_complexBatch.h5ad", verbose=False)
    run(algorithm=kruskal, data="data/preprocessed/human_pancreas_norm_complexBatch.h5ad", similarity_function=amwjmsi, verbose=False)

    print("FINISHED!")

    # run(algorithm=scanorama, similarity_function=amwjmsi, batch_normalize=scib_normalize)
    # run(algorithm=kruskal, similarity_function=experiment_amwjmsi, normalize=False)
