# imports
    # standard
import scanpy as sc
import sys
import os
sys.path.append(os.getcwd())
    # user
# from scorer_parallel import *
from src.library.scorer import *
from src.library.adata_preprocessing import *
from src.library.plotter import plotter
from src.metrics.Scanorama.scanorama import scanorama
from src.metrics.Stress.stress import kruskal

# constants
data = "data/pancreas/human_pancreas_norm_complexBatch.h5ad"
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
def run(data: str=data, 
        algorithm: Callable=algorithm, 
        dataset_name: str=dataset_name, 
        data_folder: str=data_folder, 
        batch_save_file: str=batch_save_file, 
        celltype_save_file: str=celltype_save_file,
        combiner: Callable=combiner,
        dataset_normalize: Callable=dataset_normalize,
        batch_normalize: Callable=batch_normalize,
        celltype_normalize: Callable=celltype_normalize,
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
    scorer_batch(adatas, 
                 data_folder=data_folder, 
                 algorithm=algorithm, 
                 dataset_name=dataset_name, 
                 save_file=batch_save_file,
                 batch_normalize=batch_normalize,
                 processing=processing,
                 **kwargs)
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
    print("\tFetching saved data")
    folder = f"{data_folder}/{dataset_name}/{algorithm.__name__}/{processing}/"
    with open(f"{folder}{batch_save_file}.pkl", "rb") as f:
        batch_scores = pickle.load(f)
    with open(f"{folder}{batch_save_file}-jaccard.pkl", "rb") as f:
        batch_scores_jaccard = pickle.load(f)
    with open(f"{folder}{celltype_save_file}-{combiner.__name__}.pkl", "rb") as f:
        celltype_scores = pickle.load(f)

    batch_scores = {k:v for k, v in batch_scores.items() if k[0] != k[1]}
    celltype_scores = {k:v for k, v in celltype_scores.items() if k[0] != k[1]}
    print("\tPlotting")
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
    print("1")
    # run(algorithm=scanorama, dataset_normalize=scib_normalize, verbose=False)
    # run(algorithm=scanorama, dataset_normalize=scib_normalize, combiner=amwjmsi, verbose=False)
    # print("2")
    # run(algorithm=scanorama, batch_normalize=scib_normalize, celltype_normalize=None, verbose=False)
    # run(algorithm=scanorama, batch_normalize=scib_normalize, celltype_normalize=None, combiner=amwjmsi, verbose=False)
    # print("3")
    # run(algorithm=scanorama, batch_normalize=None, celltype_normalize=scib_normalize, verbose=False)
    # print("3.5")
    # run(algorithm=scanorama, batch_normalize=None, celltype_normalize=scib_normalize, combiner=amwjmsi, verbose=False)
    # print("4")
    # run(algorithm=scanorama, data="data/preprocessed/human_pancreas_norm_complexBatch.h5ad", verbose=False)
    # run(algorithm=scanorama, data="data/preprocessed/human_pancreas_norm_complexBatch.h5ad", combiner=amwjmsi, verbose=False)

    # print("1")
    # run(algorithm=kruskal, dataset_normalize=scib_normalize, verbose=False)
    # print("2")
    # run(algorithm=kruskal, dataset_normalize=scib_normalize, combiner=amwjmsi, verbose=False)
    # print("3")
    # run(algorithm=kruskal, batch_normalize=scib_normalize, celltype_normalize=None, verbose=False)
    # print("4")
    # run(algorithm=kruskal, batch_normalize=scib_normalize, celltype_normalize=None, combiner=amwjmsi, verbose=False)
    # print("5")
    # run(algorithm=kruskal, batch_normalize=None, celltype_normalize=scib_normalize, verbose=False)
    # print("6")
    # run(algorithm=kruskal, batch_normalize=None, celltype_normalize=scib_normalize, combiner=amwjmsi, verbose=False)
    # print("7")
    # run(algorithm=kruskal, verbose=False)
    # print("8")
    # run(algorithm=kruskal, combiner=amwjmsi, verbose=False)

    print("FINISHED!")

    # run(algorithm=scanorama, combiner=amwjmsi, batch_normalize=scib_normalize)
    # run(algorithm=kruskal, combiner=experiment_amwjmsi, normalize=False)
    pass
