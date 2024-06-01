# imports
    # standard libraries
import scanpy as sc
import sys, os
sys.path.append(os.getcwd())
    # user libraries
from src.library.scorer import *
from src.library.adata_preprocessing import *
from src.library.plotter import plotter
from src.library.json_handler import *
from src.metrics.Scanorama.scanorama import scanorama
from src.metrics.Stress.stress import kruskal
from src.metrics.Seurat.seurat import seurat

# constants
DATA = "data/pancreas/human_pancreas_norm_complexBatch.h5ad"
ALGORITHM = scanorama
DATASET_NAME = "pancreas"
DATA_FOLDER = "data"
BATCH_SAVE_FILE = "batch-scores"
CELLTYPE_SAVE_FILE = "celltype-scores"
COMBINER = mwjmsi
DATASET_NORMALIZE = None
BATCH_NORMALIZE = None
CELLTYPE_NORMALIZE = None


# functions
def run(data: str=DATA, 
        algorithm: Callable=ALGORITHM, 
        dataset_name: str=DATASET_NAME, 
        data_folder: str=DATA_FOLDER, 
        batch_save_file: str=BATCH_SAVE_FILE, 
        celltype_save_file: str=CELLTYPE_SAVE_FILE,
        combiner: Callable=COMBINER,
        dataset_normalize: Callable=DATASET_NORMALIZE,
        batch_normalize: Callable=BATCH_NORMALIZE,
        celltype_normalize: Callable=CELLTYPE_NORMALIZE,
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

    score(adatas, 
          data_folder=data_folder, 
          algorithm=algorithm, 
          dataset_name=dataset_name, 
          batch_save_file=batch_save_file,
          celltype_save_file=celltype_save_file, 
          combiner=combiner,
          batch_normalize=batch_normalize,
          celltype_normalize=celltype_normalize,
          processing=processing,
          **kwargs)
    
    print("Plotting scores")
    print("\tFetching saved data")
    folder = f"{data_folder}/{dataset_name}/{algorithm.__name__}/{processing}/"
    batch_scores = json_reader(f"{folder}{batch_save_file}.json")
    batch_scores_jaccard = json_reader(f"{folder}gjsi-scores.json")
    celltype_scores = json_reader(f"{folder}{celltype_save_file}-{combiner.__name__}.json")

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
    run(algorithm=seurat, dataset_normalize=scib_normalize, verbose=False)
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
    pass
