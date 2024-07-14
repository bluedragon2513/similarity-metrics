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
from src.metrics.LCC.lcc import lcc
from src.metrics.Stress.stress import kruskal
from src.metrics.Seurat.seurat import seurat

# constants
DATA = "data/lung/Lung_atlas_public.h5ad"
ALGORITHM = scanorama
DATASET_NAME = "lung"
DATA_FOLDER = "data"
BATCH_SAVE_FILE = "batch-scores"
CELLTYPE_SAVE_FILE = "celltype-scores"
COMBINER = mwjmsi
DATASET_NORMALIZE = None
BATCH_NORMALIZE = None
CELLTYPE_NORMALIZE = None
RERUN = True


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
        **kwargs) -> None:
    
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
          rerun=RERUN,
          **kwargs)
    
    plot(
        algorithm=algorithm,
        dataset_name=dataset_name,
        data_folder=data_folder,
        batch_save_file=batch_save_file,
        celltype_save_file=celltype_save_file,
        combiner=combiner,
        processing=processing,
        **kwargs
    )

def plot(algorithm: Callable=ALGORITHM, 
        dataset_name: str=DATASET_NAME, 
        data_folder: str=DATA_FOLDER, 
        batch_save_file: str=BATCH_SAVE_FILE, 
        celltype_save_file: str=CELLTYPE_SAVE_FILE,
        combiner: Callable=COMBINER,
        processing: str="",
        **kwargs) -> None:
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
        d1=batch_scores, 
        d2=celltype_scores, 
        save_file=f"{folder}/{combiner.__name__}-new", 
        title=f"{dataset_name.title()}: batch vs. celltype ({algorithm.__name__.title()})",
        annotations=batch_scores_jaccard,
        show=False,
        random_colors=False,
        b1_colors=True
    )

    print("done")

# main
if __name__ == "__main__":
    alg = lcc
    # print("1 - dataset normalization")
    run(algorithm=seurat, dataset_normalize=scib_normalize, verbose=False)
    run(algorithm=seurat, dataset_normalize=scib_normalize, combiner=amwjmsi, verbose=False)
    # print("2 - batch normalization")
    run(algorithm=seurat, batch_normalize=scib_normalize, celltype_normalize=None, verbose=False)
    run(algorithm=seurat, batch_normalize=scib_normalize, celltype_normalize=None, combiner=amwjmsi, verbose=False)
    # print("3 - cell type normalization")
    # run(algorithm=alg, batch_normalize=None, celltype_normalize=scib_normalize, verbose=False)
    # run(algorithm=alg, batch_normalize=None, celltype_normalize=scib_normalize, combiner=amwjmsi, verbose=False)
    # print("4 - no normalization")
    # run(algorithm=alg, verbose=False)
    # run(algorithm=alg, combiner=amwjmsi, verbose=False)

    # alg2 = kruskal
    # print("1")
    # run(algorithm=alg2, dataset_normalize=scib_normalize, verbose=False)
    # print("2")
    # run(algorithm=alg2, dataset_normalize=scib_normalize, combiner=amwjmsi, verbose=False)
    # print("3")
    # run(algorithm=alg2, batch_normalize=scib_normalize, celltype_normalize=None, verbose=False)
    # print("4")
    # run(algorithm=alg2, batch_normalize=scib_normalize, celltype_normalize=None, combiner=amwjmsi, verbose=False)
    # print("5")
    # run(algorithm=alg2, batch_normalize=None, celltype_normalize=scib_normalize, verbose=False)
    # print("6")
    # run(algorithm=alg2, batch_normalize=None, celltype_normalize=scib_normalize, combiner=amwjmsi, verbose=False)
    # print("7")
    # run(algorithm=alg2, verbose=False)
    # print("8")
    # run(algorithm=alg2, combiner=amwjmsi, verbose=False)

    print("FINISHED!")
    pass
