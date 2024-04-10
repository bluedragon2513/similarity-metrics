# imports
    # standard
import scanpy as sc
from concurrent.futures import ProcessPoolExecutor
    # user
# from scorer_parallel import *
from scorer import *
from adata_preprocessing import *
from plotter import plotter

# constants
data = "data/preprocessed/human_pancreas_preprocessed.h5ad"
algorithm = "scanorama"
dataset_name = "pancreas"
folder = "data/"
batch_save_file = "batch-scores"
celltype_save_file = "celltype-scores"
similarity_function = mwjmsi

# functions
def run(data=data, 
        algorithm=algorithm, 
        dataset_name=dataset_name, 
        folder=folder, 
        batch_save_file=batch_save_file, 
        celltype_save_file=celltype_save_file,
        similarity_function=similarity_function):
    adatas = filter_batches(sc.read_h5ad(data))

    scorer_batch(adatas, folder=folder, algorithm=algorithm, dataset_name=dataset_name, save_file=batch_save_file)
    scorer_celltype(adatas, folder=folder, algorithm=algorithm, dataset_name=dataset_name, save_file=celltype_save_file, similarity_function=similarity_function)
    
    with open(f"{folder}{algorithm}/{dataset_name}/{batch_save_file}.pkl", "rb") as f:
        batch_scores = pickle.load(f)
    with open(f"{folder}{algorithm}/{dataset_name}/{celltype_save_file}-{similarity_function.__name__}.pkl", "rb") as f:
        celltype_scores = pickle.load(f)

    batch_scores = {k:v for k, v in batch_scores.items() if k[0] != k[1]}
    celltype_scores = {k:v for k, v in celltype_scores.items() if k[0] != k[1]}
    plotter(batch_scores, celltype_scores, f"{folder}{algorithm}/{dataset_name}/{similarity_function.__name__}", title=f"{dataset_name.title()}: batch vs. celltype ({algorithm.title()})")

    print("done")

# main
if __name__ == "__main__":
    run()
    run(similarity_function=amwjmsi)