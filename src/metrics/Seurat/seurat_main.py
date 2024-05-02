from seurat import read_file
from seurat import adata_to_seurat
from seurat import get_scores


def main():
    adata = read_file("/Users/carsonnannini/Research/similarity-metrics/data/preprocessed/human_pancreas_preprocessed.h5ad")
    get_scores(adata)

main()