import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from adata_preprocessing import filter_batches

def celltype_counts(adata):
    # plt.rc('')
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rcParams.update({'font.size' : 30})

    fig, axs = plt.subplots(3, 3, figsize=(30,30))
    axs = np.array(axs).flatten()

    adatas = filter_batches(adata)
    for ax, batch in zip(axs, adatas):
        values, counts = np.unique(batch.obs['celltype'], return_counts=True)
        ax.pie(counts, labels=values, autopct='%1.0f%%')
        ax.set_title(f"{batch.obs['tech'].iloc[0]}: {len(batch.obs['celltype'])} cells")
    fig.tight_layout()
    fig.subplots_adjust(top=9/10)
    fig.suptitle("Celltype Ratio per Batch")
    fig.savefig("data/pancreas/analysis/batch-celltype-ratio")
    fig.show()


if __name__ == "__main__":
    adata = sc.read_h5ad("data/preprocessed/human_pancreas_preprocessed.h5ad")
    celltype_counts(adata)
    # values, counts = np.unique(adata.obs['celltype'], return_counts=True)
    # print(values)
    # print(counts)
    # plt.pie(counts, labels=values, autopct='%1.0f%%')
    # plt.show()