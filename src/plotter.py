# imports
    # standard libraries
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pickle
    # user libraries
from scorer import *

def plotter(d1, d2, save_file, title="", annotations=None):
    x = [d1[key] for key in d1.keys() if key in d2.keys()]
    y = [d2[key] for key in d1.keys() if key in d2.keys()]
    keys = [key for key in d1.keys() if key in d2.keys()]
    # configure figure
    plt.figure(figsize=(20,20))
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rcParams.update({'font.size' : 20})

    plt.scatter(
        x, y,
        c=[np.random.choice(plt.rcParams["axes.prop_cycle"].by_key()['color']) for _ in x],
        s=1000
    )

    # calculate regression
    xr, yr = np.array(x).reshape(-1, 1), np.array(y).reshape(-1, 1)
    reg = LinearRegression().fit(xr, yr)
    x_val = [0, np.max(x)]
    y_val = [reg.intercept_[0], reg.intercept_[0] + reg.coef_[0][0]*np.max(x)]
    plt.plot(x_val, y_val)

    plt.title(title + f"\nR^2: {round(reg.score(xr, yr), 5)}\n {round(np.sum(np.abs(np.array(x)-np.array(y))), 5)}")
    plt.xlabel("Batch")
    plt.ylabel("Celltypes (Ground Truth)")
    plt.savefig(fname=save_file+".png", format="png")
    # plt.savefig(fname=save_file+".svg", format="svg") 

    # annotate it
    if annotations==None:
        for i, key in enumerate(keys):
            plt.annotate(key, (x[i], y[i]))
    else:
        for i, key in enumerate(keys):
            plt.annotate(
                f"{' '.join([str(k) for k in key])} {round(annotations[key], 2)}", (x[i], y[i]))
    plt.savefig(fname=save_file+"-annotations.png", format="png")
    # plt.savefig(fname=save_file+"-annotations.svg", format="svg")

    plt.show()




if __name__ == "__main__":
    with open("data/scanorama/pancreas/batch-scores-adjusted.pkl", "rb") as f:
        batch_scores = pickle.load(f)
    with open("data/scanorama/pancreas/batch-scores-adjusted.pkl", "rb") as f:
        batch_scores_jaccard = pickle.load(f)
    with open("data/scanorama/pancreas/celltype-scores-adjusted.pkl", "rb") as f:
        celltype_scores = pickle.load(f)
    
    batch_scores = {k:v for k, v in batch_scores.items() if k[0] != k[1]}
    celltype_scores = {k:v for k, v in celltype_scores.items() if k[0] != k[1]}
    batch_scores_jaccard = {k:v for k, v in batch_scores_jaccard.items() if k[0] != k[1]}
    plotter(batch_scores, celltype_scores, "data/pancreas/pancreas-amwjmsmi", title="Pancreas: batch vs. celltype (Scanorama)")