# imports
    # standard libraries
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import pickle
    # user libraries
from scorer import *

def plotter(d1, d2, save_file, title=""):
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

    plt.title(title + f"\n{round(reg.score(xr, yr), 5)}")
    plt.xlabel("Batch (Ground Truth)")
    plt.ylabel("Celltypes")
    plt.savefig(fname=save_file+".svg", format="svg")
    plt.savefig(fname=save_file+".png", format="png")

    # annotate it
    for i, key in enumerate(keys):
        plt.annotate(key, (x[i], y[i]))
    plt.savefig(fname=save_file+"-annotations.png", format="png")
    plt.savefig(fname=save_file+"-annotations.svg", format="svg")

    plt.show()




if __name__ == "__main__":
    with open("data/pancreas/pancreas-batch-scores-adjusted.pkl", "rb") as f:
        batch_scores = pickle.load(f)
    with open("data/pancreas/pancreas-celltype-scores-adjusted.pkl", "rb") as f:
        celltype_scores = pickle.load(f)
    # print(batch_scores)
    # print(celltype_scores)
    # print("begin")
    # for k, v in batch_scores.items():
    #     if k[0] == k[1]:
    #         print(k)
    # for k, v in celltype_scores.items():
    #     if k[0] == k[1]:
    #         print(k)
    # print("end")
    
    batch_scores = {k:v for k, v in batch_scores.items() if k[0] != k[1]}
    celltype_scores = {k:v for k, v in celltype_scores.items() if k[0] != k[1]}
    plotter(batch_scores, celltype_scores, "data/pancreas/pancreas-amwjmsmi", title="Pancreas: batch vs. celltype (Scanorama)")