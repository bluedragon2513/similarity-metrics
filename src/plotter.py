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
    
    plt.scatter(
        x, y,
        c=[np.random.choice(plt.rcParams["axes.prop_cycle"].by_key()['color']) for _ in x]
    )

    # calculate regression
    xr, yr = np.array(x).reshape(-1, 1), np.array(y).reshape(-1, 1)
    reg = LinearRegression().fit(xr, yr)
    x_val = [0, np.max(x)]
    y_val = [reg.intercept_[0], reg.intercept_[0] + reg.coef_[0][0]*np.max(x)]
    plt.plot(x_val, y_val)

    plt.title(title + f"\n{round(reg.score(xr, yr), 5)}")
    plt.savefig(fname=save_file)
    plt.show()

    # make plot with annotations
    plt.clf()
    plt.scatter(
        x, y,
        c=[np.random.choice(plt.rcParams["axes.prop_cycle"].by_key()['color']) for _ in x]
    )
    plt.plot(x_val, y_val)

    for i, key in enumerate(keys):
        plt.annotate(key, (x[i], y[i]))

    plt.title(title + f"\n{round(reg.score(xr, yr), 5)}")
    plt.savefig(fname=save_file[:-3]+"annotations.png")
    plt.show()




if __name__ == "__main__":
    with open("data/pancreas/pancreas-batch-scores.pkl", "rb") as f:
        batch_scores = pickle.load(f)
    with open("data/pancreas/pancreas-celltype-scores.pkl", "rb") as f:
        celltype_scores = pickle.load(f)
    print(celltype_scores)
    plotter(batch_scores, celltype_scores, "data/pancreas/pancreas-mwjmsmi.png", title="Pancreas: batch vs. celltype (Scanorama)")