# imports
    # standard libraries
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from typing import List, Dict, Tuple
    # user libraries
from src.library.scorer import *

def plotter(
        d1: Dict[Tuple[str, str], float], 
        d2: Dict[Tuple[str, str], float], 
        save_file: str, 
        title: str="", 
        annotations: Dict[Tuple[str, str], float]=None,
        show: bool=False) -> None:
    """
        
    """
    # Get the x, y values + batch names
    x = [d1[key] for key in d1.keys() if key in d2.keys()]
    y = [d2[key] for key in d1.keys() if key in d2.keys()]
    batch_name = [key for key in d1.keys() if key in d2.keys()]

    # configure figure
    plt.figure(figsize=(20,20))
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rcParams.update({'font.size' : 20})

    plt.scatter(x, y, c=rand_colors(x), s=1000)

    # calculate regression
    reg_x, reg_y, reg_score = calculate_regression(x, y)
    plt.plot(reg_x, reg_y)

    plt.title(title + f"\nR^2: {round(reg_score, 5)}\n {round(np.sum(np.abs(np.array(x)-np.array(y))), 5)}")
    plt.xlabel("Batch")
    plt.ylabel("Celltypes (Ground Truth)")
    plt.savefig(fname=save_file+".png", format="png")
    # plt.savefig(fname=save_file+".svg", format="svg") 

    # annotate it
    if annotations==None:
        for i, key in enumerate(batch_name):
            plt.annotate(key, (x[i], y[i]))
    else:
        for i, key in enumerate(batch_name):
            plt.annotate(
                f"{' '.join([str(k) for k in key])} {round(annotations[key], 2)}", (x[i], y[i]))
    plt.savefig(fname=save_file+"-annotations.png", format="png")
    # plt.savefig(fname=save_file+"-annotations.svg", format="svg")

    if show:
        plt.show()


def calculate_regression(
        x: List[float], 
        y: List[float]) -> Tuple[List[float], List[float], float]:
    """
        Calculates the regression statistics
        Returns the regression line defined by x_val, y_val
        Returns the R^2 (squared residuals)

        Parameter Examples:
            x = [1.0, 2.0, 3.0, ..., x_20]
            y = [1.0, 2.0, 3.0, ..., y_20]

        Return Examples:
            x_val = [0, x_20]
            y_val = [0, y_20]
            reg.score(x,y) = 1
    """
    x, y = np.array(x).reshape(-1, 1), np.array(y).reshape(-1, 1)

    reg = LinearRegression().fit(x, y)
    x_val = [0, np.max(x)]
    y_val = [reg.intercept_[0], reg.intercept_[0] + reg.coef_[0][0]*np.max(x)]

    return x_val, y_val, reg.score(x, y)

def rand_colors(x: List) -> List:
    return [np.random.choice(plt.rcParams["axes.prop_cycle"].by_key()['color']) for _ in x]

# if __name__ == "__main__":
#     with open("data/scanorama/pancreas/batch-scores-adjusted.pkl", "rb") as f:
#         batch_scores = pickle.load(f)
#     with open("data/scanorama/pancreas/batch-scores-adjusted.pkl", "rb") as f:
#         batch_scores_jaccard = pickle.load(f)
#     with open("data/scanorama/pancreas/celltype-scores-adjusted.pkl", "rb") as f:
#         celltype_scores = pickle.load(f)
    
#     batch_scores = {k:v for k, v in batch_scores.items() if k[0] != k[1]}
#     celltype_scores = {k:v for k, v in celltype_scores.items() if k[0] != k[1]}
#     batch_scores_jaccard = {k:v for k, v in batch_scores_jaccard.items() if k[0] != k[1]}
#     plotter(batch_scores, celltype_scores, "data/pancreas/pancreas-amwjmsmi", title="Pancreas: batch vs. celltype (Scanorama)")