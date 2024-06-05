# imports
    # standard libraries
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import LabelEncoder
from typing import List, Dict, Tuple
    # user libraries
from src.library.scorer import *


# functions
def plotter(
        d1: Dict[Tuple[str, str], float], 
        d2: Dict[Tuple[str, str], float], 
        save_file: str, 
        title: str="", 
        annotations: Dict[Tuple[str, str], float]=None,
        show: bool=False,
        random_colors: bool=False,
        b1_colors: bool=False) -> None:
    """
        Plots the batch vs. cell type scores
        Includes:
            Regression line
        Colors do not annotate for anything.
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

    # color and plot figure
    colors = rand_colors(x) if random_colors else batch_color(batch_name) if b1_colors else None
    plt.scatter(x, y, c=colors, s=1000)

    # calculate regression
    reg_x, reg_y, reg_score = calculate_regression(x, y)
    plt.plot(reg_x, reg_y)

    # calculate average absolute difference
    abs_diff = np.sum(np.abs(np.array(x) - np.array(y)))
    avg_abs_diff = abs_diff / np.array(x).shape[0]

    plt.title(title + f"\nR^2: {round(reg_score, 5)}\n {round(avg_abs_diff, 5)}")
    plt.xlabel("Batch")
    plt.ylabel("Cell Types (Ground Truth)")
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

def batch_color(batch_pair_name: List[Tuple[str, str]]) -> List:
    batch1 = [b1 for (b1, _) in batch_pair_name]

    # thank you to: https://stackoverflow.com/a/25730396/21168736
    color = plt.cm.rainbow(np.linspace(0, 1, len(np.unique(batch1))))
    return [color[x] for x in LabelEncoder().fit_transform(batch1)]