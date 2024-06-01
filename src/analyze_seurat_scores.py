import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from library.scorer import *
from library.adata_preprocessing import *
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler



def score_celltype(anndatas, similarity_function=amwjmsi):
    cell_types = np.unique([anndata.obs['celltype'].cat.categories for anndata in anndatas]) 
    batch_types = np.unique([anndata.obs['tech'].iloc[0] for anndata in anndatas])
    num_batch_types = len(batch_types)
    df = pd.DataFrame(index=batch_types, columns=batch_types)

    for i in range(num_batch_types):
        for j in range(i, num_batch_types):
            batch_i = batch_types[i]
            batch_j = batch_types[j]
            score = similarity_function(anndatas[i], anndatas[j], cell_types)
            df.at[batch_i, batch_j] = score
            df.at[batch_j, batch_i] = score

    return df

def plot_dataframes(batches, celltypes):
    x_vals = []
    y_vals = []
    keys = []

    for i, (x, y) in enumerate(zip(batches.values.flatten(), celltypes.values.flatten())):
        if x != 0 and y != 0:
            x_vals.append(x)
            y_vals.append(y)
            column_name = batches.columns[i % len(batches.columns)][:-4]
            row_name_index = i // len(celltypes.columns)
            row_name = celltypes.index[row_name_index]
            keys.append((column_name, row_name))

    plt.figure(figsize=(8, 6))
    plt.scatter(x_vals, y_vals, alpha=0.8, edgecolors='w')
    xr, yr = np.array(x_vals).reshape(-1, 1), np.array(y_vals).reshape(-1, 1)
    reg = LinearRegression().fit(xr, yr)
    x_val = [0, np.max(x_vals)]
    y_val = [reg.intercept_[0], reg.intercept_[0] + reg.coef_[0][0]*np.max(x_vals)]
    plt.plot(x_val, y_val)
    print(reg.score(xr, yr))
    plt.xlabel('Batch Score')
    plt.ylabel('Celltype Score')
    plt.title('Batch vs Celltype (Seurat)')
    for i, key in enumerate(keys):
        plt.annotate(key, (x_vals[i], y_vals[i]))
    plt.show()


# df = pd.read_csv('data/pancreas/seurat_scores.csv')
# data = "data/preprocessed/human_pancreas_preprocessed.h5ad"
# similarity_function = mwjmsi
# adatas = filter_batches(sc.read_h5ad(data))
# scores = score_celltype(adatas, mwjmsi)
# for column in df.columns:
#     df[column].fillna(0, inplace=True)
# df = df.iloc[:, 1:]
# df_dict = df.to_dict(orient='records')
# dict_list = [{k[:-4]:v for k,v in x.items()} for x in df_dict]
# processed_dict = {}
# for i in range(len(dict_list)):
#     store_key = list(dict_list[i].keys())[i]
#     for key in dict_list[i].keys():
#         if dict_list[i][key] != 0:
#             new_key = (store_key, key)
#             processed_dict[new_key] = dict_list[i][key]

#plot_dataframes(df, scores)
# pos = nx.bipartite_layout(graph, df['celseq2.rds'])
# nx.draw(graph, pos, with_labels=True, node_color='skyblue', node_size=2000, font_size=10)
# plt.title('Complete Bipartite Graph')
# plt.show()
# adata = sc.AnnData(df)