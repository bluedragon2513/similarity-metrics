# Similarity Metrics 

## Purpose:
The purpose of this project is to compare the effectiveness of methods used by algorithms such as Seurat and Scanorama when it comes to determining the integrability or similarity between two data batches. By comparing these methods, we can determine which techniques or procedures are more effective, which will ultimately help us formulate our own algorithm that can compute the similarity between batches with higher accuracy.

## Recommended Specs for Running: 
- R: version = 4.3.3.
- Python: version >= 3.11.4.

## Plans for Future
In the future, we will looking to move forward with testing each of our algorithms with other data other than the "pancreas" and "lung" datasets. 
We will also continue to investigate new methods and procedures for measuring similarity between batches. 

## Parameters for Execution: 
- algorithm: determines the algorithm being used to obtain the scores between batches. This parameter can be set to seurat, scanorama or kruskal. 
- data: A string that represents the relative or local path to the h5ad file that contains the data being evaluated. 
- batch_normalize: If this parameter is set to scib_normalize, each batch will individually be normalized after being filtered. 
- celltype_normalize: If this parameter is set to scib_normalize, each celltype batch will individually be normalized after being filtered.
- dataset_normalize: If this parameter is set to scib_normalize, the whole dataset will be normalized prior to being filtered.
- combiner: If this parameter is set to modified jaccard similarity (mwjmsi), it will use the modified jaccard similarity formula to calculate celltype scores. If it is set to adjusted modified jaccard similarity (amwjmsi), it will use the adjusted jaccard similarity formula to calculate the celltype score. 


Example: 
```
run(algorithm=scanorama, data="data/preprocessed/human_pancreas_norm_complexBatch.h5ad", combiner=amwjmsi)
```
Will run the pipeline with scanorama on without any normalization

## Preprocessing 
### Normalization:
- Dataset Normalization: Before filtering the batches, use scib_normalize() from scib to preprocess the entire dataset.
- Batch Normalization: Once the data has been filtered into batches, normalize the batches via scib_normalize(). 
- Celltype Normalization: Once the data has been filtered by celltype, normalize the celltype batches via scib_normalize().
- None Normalization: Continue to the next step with unprocessed data.
- Note: only one normalization method will be applied at a time. 

### Filtering:
- Extract individual batches from the dataset given from the h5ad file using filter_batches(). Put batches into a list.
- Extract individual celltypes from the dataset using filter_celltype(). 
        
## Algorithm Scoring:
### Scanorama:
- Use Scanorama's find_alignments() to find k mutual nearest neighbors (MNNs) for each batch in the pair.
- find_alignments():
    1) Call dimensionality_reduce() to apply Principal Component Analysis (PCA) on the batches.
    2) Use a table to count MNNs found between datasets with find_alignments_table(), fill_table() and nn_approx()
    3) Use Scanorama's formula/metric for similarity. This is defined by taking the minimum between the number of MNNs in dataset A over the number of cells in dataset A and the number of MNNs in dataset B over the number of cells in dataset B.

### Seurat:
- Save pairs of batches/celltypes as SeuratObjects in dataset_1.rds and dataset_2.rds using adata_to_seurat().
- Call run_seurat() in Seurat.R with appropriate parameters.
- If the number of cells in either dataset is less than 3, duplicate one of the cells so there are 3 cells. 
- Read objects and pass them into FindIntegrationAnchors() (found in Integration.R) with the default parameters. 
  FindIntegrationAnchors()'s parameters will be adjusted for batches/celltypes with cells less than the number of
  dimensions being reduced.
- FindIntegrationAnchors():
    1) Scale data in batches with ScaleData().
    2) Retrieve integration features with FindIntegrationFeatures().
    3) Determine pairwise combinations and proper offsets.
    4) Slim SeuratObjects and remove unnecessary content with DietSeurat().
    5) Run Canonical Correlation Analysis (CCA) as with RunCCA().
    6) Apply L2 normalization with L2Dim().
    7) Retrieve MNNs found in each batch.
- Apply Scanorama's scoring formula (described above) on the MNNs returned by FindIntegrationAnchors().

### Kruskal: 
- . . .

### gjsi:
- . . .

## Scoring Celltypes w/ Combiner
### Modified Jaccard Similarity (mwjmsi)
- When scoring celltypes, mwjmsi will calculate the minimum number of matching celltypes as well as the maximum
  number of matching celltypes across batch pairs. 
- Calculate the celltype scores using the mimimum number of matching celltypes, maximum number of matching celltypes
  and the scores generated by the algorithm. Ran in mwjmsi().

### Adjusted Modified Jaccard Similarity (amwjmsi)
- amwjmsi is slightly different to the degree that it will create a maximum weight bipartite graph of unmatched
  celltypes using the scoring algorithm.
- Following the graph, it will then calculate the minimum number of cells for each pair as well as the maximum
  number of cells and use those values (combined with modified jaccard similarity) to score the celltypes. 
  Ran in amwjmsi().

## Store in json File:
- The scores from the algorithm will be stored in a dictionary where the key represents a pair of batches/celltypes
and the value is the score.
- Four json files will be created. 
  1) batch-scores.json
  2) celltype-scores-amwjmsi.json
  3) celltype-scores-amwjmsi.json
  3) gjsi-scores.json
- The dictionary is then written into its appropriate json file. The location of the file will be dependent on the algorithm and 
normalization used. 

## Plot Data:
### Graph
- Read the dictionaries from the json files.
- Create a scatter plot where:
  1) Each point represents a pair of batches.
  2) X-axis: batch scores
  2) Y-axis: celltype scores
  3) R-Score (described below)
  4) Sum of Absolute Difference (described below) 
- Four plots should be generated
  1) mwjmsi with annotations (batch pair labels) 
  2) mwjmsi without annotations
  3) amwjmsi with annotations
  4) amwjmsi without annotations 
- The location the plots are saved in is dependent on the normalization and algorithm used.
- Written as a .png image. 
- Ran in plotter().

### R-Score:
- Represents the score of the squared residuals.
- Defines how accurate the algorithm's scores were in correlation to the ground truth. 
- Can be found under title.
- Calculated in calculate_regression().

### Sum of Absolute Difference:
- Took the sum of the absolute difference of each point divided by the total number of points. 
- Describes how integrable the dataset is and how it defines the degree to which the algorithm
  takes the celltype information into account. 
- Score can be located underneath the r-score.
             
