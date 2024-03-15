import numpy as np
# This is useful to normalize pseudotime
def minmax_normalise(arr):
        
    norm_arr = []
    arr = np.asarray(arr)
    arr_max = max(arr)
    arr_min = min(arr)
    for i in range(len(arr)):
        norm_arr.append((arr[i] - arr_min )/(arr_max  - arr_min )) 
    return norm_arr

    
def normStage(arr,start,newstart,end=999,ranges=0.4):
    arr_copy=arr.copy()
    largeLogic=(arr_copy>start)&(arr_copy<end)
    arrLarge=arr_copy[largeLogic]
    arrNew=newstart+minmax_normalise(arrLarge)*ranges
    arr_copy[largeLogic]=arrNew
    return(arr_copy)
    
# This function is similarity algrithom. we use random forest to predict similarity

def downsample_and_predict(reflatent_seurat,label ,scANVI, subsample_size):
    cell_ranks = reflatent_seurat.obs[label].cat.categories
    downsampled_cells = []

    for cell_rank in cell_ranks:
        current_cells = np.where(reflatent_seurat.obs[label] == cell_rank)[0]
        if min(len(current_cells), subsample_size) < subsample_size:
            print(f"The current ident {cell_rank} is smaller than the sample size, adjust to the ident size")

        downsampled_cells.extend(np.random.choice(current_cells, size=min(len(current_cells), subsample_size), replace=False))

    downsample_seurat = reflatent_seurat[downsampled_cells,: ]

    # Prepare data for RandomForest
    X = downsample_seurat.obsm[scANVI]
    y = LabelEncoder().fit_transform(downsample_seurat.obs[label])

    # Train a RandomForest classifier
    rf_model = RandomForestClassifier()
    rf_model.fit(X, y)

    # Get the class probabilities
    class_probs = rf_model.predict_proba(X)
    class_labels = rf_model.classes_

    # Create a DataFrame from the class probabilities
    votes = pd.DataFrame(data=class_probs, columns=class_labels)
    votes['ident'] = downsample_seurat.obs[label].values

    # Calculate class probabilities per cluster
    summary = votes.groupby('ident').mean()

    # Normalize the summary matrix
    summary = summary.divide(summary.max(axis=1), axis=0).fillna(0)
    np.fill_diagonal(summary.values, 0)
    summary.columns=summary.index.values


    return summary
