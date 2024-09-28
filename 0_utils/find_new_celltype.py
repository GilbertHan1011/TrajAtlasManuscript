import numpy as np
import pandas as pd

def identify_novel_cell_types(
    extended_atlas: pd.DataFrame, 
    leidenDf: pd.DataFrame, 
    uncertainty_column: str = "transf_ann_level_3_uncert", 
    percentile: int = 80
) -> dict:
    """
    Identify novel cell types based on uncertainty levels from an extended atlas.

    Parameters:
    - extended_atlas: A DataFrame containing the uncertainty annotations and core/extend classifications.
    - leidenDf: A DataFrame with cell type annotations.
    - uncertainty_column: The column name for uncertainty levels (default is "transf_ann_level_3_uncert").
    - percentile: The percentile threshold to determine novel cell types (default is 80).

    Returns:
    - novelCellType: A dictionary where keys are column names from leidenDf and values are arrays of novel cell types.
    """
    
    # Calculate the threshold for novel cell type identification
    threshold = np.percentile(
        extended_atlas.obs[uncertainty_column][extended_atlas.obs["Core_or_Extand"] == "Core"],
        percentile
    )
    
    novelCellType = {}
    
    # Iterate through each column in leidenDf
    for col in leidenDf.columns:
        # Create dummy variables for the cell types
        x = pd.get_dummies(leidenDf[col], drop_first=True)
        
        # Multiply by the uncertainty levels to assess mean uncertainty
        result = x * np.array(extended_atlas.obs[uncertainty_column]).reshape(-1, 1)
        
        # Replace zeros with NaN to avoid affecting mean calculations
        result[result == 0] = np.nan
        
        # Calculate the mean uncertainty for each cell type
        result_mean = np.nanmean(result, axis=0)
        
        # Identify novel cell types based on the mean uncertainty compared to the threshold
        novelType = np.array(result.columns[result_mean > threshold])
        
        # Store the novel cell types in the dictionary
        novelCellType[col] = novelType
    
    return novelCellType
