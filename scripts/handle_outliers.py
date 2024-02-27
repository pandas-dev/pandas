import pandas as pd
import numpy as np


"""
Detect and handle outliers in a DataFrame.

Parameters:
    data: DataFrame
    method: str, default 'z-score'. The method used for outlier detection. Options: 'z-score' or 'IQR' (Interquartile Range).
    threshold: float, default 3. The threshold for identifying outliers. Data points beyond this threshold are considered outliers.

Returns:
    DataFrame: DataFrame with outliers handled (replaced or removed).
"""
def handle_outliers(data, method='z-score', threshold=3):
    if method == 'z-score':
        z_scores = np.abs((data - data.mean()) / data.std())
        data_no_outliers = data[(z_scores < threshold).all(axis=1)]
    
    elif method == 'IQR':
        Q1 = data.quantile(0.25)
        Q3 = data.quantile(0.75)
        IQR = Q3 - Q1
        data_no_outliers = data[~((data < (Q1 - 1.5 * IQR)) | (data > (Q3 + 1.5 * IQR))).any(axis=1)]

    else:
        raise ValueError("Invalid method. Use z-score or IQR")
    
    return data_no_outliers