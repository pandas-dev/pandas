import numpy as np

def percentile_scaling(data):
    data = np.array(data)
    min_val = np.min(data)
    max_val = np.max(data)
    if max_val == min_val:
        raise ValueError("Cannot scale data with identical values.")

    scaled = 100 * (data - min_val) / (max_val - min_val)
    return scaled.tolist()
