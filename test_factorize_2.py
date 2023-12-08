import pandas as pd
import numpy as np

def custom_factorize(data, original_factorization=None):
    if original_factorization is None:
        # Ensuring the data is in a form accepted by pd.factorize
        if not isinstance(data, (pd.Series, pd.Index, np.ndarray)):
            data = np.asarray(data)
        return pd.factorize(data)
    
    original_uniques, original_codes = original_factorization
    unique_to_code = {unq: code for unq, code in zip(original_uniques, original_codes)}

    # Preparing an output array for new codes
    new_codes = np.empty(len(data), dtype=int)
    new_uniques = list(original_uniques)

    # Assigning new codes based on original factorization
    next_code = max(original_codes) + 1  # Starting from the next code after the max
    for i, item in enumerate(data):
        if item in unique_to_code:
            new_codes[i] = unique_to_code[item]
        else:
            unique_to_code[item] = next_code
            new_uniques.append(item)
            new_codes[i] = next_code
            next_code += 1

    return new_codes, new_uniques

# Example usage
original_data = ['a', 'b', 'c', 'c', 'd', 'e']
new_data = ['a', 'd', 'e', 'k']

# Original factorization
original_codes, original_uniques = custom_factorize(original_data)
print("Original Data:", original_data)
print("Original Codes:", original_codes)

# Apply same factorization to new data
new_codes, new_uniques = custom_factorize(new_data, original_factorization=(original_uniques, original_codes))
print("\nNew Data:", new_data)
print("New Codes:", new_codes)
