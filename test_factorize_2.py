import pandas as pd
import numpy as np

def custom_factorize(data, original_factorization=None):
    if original_factorization is None:
        return pd.factorize(data)
    
    original_uniques, original_codes = original_factorization
    unique_to_code = dict(zip(original_uniques, original_codes))

    # Convert data to a numpy array for efficient processing
    values = np.asarray(data)

    # Prepare containers for the new codes and uniques
    new_codes = np.empty(len(values), dtype=int)
    new_uniques = list(original_uniques)

    # Assign codes to new data and handle new uniques
    for i, item in enumerate(values):
        code = unique_to_code.get(item, None)
        if code is None:
            # Assign a new code and update the mapping
            code = len(new_uniques)
            unique_to_code[item] = code
            new_uniques.append(item)
        new_codes[i] = code

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