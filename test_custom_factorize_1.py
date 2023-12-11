import pandas as pd

def factorize1(data, original_factorization=None):
    if original_factorization is None:
        return pd.factorize(data)
    
    original_uniques, original_codes = original_factorization
    unique_to_code = dict(zip(original_uniques, range(len(original_uniques))))

    # Map existing data to original codes, assign new codes to new uniques
    new_codes = []
    new_uniques = list(original_uniques)
    next_code = len(original_uniques)

    for item in data:
        if item in unique_to_code:
            new_codes.append(unique_to_code[item])
        else:
            unique_to_code[item] = next_code
            new_uniques.append(item)
            new_codes.append(next_code)
            next_code += 1

    return new_codes, new_uniques

# Example usage
original_data = ['a', 'b', 'c', 'c', 'd', 'e']
new_data = ['a', 'd', 'e', 'k']

# Original factorization
# original_codes, original_uniques = factorize(original_data)
# print("Original Data:", original_data)
# print("Original Codes:", original_codes)

# # Apply same factorization to new data
# new_codes, new_uniques = custom_factorize(new_data, original_factorization=(original_uniques, original_codes))
# print("\nNew Data:", new_data)
# print("New Codes:", new_codes)


# Original factorization
original_codes, original_uniques = factorize1(original_data)
print("Original Data:", original_data)
print("Original Codes:", original_codes)

# Apply same factorization to new data
new_codes, new_uniques = factorize1(new_data)
new_codes1, new_uniques1 = factorize1(new_data, original_factorization=(original_uniques, original_codes))

print(new_codes, new_uniques)
print(new_codes1, new_uniques1)

# new codes will be 0123
# new codes1 will be 0345