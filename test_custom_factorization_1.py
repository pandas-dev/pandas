import pandas as pd

def custom_factorize(data, original_categories=None):
    if original_categories is None:
        return pd.factorize(data)

    # Create a mapping for the original categories
    category_to_code = {category: i for i, category in enumerate(original_categories)}

    # Assign codes, adding new categories with consecutive codes
    next_code = len(original_categories)
    coded_data = []
    for item in data:
        if item in category_to_code:
            coded_data.append(category_to_code[item])
        else:
            category_to_code[item] = next_code
            coded_data.append(next_code)
            next_code += 1

    return coded_data, list(category_to_code.keys())

# Example usage
original_data = ['a', 'b', 'c', 'c', 'd', 'e']
new_data = ['a', 'd', 'e', 'k']

# Original factorization
original_codes, original_categories = custom_factorize(original_data)
print("Original Data:", original_data)
print("Original Codes:", original_codes)

# Apply same factorization to new data
new_codes, _ = custom_factorize(new_data, original_categories=original_categories)
print("\nNew Data:", new_data)
print("New Codes:", new_codes)
