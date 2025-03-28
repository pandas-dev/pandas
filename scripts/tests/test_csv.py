import pandas as pd

print(pd.__file__)
print(pd.__version__)

import numpy as np

# # Create a DataFrame with NumPy arrays
# df = pd.DataFrame({
#     'id': [1, 2],
#     'embedding': [np.array([0.1, 0.2, 0.3]), np.array([0.4, 0.5, 0.6])]
# })

# # Save to CSV (where your custom preserve_complex logic resides)
# csv_file = "test_numpy_array.csv"
# df.to_csv(csv_file, index=False, preserve_complex=True)

# # Read back the raw CSV content (as text only)
# with open(csv_file, "r") as f:
#     csv_content = f.read()

# print(f"Saved CSV:\n{csv_content}")

# # Simple test: check that our JSON-ified arrays are present in the CSV text
# try:
#     assert "[0.1, 0.2, 0.3]" in csv_content
#     assert "[0.4, 0.5, 0.6]" in csv_content
#     print("\nTest Passed: The CSV output includes JSON-serialized arrays for 'embedding'.")
# except AssertionError:
#     print("\nTest Failed: The CSV does not appear to have JSON-serialized arrays as expected!")
#     raise



# TEST2
# Create a DataFrame with NumPy arrays
df = pd.DataFrame({
    "id": [1, 2],
    "embedding": [np.array([0.1, 0.2, 0.3]), np.array([0.4, 0.5, 0.6])]
})

# Save to CSV
csv_file = "test_numpy_array.csv"
df.to_csv(csv_file, index=False, preserve_complex=True)
print(f"Saved CSV:\n{open(csv_file).read()}")

# Read back the CSV
df_loaded = pd.read_csv(csv_file, preserve_complex=True)

# Print results
print("\nLoaded DataFrame:")
print(df_loaded)

# ✅ **Make the test fail by checking if we correctly load NumPy arrays**
try:
    assert isinstance(df_loaded["embedding"][0], np.ndarray), "Test Failed: Embeddings were not preserved as NumPy arrays!"
    print("\nTest Passed: Embeddings were correctly preserved as NumPy arrays")
except AssertionError as e:
    print("\nTest Failed: Pandas does not preserve NumPy arrays in CSV, needs improvement!")
    raise e
