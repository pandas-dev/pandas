import pandas as pd
import numpy as np

# Create a DataFrame with NumPy arrays
df = pd.DataFrame({
    'id': [1, 2],
    'embedding': [np.array([0.1, 0.2, 0.3]), np.array([0.4, 0.5, 0.6])]
})

# Save to CSV
csv_file = "test_numpy_array.csv"
df.to_csv(csv_file, index=False, preserve_complex=True)
print(f"Saved CSV:\n{open(csv_file).read()}")

# Read back the CSV
df_loaded = pd.read_csv(csv_file)

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
