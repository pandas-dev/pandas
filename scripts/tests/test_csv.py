import pandas as pd
import numpy as np

# Create a DataFrame with NumPy arrays
df = pd.DataFrame({
    'id': [1, 2],
    'embedding': [np.array([0.1, 0.2, 0.3]), np.array([0.4, 0.5, 0.6])]
})

# Save to CSV
csv_file = "test_numpy_array.csv"
df.to_csv(csv_file, index=False)
print(f"Saved CSV:\n{open(csv_file).read()}")

# Read back the CSV
df_loaded = pd.read_csv(csv_file)

# Print results
print("\nLoaded DataFrame:")
print(df_loaded)

# Check if numpy arrays were converted to strings
if isinstance(df_loaded["embedding"][0], str):
    print("\nTest Passed: NumPy arrays were converted to strings in CSV")
else:
    print("\nTest Failed: NumPy arrays were not converted as expected")
