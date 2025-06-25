import pandas as pd
import numpy as np

# Sample DataFrames
df1 = pd.DataFrame({
    "id": [1, 2, 3],
    "value1": ["A", "B", "C"]
})

df2 = pd.DataFrame({
    "id": [2, 3, 4],
    "value2": ["X", "Y", "Z"]
})

print("DataFrame 1:")
print(df1)

print("\nDataFrame 2:")
print(df2)

# Merge with coalesce_keys=True (default behavior: merge keys into one column)
merged_coalesced = pd.merge(
    df1, df2,
    on="id",
    how="outer",
    coalesce_keys=True,
    suffixes=("", "_right"),
)

# Merge with coalesce_keys=False (new behavior: keep left and right keys separately)
merged_separated = pd.merge(
    df1, df2,
    on="id",
    how="outer",
    coalesce_keys=False,
    suffixes=("", "_right"),
)

print("\nMerge result with coalesce_keys=True (merged 'id' column):")
print(merged_coalesced)

print("\nMerge result with coalesce_keys=False (separate 'id' and 'id_right' columns):")
print(merged_separated)

