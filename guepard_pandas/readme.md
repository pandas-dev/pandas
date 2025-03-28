# Guepard-Pandas Wrapper

## Introduction
The Guepard-Pandas Wrapper is an extension of the Pandas DataFrame that integrates seamlessly with Guepard’s data versioning capabilities. This wrapper allows data engineers to use DataFrames as usual while automatically tracking versions, enabling rollback, and maintaining historical snapshots without additional effort.

## Features
- Automated version tracking for DataFrames.
- Easy rollback to previous states.
- Seamless integration with Guepard, ensuring efficient storage and retrieval.

## Example Usage
```python
import pandas as pd
from guepard_pandas.guepard_dataframe import GuepardDataFrame

# Load a DataFrame
df = GuepardDataFrame(pd.read_csv("data.csv"), dataset_id="1234")

# Modify it
df["new_col"] = df["existing_col"] * 2

# Commit the changes
df.commit("Added new column")

# List versions
print(df.list_versions())

# Rollback to an older version
df.rollback(version_id="20240326_123456")
```

## Implementation Plan
1. Prototype Development
   - Extend `pd.DataFrame` with versioning methods.
   - Implement basic version storage using Parquet or Pickle.

2. Integration with Guepard API
   - Store versions directly in Guepard’s data management system.
   - Optimize performance for large DataFrames.

3. Testing & Optimization
   - Benchmark storage and retrieval performance.
   - Validate Pandas compatibility.

## Conclusion
This wrapper offers an elegant solution to integrate version control within Pandas using Guepard, enhancing data engineering workflows while maintaining full compatibility with Pandas.

Next Steps: Review feedback and develop a proof of concept.