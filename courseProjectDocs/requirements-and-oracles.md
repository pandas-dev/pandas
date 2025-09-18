# Requirements and Test Oracles

## Functional Requirements

- **FR-1**: The system shall handle missing data by representing missing values as NaN, NA or NaT in both floating-point and non-floating-point data.  
- **FR-2**: The system shall support size mutability of tabular structures, allowing columns to be inserted or deleted from a DataFrame or higher-dimensional object.  
- **FR-3**: The system shall automatically and explicitly align data when performing operations on objects, ensuring labels are aligned or allowing the user to ignore labels for automatic alignment.  
- **FR-4**: The system shall provide flexible group-by functionality to perform split-apply-combine operations for aggregating or transforming data.  
- **FR-5**: The system shall provide robust I/O tools for loading data from flat files (CSV and delimited), Excel files and databases and for saving/loading data using the ultrafast HDF5 format.  
- **FR-6**: The system shall provide time-series-specific functionality such as date-range generation, frequency conversion, moving-window statistics, and date shifting/lagging.  

## Non-Functional Requirements

- **NFR-1**: The system shall provide fast, flexible and expressive data structures designed to make working with relational or labeled data easy and intuitive.  
- **NFR-2**: The system shall be powerful and flexible, aiming to be the most powerful open-source data analysis/manipulation tool available.  
- **NFR-3**: The system shall provide robust I/O capabilities that load and save data efficiently, including the ultrafast HDF5 format.  

---

## Test Oracles

| Requirement ID | Requirement Description | Test Oracle (Expected Behavior) |
|----------------|--------------------------|----------------------------------|
| **FR-1** | Handle missing data with NaN/NA/NaT representations | When a DataFrame column contains a missing value, the system should represent it as NaN (or NA/NaT for date types) and subsequent computations should treat the value as missing. |
| **FR-2** | Support size mutability – columns can be inserted/deleted | After inserting a new column into a DataFrame, the number of columns increases and the new column is accessible by label; after deleting it, the column should no longer exist and the shape of the DataFrame reflects the removal. |
| **FR-3** | Automatic and explicit data alignment across objects | When adding two Series objects with misaligned indexes, the system should align on index labels and introduce missing values where labels do not match. |
| **FR-4** | Provide flexible group-by functionality | When grouping a DataFrame by a categorical column and applying a sum aggregation, the resulting object should contain aggregated sums for each group that equal the sum of values in the original DataFrame for that group. |
| **FR-5** | Robust I/O tools for loading and saving data | Reading a CSV file containing 100 rows and 5 columns should create a DataFrame with 100 rows and 5 columns and values that match the file; saving to HDF5 and then reloading should yield an identical DataFrame. |
| **FR-6** | Time-series-specific functionality | Generating a date range between “2023-01-01” and “2023-01-10” with a daily frequency should produce a sequence of 10 dates; shifting the resulting series by one period should move each date forward by one day. |
| **NFR-1** | Provide fast, flexible and expressive data structures | Creating and slicing a DataFrame with 10,000 rows should complete within an acceptable threshold (e.g., under 50 ms) in standard hardware, reflecting expected performance. |
| **NFR-2** | Be a powerful and flexible open-source data analysis tool | The API should allow users to chain multiple operations (e.g., filtering, grouping and aggregation) in a single fluent expression; the resulting code should remain readable and the operations should execute correctly. |
| **NFR-3** | Provide robust I/O capabilities | Loading a large CSV file (e.g., 1 GB) and saving it to HDF5 should not crash and should complete without data corruption; memory usage should remain within reasonable bounds relative to the file size. |
