# Requirements and Test Oracles

## Functional Requirements
1. The system shall support reading from a JSON file
2. The system shall support reading from an XML file
3. The system shall support reading from an HTML file
4. The system shall support reading from a SQL file
5. The system shall support reading from an Excel file
6. The system shall support reading from a CSV file
7. The system shall allow the user to create the Series data object, labeled arrays of any type
8. The system shall allow the user to create the DataFrame data object, a 2D array (table with rows and columns)
9. The system shall support a human readable output of pandas objects
10. The system shall allow the user to plot a DataFrame
11. The system shall not assert the truthiness of a pandas object
...

## Non-Functional Requirements
1. The system shall function on modern python version
2. The system shall maintain compatibility with its notable dependencies (numpy, matplotlib)
3. The system shall be deterministic
...

## Test Oracles

| Requirement ID | Requirement Description | Test Oracle (Expected Behavior) |
|-----------------------|-----------------------------------|---------------------------------------------|
| FR-1                   | The system shall ………..| After adding "Buy milk"................|
| FR-2                   | The system shal…. ……..| After deleting `"Buy milk".............|
| NFR-1                | The system shall………... | When…………..within 1 second. |
| FR-4                   | ……..                                |............                                         |

| Requirement ID | Requirement Description                                                                                       | Test Oracle (Expected Behavior)                                                                                                                                                                                |
| -------------- | ------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| FR-1           | The system shall support reading from a JSON file                                                             | After providing a valid JSON file containing tabular data, pandas shall successfully create a pandas object containing the corresponding data without raising an exception.                                    |
| FR-2           | The system shall support reading from an XML file                                                             | After providing a valid XML file containing tabular data, pandas shall successfully create a pandas object containing the corresponding data without raising an exception.                                     |
| FR-3           | The system shall support reading from an HTML file                                                            | After providing a valid HTML file containing one or more tables, pandas shall successfully extract the table data into a pandas object without raising an exception.                                           |
| FR-4           | The system shall support reading from a SQL file                                                              | After providing a valid SQL data source/query, pandas shall successfully retrieve the resulting tabular data into a pandas object without raising an exception.                                                |
| FR-5           | The system shall support reading from an Excel file                                                           | After providing a valid Excel file containing tabular data, pandas shall successfully create a pandas object containing the corresponding data without raising an exception.                                   |
| FR-6           | The system shall support reading from a CSV file                                                              | After providing a valid CSV file containing supported tabular data, pandas shall successfully read the data into a pandas object without raising an exception.                                                |
| FR-7           | The system shall allow the user to create the Series data object, labeled arrays of any type                  | After providing a collection of supported values and optional labels, pandas shall create a Series whose values and labels match the provided input.                                                           |
| FR-8           | The system shall allow the user to create the DataFrame data object, a 2D array (table with rows and columns) | After providing valid two-dimensional tabular data, pandas shall create a DataFrame whose rows, columns, and values correspond to the provided input.                                                          |
| FR-9           | The system shall support a human readable output of pandas objects                                            | After requesting the string representation of a pandas Series or DataFrame, pandas shall return a human-readable representation containing the object's relevant data and structure.                           |
| FR-10          | The system shall allow the user to plot a DataFrame                                                           | After providing a DataFrame containing plottable data and requesting a plot, pandas shall produce a plot object without raising an exception.                                                                  |
| FR-11          | The system shall not assert the truthiness of a pandas object                                                 | After attempting to evaluate the truth value of a pandas Series or DataFrame, pandas shall raise the documented ambiguity-related exception rather than implicitly evaluating the object as `True` or `False`. |
| NFR-1          | The system shall function on modern Python versions                                                           | When pandas is installed on a supported modern Python version, importing pandas and executing supported basic operations shall complete without a Python-version compatibility error.                          |
| NFR-2          | The system shall maintain compatibility with its notable dependencies (NumPy, Matplotlib)                     | When supported versions of NumPy and Matplotlib are installed, pandas shall import and perform operations that depend on these packages without dependency-related errors.                                     |
| NFR-3          | The system shall be deterministic                                                                             | When the same input, environment, and operation are provided multiple times, pandas shall produce equivalent results for each execution.                                                                       |
