import numpy as np
import pandas as pd

print("pandas version: ", pd.__version__)
print([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

array = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], order='C')

print("Numpy array is C-contiguous: ", array.flags.c_contiguous)
print(array)


dataframe = pd.DataFrame(array, index = pd.MultiIndex.from_tuples([('A', 'U'), ('A', 'V'), ('B', 'W')], names=['dim_one', 'dim_two']))
print("DataFrame is C-contiguous: ", dataframe.values.flags.c_contiguous)
print(dataframe.values)

dataframe_copy = dataframe.copy()
print("Copy of DataFrame is C-contiguous: ", dataframe_copy.values.flags.c_contiguous)
print(dataframe_copy.values)

dataframe_copy = dataframe_copy.copy()
print("Copy of Copy of DataFrame is C-contiguous: ", dataframe_copy.values.flags.c_contiguous)
print(dataframe_copy.values)

aggregated_dataframe = dataframe.groupby('dim_one').sum()
print("Aggregated copy of copy DataFrame is C-contiguous: ", aggregated_dataframe.values.flags.c_contiguous)
print(aggregated_dataframe.values)

aggregated_dataframe = dataframe.groupby('dim_one').sum()
print("Aggregated DataFrame is C-contiguous: ", aggregated_dataframe.values.flags.c_contiguous)
print(aggregated_dataframe.values)

print("===========================")


from pandas.core.layout import ArrayLayout
ArrayLayout().order = 'C'

print("pandas version: ", pd.__version__)
print([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

array = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], order='C')

print("Numpy array is C-contiguous: ", array.flags.c_contiguous)
print(array)


dataframe = pd.DataFrame(array, index = pd.MultiIndex.from_tuples([('A', 'U'), ('A', 'V'), ('B', 'W')], names=['dim_one', 'dim_two']))
print("DataFrame is C-contiguous: ", dataframe.values.flags.c_contiguous)
print(dataframe.values)

dataframe_copy = dataframe.copy()
print("Copy of DataFrame is C-contiguous: ", dataframe_copy.values.flags.c_contiguous)
print(dataframe_copy.values)

dataframe_copy = dataframe_copy.copy()
print("Copy of Copy of DataFrame is C-contiguous: ", dataframe_copy.values.flags.c_contiguous)
print(dataframe_copy.values)

aggregated_dataframe = dataframe.groupby('dim_one').sum()
print("Aggregated copy of copy DataFrame is C-contiguous: ", aggregated_dataframe.values.flags.c_contiguous)
print(aggregated_dataframe.values)

aggregated_dataframe = dataframe.groupby('dim_one').sum()
print("Aggregated DataFrame is C-contiguous: ", aggregated_dataframe.values.flags.c_contiguous)
print(aggregated_dataframe.values)

## Output in Jupyter Notebook
# pandas version:  0.23.4
# Numpy array is C-contiguous:  True
# DataFrame is C-contiguous:  True
# Copy of DataFrame is C-contiguous:  False
# Aggregated DataFrame is C-contiguous:  False