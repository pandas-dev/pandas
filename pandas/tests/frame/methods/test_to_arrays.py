from pandas.core.internals.construction import to_arrays
import numpy as np
from pandas.core.indexes.api import ensure_index
from numpy import array

def test_to_arrays():
    # GH 59717
    data = np.array([
        ('John', 25, 'New York', 50000),
        ('Jane', 30, 'San Francisco', 75000),
        ('Bob', 35, 'Chicago', 65000),
        ('Alice', 28, 'Los Angeles', 60000)
    ], dtype=[('name', 'U10'), ('age', 'i4'), ('city', 'U15'), ('salary', 'i4')])

    columns = ['name', 'salary', 'city']
    indexed_columns = ensure_index(columns)

    actual_arrays, actual_cols = to_arrays(data, indexed_columns)
    expected_arrays = [array(['John', 'Jane', 'Bob', 'Alice'], dtype='<U10'), 
                       array([50000, 75000, 65000, 60000], dtype=int),
                       array(['New York', 'San Francisco', 'Chicago', 'Los Angeles'], dtype='<U15')]
    
    for actual, expected in zip(actual_arrays, expected_arrays):
        assert np.array_equal(actual, expected), f"Arrays do not match:\nActual: {actual}\nExpected: {expected}"
    
    assert actual_cols.equals(indexed_columns), "Columns do not match"