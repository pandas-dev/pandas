import numpy as np
import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import NumpyExtensionArray

class TestMaskedArrays:
    def test_masked_array_integer(self):
        data = np.ma.array([1, 2, 3, 4], mask=[False, True, False, True]) 
        result = pd.array(data)
        expected = pd.array([1, None, 3, None], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)
        
    def test_masked_array_float(self):
        data = np.ma.array([1.1, 2.2, 3.3, 4.4], mask=[False, True, False, True]) 
        result = pd.array(data)
        expected = pd.array([1.1, None, 3.3, None], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)    
        
    def test_masked_array_string(self):
        data = np.ma.array(["foo", "bar", "baz", "qux"], mask=[False, True, False, True]) 
        result = pd.array(data)
        expected = pd.array(["foo", None, "baz", None], dtype="string")
        tm.assert_extension_array_equal(result, expected)
        
    def test_masked_array_boolean(self):
        data = np.ma.array([True, False, True, False], mask=[False, True, False, True]) 
        result = pd.array(data)
        expected = pd.array([True, None, True, None], dtype="boolean")
        tm.assert_extension_array_equal(result, expected)
    
    def test_masked_array_no_mask(self):
        data = np.ma.array([1, 2, 3, 4], mask=[False, False, False, False]) 
        result = pd.array(data)
        # Requirement: even without masked values, return nullable EA
        expected = pd.array([1, 2, 3, 4], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)
        
    def test_masked_array_all_masked(self):
        data = np.ma.array([1, 2, 3, 4], mask=[True, True, True, True]) 
        result = pd.array(data)
        # Requirement: all-masked results in NumpyExtensionArray(object)
        expected = NumpyExtensionArray(np.array([None, None, None, None], dtype=object))
        tm.assert_extension_array_equal(result, expected)
    
    def test_masked_array_with_explicit_dtype(self):
        data = np.ma.array([1, 2, 3, 4], mask=[False, True, False, True]) 
        result = pd.array(data, dtype="Float64")
        expected = pd.array([1, None, 3, None], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)