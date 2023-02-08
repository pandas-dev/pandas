import pandas._testing as tm
import numpy as np
from pandas.core.construction import array as array_test

def test_object_to_bool_array():
    array1 = array_test(np.array([True, False]), dtype="Int64")  
    array2 = array_test(np.array([True, False], dtype=object), dtype="Int64")
    tm.assert_numpy_array_equal(array1, array2)