import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

def test_get_indexer_pd_na_object_index():
    # GH#65419
    idx = pd.Index([np.nan, "b"])
    
    # Scalar lookup
    assert idx.get_loc(pd.NA) == 0
    
    # List lookup
    res_list = idx.get_indexer([pd.NA])
    tm.assert_numpy_array_equal(res_list, np.array([0], dtype=np.intp))
    
    # ndarray lookup
    res_ndarray = idx.get_indexer(np.array([pd.NA], dtype=object))
    tm.assert_numpy_array_equal(res_ndarray, np.array([0], dtype=np.intp))
    
    # Index lookup
    res_index = idx.get_indexer(pd.Index([pd.NA]))
    tm.assert_numpy_array_equal(res_index, np.array([0], dtype=np.intp))

def test_drop_pd_na_object_index():
    # GH#65419
    idx = pd.Index([np.nan, "b"])
    
    # Drop with pd.NA should work and drop the NaN-labeled entry
    result = idx.drop([pd.NA])
    expected = pd.Index(["b"])
    tm.assert_index_equal(result, expected)
    
    # Drop with np.nan should also work
    result = idx.drop([np.nan])
    tm.assert_index_equal(result, expected)
