import ctypes

import pandas as pd


def test_func():
    # TODO: this is hard coded to only work on Mac
    libc = ctypes.PyDLL('pandas/_libs/ndframe_data.cpython-37m-darwin.so')
    df = pd.DataFrame([[1, "2", 3., 4., 5], [6, "7", 8., 9, 10]])
    obj = libc.PdOrderedArrays_New(ctypes.py_object(df), ctypes.c_int(0))
    # TODO: Improve test by accessing numpy array elements sequentially
    # and testing for appropriate values
    #libc.PdOrderedArrays_Destroy(obj)
