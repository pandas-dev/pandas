import ctypes

import pandas as pd


def test_func():
    libc = ctypes.PyDLL('pandas/_libs/ndframe_data.cpython-37m-darwin.so')
    df = pd.DataFrame([[1, "2", 3., 4., 5], [6, "7", 8., 9, 10]])
    obj = libc.PdOrderedArrays_New(ctypes.py_object(df), ctypes.c_int(0))
    #libc.PdOrderedArrays_Destroy(obj)
