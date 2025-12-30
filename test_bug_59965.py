import pandas as pd
import numpy as np

s1 = pd.Series({"a": 0.0, "b": 1, "c": 1, "d": 0})
s2 = pd.Series({"a": 0.0, "b": 2, "c": 2, "d": 2})
s3 = s1/s2
print("s3 (regular float64):", s3)
print("s3.mean(skipna=True):", s3.mean(skipna=True))
print()

s4 = s1.convert_dtypes()/s2.convert_dtypes()
print("s4 (FloatingArray):", s4)
print("s4 dtype:", s4.dtype) 
print("s4 values:", s4.values)
print("s4._values:", s4._values)
print("s4._values._data:", s4._values._data)
print("s4._values._mask:", s4._values._mask)
print("s4.mean(skipna=True):", s4.mean(skipna=True))
print("Expected: 0.333...")
print()

s5 = pd.Series([None,0.5,0.5,0]).convert_dtypes()
print("s5 (FloatingArray with None):", s5)
print("s5.mean(skipna=True):", s5.mean(skipna=True))
