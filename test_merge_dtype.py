import numpy as np
import pandas as pd

def test_1():
    null_keys = pd.DataFrame(
        {
            "key": pd.Series([0, 1, pd.NA], dtype=pd.Int64Dtype()),
            "val1": pd.Series(["a", "b", "c"]),
        }
    )
    non_null_keys = pd.DataFrame(
        {
            "key": pd.Series([0, 1, 2], dtype=pd.Int64Dtype()),  # This fails!
            # "key": pd.Series([0, 1, 2], dtype=pd.UInt64Dtype()),  # This also fails!
            # "key": pd.Series([0, 1, 2], dtype=pd.Int64Dtype()),  # This works!
            "val2": pd.Series(["d", "e", "f"]),
        }
    )
    print(pd.merge(null_keys, non_null_keys, how="left", on="key").to_markdown())
    print("Test successful!")

def test_2():
    null_keys = pd.DataFrame(
        {
            "key": pd.Series([0, 1, pd.NA], dtype=pd.int64),
            "val1": pd.Series(["a", "b", "c"]),
        }
    )
    non_null_keys = pd.DataFrame(
        {
            "key": pd.Series([0, 1, 2], dtype=pd.Int64Dtype()),  # This fails!
            # "key": pd.Series([0, 1, 2], dtype=pd.UInt64Dtype()),  # This also fails!
            # "key": pd.Series([0, 1, 2], dtype=pd.Int64Dtype()),  # This works!
            "val2": pd.Series(["d", "e", "f"]),
        }
    )
    print(pd.merge(null_keys, non_null_keys, how="left", on="key").to_markdown())
    print("Test successful!")

def test_2():
    null_keys = pd.DataFrame(
        {
            "key": pd.Series([0, 1, pd.NA], dtype=pd.UInt64Dtype),
            "val1": pd.Series(["a", "b", "c"]),
        }
    )
    non_null_keys = pd.DataFrame(
        {
            "key": pd.Series([0, 1, 2], dtype=pd.Int64Dtype()),  # This fails!
            # "key": pd.Series([0, 1, 2], dtype=pd.UInt64Dtype()),  # This also fails!
            # "key": pd.Series([0, 1, 2], dtype=pd.Int64Dtype()),  # This works!
            "val2": pd.Series(["d", "e", "f"]),
        }
    )
    print(pd.merge(null_keys, non_null_keys, how="left", on="key").to_markdown())
    print("Test successful!")
