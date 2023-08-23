# Needed for new arrow string dtype

import pandas as pd

object_pyarrow_numpy = ("object", "string[pyarrow_numpy]")


def _convert_na_value(ser, expected):
    if ser.dtype != object:
        # GH#18463
        expected = expected.fillna(pd.NA)
    return expected
