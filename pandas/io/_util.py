from __future__ import annotations

import pandas as pd


def _arrow_dtype_mapping(api) -> dict:
    return {
        api.int8(): pd.Int8Dtype(),
        api.int16(): pd.Int16Dtype(),
        api.int32(): pd.Int32Dtype(),
        api.int64(): pd.Int64Dtype(),
        api.uint8(): pd.UInt8Dtype(),
        api.uint16(): pd.UInt16Dtype(),
        api.uint32(): pd.UInt32Dtype(),
        api.uint64(): pd.UInt64Dtype(),
        api.bool_(): pd.BooleanDtype(),
        api.string(): pd.StringDtype(),
        api.float32(): pd.Float32Dtype(),
        api.float64(): pd.Float64Dtype(),
    }
