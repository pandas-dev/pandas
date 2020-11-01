import pytest

from pandas.core.dtypes.common import is_bool_dtype, is_float_dtype, is_integer_dtype

import pandas as pd
from pandas import isna

from .base import BaseExtensionTests


@pytest.fixture(
    params=[
        lambda df: df.to_dict()["A"].values(),
        lambda df: [record["A"] for record in df.to_dict(orient="records")],
    ]
)
def frame_conversion(request):
    return request.param


class BaseReturnTypeTests(BaseExtensionTests):
    def get_native_dtype(self, dtype):
        if is_integer_dtype(dtype):
            return int
        elif is_float_dtype(dtype):
            return float
        elif is_bool_dtype(dtype):
            return bool
        else:
            raise ValueError("invalid dtype provided")

    @pytest.mark.parametrize(
        "func",
        [
            lambda s: s.tolist(),
            lambda s: s.to_dict().values(),
            lambda s: [v for _, v in s.iteritems()],
            lambda s: list(iter(s)),
        ],
    )
    def test_series(self, all_data, func):
        # GH 29738
        s = pd.Series(all_data)
        native_dtype = self.get_native_dtype(all_data.dtype)

        assert all(isinstance(val, native_dtype) or isna(val) for val in func(s))

    @pytest.mark.parametrize(
        "func",
        [
            lambda df: df.to_dict()["A"].values(),
            lambda df: [record["A"] for record in df.to_dict(orient="records")],
        ],
    )
    def test_frame(self, all_data, func):
        # GH 29738
        s = pd.DataFrame({"A": all_data})
        native_dtype = self.get_native_dtype(all_data.dtype)

        assert all(isinstance(val, native_dtype) or isna(val) for val in func(s))


# class TestFrameReturnTypesArePythonNative:
#     def test_boolean(self, boolean_data, frame_conversion):
#         # GH 29738
#         s = pd.DataFrame({"A": boolean_data}, dtype="boolean")
#         assert all(isinstance(val, bool) or isna(val) for val in frame_conversion(s))

#     def test_float(self, float_data, float_ea_dtype, frame_conversion):
#         # GH 29738
#         s = pd.DataFrame({"A": float_data}, dtype=float_ea_dtype)
#         assert all(isinstance(val, float) or isna(val) for val in frame_conversion(s))

#     def test_int(self, int_data, any_nullable_int_dtype, frame_conversion):
#         # GH 29738
#         s = pd.DataFrame({"A": int_data}, dtype=any_nullable_int_dtype)
#         assert all(isinstance(val, int) or isna(val) for val in frame_conversion(s))
