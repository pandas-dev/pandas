import numpy as np
import pytest

import pandas as pd
from pandas import isna


@pytest.fixture(
    params=[
        lambda s: s.tolist(),
        lambda s: s.to_dict().values(),
        lambda s: [v for _, v in s.iteritems()],
        lambda s: list(iter(s)),
    ]
)
def series_conversion(request):
    return request.param


@pytest.fixture(
    params=[
        lambda df: df.to_dict()["A"].values(),
        lambda df: [record["A"] for record in df.to_dict(orient="records")],
    ]
)
def frame_conversion(request):
    return request.param


@pytest.fixture(params=[[True, False], [False, pd.NA], [False, np.nan]])
def boolean_data(request):
    return request.param


@pytest.fixture(params=[[1.0, 2.0], [1.0, pd.NA], [1.0, np.nan]])
def float_data(request):
    return request.param


@pytest.fixture(params=[[1, 2], [1, pd.NA], [1, np.nan]])
def int_data(request):
    return request.param


class TestSeriesReturnTypesArePythonNative:
    def test_boolean(self, boolean_data, series_conversion):
        # GH 29738
        s = pd.Series(boolean_data, dtype="boolean")
        assert all(isinstance(val, bool) or isna(val) for val in series_conversion(s))

    def test_float(self, float_data, float_ea_dtype, series_conversion):
        # GH 29738
        s = pd.Series(float_data, dtype=float_ea_dtype)
        assert all(isinstance(val, float) or isna(val) for val in series_conversion(s))

    def test_int(self, int_data, any_nullable_int_dtype, series_conversion):
        # GH 29738
        s = pd.Series(int_data, dtype=any_nullable_int_dtype)
        assert all(isinstance(val, int) or isna(val) for val in series_conversion(s))


class TestFrameReturnTypesArePythonNative:
    def test_boolean(self, boolean_data, frame_conversion):
        # GH 29738
        s = pd.DataFrame({"A": boolean_data}, dtype="boolean")
        assert all(isinstance(val, bool) or isna(val) for val in frame_conversion(s))

    def test_float(self, float_data, float_ea_dtype, frame_conversion):
        # GH 29738
        s = pd.DataFrame({"A": float_data}, dtype=float_ea_dtype)
        assert all(isinstance(val, float) or isna(val) for val in frame_conversion(s))

    def test_int(self, int_data, any_nullable_int_dtype, frame_conversion):
        # GH 29738
        s = pd.DataFrame({"A": int_data}, dtype=any_nullable_int_dtype)
        assert all(isinstance(val, int) or isna(val) for val in frame_conversion(s))
