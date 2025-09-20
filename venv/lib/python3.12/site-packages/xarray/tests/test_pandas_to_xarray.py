# This file contains code vendored from pandas
# For reference, here is a copy of the pandas copyright notice:

# BSD 3-Clause License

# Copyright (c) 2008-2011, AQR Capital Management, LLC, Lambda Foundry, Inc. and PyData Development Team
# All rights reserved.

# Copyright (c) 2011-2025, Open source contributors.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import pandas as pd
import pandas._testing as tm
import pytest
from packaging.version import Version
from pandas import (
    Categorical,
    CategoricalIndex,
    DataFrame,
    Index,
    IntervalIndex,
    MultiIndex,
    RangeIndex,
    Series,
    date_range,
    period_range,
    timedelta_range,
)

indices_dict: dict[str, Index] = {
    "object": Index([f"pandas_{i}" for i in range(10)], dtype=object),
    "string": Index([f"pandas_{i}" for i in range(10)], dtype="str"),
    "datetime": date_range("2020-01-01", periods=10),
    "datetime-tz": date_range("2020-01-01", periods=10, tz="US/Pacific"),
    "period": period_range("2020-01-01", periods=10, freq="D"),
    "timedelta": timedelta_range(start="1 day", periods=10, freq="D"),
    "range": RangeIndex(10),
    "int8": Index(np.arange(10), dtype="int8"),
    "int16": Index(np.arange(10), dtype="int16"),
    "int32": Index(np.arange(10), dtype="int32"),
    "int64": Index(np.arange(10), dtype="int64"),
    "uint8": Index(np.arange(10), dtype="uint8"),
    "uint16": Index(np.arange(10), dtype="uint16"),
    "uint32": Index(np.arange(10), dtype="uint32"),
    "uint64": Index(np.arange(10), dtype="uint64"),
    "float32": Index(np.arange(10), dtype="float32"),
    "float64": Index(np.arange(10), dtype="float64"),
    "bool-object": Index([True, False] * 5, dtype=object),
    "bool-dtype": Index([True, False] * 5, dtype=bool),
    "complex64": Index(
        np.arange(10, dtype="complex64") + 1.0j * np.arange(10, dtype="complex64")
    ),
    "complex128": Index(
        np.arange(10, dtype="complex128") + 1.0j * np.arange(10, dtype="complex128")
    ),
    "categorical": CategoricalIndex(list("abcd") * 2),
    "interval": IntervalIndex.from_breaks(np.linspace(0, 100, num=11, dtype="int")),
    "empty": Index([]),
    # "tuples": MultiIndex.from_tuples(zip(["foo", "bar", "baz"], [1, 2, 3])),
    # "mi-with-dt64tz-level": _create_mi_with_dt64tz_level(),
    # "multi": _create_multiindex(),
    "repeats": Index([0, 0, 1, 1, 2, 2]),
    "nullable_int": Index(np.arange(10), dtype="Int64"),
    "nullable_uint": Index(np.arange(10), dtype="UInt16"),
    "nullable_float": Index(np.arange(10), dtype="Float32"),
    "nullable_bool": Index(np.arange(10).astype(bool), dtype="boolean"),
    "string-python": Index(
        pd.array([f"pandas_{i}" for i in range(10)], dtype="string[python]")
    ),
}


@pytest.fixture(
    params=[
        key for key, value in indices_dict.items() if not isinstance(value, MultiIndex)
    ]
)
def index_flat(request):
    """
    index fixture, but excluding MultiIndex cases.
    """
    key = request.param
    return indices_dict[key].copy()


class TestDataFrameToXArray:
    @pytest.fixture
    def df(self):
        return DataFrame(
            {
                "a": list("abcd"),
                "b": list(range(1, 5)),
                "c": np.arange(3, 7).astype("u1"),
                "d": np.arange(4.0, 8.0, dtype="float64"),
                "e": [True, False, True, False],
                "f": Categorical(list("abcd")),
                "g": date_range("20130101", periods=4),
                "h": date_range("20130101", periods=4, tz="US/Eastern"),
            }
        )

    def test_to_xarray_index_types(self, index_flat, df):
        index = index_flat
        # MultiIndex is tested in test_to_xarray_with_multiindex
        if len(index) == 0:
            pytest.skip("Test doesn't make sense for empty index")

        from xarray import Dataset

        df.index = index[:4]
        df.index.name = "foo"
        df.columns.name = "bar"
        result = df.to_xarray()
        assert result.sizes["foo"] == 4
        assert len(result.coords) == 1
        assert len(result.data_vars) == 8
        tm.assert_almost_equal(list(result.coords.keys()), ["foo"])
        assert isinstance(result, Dataset)

        # idempotency
        # datetimes w/tz are preserved
        # column names are lost
        expected = df.copy()
        expected.columns.name = None
        tm.assert_frame_equal(result.to_dataframe(), expected)

    def test_to_xarray_empty(self, df):
        from xarray import Dataset

        df.index.name = "foo"
        result = df[0:0].to_xarray()
        assert result.sizes["foo"] == 0
        assert isinstance(result, Dataset)

    def test_to_xarray_with_multiindex(self, df):
        from xarray import Dataset

        # MultiIndex
        df.index = MultiIndex.from_product([["a"], range(4)], names=["one", "two"])
        result = df.to_xarray()
        assert result.sizes["one"] == 1
        assert result.sizes["two"] == 4
        assert len(result.coords) == 2
        assert len(result.data_vars) == 8
        tm.assert_almost_equal(list(result.coords.keys()), ["one", "two"])
        assert isinstance(result, Dataset)

        result = result.to_dataframe()
        expected = df.copy()
        expected["f"] = expected["f"].astype(
            object if Version(pd.__version__) < Version("3.0.0dev0") else str
        )
        expected.columns.name = None
        tm.assert_frame_equal(result, expected)


class TestSeriesToXArray:
    def test_to_xarray_index_types(self, index_flat):
        index = index_flat
        # MultiIndex is tested in test_to_xarray_with_multiindex

        from xarray import DataArray

        ser = Series(range(len(index)), index=index, dtype="int64")
        ser.index.name = "foo"
        result = ser.to_xarray()
        repr(result)
        assert len(result) == len(index)
        assert len(result.coords) == 1
        tm.assert_almost_equal(list(result.coords.keys()), ["foo"])
        assert isinstance(result, DataArray)

        # idempotency
        tm.assert_series_equal(result.to_series(), ser)

    def test_to_xarray_empty(self):
        from xarray import DataArray

        ser = Series([], dtype=object)
        ser.index.name = "foo"
        result = ser.to_xarray()
        assert len(result) == 0
        assert len(result.coords) == 1
        tm.assert_almost_equal(list(result.coords.keys()), ["foo"])
        assert isinstance(result, DataArray)

    def test_to_xarray_with_multiindex(self):
        from xarray import DataArray

        mi = MultiIndex.from_product([["a", "b"], range(3)], names=["one", "two"])
        ser = Series(range(6), dtype="int64", index=mi)
        result = ser.to_xarray()
        assert len(result) == 2
        tm.assert_almost_equal(list(result.coords.keys()), ["one", "two"])
        assert isinstance(result, DataArray)
        res = result.to_series()
        tm.assert_series_equal(res, ser)
