import re

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd


class TestSetitemValidation:
    def _check_setitem_invalid(self, arr, invalid):
        msg = f"Invalid value '{invalid!s}' for dtype '{arr.dtype}'"
        msg = re.escape(msg)
        with pytest.raises(TypeError, match=msg):
            arr[0] = invalid

        with pytest.raises(TypeError, match=msg):
            arr[:] = invalid

        with pytest.raises(TypeError, match=msg):
            arr[[0]] = invalid

        # FIXME: don't leave commented-out
        # with pytest.raises(TypeError):
        #    arr[[0]] = [invalid]

        # with pytest.raises(TypeError):
        #    arr[[0]] = np.array([invalid], dtype=object)

        # Series non-coercion, behavior subject to change
        ser = pd.Series(arr)
        with pytest.raises(TypeError, match=msg):
            ser[0] = invalid
            # TODO: so, so many other variants of this...

    _invalid_scalars = [
        1 + 2j,
        "True",
        "1",
        "1.0",
        pd.NaT,
        np.datetime64("NaT"),
        np.timedelta64("NaT"),
    ]

    @pytest.mark.parametrize(
        "invalid", [*_invalid_scalars, 1, 1.0, np.int64(1), np.float64(1)]
    )
    def test_setitem_validation_scalar_bool(self, invalid):
        arr = pd.array([True, False, None], dtype="boolean")
        self._check_setitem_invalid(arr, invalid)

    @pytest.mark.parametrize("invalid", [*_invalid_scalars, True, 1.5, np.float64(1.5)])
    def test_setitem_validation_scalar_int(self, invalid, any_int_ea_dtype):
        arr = pd.array([1, 2, None], dtype=any_int_ea_dtype)
        self._check_setitem_invalid(arr, invalid)

    @pytest.mark.parametrize("invalid", [*_invalid_scalars, True])
    def test_setitem_validation_scalar_float(self, invalid, float_ea_dtype):
        arr = pd.array([1, 2, None], dtype=float_ea_dtype)
        self._check_setitem_invalid(arr, invalid)


@pytest.mark.parametrize(
    "dtype",
    [
        "Float64",
        pytest.param("float64[pyarrow]", marks=td.skip_if_no("pyarrow")),
    ],
)
@pytest.mark.parametrize("indexer", [1, [1], [False, True, False]])
def test_setitem_nan_in_float64_array(dtype, indexer, using_nan_is_na):
    arr = pd.array([0, pd.NA, 1], dtype=dtype)

    arr[indexer] = np.nan
    if not using_nan_is_na:
        assert np.isnan(arr[1])
    else:
        assert arr[1] is pd.NA


@pytest.mark.parametrize(
    "dtype",
    [
        "Int64",
        pytest.param("int64[pyarrow]", marks=td.skip_if_no("pyarrow")),
    ],
)
@pytest.mark.parametrize("indexer", [1, [1], [False, True, False]])
def test_setitem_nan_in_int64_array(dtype, indexer, using_nan_is_na):
    arr = pd.array([0, 1, 2], dtype=dtype)
    if not using_nan_is_na:
        err = TypeError
        msg = "Invalid value 'nan' for dtype 'Int64'"
        if dtype == "int64[pyarrow]":
            import pyarrow as pa

            err = pa.lib.ArrowInvalid
            msg = "Could not convert nan with type float"
        with pytest.raises(err, match=msg):
            arr[indexer] = np.nan
        assert arr[1] == 1
    else:
        arr[indexer] = np.nan
        assert arr[1] is pd.NA
