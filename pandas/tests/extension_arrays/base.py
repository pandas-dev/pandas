import operator

import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas.compat import StringIO
from pandas.core.internals import ExtensionBlock
from pandas.core.dtypes.common import is_extension_array_dtype
from pandas.core.dtypes.dtypes import ExtensionDtype


class BaseDtypeTests(object):
    """Base class for ExtensionDtype classes"""

    @pytest.fixture
    def dtype(self):
        """A fixture providing the ExtensionDtype to validate."""
        raise NotImplementedError

    def test_name(self, dtype):
        assert isinstance(dtype.name, str)

    def test_kind(self, dtype):
        valid = set('biufcmMOSUV')
        if dtype.kind is not None:
            assert dtype.kind in valid

    def test_construct_from_string_own_name(self, dtype):
        result = dtype.construct_from_string(dtype.name)
        assert type(result) is type(dtype)

        # check OK as classmethod
        result = type(dtype).construct_from_string(dtype.name)
        assert type(result) is type(dtype)

    def test_is_dtype_from_name(self, dtype):
        result = type(dtype).is_dtype(dtype.name)
        assert result is True

    def test_is_dtype_from_self(self, dtype):
        result = type(dtype).is_dtype(dtype)
        assert result is True


class BaseArrayTests(object):
    """Base class for extension array classes.

    Subclasses should implement the following fixtures

    * test_data
    * test_data_missing
    """

    @pytest.fixture
    def test_data(self):
        """Length-100 array for this type."""
        raise NotImplementedError

    @pytest.fixture
    def test_data_missing(self):
        """Length-2 array with [NA, Valid]"""
        raise NotImplementedError

    @pytest.fixture
    def na_cmp(self):
        """Binary operator for comparing NA values.

        Should return a function of two arguments that returns
        True if both arguments are (scalar) NA for your type.

        By defult, uses ``operator.or``
        """
        return operator.is_

    def test_len(self, test_data):
        assert len(test_data) == 100

    def test_ndim(self, test_data):
        assert test_data.ndim == 1

    def test_can_hold_na_valid(self, test_data):
        assert test_data._can_hold_na() in {True, False}

    def test_series_constructor(self, test_data):
        result = pd.Series(test_data)
        assert result.dtype == test_data.dtype
        assert len(result) == len(test_data)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    @pytest.mark.parametrize("from_series", [True, False])
    def test_dataframe_constructor(self, test_data, from_series):
        if from_series:
            test_data = pd.Series(test_data)
        result = pd.DataFrame({"A": test_data})
        assert result.dtypes['A'] == test_data.dtype
        assert result.shape == (len(test_data), 1)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    def test_concat(self, test_data):
        result = pd.concat([
            pd.Series(test_data),
            pd.Series(test_data),
        ], ignore_index=True)
        assert len(result) == len(test_data) * 2

    def test_iloc(self, test_data):
        ser = pd.Series(test_data)
        result = ser.iloc[:4]
        expected = pd.Series(test_data[:4])
        tm.assert_series_equal(result, expected)

        result = ser.iloc[[0, 1, 2, 3]]
        tm.assert_series_equal(result, expected)

    def test_loc(self, test_data):
        ser = pd.Series(test_data)
        result = ser.loc[:3]
        expected = pd.Series(test_data[:4])
        tm.assert_series_equal(result, expected)

        result = ser.loc[[0, 1, 2, 3]]
        tm.assert_series_equal(result, expected)

    def test_repr(self, test_data):
        ser = pd.Series(test_data)
        assert test_data.dtype.name in repr(ser)

        df = pd.DataFrame({"A": test_data})
        repr(df)

    def test_dtype_name_in_info(self, test_data):
        buf = StringIO()
        pd.DataFrame({"A": test_data}).info(buf=buf)
        result = buf.getvalue()
        assert test_data.dtype.name in result

    def test_memory_usage(self, test_data):
        s = pd.Series(test_data)
        result = s.memory_usage(index=False)
        assert result == s.nbytes

    def test_is_extension_array_dtype(self, test_data):
        assert is_extension_array_dtype(test_data)
        assert is_extension_array_dtype(test_data.dtype)
        assert is_extension_array_dtype(pd.Series(test_data))
        assert isinstance(test_data.dtype, ExtensionDtype)

    def test_array_interface(self, test_data):
        result = np.array(test_data)
        assert result[0] == test_data[0]

    def test_getitem_scalar(self, test_data):
        result = test_data[0]
        assert isinstance(result, test_data.dtype.type)

        result = pd.Series(test_data)[0]
        assert isinstance(result, test_data.dtype.type)

    def test_getitem_scalar_na(self, test_data_missing, na_cmp):
        result = test_data_missing[0]
        assert na_cmp(result, test_data_missing._fill_value)

    def test_getitem_mask(self, test_data):
        # Empty mask, raw array
        mask = np.zeros(len(test_data), dtype=bool)
        result = test_data[mask]
        assert len(result) == 0
        assert isinstance(result, type(test_data))

        # Empty mask, in series
        mask = np.zeros(len(test_data), dtype=bool)
        result = pd.Series(test_data)[mask]
        assert len(result) == 0
        assert result.dtype == test_data.dtype

        # non-empty mask, raw array
        mask[0] = True
        result = test_data[mask]
        assert len(result) == 1
        assert isinstance(result, type(test_data))

        # non-empty mask, in series
        result = pd.Series(test_data)[mask]
        assert len(result) == 1
        assert result.dtype == test_data.dtype

    def test_getitem_slice(self, test_data):
        # getitem[slice] should return an array
        result = test_data[slice(0)]  # empty
        assert isinstance(result, type(test_data))

        result = test_data[slice(1)]  # scalar
        assert isinstance(result, type(test_data))

    def test_take_sequence(self, test_data):
        result = pd.Series(test_data[[0, 1, 3]])
        assert result.iloc[0] == test_data[0]
        assert result.iloc[1] == test_data[1]
        assert result.iloc[2] == test_data[3]

    def test_isna(self, test_data_missing):
        if test_data_missing._can_hold_na:
            expected = np.array([True, False])
        else:
            expected = np.array([False, False])

        result = pd.isna(test_data_missing)
        tm.assert_numpy_array_equal(result, expected)

        result = pd.Series(test_data_missing).isna()
        expected = pd.Series(expected)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("method", [
        "mean", "sum", "prod", "mad", "sem", "var", "std",
        "skew", "kurt", "median"
    ])
    def test_nuisance_dropped(self, test_data, method):
        test_data = test_data[:5]
        func = operator.methodcaller(method)
        df = pd.DataFrame({"A": np.arange(len(test_data)),
                           "B": test_data})
        assert len(func(df)) == 1

    @pytest.mark.parametrize("method", [min, max])
    def test_reduction_orderable(self, test_data, method):
        test_data = test_data[:5]
        func = operator.methodcaller(method.__name__)
        df = pd.DataFrame({"A": np.arange(len(test_data)),
                           "B": test_data})
        result = func(df)
        assert len(result) == 2

        expected = method(test_data)
        assert result['B'] == expected

    @pytest.mark.parametrize("method", ['cummax', 'cummin'])
    @pytest.mark.xfail(reason="Assumes comparable to floating.")
    def test_cumulative_orderable(self, test_data, method):
        # Upcast to object
        # https://github.com/pandas-dev/pandas/issues/19296
        # assert result.dtypes['B'] == test_data.dtype
        test_data = test_data[:5]
        func = operator.methodcaller(method)
        df = pd.DataFrame({"A": np.arange(len(test_data)),
                           "B": test_data})
        result = func(df)
        assert result.shape == df.shape

    @pytest.mark.parametrize("binop", [
        operator.add,
        operator.sub,
        operator.lt,
        operator.le,
        operator.ge,
        operator.gt,
        operator.pow,
    ], ids=lambda x: x.__name__)
    def test_binops(self, test_data, binop):
        # Assert that binops work between DataFrames / Series with this type
        # if binops work between arrays of this type. Extra tests will be
        # needed for, e.g., Array + scalar
        test_data = test_data[:5]
        df = pd.DataFrame({
            "A": np.arange(len(test_data)),
            "B": test_data
        })

        try:
            expected = pd.DataFrame({
                "A": binop(df['A'], df['A']),
                "B": binop(df['B'].values, df['B'].values),
            })
        except Exception:
            msg = "{} not supported for {}".format(binop.__name__,
                                                   test_data.dtype.name)
            raise pytest.skip(msg)

        result = binop(df, df)
        tm.assert_frame_equal(result, expected)

        # homogeneous frame
        result = binop(df[['B']], df[['B']])
        tm.assert_frame_equal(result, expected[['B']])

        # series
        result = binop(df['B'], df['B'])
        tm.assert_series_equal(result, expected['B'])

    def test_as_ndarray(self, test_data):
        np.array(test_data, dtype=test_data.dtype.kind)

    def test_align(self, test_data):
        a = test_data[:3]
        b = test_data[2:5]
        r1, r2 = pd.Series(a).align(pd.Series(b, index=[1, 2, 3]))

        # TODO: assumes that the ctor can take a list of scalars of the type
        e1 = pd.Series(type(test_data)(list(a) + [test_data._fill_value]))
        e2 = pd.Series(type(test_data)([test_data._fill_value] + list(b)))
        tm.assert_series_equal(r1, e1)
        tm.assert_series_equal(r2, e2)

    @pytest.mark.xfail(reason="GH-19342")
    def test_series_given_index(self, test_data):
        result = pd.Series(test_data[:3], index=[0, 1, 2, 3, 4])
        assert result.dtype == test_data.dtype
        assert len(result) == 5
        assert len(result.values) == 5
        assert pd.isna(result.loc[[3, 4]]).all()
