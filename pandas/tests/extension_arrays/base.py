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

    * data
    * data_missing
    """

    @pytest.fixture
    def data(self):
        """Length-100 array for this type."""
        raise NotImplementedError

    @pytest.fixture
    def data_missing(self):
        """Length-2 array with [NA, Valid]"""
        raise NotImplementedError

    @pytest.fixture(params=['data', 'data_missing'])
    def all_data(self, request, data, data_missing):
        if request.param == 'data':
            return data
        elif request.param == 'data_missing':
            return data_missing

    @pytest.fixture
    def na_cmp(self):
        """Binary operator for comparing NA values.

        Should return a function of two arguments that returns
        True if both arguments are (scalar) NA for your type.

        By defult, uses ``operator.or``
        """
        return operator.is_

    def test_len(self, data):
        assert len(data) == 100

    def test_ndim(self, data):
        assert data.ndim == 1

    def test_can_hold_na_valid(self, data):
        assert data._can_hold_na() in {True, False}

    def test_series_constructor(self, data):
        result = pd.Series(data)
        assert result.dtype == data.dtype
        assert len(result) == len(data)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    @pytest.mark.parametrize("from_series", [True, False])
    def dataframe_constructor(self, data, from_series):
        if from_series:
            data = pd.Series(data)
        result = pd.DataFrame({"A": data})
        assert result.dtypes['A'] == data.dtype
        assert result.shape == (len(data), 1)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    def test_concat(self, data):
        result = pd.concat([
            pd.Series(data),
            pd.Series(data),
        ], ignore_index=True)
        assert len(result) == len(data) * 2

    def test_iloc(self, data):
        ser = pd.Series(data)
        result = ser.iloc[:4]
        expected = pd.Series(data[:4])
        tm.assert_series_equal(result, expected)

        result = ser.iloc[[0, 1, 2, 3]]
        tm.assert_series_equal(result, expected)

    def test_loc(self, data):
        ser = pd.Series(data)
        result = ser.loc[:3]
        expected = pd.Series(data[:4])
        tm.assert_series_equal(result, expected)

        result = ser.loc[[0, 1, 2, 3]]
        tm.assert_series_equal(result, expected)

    def test_repr(self, data):
        ser = pd.Series(data)
        assert data.dtype.name in repr(ser)

        df = pd.DataFrame({"A": data})
        repr(df)

    def test_dtype_name_in_info(self, data):
        buf = StringIO()
        pd.DataFrame({"A": data}).info(buf=buf)
        result = buf.getvalue()
        assert data.dtype.name in result

    def test_memory_usage(self, data):
        s = pd.Series(data)
        result = s.memory_usage(index=False)
        assert result == s.nbytes

    def test_is_extension_array_dtype(self, data):
        assert is_extension_array_dtype(data)
        assert is_extension_array_dtype(data.dtype)
        assert is_extension_array_dtype(pd.Series(data))
        assert isinstance(data.dtype, ExtensionDtype)

    def test_array_interface(self, data):
        result = np.array(data)
        assert result[0] == data[0]

    def test_getitem_scalar(self, data):
        result = data[0]
        assert isinstance(result, data.dtype.type)

        result = pd.Series(data)[0]
        assert isinstance(result, data.dtype.type)

    def test_getitem_scalar_na(self, data_missing, na_cmp):
        result = data_missing[0]
        assert na_cmp(result, data_missing._fill_value)

    def test_getitem_mask(self, data):
        # Empty mask, raw array
        mask = np.zeros(len(data), dtype=bool)
        result = data[mask]
        assert len(result) == 0
        assert isinstance(result, type(data))

        # Empty mask, in series
        mask = np.zeros(len(data), dtype=bool)
        result = pd.Series(data)[mask]
        assert len(result) == 0
        assert result.dtype == data.dtype

        # non-empty mask, raw array
        mask[0] = True
        result = data[mask]
        assert len(result) == 1
        assert isinstance(result, type(data))

        # non-empty mask, in series
        result = pd.Series(data)[mask]
        assert len(result) == 1
        assert result.dtype == data.dtype

    def test_getitem_slice(self, data):
        # getitem[slice] should return an array
        result = data[slice(0)]  # empty
        assert isinstance(result, type(data))

        result = data[slice(1)]  # scalar
        assert isinstance(result, type(data))

    def test_take_sequence(self, data):
        result = pd.Series(data)[[0, 1, 3]]
        assert result.iloc[0] == data[0]
        assert result.iloc[1] == data[1]
        assert result.iloc[2] == data[3]

    def test_isna(self, data_missing):
        if data_missing._can_hold_na:
            expected = np.array([True, False])
        else:
            expected = np.array([False, False])

        result = pd.isna(data_missing)
        tm.assert_numpy_array_equal(result, expected)

        result = pd.Series(data_missing).isna()
        expected = pd.Series(expected)
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("method", [
        "mean", "sum", "prod", "mad", "sem", "var", "std",
        "skew", "kurt", "median"
    ])
    def test_nuisance_dropped(self, data, method):
        data = data[:5]
        func = operator.methodcaller(method)
        df = pd.DataFrame({"A": np.arange(len(data)),
                           "B": data})
        obj = pd.DataFrame({"A": np.arange(len(data)),
                            "B": np.array(data, dtype=object)})

        assert len(func(df)) == len(func(obj))

    @pytest.mark.parametrize("method", [min, max])
    def test_reduction_orderable(self, data, method):
        data = data[:5]
        func = operator.methodcaller(method.__name__)
        df = pd.DataFrame({"A": np.arange(len(data)),
                           "B": data})
        result = func(df)
        assert len(result) == 2

        expected = method(data)
        assert result['B'] == expected

    @pytest.mark.parametrize("method", ['cummax', 'cummin'])
    @pytest.mark.xfail(reason="Assumes comparable to floating.")
    def test_cumulative_orderable(self, data, method):
        # Upcast to object
        # https://github.com/pandas-dev/pandas/issues/19296
        # assert result.dtypes['B'] == data.dtype
        data = data[:5]
        func = operator.methodcaller(method)
        df = pd.DataFrame({"A": np.arange(len(data)),
                           "B": data})
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
    def test_binops(self, data, binop):
        # Assert that binops work between DataFrames / Series with this type
        # if binops work between arrays of this type. Extra tests will be
        # needed for, e.g., Array + scalar
        data = data[:5]
        df = pd.DataFrame({
            "A": np.arange(len(data)),
            "B": data
        })

        try:
            expected = pd.DataFrame({
                "A": binop(df['A'], df['A']),
                "B": binop(df['B'].values, df['B'].values),
            })
        except Exception:
            msg = "{} not supported for {}".format(binop.__name__,
                                                   data.dtype.name)
            raise pytest.skip(msg)

        result = binop(df, df)
        tm.assert_frame_equal(result, expected)

        # homogeneous frame
        result = binop(df[['B']], df[['B']])
        tm.assert_frame_equal(result, expected[['B']])

        # series
        result = binop(df['B'], df['B'])
        tm.assert_series_equal(result, expected['B'])

    def test_as_ndarray(self, data):
        np.array(data, dtype=data.dtype.kind)

    def test_align(self, data):
        a = data[:3]
        b = data[2:5]
        r1, r2 = pd.Series(a).align(pd.Series(b, index=[1, 2, 3]))

        # TODO: assumes that the ctor can take a list of scalars of the type
        e1 = pd.Series(type(data)(list(a) + [data._fill_value]))
        e2 = pd.Series(type(data)([data._fill_value] + list(b)))
        tm.assert_series_equal(r1, e1)
        tm.assert_series_equal(r2, e2)

    @pytest.mark.xfail(reason="GH-19342")
    def test_series_given_index(self, data):
        result = pd.Series(data[:3], index=[0, 1, 2, 3, 4])
        assert result.dtype == data.dtype
        assert len(result) == 5
        assert len(result.values) == 5
        assert pd.isna(result.loc[[3, 4]]).all()

    @pytest.mark.parametrize('dropna', [True, False])
    def test_value_counts(self, all_data, dropna):
        all_data = all_data[:10]
        if dropna:
            other = np.array(all_data[~all_data.isna()])
        else:
            other = all_data

        result = pd.Series(all_data).value_counts(dropna=dropna).sort_index()
        expected = pd.Series(other).value_counts(dropna=dropna).sort_index()

        tm.assert_series_equal(result, expected)
