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

    # ------------------------------------------------------------------------
    # Fixtures
    # ------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------
    # Interface
    # ------------------------------------------------------------------------

    def test_len(self, data):
        assert len(data) == 100

    def test_ndim(self, data):
        assert data.ndim == 1

    def test_can_hold_na_valid(self, data):
        assert data._can_hold_na in {True, False}

    def test_memory_usage(self, data):
        s = pd.Series(data)
        result = s.memory_usage(index=False)
        assert result == s.nbytes

    def test_array_interface(self, data):
        result = np.array(data)
        assert result[0] == data[0]

    def test_as_ndarray_with_dtype_kind(self, data):
        np.array(data, dtype=data.dtype.kind)

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

    # ------------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------------

    def test_series_constructor(self, data):
        result = pd.Series(data)
        assert result.dtype == data.dtype
        assert len(result) == len(data)
        assert isinstance(result._data.blocks[0], ExtensionBlock)
        assert result._data.blocks[0].values is data

    @pytest.mark.parametrize("from_series", [True, False])
    def dataframe_constructor(self, data, from_series):
        if from_series:
            data = pd.Series(data)
        result = pd.DataFrame({"A": data})
        assert result.dtypes['A'] == data.dtype
        assert result.shape == (len(data), 1)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    @pytest.mark.xfail(reason="GH-19342")
    def test_series_given_mismatched_index_raises(self, data):
        msg = 'Wrong number of items passed 3, placement implies 4'
        with tm.assert_raises_regex(ValueError, None) as m:
            pd.Series(data[:3], index=[0, 1, 2, 3, 4])

        assert m.match(msg)

    # ------------------------------------------------------------------------
    # Reshaping
    # ------------------------------------------------------------------------

    def test_concat(self, data):
        result = pd.concat([
            pd.Series(data),
            pd.Series(data),
        ], ignore_index=True)
        assert len(result) == len(data) * 2
        assert result.dtype == data.dtype
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    # ------------------------------------------------------------------------
    # Indexing - getting
    # ------------------------------------------------------------------------

    def test_iloc_series(self, data):
        ser = pd.Series(data)
        result = ser.iloc[:4]
        expected = pd.Series(data[:4])
        tm.assert_series_equal(result, expected)

        result = ser.iloc[[0, 1, 2, 3]]
        tm.assert_series_equal(result, expected)

    def test_iloc_frame(self, data):
        df = pd.DataFrame({"A": data, 'B': np.arange(len(data))})
        expected = pd.DataFrame({"A": data[:4]})

        # slice -> frame
        result = df.iloc[:4, [0]]
        tm.assert_frame_equal(result, expected)

        # sequence -> frame
        result = df.iloc[[0, 1, 2, 3], [0]]
        tm.assert_frame_equal(result, expected)

        expected = pd.Series(data[:4], name='A')

        # slice -> series
        result = df.iloc[:4, 0]
        tm.assert_series_equal(result, expected)

        # sequence -> series
        result = df.iloc[:4, 0]
        tm.assert_series_equal(result, expected)

    def test_loc_series(self, data):
        ser = pd.Series(data)
        result = ser.loc[:3]
        expected = pd.Series(data[:4])
        tm.assert_series_equal(result, expected)

        result = ser.loc[[0, 1, 2, 3]]
        tm.assert_series_equal(result, expected)

    def test_loc_frame(self, data):
        df = pd.DataFrame({"A": data, 'B': np.arange(len(data))})
        expected = pd.DataFrame({"A": data[:4]})

        # slice -> frame
        result = df.loc[:3, ['A']]
        tm.assert_frame_equal(result, expected)

        # sequence -> frame
        result = df.loc[[0, 1, 2, 3], ['A']]
        tm.assert_frame_equal(result, expected)

        expected = pd.Series(data[:4], name='A')

        # slice -> series
        result = df.loc[:3, 'A']
        tm.assert_series_equal(result, expected)

        # sequence -> series
        result = df.loc[:3, 'A']
        tm.assert_series_equal(result, expected)

    def test_is_extension_array_dtype(self, data):
        assert is_extension_array_dtype(data)
        assert is_extension_array_dtype(data.dtype)
        assert is_extension_array_dtype(pd.Series(data))
        assert isinstance(data.dtype, ExtensionDtype)

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

    # ------------------------------------------------------------------------
    # Indexing - Setting
    # ------------------------------------------------------------------------

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_setitem_scalar(self, data):
        arr = pd.Series(data)
        arr[0] = data[1]
        assert arr[0] == data[1]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_setitem_sequence(self, data):
        arr = pd.Series(data)
        original = data.copy()

        arr[[0, 1]] = [data[1], data[0]]
        assert arr[0] == original[1]
        assert arr[1] == original[0]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_setitem_sequence_broadcasts(self, data):
        arr = pd.Series(data)

        arr[[0, 1]] = data[2]
        assert arr[0] == data[2]
        assert arr[1] == data[2]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    @pytest.mark.parametrize('setter', ['loc', 'iloc'])
    def test_set_scalar(self, data, setter):
        arr = pd.Series(data)
        setter = getattr(arr, setter)
        operator.setitem(setter, 0, data[1])
        assert arr[0] == data[1]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_set_loc_scalar_mixed(self, data):
        df = pd.DataFrame({"A": np.arange(len(data)), "B": data})
        df.loc[0, 'B'] = data[1]
        assert df.loc[0, 'B'] == data[1]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_set_loc_scalar_single(self, data):
        df = pd.DataFrame({"B": data})
        df.loc[10, 'B'] = data[1]
        assert df.loc[10, 'B'] == data[1]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_set_loc_scalar_multiple_homogoneous(self, data):
        df = pd.DataFrame({"A": data, "B": data})
        df.loc[10, 'B'] = data[1]
        assert df.loc[10, 'B'] == data[1]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_set_iloc_scalar_mixed(self, data):
        df = pd.DataFrame({"A": np.arange(len(data)), "B": data})
        df.iloc[0, 1] = data[1]
        assert df.loc[0, 'B'] == data[1]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_set_iloc_scalar_single(self, data):
        df = pd.DataFrame({"B": data})
        df.iloc[10, 0] = data[1]
        assert df.loc[10, 'B'] == data[1]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_set_iloc_scalar_multiple_homogoneous(self, data):
        df = pd.DataFrame({"A": data, "B": data})
        df.iloc[10, 1] = data[1]
        assert df.loc[10, 'B'] == data[1]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_set_mask_aligned(self, data):
        ser = pd.Series(data)
        mask = np.zeros(len(data), dtype=bool)
        mask[:2] = True

        ser[mask] = data[5:7]
        assert ser[0] == data[5]
        assert ser[1] == data[6]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_set_mask_broadcast(self, data):
        ser = pd.Series(data)
        mask = np.zeros(len(data), dtype=bool)
        mask[:2] = True

        ser[mask] = data[10]
        assert ser[0] == data[10]
        assert ser[1] == data[10]

    @pytest.mark.xfail(reason="ExtensionBlock.__setitem__ not implemented.")
    def test_setitem_expand_columns(self, data):
        df = pd.DataFrame({"A": data})
        df['B'] = 1
        assert len(df.columns) == 2

    # ------------------------------------------------------------------------
    # Methods
    # ------------------------------------------------------------------------

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

    def test_dropna(self, data_missing):
        result = pd.Series(data_missing).dropna()
        expected = pd.Series(data_missing).iloc[[1]]
        tm.assert_series_equal(result, expected)

    def test_align(self, data):
        a = data[:3]
        b = data[2:5]
        r1, r2 = pd.Series(a).align(pd.Series(b, index=[1, 2, 3]))

        # TODO: assumes that the ctor can take a list of scalars of the type
        e1 = pd.Series(type(data)(list(a) + [data._fill_value]))
        e2 = pd.Series(type(data)([data._fill_value] + list(b)))
        tm.assert_series_equal(r1, e1)
        tm.assert_series_equal(r2, e2)

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

    def test_count(self, data_missing):
        df = pd.DataFrame({"A": data_missing})
        result = df.count(axis='columns')
        expected = pd.Series([0, 1])
        tm.assert_series_equal(result, expected)

    def test_dropna_series(self, data_missing):
        ser = pd.Series(data_missing)
        result = ser.dropna()
        expected = ser.iloc[[1]]
        tm.assert_series_equal(result, expected)

    def test_dropna_frame(self, data_missing):
        df = pd.DataFrame({"A": data_missing})

        # defaults
        result = df.dropna()
        expected = df.iloc[[1]]
        tm.assert_frame_equal(result, expected)

        # axis = 1
        result = df.dropna(axis='columns')
        expected = pd.DataFrame(index=[0, 1])
        tm.assert_frame_equal(result, expected)

        # multiple
        df = pd.DataFrame({"A": data_missing,
                           "B": [1, np.nan]})
        result = df.dropna()
        expected = df.iloc[:0]
        tm.assert_frame_equal(result, expected)
