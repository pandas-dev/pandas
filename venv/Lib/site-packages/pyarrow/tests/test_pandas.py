# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

import gc
import decimal
import json
import multiprocessing as mp
import sys
import warnings

from collections import OrderedDict
from datetime import date, datetime, time, timedelta, timezone

import hypothesis as h
import hypothesis.strategies as st
import pytest
try:
    import numpy as np
    import numpy.testing as npt
    try:
        _np_VisibleDeprecationWarning = np.VisibleDeprecationWarning
    except AttributeError:
        from numpy.exceptions import (
            VisibleDeprecationWarning as _np_VisibleDeprecationWarning
        )
except ImportError:
    np = None

from pyarrow.pandas_compat import get_logical_type, _pandas_api
from pyarrow.tests.util import invoke_script, random_ascii, rands
import pyarrow.tests.strategies as past
import pyarrow.tests.util as test_util
from pyarrow.vendored.version import Version

import pyarrow as pa
try:
    from pyarrow import parquet as pq
except ImportError:
    pass

try:
    import pandas as pd
    import pandas.testing as tm
    from .pandas_examples import dataframe_with_arrays, dataframe_with_lists
except ImportError:
    pass


# Marks all of the tests in this module
pytestmark = pytest.mark.pandas


def _alltypes_example(size=100):
    return pd.DataFrame({
        'uint8': np.arange(size, dtype=np.uint8),
        'uint16': np.arange(size, dtype=np.uint16),
        'uint32': np.arange(size, dtype=np.uint32),
        'uint64': np.arange(size, dtype=np.uint64),
        'int8': np.arange(size, dtype=np.int16),
        'int16': np.arange(size, dtype=np.int16),
        'int32': np.arange(size, dtype=np.int32),
        'int64': np.arange(size, dtype=np.int64),
        'float32': np.arange(size, dtype=np.float32),
        'float64': np.arange(size, dtype=np.float64),
        'bool': np.random.randn(size) > 0,
        'datetime[s]': np.arange("2016-01-01T00:00:00.001", size,
                                 dtype='datetime64[s]'),
        'datetime[ms]': np.arange("2016-01-01T00:00:00.001", size,
                                  dtype='datetime64[ms]'),
        'datetime[us]': np.arange("2016-01-01T00:00:00.001", size,
                                  dtype='datetime64[us]'),
        'datetime[ns]': np.arange("2016-01-01T00:00:00.001", size,
                                  dtype='datetime64[ns]'),
        'timedelta64[s]': np.arange(0, size, dtype='timedelta64[s]'),
        'timedelta64[ms]': np.arange(0, size, dtype='timedelta64[ms]'),
        'timedelta64[us]': np.arange(0, size, dtype='timedelta64[us]'),
        'timedelta64[ns]': np.arange(0, size, dtype='timedelta64[ns]'),
        'str': [str(x) for x in range(size)],
        'str_with_nulls': [None] + [str(x) for x in range(size - 2)] + [None],
        'empty_str': [''] * size
    })


def _check_pandas_roundtrip(df, expected=None, use_threads=False,
                            expected_schema=None,
                            check_dtype=True, schema=None,
                            preserve_index=False,
                            as_batch=False):
    klass = pa.RecordBatch if as_batch else pa.Table
    table = klass.from_pandas(df, schema=schema,
                              preserve_index=preserve_index,
                              nthreads=2 if use_threads else 1)
    result = table.to_pandas(use_threads=use_threads)

    if expected_schema:
        # all occurrences of _check_pandas_roundtrip passes expected_schema
        # without the pandas generated key-value metadata
        assert table.schema.equals(expected_schema)

    if expected is None:
        expected = df

        for col in expected.columns:
            if expected[col].dtype == 'object':
                expected[col] = expected[col].replace({np.nan: None})

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", "elementwise comparison failed", DeprecationWarning)
        tm.assert_frame_equal(result, expected, check_dtype=check_dtype,
                              check_index_type=('equiv' if preserve_index
                                                else False))


def _check_series_roundtrip(s, type_=None, expected_pa_type=None):
    arr = pa.array(s, from_pandas=True, type=type_)

    if type_ is not None and expected_pa_type is None:
        expected_pa_type = type_

    if expected_pa_type is not None:
        assert arr.type == expected_pa_type

    result = pd.Series(arr.to_pandas(), name=s.name)
    tm.assert_series_equal(s, result)


def _check_array_roundtrip(values, expected=None, mask=None,
                           type=None):
    arr = pa.array(values, from_pandas=True, mask=mask, type=type)
    result = arr.to_pandas()

    values_nulls = pd.isnull(values)
    if mask is None:
        assert arr.null_count == values_nulls.sum()
    else:
        assert arr.null_count == (mask | values_nulls).sum()

    if expected is None:
        if mask is None:
            expected = pd.Series(values)
        else:
            expected = pd.Series(values).copy()
            expected[mask.copy()] = None

        if expected.dtype == 'object':
            expected = expected.replace({np.nan: None})

    tm.assert_series_equal(pd.Series(result), expected, check_names=False)


def _check_array_from_pandas_roundtrip(np_array, type=None):
    arr = pa.array(np_array, from_pandas=True, type=type)
    result = arr.to_pandas()
    npt.assert_array_equal(result, np_array)


class TestConvertMetadata:
    """
    Conversion tests for Pandas metadata & indices.
    """

    def test_non_string_columns(self):
        df = pd.DataFrame({0: [1, 2, 3]})
        table = pa.Table.from_pandas(df)
        assert table.field(0).name == '0'

    def test_non_string_columns_with_index(self):
        df = pd.DataFrame({0: [1.0, 2.0, 3.0], 1: [4.0, 5.0, 6.0]})
        df = df.set_index(0)

        # assert that the from_pandas raises the warning
        with pytest.warns(UserWarning):
            table = pa.Table.from_pandas(df)
            assert table.field(0).name == '1'

        expected = df.copy()
        # non-str index name will be converted to str
        expected.index.name = str(expected.index.name)
        with pytest.warns(UserWarning):
            _check_pandas_roundtrip(df, expected=expected,
                                    preserve_index=True)

    def test_from_pandas_with_columns(self):
        df = pd.DataFrame({0: [1, 2, 3], 1: [1, 3, 3], 2: [2, 4, 5]},
                          columns=[1, 0])

        table = pa.Table.from_pandas(df, columns=[0, 1])
        expected = pa.Table.from_pandas(df[[0, 1]])
        assert expected.equals(table)

        record_batch_table = pa.RecordBatch.from_pandas(df, columns=[0, 1])
        record_batch_expected = pa.RecordBatch.from_pandas(df[[0, 1]])
        assert record_batch_expected.equals(record_batch_table)

    def test_column_index_names_are_preserved(self):
        df = pd.DataFrame({'data': [1, 2, 3]})
        df.columns.names = ['a']
        _check_pandas_roundtrip(df, preserve_index=True)

    @pytest.mark.parametrize("tz", [None, "Europe/Brussels"])
    def test_column_index_names_datetime(self, tz):
        # ARROW-13756
        # Bug if index is timezone aware DataTimeIndex

        df = pd.DataFrame(
            np.random.randn(5, 3),
            columns=pd.date_range("2021-01-01", periods=3, freq="50D", tz=tz)
        )
        _check_pandas_roundtrip(df, preserve_index=True)

    def test_column_index_names_with_decimal(self):
        # GH-41503: Test valid roundtrip with decimal value in column index
        df = pd.DataFrame(
            [[decimal.Decimal(5), decimal.Decimal(6)]],
            columns=pd.MultiIndex.from_product(
                [[decimal.Decimal(1)], [decimal.Decimal(2), decimal.Decimal(3)]]
            ),
            index=[decimal.Decimal(4)],
        )
        _check_pandas_roundtrip(df, preserve_index=True)

    def test_range_index_shortcut(self):
        # ARROW-1639
        index_name = 'foo'
        df = pd.DataFrame({'a': [1, 2, 3, 4]},
                          index=pd.RangeIndex(0, 8, step=2, name=index_name))

        df2 = pd.DataFrame({'a': [4, 5, 6, 7]},
                           index=pd.RangeIndex(0, 4))

        table = pa.Table.from_pandas(df)
        table_no_index_name = pa.Table.from_pandas(df2)

        # The RangeIndex is tracked in the metadata only
        assert len(table.schema) == 1

        result = table.to_pandas()
        tm.assert_frame_equal(result, df)
        assert isinstance(result.index, pd.RangeIndex)
        assert _pandas_api.get_rangeindex_attribute(result.index, 'step') == 2
        assert result.index.name == index_name

        result2 = table_no_index_name.to_pandas()
        tm.assert_frame_equal(result2, df2)
        assert isinstance(result2.index, pd.RangeIndex)
        assert _pandas_api.get_rangeindex_attribute(result2.index, 'step') == 1
        assert result2.index.name is None

    def test_range_index_force_serialization(self):
        # ARROW-5427: preserve_index=True will force the RangeIndex to
        # be serialized as a column rather than tracked more
        # efficiently as metadata
        df = pd.DataFrame({'a': [1, 2, 3, 4]},
                          index=pd.RangeIndex(0, 8, step=2, name='foo'))

        table = pa.Table.from_pandas(df, preserve_index=True)
        assert table.num_columns == 2
        assert 'foo' in table.column_names

        restored = table.to_pandas()
        tm.assert_frame_equal(restored, df)

    def test_rangeindex_doesnt_warn(self):
        # ARROW-5606: pandas 0.25 deprecated private _start/stop/step
        # attributes -> can be removed if support < pd 0.25 is dropped
        df = pd.DataFrame(np.random.randn(4, 2), columns=['a', 'b'])

        with warnings.catch_warnings():
            warnings.simplefilter(action="error")
            # make_block deprecation in pandas, still under discussion
            # https://github.com/pandas-dev/pandas/pull/56422
            # https://github.com/pandas-dev/pandas/issues/40226
            warnings.filterwarnings(
                "ignore", "make_block is deprecated", DeprecationWarning
            )
            _check_pandas_roundtrip(df, preserve_index=True)

    def test_multiindex_columns(self):
        columns = pd.MultiIndex.from_arrays([
            ['one', 'two'], ['X', 'Y']
        ])
        df = pd.DataFrame([(1, 'a'), (2, 'b'), (3, 'c')], columns=columns)
        _check_pandas_roundtrip(df, preserve_index=True)

    def test_multiindex_columns_with_dtypes(self):
        columns = pd.MultiIndex.from_arrays(
            [
                ['one', 'two'],
                pd.DatetimeIndex(['2017-08-01', '2017-08-02']),
            ],
            names=['level_1', 'level_2'],
        )
        df = pd.DataFrame([(1, 'a'), (2, 'b'), (3, 'c')], columns=columns)
        _check_pandas_roundtrip(df, preserve_index=True)

    def test_multiindex_with_column_dtype_object(self):
        # ARROW-3651 & ARROW-9096
        # Bug when dtype of the columns is object.

        # uinderlying dtype: integer
        df = pd.DataFrame([1], columns=pd.Index([1], dtype=object))
        _check_pandas_roundtrip(df, preserve_index=True)

        # underlying dtype: floating
        df = pd.DataFrame([1], columns=pd.Index([1.1], dtype=object))
        _check_pandas_roundtrip(df, preserve_index=True)

        # underlying dtype: datetime
        # ARROW-9096: a simple roundtrip now works
        df = pd.DataFrame([1], columns=pd.Index(
            [datetime(2018, 1, 1)], dtype="object"))
        _check_pandas_roundtrip(df, preserve_index=True)

    def test_multiindex_columns_unicode(self):
        columns = pd.MultiIndex.from_arrays([['あ', 'い'], ['X', 'Y']])
        df = pd.DataFrame([(1, 'a'), (2, 'b'), (3, 'c')], columns=columns)
        _check_pandas_roundtrip(df, preserve_index=True)

    def test_multiindex_doesnt_warn(self):
        # ARROW-3953: pandas 0.24 rename of MultiIndex labels to codes
        columns = pd.MultiIndex.from_arrays([['one', 'two'], ['X', 'Y']])
        df = pd.DataFrame([(1, 'a'), (2, 'b'), (3, 'c')], columns=columns)

        with warnings.catch_warnings():
            warnings.simplefilter(action="error")
            # make_block deprecation in pandas, still under discussion
            # https://github.com/pandas-dev/pandas/pull/56422
            # https://github.com/pandas-dev/pandas/issues/40226
            warnings.filterwarnings(
                "ignore", "make_block is deprecated", DeprecationWarning
            )
            _check_pandas_roundtrip(df, preserve_index=True)

    def test_multiindex_rangeindex(self):
        # https://github.com/apache/arrow/issues/33473
        multiindex = pd.MultiIndex.from_arrays(
            [pd.RangeIndex(0, 2), pd.Index([1, 2])]
        )
        df = pd.DataFrame(pd.Series([1, 2], name="a"), index=multiindex)
        _check_pandas_roundtrip(df, preserve_index=None)

    def test_integer_index_column(self):
        df = pd.DataFrame([(1, 'a'), (2, 'b'), (3, 'c')])
        _check_pandas_roundtrip(df, preserve_index=True)

    def test_float_column_index_with_missing(self):
        df = pd.DataFrame([(1, 'a'), (2, 'b'), (3, 'c')], columns=[1.5, np.nan])
        _check_pandas_roundtrip(df, preserve_index=True)

    @pytest.mark.filterwarnings(
        "ignore:The DataFrame has column names of mixed type:UserWarning"
    )
    def test_string_column_index_with_missing(self):
        df = pd.DataFrame([(1, 'a'), (2, 'b'), (3, 'c')], columns=["A", None])
        _check_pandas_roundtrip(df, preserve_index=True)

    def test_index_metadata_field_name(self):
        # test None case, and strangely named non-index columns
        df = pd.DataFrame(
            [(1, 'a', 3.1), (2, 'b', 2.2), (3, 'c', 1.3)],
            index=pd.MultiIndex.from_arrays(
                [['c', 'b', 'a'], [3, 2, 1]],
                names=[None, 'foo']
            ),
            columns=['a', None, '__index_level_0__'],
        )
        if _pandas_api.uses_string_dtype():
            t = pa.Table.from_pandas(df, preserve_index=True)
        else:
            with pytest.warns(UserWarning):
                t = pa.Table.from_pandas(df, preserve_index=True)
        js = t.schema.pandas_metadata

        col1, col2, col3, idx0, foo = js['columns']

        assert col1['name'] == 'a'
        assert col1['name'] == col1['field_name']

        if _pandas_api.uses_string_dtype():
            assert np.isnan(col2['name'])
            assert col2['field_name'] == 'nan'
        else:
            assert col2['name'] is None
            assert col2['field_name'] == 'None'

        assert col3['name'] == '__index_level_0__'
        assert col3['name'] == col3['field_name']

        idx0_descr, foo_descr = js['index_columns']
        assert idx0_descr == '__index_level_0__'
        assert idx0['field_name'] == idx0_descr
        assert idx0['name'] is None

        assert foo_descr == 'foo'
        assert foo['field_name'] == foo_descr
        assert foo['name'] == foo_descr

    def test_categorical_column_index(self):
        df = pd.DataFrame(
            [(1, 'a', 2.0), (2, 'b', 3.0), (3, 'c', 4.0)],
            columns=pd.Index(list('def'), dtype='category')
        )
        t = pa.Table.from_pandas(df, preserve_index=True)
        js = t.schema.pandas_metadata

        column_indexes, = js['column_indexes']
        assert column_indexes['name'] is None
        assert column_indexes['pandas_type'] == 'categorical'
        assert column_indexes['numpy_type'] == 'int8'

        md = column_indexes['metadata']
        assert md['num_categories'] == 3
        assert md['ordered'] is False

    def test_string_column_index(self):
        df = pd.DataFrame(
            [(1, 'a', 2.0), (2, 'b', 3.0), (3, 'c', 4.0)],
            columns=pd.Index(list('def'), name='stringz')
        )
        t = pa.Table.from_pandas(df, preserve_index=True)
        js = t.schema.pandas_metadata

        column_indexes, = js['column_indexes']
        assert column_indexes['name'] == 'stringz'
        assert column_indexes['name'] == column_indexes['field_name']
        assert column_indexes['numpy_type'] == (
            'str' if _pandas_api.uses_string_dtype() else 'object'
        )
        assert column_indexes['pandas_type'] == 'unicode'

        md = column_indexes['metadata']

        assert len(md) == 1
        assert md['encoding'] == 'UTF-8'

    @pytest.mark.parametrize('unit', ['us', 'ns'])
    def test_datetimetz_column_index(self, unit):
        ext_kwargs = {}
        if Version(pd.__version__) >= Version("2.0.0"):
            # unit argument not supported on date_range for pandas < 2.0.0
            ext_kwargs = {'unit': unit}
        df = pd.DataFrame(
            [(1, 'a', 2.0), (2, 'b', 3.0), (3, 'c', 4.0)],
            columns=pd.date_range(
                start='2017-01-01', periods=3, tz='America/New_York', **ext_kwargs
            )
        )
        t = pa.Table.from_pandas(df, preserve_index=True)
        js = t.schema.pandas_metadata

        column_indexes, = js['column_indexes']
        assert column_indexes['name'] is None
        assert column_indexes['pandas_type'] == 'datetimetz'
        if ext_kwargs:
            assert column_indexes['numpy_type'] == f'datetime64[{unit}]'
        else:
            assert column_indexes['numpy_type'] == 'datetime64[ns]'

        md = column_indexes['metadata']
        assert md['timezone'] == 'America/New_York'

    def test_datetimetz_row_index(self):
        df = pd.DataFrame({
            'a': pd.date_range(
                start='2017-01-01', periods=3, tz='America/New_York'
            )
        })
        df = df.set_index('a')

        _check_pandas_roundtrip(df, preserve_index=True)

    def test_categorical_row_index(self):
        df = pd.DataFrame({'a': [1, 2, 3], 'b': [1, 2, 3]})
        df['a'] = df.a.astype('category')
        df = df.set_index('a')

        _check_pandas_roundtrip(df, preserve_index=True)

    def test_duplicate_column_names_does_not_crash(self):
        df = pd.DataFrame([(1, 'a'), (2, 'b')], columns=list('aa'))
        with pytest.raises(ValueError):
            pa.Table.from_pandas(df)

    def test_dictionary_indices_boundscheck(self):
        # ARROW-1658. No validation of indices leads to segfaults in pandas
        indices = [[0, 1], [0, -1]]

        for inds in indices:
            arr = pa.DictionaryArray.from_arrays(inds, ['a'], safe=False)
            batch = pa.RecordBatch.from_arrays([arr], ['foo'])
            table = pa.Table.from_batches([batch, batch, batch])

            with pytest.raises(IndexError):
                arr.to_pandas()

            with pytest.raises(IndexError):
                table.to_pandas()

    def test_unicode_with_unicode_column_and_index(self):
        df = pd.DataFrame({'あ': ['い']}, index=['う'])

        _check_pandas_roundtrip(df, preserve_index=True)

    def test_mixed_column_names(self):
        # mixed type column names are not reconstructed exactly
        df = pd.DataFrame({'a': [1, 2], 'b': [3, 4]})

        for cols in [['あ', b'a'], [1, '2'], [1, 1.5]]:
            df.columns = pd.Index(cols, dtype=object)

            # assert that the from_pandas raises the warning
            with pytest.warns(UserWarning):
                pa.Table.from_pandas(df)

            expected = df.copy()
            expected.columns = df.columns.values.astype(str)
            with pytest.warns(UserWarning):
                _check_pandas_roundtrip(df, expected=expected,
                                        preserve_index=True)

    def test_binary_column_name(self):
        if Version("2.0.0") <= Version(pd.__version__) < Version("3.0.0"):
            # TODO: regression in pandas, hopefully fixed in next version
            # https://issues.apache.org/jira/browse/ARROW-18394
            # https://github.com/pandas-dev/pandas/issues/50127
            pytest.skip("Regression in pandas 2.0.0")
        column_data = ['い']
        key = 'あ'.encode()
        data = {key: column_data}
        df = pd.DataFrame(data)

        # we can't use _check_pandas_roundtrip here because our metadata
        # is always decoded as utf8: even if binary goes in, utf8 comes out
        t = pa.Table.from_pandas(df, preserve_index=True)
        df2 = t.to_pandas()
        assert df.values[0] == df2.values[0]
        assert df.index.values[0] == df2.index.values[0]
        assert df.columns[0] == key

    def test_multiindex_duplicate_values(self):
        num_rows = 3
        numbers = list(range(num_rows))
        index = pd.MultiIndex.from_arrays(
            [['foo', 'foo', 'bar'], numbers],
            names=['foobar', 'some_numbers'],
        )

        df = pd.DataFrame({'numbers': numbers}, index=index)

        _check_pandas_roundtrip(df, preserve_index=True)

    def test_metadata_with_mixed_types(self):
        df = pd.DataFrame({'data': [b'some_bytes', 'some_unicode']})
        table = pa.Table.from_pandas(df)
        js = table.schema.pandas_metadata
        assert 'mixed' not in js
        data_column = js['columns'][0]
        assert data_column['pandas_type'] == 'bytes'
        assert data_column['numpy_type'] == 'object'

    def test_ignore_metadata(self):
        df = pd.DataFrame({'a': [1, 2, 3], 'b': ['foo', 'bar', 'baz']},
                          index=['one', 'two', 'three'])
        table = pa.Table.from_pandas(df)

        result = table.to_pandas(ignore_metadata=True)
        expected = (table.cast(table.schema.remove_metadata())
                    .to_pandas())

        tm.assert_frame_equal(result, expected)

    def test_list_metadata(self):
        df = pd.DataFrame({'data': [[1], [2, 3, 4], [5] * 7]})
        schema = pa.schema([pa.field('data', type=pa.list_(pa.int64()))])
        table = pa.Table.from_pandas(df, schema=schema)
        js = table.schema.pandas_metadata
        assert 'mixed' not in js
        data_column = js['columns'][0]
        assert data_column['pandas_type'] == 'list[int64]'
        assert data_column['numpy_type'] == 'object'

    def test_struct_metadata(self):
        df = pd.DataFrame({'dicts': [{'a': 1, 'b': 2}, {'a': 3, 'b': 4}]})
        table = pa.Table.from_pandas(df)
        pandas_metadata = table.schema.pandas_metadata
        assert pandas_metadata['columns'][0]['pandas_type'] == 'object'

    def test_decimal_metadata(self):
        expected = pd.DataFrame({
            'decimals': [
                decimal.Decimal('394092382910493.12341234678'),
                -decimal.Decimal('314292388910493.12343437128'),
            ]
        })
        table = pa.Table.from_pandas(expected)
        js = table.schema.pandas_metadata
        assert 'mixed' not in js
        data_column = js['columns'][0]
        assert data_column['pandas_type'] == 'decimal'
        assert data_column['numpy_type'] == 'object'
        assert data_column['metadata'] == {'precision': 26, 'scale': 11}

    @pytest.mark.parametrize('typ', [
        pa.decimal32,
        pa.decimal64,
        pa.decimal128,
        pa.decimal256,
    ])
    def test_decimal_other_bitwidts(self, typ):
        df = pd.DataFrame({'a': [decimal.Decimal('3.14')]})
        schema = pa.schema([pa.field('a', type=typ(4, 2))])
        table = pa.Table.from_pandas(df, schema=schema)
        col_meta = table.schema.pandas_metadata['columns'][0]
        assert col_meta['pandas_type'] == 'decimal'
        assert col_meta['metadata'] == {'precision': 4, 'scale': 2}

    def test_table_column_subset_metadata(self):
        # ARROW-1883
        # non-default index
        for index in [
                pd.Index(['a', 'b', 'c'], name='index'),
                pd.date_range("2017-01-01", periods=3, tz='Europe/Brussels')]:
            df = pd.DataFrame({'a': [1, 2, 3],
                               'b': [.1, .2, .3]}, index=index)
            table = pa.Table.from_pandas(df)

            table_subset = table.remove_column(1)
            result = table_subset.to_pandas()
            expected = df[['a']]
            if isinstance(df.index, pd.DatetimeIndex):
                df.index.freq = None
            tm.assert_frame_equal(result, expected)

            table_subset2 = table_subset.remove_column(1)
            result = table_subset2.to_pandas()
            tm.assert_frame_equal(result, df[['a']].reset_index(drop=True))

    def test_to_pandas_column_subset_multiindex(self):
        # ARROW-10122
        df = pd.DataFrame(
            {"first": list(range(5)),
             "second": list(range(5)),
             "value": np.arange(5)}
        )
        table = pa.Table.from_pandas(df.set_index(["first", "second"]))

        subset = table.select(["first", "value"])
        result = subset.to_pandas()
        expected = df[["first", "value"]].set_index("first")
        tm.assert_frame_equal(result, expected)

    def test_empty_list_metadata(self):
        # Create table with array of empty lists, forced to have type
        # list(string) in pyarrow
        c1 = [["test"], ["a", "b"], None]
        c2 = [[], [], []]
        arrays = OrderedDict([
            ('c1', pa.array(c1, type=pa.list_(pa.string()))),
            ('c2', pa.array(c2, type=pa.list_(pa.string()))),
        ])
        rb = pa.RecordBatch.from_arrays(
            list(arrays.values()),
            list(arrays.keys())
        )
        tbl = pa.Table.from_batches([rb])

        # First roundtrip changes schema, because pandas cannot preserve the
        # type of empty lists
        df = tbl.to_pandas()
        tbl2 = pa.Table.from_pandas(df)
        md2 = tbl2.schema.pandas_metadata

        # Second roundtrip
        df2 = tbl2.to_pandas()
        expected = pd.DataFrame(OrderedDict([('c1', c1), ('c2', c2)]))

        tm.assert_frame_equal(df2, expected)

        assert md2['columns'] == [
            {
                'name': 'c1',
                'field_name': 'c1',
                'metadata': None,
                'numpy_type': 'object',
                'pandas_type': 'list[unicode]',
            },
            {
                'name': 'c2',
                'field_name': 'c2',
                'metadata': None,
                'numpy_type': 'object',
                'pandas_type': 'list[empty]',
            }
        ]

    def test_metadata_pandas_version(self):
        df = pd.DataFrame({'a': [1, 2, 3], 'b': [1, 2, 3]})
        table = pa.Table.from_pandas(df)
        assert table.schema.pandas_metadata['pandas_version'] is not None

    def test_mismatch_metadata_schema(self):
        # ARROW-10511
        # It is possible that the metadata and actual schema is not fully
        # matching (eg no timezone information for tz-aware column)
        # -> to_pandas() conversion should not fail on that
        ext_kwargs = {}
        if Version(pd.__version__) >= Version("2.0.0"):
            # unit argument not supported on date_range for pandas < 2.0.0
            ext_kwargs = {'unit': 'ns'}
        df = pd.DataFrame({"datetime": pd.date_range(
            "2020-01-01", periods=3, **ext_kwargs)})

        # OPTION 1: casting after conversion
        table = pa.Table.from_pandas(df)
        # cast the "datetime" column to be tz-aware
        new_col = table["datetime"].cast(pa.timestamp('ns', tz="UTC"))
        new_table1 = table.set_column(
            0, pa.field("datetime", new_col.type), new_col
        )

        # OPTION 2: specify schema during conversion
        schema = pa.schema([("datetime", pa.timestamp('ns', tz="UTC"))])
        new_table2 = pa.Table.from_pandas(df, schema=schema)

        expected = df.copy()
        expected["datetime"] = expected["datetime"].dt.tz_localize("UTC")

        for new_table in [new_table1, new_table2]:
            # ensure the new table still has the pandas metadata
            assert new_table.schema.pandas_metadata is not None
            # convert to pandas
            result = new_table.to_pandas()
            tm.assert_frame_equal(result, expected)


class TestConvertPrimitiveTypes:
    """
    Conversion tests for primitive (e.g. numeric) types.
    """

    def test_float_no_nulls(self):
        data = {}
        fields = []
        dtypes = [('f2', pa.float16()),
                  ('f4', pa.float32()),
                  ('f8', pa.float64())]
        num_values = 100

        for numpy_dtype, arrow_dtype in dtypes:
            values = np.random.randn(num_values)
            data[numpy_dtype] = values.astype(numpy_dtype)
            fields.append(pa.field(numpy_dtype, arrow_dtype))

        df = pd.DataFrame(data)
        schema = pa.schema(fields)
        _check_pandas_roundtrip(df, expected_schema=schema)

    def test_float_nulls(self):
        num_values = 100

        null_mask = np.random.randint(0, 10, size=num_values) < 3
        dtypes = [('f2', pa.float16()),
                  ('f4', pa.float32()),
                  ('f8', pa.float64())]
        names = ['f2', 'f4', 'f8']
        expected_cols = []

        arrays = []
        fields = []
        for name, arrow_dtype in dtypes:
            values = np.random.randn(num_values).astype(name)

            arr = pa.array(values, from_pandas=True, mask=null_mask)
            arrays.append(arr)
            fields.append(pa.field(name, arrow_dtype))
            values[null_mask] = np.nan

            expected_cols.append(values)

        ex_frame = pd.DataFrame(dict(zip(names, expected_cols)),
                                columns=names)

        table = pa.Table.from_arrays(arrays, names)
        assert table.schema.equals(pa.schema(fields))
        result = table.to_pandas()
        tm.assert_frame_equal(result, ex_frame)

    def test_float_nulls_to_ints(self):
        # ARROW-2135
        df = pd.DataFrame({"a": [1.0, 2.0, np.nan]})
        schema = pa.schema([pa.field("a", pa.int16(), nullable=True)])
        table = pa.Table.from_pandas(df, schema=schema, safe=False)
        assert table[0].to_pylist() == [1, 2, None]
        tm.assert_frame_equal(df, table.to_pandas())

    def test_float_nulls_to_boolean(self):
        s = pd.Series([0.0, 1.0, 2.0, None, -3.0])
        expected = pd.Series([False, True, True, None, True])
        _check_array_roundtrip(s, expected=expected, type=pa.bool_())

    def test_series_from_pandas_false_respected(self):
        # Check that explicit from_pandas=False is respected
        s = pd.Series([0.0, np.nan])
        arr = pa.array(s, from_pandas=False)
        assert arr.null_count == 0
        assert np.isnan(arr[1].as_py())

    def test_integer_no_nulls(self):
        data = OrderedDict()
        fields = []

        numpy_dtypes = [
            ('i1', pa.int8()), ('i2', pa.int16()),
            ('i4', pa.int32()), ('i8', pa.int64()),
            ('u1', pa.uint8()), ('u2', pa.uint16()),
            ('u4', pa.uint32()), ('u8', pa.uint64()),
            ('longlong', pa.int64()), ('ulonglong', pa.uint64())
        ]
        num_values = 100

        for dtype, arrow_dtype in numpy_dtypes:
            info = np.iinfo(dtype)
            values = np.random.randint(max(info.min, np.iinfo(np.int_).min),
                                       min(info.max, np.iinfo(np.int_).max),
                                       size=num_values, dtype=dtype)
            data[dtype] = values.astype(dtype)
            fields.append(pa.field(dtype, arrow_dtype))

        df = pd.DataFrame(data)
        schema = pa.schema(fields)
        _check_pandas_roundtrip(df, expected_schema=schema)

    def test_all_integer_types(self):
        # Test all Numpy integer aliases
        data = OrderedDict()
        numpy_dtypes = ['i1', 'i2', 'i4', 'i8', 'u1', 'u2', 'u4', 'u8',
                        'byte', 'ubyte', 'short', 'ushort', 'intc', 'uintc',
                        'int_', 'uint', 'longlong', 'ulonglong']
        for dtype in numpy_dtypes:
            data[dtype] = np.arange(12, dtype=dtype)
        df = pd.DataFrame(data)
        _check_pandas_roundtrip(df)

        # Do the same with pa.array()
        # (for some reason, it doesn't use the same code paths at all)
        for np_arr in data.values():
            arr = pa.array(np_arr)
            assert arr.to_pylist() == np_arr.tolist()

    def test_integer_byteorder(self):
        # Byteswapped arrays are not supported yet
        int_dtypes = ['i1', 'i2', 'i4', 'i8', 'u1', 'u2', 'u4', 'u8']
        for dt in int_dtypes:
            for order in '=<>':
                data = np.array([1, 2, 42], dtype=order + dt)
                for np_arr in (data, data[::2]):
                    if data.dtype.isnative:
                        arr = pa.array(data)
                        assert arr.to_pylist() == data.tolist()
                    else:
                        with pytest.raises(NotImplementedError):
                            arr = pa.array(data)

    def test_integer_with_nulls(self):
        # pandas requires upcast to float dtype

        int_dtypes = ['i1', 'i2', 'i4', 'i8', 'u1', 'u2', 'u4', 'u8']
        num_values = 100

        null_mask = np.random.randint(0, 10, size=num_values) < 3

        expected_cols = []
        arrays = []
        for name in int_dtypes:
            values = np.random.randint(0, 100, size=num_values)

            arr = pa.array(values, mask=null_mask)
            arrays.append(arr)

            expected = values.astype('f8')
            expected[null_mask] = np.nan

            expected_cols.append(expected)

        ex_frame = pd.DataFrame(dict(zip(int_dtypes, expected_cols)),
                                columns=int_dtypes)

        table = pa.Table.from_arrays(arrays, int_dtypes)
        result = table.to_pandas()

        tm.assert_frame_equal(result, ex_frame)

    def test_array_from_pandas_type_cast(self):
        arr = np.arange(10, dtype='int64')

        target_type = pa.int8()

        result = pa.array(arr, type=target_type)
        expected = pa.array(arr.astype('int8'))
        assert result.equals(expected)

    def test_boolean_no_nulls(self):
        num_values = 100

        np.random.seed(0)

        df = pd.DataFrame({'bools': np.random.randn(num_values) > 0})
        field = pa.field('bools', pa.bool_())
        schema = pa.schema([field])
        _check_pandas_roundtrip(df, expected_schema=schema)

    def test_boolean_nulls(self):
        # pandas requires upcast to object dtype
        num_values = 100
        np.random.seed(0)

        mask = np.random.randint(0, 10, size=num_values) < 3
        values = np.random.randint(0, 10, size=num_values) < 5

        arr = pa.array(values, mask=mask)

        expected = values.astype(object)
        expected[mask] = None

        field = pa.field('bools', pa.bool_())
        schema = pa.schema([field])
        ex_frame = pd.DataFrame({'bools': expected})

        table = pa.Table.from_arrays([arr], ['bools'])
        assert table.schema.equals(schema)
        result = table.to_pandas()

        tm.assert_frame_equal(result, ex_frame)

    def test_boolean_to_int(self):
        # test from dtype=bool
        s = pd.Series([True, True, False, True, True] * 2)
        expected = pd.Series([1, 1, 0, 1, 1] * 2)
        _check_array_roundtrip(s, expected=expected, type=pa.int64())

    def test_boolean_objects_to_int(self):
        # test from dtype=object
        s = pd.Series([True, True, False, True, True] * 2, dtype=object)
        expected = pd.Series([1, 1, 0, 1, 1] * 2)
        expected_msg = 'Expected integer, got bool'
        with pytest.raises(pa.ArrowTypeError, match=expected_msg):
            _check_array_roundtrip(s, expected=expected, type=pa.int64())

    def test_boolean_nulls_to_float(self):
        # test from dtype=object
        s = pd.Series([True, True, False, None, True] * 2)
        expected = pd.Series([1.0, 1.0, 0.0, None, 1.0] * 2)
        _check_array_roundtrip(s, expected=expected, type=pa.float64())

    def test_boolean_multiple_columns(self):
        # ARROW-6325 (multiple columns resulting in strided conversion)
        df = pd.DataFrame(np.ones((3, 2), dtype='bool'), columns=['a', 'b'])
        _check_pandas_roundtrip(df)

    def test_float_object_nulls(self):
        arr = np.array([None, 1.5, np.float64(3.5)] * 5, dtype=object)
        df = pd.DataFrame({'floats': arr})
        expected = pd.DataFrame({'floats': pd.to_numeric(arr)})
        field = pa.field('floats', pa.float64())
        schema = pa.schema([field])
        _check_pandas_roundtrip(df, expected=expected,
                                expected_schema=schema)

    def test_float_with_null_as_integer(self):
        # ARROW-2298
        s = pd.Series([np.nan, 1., 2., np.nan])

        types = [pa.int8(), pa.int16(), pa.int32(), pa.int64(),
                 pa.uint8(), pa.uint16(), pa.uint32(), pa.uint64()]
        for ty in types:
            result = pa.array(s, type=ty)
            expected = pa.array([None, 1, 2, None], type=ty)
            assert result.equals(expected)

            df = pd.DataFrame({'has_nulls': s})
            schema = pa.schema([pa.field('has_nulls', ty)])
            result = pa.Table.from_pandas(df, schema=schema,
                                          preserve_index=False)
            assert result[0].chunk(0).equals(expected)

    def test_int_object_nulls(self):
        arr = np.array([None, 1, np.int64(3)] * 5, dtype=object)
        df = pd.DataFrame({'ints': arr})
        expected = pd.DataFrame({'ints': pd.to_numeric(arr)})
        field = pa.field('ints', pa.int64())
        schema = pa.schema([field])
        _check_pandas_roundtrip(df, expected=expected,
                                expected_schema=schema)

    def test_boolean_object_nulls(self):
        arr = np.array([False, None, True] * 100, dtype=object)
        df = pd.DataFrame({'bools': arr})
        field = pa.field('bools', pa.bool_())
        schema = pa.schema([field])
        _check_pandas_roundtrip(df, expected_schema=schema)

    def test_all_nulls_cast_numeric(self):
        arr = np.array([None], dtype=object)

        def _check_type(t):
            a2 = pa.array(arr, type=t)
            assert a2.type == t
            assert a2[0].as_py() is None

        _check_type(pa.int32())
        _check_type(pa.float64())

    def test_half_floats_from_numpy(self):
        arr = np.array([1.5, np.nan], dtype=np.float16)
        a = pa.array(arr, type=pa.float16())
        x, y = a.to_pylist()
        assert isinstance(x, float)
        assert x == 1.5
        assert isinstance(y, float)
        assert np.isnan(y)

        a = pa.array(arr, type=pa.float16(), from_pandas=True)
        x, y = a.to_pylist()
        assert isinstance(x, float)
        assert x == 1.5
        assert y is None


@pytest.mark.parametrize('dtype',
                         ['i1', 'i2', 'i4', 'i8', 'u1', 'u2', 'u4', 'u8'])
def test_array_integer_object_nulls_option(dtype):
    num_values = 100

    null_mask = np.random.randint(0, 10, size=num_values) < 3
    values = np.random.randint(0, 100, size=num_values, dtype=dtype)

    array = pa.array(values, mask=null_mask)

    if null_mask.any():
        expected = values.astype('O')
        expected[null_mask] = None
    else:
        expected = values

    result = array.to_pandas(integer_object_nulls=True)

    np.testing.assert_equal(result, expected)


@pytest.mark.parametrize('dtype',
                         ['i1', 'i2', 'i4', 'i8', 'u1', 'u2', 'u4', 'u8'])
def test_table_integer_object_nulls_option(dtype):
    num_values = 100

    null_mask = np.random.randint(0, 10, size=num_values) < 3
    values = np.random.randint(0, 100, size=num_values, dtype=dtype)

    array = pa.array(values, mask=null_mask)

    if null_mask.any():
        expected = values.astype('O')
        expected[null_mask] = None
    else:
        expected = values

    expected = pd.DataFrame({dtype: expected})

    table = pa.Table.from_arrays([array], [dtype])
    result = table.to_pandas(integer_object_nulls=True)

    tm.assert_frame_equal(result, expected)


class TestConvertDateTimeLikeTypes:
    """
    Conversion tests for datetime- and timestamp-like types (date64, etc.).
    """

    def test_timestamps_notimezone_no_nulls(self):
        df = pd.DataFrame({
            'datetime64': np.array([
                '2007-07-13T01:23:34.123456789',
                '2006-01-13T12:34:56.432539784',
                '2010-08-13T05:46:57.437699912'],
                dtype='datetime64[ns]')
        })
        field = pa.field('datetime64', pa.timestamp('ns'))
        schema = pa.schema([field])
        _check_pandas_roundtrip(
            df,
            expected_schema=schema,
        )

    def test_timestamps_notimezone_nulls(self):
        df = pd.DataFrame({
            'datetime64': np.array([
                '2007-07-13T01:23:34.123456789',
                None,
                '2010-08-13T05:46:57.437699912'],
                dtype='datetime64[ns]')
        })
        field = pa.field('datetime64', pa.timestamp('ns'))
        schema = pa.schema([field])
        _check_pandas_roundtrip(
            df,
            expected_schema=schema,
        )

    @pytest.mark.parametrize('unit', ['s', 'ms', 'us', 'ns'])
    def test_timestamps_with_timezone(self, unit):
        if Version(pd.__version__) < Version("2.0.0") and unit != 'ns':
            pytest.skip("pandas < 2.0 only supports nanosecond datetime64")
        df = pd.DataFrame({
            'datetime64': np.array([
                '2007-07-13T01:23:34.123',
                '2006-01-13T12:34:56.432',
                '2010-08-13T05:46:57.437'],
                dtype=f'datetime64[{unit}]')
        })
        df['datetime64'] = df['datetime64'].dt.tz_localize('US/Eastern')
        _check_pandas_roundtrip(df)

        _check_series_roundtrip(df['datetime64'])

        # drop-in a null
        df = pd.DataFrame({
            'datetime64': np.array([
                '2007-07-13T01:23:34.123456789',
                None,
                '2006-01-13T12:34:56.432539784',
                '2010-08-13T05:46:57.437699912'],
                dtype=f'datetime64[{unit}]')
        })
        df['datetime64'] = df['datetime64'].dt.tz_localize('US/Eastern')

        _check_pandas_roundtrip(df)

    def test_python_datetime(self):
        # ARROW-2106
        date_array = [datetime.today() + timedelta(days=x) for x in range(10)]
        df = pd.DataFrame({
            'datetime': pd.Series(date_array, dtype=object)
        })

        table = pa.Table.from_pandas(df)
        assert isinstance(table[0].chunk(0), pa.TimestampArray)

        result = table.to_pandas()
        # Pandas v2 defaults to [ns], but Arrow defaults to [us] time units
        # so we need to cast the pandas dtype. Pandas v1 will always silently
        # coerce to [ns] due to lack of non-[ns] support.
        expected_df = pd.DataFrame({
            'datetime': pd.Series(date_array, dtype='datetime64[us]')
        })
        tm.assert_frame_equal(expected_df, result)

    def test_python_datetime_with_pytz_tzinfo(self):
        pytz = pytest.importorskip("pytz")

        for tz in [pytz.utc, pytz.timezone('US/Eastern'), pytz.FixedOffset(1)]:
            values = [datetime(2018, 1, 1, 12, 23, 45, tzinfo=tz)]
            df = pd.DataFrame({'datetime': values})
            _check_pandas_roundtrip(df)

    @h.given(st.none() | past.timezones)
    @h.settings(deadline=None)
    def test_python_datetime_with_pytz_timezone(self, tz):
        if str(tz) in ["build/etc/localtime", "Factory"]:
            pytest.skip("Localtime timezone not supported")
        values = [datetime(2018, 1, 1, 12, 23, 45, tzinfo=tz)]
        df = pd.DataFrame({'datetime': values})
        _check_pandas_roundtrip(df, check_dtype=False)

    def test_python_datetime_with_timezone_tzinfo(self):
        pytz = pytest.importorskip("pytz")
        from datetime import timezone

        values = [datetime(2018, 1, 1, 12, 23, 45, tzinfo=timezone.utc)]
        # also test with index to ensure both paths roundtrip (ARROW-9962)
        df = pd.DataFrame({'datetime': values}, index=values)
        _check_pandas_roundtrip(df, preserve_index=True)

        # datetime.timezone is going to be pytz.FixedOffset
        hours = 1
        tz_timezone = timezone(timedelta(hours=hours))
        tz_pytz = pytz.FixedOffset(hours * 60)
        values = [datetime(2018, 1, 1, 12, 23, 45, tzinfo=tz_timezone)]
        values_exp = [datetime(2018, 1, 1, 12, 23, 45, tzinfo=tz_pytz)]
        df = pd.DataFrame({'datetime': values}, index=values)
        df_exp = pd.DataFrame({'datetime': values_exp}, index=values_exp)
        _check_pandas_roundtrip(df, expected=df_exp, preserve_index=True)

    def test_python_datetime_subclass(self):

        class MyDatetime(datetime):
            # see https://github.com/pandas-dev/pandas/issues/21142
            nanosecond = 0.0

        date_array = [MyDatetime(2000, 1, 1, 1, 1, 1)]
        df = pd.DataFrame({"datetime": pd.Series(date_array, dtype=object)})

        table = pa.Table.from_pandas(df)
        assert isinstance(table[0].chunk(0), pa.TimestampArray)

        result = table.to_pandas()

        # Pandas v2 defaults to [ns], but Arrow defaults to [us] time units
        # so we need to cast the pandas dtype. Pandas v1 will always silently
        # coerce to [ns] due to lack of non-[ns] support.
        expected_df = pd.DataFrame(
            {"datetime": pd.Series(date_array, dtype='datetime64[us]')})

        # https://github.com/pandas-dev/pandas/issues/21142
        expected_df["datetime"] = pd.to_datetime(expected_df["datetime"])

        tm.assert_frame_equal(expected_df, result)

    def test_python_date_subclass(self):

        class MyDate(date):
            pass

        date_array = [MyDate(2000, 1, 1)]
        df = pd.DataFrame({"date": pd.Series(date_array, dtype=object)})

        table = pa.Table.from_pandas(df)
        assert isinstance(table[0].chunk(0), pa.Date32Array)

        result = table.to_pandas()
        expected_df = pd.DataFrame(
            {"date": np.array([date(2000, 1, 1)], dtype=object)}
        )
        tm.assert_frame_equal(expected_df, result)

    def test_datetime64_to_date32(self):
        # ARROW-1718
        arr = pa.array([date(2017, 10, 23), None])
        c = pa.chunked_array([arr])
        s = c.to_pandas()

        arr2 = pa.Array.from_pandas(s, type=pa.date32())

        assert arr2.equals(arr.cast('date32'))

    @pytest.mark.parametrize('mask', [
        None,
        [True, False, False, True, False, False],
    ])
    def test_pandas_datetime_to_date64(self, mask):
        if mask:
            mask = np.array(mask)
        s = pd.to_datetime([
            '2018-05-10T00:00:00',
            '2018-05-11T00:00:00',
            '2018-05-12T00:00:00',
            '2018-05-10T10:24:01',
            '2018-05-11T10:24:01',
            '2018-05-12T10:24:01',
        ])
        arr = pa.Array.from_pandas(s, type=pa.date64(), mask=mask)

        data = np.array([
            date(2018, 5, 10),
            date(2018, 5, 11),
            date(2018, 5, 12),
            date(2018, 5, 10),
            date(2018, 5, 11),
            date(2018, 5, 12),
        ])
        expected = pa.array(data, mask=mask, type=pa.date64())

        assert arr.equals(expected)

    @pytest.mark.parametrize("coerce_to_ns,expected_dtype",
                             [(False, 'datetime64[ms]'),
                              (True, 'datetime64[ns]')])
    def test_array_types_date_as_object(self, coerce_to_ns, expected_dtype):
        data = [date(2000, 1, 1),
                None,
                date(1970, 1, 1),
                date(2040, 2, 26)]
        expected_days = np.array(['2000-01-01', None, '1970-01-01',
                                  '2040-02-26'], dtype='datetime64[D]')

        if Version(pd.__version__) < Version("2.0.0"):
            # ARROW-3789: Coerce date/timestamp types to datetime64[ns]
            expected_dtype = 'datetime64[ns]'

        expected = np.array(['2000-01-01', None, '1970-01-01',
                             '2040-02-26'], dtype=expected_dtype)

        objects = [pa.array(data),
                   pa.chunked_array([data])]

        for obj in objects:
            result = obj.to_pandas(coerce_temporal_nanoseconds=coerce_to_ns)
            expected_obj = expected_days.astype(object)
            assert result.dtype == expected_obj.dtype
            npt.assert_array_equal(result, expected_obj)

            result = obj.to_pandas(date_as_object=False,
                                   coerce_temporal_nanoseconds=coerce_to_ns)
            assert result.dtype == expected.dtype
            npt.assert_array_equal(result, expected)

    @pytest.mark.parametrize("coerce_to_ns,expected_type",
                             [(False, 'datetime64[ms]'),
                              (True, 'datetime64[ns]')])
    def test_table_convert_date_as_object(self, coerce_to_ns, expected_type):
        df = pd.DataFrame({
            'date': [date(2000, 1, 1),
                     None,
                     date(1970, 1, 1),
                     date(2040, 2, 26)]})

        table = pa.Table.from_pandas(df, preserve_index=False)

        df_datetime = table.to_pandas(date_as_object=False,
                                      coerce_temporal_nanoseconds=coerce_to_ns)
        df_object = table.to_pandas()

        tm.assert_frame_equal(df.astype(expected_type), df_datetime,
                              check_dtype=True)
        tm.assert_frame_equal(df, df_object, check_dtype=True)

    @pytest.mark.parametrize("arrow_type",
                             [pa.date32(), pa.date64(), pa.timestamp('s'),
                              pa.timestamp('ms'), pa.timestamp('us'),
                              pa.timestamp('ns'), pa.timestamp('s', 'UTC'),
                              pa.timestamp('ms', 'UTC'), pa.timestamp('us', 'UTC'),
                              pa.timestamp('ns', 'UTC')])
    def test_array_coerce_temporal_nanoseconds(self, arrow_type):
        data = [date(2000, 1, 1), datetime(2001, 1, 1)]
        expected = pd.Series(data)
        arr = pa.array(data).cast(arrow_type)
        result = arr.to_pandas(
            coerce_temporal_nanoseconds=True, date_as_object=False)
        expected_tz = None
        if hasattr(arrow_type, 'tz') and arrow_type.tz is not None:
            expected_tz = 'UTC'
        expected_type = pa.timestamp('ns', expected_tz).to_pandas_dtype()
        tm.assert_series_equal(result, expected.astype(expected_type))

    @pytest.mark.parametrize("arrow_type",
                             [pa.date32(), pa.date64(), pa.timestamp('s'),
                              pa.timestamp('ms'), pa.timestamp('us'),
                              pa.timestamp('ns'), pa.timestamp('s', 'UTC'),
                              pa.timestamp('ms', 'UTC'), pa.timestamp('us', 'UTC'),
                              pa.timestamp('ns', 'UTC')])
    def test_table_coerce_temporal_nanoseconds(self, arrow_type):
        data = [date(2000, 1, 1), datetime(2001, 1, 1)]
        schema = pa.schema([pa.field('date', arrow_type)])
        expected_df = pd.DataFrame({'date': data})
        table = pa.table([pa.array(data)], schema=schema)
        result_df = table.to_pandas(
            coerce_temporal_nanoseconds=True, date_as_object=False)
        expected_tz = None
        if hasattr(arrow_type, 'tz') and arrow_type.tz is not None:
            expected_tz = 'UTC'
        expected_type = pa.timestamp('ns', expected_tz).to_pandas_dtype()
        tm.assert_frame_equal(result_df, expected_df.astype(expected_type))

    def test_date_infer(self):
        df = pd.DataFrame({
            'date': [date(2000, 1, 1),
                     None,
                     date(1970, 1, 1),
                     date(2040, 2, 26)]})
        table = pa.Table.from_pandas(df, preserve_index=False)
        field = pa.field('date', pa.date32())

        # schema's metadata is generated by from_pandas conversion
        expected_schema = pa.schema([field], metadata=table.schema.metadata)
        assert table.schema.equals(expected_schema)

        result = table.to_pandas()
        tm.assert_frame_equal(result, df)

    def test_date_mask(self):
        arr = np.array([date(2017, 4, 3), date(2017, 4, 4)],
                       dtype='datetime64[D]')
        mask = [True, False]
        result = pa.array(arr, mask=np.array(mask))
        expected = np.array([None, date(2017, 4, 4)], dtype='datetime64[D]')
        expected = pa.array(expected, from_pandas=True)
        assert expected.equals(result)

    def test_date_objects_typed(self):
        arr = np.array([
            date(2017, 4, 3),
            None,
            date(2017, 4, 4),
            date(2017, 4, 5)], dtype=object)

        arr_i4 = np.array([17259, -1, 17260, 17261], dtype='int32')
        arr_i8 = arr_i4.astype('int64') * 86400000
        mask = np.array([False, True, False, False])

        t32 = pa.date32()
        t64 = pa.date64()

        a32 = pa.array(arr, type=t32)
        a64 = pa.array(arr, type=t64)

        a32_expected = pa.array(arr_i4, mask=mask, type=t32)
        a64_expected = pa.array(arr_i8, mask=mask, type=t64)

        assert a32.equals(a32_expected)
        assert a64.equals(a64_expected)

        # Test converting back to pandas
        colnames = ['date32', 'date64']
        table = pa.Table.from_arrays([a32, a64], colnames)

        ex_values = (np.array(['2017-04-03', '2017-04-04', '2017-04-04',
                               '2017-04-05'],
                              dtype='datetime64[D]'))
        ex_values[1] = pd.NaT.value

        # date32 and date64 convert to [ms] in pandas v2, but
        # in pandas v1 they are silently coerced to [ns]
        ex_datetime64ms = ex_values.astype('datetime64[ms]')
        expected_pandas = pd.DataFrame({'date32': ex_datetime64ms,
                                        'date64': ex_datetime64ms},
                                       columns=colnames)
        table_pandas = table.to_pandas(date_as_object=False)
        tm.assert_frame_equal(table_pandas, expected_pandas)

        table_pandas_objects = table.to_pandas()
        ex_objects = ex_values.astype('object')
        expected_pandas_objects = pd.DataFrame({'date32': ex_objects,
                                                'date64': ex_objects},
                                               columns=colnames)
        tm.assert_frame_equal(table_pandas_objects,
                              expected_pandas_objects)

    def test_pandas_null_values(self):
        # ARROW-842
        pd_NA = getattr(pd, 'NA', None)
        values = np.array([datetime(2000, 1, 1), pd.NaT, pd_NA], dtype=object)
        values_with_none = np.array([datetime(2000, 1, 1), None, None],
                                    dtype=object)
        result = pa.array(values, from_pandas=True)
        expected = pa.array(values_with_none, from_pandas=True)
        assert result.equals(expected)
        assert result.null_count == 2

        # ARROW-9407
        assert pa.array([pd.NaT], from_pandas=True).type == pa.null()
        assert pa.array([pd_NA], from_pandas=True).type == pa.null()

    def test_dates_from_integers(self):
        t1 = pa.date32()
        t2 = pa.date64()

        arr = np.array([17259, 17260, 17261], dtype='int32')
        arr2 = arr.astype('int64') * 86400000

        a1 = pa.array(arr, type=t1)
        a2 = pa.array(arr2, type=t2)

        expected = date(2017, 4, 3)
        assert a1[0].as_py() == expected
        assert a2[0].as_py() == expected

    def test_pytime_from_pandas(self):
        pytimes = [time(1, 2, 3, 1356),
                   time(4, 5, 6, 1356)]

        # microseconds
        t1 = pa.time64('us')

        aobjs = np.array(pytimes + [None], dtype=object)
        parr = pa.array(aobjs)
        assert parr.type == t1
        assert parr[0].as_py() == pytimes[0]
        assert parr[1].as_py() == pytimes[1]
        assert parr[2].as_py() is None

        # DataFrame
        df = pd.DataFrame({'times': aobjs})
        batch = pa.RecordBatch.from_pandas(df)
        assert batch[0].equals(parr)

        # Test ndarray of int64 values
        arr = np.array([_pytime_to_micros(v) for v in pytimes],
                       dtype='int64')

        a1 = pa.array(arr, type=pa.time64('us'))
        assert a1[0].as_py() == pytimes[0]

        a2 = pa.array(arr * 1000, type=pa.time64('ns'))
        assert a2[0].as_py() == pytimes[0]

        a3 = pa.array((arr / 1000).astype('i4'),
                      type=pa.time32('ms'))
        assert a3[0].as_py() == pytimes[0].replace(microsecond=1000)

        a4 = pa.array((arr / 1000000).astype('i4'),
                      type=pa.time32('s'))
        assert a4[0].as_py() == pytimes[0].replace(microsecond=0)

    def test_arrow_time_to_pandas(self):
        pytimes = [time(1, 2, 3, 1356),
                   time(4, 5, 6, 1356),
                   time(0, 0, 0)]

        expected = np.array(pytimes[:2] + [None])
        expected_ms = np.array([x.replace(microsecond=1000)
                                for x in pytimes[:2]] +
                               [None])
        expected_s = np.array([x.replace(microsecond=0)
                               for x in pytimes[:2]] +
                              [None])

        arr = np.array([_pytime_to_micros(v) for v in pytimes],
                       dtype='int64')
        arr = np.array([_pytime_to_micros(v) for v in pytimes],
                       dtype='int64')

        null_mask = np.array([False, False, True], dtype=bool)

        a1 = pa.array(arr, mask=null_mask, type=pa.time64('us'))
        a2 = pa.array(arr * 1000, mask=null_mask,
                      type=pa.time64('ns'))

        a3 = pa.array((arr / 1000).astype('i4'), mask=null_mask,
                      type=pa.time32('ms'))
        a4 = pa.array((arr / 1000000).astype('i4'), mask=null_mask,
                      type=pa.time32('s'))

        names = ['time64[us]', 'time64[ns]', 'time32[ms]', 'time32[s]']
        batch = pa.RecordBatch.from_arrays([a1, a2, a3, a4], names)

        for arr, expected_values in [(a1, expected),
                                     (a2, expected),
                                     (a3, expected_ms),
                                     (a4, expected_s)]:
            result_pandas = arr.to_pandas()
            assert (result_pandas.values == expected_values).all()

        df = batch.to_pandas()
        expected_df = pd.DataFrame({'time64[us]': expected,
                                    'time64[ns]': expected,
                                    'time32[ms]': expected_ms,
                                    'time32[s]': expected_s},
                                   columns=names)

        tm.assert_frame_equal(df, expected_df)

    def test_numpy_datetime64_columns(self):
        datetime64_ns = np.array([
            '2007-07-13T01:23:34.123456789',
            None,
            '2006-01-13T12:34:56.432539784',
            '2010-08-13T05:46:57.437699912'],
            dtype='datetime64[ns]')
        _check_array_from_pandas_roundtrip(datetime64_ns)

        datetime64_us = np.array([
            '2007-07-13T01:23:34.123456',
            None,
            '2006-01-13T12:34:56.432539',
            '2010-08-13T05:46:57.437699'],
            dtype='datetime64[us]')
        _check_array_from_pandas_roundtrip(datetime64_us)

        datetime64_ms = np.array([
            '2007-07-13T01:23:34.123',
            None,
            '2006-01-13T12:34:56.432',
            '2010-08-13T05:46:57.437'],
            dtype='datetime64[ms]')
        _check_array_from_pandas_roundtrip(datetime64_ms)

        datetime64_s = np.array([
            '2007-07-13T01:23:34',
            None,
            '2006-01-13T12:34:56',
            '2010-08-13T05:46:57'],
            dtype='datetime64[s]')
        _check_array_from_pandas_roundtrip(datetime64_s)

    def test_timestamp_to_pandas_coerces_to_ns(self):
        # non-ns timestamp gets cast to ns on conversion to pandas
        if Version(pd.__version__) >= Version("2.0.0"):
            pytest.skip("pandas >= 2.0 supports non-nanosecond datetime64")

        arr = pa.array([1, 2, 3], pa.timestamp('ms'))
        expected = pd.Series(pd.to_datetime([1, 2, 3], unit='ms'))
        s = arr.to_pandas()
        tm.assert_series_equal(s, expected)
        arr = pa.chunked_array([arr])
        s = arr.to_pandas()
        tm.assert_series_equal(s, expected)

    def test_timestamp_to_pandas_out_of_bounds(self):
        # ARROW-7758 check for out of bounds timestamps for non-ns timestamps
        # that end up getting coerced into ns timestamps.

        for unit in ['s', 'ms', 'us']:
            for tz in [None, 'America/New_York']:
                arr = pa.array([datetime(1, 1, 1)], pa.timestamp(unit, tz=tz))
                table = pa.table({'a': arr})

                msg = "would result in out of bounds timestamp"
                with pytest.raises(ValueError, match=msg):
                    arr.to_pandas(coerce_temporal_nanoseconds=True)

                with pytest.raises(ValueError, match=msg):
                    table.to_pandas(coerce_temporal_nanoseconds=True)

                with pytest.raises(ValueError, match=msg):
                    # chunked array
                    table.column('a').to_pandas(coerce_temporal_nanoseconds=True)

                # just ensure those don't give an error, but do not
                # check actual garbage output
                arr.to_pandas(safe=False, coerce_temporal_nanoseconds=True)
                table.to_pandas(safe=False, coerce_temporal_nanoseconds=True)
                table.column('a').to_pandas(
                    safe=False, coerce_temporal_nanoseconds=True)

    def test_timestamp_to_pandas_empty_chunked(self):
        # ARROW-7907 table with chunked array with 0 chunks
        table = pa.table({'a': pa.chunked_array([], type=pa.timestamp('us'))})
        result = table.to_pandas()
        expected = pd.DataFrame({'a': pd.Series([], dtype="datetime64[us]")})
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize('dtype', [pa.date32(), pa.date64()])
    def test_numpy_datetime64_day_unit(self, dtype):
        datetime64_d = np.array([
            '2007-07-13',
            None,
            '2006-01-15',
            '2010-08-19'],
            dtype='datetime64[D]')
        _check_array_from_pandas_roundtrip(datetime64_d, type=dtype)

    def test_array_from_pandas_date_with_mask(self):
        m = np.array([True, False, True])
        data = pd.Series([
            date(1990, 1, 1),
            date(1991, 1, 1),
            date(1992, 1, 1)
        ])

        result = pa.Array.from_pandas(data, mask=m)

        expected = pd.Series([None, date(1991, 1, 1), None])
        assert pa.Array.from_pandas(expected).equals(result)

    @pytest.mark.skipif(
        np is not None and Version('1.16.0') <= Version(
            np.__version__) < Version('1.16.1'),
        reason='Until numpy/numpy#12745 is resolved')
    def test_fixed_offset_timezone(self):
        df = pd.DataFrame({
            'a': [
                pd.Timestamp('2012-11-11 00:00:00+01:00'),
                pd.NaT
            ]
        })
        # 'check_dtype=False' because pandas >= 2 uses datetime.timezone
        # instead of pytz.FixedOffset, and thus the dtype is not exactly
        # identical (pyarrow still defaults to pytz)
        # TODO remove if https://github.com/apache/arrow/issues/15047 is fixed
        _check_pandas_roundtrip(df, check_dtype=False)

    @pytest.mark.parametrize("unit", ['s', 'ms', 'us', 'ns'])
    def test_timedeltas_no_nulls(self, unit):
        if Version(pd.__version__) < Version("2.0.0"):
            unit = 'ns'
        df = pd.DataFrame({
            'timedelta64': np.array([0, 3600000000000, 7200000000000],
                                    dtype=f'timedelta64[{unit}]')
        })
        field = pa.field('timedelta64', pa.duration(unit))
        schema = pa.schema([field])
        _check_pandas_roundtrip(
            df,
            expected_schema=schema,
        )

    @pytest.mark.parametrize("unit", ['s', 'ms', 'us', 'ns'])
    def test_timedeltas_nulls(self, unit):
        if Version(pd.__version__) < Version("2.0.0"):
            unit = 'ns'
        df = pd.DataFrame({
            'timedelta64': np.array([0, None, 7200000000000],
                                    dtype=f'timedelta64[{unit}]')
        })
        field = pa.field('timedelta64', pa.duration(unit))
        schema = pa.schema([field])
        _check_pandas_roundtrip(
            df,
            expected_schema=schema,
        )

    def test_month_day_nano_interval(self):
        from pandas.tseries.offsets import DateOffset
        df = pd.DataFrame({
            'date_offset': [None,
                            DateOffset(days=3600, months=3600, microseconds=3,
                                       nanoseconds=600)]
        })
        schema = pa.schema([('date_offset', pa.month_day_nano_interval())])
        _check_pandas_roundtrip(
            df,
            expected_schema=schema)


# ----------------------------------------------------------------------
# Conversion tests for string and binary types.


class TestConvertStringLikeTypes:

    def test_pandas_unicode(self):
        repeats = 1000
        values = ['foo', None, 'bar', 'mañana', np.nan]
        df = pd.DataFrame({'strings': values * repeats})
        field = pa.field(
            'strings',
            pa.large_string() if _pandas_api.uses_string_dtype() else pa.string()
        )
        schema = pa.schema([field])
        ex_values = ['foo', None, 'bar', 'mañana', None]
        expected = pd.DataFrame({'strings': ex_values * repeats})

        _check_pandas_roundtrip(df, expected=expected, expected_schema=schema)

    def test_bytes_to_binary(self):
        values = ['qux', b'foo', None, bytearray(b'barz'), 'qux', np.nan]
        df = pd.DataFrame({'strings': values})

        table = pa.Table.from_pandas(df)
        assert table[0].type == pa.binary()

        values2 = [b'qux', b'foo', None, b'barz', b'qux', None]
        expected = pd.DataFrame({'strings': values2})
        _check_pandas_roundtrip(df, expected)

    @pytest.mark.large_memory
    def test_bytes_exceed_2gb(self):
        v1 = b'x' * 100000000
        v2 = b'x' * 147483646

        # ARROW-2227, hit exactly 2GB on the nose
        df = pd.DataFrame({
            'strings': [v1] * 20 + [v2] + ['x'] * 20
        })
        arr = pa.array(df['strings'])
        assert isinstance(arr, pa.ChunkedArray)
        assert arr.num_chunks == 2
        arr = None

        table = pa.Table.from_pandas(df)
        assert table[0].num_chunks == 2

    @pytest.mark.large_memory
    @pytest.mark.parametrize('char', ['x', b'x'])
    def test_auto_chunking_pandas_series_of_strings(self, char):
        # ARROW-2367
        v1 = char * 100000000
        v2 = char * 147483646

        df = pd.DataFrame({
            'strings': [[v1]] * 20 + [[v2]] + [[b'x']]
        })
        arr = pa.array(df['strings'], from_pandas=True)
        arr.validate(full=True)
        assert isinstance(arr, pa.ChunkedArray)
        assert arr.num_chunks == 2
        assert len(arr.chunk(0)) == 21
        assert len(arr.chunk(1)) == 1

    def test_fixed_size_bytes(self):
        values = [b'foo', None, bytearray(b'bar'), None, None, b'hey']
        df = pd.DataFrame({'strings': values})
        schema = pa.schema([pa.field('strings', pa.binary(3))])
        table = pa.Table.from_pandas(df, schema=schema)
        assert table.schema[0].type == schema[0].type
        assert table.schema[0].name == schema[0].name
        result = table.to_pandas()
        tm.assert_frame_equal(result, df)

    def test_fixed_size_bytes_does_not_accept_varying_lengths(self):
        values = [b'foo', None, b'ba', None, None, b'hey']
        df = pd.DataFrame({'strings': values})
        schema = pa.schema([pa.field('strings', pa.binary(3))])
        with pytest.raises(pa.ArrowInvalid):
            pa.Table.from_pandas(df, schema=schema)

    def test_variable_size_bytes(self):
        s = pd.Series([b'123', b'', b'a', None])
        _check_series_roundtrip(s, type_=pa.binary())

    def test_binary_from_bytearray(self):
        s = pd.Series([bytearray(b'123'), bytearray(b''), bytearray(b'a'),
                       None])
        # Explicitly set type
        _check_series_roundtrip(s, type_=pa.binary())
        # Infer type from bytearrays
        _check_series_roundtrip(s, expected_pa_type=pa.binary())

    def test_large_binary(self):
        s = pd.Series([b'123', b'', b'a', None])
        _check_series_roundtrip(s, type_=pa.large_binary())
        df = pd.DataFrame({'a': s})
        _check_pandas_roundtrip(
            df, schema=pa.schema([('a', pa.large_binary())]))

    def test_large_string(self):
        s = pd.Series(['123', '', 'a', None])
        _check_series_roundtrip(s, type_=pa.large_string())
        df = pd.DataFrame({'a': s})
        _check_pandas_roundtrip(
            df, schema=pa.schema([('a', pa.large_string())]))

    def test_binary_view(self):
        s = pd.Series([b'123', b'', b'a', None])
        _check_series_roundtrip(s, type_=pa.binary_view())
        df = pd.DataFrame({'a': s})
        _check_pandas_roundtrip(
            df, schema=pa.schema([('a', pa.binary_view())]))

    def test_string_view(self):
        s = pd.Series(['123', '', 'a', None])
        _check_series_roundtrip(s, type_=pa.string_view())
        df = pd.DataFrame({'a': s})
        _check_pandas_roundtrip(
            df, schema=pa.schema([('a', pa.string_view())]))

    def test_table_empty_str(self):
        values = ['', '', '', '', '']
        df = pd.DataFrame({'strings': values})
        field = pa.field('strings', pa.string())
        schema = pa.schema([field])
        table = pa.Table.from_pandas(df, schema=schema)

        result1 = table.to_pandas(strings_to_categorical=False)
        expected1 = pd.DataFrame({'strings': values})
        tm.assert_frame_equal(result1, expected1, check_dtype=True)

        result2 = table.to_pandas(strings_to_categorical=True)
        expected2 = pd.DataFrame({'strings': pd.Categorical(values)})
        tm.assert_frame_equal(result2, expected2, check_dtype=True)

    def test_selective_categoricals(self):
        values = ['', '', '', '', '']
        df = pd.DataFrame({'strings': values})
        field = pa.field('strings', pa.string())
        schema = pa.schema([field])
        table = pa.Table.from_pandas(df, schema=schema)
        expected_str = pd.DataFrame({'strings': values})
        expected_cat = pd.DataFrame({'strings': pd.Categorical(values)})

        result1 = table.to_pandas(categories=['strings'])
        tm.assert_frame_equal(result1, expected_cat, check_dtype=True)
        result2 = table.to_pandas(categories=[])
        tm.assert_frame_equal(result2, expected_str, check_dtype=True)
        result3 = table.to_pandas(categories=('strings',))
        tm.assert_frame_equal(result3, expected_cat, check_dtype=True)
        result4 = table.to_pandas(categories=tuple())
        tm.assert_frame_equal(result4, expected_str, check_dtype=True)

    def test_to_pandas_categorical_zero_length(self):
        # ARROW-3586
        array = pa.array([], type=pa.int32())
        table = pa.Table.from_arrays(arrays=[array], names=['col'])
        # This would segfault under 0.11.0
        table.to_pandas(categories=['col'])

    def test_to_pandas_categories_already_dictionary(self):
        # Showed up in ARROW-6434, ARROW-6435
        array = pa.array(['foo', 'foo', 'foo', 'bar']).dictionary_encode()
        table = pa.Table.from_arrays(arrays=[array], names=['col'])
        result = table.to_pandas(categories=['col'])
        assert table.to_pandas().equals(result)

    @pytest.mark.parametrize(
        "string_type", [pa.string(), pa.large_string(), pa.string_view()]
    )
    def test_table_str_to_categorical_without_na(self, string_type):
        values = ['a', 'a', 'b', 'b', 'c']
        df = pd.DataFrame({'strings': values})
        field = pa.field('strings', string_type)
        schema = pa.schema([field])
        table = pa.Table.from_pandas(df, schema=schema)

        result = table.to_pandas(strings_to_categorical=True)
        expected = pd.DataFrame({'strings': pd.Categorical(values)})
        tm.assert_frame_equal(result, expected, check_dtype=True)

        with pytest.raises(pa.ArrowInvalid):
            table.to_pandas(strings_to_categorical=True,
                            zero_copy_only=True)

        # chunked array
        result = table["strings"].to_pandas(strings_to_categorical=True)
        expected = pd.Series(pd.Categorical(values), name="strings")
        tm.assert_series_equal(result, expected)

        with pytest.raises(pa.ArrowInvalid):
            table["strings"].to_pandas(strings_to_categorical=True,
                                       zero_copy_only=True)

    @pytest.mark.parametrize(
        "string_type", [pa.string(), pa.large_string(), pa.string_view()]
    )
    def test_table_str_to_categorical_with_na(self, string_type):
        values = [None, 'a', 'b', np.nan]
        df = pd.DataFrame({'strings': values})
        field = pa.field('strings', string_type)
        schema = pa.schema([field])
        table = pa.Table.from_pandas(df, schema=schema)

        result = table.to_pandas(strings_to_categorical=True)
        expected = pd.DataFrame({'strings': pd.Categorical(values)})
        tm.assert_frame_equal(result, expected, check_dtype=True)

        with pytest.raises(pa.ArrowInvalid):
            table.to_pandas(strings_to_categorical=True,
                            zero_copy_only=True)

        # chunked array
        result = table["strings"].to_pandas(strings_to_categorical=True)
        expected = pd.Series(pd.Categorical(values), name="strings")
        tm.assert_series_equal(result, expected)

        with pytest.raises(pa.ArrowInvalid):
            table["strings"].to_pandas(strings_to_categorical=True,
                                       zero_copy_only=True)

    # Regression test for ARROW-2101
    def test_array_of_bytes_to_strings(self):
        converted = pa.array(np.array([b'x'], dtype=object), pa.string())
        assert converted.type == pa.string()

    # Make sure that if an ndarray of bytes is passed to the array
    # constructor and the type is string, it will fail if those bytes
    # cannot be converted to utf-8
    def test_array_of_bytes_to_strings_bad_data(self):
        with pytest.raises(
                pa.lib.ArrowInvalid,
                match="was not a utf8 string"):
            pa.array(np.array([b'\x80\x81'], dtype=object), pa.string())

    def test_numpy_string_array_to_fixed_size_binary(self):
        arr = np.array([b'foo', b'bar', b'baz'], dtype='|S3')

        converted = pa.array(arr, type=pa.binary(3))
        expected = pa.array(list(arr), type=pa.binary(3))
        assert converted.equals(expected)

        mask = np.array([False, True, False])
        converted = pa.array(arr, type=pa.binary(3), mask=mask)
        expected = pa.array([b'foo', None, b'baz'], type=pa.binary(3))
        assert converted.equals(expected)

        with pytest.raises(pa.lib.ArrowInvalid,
                           match=r'Got bytestring of length 3 \(expected 4\)'):
            arr = np.array([b'foo', b'bar', b'baz'], dtype='|S3')
            pa.array(arr, type=pa.binary(4))

        with pytest.raises(
                pa.lib.ArrowInvalid,
                match=r'Got bytestring of length 12 \(expected 3\)'):
            arr = np.array([b'foo', b'bar', b'baz'], dtype='|U3')
            pa.array(arr, type=pa.binary(3))


class TestConvertDecimalTypes:
    """
    Conversion test for decimal types.
    """
    decimal32 = [
        decimal.Decimal('-1234.123'),
        decimal.Decimal('1234.439')
    ]
    decimal64 = [
        decimal.Decimal('-129934.123331'),
        decimal.Decimal('129534.123731')
    ]
    decimal128 = [
        decimal.Decimal('394092382910493.12341234678'),
        decimal.Decimal('-314292388910493.12343437128')
    ]

    @pytest.mark.parametrize(('values', 'expected_type'), [
        pytest.param(decimal32, pa.decimal128(7, 3), id='decimal32'),
        pytest.param(decimal64, pa.decimal128(12, 6), id='decimal64'),
        pytest.param(decimal128, pa.decimal128(26, 11), id='decimal128')
    ])
    def test_decimal_from_pandas(self, values, expected_type):
        expected = pd.DataFrame({'decimals': values})
        table = pa.Table.from_pandas(expected, preserve_index=False)
        field = pa.field('decimals', expected_type)

        # schema's metadata is generated by from_pandas conversion
        expected_schema = pa.schema([field], metadata=table.schema.metadata)
        assert table.schema.equals(expected_schema)

    @pytest.mark.parametrize('values', [
        pytest.param(decimal32, id='decimal32'),
        pytest.param(decimal64, id='decimal64'),
        pytest.param(decimal128, id='decimal128')
    ])
    def test_decimal_to_pandas(self, values):
        expected = pd.DataFrame({'decimals': values})
        converted = pa.Table.from_pandas(expected)
        df = converted.to_pandas()
        tm.assert_frame_equal(df, expected)

    def test_decimal_fails_with_truncation(self):
        data1 = [decimal.Decimal('1.234')]
        type1 = pa.decimal128(10, 2)
        with pytest.raises(pa.ArrowInvalid):
            pa.array(data1, type=type1)

        data2 = [decimal.Decimal('1.2345')]
        type2 = pa.decimal128(10, 3)
        with pytest.raises(pa.ArrowInvalid):
            pa.array(data2, type=type2)

    def test_decimal_with_different_precisions(self):
        data = [
            decimal.Decimal('0.01'),
            decimal.Decimal('0.001'),
        ]
        series = pd.Series(data)
        array = pa.array(series)
        assert array.to_pylist() == data
        assert array.type == pa.decimal128(3, 3)

        array = pa.array(data, type=pa.decimal128(12, 5))
        expected = [decimal.Decimal('0.01000'), decimal.Decimal('0.00100')]
        assert array.to_pylist() == expected

    def test_decimal_with_None_explicit_type(self):
        series = pd.Series([decimal.Decimal('3.14'), None])
        _check_series_roundtrip(series, type_=pa.decimal128(12, 5))

        # Test that having all None values still produces decimal array
        series = pd.Series([None] * 2)
        _check_series_roundtrip(series, type_=pa.decimal128(12, 5))

    def test_decimal_with_None_infer_type(self):
        series = pd.Series([decimal.Decimal('3.14'), None])
        _check_series_roundtrip(series, expected_pa_type=pa.decimal128(3, 2))

    def test_strided_objects(self, tmpdir):
        # see ARROW-3053
        data = {
            'a': {0: 'a'},
            'b': {0: decimal.Decimal('0.0')}
        }

        # This yields strided objects
        df = pd.DataFrame.from_dict(data)
        _check_pandas_roundtrip(df)

    @pytest.mark.parametrize("typ", [
        pa.decimal32,
        pa.decimal64,
        pa.decimal128,
        pa.decimal256,
    ])
    def test_decimal_array_to_pandas(self, typ):
        data = [decimal.Decimal('3.14'), None]
        arr = pa.array(data, type=typ(3, 2))
        result = arr.to_pandas()
        expected = pd.Series(data)
        tm.assert_series_equal(result, expected)


class TestConvertListTypes:
    """
    Conversion tests for list<> types.
    """

    def test_column_of_arrays(self):
        df, schema = dataframe_with_arrays()
        _check_pandas_roundtrip(df, schema=schema, expected_schema=schema)
        table = pa.Table.from_pandas(df, schema=schema, preserve_index=False)

        # schema's metadata is generated by from_pandas conversion
        expected_schema = schema.with_metadata(table.schema.metadata)
        assert table.schema.equals(expected_schema)

        for column in df.columns:
            field = schema.field(column)
            _check_array_roundtrip(df[column], type=field.type)

    def test_column_of_arrays_to_py(self):
        # Test regression in ARROW-1199 not caught in above test
        dtype = 'i1'
        arr = np.array([
            np.arange(10, dtype=dtype),
            np.arange(5, dtype=dtype),
            None,
            np.arange(1, dtype=dtype)
        ], dtype=object)
        type_ = pa.list_(pa.int8())
        parr = pa.array(arr, type=type_)

        assert parr[0].as_py() == list(range(10))
        assert parr[1].as_py() == list(range(5))
        assert parr[2].as_py() is None
        assert parr[3].as_py() == [0]

    def test_column_of_boolean_list(self):
        # ARROW-4370: Table to pandas conversion fails for list of bool
        array = pa.array([[True, False], [True]], type=pa.list_(pa.bool_()))
        table = pa.Table.from_arrays([array], names=['col1'])
        df = table.to_pandas()

        expected_df = pd.DataFrame({'col1': [[True, False], [True]]})
        tm.assert_frame_equal(df, expected_df)

        s = table[0].to_pandas()
        tm.assert_series_equal(pd.Series(s), df['col1'], check_names=False)

    def test_column_of_decimal_list(self):
        array = pa.array([[decimal.Decimal('1'), decimal.Decimal('2')],
                          [decimal.Decimal('3.3')]],
                         type=pa.list_(pa.decimal128(2, 1)))
        table = pa.Table.from_arrays([array], names=['col1'])
        df = table.to_pandas()

        expected_df = pd.DataFrame(
            {'col1': [[decimal.Decimal('1'), decimal.Decimal('2')],
                      [decimal.Decimal('3.3')]]})
        tm.assert_frame_equal(df, expected_df)

    def test_nested_types_from_ndarray_null_entries(self):
        # Root cause of ARROW-6435
        s = pd.Series(np.array([np.nan, np.nan], dtype=object))

        for ty in [pa.list_(pa.int64()),
                   pa.large_list(pa.int64()),
                   pa.struct([pa.field('f0', 'int32')])]:
            result = pa.array(s, type=ty)
            expected = pa.array([None, None], type=ty)
            assert result.equals(expected)

            with pytest.raises(TypeError):
                pa.array(s.values, type=ty)

    def test_column_of_lists(self):
        df, schema = dataframe_with_lists()
        _check_pandas_roundtrip(df, schema=schema, expected_schema=schema)
        table = pa.Table.from_pandas(df, schema=schema, preserve_index=False)

        # schema's metadata is generated by from_pandas conversion
        expected_schema = schema.with_metadata(table.schema.metadata)
        assert table.schema.equals(expected_schema)

        for column in df.columns:
            field = schema.field(column)
            _check_array_roundtrip(df[column], type=field.type)

    def test_column_of_lists_first_empty(self):
        # ARROW-2124
        num_lists = [[], [2, 3, 4], [3, 6, 7, 8], [], [2]]
        series = pd.Series([np.array(s, dtype=float) for s in num_lists])
        arr = pa.array(series)
        result = pd.Series(arr.to_pandas())
        tm.assert_series_equal(result, series)

    def test_column_of_lists_chunked(self):
        # ARROW-1357
        df = pd.DataFrame({
            'lists': np.array([
                [1, 2],
                None,
                [2, 3],
                [4, 5],
                [6, 7],
                [8, 9]
            ], dtype=object)
        })

        schema = pa.schema([
            pa.field('lists', pa.list_(pa.int64()))
        ])

        t1 = pa.Table.from_pandas(df[:2], schema=schema)
        t2 = pa.Table.from_pandas(df[2:], schema=schema)

        table = pa.concat_tables([t1, t2])
        result = table.to_pandas()

        tm.assert_frame_equal(result, df)

    def test_empty_column_of_lists_chunked(self):
        df = pd.DataFrame({
            'lists': np.array([], dtype=object)
        })

        schema = pa.schema([
            pa.field('lists', pa.list_(pa.int64()))
        ])

        table = pa.Table.from_pandas(df, schema=schema)
        result = table.to_pandas()

        tm.assert_frame_equal(result, df)

    def test_column_of_lists_chunked2(self):
        data1 = [[0, 1], [2, 3], [4, 5], [6, 7], [10, 11],
                 [12, 13], [14, 15], [16, 17]]
        data2 = [[8, 9], [18, 19]]

        a1 = pa.array(data1)
        a2 = pa.array(data2)

        t1 = pa.Table.from_arrays([a1], names=['a'])
        t2 = pa.Table.from_arrays([a2], names=['a'])

        concatenated = pa.concat_tables([t1, t2])

        result = concatenated.to_pandas()
        expected = pd.DataFrame({'a': data1 + data2})

        tm.assert_frame_equal(result, expected)

    def test_column_of_lists_strided(self):
        df, schema = dataframe_with_lists()
        df = pd.concat([df] * 6, ignore_index=True)

        arr = df['int64'].values[::3]
        assert arr.strides[0] != 8

        _check_array_roundtrip(arr)

    def test_nested_lists_all_none(self):
        data = np.array([[None, None], None], dtype=object)

        arr = pa.array(data)
        expected = pa.array(list(data))
        assert arr.equals(expected)
        assert arr.type == pa.list_(pa.null())

        data2 = np.array([None, None, [None, None],
                          np.array([None, None], dtype=object)],
                         dtype=object)
        arr = pa.array(data2)
        expected = pa.array([None, None, [None, None], [None, None]])
        assert arr.equals(expected)

    def test_nested_lists_all_empty(self):
        # ARROW-2128
        data = pd.Series([[], [], []])
        arr = pa.array(data)
        expected = pa.array(list(data))
        assert arr.equals(expected)
        assert arr.type == pa.list_(pa.null())

    def test_nested_list_first_empty(self):
        # ARROW-2711
        data = pd.Series([[], ["a"]])
        arr = pa.array(data)
        expected = pa.array(list(data))
        assert arr.equals(expected)
        assert arr.type == pa.list_(pa.string())

    def test_nested_smaller_ints(self):
        # ARROW-1345, ARROW-2008, there were some type inference bugs happening
        # before
        data = pd.Series([np.array([1, 2, 3], dtype='i1'), None])
        result = pa.array(data)
        result2 = pa.array(data.values)
        expected = pa.array([[1, 2, 3], None], type=pa.list_(pa.int8()))
        assert result.equals(expected)
        assert result2.equals(expected)

        data3 = pd.Series([np.array([1, 2, 3], dtype='f4'), None])
        result3 = pa.array(data3)
        expected3 = pa.array([[1, 2, 3], None], type=pa.list_(pa.float32()))
        assert result3.equals(expected3)

    def test_infer_lists(self):
        data = OrderedDict([
            ('nan_ints', [[np.nan, 1], [2, 3]]),
            ('ints', [[0, 1], [2, 3]]),
            ('strs', [[None, 'b'], ['c', 'd']]),
            ('nested_strs', [[[None, 'b'], ['c', 'd']], None])
        ])
        df = pd.DataFrame(data)

        expected_schema = pa.schema([
            pa.field('nan_ints', pa.list_(pa.int64())),
            pa.field('ints', pa.list_(pa.int64())),
            pa.field('strs', pa.list_(pa.string())),
            pa.field('nested_strs', pa.list_(pa.list_(pa.string())))
        ])

        _check_pandas_roundtrip(df, expected_schema=expected_schema)

    def test_fixed_size_list(self):
        # ARROW-7365
        fixed_ty = pa.list_(pa.int64(), list_size=4)
        variable_ty = pa.list_(pa.int64())

        data = [[0, 1, 2, 3], None, [4, 5, 6, 7], [8, 9, 10, 11]]
        fixed_arr = pa.array(data, type=fixed_ty)
        variable_arr = pa.array(data, type=variable_ty)

        result = fixed_arr.to_pandas()
        expected = variable_arr.to_pandas()

        for left, right in zip(result, expected):
            if left is None:
                assert right is None
            npt.assert_array_equal(left, right)

    def test_infer_numpy_array(self):
        data = OrderedDict([
            ('ints', [
                np.array([0, 1], dtype=np.int64),
                np.array([2, 3], dtype=np.int64)
            ])
        ])
        df = pd.DataFrame(data)
        expected_schema = pa.schema([
            pa.field('ints', pa.list_(pa.int64()))
        ])

        _check_pandas_roundtrip(df, expected_schema=expected_schema)

    def test_to_list_of_structs_pandas(self):
        ints = pa.array([1, 2, 3], pa.int32())
        strings = pa.array([['a', 'b'], ['c', 'd'], ['e', 'f']],
                           pa.list_(pa.string()))
        structs = pa.StructArray.from_arrays([ints, strings], ['f1', 'f2'])
        data = pa.ListArray.from_arrays([0, 1, 3], structs)

        expected = pd.Series([
            [{'f1': 1, 'f2': ['a', 'b']}],
            [{'f1': 2, 'f2': ['c', 'd']},
             {'f1': 3, 'f2': ['e', 'f']}]
        ])

        series = pd.Series(data.to_pandas())

        # pandas.testing generates a
        # DeprecationWarning: elementwise comparison failed
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "elementwise comparison failed",
                                    DeprecationWarning)
            tm.assert_series_equal(series, expected)

    def test_to_list_of_maps_pandas(self):
        if ((Version(np.__version__) >= Version("1.25.0.dev0")) and
                (Version(pd.__version__) < Version("2.0.0"))):
            # TODO: regression in pandas with numpy 1.25dev
            # https://github.com/pandas-dev/pandas/issues/50360
            pytest.skip("Regression in pandas with numpy 1.25")
        data = [
            [[('foo', ['a', 'b']), ('bar', ['c', 'd'])]],
            [[('baz', []), ('qux', None), ('quux', [None, 'e'])], [('quz', ['f', 'g'])]]
        ]
        arr = pa.array(data, pa.list_(pa.map_(pa.utf8(), pa.list_(pa.utf8()))))
        series = arr.to_pandas()
        expected = pd.Series(data)

        # pandas.testing generates a
        # DeprecationWarning: elementwise comparison failed
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "elementwise comparison failed",
                                    DeprecationWarning)
            tm.assert_series_equal(series, expected)

    def test_to_list_of_maps_pandas_sliced(self):
        """
        A slightly more rigorous test for chunk/slice combinations
        """

        if ((Version(np.__version__) >= Version("1.25.0.dev0")) and
                (Version(pd.__version__) < Version("2.0.0"))):
            # TODO: regression in pandas with numpy 1.25dev
            # https://github.com/pandas-dev/pandas/issues/50360
            pytest.skip("Regression in pandas with numpy 1.25")

        keys = pa.array(['ignore', 'foo', 'bar', 'baz',
                         'qux', 'quux', 'ignore']).slice(1, 5)
        items = pa.array(
            [['ignore'], ['ignore'], ['a', 'b'], ['c', 'd'], [], None, [None, 'e']],
            pa.list_(pa.string()),
        ).slice(2, 5)
        map = pa.MapArray.from_arrays([0, 2, 4], keys, items)
        arr = pa.ListArray.from_arrays([0, 1, 2], map)

        series = arr.to_pandas()
        expected = pd.Series([
            [[('foo', ['a', 'b']), ('bar', ['c', 'd'])]],
            [[('baz', []), ('qux', None)]],
        ])

        series_sliced = arr.slice(1, 2).to_pandas()
        expected_sliced = pd.Series([
            [[('baz', []), ('qux', None)]],
        ])

        # pandas.testing generates a
        # DeprecationWarning: elementwise comparison failed
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", "elementwise comparison failed",
                                    DeprecationWarning)
            tm.assert_series_equal(series, expected)
            tm.assert_series_equal(series_sliced, expected_sliced)

    @pytest.mark.parametrize('t,data,expected', [
        (
            pa.int64,
            [[1, 2], [3], None],
            [None, [3], None]
        ),
        (
            pa.string,
            [['aaa', 'bb'], ['c'], None],
            [None, ['c'], None]
        ),
        (
            pa.null,
            [[None, None], [None], None],
            [None, [None], None]
        )
    ])
    def test_array_from_pandas_typed_array_with_mask(self, t, data, expected):
        m = np.array([True, False, True])

        s = pd.Series(data)
        result = pa.Array.from_pandas(s, mask=m, type=pa.list_(t()))

        assert pa.Array.from_pandas(expected,
                                    type=pa.list_(t())).equals(result)

    def test_empty_list_roundtrip(self):
        empty_list_array = np.empty((3,), dtype=object)
        empty_list_array.fill([])

        df = pd.DataFrame({'a': np.array(['1', '2', '3']),
                           'b': empty_list_array})
        tbl = pa.Table.from_pandas(df)

        result = tbl.to_pandas()

        tm.assert_frame_equal(result, df)

    def test_array_from_nested_arrays(self):
        df, schema = dataframe_with_arrays()
        for field in schema:
            arr = df[field.name].values
            expected = pa.array(list(arr), type=field.type)
            result = pa.array(arr)
            assert result.type == field.type  # == list<scalar>
            assert result.equals(expected)

    def test_nested_large_list(self):
        s = (pa.array([[[1, 2, 3], [4]], None],
                      type=pa.large_list(pa.large_list(pa.int64())))
             .to_pandas())

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",
                                    "Creating an ndarray from ragged nested",
                                    _np_VisibleDeprecationWarning)
            warnings.filterwarnings("ignore", "elementwise comparison failed",
                                    DeprecationWarning)
            tm.assert_series_equal(
                s, pd.Series([[[1, 2, 3], [4]], None], dtype=object),
                check_names=False)

    def test_large_binary_list(self):
        for list_type_factory in (pa.list_, pa.large_list):
            s = (pa.array([["aa", "bb"], None, ["cc"], []],
                          type=list_type_factory(pa.large_binary()))
                 .to_pandas())
            tm.assert_series_equal(
                s, pd.Series([[b"aa", b"bb"], None, [b"cc"], []]),
                check_names=False)
            s = (pa.array([["aa", "bb"], None, ["cc"], []],
                          type=list_type_factory(pa.large_string()))
                 .to_pandas())
            tm.assert_series_equal(
                s, pd.Series([["aa", "bb"], None, ["cc"], []]),
                check_names=False)

    def test_list_of_dictionary(self):
        child = pa.array(["foo", "bar", None, "foo"]).dictionary_encode()
        arr = pa.ListArray.from_arrays([0, 1, 3, 3, 4], child)

        # Expected a Series of lists
        expected = pd.Series(arr.to_pylist())
        tm.assert_series_equal(arr.to_pandas(), expected)

        # Same but with nulls
        arr = arr.take([0, 1, None, 3])
        expected[2] = None
        tm.assert_series_equal(arr.to_pandas(), expected)

    @pytest.mark.large_memory
    def test_auto_chunking_on_list_overflow(self):
        # ARROW-9976
        n = 2**21
        df = pd.DataFrame.from_dict({
            "a": list(np.zeros((n, 2**10), dtype='uint8')),
            "b": range(n)
        })
        table = pa.Table.from_pandas(df)
        table.validate(full=True)

        column_a = table[0]
        assert column_a.num_chunks == 2
        assert len(column_a.chunk(0)) == 2**21 - 1
        assert len(column_a.chunk(1)) == 1

    def test_map_array_roundtrip(self):
        data = [[(b'a', 1), (b'b', 2)],
                [(b'c', 3)],
                [(b'd', 4), (b'e', 5), (b'f', 6)],
                [(b'g', 7)]]

        df = pd.DataFrame({"map": data})
        schema = pa.schema([("map", pa.map_(pa.binary(), pa.int32()))])

        _check_pandas_roundtrip(df, schema=schema)

    def test_map_array_chunked(self):
        data1 = [[(b'a', 1), (b'b', 2)],
                 [(b'c', 3)],
                 [(b'd', 4), (b'e', 5), (b'f', 6)],
                 [(b'g', 7)]]
        data2 = [[(k, v * 2) for k, v in row] for row in data1]

        arr1 = pa.array(data1, type=pa.map_(pa.binary(), pa.int32()))
        arr2 = pa.array(data2, type=pa.map_(pa.binary(), pa.int32()))
        arr = pa.chunked_array([arr1, arr2])

        expected = pd.Series(data1 + data2)
        actual = arr.to_pandas()
        tm.assert_series_equal(actual, expected, check_names=False)

    def test_map_array_with_nulls(self):
        data = [[(b'a', 1), (b'b', 2)],
                None,
                [(b'd', 4), (b'e', 5), (b'f', None)],
                [(b'g', 7)]]

        # None value in item array causes upcast to float
        expected = [[(k, float(v) if v is not None else None) for k, v in row]
                    if row is not None else None for row in data]
        expected = pd.Series(expected)

        arr = pa.array(data, type=pa.map_(pa.binary(), pa.int32()))
        actual = arr.to_pandas()
        tm.assert_series_equal(actual, expected, check_names=False)

    def test_map_array_dictionary_encoded(self):
        offsets = pa.array([0, 3, 5])
        items = pa.array(['a', 'b', 'c', 'a', 'd']).dictionary_encode()
        keys = pa.array(list(range(len(items))))
        arr = pa.MapArray.from_arrays(offsets, keys, items)

        # Dictionary encoded values converted to dense
        expected = pd.Series(
            [[(0, 'a'), (1, 'b'), (2, 'c')], [(3, 'a'), (4, 'd')]])

        actual = arr.to_pandas()
        tm.assert_series_equal(actual, expected, check_names=False)

    def test_list_no_duplicate_base(self):
        # ARROW-18400
        arr = pa.array([[1, 2], [3, 4, 5], None, [6, None], [7, 8]])
        chunked_arr = pa.chunked_array([arr.slice(0, 3), arr.slice(3, 1)])

        np_arr = chunked_arr.to_numpy()

        expected = np.array([[1., 2.], [3., 4., 5.], None,
                            [6., np.nan]], dtype="object")
        for left, right in zip(np_arr, expected):
            if right is None:
                assert left == right
            else:
                npt.assert_array_equal(left, right)

        expected_base = np.array([[1., 2., 3., 4., 5., 6., np.nan]])
        npt.assert_array_equal(np_arr[0].base, expected_base)

        np_arr_sliced = chunked_arr.slice(1, 3).to_numpy()

        expected = np.array([[3, 4, 5], None, [6, np.nan]], dtype="object")
        for left, right in zip(np_arr_sliced, expected):
            if right is None:
                assert left == right
            else:
                npt.assert_array_equal(left, right)

        expected_base = np.array([[3., 4., 5., 6., np.nan]])
        npt.assert_array_equal(np_arr_sliced[0].base, expected_base)

    def test_list_values_behind_null(self):
        arr = pa.ListArray.from_arrays(
            offsets=pa.array([0, 2, 4, 6]),
            values=pa.array([1, 2, 99, 99, 3, None]),
            mask=pa.array([False, True, False])
        )
        np_arr = arr.to_numpy(zero_copy_only=False)

        expected = np.array([[1., 2.], None, [3., np.nan]], dtype="object")
        for left, right in zip(np_arr, expected):
            if right is None:
                assert left == right
            else:
                npt.assert_array_equal(left, right)

    @pytest.mark.parametrize("klass", [pa.ListViewArray, pa.LargeListViewArray])
    def test_list_view_to_pandas_with_in_order_offsets(self, klass):
        arr = klass.from_arrays(
            offsets=pa.array([0, 2, 4]),
            sizes=pa.array([2, 2, 2]),
            values=pa.array([1, 2, 3, 4, 5, 6]),
        )

        actual = arr.to_pandas()
        expected = pd.Series([[1, 2], [3, 4], [5, 6]])

        tm.assert_series_equal(actual, expected)

    @pytest.mark.parametrize("klass", [pa.ListViewArray, pa.LargeListViewArray])
    def test_list_view_to_pandas_with_out_of_order_offsets(self, klass):
        arr = klass.from_arrays(
            offsets=pa.array([2, 4, 0]),
            sizes=pa.array([2, 2, 2]),
            values=pa.array([1, 2, 3, 4, 5, 6]),
        )

        actual = arr.to_pandas()
        expected = pd.Series([[3, 4], [5, 6], [1, 2]])

        tm.assert_series_equal(actual, expected)

    @pytest.mark.parametrize("klass", [pa.ListViewArray, pa.LargeListViewArray])
    def test_list_view_to_pandas_with_overlapping_offsets(self, klass):
        arr = klass.from_arrays(
            offsets=pa.array([0, 1, 2]),
            sizes=pa.array([4, 4, 4]),
            values=pa.array([1, 2, 3, 4, 5, 6]),
        )

        actual = arr.to_pandas()
        expected = pd.Series([[1, 2, 3, 4], [2, 3, 4, 5], [3, 4, 5, 6]])

        tm.assert_series_equal(actual, expected)

    @pytest.mark.parametrize("klass", [pa.ListViewArray, pa.LargeListViewArray])
    def test_list_view_to_pandas_with_null_values(self, klass):
        arr = klass.from_arrays(
            offsets=pa.array([0, 2, 2]),
            sizes=pa.array([2, 0, 0]),
            values=pa.array([1, None]),
            mask=pa.array([False, False, True])
        )

        actual = arr.to_pandas()
        expected = pd.Series([[1, np.nan], [], None])

        tm.assert_series_equal(actual, expected)

    @pytest.mark.parametrize("klass", [pa.ListViewArray, pa.LargeListViewArray])
    def test_list_view_to_pandas_multiple_chunks(self, klass):
        gc.collect()
        bytes_start = pa.total_allocated_bytes()
        arr1 = klass.from_arrays(
            offsets=pa.array([2, 1, 0]),
            sizes=pa.array([2, 2, 2]),
            values=pa.array([1, 2, 3, 4])
        )
        arr2 = klass.from_arrays(
            offsets=pa.array([0, 1, 1]),
            sizes=pa.array([3, 3, 0]),
            values=pa.array([5, 6, 7, None]),
            mask=pa.array([False, False, True])
        )
        arr = pa.chunked_array([arr1, arr2])

        actual = arr.to_pandas()
        expected = pd.Series([[3, 4], [2, 3], [1, 2], [5, 6, 7], [6, 7, np.nan], None])

        tm.assert_series_equal(actual, expected)

        del actual
        del arr
        del arr1
        del arr2
        bytes_end = pa.total_allocated_bytes()
        assert bytes_end == bytes_start


class TestConvertStructTypes:
    """
    Conversion tests for struct types.
    """

    def test_pandas_roundtrip(self):
        df = pd.DataFrame({'dicts': [{'a': 1, 'b': 2}, {'a': 3, 'b': 4}]})

        expected_schema = pa.schema([
            ('dicts', pa.struct([('a', pa.int64()), ('b', pa.int64())])),
        ])

        _check_pandas_roundtrip(df, expected_schema=expected_schema)

        # specifying schema explicitly in from_pandas
        _check_pandas_roundtrip(
            df, schema=expected_schema, expected_schema=expected_schema)

    def test_to_pandas(self):
        ints = pa.array([None, 2, 3], type=pa.int64())
        strs = pa.array(['a', None, 'c'], type=pa.string())
        bools = pa.array([True, False, None], type=pa.bool_())
        arr = pa.StructArray.from_arrays(
            [ints, strs, bools],
            ['ints', 'strs', 'bools'])

        expected = pd.Series([
            {'ints': None, 'strs': 'a', 'bools': True},
            {'ints': 2, 'strs': None, 'bools': False},
            {'ints': 3, 'strs': 'c', 'bools': None},
        ])

        series = pd.Series(arr.to_pandas())
        tm.assert_series_equal(series, expected)

    def test_to_pandas_multiple_chunks(self):
        # ARROW-11855
        gc.collect()
        bytes_start = pa.total_allocated_bytes()
        ints1 = pa.array([1], type=pa.int64())
        ints2 = pa.array([2], type=pa.int64())
        arr1 = pa.StructArray.from_arrays([ints1], ['ints'])
        arr2 = pa.StructArray.from_arrays([ints2], ['ints'])
        arr = pa.chunked_array([arr1, arr2])

        expected = pd.Series([
            {'ints': 1},
            {'ints': 2}
        ])

        series = pd.Series(arr.to_pandas())
        tm.assert_series_equal(series, expected)

        del series
        del arr
        del arr1
        del arr2
        del ints1
        del ints2
        bytes_end = pa.total_allocated_bytes()
        assert bytes_end == bytes_start

    def test_from_numpy(self):
        dt = np.dtype([('x', np.int32),
                       (('y_title', 'y'), np.bool_)])
        ty = pa.struct([pa.field('x', pa.int32()),
                        pa.field('y', pa.bool_())])

        data = np.array([], dtype=dt)
        arr = pa.array(data, type=ty)
        assert arr.to_pylist() == []

        data = np.array([(42, True), (43, False)], dtype=dt)
        arr = pa.array(data, type=ty)
        assert arr.to_pylist() == [{'x': 42, 'y': True},
                                   {'x': 43, 'y': False}]

        # With mask
        arr = pa.array(data, mask=np.bool_([False, True]), type=ty)
        assert arr.to_pylist() == [{'x': 42, 'y': True}, None]

        # Trivial struct type
        dt = np.dtype([])
        ty = pa.struct([])

        data = np.array([], dtype=dt)
        arr = pa.array(data, type=ty)
        assert arr.to_pylist() == []

        data = np.array([(), ()], dtype=dt)
        arr = pa.array(data, type=ty)
        assert arr.to_pylist() == [{}, {}]

    def test_from_numpy_nested(self):
        # Note: an object field inside a struct
        dt = np.dtype([('x', np.dtype([('xx', np.int8),
                                       ('yy', np.bool_)])),
                       ('y', np.int16),
                       ('z', np.object_)])
        # Note: itemsize is not necessarily a multiple of sizeof(object)
        # object_ is 8 bytes on 64-bit systems, 4 bytes on 32-bit systems
        assert dt.itemsize == (12 if sys.maxsize > 2**32 else 8)
        ty = pa.struct([pa.field('x', pa.struct([pa.field('xx', pa.int8()),
                                                 pa.field('yy', pa.bool_())])),
                        pa.field('y', pa.int16()),
                        pa.field('z', pa.string())])

        data = np.array([], dtype=dt)
        arr = pa.array(data, type=ty)
        assert arr.to_pylist() == []

        data = np.array([
            ((1, True), 2, 'foo'),
            ((3, False), 4, 'bar')], dtype=dt)
        arr = pa.array(data, type=ty)
        assert arr.to_pylist() == [
            {'x': {'xx': 1, 'yy': True}, 'y': 2, 'z': 'foo'},
            {'x': {'xx': 3, 'yy': False}, 'y': 4, 'z': 'bar'}]

    @pytest.mark.slow
    @pytest.mark.large_memory
    def test_from_numpy_large(self):
        # Exercise rechunking + nulls
        target_size = 3 * 1024**3  # 4GB
        dt = np.dtype([('x', np.float64), ('y', 'object')])
        bs = 65536 - dt.itemsize
        block = b'.' * bs
        n = target_size // (bs + dt.itemsize)
        data = np.zeros(n, dtype=dt)
        data['x'] = np.random.random_sample(n)
        data['y'] = block
        # Add implicit nulls
        data['x'][data['x'] < 0.2] = np.nan

        ty = pa.struct([pa.field('x', pa.float64()),
                        pa.field('y', pa.binary())])
        arr = pa.array(data, type=ty, from_pandas=True)
        arr.validate(full=True)
        assert arr.num_chunks == 2

        def iter_chunked_array(arr):
            for chunk in arr.iterchunks():
                yield from chunk

        def check(arr, data, mask=None):
            assert len(arr) == len(data)
            xs = data['x']
            ys = data['y']
            for i, obj in enumerate(iter_chunked_array(arr)):
                try:
                    d = obj.as_py()
                    if mask is not None and mask[i]:
                        assert d is None
                    else:
                        x = xs[i]
                        if np.isnan(x):
                            assert d['x'] is None
                        else:
                            assert d['x'] == x
                        assert d['y'] == ys[i]
                except Exception:
                    print("Failed at index", i)
                    raise

        check(arr, data)
        del arr

        # Now with explicit mask
        mask = np.random.random_sample(n) < 0.2
        arr = pa.array(data, type=ty, mask=mask, from_pandas=True)
        arr.validate(full=True)
        assert arr.num_chunks == 2

        check(arr, data, mask)
        del arr

    def test_from_numpy_bad_input(self):
        ty = pa.struct([pa.field('x', pa.int32()),
                        pa.field('y', pa.bool_())])
        dt = np.dtype([('x', np.int32),
                       ('z', np.bool_)])

        data = np.array([], dtype=dt)
        with pytest.raises(ValueError,
                           match="Missing field 'y'"):
            pa.array(data, type=ty)
        data = np.int32([])
        with pytest.raises(TypeError,
                           match="Expected struct array"):
            pa.array(data, type=ty)

    def test_from_tuples(self):
        df = pd.DataFrame({'tuples': [(1, 2), (3, 4)]})
        expected_df = pd.DataFrame(
            {'tuples': [{'a': 1, 'b': 2}, {'a': 3, 'b': 4}]})

        # conversion from tuples works when specifying expected struct type
        struct_type = pa.struct([('a', pa.int64()), ('b', pa.int64())])

        arr = np.asarray(df['tuples'])
        _check_array_roundtrip(
            arr, expected=expected_df['tuples'], type=struct_type)

        expected_schema = pa.schema([('tuples', struct_type)])
        _check_pandas_roundtrip(
            df, expected=expected_df, schema=expected_schema,
            expected_schema=expected_schema)

    def test_struct_of_dictionary(self):
        names = ['ints', 'strs']
        children = [pa.array([456, 789, 456]).dictionary_encode(),
                    pa.array(["foo", "foo", None]).dictionary_encode()]
        arr = pa.StructArray.from_arrays(children, names=names)

        # Expected a Series of {field name: field value} dicts
        rows_as_tuples = zip(*(child.to_pylist() for child in children))
        rows_as_dicts = [dict(zip(names, row)) for row in rows_as_tuples]

        expected = pd.Series(rows_as_dicts)
        tm.assert_series_equal(arr.to_pandas(), expected)

        # Same but with nulls
        arr = arr.take([0, None, 2])
        expected[1] = None
        tm.assert_series_equal(arr.to_pandas(), expected)


class TestZeroCopyConversion:
    """
    Tests that zero-copy conversion works with some types.
    """

    def test_zero_copy_success(self):
        result = pa.array([0, 1, 2]).to_pandas(zero_copy_only=True)
        npt.assert_array_equal(result, [0, 1, 2])

    def test_zero_copy_dictionaries(self):
        arr = pa.DictionaryArray.from_arrays(
            np.array([0, 0]),
            np.array([5], dtype="int64"),
        )

        result = arr.to_pandas(zero_copy_only=True)
        values = pd.Categorical([5, 5])

        tm.assert_series_equal(pd.Series(result), pd.Series(values),
                               check_names=False)

    def test_zero_copy_timestamp(self):
        arr = np.array(['2007-07-13'], dtype='datetime64[ns]')
        result = pa.array(arr).to_pandas(zero_copy_only=True)
        npt.assert_array_equal(result, arr)

    def test_zero_copy_duration(self):
        arr = np.array([1], dtype='timedelta64[ns]')
        result = pa.array(arr).to_pandas(zero_copy_only=True)
        npt.assert_array_equal(result, arr)

    def check_zero_copy_failure(self, arr):
        with pytest.raises(pa.ArrowInvalid):
            arr.to_pandas(zero_copy_only=True)

    def test_zero_copy_failure_on_object_types(self):
        self.check_zero_copy_failure(pa.array(['A', 'B', 'C']))

    def test_zero_copy_failure_with_int_when_nulls(self):
        self.check_zero_copy_failure(pa.array([0, 1, None]))

    def test_zero_copy_failure_with_float_when_nulls(self):
        self.check_zero_copy_failure(pa.array([0.0, 1.0, None]))

    def test_zero_copy_failure_on_bool_types(self):
        self.check_zero_copy_failure(pa.array([True, False]))

    def test_zero_copy_failure_on_list_types(self):
        arr = pa.array([[1, 2], [8, 9]], type=pa.list_(pa.int64()))
        self.check_zero_copy_failure(arr)

    def test_zero_copy_failure_on_timestamp_with_nulls(self):
        arr = np.array([1, None], dtype='datetime64[ns]')
        self.check_zero_copy_failure(pa.array(arr))

    def test_zero_copy_failure_on_duration_with_nulls(self):
        arr = np.array([1, None], dtype='timedelta64[ns]')
        self.check_zero_copy_failure(pa.array(arr))


def _non_threaded_conversion():
    df = _alltypes_example()
    _check_pandas_roundtrip(df, use_threads=False)
    _check_pandas_roundtrip(df, use_threads=False, as_batch=True)


def _threaded_conversion():
    df = _alltypes_example()
    _check_pandas_roundtrip(df, use_threads=True)
    _check_pandas_roundtrip(df, use_threads=True, as_batch=True)


class TestConvertMisc:
    """
    Miscellaneous conversion tests.
    """

    type_pairs = [
        ("int8", pa.int8()),
        ("int16", pa.int16()),
        ("int32", pa.int32()),
        ("int64", pa.int64()),
        ("uint8", pa.uint8()),
        ("uint16", pa.uint16()),
        ("uint32", pa.uint32()),
        ("uint64", pa.uint64()),
        ("float16", pa.float16()),
        ("float32", pa.float32()),
        ("float64", pa.float64()),
        # XXX unsupported
        # (np.dtype([('a', 'i2')]), pa.struct([pa.field('a', pa.int16())])),
        ("object", pa.string()),
        ("object", pa.binary()),
        ("object", pa.binary(10)),
        ("object", pa.list_(pa.int64())),
    ]

    def test_all_none_objects(self):
        df = pd.DataFrame({'a': [None, None, None]})
        _check_pandas_roundtrip(df)

    def test_all_none_category(self):
        df = pd.DataFrame({'a': [None, None, None]})
        df['a'] = df['a'].astype('category')
        _check_pandas_roundtrip(df)

    def test_empty_arrays(self):
        for dtype_str, pa_type in self.type_pairs:
            arr = np.array([], dtype=np.dtype(dtype_str))
            _check_array_roundtrip(arr, type=pa_type)

    def test_non_threaded_conversion(self):
        _non_threaded_conversion()

    @pytest.mark.processes
    @pytest.mark.threading
    def test_threaded_conversion_multiprocess(self):
        # Parallel conversion should work from child processes too (ARROW-2963)
        pool = mp.Pool(2)
        try:
            pool.apply(_threaded_conversion)
        finally:
            pool.close()
            pool.join()

    def test_category(self):
        repeats = 5
        v1 = ['foo', None, 'bar', 'qux', np.nan]
        v2 = [4, 5, 6, 7, 8]
        v3 = [b'foo', None, b'bar', b'qux', np.nan]

        arrays = {
            'cat_strings': pd.Categorical(v1 * repeats),
            'cat_strings_with_na': pd.Categorical(v1 * repeats,
                                                  categories=['foo', 'bar']),
            'cat_ints': pd.Categorical(v2 * repeats),
            'cat_binary': pd.Categorical(v3 * repeats),
            'cat_strings_ordered': pd.Categorical(
                v1 * repeats, categories=['bar', 'qux', 'foo'],
                ordered=True),
            'ints': v2 * repeats,
            'ints2': v2 * repeats,
            'strings': v1 * repeats,
            'strings2': v1 * repeats,
            'strings3': v3 * repeats}
        df = pd.DataFrame(arrays)
        _check_pandas_roundtrip(df)

        for k in arrays:
            _check_array_roundtrip(arrays[k])

    def test_category_implicit_from_pandas(self):
        # ARROW-3374
        def _check(v):
            arr = pa.array(v)
            result = arr.to_pandas()
            tm.assert_series_equal(pd.Series(result), pd.Series(v))

        arrays = [
            pd.Categorical(['a', 'b', 'c'], categories=['a', 'b']),
            pd.Categorical(['a', 'b', 'c'], categories=['a', 'b'],
                           ordered=True)
        ]
        for arr in arrays:
            _check(arr)

    def test_empty_category(self):
        # ARROW-2443
        df = pd.DataFrame({'cat': pd.Categorical([])})
        _check_pandas_roundtrip(df)

    def test_category_zero_chunks(self):
        # ARROW-5952
        for pa_type, dtype in [(pa.string(), 'object'), (pa.int64(), 'int64')]:
            a = pa.chunked_array([], pa.dictionary(pa.int8(), pa_type))
            result = a.to_pandas()
            expected = pd.Categorical([], categories=np.array([], dtype=dtype))
            tm.assert_series_equal(pd.Series(result), pd.Series(expected))

            table = pa.table({'a': a})
            result = table.to_pandas()
            expected = pd.DataFrame({'a': expected})
            tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "data,error_type",
        [
            ({"a": ["a", 1, 2.0]}, pa.ArrowTypeError),
            ({"a": ["a", 1, 2.0]}, pa.ArrowTypeError),
            ({"a": [1, True]}, pa.ArrowTypeError),
            ({"a": [True, "a"]}, pa.ArrowInvalid),
            ({"a": [1, "a"]}, pa.ArrowInvalid),
            ({"a": [1.0, "a"]}, pa.ArrowInvalid),
        ],
    )
    def test_mixed_types_fails(self, data, error_type):
        df = pd.DataFrame(data)
        msg = "Conversion failed for column a with type object"
        with pytest.raises(error_type, match=msg):
            pa.Table.from_pandas(df)

    def test_strided_data_import(self):
        cases = []

        columns = ['a', 'b', 'c']
        N, K = 100, 3
        random_numbers = np.random.randn(N, K).copy() * 100

        numeric_dtypes = ['i1', 'i2', 'i4', 'i8', 'u1', 'u2', 'u4', 'u8',
                          'f4', 'f8']

        for type_name in numeric_dtypes:
            # Casting np.float64 -> uint32 or uint64 throws a RuntimeWarning
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                cases.append(random_numbers.astype(type_name))

        # strings
        cases.append(np.array([random_ascii(10) for i in range(N * K)],
                              dtype=object)
                     .reshape(N, K).copy())

        # booleans
        boolean_objects = (np.array([True, False, True] * N, dtype=object)
                           .reshape(N, K).copy())

        # add some nulls, so dtype comes back as objects
        boolean_objects[5] = None
        cases.append(boolean_objects)

        cases.append(np.arange("2016-01-01T00:00:00.001", N * K,
                               dtype='datetime64[ms]')
                     .reshape(N, K).copy())

        strided_mask = (random_numbers > 0).astype(bool)[:, 0]

        for case in cases:
            df = pd.DataFrame(case, columns=columns)
            col = df['a']

            _check_pandas_roundtrip(df)
            _check_array_roundtrip(col)
            _check_array_roundtrip(col, mask=strided_mask)

    def test_all_nones(self):
        def _check_series(s):
            converted = pa.array(s)
            assert isinstance(converted, pa.NullArray)
            assert len(converted) == 3
            assert converted.null_count == 3
            for item in converted:
                assert item is pa.NA

        _check_series(pd.Series([None] * 3, dtype=object))
        _check_series(pd.Series([np.nan] * 3, dtype=object))
        _check_series(pd.Series([None, np.nan, None], dtype=object))

    def test_partial_schema(self):
        data = OrderedDict([
            ('a', [0, 1, 2, 3, 4]),
            ('b', np.array([-10, -5, 0, 5, 10], dtype=np.int32)),
            ('c', [-10, -5, 0, 5, 10])
        ])
        df = pd.DataFrame(data)

        partial_schema = pa.schema([
            pa.field('c', pa.int64()),
            pa.field('a', pa.int64())
        ])

        _check_pandas_roundtrip(df, schema=partial_schema,
                                expected=df[['c', 'a']],
                                expected_schema=partial_schema)

    def test_table_batch_empty_dataframe(self):
        df = pd.DataFrame({})
        _check_pandas_roundtrip(df, preserve_index=None)
        _check_pandas_roundtrip(df, preserve_index=None, as_batch=True)

        expected = pd.DataFrame(columns=pd.Index([]))
        _check_pandas_roundtrip(df, expected, preserve_index=False)
        _check_pandas_roundtrip(df, expected, preserve_index=False, as_batch=True)

        df2 = pd.DataFrame({}, index=[0, 1, 2])
        _check_pandas_roundtrip(df2, preserve_index=True)
        _check_pandas_roundtrip(df2, as_batch=True, preserve_index=True)

    def test_convert_empty_table(self):
        arr = pa.array([], type=pa.int64())
        empty_objects = pd.Series(np.array([], dtype=object))
        tm.assert_series_equal(arr.to_pandas(),
                               pd.Series(np.array([], dtype=np.int64)))
        arr = pa.array([], type=pa.string())
        tm.assert_series_equal(arr.to_pandas(), empty_objects)
        arr = pa.array([], type=pa.list_(pa.int64()))
        tm.assert_series_equal(arr.to_pandas(), empty_objects)
        arr = pa.array([], type=pa.struct([pa.field('a', pa.int64())]))
        tm.assert_series_equal(arr.to_pandas(), empty_objects)

    def test_non_natural_stride(self):
        """
        ARROW-2172: converting from a Numpy array with a stride that's
        not a multiple of itemsize.
        """
        dtype = np.dtype([('x', np.int32), ('y', np.int16)])
        data = np.array([(42, -1), (-43, 2)], dtype=dtype)
        assert data.strides == (6,)
        arr = pa.array(data['x'], type=pa.int32())
        assert arr.to_pylist() == [42, -43]
        arr = pa.array(data['y'], type=pa.int16())
        assert arr.to_pylist() == [-1, 2]

    def test_array_from_strided_numpy_array(self):
        # ARROW-5651
        np_arr = np.arange(0, 10, dtype=np.float32)[1:-1:2]
        pa_arr = pa.array(np_arr, type=pa.float64())
        expected = pa.array([1.0, 3.0, 5.0, 7.0], type=pa.float64())
        pa_arr.equals(expected)

    def test_safe_unsafe_casts(self):
        # ARROW-2799
        df = pd.DataFrame({
            'A': list('abc'),
            'B': np.linspace(0, 1, 3)
        })

        schema = pa.schema([
            pa.field('A', pa.string()),
            pa.field('B', pa.int32())
        ])

        with pytest.raises(ValueError):
            pa.Table.from_pandas(df, schema=schema)

        table = pa.Table.from_pandas(df, schema=schema, safe=False)
        assert table.column('B').type == pa.int32()

    def test_error_sparse(self):
        # ARROW-2818
        try:
            df = pd.DataFrame({'a': pd.arrays.SparseArray([1, np.nan, 3])})
        except AttributeError:
            # pandas.arrays module introduced in pandas 0.24
            df = pd.DataFrame({'a': pd.SparseArray([1, np.nan, 3])})
        with pytest.raises(TypeError, match="Sparse pandas data"):
            pa.Table.from_pandas(df)


def test_safe_cast_from_float_with_nans_to_int():
    values = pd.Series([1, 2, None, 4])
    arr = pa.Array.from_pandas(values, type=pa.int32(), safe=True)
    expected = pa.array([1, 2, None, 4], type=pa.int32())
    assert arr.equals(expected)


def _fully_loaded_dataframe_example():
    index = pd.MultiIndex.from_arrays([
        pd.date_range('2000-01-01', periods=5).repeat(2),
        np.tile(np.array(['foo', 'bar'], dtype=object), 5)
    ])

    c1 = pd.date_range('2000-01-01', periods=10)
    data = {
        0: c1,
        1: c1.tz_localize('utc'),
        2: c1.tz_localize('US/Eastern'),
        3: c1[::2].tz_localize('utc').repeat(2).astype('category'),
        4: ['foo', 'bar'] * 5,
        5: pd.Series(['foo', 'bar'] * 5).astype('category').values,
        6: [True, False] * 5,
        7: np.random.randn(10),
        8: np.random.randint(0, 100, size=10),
        9: pd.period_range('2013', periods=10, freq='M'),
        10: pd.interval_range(start=1, freq=1, periods=10),
    }
    return pd.DataFrame(data, index=index)


@pytest.mark.parametrize('columns', ([b'foo'], ['foo']))
def test_roundtrip_with_bytes_unicode(columns):
    if Version("2.0.0") <= Version(pd.__version__) < Version("3.0.0"):
        # TODO: regression in pandas, hopefully fixed in next version
        # https://issues.apache.org/jira/browse/ARROW-18394
        # https://github.com/pandas-dev/pandas/issues/50127
        pytest.skip("Regression in pandas 2.0.0")

    df = pd.DataFrame(columns=columns)
    table1 = pa.Table.from_pandas(df)
    table2 = pa.Table.from_pandas(table1.to_pandas())
    assert table1.equals(table2)
    assert table1.schema.equals(table2.schema)
    assert table1.schema.metadata == table2.schema.metadata


def _pytime_from_micros(val):
    microseconds = val % 1000000
    val //= 1000000
    seconds = val % 60
    val //= 60
    minutes = val % 60
    hours = val // 60
    return time(hours, minutes, seconds, microseconds)


def _pytime_to_micros(pytime):
    return (pytime.hour * 3600000000 +
            pytime.minute * 60000000 +
            pytime.second * 1000000 +
            pytime.microsecond)


def test_convert_unsupported_type_error_message():
    # ARROW-1454

    # custom python objects
    class A:
        pass

    df = pd.DataFrame({'a': [A(), A()]})

    msg = 'Conversion failed for column a with type object'
    with pytest.raises(ValueError, match=msg):
        pa.Table.from_pandas(df)


# ----------------------------------------------------------------------
# Hypothesis tests


@h.given(past.arrays(past.pandas_compatible_types))
def test_array_to_pandas_roundtrip(arr):
    s = arr.to_pandas()
    restored = pa.array(s, type=arr.type, from_pandas=True)
    assert restored.equals(arr)


# ----------------------------------------------------------------------
# Test object deduplication in to_pandas


def _generate_dedup_example(nunique, repeats):
    unique_values = [rands(10) for i in range(nunique)]
    return unique_values * repeats


def _assert_nunique(obj, expected):
    assert len({id(x) for x in obj}) == expected


def test_to_pandas_deduplicate_strings_array_types():
    if _pandas_api.uses_string_dtype():
        pytest.skip(
            "pandas uses string dtype and not object dtype, keyword has no effect"
        )
    nunique = 100
    repeats = 10
    values = _generate_dedup_example(nunique, repeats)

    for arr in [pa.array(values, type=pa.binary()),
                pa.array(values, type=pa.utf8()),
                pa.chunked_array([values, values])]:
        _assert_nunique(arr.to_pandas(), nunique)
        _assert_nunique(arr.to_pandas(deduplicate_objects=False), len(arr))


def test_to_pandas_deduplicate_strings_table_types():
    if _pandas_api.uses_string_dtype():
        pytest.skip(
            "pandas uses string dtype and not object dtype, keyword has no effect"
        )
    nunique = 100
    repeats = 10
    values = _generate_dedup_example(nunique, repeats)

    arr = pa.array(values)
    rb = pa.RecordBatch.from_arrays([arr], ['foo'])
    tbl = pa.Table.from_batches([rb])

    for obj in [rb, tbl]:
        _assert_nunique(obj.to_pandas()['foo'], nunique)
        _assert_nunique(obj.to_pandas(deduplicate_objects=False)['foo'],
                        len(obj))


def test_to_pandas_deduplicate_integers_as_objects():
    nunique = 100
    repeats = 10

    # Python automatically interns smaller integers
    unique_values = list(np.random.randint(10000000, 1000000000, size=nunique))
    unique_values[nunique // 2] = None

    arr = pa.array(unique_values * repeats)

    _assert_nunique(arr.to_pandas(integer_object_nulls=True), nunique)
    _assert_nunique(arr.to_pandas(integer_object_nulls=True,
                                  deduplicate_objects=False),
                    # Account for None
                    (nunique - 1) * repeats + 1)


def test_to_pandas_deduplicate_date_time():
    nunique = 100
    repeats = 10

    unique_values = list(range(nunique))

    cases = [
        # raw type, array type, to_pandas options
        ('int32', 'date32', {'date_as_object': True}),
        ('int64', 'date64', {'date_as_object': True}),
        ('int32', 'time32[ms]', {}),
        ('int64', 'time64[us]', {})
    ]

    for raw_type, array_type, pandas_options in cases:
        raw_arr = pa.array(unique_values * repeats, type=raw_type)
        casted_arr = raw_arr.cast(array_type)

        _assert_nunique(casted_arr.to_pandas(**pandas_options),
                        nunique)
        _assert_nunique(casted_arr.to_pandas(deduplicate_objects=False,
                                             **pandas_options),
                        len(casted_arr))


# ---------------------------------------------------------------------

def test_table_from_pandas_checks_field_nullability():
    # ARROW-2136
    df = pd.DataFrame({'a': [1.2, 2.1, 3.1],
                       'b': [np.nan, 'string', 'foo']})
    schema = pa.schema([pa.field('a', pa.float64(), nullable=False),
                        pa.field('b', pa.utf8(), nullable=False)])

    with pytest.raises(ValueError):
        pa.Table.from_pandas(df, schema=schema)


def test_table_from_pandas_keeps_column_order_of_dataframe():
    df1 = pd.DataFrame(OrderedDict([
        ('partition', [0, 0, 1, 1]),
        ('arrays', [[0, 1, 2], [3, 4], None, None]),
        ('floats', [None, None, 1.1, 3.3])
    ]))
    df2 = df1[['floats', 'partition', 'arrays']]

    schema1 = pa.schema([
        ('partition', pa.int64()),
        ('arrays', pa.list_(pa.int64())),
        ('floats', pa.float64()),
    ])
    schema2 = pa.schema([
        ('floats', pa.float64()),
        ('partition', pa.int64()),
        ('arrays', pa.list_(pa.int64()))
    ])

    table1 = pa.Table.from_pandas(df1, preserve_index=False)
    table2 = pa.Table.from_pandas(df2, preserve_index=False)

    assert table1.schema.equals(schema1)
    assert table2.schema.equals(schema2)


def test_table_from_pandas_keeps_column_order_of_schema():
    # ARROW-3766
    df = pd.DataFrame(OrderedDict([
        ('partition', [0, 0, 1, 1]),
        ('arrays', [[0, 1, 2], [3, 4], None, None]),
        ('floats', [None, None, 1.1, 3.3])
    ]))

    schema = pa.schema([
        ('floats', pa.float64()),
        ('arrays', pa.list_(pa.int32())),
        ('partition', pa.int32())
    ])

    df1 = df[df.partition == 0]
    df2 = df[df.partition == 1][['floats', 'partition', 'arrays']]

    table1 = pa.Table.from_pandas(df1, schema=schema, preserve_index=False)
    table2 = pa.Table.from_pandas(df2, schema=schema, preserve_index=False)

    assert table1.schema.equals(schema)
    assert table1.schema.equals(table2.schema)


def test_table_from_pandas_columns_argument_only_does_filtering():
    df = pd.DataFrame(OrderedDict([
        ('partition', [0, 0, 1, 1]),
        ('arrays', [[0, 1, 2], [3, 4], None, None]),
        ('floats', [None, None, 1.1, 3.3])
    ]))

    columns1 = ['arrays', 'floats', 'partition']
    schema1 = pa.schema([
        ('arrays', pa.list_(pa.int64())),
        ('floats', pa.float64()),
        ('partition', pa.int64())
    ])

    columns2 = ['floats', 'partition']
    schema2 = pa.schema([
        ('floats', pa.float64()),
        ('partition', pa.int64())
    ])

    table1 = pa.Table.from_pandas(df, columns=columns1, preserve_index=False)
    table2 = pa.Table.from_pandas(df, columns=columns2, preserve_index=False)

    assert table1.schema.equals(schema1)
    assert table2.schema.equals(schema2)


def test_table_from_pandas_columns_and_schema_are_mutually_exclusive():
    df = pd.DataFrame(OrderedDict([
        ('partition', [0, 0, 1, 1]),
        ('arrays', [[0, 1, 2], [3, 4], None, None]),
        ('floats', [None, None, 1.1, 3.3])
    ]))
    schema = pa.schema([
        ('partition', pa.int32()),
        ('arrays', pa.list_(pa.int32())),
        ('floats', pa.float64()),
    ])
    columns = ['arrays', 'floats']

    with pytest.raises(ValueError):
        pa.Table.from_pandas(df, schema=schema, columns=columns)


def test_table_from_pandas_keeps_schema_nullability():
    # ARROW-5169
    df = pd.DataFrame({'a': [1, 2, 3, 4]})

    schema = pa.schema([
        pa.field('a', pa.int64(), nullable=False),
    ])

    table = pa.Table.from_pandas(df)
    assert table.schema.field('a').nullable is True
    table = pa.Table.from_pandas(df, schema=schema)
    assert table.schema.field('a').nullable is False


def test_table_from_pandas_schema_index_columns():
    # ARROW-5220
    df = pd.DataFrame({'a': [1, 2, 3], 'b': [0.1, 0.2, 0.3]})

    schema = pa.schema([
        ('a', pa.int64()),
        ('b', pa.float64()),
        ('index', pa.int64()),
    ])

    # schema includes index with name not in dataframe
    with pytest.raises(KeyError, match="name 'index' present in the"):
        pa.Table.from_pandas(df, schema=schema)

    df.index.name = 'index'

    # schema includes correct index name -> roundtrip works
    _check_pandas_roundtrip(df, schema=schema, preserve_index=True,
                            expected_schema=schema)

    # schema includes correct index name but preserve_index=False
    with pytest.raises(ValueError, match="'preserve_index=False' was"):
        pa.Table.from_pandas(df, schema=schema, preserve_index=False)

    # in case of preserve_index=None -> RangeIndex serialized as metadata
    # clashes with the index in the schema
    with pytest.raises(ValueError, match="name 'index' is present in the "
                                         "schema, but it is a RangeIndex"):
        pa.Table.from_pandas(df, schema=schema, preserve_index=None)

    df.index = pd.Index([0, 1, 2], name='index')

    # for non-RangeIndex, both preserve_index=None and True work
    _check_pandas_roundtrip(df, schema=schema, preserve_index=None,
                            expected_schema=schema)
    _check_pandas_roundtrip(df, schema=schema, preserve_index=True,
                            expected_schema=schema)

    # schema has different order (index column not at the end)
    schema = pa.schema([
        ('index', pa.int64()),
        ('a', pa.int64()),
        ('b', pa.float64()),
    ])
    _check_pandas_roundtrip(df, schema=schema, preserve_index=None,
                            expected_schema=schema)
    _check_pandas_roundtrip(df, schema=schema, preserve_index=True,
                            expected_schema=schema)

    # schema does not include the index -> index is not included as column
    # even though preserve_index=True/None
    schema = pa.schema([
        ('a', pa.int64()),
        ('b', pa.float64()),
    ])
    expected = df.copy()
    expected = expected.reset_index(drop=True)
    _check_pandas_roundtrip(df, schema=schema, preserve_index=None,
                            expected_schema=schema, expected=expected)
    _check_pandas_roundtrip(df, schema=schema, preserve_index=True,
                            expected_schema=schema, expected=expected)

    # dataframe with a MultiIndex
    df.index = pd.MultiIndex.from_tuples([('a', 1), ('a', 2), ('b', 1)],
                                         names=['level1', 'level2'])
    schema = pa.schema([
        ('level1', pa.string()),
        ('level2', pa.int64()),
        ('a', pa.int64()),
        ('b', pa.float64()),
    ])
    _check_pandas_roundtrip(df, schema=schema, preserve_index=True,
                            expected_schema=schema)
    _check_pandas_roundtrip(df, schema=schema, preserve_index=None,
                            expected_schema=schema)

    # only one of the levels of the MultiIndex is included
    schema = pa.schema([
        ('level2', pa.int64()),
        ('a', pa.int64()),
        ('b', pa.float64()),
    ])
    expected = df.copy()
    expected = expected.reset_index('level1', drop=True)
    _check_pandas_roundtrip(df, schema=schema, preserve_index=True,
                            expected_schema=schema, expected=expected)
    _check_pandas_roundtrip(df, schema=schema, preserve_index=None,
                            expected_schema=schema, expected=expected)


def test_table_from_pandas_schema_index_columns__unnamed_index():
    # ARROW-6999 - unnamed indices in specified schema
    df = pd.DataFrame({'a': [1, 2, 3], 'b': [0.1, 0.2, 0.3]})

    expected_schema = pa.schema([
        ('a', pa.int64()),
        ('b', pa.float64()),
        ('__index_level_0__', pa.int64()),
    ])

    schema = pa.Schema.from_pandas(df, preserve_index=True)
    table = pa.Table.from_pandas(df, preserve_index=True, schema=schema)
    assert table.schema.remove_metadata().equals(expected_schema)

    # non-RangeIndex (preserved by default)
    df = pd.DataFrame({'a': [1, 2, 3], 'b': [0.1, 0.2, 0.3]}, index=[0, 1, 2])
    schema = pa.Schema.from_pandas(df)
    table = pa.Table.from_pandas(df, schema=schema)
    assert table.schema.remove_metadata().equals(expected_schema)


def test_table_from_pandas_schema_with_custom_metadata():
    # ARROW-7087 - metadata disappear from pandas
    df = pd.DataFrame()
    schema = pa.Schema.from_pandas(df).with_metadata({'meta': 'True'})
    table = pa.Table.from_pandas(df, schema=schema)
    assert table.schema.metadata.get(b'meta') == b'True'


def test_table_from_pandas_schema_field_order_metadata():
    # ARROW-10532
    # ensure that a different field order in specified schema doesn't
    # mangle metadata
    df = pd.DataFrame({
        "datetime": pd.date_range("2020-01-01T00:00:00Z", freq="h", periods=2),
        "float": np.random.randn(2)
    })

    schema = pa.schema([
        pa.field("float", pa.float32(), nullable=True),
        pa.field("datetime", pa.timestamp("s", tz="UTC"), nullable=False)
    ])

    table = pa.Table.from_pandas(df, schema=schema)
    assert table.schema.equals(schema)
    metadata_float = table.schema.pandas_metadata["columns"][0]
    assert metadata_float["name"] == "float"
    assert metadata_float["metadata"] is None
    metadata_datetime = table.schema.pandas_metadata["columns"][1]
    assert metadata_datetime["name"] == "datetime"
    assert metadata_datetime["metadata"] == {'timezone': 'UTC'}

    result = table.to_pandas()
    coerce_cols_to_types = {"float": "float32"}
    if Version(pd.__version__) >= Version("2.0.0"):
        # Pandas v2 now support non-nanosecond time units
        coerce_cols_to_types["datetime"] = "datetime64[s, UTC]"
    expected = df[["float", "datetime"]].astype(coerce_cols_to_types)

    tm.assert_frame_equal(result, expected)


# ----------------------------------------------------------------------
# RecordBatch, Table


def test_recordbatch_from_to_pandas():
    data = pd.DataFrame({
        'c1': np.array([1, 2, 3, 4, 5], dtype='int64'),
        'c2': np.array([1, 2, 3, 4, 5], dtype='uint32'),
        'c3': np.random.randn(5),
        'c4': ['foo', 'bar', None, 'baz', 'qux'],
        'c5': [False, True, False, True, False]
    })

    batch = pa.RecordBatch.from_pandas(data)
    result = batch.to_pandas()
    tm.assert_frame_equal(data, result)


def test_recordbatchlist_to_pandas():
    data1 = pd.DataFrame({
        'c1': np.array([1, 1, 2], dtype='uint32'),
        'c2': np.array([1.0, 2.0, 3.0], dtype='float64'),
        'c3': [True, None, False],
        'c4': ['foo', 'bar', None]
    })

    data2 = pd.DataFrame({
        'c1': np.array([3, 5], dtype='uint32'),
        'c2': np.array([4.0, 5.0], dtype='float64'),
        'c3': [True, True],
        'c4': ['baz', 'qux']
    })

    batch1 = pa.RecordBatch.from_pandas(data1)
    batch2 = pa.RecordBatch.from_pandas(data2)

    table = pa.Table.from_batches([batch1, batch2])
    result = table.to_pandas()
    data = pd.concat([data1, data2]).reset_index(drop=True)
    tm.assert_frame_equal(data, result)


def test_recordbatch_table_pass_name_to_pandas():
    rb = pa.record_batch([pa.array([1, 2, 3, 4])], names=['a0'])
    t = pa.table([pa.array([1, 2, 3, 4])], names=['a0'])
    assert rb[0].to_pandas().name == 'a0'
    assert t[0].to_pandas().name == 'a0'


# ----------------------------------------------------------------------
# Metadata serialization


@pytest.mark.parametrize(
    ('type', 'expected'),
    [
        (pa.null(), 'empty'),
        (pa.bool_(), 'bool'),
        (pa.int8(), 'int8'),
        (pa.int16(), 'int16'),
        (pa.int32(), 'int32'),
        (pa.int64(), 'int64'),
        (pa.uint8(), 'uint8'),
        (pa.uint16(), 'uint16'),
        (pa.uint32(), 'uint32'),
        (pa.uint64(), 'uint64'),
        (pa.float16(), 'float16'),
        (pa.float32(), 'float32'),
        (pa.float64(), 'float64'),
        (pa.date32(), 'date'),
        (pa.date64(), 'date'),
        (pa.binary(), 'bytes'),
        (pa.binary(length=4), 'bytes'),
        (pa.string(), 'unicode'),
        (pa.list_(pa.list_(pa.int16())), 'list[list[int16]]'),
        (pa.decimal128(18, 3), 'decimal'),
        (pa.timestamp('ms'), 'datetime'),
        (pa.timestamp('us', 'UTC'), 'datetimetz'),
        (pa.time32('s'), 'time'),
        (pa.time64('us'), 'time')
    ]
)
def test_logical_type(type, expected):
    assert get_logical_type(type) == expected


# ----------------------------------------------------------------------
# to_pandas uses MemoryPool

def test_array_uses_memory_pool():
    # ARROW-6570
    N = 10000
    arr = pa.array(np.arange(N, dtype=np.int64),
                   mask=np.random.randint(0, 2, size=N).astype(np.bool_))

    # In the case the gc is caught loading
    gc.collect()

    prior_allocation = pa.total_allocated_bytes()

    x = arr.to_pandas()
    assert pa.total_allocated_bytes() == (prior_allocation + N * 8)
    x = None  # noqa
    gc.collect()

    assert pa.total_allocated_bytes() == prior_allocation

    # zero copy does not allocate memory
    arr = pa.array(np.arange(N, dtype=np.int64))

    prior_allocation = pa.total_allocated_bytes()
    x = arr.to_pandas()  # noqa
    assert pa.total_allocated_bytes() == prior_allocation


def test_singleton_blocks_zero_copy():
    # Part of ARROW-3789
    t = pa.table([pa.array(np.arange(1000, dtype=np.int64))], ['f0'])

    # Zero copy if split_blocks=True
    _check_to_pandas_memory_unchanged(t, split_blocks=True)

    prior_allocation = pa.total_allocated_bytes()
    result = t.to_pandas()
    # access private `_values` because the public `values` is made read-only by pandas
    assert result['f0']._values.flags.writeable
    assert pa.total_allocated_bytes() > prior_allocation


def _check_to_pandas_memory_unchanged(obj, **kwargs):
    prior_allocation = pa.total_allocated_bytes()
    x = obj.to_pandas(**kwargs)  # noqa

    # Memory allocation unchanged -- either zero copy or self-destructing
    if _pandas_api.uses_string_dtype():
        # for the string array of the columns Index
        # -> increase the size to account for overallocation for small arrays
        max_index_allocation = max(192, x.columns.nbytes * 2)
        assert pa.total_allocated_bytes() <= (prior_allocation + max_index_allocation)
    else:
        assert pa.total_allocated_bytes() == prior_allocation


def test_to_pandas_split_blocks():
    # ARROW-3789
    t = pa.table([
        pa.array([1, 2, 3, 4, 5]*100, type='i1'),
        pa.array([1, 2, 3, 4, 5]*100, type='i4'),
        pa.array([1, 2, 3, 4, 5]*100, type='i8'),
        pa.array([1, 2, 3, 4, 5]*100, type='f4'),
        pa.array([1, 2, 3, 4, 5]*100, type='f8'),
        pa.array([1, 2, 3, 4, 5]*100, type='f8'),
        pa.array([1, 2, 3, 4, 5]*100, type='f8'),
        pa.array([1, 2, 3, 4, 5]*100, type='f8'),
    ], [f'f{i}' for i in range(8)])

    _check_blocks_created(t, 8)
    _check_to_pandas_memory_unchanged(t, split_blocks=True)


def _get_mgr(df):
    if Version(pd.__version__) < Version("1.1.0"):
        return df._data
    else:
        return df._mgr


def _check_blocks_created(t, number):
    x = t.to_pandas(split_blocks=True)
    assert len(_get_mgr(x).blocks) == number


def test_to_pandas_self_destruct():
    K = 50

    def _make_table():
        return pa.table([
            # Slice to force a copy
            pa.array(np.random.randn(10000)[::2])
            for i in range(K)
        ], [f'f{i}' for i in range(K)])

    t = _make_table()
    _check_to_pandas_memory_unchanged(t, split_blocks=True, self_destruct=True)

    # Check non-split-block behavior
    t = _make_table()
    _check_to_pandas_memory_unchanged(t, self_destruct=True)


def test_table_uses_memory_pool():
    N = 10000
    arr = pa.array(np.arange(N, dtype=np.int64))
    t = pa.table([arr, arr, arr], ['f0', 'f1', 'f2'])

    prior_allocation = pa.total_allocated_bytes()
    x = t.to_pandas()

    new_allocation = 3 * N * 8
    if _pandas_api.uses_string_dtype():
        # for the small columns Index
        new_allocation += 128

    assert pa.total_allocated_bytes() == (prior_allocation + new_allocation)

    # Check successful garbage collection
    x = None  # noqa
    gc.collect()
    assert pa.total_allocated_bytes() == prior_allocation


def test_object_leak_in_numpy_array():
    # ARROW-6876
    arr = pa.array([{'a': 1}])
    np_arr = arr.to_pandas()
    assert np_arr.dtype == np.dtype('object')
    obj = np_arr[0]
    refcount = sys.getrefcount(obj)
    assert sys.getrefcount(obj) == refcount
    del np_arr
    assert sys.getrefcount(obj) == refcount - 1


def test_object_leak_in_dataframe():
    # ARROW-6876
    arr = pa.array([{'a': 1}])
    table = pa.table([arr], ['f0'])
    col = table.to_pandas()['f0']
    assert col.dtype == np.dtype('object')
    obj = col[0]
    refcount = sys.getrefcount(obj)
    assert sys.getrefcount(obj) == refcount
    del col
    assert sys.getrefcount(obj) == refcount - 1


# ----------------------------------------------------------------------
# Some nested array tests array tests


def test_array_from_py_float32():
    data = [[1.2, 3.4], [9.0, 42.0]]

    t = pa.float32()

    arr1 = pa.array(data[0], type=t)
    arr2 = pa.array(data, type=pa.list_(t))

    expected1 = np.array(data[0], dtype=np.float32)
    expected2 = pd.Series([np.array(data[0], dtype=np.float32),
                           np.array(data[1], dtype=np.float32)])

    assert arr1.type == t
    assert arr1.equals(pa.array(expected1))
    assert arr2.equals(pa.array(expected2))


# ----------------------------------------------------------------------
# Timestamp tests


def test_cast_timestamp_unit():
    # ARROW-1680
    val = datetime.now()
    s = pd.Series([val])
    s_nyc = s.dt.tz_localize('tzlocal()').dt.tz_convert('America/New_York')

    us_with_tz = pa.timestamp('us', tz='America/New_York')

    arr = pa.Array.from_pandas(s_nyc, type=us_with_tz)

    # ARROW-1906
    assert arr.type == us_with_tz

    arr2 = pa.Array.from_pandas(s, type=pa.timestamp('us'))

    assert arr[0].as_py() == s_nyc[0].to_pydatetime()
    assert arr2[0].as_py() == s[0].to_pydatetime()

    # Disallow truncation
    arr = pa.array([123123], type='int64').cast(pa.timestamp('ms'))
    expected = pa.array([123], type='int64').cast(pa.timestamp('s'))

    # sanity check that the cast worked right
    assert arr.type == pa.timestamp('ms')

    target = pa.timestamp('s')
    with pytest.raises(ValueError):
        arr.cast(target)

    result = arr.cast(target, safe=False)
    assert result.equals(expected)

    # ARROW-1949
    series = pd.Series([pd.Timestamp(1), pd.Timestamp(10), pd.Timestamp(1000)])
    expected = pa.array([0, 0, 1], type=pa.timestamp('us'))

    with pytest.raises(ValueError):
        pa.array(series, type=pa.timestamp('us'))

    with pytest.raises(ValueError):
        pa.Array.from_pandas(series, type=pa.timestamp('us'))

    result = pa.Array.from_pandas(series, type=pa.timestamp('us'), safe=False)
    assert result.equals(expected)

    result = pa.array(series, type=pa.timestamp('us'), safe=False)
    assert result.equals(expected)


def test_nested_with_timestamp_tz_round_trip():
    ts = pd.Timestamp.now()
    ts_dt = ts.to_pydatetime()
    arr = pa.array([ts_dt], type=pa.timestamp('us', tz='America/New_York'))
    struct = pa.StructArray.from_arrays([arr, arr], ['start', 'stop'])

    result = struct.to_pandas()
    restored = pa.array(result)
    assert restored.equals(struct)


def test_nested_with_timestamp_tz():
    # ARROW-7723
    ts = pd.Timestamp.now()
    ts_dt = ts.to_pydatetime()

    # XXX: Ensure that this data does not get promoted to nanoseconds (and thus
    # integers) to preserve behavior in 0.15.1
    for unit in ['s', 'ms', 'us']:
        if unit in ['s', 'ms']:
            # This is used for verifying timezone conversion to micros are not
            # important
            def truncate(x): return x.replace(microsecond=0)
        else:
            def truncate(x): return x
        arr = pa.array([ts], type=pa.timestamp(unit))
        arr2 = pa.array([ts], type=pa.timestamp(unit, tz='America/New_York'))

        arr3 = pa.StructArray.from_arrays([arr, arr], ['start', 'stop'])
        arr4 = pa.StructArray.from_arrays([arr2, arr2], ['start', 'stop'])

        result = arr3.to_pandas()
        assert isinstance(result[0]['start'], datetime)
        assert result[0]['start'].tzinfo is None
        assert isinstance(result[0]['stop'], datetime)
        assert result[0]['stop'].tzinfo is None

        result = arr4.to_pandas()
        assert isinstance(result[0]['start'], datetime)
        assert result[0]['start'].tzinfo is not None
        utc_dt = result[0]['start'].astimezone(timezone.utc)
        assert truncate(utc_dt).replace(tzinfo=None) == truncate(ts_dt)
        assert isinstance(result[0]['stop'], datetime)
        assert result[0]['stop'].tzinfo is not None

        # same conversion for table
        result = pa.table({'a': arr3}).to_pandas()
        assert isinstance(result['a'][0]['start'], datetime)
        assert result['a'][0]['start'].tzinfo is None
        assert isinstance(result['a'][0]['stop'], datetime)
        assert result['a'][0]['stop'].tzinfo is None

        result = pa.table({'a': arr4}).to_pandas()
        assert isinstance(result['a'][0]['start'], datetime)
        assert result['a'][0]['start'].tzinfo is not None
        assert isinstance(result['a'][0]['stop'], datetime)
        assert result['a'][0]['stop'].tzinfo is not None


# ----------------------------------------------------------------------
# DictionaryArray tests


def test_dictionary_with_pandas():
    src_indices = np.repeat([0, 1, 2], 2)
    dictionary = np.array(['foo', 'bar', 'baz'], dtype=object)
    mask = np.array([False, False, True, False, False, False])

    for index_type in ['uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32',
                       'uint64', 'int64']:
        indices = src_indices.astype(index_type)
        d1 = pa.DictionaryArray.from_arrays(indices, dictionary)
        d2 = pa.DictionaryArray.from_arrays(indices, dictionary, mask=mask)

        if index_type == 'uint64':
            # uint64 is not supported due to overflow risk (values > 2^63-1)
            with pytest.raises(TypeError,
                               match="UInt64 dictionary indices"):
                d1.to_pandas()
            continue

        pandas1 = d1.to_pandas()
        # Pandas Categorical uses signed int codes. Arrow converts:
        # uint8 to int16, uint16 to int32, uint32 to int64, signed types unchanged
        if index_type == 'uint8':
            compare_indices = indices.astype('int16')
        elif index_type == 'uint16':
            compare_indices = indices.astype('int32')
        elif index_type == 'uint32':
            compare_indices = indices.astype('int64')
        else:
            compare_indices = indices
        ex_pandas1 = pd.Categorical.from_codes(compare_indices, categories=dictionary)

        tm.assert_series_equal(pd.Series(pandas1), pd.Series(ex_pandas1))

        pandas2 = d2.to_pandas()
        assert pandas2.isnull().sum() == 1

        # Use same conversion as above for comparison
        if index_type == 'uint8':
            signed_indices = indices.astype('int16')
        elif index_type == 'uint16':
            signed_indices = indices.astype('int32')
        elif index_type == 'uint32':
            signed_indices = indices.astype('int64')
        else:
            signed_indices = indices
        ex_pandas2 = pd.Categorical.from_codes(np.where(mask, -1,
                                                        signed_indices),
                                               categories=dictionary)

        tm.assert_series_equal(pd.Series(pandas2), pd.Series(ex_pandas2))


def random_strings(n, item_size, pct_null=0, dictionary=None):
    if dictionary is not None:
        result = dictionary[np.random.randint(0, len(dictionary), size=n)]
    else:
        result = np.array([random_ascii(item_size) for i in range(n)],
                          dtype=object)

    if pct_null > 0:
        result[np.random.rand(n) < pct_null] = None

    return result


def test_variable_dictionary_to_pandas():
    np.random.seed(12345)

    d1 = pa.array(random_strings(100, 32), type='string')
    d2 = pa.array(random_strings(100, 16), type='string')
    d3 = pa.array(random_strings(10000, 10), type='string')

    a1 = pa.DictionaryArray.from_arrays(
        np.random.randint(0, len(d1), size=1000, dtype='i4'),
        d1
    )
    a2 = pa.DictionaryArray.from_arrays(
        np.random.randint(0, len(d2), size=1000, dtype='i4'),
        d2
    )

    # With some nulls
    a3 = pa.DictionaryArray.from_arrays(
        np.random.randint(0, len(d3), size=1000, dtype='i4'), d3)

    i4 = pa.array(
        np.random.randint(0, len(d3), size=1000, dtype='i4'),
        mask=np.random.rand(1000) < 0.1
    )
    a4 = pa.DictionaryArray.from_arrays(i4, d3)

    expected_dict = pa.concat_arrays([d1, d2, d3])

    a = pa.chunked_array([a1, a2, a3, a4])
    a_dense = pa.chunked_array([a1.cast('string'),
                                a2.cast('string'),
                                a3.cast('string'),
                                a4.cast('string')])

    result = a.to_pandas()
    result_dense = a_dense.to_pandas()

    assert (result.cat.categories == expected_dict.to_pandas()).all()

    expected_dense = result.astype('str')
    expected_dense[result_dense.isnull()] = None
    tm.assert_series_equal(result_dense, expected_dense)


def test_dictionary_encoded_nested_to_pandas():
    # ARROW-6899
    child = pa.array(['a', 'a', 'a', 'b', 'b']).dictionary_encode()

    arr = pa.ListArray.from_arrays([0, 3, 5], child)

    result = arr.to_pandas()
    expected = pd.Series([np.array(['a', 'a', 'a'], dtype=object),
                          np.array(['b', 'b'], dtype=object)])

    tm.assert_series_equal(result, expected)


def test_dictionary_from_pandas():
    cat = pd.Categorical(['a', 'b', 'a'])
    expected_type = pa.dictionary(
        pa.int8(),
        pa.large_string() if _pandas_api.uses_string_dtype() else pa.string()
    )

    result = pa.array(cat)
    assert result.to_pylist() == ['a', 'b', 'a']
    assert result.type.equals(expected_type)

    # with missing values in categorical
    cat = pd.Categorical(['a', 'b', None, 'a'])

    result = pa.array(cat)
    assert result.to_pylist() == ['a', 'b', None, 'a']
    assert result.type.equals(expected_type)

    # with additional mask
    result = pa.array(cat, mask=np.array([False, False, False, True]))
    assert result.to_pylist() == ['a', 'b', None, None]
    assert result.type.equals(expected_type)


def test_dictionary_from_pandas_specified_type():
    # ARROW-7168 - ensure specified type is always respected

    # the same as cat = pd.Categorical(['a', 'b']) but explicit about dtypes
    cat = pd.Categorical.from_codes(
        np.array([0, 1], dtype='int8'), np.array(['a', 'b'], dtype=object))

    # different index type -> allow this
    # (the type of the 'codes' in pandas is not part of the data type)
    typ = pa.dictionary(index_type=pa.int16(), value_type=pa.string())
    result = pa.array(cat, type=typ)
    assert result.type.equals(typ)
    assert result.to_pylist() == ['a', 'b']

    # mismatching values type -> raise error
    typ = pa.dictionary(index_type=pa.int8(), value_type=pa.int64())
    with pytest.raises(pa.ArrowInvalid):
        result = pa.array(cat, type=typ)

    # mismatching order -> raise error
    typ = pa.dictionary(
        index_type=pa.int8(), value_type=pa.string(), ordered=True)
    msg = "The 'ordered' flag of the passed categorical values "
    with pytest.raises(ValueError, match=msg):
        result = pa.array(cat, type=typ)
    assert result.to_pylist() == ['a', 'b']

    # with mask
    typ = pa.dictionary(index_type=pa.int16(), value_type=pa.string())
    result = pa.array(cat, type=typ, mask=np.array([False, True]))
    assert result.type.equals(typ)
    assert result.to_pylist() == ['a', None]

    # empty categorical -> be flexible in values type to allow
    cat = pd.Categorical([])

    typ = pa.dictionary(index_type=pa.int8(), value_type=pa.string())
    result = pa.array(cat, type=typ)
    assert result.type.equals(typ)
    assert result.to_pylist() == []
    typ = pa.dictionary(index_type=pa.int8(), value_type=pa.int64())
    result = pa.array(cat, type=typ)
    assert result.type.equals(typ)
    assert result.to_pylist() == []

    # passing non-dictionary type
    cat = pd.Categorical(['a', 'b'])
    result = pa.array(cat, type=pa.string())
    expected = pa.array(['a', 'b'], type=pa.string())
    assert result.equals(expected)
    assert result.to_pylist() == ['a', 'b']


def test_convert_categories_to_array_with_string_pyarrow_dtype():
    # gh-33727: categories should be converted to pa.Array
    if Version(pd.__version__) < Version("1.3.0"):
        pytest.skip("PyArrow backed string data type introduced in pandas 1.3.0")

    df = pd.DataFrame({"x": ["foo", "bar", "foo"]}, dtype="string[pyarrow]")
    df = df.astype("category")
    indices = pa.array(df['x'].cat.codes)
    dictionary = pa.array(df["x"].cat.categories.values)
    assert isinstance(dictionary, pa.Array)

    expected = pa.Array.from_pandas(df['x'])
    result = pa.DictionaryArray.from_arrays(indices, dictionary)
    assert result == expected


# ----------------------------------------------------------------------
# Array protocol in pandas conversions tests


def test_array_protocol():
    df = pd.DataFrame({'a': pd.Series([1, 2, None], dtype='Int64')})

    # __arrow_array__ added to pandas IntegerArray in 0.26.0.dev

    # default conversion
    result = pa.table(df)
    expected = pa.array([1, 2, None], pa.int64())
    assert result[0].chunk(0).equals(expected)

    # with specifying schema
    schema = pa.schema([('a', pa.float64())])
    result = pa.table(df, schema=schema)
    expected2 = pa.array([1, 2, None], pa.float64())
    assert result[0].chunk(0).equals(expected2)

    # pass Series to pa.array
    result = pa.array(df['a'])
    assert result.equals(expected)
    result = pa.array(df['a'], type=pa.float64())
    assert result.equals(expected2)

    # pass actual ExtensionArray to pa.array
    result = pa.array(df['a'].values)
    assert result.equals(expected)
    result = pa.array(df['a'].values, type=pa.float64())
    assert result.equals(expected2)


class DummyExtensionType(pa.ExtensionType):

    def __init__(self):
        super().__init__(pa.int64(),
                         'pyarrow.tests.test_pandas.DummyExtensionType')

    def __arrow_ext_serialize__(self):
        return b''

    @classmethod
    def __arrow_ext_deserialize__(cls, storage_type, serialized):
        assert serialized == b''
        assert storage_type == pa.int64()
        return cls()


def PandasArray__arrow_array__(self, type=None):
    # hardcode dummy return regardless of self - we only want to check that
    # this method is correctly called
    storage = pa.array([1, 2, 3], type=pa.int64())
    return pa.ExtensionArray.from_storage(DummyExtensionType(), storage)


def test_array_protocol_pandas_extension_types(monkeypatch):
    # ARROW-7022 - ensure protocol works for Period / Interval extension dtypes

    storage = pa.array([1, 2, 3], type=pa.int64())
    expected = pa.ExtensionArray.from_storage(DummyExtensionType(), storage)

    monkeypatch.setattr(pd.arrays.PeriodArray, "__arrow_array__",
                        PandasArray__arrow_array__, raising=False)
    monkeypatch.setattr(pd.arrays.IntervalArray, "__arrow_array__",
                        PandasArray__arrow_array__, raising=False)
    for arr in [pd.period_range("2012-01-01", periods=3, freq="D").array,
                pd.interval_range(1, 4).array]:
        result = pa.array(arr)
        assert result.equals(expected)
        result = pa.array(pd.Series(arr))
        assert result.equals(expected)
        result = pa.array(pd.Index(arr))
        assert result.equals(expected)
        result = pa.table(pd.DataFrame({'a': arr})).column('a').chunk(0)
        assert result.equals(expected)


# ----------------------------------------------------------------------
# Pandas ExtensionArray support


def _Int64Dtype__from_arrow__(self, array):
    # for test only deal with single chunk for now
    # TODO: do we require handling of chunked arrays in the protocol?
    if isinstance(array, pa.Array):
        arr = array
    else:
        # ChunkedArray - here only deal with a single chunk for the test
        arr = array.chunk(0)
    buflist = arr.buffers()
    data = np.frombuffer(buflist[-1], dtype='int64')[
        arr.offset:arr.offset + len(arr)]
    bitmask = buflist[0]
    if bitmask is not None:
        mask = pa.BooleanArray.from_buffers(
            pa.bool_(), len(arr), [None, bitmask])
        mask = np.asarray(mask)
    else:
        mask = np.ones(len(arr), dtype=bool)
    int_arr = pd.arrays.IntegerArray(data.copy(), ~mask, copy=False)
    return int_arr


def test_convert_to_extension_array(monkeypatch):
    # table converted from dataframe with extension types (so pandas_metadata
    # has this information)
    df = pd.DataFrame(
        {'a': [1, 2, 3], 'b': pd.array([2, 3, 4], dtype='Int64'),
         'c': [4, 5, 6]})
    table = pa.table(df)

    # Int64Dtype is recognized -> convert to extension block by default
    # for a proper roundtrip
    result = table.to_pandas()
    assert _get_mgr(result).blocks[0].values.dtype == np.dtype("int64")
    assert _get_mgr(result).blocks[1].values.dtype == pd.Int64Dtype()
    tm.assert_frame_equal(result, df)

    # test with missing values
    df2 = pd.DataFrame({'a': pd.array([1, 2, None], dtype='Int64')})
    table2 = pa.table(df2)
    result = table2.to_pandas()
    assert _get_mgr(result).blocks[0].values.dtype == pd.Int64Dtype()
    tm.assert_frame_equal(result, df2)

    # monkeypatch pandas Int64Dtype to *not* have the protocol method
    if Version(pd.__version__) < Version("1.3.0.dev"):
        monkeypatch.delattr(
            pd.core.arrays.integer._IntegerDtype, "__from_arrow__")
    else:
        monkeypatch.delattr(
            pd.core.arrays.integer.NumericDtype, "__from_arrow__")
    # Int64Dtype has no __from_arrow__ -> use normal conversion
    result = table.to_pandas()
    assert len(_get_mgr(result).blocks) == 1
    assert _get_mgr(result).blocks[0].values.dtype == np.dtype("int64")


class MyCustomIntegerType(pa.ExtensionType):

    def __init__(self):
        super().__init__(pa.int64(),
                         'pyarrow.tests.test_pandas.MyCustomIntegerType')

    def __arrow_ext_serialize__(self):
        return b''

    def to_pandas_dtype(self):
        return pd.Int64Dtype()


def test_conversion_extensiontype_to_extensionarray(monkeypatch):
    # converting extension type to linked pandas ExtensionDtype/Array
    storage = pa.array([1, 2, 3, 4], pa.int64())
    arr = pa.ExtensionArray.from_storage(MyCustomIntegerType(), storage)
    table = pa.table({'a': arr})

    # extension type points to Int64Dtype, which knows how to create a
    # pandas ExtensionArray
    result = arr.to_pandas()
    assert _get_mgr(result).blocks[0].values.dtype == pd.Int64Dtype()
    expected = pd.Series([1, 2, 3, 4], dtype='Int64')
    tm.assert_series_equal(result, expected)

    result = table.to_pandas()
    assert _get_mgr(result).blocks[0].values.dtype == pd.Int64Dtype()
    expected = pd.DataFrame({'a': pd.array([1, 2, 3, 4], dtype='Int64')})
    tm.assert_frame_equal(result, expected)

    # monkeypatch pandas Int64Dtype to *not* have the protocol method
    # (remove the version added above and the actual version for recent pandas)
    if Version(pd.__version__) < Version("1.3.0.dev"):
        monkeypatch.delattr(
            pd.core.arrays.integer._IntegerDtype, "__from_arrow__")
    else:
        monkeypatch.delattr(
            pd.core.arrays.integer.NumericDtype, "__from_arrow__")

    result = arr.to_pandas()
    assert _get_mgr(result).blocks[0].values.dtype == np.dtype("int64")
    expected = pd.Series([1, 2, 3, 4])
    tm.assert_series_equal(result, expected)

    with pytest.raises(ValueError):
        table.to_pandas()


def test_to_pandas_extension_dtypes_mapping():
    table = pa.table({'a': pa.array([1, 2, 3], pa.int64())})

    # default use numpy dtype
    result = table.to_pandas()
    assert result['a'].dtype == np.dtype('int64')

    # specify to override the default
    result = table.to_pandas(types_mapper={pa.int64(): pd.Int64Dtype()}.get)
    assert isinstance(result['a'].dtype, pd.Int64Dtype)

    # types that return None in function get normal conversion
    table = pa.table({'a': pa.array([1, 2, 3], pa.int32())})
    result = table.to_pandas(types_mapper={pa.int64(): pd.Int64Dtype()}.get)
    assert result['a'].dtype == np.dtype('int32')

    # `types_mapper` overrules the pandas metadata
    table = pa.table(pd.DataFrame({'a': pd.array([1, 2, 3], dtype="Int64")}))
    result = table.to_pandas()
    assert isinstance(result['a'].dtype, pd.Int64Dtype)
    result = table.to_pandas(
        types_mapper={pa.int64(): pd.PeriodDtype('D')}.get)
    assert isinstance(result['a'].dtype, pd.PeriodDtype)


def test_to_pandas_extension_dtypes_mapping_complex_type():
    # https://github.com/apache/arrow/pull/44720
    if Version(pd.__version__) < Version("1.5.2"):
        pytest.skip("Test relies on pd.ArrowDtype")
    pa_type = pa.struct(
        [
            pa.field("bar", pa.bool_(), nullable=False),
            pa.field("baz", pa.float32(), nullable=True),
        ],
    )
    pd_type = pd.ArrowDtype(pa_type)
    schema = pa.schema([pa.field("foo", pa_type)])
    df0 = pd.DataFrame(
        [
            {"foo": {"bar": True, "baz": np.float32(1)}},
            {"foo": {"bar": True, "baz": None}},
        ],
    ).astype({"foo": pd_type})

    # Round trip df0 into df1
    table = pa.Table.from_pandas(df0, schema=schema)
    df1 = table.to_pandas(types_mapper=pd.ArrowDtype)
    pd.testing.assert_frame_equal(df0, df1)


def test_array_to_pandas():
    if Version(pd.__version__) < Version("1.1"):
        pytest.skip("ExtensionDtype to_pandas method missing")

    for arr in [pd.period_range("2012-01-01", periods=3, freq="D").array,
                pd.interval_range(1, 4).array]:
        result = pa.array(arr).to_pandas()
        expected = pd.Series(arr)
        tm.assert_series_equal(result, expected)

        result = pa.table({"col": arr})["col"].to_pandas()
        expected = pd.Series(arr, name="col")
        tm.assert_series_equal(result, expected)


def test_roundtrip_empty_table_with_extension_dtype_index():
    df = pd.DataFrame(index=pd.interval_range(start=0, end=3))
    table = pa.table(df)
    if Version(pd.__version__) > Version("1.0"):
        tm.assert_index_equal(table.to_pandas().index, df.index)
    else:
        tm.assert_index_equal(table.to_pandas().index,
                              pd.Index([{'left': 0, 'right': 1},
                                        {'left': 1, 'right': 2},
                                        {'left': 2, 'right': 3}],
                                       dtype='object'))


@pytest.mark.parametrize("index", ["a", ["a", "b"]])
def test_to_pandas_types_mapper_index(index):
    if Version(pd.__version__) < Version("1.5.0"):
        pytest.skip("ArrowDtype missing")
    df = pd.DataFrame(
        {
            "a": [1, 2],
            "b": [3, 4],
            "c": [5, 6],
        },
        dtype=pd.ArrowDtype(pa.int64()),
    ).set_index(index)
    expected = df.copy()
    table = pa.table(df)
    result = table.to_pandas(types_mapper=pd.ArrowDtype)
    tm.assert_frame_equal(result, expected)


def test_array_to_pandas_types_mapper():
    # https://issues.apache.org/jira/browse/ARROW-9664
    if Version(pd.__version__) < Version("1.2.0"):
        pytest.skip("Float64Dtype extension dtype missing")

    data = pa.array([1, 2, 3], pa.int64())

    # Test with mapper function
    types_mapper = {pa.int64(): pd.Int64Dtype()}.get
    result = data.to_pandas(types_mapper=types_mapper)
    assert result.dtype == pd.Int64Dtype()

    # Test mapper function returning None
    types_mapper = {pa.int64(): None}.get
    result = data.to_pandas(types_mapper=types_mapper)
    assert result.dtype == np.dtype("int64")

    # Test mapper function not containing the dtype
    types_mapper = {pa.float64(): pd.Float64Dtype()}.get
    result = data.to_pandas(types_mapper=types_mapper)
    assert result.dtype == np.dtype("int64")


@pytest.mark.pandas
def test_chunked_array_to_pandas_types_mapper():
    # https://issues.apache.org/jira/browse/ARROW-9664
    if Version(pd.__version__) < Version("1.2.0"):
        pytest.skip("Float64Dtype extension dtype missing")

    data = pa.chunked_array([pa.array([1, 2, 3], pa.int64())])
    assert isinstance(data, pa.ChunkedArray)

    # Test with mapper function
    types_mapper = {pa.int64(): pd.Int64Dtype()}.get
    result = data.to_pandas(types_mapper=types_mapper)
    assert result.dtype == pd.Int64Dtype()

    # Test mapper function returning None
    types_mapper = {pa.int64(): None}.get
    result = data.to_pandas(types_mapper=types_mapper)
    assert result.dtype == np.dtype("int64")

    # Test mapper function not containing the dtype
    types_mapper = {pa.float64(): pd.Float64Dtype()}.get
    result = data.to_pandas(types_mapper=types_mapper)
    assert result.dtype == np.dtype("int64")


# ----------------------------------------------------------------------
# Legacy metadata compatibility tests


def test_metadata_compat_range_index_pre_0_12():
    # Forward compatibility for metadata created from pandas.RangeIndex
    # prior to pyarrow 0.13.0
    a_values = ['foo', 'bar', None, 'baz']
    b_values = ['a', 'a', 'b', 'b']
    a_arrow = pa.array(a_values, type='utf8')
    b_arrow = pa.array(b_values, type='utf8')

    rng_index_arrow = pa.array([0, 2, 4, 6], type='int64')

    gen_name_0 = '__index_level_0__'
    gen_name_1 = '__index_level_1__'

    # Case 1: named RangeIndex
    e1 = pd.DataFrame({
        'a': a_values
    }, index=pd.RangeIndex(0, 8, step=2, name='qux'))
    t1 = pa.Table.from_arrays([a_arrow, rng_index_arrow],
                              names=['a', 'qux'])
    t1 = t1.replace_schema_metadata({
        b'pandas': json.dumps(
            {'index_columns': ['qux'],
             'column_indexes': [{'name': None,
                                 'field_name': None,
                                 'pandas_type': 'unicode',
                                 'numpy_type': 'object',
                                 'metadata': {'encoding': 'UTF-8'}}],
             'columns': [{'name': 'a',
                          'field_name': 'a',
                          'pandas_type': 'unicode',
                          'numpy_type': 'object',
                          'metadata': None},
                         {'name': 'qux',
                          'field_name': 'qux',
                          'pandas_type': 'int64',
                          'numpy_type': 'int64',
                          'metadata': None}],
             'pandas_version': '0.23.4'}
        )})
    r1 = t1.to_pandas()
    tm.assert_frame_equal(r1, e1)

    # Case 2: named RangeIndex, but conflicts with an actual column
    e2 = pd.DataFrame({
        'qux': a_values
    }, index=pd.RangeIndex(0, 8, step=2, name='qux'))
    t2 = pa.Table.from_arrays([a_arrow, rng_index_arrow],
                              names=['qux', gen_name_0])
    t2 = t2.replace_schema_metadata({
        b'pandas': json.dumps(
            {'index_columns': [gen_name_0],
             'column_indexes': [{'name': None,
                                 'field_name': None,
                                 'pandas_type': 'unicode',
                                 'numpy_type': 'object',
                                 'metadata': {'encoding': 'UTF-8'}}],
             'columns': [{'name': 'a',
                          'field_name': 'a',
                          'pandas_type': 'unicode',
                          'numpy_type': 'object',
                          'metadata': None},
                         {'name': 'qux',
                          'field_name': gen_name_0,
                          'pandas_type': 'int64',
                          'numpy_type': 'int64',
                          'metadata': None}],
             'pandas_version': '0.23.4'}
        )})
    r2 = t2.to_pandas()
    tm.assert_frame_equal(r2, e2)

    # Case 3: unnamed RangeIndex
    e3 = pd.DataFrame({
        'a': a_values
    }, index=pd.RangeIndex(0, 8, step=2, name=None))
    t3 = pa.Table.from_arrays([a_arrow, rng_index_arrow],
                              names=['a', gen_name_0])
    t3 = t3.replace_schema_metadata({
        b'pandas': json.dumps(
            {'index_columns': [gen_name_0],
             'column_indexes': [{'name': None,
                                 'field_name': None,
                                 'pandas_type': 'unicode',
                                 'numpy_type': 'object',
                                 'metadata': {'encoding': 'UTF-8'}}],
             'columns': [{'name': 'a',
                          'field_name': 'a',
                          'pandas_type': 'unicode',
                          'numpy_type': 'object',
                          'metadata': None},
                         {'name': None,
                          'field_name': gen_name_0,
                          'pandas_type': 'int64',
                          'numpy_type': 'int64',
                          'metadata': None}],
             'pandas_version': '0.23.4'}
        )})
    r3 = t3.to_pandas()
    tm.assert_frame_equal(r3, e3)

    # Case 4: MultiIndex with named RangeIndex
    e4 = pd.DataFrame({
        'a': a_values
    }, index=[pd.RangeIndex(0, 8, step=2, name='qux'), b_values])
    t4 = pa.Table.from_arrays([a_arrow, rng_index_arrow, b_arrow],
                              names=['a', 'qux', gen_name_1])
    t4 = t4.replace_schema_metadata({
        b'pandas': json.dumps(
            {'index_columns': ['qux', gen_name_1],
             'column_indexes': [{'name': None,
                                 'field_name': None,
                                 'pandas_type': 'unicode',
                                 'numpy_type': 'object',
                                 'metadata': {'encoding': 'UTF-8'}}],
             'columns': [{'name': 'a',
                          'field_name': 'a',
                          'pandas_type': 'unicode',
                          'numpy_type': 'object',
                          'metadata': None},
                         {'name': 'qux',
                          'field_name': 'qux',
                          'pandas_type': 'int64',
                          'numpy_type': 'int64',
                          'metadata': None},
                         {'name': None,
                          'field_name': gen_name_1,
                          'pandas_type': 'unicode',
                          'numpy_type': 'object',
                          'metadata': None}],
             'pandas_version': '0.23.4'}
        )})
    r4 = t4.to_pandas()
    tm.assert_frame_equal(r4, e4)

    # Case 4: MultiIndex with unnamed RangeIndex
    e5 = pd.DataFrame({
        'a': a_values
    }, index=[pd.RangeIndex(0, 8, step=2, name=None), b_values])
    t5 = pa.Table.from_arrays([a_arrow, rng_index_arrow, b_arrow],
                              names=['a', gen_name_0, gen_name_1])
    t5 = t5.replace_schema_metadata({
        b'pandas': json.dumps(
            {'index_columns': [gen_name_0, gen_name_1],
             'column_indexes': [{'name': None,
                                 'field_name': None,
                                 'pandas_type': 'unicode',
                                 'numpy_type': 'object',
                                 'metadata': {'encoding': 'UTF-8'}}],
             'columns': [{'name': 'a',
                          'field_name': 'a',
                          'pandas_type': 'unicode',
                          'numpy_type': 'object',
                          'metadata': None},
                         {'name': None,
                          'field_name': gen_name_0,
                          'pandas_type': 'int64',
                          'numpy_type': 'int64',
                          'metadata': None},
                         {'name': None,
                          'field_name': gen_name_1,
                          'pandas_type': 'unicode',
                          'numpy_type': 'object',
                          'metadata': None}],
             'pandas_version': '0.23.4'}
        )})
    r5 = t5.to_pandas()
    tm.assert_frame_equal(r5, e5)


def test_metadata_compat_missing_field_name():
    # Combination of missing field name but with index column as metadata.
    # This combo occurs in the latest versions of fastparquet (0.3.2), but not
    # in pyarrow itself (since field_name was added in 0.8, index as metadata
    # only added later)

    a_values = [1, 2, 3, 4]
    b_values = ['a', 'b', 'c', 'd']
    a_arrow = pa.array(a_values, type='int64')
    b_arrow = pa.array(b_values, type='utf8')

    expected = pd.DataFrame({
        'a': a_values,
        'b': b_values,
    }, index=pd.RangeIndex(0, 8, step=2, name='qux'))
    table = pa.table({'a': a_arrow, 'b': b_arrow})

    # metadata generated by fastparquet 0.3.2 with missing field_names
    table = table.replace_schema_metadata({
        b'pandas': json.dumps({
            'column_indexes': [
                {'field_name': None,
                 'metadata': None,
                 'name': None,
                 'numpy_type': 'object',
                 'pandas_type': 'mixed-integer'}
            ],
            'columns': [
                {'metadata': None,
                 'name': 'a',
                 'numpy_type': 'int64',
                 'pandas_type': 'int64'},
                {'metadata': None,
                 'name': 'b',
                 'numpy_type': 'object',
                 'pandas_type': 'unicode'}
            ],
            'index_columns': [
                {'kind': 'range',
                 'name': 'qux',
                 'start': 0,
                 'step': 2,
                 'stop': 8}
            ],
            'pandas_version': '0.25.0'}

        )})
    result = table.to_pandas()
    tm.assert_frame_equal(result, expected)


def test_metadata_index_name_not_json_serializable():
    name = np.int64(6)  # not json serializable by default
    table = pa.table(pd.DataFrame(index=pd.RangeIndex(0, 4, name=name)))
    metadata = table.schema.pandas_metadata
    assert metadata['index_columns'][0]['name'] == '6'


def test_metadata_index_name_is_json_serializable():
    name = 6  # json serializable by default
    table = pa.table(pd.DataFrame(index=pd.RangeIndex(0, 4, name=name)))
    metadata = table.schema.pandas_metadata
    assert metadata['index_columns'][0]['name'] == 6


def make_df_with_timestamps():
    # Some of the milliseconds timestamps deliberately don't fit in the range
    # that is possible with nanosecond timestamps.
    df = pd.DataFrame({
        'dateTimeMs': [
            np.datetime64('0001-01-01 00:00', 'ms'),
            np.datetime64('2012-05-02 12:35', 'ms'),
            np.datetime64('2012-05-03 15:42', 'ms'),
            np.datetime64('3000-05-03 15:42', 'ms'),
        ],
        'dateTimeNs': [
            np.datetime64('1991-01-01 00:00', 'ns'),
            np.datetime64('2012-05-02 12:35', 'ns'),
            np.datetime64('2012-05-03 15:42', 'ns'),
            np.datetime64('2050-05-03 15:42', 'ns'),
        ],
    })
    df['dateTimeMs'] = df['dateTimeMs'].astype('object')
    # Not part of what we're testing, just ensuring that the inputs are what we
    # expect.
    assert (df.dateTimeMs.dtype, df.dateTimeNs.dtype) == (
        # O == object, M8[ns] == timestamp64[ns]
        np.dtype("O"), np.dtype("M8[ns]")
    )
    return df


@pytest.mark.parquet
def test_timestamp_as_object_parquet(tempdir):
    # Timestamps can be stored as Parquet and reloaded into Pandas with no loss
    # of information if the timestamp_as_object option is True.
    df = make_df_with_timestamps()
    table = pa.Table.from_pandas(df)
    filename = tempdir / "timestamps_from_pandas.parquet"
    pq.write_table(table, filename)
    result = pq.read_table(filename)
    df2 = result.to_pandas(timestamp_as_object=True)
    tm.assert_frame_equal(df, df2)


def test_timestamp_as_object_out_of_range():
    # Out of range timestamps can be converted Arrow and reloaded into Pandas
    # with no loss of information if the timestamp_as_object option is True.
    df = make_df_with_timestamps()
    table = pa.Table.from_pandas(df)
    df2 = table.to_pandas(timestamp_as_object=True)
    tm.assert_frame_equal(df, df2)


@pytest.mark.parametrize("resolution", ["s", "ms", "us"])
@pytest.mark.parametrize("tz", [None, "America/New_York"])
# One datetime outside nanosecond range, one inside nanosecond range:
@pytest.mark.parametrize("dt", [datetime(1553, 1, 1), datetime(2020, 1, 1)])
def test_timestamp_as_object_non_nanosecond(resolution, tz, dt):
    # Timestamps can be converted Arrow and reloaded into Pandas with no loss
    # of information if the timestamp_as_object option is True.
    arr = pa.array([dt], type=pa.timestamp(resolution, tz=tz))
    table = pa.table({'a': arr})

    for result in [
        arr.to_pandas(timestamp_as_object=True),
        table.to_pandas(timestamp_as_object=True)['a']
    ]:
        assert result.dtype == object
        assert isinstance(result[0], datetime)
        if tz:
            assert result[0].tzinfo is not None
            expected = result[0].tzinfo.fromutc(dt)
        else:
            assert result[0].tzinfo is None
            expected = dt
        assert result[0] == expected


def test_timestamp_as_object_fixed_offset():
    # ARROW-16547 to_pandas with timestamp_as_object=True and FixedOffset
    pytz = pytest.importorskip("pytz")
    import datetime
    timezone = pytz.FixedOffset(120)
    dt = timezone.localize(datetime.datetime(2022, 5, 12, 16, 57))

    table = pa.table({"timestamp_col": pa.array([dt])})
    result = table.to_pandas(timestamp_as_object=True)
    assert pa.table(result) == table


@pytest.mark.processes
def test_threaded_pandas_import():
    invoke_script("pandas_threaded_import.py")


def test_does_not_mutate_timedelta_dtype():
    expected = np.dtype('m8')

    assert np.dtype(np.timedelta64) == expected

    df = pd.DataFrame({"a": [np.timedelta64()]})
    t = pa.Table.from_pandas(df)
    t.to_pandas()

    assert np.dtype(np.timedelta64) == expected


def test_does_not_mutate_timedelta_nested():
    # ARROW-17893: dataframe with timedelta and a list of dictionary
    # also with timedelta produces wrong result with to_pandas

    from datetime import timedelta
    timedelta_1 = [{"timedelta_1": timedelta(seconds=12, microseconds=1)}]
    timedelta_2 = [timedelta(hours=3, minutes=40, seconds=23)]
    table = pa.table({"timedelta_1": timedelta_1, "timedelta_2": timedelta_2})
    df = table.to_pandas()

    assert df["timedelta_2"][0].to_pytimedelta() == timedelta_2[0]


def test_roundtrip_nested_map_table_with_pydicts():
    schema = pa.schema([
        pa.field(
            "a",
            pa.list_(
                pa.map_(pa.int8(), pa.struct([pa.field("b", pa.binary())]))
            )
        )
    ])
    table = pa.table([[
        [[(1, None)]],
        None,
        [
            [(2, {"b": b"abc"})],
            [(3, {"b": None}), (4, {"b": b"def"})],
        ]
    ]],
        schema=schema,
    )

    expected_default_df = pd.DataFrame(
        {"a": [[[(1, None)]], None, [[(2, {"b": b"abc"})],
                                     [(3, {"b": None}), (4, {"b": b"def"})]]]}
    )
    expected_as_pydicts_df = pd.DataFrame(
        {"a": [
            [{1: None}],
            None,
            [{2: {"b": b"abc"}}, {3: {"b": None}, 4: {"b": b"def"}}],
        ]}
    )

    default_df = table.to_pandas()
    as_pydicts_df = table.to_pandas(maps_as_pydicts="strict")

    tm.assert_frame_equal(default_df, expected_default_df)
    tm.assert_frame_equal(as_pydicts_df, expected_as_pydicts_df)

    table_default_roundtrip = pa.Table.from_pandas(default_df, schema=schema)
    assert table.equals(table_default_roundtrip)

    table_as_pydicts_roundtrip = pa.Table.from_pandas(as_pydicts_df, schema=schema)
    assert table.equals(table_as_pydicts_roundtrip)


def test_roundtrip_nested_map_array_with_pydicts_sliced():
    """
    Slightly more robust test with chunking and slicing
    """
    keys_1 = pa.array(['foo', 'bar'])
    keys_2 = pa.array(['baz', 'qux', 'quux', 'quz'])
    keys_3 = pa.array([], pa.string())

    items_1 = pa.array(
        [['a', 'b'], ['c', 'd']],
        pa.list_(pa.string()),
    )
    items_2 = pa.array(
        [[], None, [None, 'e'], ['f', 'g']],
        pa.list_(pa.string()),
    )
    items_3 = pa.array(
        [],
        pa.list_(pa.string()),
    )

    map_chunk_1 = pa.MapArray.from_arrays([0, 2], keys_1, items_1)
    map_chunk_2 = pa.MapArray.from_arrays([0, 3, 4], keys_2, items_2)
    map_chunk_3 = pa.MapArray.from_arrays([0, 0], keys_3, items_3)
    chunked_array = pa.chunked_array([
        pa.ListArray.from_arrays([0, 1], map_chunk_1).slice(0),
        pa.ListArray.from_arrays([0, 1], map_chunk_2.slice(1)).slice(0),
        pa.ListArray.from_arrays([0, 0], map_chunk_3).slice(0),
    ])

    series_default = chunked_array.to_pandas()
    expected_series_default = pd.Series([
        [[('foo', ['a', 'b']), ('bar', ['c', 'd'])]],
        [[('quz', ['f', 'g'])]],
        [],
    ])

    series_pydicts = chunked_array.to_pandas(maps_as_pydicts="strict")
    expected_series_pydicts = pd.Series([
        [{'foo': ['a', 'b'], 'bar': ['c', 'd']}],
        [{'quz': ['f', 'g']}],
        [],
    ])

    sliced = chunked_array.slice(1, 3)
    series_default_sliced = sliced.to_pandas()
    expected_series_default_sliced = pd.Series([
        [[('quz', ['f', 'g'])]],
        [],
    ])

    series_pydicts_sliced = sliced.to_pandas(maps_as_pydicts="strict")
    expected_series_pydicts_sliced = pd.Series([
        [{'quz': ['f', 'g']}],
        [],
    ])

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "elementwise comparison failed",
                                DeprecationWarning)
        tm.assert_series_equal(series_default, expected_series_default)
        tm.assert_series_equal(series_pydicts, expected_series_pydicts)
        tm.assert_series_equal(series_default_sliced, expected_series_default_sliced)
        tm.assert_series_equal(series_pydicts_sliced, expected_series_pydicts_sliced)

    ty = pa.list_(pa.map_(pa.string(), pa.list_(pa.string())))

    def assert_roundtrip(series: pd.Series, data) -> None:
        array_roundtrip = pa.chunked_array(pa.Array.from_pandas(series, type=ty))
        array_roundtrip.validate(full=True)
        assert data.equals(array_roundtrip)

    assert_roundtrip(series_default, chunked_array)
    assert_roundtrip(series_pydicts, chunked_array)
    assert_roundtrip(series_default_sliced, sliced)
    assert_roundtrip(series_pydicts_sliced, sliced)


def test_roundtrip_map_array_with_pydicts_duplicate_keys():
    keys = pa.array(['foo', 'bar', 'foo'])
    items = pa.array(
        [['a', 'b'], ['c', 'd'], ['1', '2']],
        pa.list_(pa.string()),
    )
    offsets = [0, 3]
    maps = pa.MapArray.from_arrays(offsets, keys, items)
    ty = pa.map_(pa.string(), pa.list_(pa.string()))

    # ------------------------
    # With maps as pydicts
    with pytest.raises(pa.lib.ArrowException):
        # raises because of duplicate keys
        maps.to_pandas(maps_as_pydicts="strict")
    series_pydicts = maps.to_pandas(maps_as_pydicts="lossy")
    # some data loss occurs for duplicate keys
    expected_series_pydicts = pd.Series([
        {'foo': ['1', '2'], 'bar': ['c', 'd']},
    ])
    # roundtrip is not possible because of data loss
    assert not maps.equals(pa.Array.from_pandas(series_pydicts, type=ty))

    # ------------------------
    # With default assoc list of tuples
    series_default = maps.to_pandas()
    expected_series_default = pd.Series([
        [('foo', ['a', 'b']), ('bar', ['c', 'd']), ('foo', ['1', '2'])],
    ])
    assert maps.equals(pa.Array.from_pandas(series_default, type=ty))

    # custom comparison for compatibility w/ Pandas 1.0.0
    # would otherwise run:
    #   tm.assert_series_equal(series_pydicts, expected_series_pydicts)
    assert len(series_pydicts) == len(expected_series_pydicts)
    for row1, row2 in zip(series_pydicts, expected_series_pydicts):
        assert len(row1) == len(row2)
        for tup1, tup2 in zip(row1.items(), row2.items()):
            assert tup1[0] == tup2[0]
            assert np.array_equal(tup1[1], tup2[1])

    # custom comparison for compatibility w/ Pandas 1.0.0
    # would otherwise run:
    #   tm.assert_series_equal(series_default, expected_series_default)
    assert len(series_default) == len(expected_series_default)
    for row1, row2 in zip(series_default, expected_series_default):
        assert len(row1) == len(row2)
        for tup1, tup2 in zip(row1, row2):
            assert tup1[0] == tup2[0]
            assert np.array_equal(tup1[1], tup2[1])


def test_unhashable_map_keys_with_pydicts():
    keys = pa.array(
        [['a', 'b'], ['c', 'd'], [], ['e'], [None, 'f'], ['g', 'h']],
        pa.list_(pa.string()),
    )
    items = pa.array(['foo', 'bar', 'baz', 'qux', 'quux', 'quz'])
    offsets = [0, 2, 6]
    maps = pa.MapArray.from_arrays(offsets, keys, items)

    # ------------------------
    # With maps as pydicts
    with pytest.raises(TypeError):
        maps.to_pandas(maps_as_pydicts="lossy")

    # ------------------------
    # With default assoc list of tuples
    series = maps.to_pandas()
    expected_series_default = pd.Series([
        [(['a', 'b'], 'foo'), (['c', 'd'], 'bar')],
        [([], 'baz'), (['e'], 'qux'), ([None, 'f'], 'quux'), (['g', 'h'], 'quz')],
    ])

    # custom comparison for compatibility w/ Pandas 1.0.0
    # would otherwise run:
    #   tm.assert_series_equal(series, expected_series_default)
    assert len(series) == len(expected_series_default)
    for row1, row2 in zip(series, expected_series_default):
        assert len(row1) == len(row2)
        for tup1, tup2 in zip(row1, row2):
            assert np.array_equal(tup1[0], tup2[0])
            assert tup1[1] == tup2[1]


def test_table_column_conversion_for_datetime():
    # GH-35235
    # pandas implemented __from_arrow__ for DatetimeTZDtype,
    # but we choose to do the conversion in Arrow instead.
    # https://github.com/pandas-dev/pandas/pull/52201
    series = pd.Series(pd.date_range("2012", periods=2, tz="Europe/Brussels"),
                       name="datetime_column")
    table = pa.table({"datetime_column": pa.array(series)})
    table_col = table.column("datetime_column")

    result = table_col.to_pandas()
    assert result.name == "datetime_column"
    tm.assert_series_equal(result, series)


def test_array_conversion_for_datetime():
    # GH-35235
    # pandas implemented __from_arrow__ for DatetimeTZDtype,
    # but we choose to do the conversion in Arrow instead.
    # https://github.com/pandas-dev/pandas/pull/52201
    series = pd.Series(pd.date_range("2012", periods=2, tz="Europe/Brussels"))
    arr = pa.array(series)

    result = arr.to_pandas()
    tm.assert_series_equal(result, series)


@pytest.mark.large_memory
def test_nested_chunking_valid():
    # GH-32439: Chunking can cause arrays to be in invalid state
    # when nested types are involved.
    # Here we simply ensure we validate correctly.

    def roundtrip(df, schema=None):
        tab = pa.Table.from_pandas(df, schema=schema)
        tab.validate(full=True)
        # we expect to trigger chunking internally
        # an assertion failure here may just mean this threshold has changed
        num_chunks = tab.column(0).num_chunks
        assert num_chunks > 1
        tm.assert_frame_equal(tab.to_pandas(self_destruct=True,
                                            maps_as_pydicts="strict"), df)

    x = b"0" * 720000000
    roundtrip(pd.DataFrame({"strings": [x, x, x]}))

    struct = {"struct_field": x}
    roundtrip(pd.DataFrame({"structs": [struct, struct, struct]}))

    lists = [x]
    roundtrip(pd.DataFrame({"lists": [lists, lists, lists]}))

    los = [struct]
    roundtrip(pd.DataFrame({"los": [los, los, los]}))

    sol = {"struct_field": lists}
    roundtrip(pd.DataFrame({"sol": [sol, sol, sol]}))

    map_of_los = {"a": los}
    map_type = pa.map_(pa.string(),
                       pa.list_(pa.struct([("struct_field", pa.binary())])))
    schema = pa.schema([("maps", map_type)])
    roundtrip(pd.DataFrame({"maps": [map_of_los, map_of_los, map_of_los]}),
              schema=schema)


def test_bytes_column_name_to_pandas():
    df = pd.DataFrame([[0.1, 0.2], [0.3, 0.4]], columns=[b'col1', b'col2'])
    table = pa.Table.from_pandas(df)
    assert table.column_names == ['col1', 'col2']
    assert table.to_pandas().equals(df)


@pytest.mark.processes
def test_is_data_frame_race_condition():
    # See https://github.com/apache/arrow/issues/39313
    test_util.invoke_script('arrow_39313.py')


def test_json_unserializable_pd_df_attrs():
    df = pd.DataFrame({"x": [1, 2, 3]})

    df.attrs["timestamp"] = datetime.fromisoformat("2025-10-28T14:20:42")

    with pytest.warns(
        UserWarning,
        match="Could not serialize pd.DataFrame.attrs:",
    ):
        df_table = pa.table(df)

    pd_metadata = json.loads(df_table.schema.metadata[b"pandas"])

    assert not pd_metadata["attributes"]
