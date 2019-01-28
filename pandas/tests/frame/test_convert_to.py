# -*- coding: utf-8 -*-

import collections
from collections import OrderedDict, defaultdict
from datetime import datetime

import numpy as np
import pytest
import pytz

from pandas.compat import long

from pandas import DataFrame, MultiIndex, Series, Timestamp, compat, date_range
from pandas.tests.frame.common import TestData
import pandas.util.testing as tm


class TestDataFrameConvertTo(TestData):

    def test_to_dict_timestamp(self):

        # GH11247
        # split/records producing np.datetime64 rather than Timestamps
        # on datetime64[ns] dtypes only

        tsmp = Timestamp('20130101')
        test_data = DataFrame({'A': [tsmp, tsmp], 'B': [tsmp, tsmp]})
        test_data_mixed = DataFrame({'A': [tsmp, tsmp], 'B': [1, 2]})

        expected_records = [{'A': tsmp, 'B': tsmp},
                            {'A': tsmp, 'B': tsmp}]
        expected_records_mixed = [{'A': tsmp, 'B': 1},
                                  {'A': tsmp, 'B': 2}]

        assert (test_data.to_dict(orient='records') ==
                expected_records)
        assert (test_data_mixed.to_dict(orient='records') ==
                expected_records_mixed)

        expected_series = {
            'A': Series([tsmp, tsmp], name='A'),
            'B': Series([tsmp, tsmp], name='B'),
        }
        expected_series_mixed = {
            'A': Series([tsmp, tsmp], name='A'),
            'B': Series([1, 2], name='B'),
        }

        tm.assert_dict_equal(test_data.to_dict(orient='series'),
                             expected_series)
        tm.assert_dict_equal(test_data_mixed.to_dict(orient='series'),
                             expected_series_mixed)

        expected_split = {
            'index': [0, 1],
            'data': [[tsmp, tsmp],
                     [tsmp, tsmp]],
            'columns': ['A', 'B']
        }
        expected_split_mixed = {
            'index': [0, 1],
            'data': [[tsmp, 1],
                     [tsmp, 2]],
            'columns': ['A', 'B']
        }

        tm.assert_dict_equal(test_data.to_dict(orient='split'),
                             expected_split)
        tm.assert_dict_equal(test_data_mixed.to_dict(orient='split'),
                             expected_split_mixed)

    def test_to_dict_index_not_unique_with_index_orient(self):
        # GH22801
        # Data loss when indexes are not unique. Raise ValueError.
        df = DataFrame({'a': [1, 2], 'b': [0.5, 0.75]}, index=['A', 'A'])
        pytest.raises(ValueError, df.to_dict, orient='index')

    def test_to_dict_invalid_orient(self):
        df = DataFrame({'A': [0, 1]})
        pytest.raises(ValueError, df.to_dict, orient='xinvalid')

    def test_to_records_dt64(self):
        df = DataFrame([["one", "two", "three"],
                        ["four", "five", "six"]],
                       index=date_range("2012-01-01", "2012-01-02"))

        # convert_datetime64 defaults to None
        expected = df.index.values[0]
        result = df.to_records()['index'][0]
        assert expected == result

        # check for FutureWarning if convert_datetime64=False is passed
        with tm.assert_produces_warning(FutureWarning):
            expected = df.index.values[0]
            result = df.to_records(convert_datetime64=False)['index'][0]
            assert expected == result

        # check for FutureWarning if convert_datetime64=True is passed
        with tm.assert_produces_warning(FutureWarning):
            expected = df.index[0]
            result = df.to_records(convert_datetime64=True)['index'][0]
            assert expected == result

    def test_to_records_with_multindex(self):
        # GH3189
        index = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                 ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        data = np.zeros((8, 4))
        df = DataFrame(data, index=index)
        r = df.to_records(index=True)['level_0']
        assert 'bar' in r
        assert 'one' not in r

    def test_to_records_with_Mapping_type(self):
        import email
        from email.parser import Parser

        compat.Mapping.register(email.message.Message)

        headers = Parser().parsestr('From: <user@example.com>\n'
                                    'To: <someone_else@example.com>\n'
                                    'Subject: Test message\n'
                                    '\n'
                                    'Body would go here\n')

        frame = DataFrame.from_records([headers])
        all(x in frame for x in ['Type', 'Subject', 'From'])

    def test_to_records_floats(self):
        df = DataFrame(np.random.rand(10, 10))
        df.to_records()

    def test_to_records_index_name(self):
        df = DataFrame(np.random.randn(3, 3))
        df.index.name = 'X'
        rs = df.to_records()
        assert 'X' in rs.dtype.fields

        df = DataFrame(np.random.randn(3, 3))
        rs = df.to_records()
        assert 'index' in rs.dtype.fields

        df.index = MultiIndex.from_tuples([('a', 'x'), ('a', 'y'), ('b', 'z')])
        df.index.names = ['A', None]
        rs = df.to_records()
        assert 'level_0' in rs.dtype.fields

    def test_to_records_with_unicode_index(self):
        # GH13172
        # unicode_literals conflict with to_records
        result = DataFrame([{u'a': u'x', u'b': 'y'}]).set_index(u'a') \
            .to_records()
        expected = np.rec.array([('x', 'y')], dtype=[('a', 'O'), ('b', 'O')])
        tm.assert_almost_equal(result, expected)

    def test_to_records_with_unicode_column_names(self):
        # xref issue: https://github.com/numpy/numpy/issues/2407
        # Issue #11879. to_records used to raise an exception when used
        # with column names containing non-ascii characters in Python 2
        result = DataFrame(data={u"accented_name_é": [1.0]}).to_records()

        # Note that numpy allows for unicode field names but dtypes need
        # to be specified using dictionary instead of list of tuples.
        expected = np.rec.array(
            [(0, 1.0)],
            dtype={"names": ["index", u"accented_name_é"],
                   "formats": ['=i8', '=f8']}
        )
        tm.assert_almost_equal(result, expected)

    def test_to_records_with_categorical(self):

        # GH8626

        # dict creation
        df = DataFrame({'A': list('abc')}, dtype='category')
        expected = Series(list('abc'), dtype='category', name='A')
        tm.assert_series_equal(df['A'], expected)

        # list-like creation
        df = DataFrame(list('abc'), dtype='category')
        expected = Series(list('abc'), dtype='category', name=0)
        tm.assert_series_equal(df[0], expected)

        # to record array
        # this coerces
        result = df.to_records()
        expected = np.rec.array([(0, 'a'), (1, 'b'), (2, 'c')],
                                dtype=[('index', '=i8'), ('0', 'O')])
        tm.assert_almost_equal(result, expected)

    @pytest.mark.parametrize("kwargs,expected", [
        # No dtypes --> default to array dtypes.
        (dict(),
         np.rec.array([(0, 1, 0.2, "a"), (1, 2, 1.5, "bc")],
                      dtype=[("index", "<i8"), ("A", "<i8"),
                             ("B", "<f8"), ("C", "O")])),

        # Should have no effect in this case.
        (dict(index=True),
         np.rec.array([(0, 1, 0.2, "a"), (1, 2, 1.5, "bc")],
                      dtype=[("index", "<i8"), ("A", "<i8"),
                             ("B", "<f8"), ("C", "O")])),

        # Column dtype applied across the board. Index unaffected.
        (dict(column_dtypes="<U4"),
         np.rec.array([("0", "1", "0.2", "a"), ("1", "2", "1.5", "bc")],
                      dtype=[("index", "<i8"), ("A", "<U4"),
                             ("B", "<U4"), ("C", "<U4")])),

        # Index dtype applied across the board. Columns unaffected.
        (dict(index_dtypes="<U1"),
         np.rec.array([("0", 1, 0.2, "a"), ("1", 2, 1.5, "bc")],
                      dtype=[("index", "<U1"), ("A", "<i8"),
                             ("B", "<f8"), ("C", "O")])),

        # Pass in a type instance.
        (dict(column_dtypes=np.unicode),
         np.rec.array([("0", "1", "0.2", "a"), ("1", "2", "1.5", "bc")],
                      dtype=[("index", "<i8"), ("A", "<U"),
                             ("B", "<U"), ("C", "<U")])),

        # Pass in a dictionary (name-only).
        (dict(column_dtypes={"A": np.int8, "B": np.float32, "C": "<U2"}),
         np.rec.array([("0", "1", "0.2", "a"), ("1", "2", "1.5", "bc")],
                      dtype=[("index", "<i8"), ("A", "i1"),
                             ("B", "<f4"), ("C", "<U2")])),

        # Pass in a dictionary (indices-only).
        (dict(index_dtypes={0: "int16"}),
         np.rec.array([(0, 1, 0.2, "a"), (1, 2, 1.5, "bc")],
                      dtype=[("index", "i2"), ("A", "<i8"),
                             ("B", "<f8"), ("C", "O")])),

        # Ignore index mappings if index is not True.
        (dict(index=False, index_dtypes="<U2"),
         np.rec.array([(1, 0.2, "a"), (2, 1.5, "bc")],
                      dtype=[("A", "<i8"), ("B", "<f8"), ("C", "O")])),

        # Non-existent names / indices in mapping should not error.
        (dict(index_dtypes={0: "int16", "not-there": "float32"}),
         np.rec.array([(0, 1, 0.2, "a"), (1, 2, 1.5, "bc")],
                      dtype=[("index", "i2"), ("A", "<i8"),
                             ("B", "<f8"), ("C", "O")])),

        # Names / indices not in mapping default to array dtype.
        (dict(column_dtypes={"A": np.int8, "B": np.float32}),
         np.rec.array([("0", "1", "0.2", "a"), ("1", "2", "1.5", "bc")],
                      dtype=[("index", "<i8"), ("A", "i1"),
                             ("B", "<f4"), ("C", "O")])),

        # Mixture of everything.
        (dict(column_dtypes={"A": np.int8, "B": np.float32},
              index_dtypes="<U2"),
         np.rec.array([("0", "1", "0.2", "a"), ("1", "2", "1.5", "bc")],
                      dtype=[("index", "<U2"), ("A", "i1"),
                             ("B", "<f4"), ("C", "O")])),

        # Invalid dype values.
        (dict(index=False, column_dtypes=list()),
         "Invalid dtype \\[\\] specified for column A"),

        (dict(index=False, column_dtypes={"A": "int32", "B": 5}),
         "Invalid dtype 5 specified for column B"),
    ])
    def test_to_records_dtype(self, kwargs, expected):
        # see gh-18146
        df = DataFrame({"A": [1, 2], "B": [0.2, 1.5], "C": ["a", "bc"]})

        if isinstance(expected, str):
            with pytest.raises(ValueError, match=expected):
                df.to_records(**kwargs)
        else:
            result = df.to_records(**kwargs)
            tm.assert_almost_equal(result, expected)

    @pytest.mark.parametrize("df,kwargs,expected", [
        # MultiIndex in the index.
        (DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                   columns=list("abc")).set_index(["a", "b"]),
         dict(column_dtypes="float64", index_dtypes={0: "int32", 1: "int8"}),
         np.rec.array([(1, 2, 3.), (4, 5, 6.), (7, 8, 9.)],
                      dtype=[("a", "<i4"), ("b", "i1"), ("c", "<f8")])),

        # MultiIndex in the columns.
        (DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                   columns=MultiIndex.from_tuples([("a", "d"), ("b", "e"),
                                                   ("c", "f")])),
         dict(column_dtypes={0: "<U1", 2: "float32"}, index_dtypes="float32"),
         np.rec.array([(0., u"1", 2, 3.), (1., u"4", 5, 6.),
                       (2., u"7", 8, 9.)],
                      dtype=[("index", "<f4"),
                             ("('a', 'd')", "<U1"),
                             ("('b', 'e')", "<i8"),
                             ("('c', 'f')", "<f4")])),

        # MultiIndex in both the columns and index.
        (DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                   columns=MultiIndex.from_tuples([
                       ("a", "d"), ("b", "e"), ("c", "f")], names=list("ab")),
                   index=MultiIndex.from_tuples([
                       ("d", -4), ("d", -5), ("f", -6)], names=list("cd"))),
         dict(column_dtypes="float64", index_dtypes={0: "<U2", 1: "int8"}),
         np.rec.array([("d", -4, 1., 2., 3.), ("d", -5, 4., 5., 6.),
                       ("f", -6, 7, 8, 9.)],
                      dtype=[("c", "<U2"), ("d", "i1"),
                             ("('a', 'd')", "<f8"), ("('b', 'e')", "<f8"),
                             ("('c', 'f')", "<f8")]))
    ])
    def test_to_records_dtype_mi(self, df, kwargs, expected):
        # see gh-18146
        result = df.to_records(**kwargs)
        tm.assert_almost_equal(result, expected)

    def test_to_records_dict_like(self):
        # see gh-18146
        class DictLike(object):
            def __init__(self, **kwargs):
                self.d = kwargs.copy()

            def __getitem__(self, key):
                return self.d.__getitem__(key)

            def __contains__(self, key):
                return key in self.d

            def keys(self):
                return self.d.keys()

        df = DataFrame({"A": [1, 2], "B": [0.2, 1.5], "C": ["a", "bc"]})

        dtype_mappings = dict(column_dtypes=DictLike(**{"A": np.int8,
                                                        "B": np.float32}),
                              index_dtypes="<U2")

        result = df.to_records(**dtype_mappings)
        expected = np.rec.array([("0", "1", "0.2", "a"),
                                 ("1", "2", "1.5", "bc")],
                                dtype=[("index", "<U2"), ("A", "i1"),
                                       ("B", "<f4"), ("C", "O")])
        tm.assert_almost_equal(result, expected)

    @pytest.mark.parametrize('mapping', [
        dict,
        collections.defaultdict(list),
        collections.OrderedDict])
    def test_to_dict(self, mapping):
        test_data = {
            'A': {'1': 1, '2': 2},
            'B': {'1': '1', '2': '2', '3': '3'},
        }

        # GH16122
        recons_data = DataFrame(test_data).to_dict(into=mapping)

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                assert (v2 == recons_data[k][k2])

        recons_data = DataFrame(test_data).to_dict("l", mapping)

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                assert (v2 == recons_data[k][int(k2) - 1])

        recons_data = DataFrame(test_data).to_dict("s", mapping)

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                assert (v2 == recons_data[k][k2])

        recons_data = DataFrame(test_data).to_dict("sp", mapping)
        expected_split = {'columns': ['A', 'B'], 'index': ['1', '2', '3'],
                          'data': [[1.0, '1'], [2.0, '2'], [np.nan, '3']]}
        tm.assert_dict_equal(recons_data, expected_split)

        recons_data = DataFrame(test_data).to_dict("r", mapping)
        expected_records = [{'A': 1.0, 'B': '1'},
                            {'A': 2.0, 'B': '2'},
                            {'A': np.nan, 'B': '3'}]
        assert isinstance(recons_data, list)
        assert (len(recons_data) == 3)
        for l, r in zip(recons_data, expected_records):
            tm.assert_dict_equal(l, r)

        # GH10844
        recons_data = DataFrame(test_data).to_dict("i")

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                assert (v2 == recons_data[k2][k])

        df = DataFrame(test_data)
        df['duped'] = df[df.columns[0]]
        recons_data = df.to_dict("i")
        comp_data = test_data.copy()
        comp_data['duped'] = comp_data[df.columns[0]]
        for k, v in compat.iteritems(comp_data):
            for k2, v2 in compat.iteritems(v):
                assert (v2 == recons_data[k2][k])

    @pytest.mark.parametrize('mapping', [
        list,
        collections.defaultdict,
        []])
    def test_to_dict_errors(self, mapping):
        # GH16122
        df = DataFrame(np.random.randn(3, 3))
        with pytest.raises(TypeError):
            df.to_dict(into=mapping)

    def test_to_dict_not_unique_warning(self):
        # GH16927: When converting to a dict, if a column has a non-unique name
        # it will be dropped, throwing a warning.
        df = DataFrame([[1, 2, 3]], columns=['a', 'a', 'b'])
        with tm.assert_produces_warning(UserWarning):
            df.to_dict()

    @pytest.mark.parametrize('tz', ['UTC', 'GMT', 'US/Eastern'])
    def test_to_records_datetimeindex_with_tz(self, tz):
        # GH13937
        dr = date_range('2016-01-01', periods=10,
                        freq='S', tz=tz)

        df = DataFrame({'datetime': dr}, index=dr)

        expected = df.to_records()
        result = df.tz_convert("UTC").to_records()

        # both converted to UTC, so they are equal
        tm.assert_numpy_array_equal(result, expected)

    # orient - orient argument to to_dict function
    # item_getter - function for extracting value from
    # the resulting dict using column name and index
    @pytest.mark.parametrize('orient,item_getter', [
        ('dict', lambda d, col, idx: d[col][idx]),
        ('records', lambda d, col, idx: d[idx][col]),
        ('list', lambda d, col, idx: d[col][idx]),
        ('split', lambda d, col, idx: d['data'][idx][d['columns'].index(col)]),
        ('index', lambda d, col, idx: d[idx][col])
    ])
    def test_to_dict_box_scalars(self, orient, item_getter):
        # 14216, 23753
        # make sure that we are boxing properly
        df = DataFrame({'a': [1, 2], 'b': [.1, .2]})
        result = df.to_dict(orient=orient)
        assert isinstance(item_getter(result, 'a', 0), (int, long))
        assert isinstance(item_getter(result, 'b', 0), float)

    def test_frame_to_dict_tz(self):
        # GH18372 When converting to dict with orient='records' columns of
        # datetime that are tz-aware were not converted to required arrays
        data = [(datetime(2017, 11, 18, 21, 53, 0, 219225, tzinfo=pytz.utc),),
                (datetime(2017, 11, 18, 22, 6, 30, 61810, tzinfo=pytz.utc,),)]
        df = DataFrame(list(data), columns=["d", ])

        result = df.to_dict(orient='records')
        expected = [
            {'d': Timestamp('2017-11-18 21:53:00.219225+0000', tz=pytz.utc)},
            {'d': Timestamp('2017-11-18 22:06:30.061810+0000', tz=pytz.utc)},
        ]
        tm.assert_dict_equal(result[0], expected[0])
        tm.assert_dict_equal(result[1], expected[1])

    @pytest.mark.parametrize('into, expected', [
        (dict, {0: {'int_col': 1, 'float_col': 1.0},
                1: {'int_col': 2, 'float_col': 2.0},
                2: {'int_col': 3, 'float_col': 3.0}}),
        (OrderedDict, OrderedDict([(0, {'int_col': 1, 'float_col': 1.0}),
                                   (1, {'int_col': 2, 'float_col': 2.0}),
                                   (2, {'int_col': 3, 'float_col': 3.0})])),
        (defaultdict(list), defaultdict(list,
                                        {0: {'int_col': 1, 'float_col': 1.0},
                                         1: {'int_col': 2, 'float_col': 2.0},
                                         2: {'int_col': 3, 'float_col': 3.0}}))
    ])
    def test_to_dict_index_dtypes(self, into, expected):
        # GH 18580
        # When using to_dict(orient='index') on a dataframe with int
        # and float columns only the int columns were cast to float

        df = DataFrame({'int_col': [1, 2, 3],
                        'float_col': [1.0, 2.0, 3.0]})

        result = df.to_dict(orient='index', into=into)
        cols = ['int_col', 'float_col']
        result = DataFrame.from_dict(result, orient='index')[cols]
        expected = DataFrame.from_dict(expected, orient='index')[cols]
        tm.assert_frame_equal(result, expected)

    def test_to_dict_numeric_names(self):
        # https://github.com/pandas-dev/pandas/issues/24940
        df = DataFrame({str(i): [i] for i in range(5)})
        result = set(df.to_dict('records')[0].keys())
        expected = set(df.columns)
        assert result == expected

    def test_to_dict_wide(self):
        # https://github.com/pandas-dev/pandas/issues/24939
        df = DataFrame({('A_{:d}'.format(i)): [i] for i in range(256)})
        result = df.to_dict('records')[0]
        expected = {'A_{:d}'.format(i): i for i in range(256)}
        assert result == expected
