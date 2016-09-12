# -*- coding: utf-8 -*-

from __future__ import print_function

from numpy import nan
import numpy as np

from pandas import compat
from pandas import (DataFrame, Series, MultiIndex, Timestamp,
                    date_range)

import pandas.util.testing as tm

from pandas.tests.frame.common import TestData


class TestDataFrameConvertTo(tm.TestCase, TestData):

    _multiprocess_can_split_ = True

    def test_to_dict(self):
        test_data = {
            'A': {'1': 1, '2': 2},
            'B': {'1': '1', '2': '2', '3': '3'},
        }
        recons_data = DataFrame(test_data).to_dict()

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                self.assertEqual(v2, recons_data[k][k2])

        recons_data = DataFrame(test_data).to_dict("l")

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                self.assertEqual(v2, recons_data[k][int(k2) - 1])

        recons_data = DataFrame(test_data).to_dict("s")

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                self.assertEqual(v2, recons_data[k][k2])

        recons_data = DataFrame(test_data).to_dict("sp")
        expected_split = {'columns': ['A', 'B'], 'index': ['1', '2', '3'],
                          'data': [[1.0, '1'], [2.0, '2'], [nan, '3']]}
        tm.assert_dict_equal(recons_data, expected_split)

        recons_data = DataFrame(test_data).to_dict("r")
        expected_records = [{'A': 1.0, 'B': '1'},
                            {'A': 2.0, 'B': '2'},
                            {'A': nan, 'B': '3'}]
        tm.assertIsInstance(recons_data, list)
        self.assertEqual(len(recons_data), 3)
        for l, r in zip(recons_data, expected_records):
            tm.assert_dict_equal(l, r)

        # GH10844
        recons_data = DataFrame(test_data).to_dict("i")

        for k, v in compat.iteritems(test_data):
            for k2, v2 in compat.iteritems(v):
                self.assertEqual(v2, recons_data[k2][k])

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

        self.assertEqual(test_data.to_dict(orient='records'),
                         expected_records)
        self.assertEqual(test_data_mixed.to_dict(orient='records'),
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

    def test_to_dict_invalid_orient(self):
        df = DataFrame({'A': [0, 1]})
        self.assertRaises(ValueError, df.to_dict, orient='xinvalid')

    def test_to_records_dt64(self):
        df = DataFrame([["one", "two", "three"],
                        ["four", "five", "six"]],
                       index=date_range("2012-01-01", "2012-01-02"))
        self.assertEqual(df.to_records()['index'][0], df.index[0])

        rs = df.to_records(convert_datetime64=False)
        self.assertEqual(rs['index'][0], df.index.values[0])

    def test_to_records_with_multindex(self):
        # GH3189
        index = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
                 ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
        data = np.zeros((8, 4))
        df = DataFrame(data, index=index)
        r = df.to_records(index=True)['level_0']
        self.assertTrue('bar' in r)
        self.assertTrue('one' not in r)

    def test_to_records_with_Mapping_type(self):
        import email
        from email.parser import Parser
        import collections

        collections.Mapping.register(email.message.Message)

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
        self.assertIn('X', rs.dtype.fields)

        df = DataFrame(np.random.randn(3, 3))
        rs = df.to_records()
        self.assertIn('index', rs.dtype.fields)

        df.index = MultiIndex.from_tuples([('a', 'x'), ('a', 'y'), ('b', 'z')])
        df.index.names = ['A', None]
        rs = df.to_records()
        self.assertIn('level_0', rs.dtype.fields)

    def test_to_records_with_unicode_index(self):
        # GH13172
        # unicode_literals conflict with to_records
        result = DataFrame([{u'a': u'x', u'b': 'y'}]).set_index(u'a')\
            .to_records()
        expected = np.rec.array([('x', 'y')], dtype=[('a', 'O'), ('b', 'O')])
        tm.assert_almost_equal(result, expected)
