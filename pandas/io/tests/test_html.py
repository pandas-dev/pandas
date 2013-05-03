import os
import re
from cStringIO import StringIO
from unittest import TestCase

import nose

import numpy as np
from numpy.testing.decorators import slow

from pandas.io.html import read_html, import_module
from pandas import DataFrame, MultiIndex
from pandas.util.testing import assert_frame_equal, network


def _skip_if_no_parser():
    try:
        import_module('lxml')
    except ImportError:
        try:
            import_module('bs4')
        except ImportError:
            raise nose.SkipTest


DATA_PATH = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data')


def _run_read_html(*args, **kwargs):
    _skip_if_no_parser()
    return read_html(*args, **kwargs)


def isframe(x):
    return isinstance(x, DataFrame)


def assert_framelist_equal(list1, list2):
    assert len(list1) == len(list2), ('lists are not of equal size '
                                      'len(list1) == {0}, '
                                      'len(list2) == {1}'.format(len(list1),
                                                                 len(list2)))
    assert all(map(lambda x, y: isframe(x) and isframe(y), list1, list2)), \
        'not all list elements are DataFrames'
    for frame_i, frame_j in zip(list1, list2):
        assert_frame_equal(frame_i, frame_j)
        assert not frame_i.empty, 'frames are both empty'


class TestLxmlReadHtml(TestCase):
    def setUp(self):
        self.spam_data = os.path.join(DATA_PATH, 'spam.html')
        self.banklist_data = os.path.join(DATA_PATH, 'failed_banklist.html')

    def run_read_html(self, *args, **kwargs):
        kwargs['flavor'] = 'lxml'
        return _run_read_html(*args, **kwargs)

    @network
    def test_banklist_url(self):
        url = 'http://www.fdic.gov/bank/individual/failed/banklist.html'
        df1 = self.run_read_html(url, 'First Federal Bank of Florida',
                                 attrs={"id": 'table'})
        df2 = self.run_read_html(url, 'Metcalf Bank', attrs={'id': 'table'})

        assert_framelist_equal(df1, df2)

    @network
    def test_spam_url(self):
        url = ('http://ndb.nal.usda.gov/ndb/foods/show/1732?fg=&man=&'
               'lfacet=&format=&count=&max=25&offset=&sort=&qlookup=spam')
        df1 = self.run_read_html(url, '.*Water.*')
        df2 = self.run_read_html(url, 'Unit')

        assert_framelist_equal(df1, df2)

    @slow
    def test_banklist(self):
        df1 = self.run_read_html(self.banklist_data, '.*Florida.*',
                                 attrs={'id': 'table'})
        df2 = self.run_read_html(self.banklist_data, 'Metcalf Bank',
                                 attrs={'id': 'table'})

        assert_framelist_equal(df1, df2)

    @slow
    def test_banklist_header(self):
        df = self.run_read_html(self.banklist_data, 'Metcalf',
                                attrs={'id': 'table'}, header=0, skiprows=1)[0]
        self.assertFalse(df.empty)
        cols = ['Bank Name', 'City', 'State', 'CERT #',
                'Acquiring Institution', 'Closing Date', 'Updated Date']
        self.assertListEqual(df.columns.values.tolist(), cols)
        self.assertEqual(df.shape[0], 499)

    def test_spam(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*',
                                 infer_types=False)
        df2 = self.run_read_html(self.spam_data, 'Unit', infer_types=False)

        assert_framelist_equal(df1, df2)

        self.assertEqual(df1[0].ix[0, 0], 'Nutrient')

    def test_spam_no_match(self):
        dfs = self.run_read_html(self.spam_data)
        for df in dfs:
            self.assertIsInstance(df, DataFrame)

    def test_banklist_no_match(self):
        dfs = self.run_read_html(self.banklist_data, attrs={'id': 'table'})
        for df in dfs:
            self.assertIsInstance(df, DataFrame)

    def test_spam_header(self):
        df = self.run_read_html(self.spam_data, '.*Water.*', header=0)[0]
        self.assertEqual(df.columns[0], 'Nutrient')
        self.assertFalse(df.empty)

    def test_skiprows_int(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*', skiprows=1)
        df2 = self.run_read_html(self.spam_data, 'Unit', skiprows=1)

        assert_framelist_equal(df1, df2)

    def test_skiprows_xrange(self):
        df1 = [self.run_read_html(self.spam_data, '.*Water.*').pop()[2:]]
        df2 = self.run_read_html(self.spam_data, 'Unit', skiprows=xrange(2))

        assert_framelist_equal(df1, df2)

    def test_skiprows_list(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*', skiprows=[1, 2])
        df2 = self.run_read_html(self.spam_data, 'Unit', skiprows=[2, 1])

        assert_framelist_equal(df1, df2)

    def test_skiprows_set(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*',
                                 skiprows=set([1, 2]))
        df2 = self.run_read_html(self.spam_data, 'Unit', skiprows=set([2, 1]))

        assert_framelist_equal(df1, df2)

    def test_skiprows_slice(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*', skiprows=1)
        df2 = self.run_read_html(self.spam_data, 'Unit', skiprows=1)

        assert_framelist_equal(df1, df2)

    def test_skiprows_slice_short(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*',
                                 skiprows=slice(2))
        df2 = self.run_read_html(self.spam_data, 'Unit', skiprows=slice(2))

        assert_framelist_equal(df1, df2)

    def test_skiprows_slice_long(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*',
                                 skiprows=slice(2, 5))
        df2 = self.run_read_html(self.spam_data, 'Unit',
                                 skiprows=slice(4, 1, -1))

        assert_framelist_equal(df1, df2)

    def test_skiprows_ndarray(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*',
                                 skiprows=np.arange(2))
        df2 = self.run_read_html(self.spam_data, 'Unit', skiprows=np.arange(2))

        assert_framelist_equal(df1, df2)

    def test_skiprows_invalid(self):
        self.assertRaises(ValueError, self.run_read_html, self.spam_data,
                          '.*Water.*', skiprows='asdf')

    def test_index(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*', index_col=0)
        df2 = self.run_read_html(self.spam_data, 'Unit', index_col=0)
        assert_framelist_equal(df1, df2)

    def test_header(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*', header=0)
        df2 = self.run_read_html(self.spam_data, 'Unit', header=0)
        assert_framelist_equal(df1, df2)
        self.assertEqual(df1[0].columns[0], 'Nutrient')

    def test_header_and_index(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*', header=0,
                                 index_col=0)
        df2 = self.run_read_html(self.spam_data, 'Unit', header=0, index_col=0)
        assert_framelist_equal(df1, df2)

    def test_infer_types(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*', header=0,
                                 index_col=0, infer_types=False)
        df2 = self.run_read_html(self.spam_data, 'Unit', header=0, index_col=0,
                                 infer_types=False)
        assert_framelist_equal(df1, df2)

        df2 = self.run_read_html(self.spam_data, 'Unit', header=0, index_col=0,
                                 infer_types=True)

        self.assertRaises(AssertionError, assert_framelist_equal, df1, df2)

    def test_string_io(self):
        with open(self.spam_data) as f:
            data1 = StringIO(f.read())

        with open(self.spam_data) as f:
            data2 = StringIO(f.read())

        df1 = self.run_read_html(data1, '.*Water.*', infer_types=False)
        df2 = self.run_read_html(data2, 'Unit', infer_types=False)
        assert_framelist_equal(df1, df2)

    def test_string(self):
        with open(self.spam_data) as f:
            data = f.read()

        df1 = self.run_read_html(data, '.*Water.*', infer_types=False)
        df2 = self.run_read_html(data, 'Unit', infer_types=False)

        assert_framelist_equal(df1, df2)

    def test_file_like(self):
        with open(self.spam_data) as f:
            df1 = self.run_read_html(f, '.*Water.*', infer_types=False)

        with open(self.spam_data) as f:
            df2 = self.run_read_html(f, 'Unit', infer_types=False)

        assert_framelist_equal(df1, df2)

    def test_bad_url_protocol(self):
        self.assertRaises(ValueError, self.run_read_html, 'git://github.com',
                          '.*Water.*')

    @slow
    def test_file_url(self):
        url = self.banklist_data
        dfs = self.run_read_html('file://' + url, 'First',
                                 attrs={'id': 'table'})
        self.assertIsInstance(dfs, list)
        for df in dfs:
            self.assertIsInstance(df, DataFrame)

    @slow
    def test_invalid_table_attrs(self):
        url = self.banklist_data
        self.assertRaises(AssertionError, self.run_read_html, url,
                          'First Federal Bank of Florida',
                          attrs={'id': 'tasdfable'})

    def _bank_data(self, *args, **kwargs):
        return self.run_read_html(self.banklist_data, 'Metcalf',
                                  attrs={'id': 'table'}, *args, **kwargs)

    @slow
    def test_multiindex_header(self):
        df = self._bank_data(header=[0, 1])[0]
        self.assertIsInstance(df.columns, MultiIndex)

    @slow
    def test_multiindex_index(self):
        df = self._bank_data(index_col=[0, 1])[0]
        self.assertIsInstance(df.index, MultiIndex)

    @slow
    def test_multiindex_header_index(self):
        df = self._bank_data(header=[0, 1], index_col=[0, 1])[0]
        self.assertIsInstance(df.columns, MultiIndex)
        self.assertIsInstance(df.index, MultiIndex)

    @slow
    def test_multiindex_header_skiprows(self):
        df = self._bank_data(header=[0, 1], skiprows=1)[0]
        self.assertIsInstance(df.columns, MultiIndex)

    @slow
    def test_multiindex_header_index_skiprows(self):
        df = self._bank_data(header=[0, 1], index_col=[0, 1], skiprows=1)[0]
        self.assertIsInstance(df.index, MultiIndex)

    @slow
    def test_regex_idempotency(self):
        url = self.banklist_data
        dfs = self.run_read_html('file://' + url,
                                 match=re.compile(re.compile('Florida')),
                                 attrs={'id': 'table'})
        self.assertIsInstance(dfs, list)
        for df in dfs:
            self.assertIsInstance(df, DataFrame)

    def test_negative_skiprows_spam(self):
        url = self.spam_data
        self.assertRaises(AssertionError, self.run_read_html, url, 'Water',
                          skiprows=-1)

    def test_negative_skiprows_banklist(self):
        url = self.banklist_data
        self.assertRaises(AssertionError, self.run_read_html, url, 'Florida',
                          skiprows=-1)

    @slow
    def test_multiple_matches(self):
        url = self.banklist_data
        dfs = self.run_read_html(url, match=r'Florida')
        self.assertIsInstance(dfs, list)
        self.assertGreater(len(dfs), 1)
        for df in dfs:
            self.assertIsInstance(df, DataFrame)


def test_invalid_flavor():
    url = 'google.com'
    nose.tools.assert_raises(AssertionError, _run_read_html, url, 'google',
                             flavor='not a* valid**++ flaver')


class TestBs4ReadHtml(TestLxmlReadHtml):
    def run_read_html(self, *args, **kwargs):
        kwargs['flavor'] = 'bs4'
        return _run_read_html(*args, **kwargs)
