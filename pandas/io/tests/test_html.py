import os
import re
from cStringIO import StringIO
from unittest import TestCase
import numbers
from urllib2 import urlopen
from contextlib import closing
import warnings

import nose

import numpy as np
from numpy.random import rand
from numpy.testing.decorators import slow

from pandas.io.html import read_html, import_module, _parse, _LxmlFrameParser
from pandas.io.html import _BeautifulSoupHtml5LibFrameParser
from pandas.io.html import _BeautifulSoupLxmlFrameParser, _remove_whitespace
from pandas import DataFrame, MultiIndex, read_csv, Timestamp
from pandas.util.testing import (assert_frame_equal, network,
                                 get_data_path)
from numpy.testing.decorators import slow

from pandas.util.testing import makeCustomDataframe as mkdf


def _have_module(module_name):
    try:
        import_module(module_name)
        return True
    except ImportError:
        return False


def _skip_if_no(module_name):
    if not _have_module(module_name):
        raise nose.SkipTest


def _skip_if_none(module_names):
    if isinstance(module_names, basestring):
        _skip_if_no(module_names)
    else:
        if not all(_have_module(module_name) for module_name in module_names):
            raise nose.SkipTest


DATA_PATH = get_data_path()


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


def _run_read_html(parser, io, match='.+', flavor='bs4', header=None,
                   index_col=None, skiprows=None, infer_types=False,
                   attrs=None):
    if isinstance(skiprows, numbers.Integral) and skiprows < 0:
        raise AssertionError('cannot skip rows starting from the end of the '
                             'data (you passed a negative value)')
    return _parse(parser, io, match, flavor, header, index_col, skiprows,
                  infer_types, attrs)


class TestLxmlReadHtml(TestCase):
    def test_to_html_compat(self):
        df = mkdf(4, 3, data_gen_f=lambda *args: rand(), c_idx_names=False,
                  r_idx_names=False).applymap('{0:.3f}'.format)
        out = df.to_html()
        res = self.run_read_html(out, attrs={'class': 'dataframe'},
                                 index_col=0)[0]
        print df.dtypes
        print res.dtypes
        assert_frame_equal(res, df)

    def setUp(self):
        self.spam_data = os.path.join(DATA_PATH, 'spam.html')
        self.banklist_data = os.path.join(DATA_PATH, 'banklist.html')

    def run_read_html(self, *args, **kwargs):
        kwargs['flavor'] = 'lxml'
        _skip_if_no('lxml')
        parser = _LxmlFrameParser
        return _run_read_html(parser, *args, **kwargs)

    @network
    @slow
    def test_banklist_url(self):
        url = 'http://www.fdic.gov/bank/individual/failed/banklist.html'
        df1 = self.run_read_html(url, 'First Federal Bank of Florida',
                                 attrs={"id": 'table'})
        df2 = self.run_read_html(url, 'Metcalf Bank', attrs={'id': 'table'})

        assert_framelist_equal(df1, df2)

    @network
    @slow
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
        def try_remove_ws(x):
            try:
                return _remove_whitespace(x)
            except AttributeError:
                return x

        df = self.run_read_html(self.banklist_data, 'Metcalf',
                                attrs={'id': 'table'}, infer_types=False)[0]
        ground_truth = read_csv(os.path.join(DATA_PATH, 'banklist.csv'),
                                converters={'Closing Date': Timestamp,
                                            'Updated Date': Timestamp})
        self.assertNotEqual(df.shape, ground_truth.shape)
        self.assertRaises(AssertionError, assert_frame_equal, df,
                          ground_truth.applymap(try_remove_ws))

    @slow
    def test_gold_canyon(self):
        gc = 'Gold Canyon'
        with open(self.banklist_data, 'r') as f:
            raw_text = f.read()

        self.assertIn(gc, raw_text)
        df = self.run_read_html(self.banklist_data, 'Gold Canyon',
                                attrs={'id': 'table'}, infer_types=False)[0]
        self.assertNotIn(gc, df.to_string())

    def test_spam(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*',
                                 infer_types=False)
        df2 = self.run_read_html(self.spam_data, 'Unit', infer_types=False)

        assert_framelist_equal(df1, df2)
        print df1[0]

        self.assertEqual(df1[0].ix[0, 0], 'Proximates')
        self.assertEqual(df1[0].columns[0], 'Nutrient')

    def test_spam_no_match(self):
        dfs = self.run_read_html(self.spam_data)
        for df in dfs:
            self.assertIsInstance(df, DataFrame)

    def test_banklist_no_match(self):
        dfs = self.run_read_html(self.banklist_data, attrs={'id': 'table'})
        for df in dfs:
            self.assertIsInstance(df, DataFrame)

    def test_spam_header(self):
        df = self.run_read_html(self.spam_data, '.*Water.*', header=0)
        df = self.run_read_html(self.spam_data, '.*Water.*', header=1)[0]
        self.assertEqual(df.columns[0], 'Water')
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

    def test_header_and_index(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*', header=1,
                                 index_col=0)
        df2 = self.run_read_html(self.spam_data, 'Unit', header=1, index_col=0)
        assert_framelist_equal(df1, df2)

    def test_infer_types(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*', index_col=0,
                                 infer_types=False)
        df2 = self.run_read_html(self.spam_data, 'Unit', index_col=0,
                                 infer_types=False)
        assert_framelist_equal(df1, df2)

        df2 = self.run_read_html(self.spam_data, 'Unit', index_col=0,
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
        url = 'http://code.google.com/p/pythonxy/wiki/StandardPlugins'
        dfs = self.run_read_html(url, match='Python',
                                 attrs={'class': 'wikitable'})
        self.assertGreater(len(dfs), 1)

    @network
    @slow
    def test_pythonxy_plugins_table(self):
        url = 'http://code.google.com/p/pythonxy/wiki/StandardPlugins'
        dfs = self.run_read_html(url, match='Python',
                                 attrs={'class': 'wikitable'})
        zz = [df.iloc[0, 0] for df in dfs]
        self.assertListEqual(sorted(zz), sorted(['Python', 'SciTE']))


def test_invalid_flavor():
    url = 'google.com'
    nose.tools.assert_raises(AssertionError, read_html, url, 'google',
                             flavor='not a* valid**++ flaver')


@slow
class TestBs4LxmlParser(TestLxmlReadHtml):
    def test(self):
        pass

    def run_read_html(self, *args, **kwargs):
        kwargs['flavor'] = 'bs4'
        _skip_if_none(('lxml', 'bs4'))
        parser = _BeautifulSoupLxmlFrameParser
        return _run_read_html(parser, *args, **kwargs)


@slow
class TestBs4Html5LibParser(TestBs4LxmlParser):
    def test(self):
        pass

    def run_read_html(self, *args, **kwargs):
        kwargs['flavor'] = 'bs4'
        _skip_if_none(('html5lib', 'bs4'))
        parser = _BeautifulSoupHtml5LibFrameParser
        return _run_read_html(parser, *args, **kwargs)

    @slow
    def test_banklist_header(self):
        def try_remove_ws(x):
            try:
                return _remove_whitespace(x)
            except AttributeError:
                return x

        df = self.run_read_html(self.banklist_data, 'Metcalf',
                                attrs={'id': 'table'})[0]
        ground_truth = read_csv(os.path.join(DATA_PATH, 'banklist.csv'),
                                converters={'Updated Date': Timestamp,
                                            'Closing Date': Timestamp})
        # these will not
        self.assertTupleEqual(df.shape, ground_truth.shape)
        old = ['First Vietnamese American Bank In Vietnamese',
               'Westernbank Puerto Rico En Espanol',
               'R-G Premier Bank of Puerto Rico En Espanol',
               'Eurobank En Espanol', 'Sanderson State Bank En Espanol',
               'Washington Mutual Bank (Including its subsidiary Washington '
               'Mutual Bank FSB)',
               'Silver State Bank En Espanol',
               'AmTrade International BankEn Espanol',
               'Hamilton Bank, NA En Espanol',
               'The Citizens Savings BankPioneer Community Bank, Inc.']
        new = ['First Vietnamese American Bank', 'Westernbank Puerto Rico',
               'R-G Premier Bank of Puerto Rico', 'Eurobank',
               'Sanderson State Bank', 'Washington Mutual Bank',
               'Silver State Bank', 'AmTrade International Bank',
               'Hamilton Bank, NA', 'The Citizens Savings Bank']
        dfnew = df.applymap(try_remove_ws).replace(old, new)
        gtnew = ground_truth.applymap(try_remove_ws)
        converted = dfnew.convert_objects(convert_numeric=True)
        assert_frame_equal(converted.convert_objects(convert_dates='coerce'),
                           gtnew)

    @slow
    def test_gold_canyon(self):
        gc = 'Gold Canyon'
        with open(self.banklist_data, 'r') as f:
            raw_text = f.read()

        self.assertIn(gc, raw_text)
        df = self.run_read_html(self.banklist_data, 'Gold Canyon',
                                attrs={'id': 'table'}, infer_types=False)[0]
        self.assertIn(gc, df.to_string())


def get_elements_from_url(url, flavor, element='table'):
    _skip_if_no('bs4')
    _skip_if_no(flavor)
    from bs4 import BeautifulSoup, SoupStrainer
    strainer = SoupStrainer(element)
    with closing(urlopen(url)) as f:
        soup = BeautifulSoup(f, features=flavor, parse_only=strainer)
    return soup.find_all(element)


@slow
def test_bs4_finds_tables():
    url = ('http://ndb.nal.usda.gov/ndb/foods/show/1732?fg=&man=&'
           'lfacet=&format=&count=&max=25&offset=&sort=&qlookup=spam')
    flavors = 'lxml', 'html5lib'
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')

        for flavor in flavors:
            assert get_elements_from_url(url, flavor, 'table')


def get_lxml_elements(url, element):

    _skip_if_no('lxml')
    from lxml.html import parse
    doc = parse(url)
    return doc.xpath('.//{0}'.format(element))


@slow
def test_lxml_finds_tables():
    url = ('http://ndb.nal.usda.gov/ndb/foods/show/1732?fg=&man=&'
           'lfacet=&format=&count=&max=25&offset=&sort=&qlookup=spam')
    assert get_lxml_elements(url, 'table')


@slow
def test_lxml_finds_tbody():
    url = ('http://ndb.nal.usda.gov/ndb/foods/show/1732?fg=&man=&'
           'lfacet=&format=&count=&max=25&offset=&sort=&qlookup=spam')
    assert get_lxml_elements(url, 'tbody')
