import os
import re
from cStringIO import StringIO
from unittest import TestCase
import warnings
from distutils.version import LooseVersion

import nose
from nose.tools import assert_raises

import numpy as np
from numpy.random import rand
from numpy.testing.decorators import slow

try:
    from importlib import import_module
except ImportError:
    import_module = __import__

from pandas.io.html import read_html
from pandas.io.common import urlopen

from pandas import DataFrame, MultiIndex, read_csv, Timestamp
from pandas.util.testing import (assert_frame_equal, network,
                                 get_data_path)

from pandas.util.testing import makeCustomDataframe as mkdf


def _have_module(module_name):
    try:
        import_module(module_name)
        return True
    except ImportError:
        return False


def _skip_if_no(module_name):
    if not _have_module(module_name):
        raise nose.SkipTest("{0} not found".format(module_name))


def _skip_if_none_of(module_names):
    if isinstance(module_names, basestring):
        _skip_if_no(module_names)
        if module_names == 'bs4':
            import bs4
            if bs4.__version__ == LooseVersion('4.2.0'):
                raise nose.SkipTest("Bad version of bs4: 4.2.0")
    else:
        not_found = [module_name for module_name in module_names if not
                     _have_module(module_name)]
        if set(not_found) & set(module_names):
            raise nose.SkipTest("{0} not found".format(not_found))
        if 'bs4' in module_names:
            import bs4
            if bs4.__version__ == LooseVersion('4.2.0'):
                raise nose.SkipTest("Bad version of bs4: 4.2.0")


DATA_PATH = get_data_path()

def isframe(x):
    return isinstance(x, DataFrame)


def assert_framelist_equal(list1, list2, *args, **kwargs):
    assert len(list1) == len(list2), ('lists are not of equal size '
                                      'len(list1) == {0}, '
                                      'len(list2) == {1}'.format(len(list1),
                                                                 len(list2)))
    assert all(map(lambda x, y: isframe(x) and isframe(y), list1, list2)), \
        'not all list elements are DataFrames'
    for frame_i, frame_j in zip(list1, list2):
        assert_frame_equal(frame_i, frame_j, *args, **kwargs)
        assert not frame_i.empty, 'frames are both empty'


def test_bs4_version_fails():
    _skip_if_none_of(('bs4', 'html5lib'))
    import bs4
    if bs4.__version__ == LooseVersion('4.2.0'):
        assert_raises(AssertionError, read_html, os.path.join(DATA_PATH,
                                                              "spam.html"),
                      flavor='bs4')


class TestReadHtmlBase(TestCase):
    def run_read_html(self, *args, **kwargs):
        kwargs['flavor'] = kwargs.get('flavor', self.flavor)
        return read_html(*args, **kwargs)

    def try_skip(self):
        _skip_if_none_of(('bs4', 'html5lib'))

    def setup_data(self):
        self.spam_data = os.path.join(DATA_PATH, 'spam.html')
        self.banklist_data = os.path.join(DATA_PATH, 'banklist.html')

    def setup_flavor(self):
        self.flavor = 'bs4'

    def setUp(self):
        self.try_skip()
        self.setup_data()
        self.setup_flavor()

    def test_to_html_compat(self):
        df = mkdf(4, 3, data_gen_f=lambda *args: rand(), c_idx_names=False,
                  r_idx_names=False).applymap('{0:.3f}'.format).astype(float)
        out = df.to_html()
        res = self.run_read_html(out, attrs={'class': 'dataframe'},
                                 index_col=0)[0]
        print (df.dtypes)
        print (res.dtypes)
        assert_frame_equal(res, df)

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

    def test_spam(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*',
                                 infer_types=False)
        df2 = self.run_read_html(self.spam_data, 'Unit', infer_types=False)

        assert_framelist_equal(df1, df2)
        print (df1[0])

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

    def test_header_and_index_no_types(self):
        df1 = self.run_read_html(self.spam_data, '.*Water.*', header=1,
                                 index_col=0, infer_types=False)
        df2 = self.run_read_html(self.spam_data, 'Unit', header=1, index_col=0,
                                 infer_types=False)
        assert_framelist_equal(df1, df2)

    def test_header_and_index_with_types(self):
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

    @network
    def test_multiple_matches(self):
        url = 'http://code.google.com/p/pythonxy/wiki/StandardPlugins'
        dfs = self.run_read_html(url, match='Python',
                                 attrs={'class': 'wikitable'})
        self.assertGreater(len(dfs), 1)

    @network
    def test_pythonxy_plugins_table(self):
        url = 'http://code.google.com/p/pythonxy/wiki/StandardPlugins'
        dfs = self.run_read_html(url, match='Python',
                                 attrs={'class': 'wikitable'})
        zz = [df.iloc[0, 0] for df in dfs]
        self.assertListEqual(sorted(zz), sorted(['Python', 'SciTE']))

    @slow
    def test_banklist_header(self):
        from pandas.io.html import _remove_whitespace
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
        self.assertTupleEqual(df.shape, ground_truth.shape)
        old = ['First Vietnamese American BankIn Vietnamese',
               'Westernbank Puerto RicoEn Espanol',
               'R-G Premier Bank of Puerto RicoEn Espanol',
               'EurobankEn Espanol', 'Sanderson State BankEn Espanol',
               'Washington Mutual Bank(Including its subsidiary Washington '
               'Mutual Bank FSB)',
               'Silver State BankEn Espanol',
               'AmTrade International BankEn Espanol',
               'Hamilton Bank, NAEn Espanol',
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

        self.assert_(gc in raw_text)
        df = self.run_read_html(self.banklist_data, 'Gold Canyon',
                                attrs={'id': 'table'}, infer_types=False)[0]
        self.assertIn(gc, df.to_string())


class TestReadHtmlLxml(TestCase):
    def run_read_html(self, *args, **kwargs):
        self.flavor = ['lxml']
        self.try_skip()
        kwargs['flavor'] = kwargs.get('flavor', self.flavor)
        return read_html(*args, **kwargs)

    def try_skip(self):
        _skip_if_no('lxml')

    def test_spam_data_fail(self):
        from lxml.etree import XMLSyntaxError
        spam_data = os.path.join(DATA_PATH, 'spam.html')
        self.assertRaises(XMLSyntaxError, self.run_read_html, spam_data,
                          flavor=['lxml'])

    def test_banklist_data_fail(self):
        from lxml.etree import XMLSyntaxError
        banklist_data = os.path.join(DATA_PATH, 'banklist.html')
        self.assertRaises(XMLSyntaxError, self.run_read_html, banklist_data, flavor=['lxml'])

    def test_works_on_valid_markup(self):
        filename = os.path.join(DATA_PATH, 'valid_markup.html')
        dfs = self.run_read_html(filename, index_col=0, flavor=['lxml'])
        self.assertIsInstance(dfs, list)
        self.assertIsInstance(dfs[0], DataFrame)

    def setUp(self):
        self.try_skip()

    @slow
    def test_fallback_success(self):
        _skip_if_none_of(('bs4', 'html5lib'))
        banklist_data = os.path.join(DATA_PATH, 'banklist.html')
        self.run_read_html(banklist_data, '.*Water.*', flavor=['lxml',
                                                               'html5lib'])


def test_invalid_flavor():
    url = 'google.com'
    nose.tools.assert_raises(ValueError, read_html, url, 'google',
                             flavor='not a* valid**++ flaver')


def get_elements_from_url(url, element='table', base_url="file://"):
    _skip_if_none_of(('bs4', 'html5lib'))
    url = "".join([base_url, url])
    from bs4 import BeautifulSoup
    with urlopen(url) as f:
        soup = BeautifulSoup(f, features='html5lib')
    return soup.find_all(element)


@slow
def test_bs4_finds_tables():
    filepath = os.path.join(DATA_PATH, "spam.html")
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        assert get_elements_from_url(filepath, 'table')


def get_lxml_elements(url, element):
    _skip_if_no('lxml')
    from lxml.html import parse
    doc = parse(url)
    return doc.xpath('.//{0}'.format(element))


@slow
def test_lxml_finds_tables():
    filepath = os.path.join(DATA_PATH, "spam.html")
    assert get_lxml_elements(filepath, 'table')


@slow
def test_lxml_finds_tbody():
    filepath = os.path.join(DATA_PATH, "spam.html")
    assert get_lxml_elements(filepath, 'tbody')
