# -*- coding: utf-8 -*-
# pylint: disable-msg=E1101,W0612

from datetime import datetime, timedelta, date
import os
import operator
import re
import warnings

import nose

from numpy import nan as NA
import numpy as np
from numpy.testing import assert_array_equal
from numpy.random import randint

from pandas.compat import range, lrange, u
import pandas.compat as compat
from pandas import (Index, Series, TimeSeries, DataFrame, isnull, notnull,
                    bdate_range, date_range, MultiIndex)
import pandas.core.common as com

from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm

import pandas.core.strings as strings


class TestStringMethods(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_api(self):

        # GH 6106
        self.assertIsNone(Series.str)

    def test_iter(self):
        # GH3638
        strs = 'google', 'wikimedia', 'wikipedia', 'wikitravel'
        ds = Series(strs)

        for s in ds.str:
            # iter must yield a Series
            tm.assert_isinstance(s, Series)

            # indices of each yielded Series should be equal to the index of
            # the original Series
            assert_array_equal(s.index, ds.index)

            for el in s:
                # each element of the series is either a basestring/str or nan
                self.assertTrue(isinstance(el, compat.string_types) or isnull(el))

        # desired behavior is to iterate until everything would be nan on the
        # next iter so make sure the last element of the iterator was 'l' in
        # this case since 'wikitravel' is the longest string
        self.assertEqual(s.dropna().values.item(), 'l')

    def test_iter_empty(self):
        ds = Series([], dtype=object)

        i, s = 100, 1

        for i, s in enumerate(ds.str):
            pass

        # nothing to iterate over so nothing defined values should remain
        # unchanged
        self.assertEqual(i, 100)
        self.assertEqual(s, 1)

    def test_iter_single_element(self):
        ds = Series(['a'])

        for i, s in enumerate(ds.str):
            pass

        self.assertFalse(i)
        assert_series_equal(ds, s)

    def test_iter_numeric_try_string(self):
        # behavior identical to empty series
        dsi = Series(lrange(4))

        i, s = 100, 'h'

        for i, s in enumerate(dsi.str):
            pass

        self.assertEqual(i, 100)
        self.assertEqual(s, 'h')

        dsf = Series(np.arange(4.))

        for i, s in enumerate(dsf.str):
            pass

        self.assertEqual(i, 100)
        self.assertEqual(s, 'h')

    def test_iter_object_try_string(self):
        ds = Series([slice(None, randint(10), randint(10, 20))
                     for _ in range(4)])

        i, s = 100, 'h'

        for i, s in enumerate(ds.str):
            pass

        self.assertEqual(i, 100)
        self.assertEqual(s, 'h')

    def test_cat(self):
        one = ['a', 'a', 'b', 'b', 'c', NA]
        two = ['a', NA, 'b', 'd', 'foo', NA]

        # single array
        result = strings.str_cat(one)
        self.assertTrue(isnull(result))

        result = strings.str_cat(one, na_rep='NA')
        exp = 'aabbcNA'
        self.assertEqual(result, exp)

        result = strings.str_cat(one, na_rep='-')
        exp = 'aabbc-'
        self.assertEqual(result, exp)

        result = strings.str_cat(one, sep='_', na_rep='NA')
        exp = 'a_a_b_b_c_NA'
        self.assertEqual(result, exp)

        # Multiple arrays
        result = strings.str_cat(one, [two], na_rep='NA')
        exp = ['aa', 'aNA', 'bb', 'bd', 'cfoo', 'NANA']
        self.assert_numpy_array_equal(result, exp)

        result = strings.str_cat(one, two)
        exp = ['aa', NA, 'bb', 'bd', 'cfoo', NA]
        tm.assert_almost_equal(result, exp)

    def test_count(self):
        values = ['foo', 'foofoo', NA, 'foooofooofommmfoo']

        result = strings.str_count(values, 'f[o]+')
        exp = [1, 2, NA, 4]
        tm.assert_almost_equal(result, exp)

        result = Series(values).str.count('f[o]+')
        tm.assert_isinstance(result, Series)
        tm.assert_almost_equal(result, exp)

        # mixed
        mixed = ['a', NA, 'b', True, datetime.today(), 'foo', None, 1, 2.]
        rs = strings.str_count(mixed, 'a')
        xp = [1, NA, 0, NA, NA, 0, NA, NA, NA]
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.count('a')
        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = [u('foo'), u('foofoo'), NA, u('foooofooofommmfoo')]

        result = strings.str_count(values, 'f[o]+')
        exp = [1, 2, NA, 4]
        tm.assert_almost_equal(result, exp)

        result = Series(values).str.count('f[o]+')
        tm.assert_isinstance(result, Series)
        tm.assert_almost_equal(result, exp)

    def test_contains(self):
        values = ['foo', NA, 'fooommm__foo', 'mmm_', 'foommm[_]+bar']
        pat = 'mmm[_]+'

        result = strings.str_contains(values, pat)
        expected = [False, NA, True, True, False]
        tm.assert_almost_equal(result, expected)

        result = strings.str_contains(values, pat, regex=False)
        expected = [False, NA, False, False, True]
        tm.assert_almost_equal(result, expected)

        values = ['foo', 'xyz', 'fooommm__foo', 'mmm_']
        result = strings.str_contains(values, pat)
        expected = [False, False, True, True]
        self.assertEqual(result.dtype, np.bool_)
        tm.assert_almost_equal(result, expected)

        # case insensitive using regex
        values = ['Foo', 'xYz', 'fOOomMm__fOo', 'MMM_']
        result = strings.str_contains(values, 'FOO|mmm', case=False)
        expected = [True, False, True, True]
        tm.assert_almost_equal(result, expected)

        # case insensitive without regex
        result = strings.str_contains(values, 'foo', regex=False, case=False)
        expected = [True, False, True, False]
        tm.assert_almost_equal(result, expected)

        # mixed
        mixed = ['a', NA, 'b', True, datetime.today(), 'foo', None, 1, 2.]
        rs = strings.str_contains(mixed, 'o')
        xp = [False, NA, False, NA, NA, True, NA, NA, NA]
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.contains('o')
        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = [u('foo'), NA, u('fooommm__foo'), u('mmm_')]
        pat = 'mmm[_]+'

        result = strings.str_contains(values, pat)
        expected = [False, np.nan, True, True]
        tm.assert_almost_equal(result, expected)

        result = strings.str_contains(values, pat, na=False)
        expected = [False, False, True, True]
        tm.assert_almost_equal(result, expected)

        values = ['foo', 'xyz', 'fooommm__foo', 'mmm_']
        result = strings.str_contains(values, pat)
        expected = [False, False, True, True]
        self.assertEqual(result.dtype, np.bool_)
        tm.assert_almost_equal(result, expected)

        # na
        values = Series(['om', 'foo',np.nan])
        res = values.str.contains('foo', na="foo")
        self.assertEqual (res.ix[2], "foo")

    def test_startswith(self):
        values = Series(['om', NA, 'foo_nom', 'nom', 'bar_foo', NA, 'foo'])

        result = values.str.startswith('foo')
        exp = Series([False, NA, True, False, False, NA, True])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = ['a', NA, 'b', True, datetime.today(), 'foo', None, 1, 2.]
        rs = strings.str_startswith(mixed, 'f')
        xp = [False, NA, False, NA, NA, True, NA, NA, NA]
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.startswith('f')
        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('om'), NA, u('foo_nom'), u('nom'), u('bar_foo'), NA,
                         u('foo')])

        result = values.str.startswith('foo')
        exp = Series([False, NA, True, False, False, NA, True])
        tm.assert_series_equal(result, exp)

        result = values.str.startswith('foo', na=True)
        tm.assert_series_equal(result, exp.fillna(True).astype(bool))

    def test_endswith(self):
        values = Series(['om', NA, 'foo_nom', 'nom', 'bar_foo', NA, 'foo'])

        result = values.str.endswith('foo')
        exp = Series([False, NA, False, False, True, NA, True])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = ['a', NA, 'b', True, datetime.today(), 'foo', None, 1, 2.]
        rs = strings.str_endswith(mixed, 'f')
        xp = [False, NA, False, NA, NA, False, NA, NA, NA]
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.endswith('f')
        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('om'), NA, u('foo_nom'), u('nom'), u('bar_foo'), NA,
                         u('foo')])

        result = values.str.endswith('foo')
        exp = Series([False, NA, False, False, True, NA, True])
        tm.assert_series_equal(result, exp)

        result = values.str.endswith('foo', na=False)
        tm.assert_series_equal(result, exp.fillna(False).astype(bool))

    def test_title(self):
        values = Series(["FOO", "BAR", NA, "Blah", "blurg"])

        result = values.str.title()
        exp = Series(["Foo", "Bar", NA, "Blah", "Blurg"])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(["FOO", NA, "bar", True, datetime.today(),
                        "blah", None, 1, 2.])
        mixed = mixed.str.title()
        exp = Series(["Foo", NA, "Bar", NA, NA, "Blah", NA, NA, NA])
        tm.assert_almost_equal(mixed, exp)

        # unicode
        values = Series([u("FOO"), NA, u("bar"), u("Blurg")])

        results = values.str.title()
        exp = Series([u("Foo"), NA, u("Bar"), u("Blurg")])

        tm.assert_series_equal(results, exp)

    def test_lower_upper(self):
        values = Series(['om', NA, 'nom', 'nom'])

        result = values.str.upper()
        exp = Series(['OM', NA, 'NOM', 'NOM'])
        tm.assert_series_equal(result, exp)

        result = result.str.lower()
        tm.assert_series_equal(result, values)

        # mixed
        mixed = Series(['a', NA, 'b', True, datetime.today(), 'foo', None,
                        1, 2.])
        mixed = mixed.str.upper()
        rs = Series(mixed).str.lower()
        xp = ['a', NA, 'b', NA, NA, 'foo', NA, NA, NA]
        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('om'), NA, u('nom'), u('nom')])

        result = values.str.upper()
        exp = Series([u('OM'), NA, u('NOM'), u('NOM')])
        tm.assert_series_equal(result, exp)

        result = result.str.lower()
        tm.assert_series_equal(result, values)

    def test_replace(self):
        values = Series(['fooBAD__barBAD', NA])

        result = values.str.replace('BAD[_]*', '')
        exp = Series(['foobar', NA])
        tm.assert_series_equal(result, exp)

        result = values.str.replace('BAD[_]*', '', n=1)
        exp = Series(['foobarBAD', NA])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(['aBAD', NA, 'bBAD', True, datetime.today(), 'fooBAD',
                        None, 1, 2.])

        rs = Series(mixed).str.replace('BAD[_]*', '')
        xp = ['a', NA, 'b', NA, NA, 'foo', NA, NA, NA]
        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('fooBAD__barBAD'), NA])

        result = values.str.replace('BAD[_]*', '')
        exp = Series([u('foobar'), NA])
        tm.assert_series_equal(result, exp)

        result = values.str.replace('BAD[_]*', '', n=1)
        exp = Series([u('foobarBAD'), NA])
        tm.assert_series_equal(result, exp)

        #flags + unicode
        values = Series([b"abcd,\xc3\xa0".decode("utf-8")])
        exp = Series([b"abcd, \xc3\xa0".decode("utf-8")])
        result = values.str.replace("(?<=\w),(?=\w)", ", ", flags=re.UNICODE)
        tm.assert_series_equal(result, exp)

    def test_repeat(self):
        values = Series(['a', 'b', NA, 'c', NA, 'd'])

        result = values.str.repeat(3)
        exp = Series(['aaa', 'bbb', NA, 'ccc', NA, 'ddd'])
        tm.assert_series_equal(result, exp)

        result = values.str.repeat([1, 2, 3, 4, 5, 6])
        exp = Series(['a', 'bb', NA, 'cccc', NA, 'dddddd'])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(['a', NA, 'b', True, datetime.today(), 'foo',
                        None, 1, 2.])

        rs = Series(mixed).str.repeat(3)
        xp = ['aaa', NA, 'bbb', NA, NA, 'foofoofoo', NA, NA, NA]
        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('a'), u('b'), NA, u('c'), NA,
                         u('d')])

        result = values.str.repeat(3)
        exp = Series([u('aaa'), u('bbb'), NA, u('ccc'), NA,
                      u('ddd')])
        tm.assert_series_equal(result, exp)

        result = values.str.repeat([1, 2, 3, 4, 5, 6])
        exp = Series([u('a'), u('bb'), NA, u('cccc'), NA,
                      u('dddddd')])
        tm.assert_series_equal(result, exp)

    def test_deprecated_match(self):
        # Old match behavior, deprecated (but still default) in 0.13
        values = Series(['fooBAD__barBAD', NA, 'foo'])

        with tm.assert_produces_warning():
            result = values.str.match('.*(BAD[_]+).*(BAD)')
        exp = Series([('BAD__', 'BAD'), NA, []])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(['aBAD_BAD', NA, 'BAD_b_BAD', True, datetime.today(),
                        'foo', None, 1, 2.])

        with tm.assert_produces_warning():
            rs = Series(mixed).str.match('.*(BAD[_]+).*(BAD)')
        xp = [('BAD_', 'BAD'), NA, ('BAD_', 'BAD'), NA, NA, [], NA, NA, NA]
        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('fooBAD__barBAD'), NA, u('foo')])

        with tm.assert_produces_warning():
            result = values.str.match('.*(BAD[_]+).*(BAD)')
        exp = Series([(u('BAD__'), u('BAD')), NA, []])
        tm.assert_series_equal(result, exp)

    def test_match(self):
        # New match behavior introduced in 0.13
        values = Series(['fooBAD__barBAD', NA, 'foo'])
        with tm.assert_produces_warning():
            result = values.str.match('.*(BAD[_]+).*(BAD)', as_indexer=True)
        exp = Series([True, NA, False])
        tm.assert_series_equal(result, exp)

        # If no groups, use new behavior even when as_indexer is False.
        # (Old behavior is pretty much useless in this case.)
        values = Series(['fooBAD__barBAD', NA, 'foo'])
        result = values.str.match('.*BAD[_]+.*BAD', as_indexer=False)
        exp = Series([True, NA, False])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(['aBAD_BAD', NA, 'BAD_b_BAD', True, datetime.today(),
                        'foo', None, 1, 2.])

        with tm.assert_produces_warning():
            rs = Series(mixed).str.match('.*(BAD[_]+).*(BAD)', as_indexer=True)
        xp = [True, NA, True, NA, NA, False, NA, NA, NA]
        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('fooBAD__barBAD'), NA, u('foo')])

        with tm.assert_produces_warning():
            result = values.str.match('.*(BAD[_]+).*(BAD)', as_indexer=True)
        exp = Series([True, NA, False])
        tm.assert_series_equal(result, exp)

        # na GH #6609
        res = Series(['a', 0, np.nan]).str.match('a', na=False)
        exp = Series([True, False, False])
        assert_series_equal(exp, res)
        res = Series(['a', 0, np.nan]).str.match('a')
        exp = Series([True, np.nan, np.nan])
        assert_series_equal(exp, res)

    def test_extract(self):
        # Contains tests like those in test_match and some others.

        values = Series(['fooBAD__barBAD', NA, 'foo'])
        er = [NA, NA]  # empty row

        result = values.str.extract('.*(BAD[_]+).*(BAD)')
        exp = DataFrame([['BAD__', 'BAD'], er, er])
        tm.assert_frame_equal(result, exp)

        # mixed
        mixed = Series(['aBAD_BAD', NA, 'BAD_b_BAD', True, datetime.today(),
                        'foo', None, 1, 2.])

        rs = Series(mixed).str.extract('.*(BAD[_]+).*(BAD)')
        exp = DataFrame([['BAD_', 'BAD'], er, ['BAD_', 'BAD'], er, er,
                         er, er, er, er])
        tm.assert_frame_equal(rs, exp)

        # unicode
        values = Series([u('fooBAD__barBAD'), NA, u('foo')])

        result = values.str.extract('.*(BAD[_]+).*(BAD)')
        exp = DataFrame([[u('BAD__'), u('BAD')], er, er])
        tm.assert_frame_equal(result, exp)

        # no groups
        s = Series(['A1', 'B2', 'C3'])
        f = lambda: s.str.extract('[ABC][123]')
        self.assertRaises(ValueError, f)

        # only non-capturing groups
        f = lambda: s.str.extract('(?:[AB]).*')
        self.assertRaises(ValueError, f)

        # one group, no matches
        result = s.str.extract('(_)')
        exp = Series([NA, NA, NA], dtype=object)
        tm.assert_series_equal(result, exp)

        # two groups, no matches
        result = s.str.extract('(_)(_)')
        exp = DataFrame([[NA, NA], [NA, NA], [NA, NA]], dtype=object)
        tm.assert_frame_equal(result, exp)

        # one group, some matches
        result = s.str.extract('([AB])[123]')
        exp = Series(['A', 'B', NA])
        tm.assert_series_equal(result, exp)

        # two groups, some matches
        result = s.str.extract('([AB])([123])')
        exp = DataFrame([['A', '1'], ['B', '2'], [NA, NA]])
        tm.assert_frame_equal(result, exp)

        # named group/groups
        result = s.str.extract('(?P<letter>[AB])(?P<number>[123])')
        exp = DataFrame([['A', '1'], ['B', '2'], [NA, NA]], columns=['letter', 'number'])
        tm.assert_frame_equal(result, exp)
        result = s.str.extract('(?P<letter>[AB])')
        exp = Series(['A', 'B', NA], name='letter')
        tm.assert_series_equal(result, exp)

        # mix named and unnamed groups
        result = s.str.extract('([AB])(?P<number>[123])')
        exp = DataFrame([['A', '1'], ['B', '2'], [NA, NA]], columns=[0, 'number'])
        tm.assert_frame_equal(result, exp)

        # one normal group, one non-capturing group
        result = s.str.extract('([AB])(?:[123])')
        exp = Series(['A', 'B', NA])
        tm.assert_series_equal(result, exp)

        # two normal groups, one non-capturing group
        result = Series(['A11', 'B22', 'C33']).str.extract('([AB])([123])(?:[123])')
        exp = DataFrame([['A', '1'], ['B', '2'], [NA, NA]])
        tm.assert_frame_equal(result, exp)

        # one optional group followed by one normal group
        result = Series(['A1', 'B2', '3']).str.extract('(?P<letter>[AB])?(?P<number>[123])')
        exp = DataFrame([['A', '1'], ['B', '2'], [NA, '3']], columns=['letter', 'number'])
        tm.assert_frame_equal(result, exp)

        # one normal group followed by one optional group
        result = Series(['A1', 'B2', 'C']).str.extract('(?P<letter>[ABC])(?P<number>[123])?')
        exp = DataFrame([['A', '1'], ['B', '2'], ['C', NA]], columns=['letter', 'number'])
        tm.assert_frame_equal(result, exp)

        # single group renames series properly
        s = Series(['A1', 'A2'])
        result = s.str.extract(r'(?P<uno>A)\d')
        tm.assert_equal(result.name, 'uno')

        # GH6348
        # not passing index to the extractor
        def check_index(index):
            data = ['A1', 'B2', 'C']
            index = index[:len(data)]
            result = Series(data, index=index).str.extract('(\d)')
            exp = Series(['1', '2', NA], index=index)
            tm.assert_series_equal(result, exp)

            result = Series(data, index=index).str.extract('(?P<letter>\D)(?P<number>\d)?')
            exp = DataFrame([['A', '1'], ['B', '2'], ['C', NA]], columns=['letter', 'number'], index=index)
            tm.assert_frame_equal(result, exp)

        for index in [ tm.makeStringIndex, tm.makeUnicodeIndex, tm.makeIntIndex,
                       tm.makeDateIndex, tm.makePeriodIndex ]:
            check_index(index())

    def test_extract_single_series_name_is_preserved(self):
        s = Series(['a3', 'b3', 'c2'], name='bob')
        r = s.str.extract(r'(?P<sue>[a-z])')
        e = Series(['a', 'b', 'c'], name='sue')
        tm.assert_series_equal(r, e)
        self.assertEqual(r.name, e.name)

    def test_empty_str_methods(self):
        empty_str = empty = Series(dtype=str)
        empty_int = Series(dtype=int)
        empty_bool = Series(dtype=bool)
        empty_list = Series(dtype=list)
        empty_bytes = Series(dtype=object)

        # GH7241
        # (extract) on empty series

        tm.assert_series_equal(empty_str, empty.str.cat(empty))
        tm.assert_equal('', empty.str.cat())
        tm.assert_series_equal(empty_str, empty.str.title())
        tm.assert_series_equal(empty_int, empty.str.count('a'))
        tm.assert_series_equal(empty_bool, empty.str.contains('a'))
        tm.assert_series_equal(empty_bool, empty.str.startswith('a'))
        tm.assert_series_equal(empty_bool, empty.str.endswith('a'))
        tm.assert_series_equal(empty_str, empty.str.lower())
        tm.assert_series_equal(empty_str, empty.str.upper())
        tm.assert_series_equal(empty_str, empty.str.replace('a','b'))
        tm.assert_series_equal(empty_str, empty.str.repeat(3))
        tm.assert_series_equal(empty_bool, empty.str.match('^a'))
        tm.assert_series_equal(empty_str, empty.str.extract('()'))
        tm.assert_frame_equal(DataFrame(columns=[0,1], dtype=str), empty.str.extract('()()'))
        tm.assert_frame_equal(DataFrame(dtype=str), empty.str.get_dummies())
        tm.assert_series_equal(empty_str, empty_list.str.join(''))
        tm.assert_series_equal(empty_int, empty.str.len())
        tm.assert_series_equal(empty_list, empty_list.str.findall('a'))
        tm.assert_series_equal(empty_str, empty.str.pad(42))
        tm.assert_series_equal(empty_str, empty.str.center(42))
        tm.assert_series_equal(empty_list, empty.str.split('a'))
        tm.assert_series_equal(empty_str, empty.str.slice(stop=1))
        tm.assert_series_equal(empty_str, empty.str.slice(step=1))
        tm.assert_series_equal(empty_str, empty.str.strip())
        tm.assert_series_equal(empty_str, empty.str.lstrip())
        tm.assert_series_equal(empty_str, empty.str.rstrip())
        tm.assert_series_equal(empty_str, empty.str.rstrip())
        tm.assert_series_equal(empty_str, empty.str.wrap(42))
        tm.assert_series_equal(empty_str, empty.str.get(0))
        tm.assert_series_equal(empty_str, empty_bytes.str.decode('ascii'))
        tm.assert_series_equal(empty_bytes, empty.str.encode('ascii'))

    def test_get_dummies(self):
        s = Series(['a|b', 'a|c', np.nan])
        result = s.str.get_dummies('|')
        expected = DataFrame([[1, 1, 0], [1, 0, 1], [0, 0, 0]],
                             columns=list('abc'))
        tm.assert_frame_equal(result, expected)

        s = Series(['a;b', 'a', 7])
        result = s.str.get_dummies(';')
        expected = DataFrame([[0, 1, 1], [0, 1, 0], [1, 0, 0]],
                             columns=list('7ab'))
        tm.assert_frame_equal(result, expected)

    def test_join(self):
        values = Series(['a_b_c', 'c_d_e', np.nan, 'f_g_h'])
        result = values.str.split('_').str.join('_')
        tm.assert_series_equal(values, result)

        # mixed
        mixed = Series(['a_b', NA, 'asdf_cas_asdf', True, datetime.today(),
                        'foo', None, 1, 2.])

        rs = Series(mixed).str.split('_').str.join('_')
        xp = Series(['a_b', NA, 'asdf_cas_asdf', NA, NA, 'foo', NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('a_b_c'), u('c_d_e'), np.nan,
                         u('f_g_h')])
        result = values.str.split('_').str.join('_')
        tm.assert_series_equal(values, result)

    def test_len(self):
        values = Series(['foo', 'fooo', 'fooooo', np.nan, 'fooooooo'])

        result = values.str.len()
        exp = values.map(lambda x: len(x) if com.notnull(x) else NA)
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(['a_b', NA, 'asdf_cas_asdf', True, datetime.today(),
                        'foo', None, 1, 2.])

        rs = Series(mixed).str.len()
        xp = Series([3, NA, 13, NA, NA, 3, NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('foo'), u('fooo'), u('fooooo'), np.nan,
                         u('fooooooo')])

        result = values.str.len()
        exp = values.map(lambda x: len(x) if com.notnull(x) else NA)
        tm.assert_series_equal(result, exp)

    def test_findall(self):
        values = Series(['fooBAD__barBAD', NA, 'foo', 'BAD'])

        result = values.str.findall('BAD[_]*')
        exp = Series([['BAD__', 'BAD'], NA, [], ['BAD']])
        tm.assert_almost_equal(result, exp)

        # mixed
        mixed = Series(['fooBAD__barBAD', NA, 'foo', True, datetime.today(),
                        'BAD', None, 1, 2.])

        rs = Series(mixed).str.findall('BAD[_]*')
        xp = Series([['BAD__', 'BAD'], NA, [], NA, NA, ['BAD'], NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('fooBAD__barBAD'), NA, u('foo'),
                         u('BAD')])

        result = values.str.findall('BAD[_]*')
        exp = Series([[u('BAD__'), u('BAD')], NA, [], [u('BAD')]])
        tm.assert_almost_equal(result, exp)

    def test_pad(self):
        values = Series(['a', 'b', NA, 'c', NA, 'eeeeee'])

        result = values.str.pad(5, side='left')
        exp = Series(['    a', '    b', NA, '    c', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

        result = values.str.pad(5, side='right')
        exp = Series(['a    ', 'b    ', NA, 'c    ', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

        result = values.str.pad(5, side='both')
        exp = Series(['  a  ', '  b  ', NA, '  c  ', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

        # mixed
        mixed = Series(['a', NA, 'b', True, datetime.today(),
                        'ee', None, 1, 2.])

        rs = Series(mixed).str.pad(5, side='left')
        xp = Series(['    a', NA, '    b', NA, NA, '   ee', NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        mixed = Series(['a', NA, 'b', True, datetime.today(),
                        'ee', None, 1, 2.])

        rs = Series(mixed).str.pad(5, side='right')
        xp = Series(['a    ', NA, 'b    ', NA, NA, 'ee   ', NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        mixed = Series(['a', NA, 'b', True, datetime.today(),
                        'ee', None, 1, 2.])

        rs = Series(mixed).str.pad(5, side='both')
        xp = Series(['  a  ', NA, '  b  ', NA, NA, '  ee ', NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('a'), u('b'), NA, u('c'), NA,
                         u('eeeeee')])

        result = values.str.pad(5, side='left')
        exp = Series([u('    a'), u('    b'), NA, u('    c'), NA,
                      u('eeeeee')])
        tm.assert_almost_equal(result, exp)

        result = values.str.pad(5, side='right')
        exp = Series([u('a    '), u('b    '), NA, u('c    '), NA,
                      u('eeeeee')])
        tm.assert_almost_equal(result, exp)

        result = values.str.pad(5, side='both')
        exp = Series([u('  a  '), u('  b  '), NA, u('  c  '), NA,
                      u('eeeeee')])
        tm.assert_almost_equal(result, exp)

    def test_center(self):
        values = Series(['a', 'b', NA, 'c', NA, 'eeeeee'])

        result = values.str.center(5)
        exp = Series(['  a  ', '  b  ', NA, '  c  ', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

        # mixed
        mixed = Series(['a', NA, 'b', True, datetime.today(),
                        'c', 'eee', None, 1, 2.])

        rs = Series(mixed).str.center(5)
        xp = Series(['  a  ', NA, '  b  ', NA, NA, '  c  ', ' eee ', NA, NA,
                     NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('a'), u('b'), NA, u('c'), NA,
                         u('eeeeee')])

        result = values.str.center(5)
        exp = Series([u('  a  '), u('  b  '), NA, u('  c  '), NA,
                      u('eeeeee')])
        tm.assert_almost_equal(result, exp)

    def test_split(self):
        values = Series(['a_b_c', 'c_d_e', NA, 'f_g_h'])

        result = values.str.split('_')
        exp = Series([['a', 'b', 'c'], ['c', 'd', 'e'], NA, ['f', 'g', 'h']])
        tm.assert_series_equal(result, exp)

        # more than one char
        values = Series(['a__b__c', 'c__d__e', NA, 'f__g__h'])
        result = values.str.split('__')
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(['a_b_c', NA, 'd_e_f', True, datetime.today(),
                        None, 1, 2.])

        rs = Series(mixed).str.split('_')
        xp = Series([['a', 'b', 'c'], NA, ['d', 'e', 'f'], NA, NA,
                     NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('a_b_c'), u('c_d_e'), NA, u('f_g_h')])

        result = values.str.split('_')
        exp = Series([[u('a'), u('b'), u('c')],
                      [u('c'), u('d'), u('e')], NA,
                      [u('f'), u('g'), u('h')]])
        tm.assert_series_equal(result, exp)

    def test_split_noargs(self):
        # #1859
        s = Series(['Wes McKinney', 'Travis  Oliphant'])

        result = s.str.split()
        self.assertEqual(result[1], ['Travis', 'Oliphant'])

    def test_split_maxsplit(self):
        # re.split 0, str.split -1
        s = Series(['bd asdf jfg', 'kjasdflqw asdfnfk'])

        result = s.str.split(n=-1)
        xp = s.str.split()
        tm.assert_series_equal(result, xp)

        result = s.str.split(n=0)
        tm.assert_series_equal(result, xp)

        xp = s.str.split('asdf')
        result = s.str.split('asdf', n=0)
        tm.assert_series_equal(result, xp)

        result = s.str.split('asdf', n=-1)
        tm.assert_series_equal(result, xp)

    def test_split_no_pat_with_nonzero_n(self):
        s = Series(['split once', 'split once too!'])
        result = s.str.split(n=1)
        expected = Series({0: ['split', 'once'], 1: ['split', 'once too!']})
        tm.assert_series_equal(expected, result)

    def test_split_to_dataframe(self):
        s = Series(['nosplit', 'alsonosplit'])
        result = s.str.split('_', return_type='frame')
        exp = DataFrame({0: Series(['nosplit', 'alsonosplit'])})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_equal_splits', 'with_no_nans'])
        result = s.str.split('_', return_type='frame')
        exp = DataFrame({0: ['some', 'with'], 1: ['equal', 'no'],
                         2: ['splits', 'nans']})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_unequal_splits', 'one_of_these_things_is_not'])
        result = s.str.split('_', return_type='frame')
        exp = DataFrame({0: ['some', 'one'], 1: ['unequal', 'of'],
                         2: ['splits', 'these'], 3: [NA, 'things'],
                         4: [NA, 'is'], 5: [NA, 'not']})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_splits', 'with_index'], index=['preserve', 'me'])
        result = s.str.split('_', return_type='frame')
        exp = DataFrame({0: ['some', 'with'], 1: ['splits', 'index']},
                        index=['preserve', 'me'])
        tm.assert_frame_equal(result, exp)

        with tm.assertRaisesRegexp(ValueError, "return_type must be"):
            s.str.split('_', return_type="some_invalid_type")

    def test_pipe_failures(self):
        # #2119
        s = Series(['A|B|C'])

        result = s.str.split('|')
        exp = Series([['A', 'B', 'C']])

        tm.assert_series_equal(result, exp)

        result = s.str.replace('|', ' ')
        exp = Series(['A B C'])

        tm.assert_series_equal(result, exp)

    def test_slice(self):
        values = Series(['aafootwo', 'aabartwo', NA, 'aabazqux'])

        result = values.str.slice(2, 5)
        exp = Series(['foo', 'bar', NA, 'baz'])
        tm.assert_series_equal(result, exp)

        for start, stop, step in [(0, 3, -1), (None, None, -1),
                                  (3, 10, 2), (3, 0, -1)]:
            try:
                result = values.str.slice(start, stop, step)
                expected = Series([s[start:stop:step] if not isnull(s) else NA for s in
                                   values])
                tm.assert_series_equal(result, expected)
            except:
                print('failed on %s:%s:%s' % (start, stop, step))
                raise

        # mixed
        mixed = Series(['aafootwo', NA, 'aabartwo', True, datetime.today(),
                        None, 1, 2.])

        rs = Series(mixed).str.slice(2, 5)
        xp = Series(['foo', NA, 'bar', NA, NA,
                     NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.slice(2, 5, -1)
        xp = Series(['oof', NA, 'rab', NA, NA,
                     NA, NA, NA])

        # unicode
        values = Series([u('aafootwo'), u('aabartwo'), NA,
                         u('aabazqux')])

        result = values.str.slice(2, 5)
        exp = Series([u('foo'), u('bar'), NA, u('baz')])
        tm.assert_series_equal(result, exp)

        result = values.str.slice(0, -1, 2)
        exp = Series([u('afow'), u('abrw'), NA, u('abzu')])
        tm.assert_series_equal(result, exp)

    def test_slice_replace(self):
        pass

    def test_strip_lstrip_rstrip(self):
        values = Series(['  aa   ', ' bb \n', NA, 'cc  '])

        result = values.str.strip()
        exp = Series(['aa', 'bb', NA, 'cc'])
        tm.assert_series_equal(result, exp)

        result = values.str.lstrip()
        exp = Series(['aa   ', 'bb \n', NA, 'cc  '])
        tm.assert_series_equal(result, exp)

        result = values.str.rstrip()
        exp = Series(['  aa', ' bb', NA, 'cc'])
        tm.assert_series_equal(result, exp)

    def test_strip_lstrip_rstrip_mixed(self):
        # mixed
        mixed = Series(['  aa  ', NA, ' bb \t\n', True, datetime.today(),
                        None, 1, 2.])

        rs = Series(mixed).str.strip()
        xp = Series(['aa', NA, 'bb', NA, NA,
                     NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.lstrip()
        xp = Series(['aa  ', NA, 'bb \t\n', NA, NA,
                     NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.rstrip()
        xp = Series(['  aa', NA, ' bb', NA, NA,
                     NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

    def test_strip_lstrip_rstrip_unicode(self):
        # unicode
        values = Series([u('  aa   '), u(' bb \n'), NA,
                         u('cc  ')])

        result = values.str.strip()
        exp = Series([u('aa'), u('bb'), NA, u('cc')])
        tm.assert_series_equal(result, exp)

        result = values.str.lstrip()
        exp = Series([u('aa   '), u('bb \n'), NA, u('cc  ')])
        tm.assert_series_equal(result, exp)

        result = values.str.rstrip()
        exp = Series([u('  aa'), u(' bb'), NA, u('cc')])
        tm.assert_series_equal(result, exp)

    def test_strip_lstrip_rstrip_args(self):
        values = Series(['xxABCxx', 'xx BNSD', 'LDFJH xx'])

        rs = values.str.strip('x')
        xp = Series(['ABC', ' BNSD', 'LDFJH '])
        assert_series_equal(rs, xp)

        rs = values.str.lstrip('x')
        xp = Series(['ABCxx', ' BNSD', 'LDFJH xx'])
        assert_series_equal(rs, xp)

        rs = values.str.rstrip('x')
        xp = Series(['xxABC', 'xx BNSD', 'LDFJH '])
        assert_series_equal(rs, xp)

    def test_strip_lstrip_rstrip_args_unicode(self):
        values = Series([u('xxABCxx'), u('xx BNSD'),
                         u('LDFJH xx')])

        rs = values.str.strip(u('x'))
        xp = Series(['ABC', ' BNSD', 'LDFJH '])
        assert_series_equal(rs, xp)

        rs = values.str.lstrip(u('x'))
        xp = Series(['ABCxx', ' BNSD', 'LDFJH xx'])
        assert_series_equal(rs, xp)

        rs = values.str.rstrip(u('x'))
        xp = Series(['xxABC', 'xx BNSD', 'LDFJH '])
        assert_series_equal(rs, xp)

    def test_wrap(self):
        # test values are: two words less than width, two words equal to width,
        # two words greater than width, one word less than width, one word
        # equal to width, one word greater than width, multiple tokens with trailing
        # whitespace equal to width
        values = Series([u('hello world'), u('hello world!'),
                         u('hello world!!'), u('abcdefabcde'),
                         u('abcdefabcdef'), u('abcdefabcdefa'),
                         u('ab ab ab ab '), u('ab ab ab ab a'),
                         u('\t')])

        # expected values
        xp = Series([u('hello world'), u('hello world!'),
                     u('hello\nworld!!'), u('abcdefabcde'),
                     u('abcdefabcdef'), u('abcdefabcdef\na'),
                     u('ab ab ab ab'), u('ab ab ab ab\na'),
                     u('')])

        rs = values.str.wrap(12, break_long_words=True)
        assert_series_equal(rs, xp)

        # test with pre and post whitespace (non-unicode), NaN, and non-ascii Unicode
        values = Series(['  pre  ', np.nan, u('\xac\u20ac\U00008000 abadcafe')])
        xp = Series(['  pre', NA, u('\xac\u20ac\U00008000 ab\nadcafe')])
        rs = values.str.wrap(6)
        assert_series_equal(rs, xp)

    def test_get(self):
        values = Series(['a_b_c', 'c_d_e', np.nan, 'f_g_h'])

        result = values.str.split('_').str.get(1)
        expected = Series(['b', 'd', np.nan, 'g'])
        tm.assert_series_equal(result, expected)

        # mixed
        mixed = Series(['a_b_c', NA, 'c_d_e', True, datetime.today(),
                        None, 1, 2.])

        rs = Series(mixed).str.split('_').str.get(1)
        xp = Series(['b', NA, 'd', NA, NA,
                     NA, NA, NA])

        tm.assert_isinstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('a_b_c'), u('c_d_e'), np.nan,
                         u('f_g_h')])

        result = values.str.split('_').str.get(1)
        expected = Series([u('b'), u('d'), np.nan, u('g')])
        tm.assert_series_equal(result, expected)

    def test_more_contains(self):
        # PR #1179
        import re

        s = Series(['A', 'B', 'C', 'Aaba', 'Baca', '', NA,
                    'CABA', 'dog', 'cat'])

        result = s.str.contains('a')
        expected = Series([False, False, False, True, True, False, np.nan,
                           False, False, True])
        assert_series_equal(result, expected)

        result = s.str.contains('a', case=False)
        expected = Series([True, False, False, True, True, False, np.nan,
                           True, False, True])
        assert_series_equal(result, expected)

        result = s.str.contains('Aa')
        expected = Series([False, False, False, True, False, False, np.nan,
                           False, False, False])
        assert_series_equal(result, expected)

        result = s.str.contains('ba')
        expected = Series([False, False, False, True, False, False, np.nan,
                           False, False, False])
        assert_series_equal(result, expected)

        result = s.str.contains('ba', case=False)
        expected = Series([False, False, False, True, True, False, np.nan,
                           True, False, False])
        assert_series_equal(result, expected)

    def test_more_replace(self):
        # PR #1179
        import re
        s = Series(['A', 'B', 'C', 'Aaba', 'Baca',
                    '', NA, 'CABA', 'dog', 'cat'])

        result = s.str.replace('A', 'YYY')
        expected = Series(['YYY', 'B', 'C', 'YYYaba', 'Baca', '', NA,
                           'CYYYBYYY', 'dog', 'cat'])
        assert_series_equal(result, expected)

        result = s.str.replace('A', 'YYY', case=False)
        expected = Series(['YYY', 'B', 'C', 'YYYYYYbYYY', 'BYYYcYYY', '', NA,
                           'CYYYBYYY', 'dog', 'cYYYt'])
        assert_series_equal(result, expected)

        result = s.str.replace('^.a|dog', 'XX-XX ', case=False)
        expected = Series(['A', 'B', 'C', 'XX-XX ba', 'XX-XX ca', '', NA,
                           'XX-XX BA', 'XX-XX ', 'XX-XX t'])
        assert_series_equal(result, expected)

    def test_string_slice_get_syntax(self):
        s = Series(['YYY', 'B', 'C', 'YYYYYYbYYY', 'BYYYcYYY', NA,
                    'CYYYBYYY', 'dog', 'cYYYt'])

        result = s.str[0]
        expected = s.str.get(0)
        assert_series_equal(result, expected)

        result = s.str[:3]
        expected = s.str.slice(stop=3)
        assert_series_equal(result, expected)

        result = s.str[2::-1]
        expected = s.str.slice(start=2, step=-1)
        assert_series_equal(result, expected)

    def test_string_slice_out_of_bounds(self):
        s = Series([(1, 2), (1,), (3,4,5)])

        result = s.str[1]
        expected = Series([2, np.nan, 4])

        assert_series_equal(result, expected)

        s = Series(['foo', 'b', 'ba'])
        result = s.str[1]
        expected = Series(['o', np.nan, 'a'])
        assert_series_equal(result, expected)

    def test_match_findall_flags(self):
        data = {'Dave': 'dave@google.com', 'Steve': 'steve@gmail.com',
                'Rob': 'rob@gmail.com', 'Wes': np.nan}
        data = Series(data)

        pat = pattern = r'([A-Z0-9._%+-]+)@([A-Z0-9.-]+)\.([A-Z]{2,4})'

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            result = data.str.match(pat, flags=re.IGNORECASE)
            assert issubclass(w[-1].category, UserWarning)
        self.assertEqual(result[0], ('dave', 'google', 'com'))

        result = data.str.findall(pat, flags=re.IGNORECASE)
        self.assertEqual(result[0][0], ('dave', 'google', 'com'))

        result = data.str.count(pat, flags=re.IGNORECASE)
        self.assertEqual(result[0], 1)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            result = data.str.contains(pat, flags=re.IGNORECASE)
            assert issubclass(w[-1].category, UserWarning)
        self.assertEqual(result[0], True)

    def test_encode_decode(self):
        base = Series([u('a'), u('b'), u('a\xe4')])
        series = base.str.encode('utf-8')

        f = lambda x: x.decode('utf-8')
        result = series.str.decode('utf-8')
        exp = series.map(f)

        tm.assert_series_equal(result, exp)

    def test_encode_decode_errors(self):
        encodeBase = Series([u('a'), u('b'), u('a\x9d')])

        self.assertRaises(UnicodeEncodeError,
                          encodeBase.str.encode, 'cp1252')

        f = lambda x: x.encode('cp1252', 'ignore')
        result = encodeBase.str.encode('cp1252', 'ignore')
        exp = encodeBase.map(f)
        tm.assert_series_equal(result, exp)

        decodeBase = Series([b'a', b'b', b'a\x9d'])

        self.assertRaises(UnicodeDecodeError,
                          decodeBase.str.decode, 'cp1252')

        f = lambda x: x.decode('cp1252', 'ignore')
        result = decodeBase.str.decode('cp1252', 'ignore')
        exp = decodeBase.map(f)

        tm.assert_series_equal(result, exp)

    def test_cat_on_filtered_index(self):
        df = DataFrame(index=MultiIndex.from_product([[2011, 2012], [1,2,3]],
                                                     names=['year', 'month']))

        df = df.reset_index()
        df = df[df.month > 1]

        str_year = df.year.astype('str')
        str_month = df.month.astype('str')
        str_both = str_year.str.cat(str_month, sep=' ')

        self.assertEqual(str_both.loc[1], '2011 2')

        str_multiple = str_year.str.cat([str_month, str_month], sep=' ')

        self.assertEqual(str_multiple.loc[1], '2011 2 2')


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
