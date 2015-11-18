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
from numpy.random import randint

from pandas.compat import range, lrange, u, unichr
import pandas.compat as compat
from pandas import (Index, Series, DataFrame, isnull, notnull,
                    bdate_range, date_range, MultiIndex)
import pandas.core.common as com

from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm

import pandas.core.strings as strings


class TestStringMethods(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_api(self):

        # GH 6106, GH 9322
        self.assertIs(Series.str, strings.StringMethods)
        self.assertIsInstance(Series(['']).str, strings.StringMethods)

        # GH 9184
        invalid = Series([1])
        with tm.assertRaisesRegexp(AttributeError, "only use .str accessor"):
            invalid.str
        self.assertFalse(hasattr(invalid, 'str'))

    def test_iter(self):
        # GH3638
        strs = 'google', 'wikimedia', 'wikipedia', 'wikitravel'
        ds = Series(strs)

        for s in ds.str:
            # iter must yield a Series
            tm.assertIsInstance(s, Series)

            # indices of each yielded Series should be equal to the index of
            # the original Series
            tm.assert_numpy_array_equal(s.index, ds.index)

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
        tm.assertIsInstance(result, Series)
        tm.assert_almost_equal(result, exp)

        # mixed
        mixed = ['a', NA, 'b', True, datetime.today(), 'foo', None, 1, 2.]
        rs = strings.str_count(mixed, 'a')
        xp = [1, NA, 0, NA, NA, 0, NA, NA, NA]
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.count('a')
        tm.assertIsInstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = [u('foo'), u('foofoo'), NA, u('foooofooofommmfoo')]

        result = strings.str_count(values, 'f[o]+')
        exp = [1, 2, NA, 4]
        tm.assert_almost_equal(result, exp)

        result = Series(values).str.count('f[o]+')
        tm.assertIsInstance(result, Series)
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
        tm.assertIsInstance(rs, Series)
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
        tm.assertIsInstance(rs, Series)
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
        tm.assertIsInstance(rs, Series)
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
        tm.assertIsInstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('om'), NA, u('nom'), u('nom')])

        result = values.str.upper()
        exp = Series([u('OM'), NA, u('NOM'), u('NOM')])
        tm.assert_series_equal(result, exp)

        result = result.str.lower()
        tm.assert_series_equal(result, values)

    def test_capitalize(self):
        values = Series(["FOO", "BAR", NA, "Blah", "blurg"])
        result = values.str.capitalize()
        exp = Series(["Foo", "Bar", NA, "Blah", "Blurg"])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(["FOO", NA, "bar", True, datetime.today(),
                        "blah", None, 1, 2.])
        mixed = mixed.str.capitalize()
        exp = Series(["Foo", NA, "Bar", NA, NA, "Blah", NA, NA, NA])
        tm.assert_almost_equal(mixed, exp)

        # unicode
        values = Series([u("FOO"), NA, u("bar"), u("Blurg")])
        results = values.str.capitalize()
        exp = Series([u("Foo"), NA, u("Bar"), u("Blurg")])
        tm.assert_series_equal(results, exp)

    def test_swapcase(self):
        values = Series(["FOO", "BAR", NA, "Blah", "blurg"])
        result = values.str.swapcase()
        exp = Series(["foo", "bar", NA, "bLAH", "BLURG"])
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(["FOO", NA, "bar", True, datetime.today(),
                        "Blah", None, 1, 2.])
        mixed = mixed.str.swapcase()
        exp = Series(["foo", NA, "BAR", NA, NA, "bLAH", NA, NA, NA])
        tm.assert_almost_equal(mixed, exp)

        # unicode
        values = Series([u("FOO"), NA, u("bar"), u("Blurg")])
        results = values.str.swapcase()
        exp = Series([u("foo"), NA, u("BAR"), u("bLURG")])
        tm.assert_series_equal(results, exp)

    def test_casemethods(self):
        values = ['aaa', 'bbb', 'CCC', 'Dddd', 'eEEE']
        s = Series(values)
        self.assertEqual(s.str.lower().tolist(), [v.lower() for v in values])
        self.assertEqual(s.str.upper().tolist(), [v.upper() for v in values])
        self.assertEqual(s.str.title().tolist(), [v.title() for v in values])
        self.assertEqual(s.str.capitalize().tolist(), [v.capitalize() for v in values])
        self.assertEqual(s.str.swapcase().tolist(), [v.swapcase() for v in values])

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
        tm.assertIsInstance(rs, Series)
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
        tm.assertIsInstance(rs, Series)
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
        tm.assertIsInstance(rs, Series)
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
        tm.assertIsInstance(rs, Series)
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

        # GH9980
        # Index only works with one regex group since
        # multi-group would expand to a frame
        idx = Index(['A1', 'A2', 'A3', 'A4', 'B5'])
        with tm.assertRaisesRegexp(ValueError, "supported"):
            idx.str.extract('([AB])([123])')

        # these should work for both Series and Index
        for klass in [Series, Index]:
            # no groups
            s_or_idx = klass(['A1', 'B2', 'C3'])
            f = lambda: s_or_idx.str.extract('[ABC][123]')
            self.assertRaises(ValueError, f)

            # only non-capturing groups
            f = lambda: s_or_idx.str.extract('(?:[AB]).*')
            self.assertRaises(ValueError, f)

            # single group renames series/index properly
            s_or_idx = klass(['A1', 'A2'])
            result = s_or_idx.str.extract(r'(?P<uno>A)\d')
            tm.assert_equal(result.name, 'uno')
            tm.assert_numpy_array_equal(result, klass(['A', 'A']))

        s = Series(['A1', 'B2', 'C3'])
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

        # one named group
        result = s.str.extract('(?P<letter>[AB])')
        exp = Series(['A', 'B', NA], name='letter')
        tm.assert_series_equal(result, exp)

        # two named groups
        result = s.str.extract('(?P<letter>[AB])(?P<number>[123])')
        exp = DataFrame([['A', '1'], ['B', '2'], [NA, NA]], columns=['letter', 'number'])
        tm.assert_frame_equal(result, exp)

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
        tm.assert_series_equal(empty_int, empty.str.find('a'))
        tm.assert_series_equal(empty_int, empty.str.rfind('a'))
        tm.assert_series_equal(empty_str, empty.str.pad(42))
        tm.assert_series_equal(empty_str, empty.str.center(42))
        tm.assert_series_equal(empty_list, empty.str.split('a'))
        tm.assert_series_equal(empty_list, empty.str.rsplit('a'))
        tm.assert_series_equal(empty_list, empty.str.partition('a', expand=False))
        tm.assert_series_equal(empty_list, empty.str.rpartition('a', expand=False))
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
        tm.assert_series_equal(empty_str, empty.str.isalnum())
        tm.assert_series_equal(empty_str, empty.str.isalpha())
        tm.assert_series_equal(empty_str, empty.str.isdigit())
        tm.assert_series_equal(empty_str, empty.str.isspace())
        tm.assert_series_equal(empty_str, empty.str.islower())
        tm.assert_series_equal(empty_str, empty.str.isupper())
        tm.assert_series_equal(empty_str, empty.str.istitle())
        tm.assert_series_equal(empty_str, empty.str.isnumeric())
        tm.assert_series_equal(empty_str, empty.str.isdecimal())
        tm.assert_series_equal(empty_str, empty.str.capitalize())
        tm.assert_series_equal(empty_str, empty.str.swapcase())
        tm.assert_series_equal(empty_str, empty.str.normalize('NFC'))
        if compat.PY3:
            table = str.maketrans('a', 'b')
        else:
            import string
            table = string.maketrans('a', 'b')
        tm.assert_series_equal(empty_str, empty.str.translate(table))

    def test_empty_str_methods_to_frame(self):
        empty_str = empty = Series(dtype=str)
        empty_df = DataFrame([])
        tm.assert_frame_equal(empty_df, empty.str.partition('a'))
        tm.assert_frame_equal(empty_df, empty.str.rpartition('a'))

    def test_ismethods(self):
        values = ['A', 'b', 'Xy', '4', '3A', '', 'TT', '55', '-', '  ']
        str_s = Series(values)
        alnum_e = [True, True, True, True, True, False, True, True, False, False]
        alpha_e = [True, True, True, False, False, False, True, False, False, False]
        digit_e = [False, False, False, True, False, False, False, True, False, False]
        num_e = [False, False, False, True, False, False, False, True, False, False]
        space_e = [False, False, False, False, False, False, False, False, False, True]
        lower_e = [False, True, False, False, False, False, False, False, False, False]
        upper_e = [True, False, False, False, True, False, True, False, False, False]
        title_e = [True, False, True, False, True, False, False, False, False, False]

        tm.assert_series_equal(str_s.str.isalnum(), Series(alnum_e))
        tm.assert_series_equal(str_s.str.isalpha(), Series(alpha_e))
        tm.assert_series_equal(str_s.str.isdigit(), Series(digit_e))
        tm.assert_series_equal(str_s.str.isspace(), Series(space_e))
        tm.assert_series_equal(str_s.str.islower(), Series(lower_e))
        tm.assert_series_equal(str_s.str.isupper(), Series(upper_e))
        tm.assert_series_equal(str_s.str.istitle(), Series(title_e))

        self.assertEqual(str_s.str.isalnum().tolist(), [v.isalnum() for v in values])
        self.assertEqual(str_s.str.isalpha().tolist(), [v.isalpha() for v in values])
        self.assertEqual(str_s.str.isdigit().tolist(), [v.isdigit() for v in values])
        self.assertEqual(str_s.str.isspace().tolist(), [v.isspace() for v in values])
        self.assertEqual(str_s.str.islower().tolist(), [v.islower() for v in values])
        self.assertEqual(str_s.str.isupper().tolist(), [v.isupper() for v in values])
        self.assertEqual(str_s.str.istitle().tolist(), [v.istitle() for v in values])

    def test_isnumeric(self):
        # 0x00bc: ¼ VULGAR FRACTION ONE QUARTER
        # 0x2605: ★ not number
        # 0x1378: ፸ ETHIOPIC NUMBER SEVENTY
        # 0xFF13: ３ Em 3
        values = ['A', '3', u'¼', u'★', u'፸', u'３', 'four']
        s = Series(values)
        numeric_e = [False, True, True, False, True, True, False]
        decimal_e = [False, True, False, False, False, True, False]
        tm.assert_series_equal(s.str.isnumeric(), Series(numeric_e))
        tm.assert_series_equal(s.str.isdecimal(), Series(decimal_e))

        unicodes = [u'A', u'3', u'¼', u'★', u'፸', u'３', u'four']
        self.assertEqual(s.str.isnumeric().tolist(), [v.isnumeric() for v in unicodes])
        self.assertEqual(s.str.isdecimal().tolist(), [v.isdecimal() for v in unicodes])

        values = ['A', np.nan, u'¼', u'★', np.nan, u'３', 'four']
        s = Series(values)
        numeric_e = [False, np.nan, True, False, np.nan, True, False]
        decimal_e = [False, np.nan, False, False, np.nan, True, False]
        tm.assert_series_equal(s.str.isnumeric(), Series(numeric_e))
        tm.assert_series_equal(s.str.isdecimal(), Series(decimal_e))

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

        # GH9980
        # Index.str does not support get_dummies() as it returns a frame
        with tm.assertRaisesRegexp(TypeError, "not supported"):
            idx = Index(['a|b', 'a|c', 'b|c'])
            idx.str.get_dummies('|')

    def test_join(self):
        values = Series(['a_b_c', 'c_d_e', np.nan, 'f_g_h'])
        result = values.str.split('_').str.join('_')
        tm.assert_series_equal(values, result)

        # mixed
        mixed = Series(['a_b', NA, 'asdf_cas_asdf', True, datetime.today(),
                        'foo', None, 1, 2.])

        rs = Series(mixed).str.split('_').str.join('_')
        xp = Series(['a_b', NA, 'asdf_cas_asdf', NA, NA, 'foo', NA, NA, NA])

        tm.assertIsInstance(rs, Series)
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

        tm.assertIsInstance(rs, Series)
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

        tm.assertIsInstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('fooBAD__barBAD'), NA, u('foo'),
                         u('BAD')])

        result = values.str.findall('BAD[_]*')
        exp = Series([[u('BAD__'), u('BAD')], NA, [], [u('BAD')]])
        tm.assert_almost_equal(result, exp)

    def test_find(self):
        values = Series(['ABCDEFG', 'BCDEFEF', 'DEFGHIJEF', 'EFGHEF', 'XXXX'])
        result = values.str.find('EF')
        tm.assert_series_equal(result, Series([4, 3, 1, 0, -1]))
        expected = np.array([v.find('EF') for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.rfind('EF')
        tm.assert_series_equal(result, Series([4, 5, 7, 4, -1]))
        expected = np.array([v.rfind('EF') for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.find('EF', 3)
        tm.assert_series_equal(result, Series([4, 3, 7, 4, -1]))
        expected = np.array([v.find('EF', 3) for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.rfind('EF', 3)
        tm.assert_series_equal(result, Series([4, 5, 7, 4, -1]))
        expected = np.array([v.rfind('EF', 3) for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.find('EF', 3, 6)
        tm.assert_series_equal(result, Series([4, 3, -1, 4, -1]))
        expected = np.array([v.find('EF', 3, 6) for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.rfind('EF', 3, 6)
        tm.assert_series_equal(result, Series([4, 3, -1, 4, -1]))
        expected = np.array([v.rfind('EF', 3, 6) for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        with tm.assertRaisesRegexp(TypeError, "expected a string object, not int"):
            result = values.str.find(0)

        with tm.assertRaisesRegexp(TypeError, "expected a string object, not int"):
            result = values.str.rfind(0)

    def test_find_nan(self):
        values = Series(['ABCDEFG', np.nan, 'DEFGHIJEF', np.nan, 'XXXX'])
        result = values.str.find('EF')
        tm.assert_series_equal(result, Series([4, np.nan, 1, np.nan, -1]))

        result = values.str.rfind('EF')
        tm.assert_series_equal(result, Series([4, np.nan, 7, np.nan, -1]))

        result = values.str.find('EF', 3)
        tm.assert_series_equal(result, Series([4, np.nan, 7, np.nan, -1]))

        result = values.str.rfind('EF', 3)
        tm.assert_series_equal(result, Series([4, np.nan, 7, np.nan, -1]))

        result = values.str.find('EF', 3, 6)
        tm.assert_series_equal(result, Series([4, np.nan, -1, np.nan, -1]))

        result = values.str.rfind('EF', 3, 6)
        tm.assert_series_equal(result, Series([4, np.nan, -1, np.nan, -1]))

    def test_index(self):
        for klass in [Series, Index]:
            s = klass(['ABCDEFG', 'BCDEFEF', 'DEFGHIJEF', 'EFGHEF'])

            result = s.str.index('EF')
            tm.assert_numpy_array_equal(result, klass([4, 3, 1, 0]))
            expected = np.array([v.index('EF') for v in s.values])
            tm.assert_numpy_array_equal(result.values, expected)

            result = s.str.rindex('EF')
            tm.assert_numpy_array_equal(result, klass([4, 5, 7, 4]))
            expected = np.array([v.rindex('EF') for v in s.values])
            tm.assert_numpy_array_equal(result.values, expected)

            result = s.str.index('EF', 3)
            tm.assert_numpy_array_equal(result, klass([4, 3, 7, 4]))
            expected = np.array([v.index('EF', 3) for v in s.values])
            tm.assert_numpy_array_equal(result.values, expected)

            result = s.str.rindex('EF', 3)
            tm.assert_numpy_array_equal(result, klass([4, 5, 7, 4]))
            expected = np.array([v.rindex('EF', 3) for v in s.values])
            tm.assert_numpy_array_equal(result.values, expected)

            result = s.str.index('E', 4, 8)
            tm.assert_numpy_array_equal(result, klass([4, 5, 7, 4]))
            expected = np.array([v.index('E', 4, 8) for v in s.values])
            tm.assert_numpy_array_equal(result.values, expected)

            result = s.str.rindex('E', 0, 5)
            tm.assert_numpy_array_equal(result, klass([4, 3, 1, 4]))
            expected = np.array([v.rindex('E', 0, 5) for v in s.values])
            tm.assert_numpy_array_equal(result.values, expected)

            with tm.assertRaisesRegexp(ValueError, "substring not found"):
                result = s.str.index('DE')

            with tm.assertRaisesRegexp(TypeError, "expected a string object, not int"):
                result = s.str.index(0)

        # test with nan
        s = Series(['abcb', 'ab', 'bcbe', np.nan])
        result = s.str.index('b')
        tm.assert_numpy_array_equal(result, Series([1, 1, 0, np.nan]))
        result = s.str.rindex('b')
        tm.assert_numpy_array_equal(result, Series([3, 1, 2, np.nan]))

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

        tm.assertIsInstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        mixed = Series(['a', NA, 'b', True, datetime.today(),
                        'ee', None, 1, 2.])

        rs = Series(mixed).str.pad(5, side='right')
        xp = Series(['a    ', NA, 'b    ', NA, NA, 'ee   ', NA, NA, NA])

        tm.assertIsInstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        mixed = Series(['a', NA, 'b', True, datetime.today(),
                        'ee', None, 1, 2.])

        rs = Series(mixed).str.pad(5, side='both')
        xp = Series(['  a  ', NA, '  b  ', NA, NA, '  ee ', NA, NA, NA])

        tm.assertIsInstance(rs, Series)
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

    def test_pad_fillchar(self):

        values = Series(['a', 'b', NA, 'c', NA, 'eeeeee'])

        result = values.str.pad(5, side='left', fillchar='X')
        exp = Series(['XXXXa', 'XXXXb', NA, 'XXXXc', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

        result = values.str.pad(5, side='right', fillchar='X')
        exp = Series(['aXXXX', 'bXXXX', NA, 'cXXXX', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

        result = values.str.pad(5, side='both', fillchar='X')
        exp = Series(['XXaXX', 'XXbXX', NA, 'XXcXX', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

        with tm.assertRaisesRegexp(TypeError, "fillchar must be a character, not str"):
            result = values.str.pad(5, fillchar='XY')

        with tm.assertRaisesRegexp(TypeError, "fillchar must be a character, not int"):
            result = values.str.pad(5, fillchar=5)

    def test_translate(self):
        for klass in [Series, Index]:
            s = klass(['abcdefg', 'abcc', 'cdddfg', 'cdefggg'])
            if not compat.PY3:
                import string
                table = string.maketrans('abc', 'cde')
            else:
                table = str.maketrans('abc', 'cde')
            result = s.str.translate(table)
            expected = klass(['cdedefg', 'cdee', 'edddfg', 'edefggg'])
            tm.assert_numpy_array_equal(result, expected)

            # use of deletechars is python 2 only
            if not compat.PY3:
                result = s.str.translate(table, deletechars='fg')
                expected = klass(['cdede', 'cdee', 'eddd', 'ede'])
                tm.assert_numpy_array_equal(result, expected)

                result = s.str.translate(None, deletechars='fg')
                expected = klass(['abcde', 'abcc', 'cddd', 'cde'])
                tm.assert_numpy_array_equal(result, expected)
            else:
                with tm.assertRaisesRegexp(ValueError, "deletechars is not a valid argument"):
                    result = s.str.translate(table, deletechars='fg')

        # Series with non-string values
        s = Series(['a', 'b', 'c', 1.2])
        expected = Series(['c', 'd', 'e', np.nan])
        result = s.str.translate(table)
        tm.assert_numpy_array_equal(result, expected)

    def test_center_ljust_rjust(self):
        values = Series(['a', 'b', NA, 'c', NA, 'eeeeee'])

        result = values.str.center(5)
        exp = Series(['  a  ', '  b  ', NA, '  c  ', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

        result = values.str.ljust(5)
        exp = Series(['a    ', 'b    ', NA, 'c    ', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

        result = values.str.rjust(5)
        exp = Series(['    a', '    b', NA, '    c', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

        # mixed
        mixed = Series(['a', NA, 'b', True, datetime.today(),
                        'c', 'eee', None, 1, 2.])

        rs = Series(mixed).str.center(5)
        xp = Series(['  a  ', NA, '  b  ', NA, NA, '  c  ', ' eee ', NA, NA,
                     NA])
        tm.assertIsInstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.ljust(5)
        xp = Series(['a    ', NA, 'b    ', NA, NA, 'c    ', 'eee  ', NA, NA,
                     NA])
        tm.assertIsInstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.rjust(5)
        xp = Series(['    a', NA, '    b', NA, NA, '    c', '  eee', NA, NA,
                     NA])
        tm.assertIsInstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        # unicode
        values = Series([u('a'), u('b'), NA, u('c'), NA,
                         u('eeeeee')])

        result = values.str.center(5)
        exp = Series([u('  a  '), u('  b  '), NA, u('  c  '), NA,
                      u('eeeeee')])
        tm.assert_almost_equal(result, exp)

        result = values.str.ljust(5)
        exp = Series([u('a    '), u('b    '), NA, u('c    '), NA,
                      u('eeeeee')])
        tm.assert_almost_equal(result, exp)

        result = values.str.rjust(5)
        exp = Series([u('    a'), u('    b'), NA, u('    c'), NA,
                      u('eeeeee')])
        tm.assert_almost_equal(result, exp)

    def test_center_ljust_rjust_fillchar(self):
        values = Series(['a', 'bb', 'cccc', 'ddddd', 'eeeeee'])

        result = values.str.center(5, fillchar='X')
        expected = Series(['XXaXX', 'XXbbX', 'Xcccc', 'ddddd', 'eeeeee'])
        tm.assert_series_equal(result, expected)
        expected = np.array([v.center(5, 'X') for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.ljust(5, fillchar='X')
        expected = Series(['aXXXX', 'bbXXX', 'ccccX', 'ddddd', 'eeeeee'])
        tm.assert_series_equal(result, expected)
        expected = np.array([v.ljust(5, 'X') for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.rjust(5, fillchar='X')
        expected = Series(['XXXXa', 'XXXbb', 'Xcccc', 'ddddd', 'eeeeee'])
        tm.assert_series_equal(result, expected)
        expected = np.array([v.rjust(5, 'X') for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        # If fillchar is not a charatter, normal str raises TypeError
        # 'aaa'.ljust(5, 'XY')
        # TypeError: must be char, not str
        with tm.assertRaisesRegexp(TypeError, "fillchar must be a character, not str"):
            result = values.str.center(5, fillchar='XY')

        with tm.assertRaisesRegexp(TypeError, "fillchar must be a character, not str"):
            result = values.str.ljust(5, fillchar='XY')

        with tm.assertRaisesRegexp(TypeError, "fillchar must be a character, not str"):
            result = values.str.rjust(5, fillchar='XY')

        with tm.assertRaisesRegexp(TypeError, "fillchar must be a character, not int"):
            result = values.str.center(5, fillchar=1)

        with tm.assertRaisesRegexp(TypeError, "fillchar must be a character, not int"):
            result = values.str.ljust(5, fillchar=1)

        with tm.assertRaisesRegexp(TypeError, "fillchar must be a character, not int"):
            result = values.str.rjust(5, fillchar=1)

    def test_zfill(self):
        values = Series(['1', '22', 'aaa', '333', '45678'])

        result = values.str.zfill(5)
        expected = Series(['00001', '00022', '00aaa', '00333', '45678'])
        tm.assert_series_equal(result, expected)
        expected = np.array([v.zfill(5) for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        result = values.str.zfill(3)
        expected = Series(['001', '022', 'aaa', '333', '45678'])
        tm.assert_series_equal(result, expected)
        expected = np.array([v.zfill(3) for v in values.values])
        tm.assert_numpy_array_equal(result.values, expected)

        values = Series(['1', np.nan, 'aaa', np.nan, '45678'])
        result = values.str.zfill(5)
        expected = Series(['00001', np.nan, '00aaa', np.nan, '45678'])
        tm.assert_series_equal(result, expected)

    def test_split(self):
        values = Series(['a_b_c', 'c_d_e', NA, 'f_g_h'])

        result = values.str.split('_')
        exp = Series([['a', 'b', 'c'], ['c', 'd', 'e'], NA, ['f', 'g', 'h']])
        tm.assert_series_equal(result, exp)

        # more than one char
        values = Series(['a__b__c', 'c__d__e', NA, 'f__g__h'])
        result = values.str.split('__')
        tm.assert_series_equal(result, exp)

        result = values.str.split('__', expand=False)
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(['a_b_c', NA, 'd_e_f', True, datetime.today(),
                        None, 1, 2.])
        result = mixed.str.split('_')
        exp = Series([['a', 'b', 'c'], NA, ['d', 'e', 'f'], NA, NA,
                     NA, NA, NA])
        tm.assertIsInstance(result, Series)
        tm.assert_almost_equal(result, exp)

        result = mixed.str.split('_', expand=False)
        tm.assertIsInstance(result, Series)
        tm.assert_almost_equal(result, exp)

        # unicode
        values = Series([u('a_b_c'), u('c_d_e'), NA, u('f_g_h')])

        result = values.str.split('_')
        exp = Series([[u('a'), u('b'), u('c')],
                      [u('c'), u('d'), u('e')], NA,
                      [u('f'), u('g'), u('h')]])
        tm.assert_series_equal(result, exp)

        result = values.str.split('_', expand=False)
        tm.assert_series_equal(result, exp)

        # regex split
        values = Series([u('a,b_c'), u('c_d,e'), NA, u('f,g,h')])
        result = values.str.split('[,_]')
        exp = Series([[u('a'), u('b'), u('c')],
                      [u('c'), u('d'), u('e')], NA,
                      [u('f'), u('g'), u('h')]])
        tm.assert_series_equal(result, exp)

    def test_rsplit(self):
        values = Series(['a_b_c', 'c_d_e', NA, 'f_g_h'])
        result = values.str.rsplit('_')
        exp = Series([['a', 'b', 'c'], ['c', 'd', 'e'], NA, ['f', 'g', 'h']])
        tm.assert_series_equal(result, exp)

        # more than one char
        values = Series(['a__b__c', 'c__d__e', NA, 'f__g__h'])
        result = values.str.rsplit('__')
        tm.assert_series_equal(result, exp)

        result = values.str.rsplit('__', expand=False)
        tm.assert_series_equal(result, exp)

        # mixed
        mixed = Series(['a_b_c', NA, 'd_e_f', True, datetime.today(),
                        None, 1, 2.])
        result = mixed.str.rsplit('_')
        exp = Series([['a', 'b', 'c'], NA, ['d', 'e', 'f'], NA, NA,
                     NA, NA, NA])
        tm.assertIsInstance(result, Series)
        tm.assert_almost_equal(result, exp)

        result = mixed.str.rsplit('_', expand=False)
        tm.assertIsInstance(result, Series)
        tm.assert_almost_equal(result, exp)

        # unicode
        values = Series([u('a_b_c'), u('c_d_e'), NA, u('f_g_h')])
        result = values.str.rsplit('_')
        exp = Series([[u('a'), u('b'), u('c')],
                      [u('c'), u('d'), u('e')], NA,
                      [u('f'), u('g'), u('h')]])
        tm.assert_series_equal(result, exp)

        result = values.str.rsplit('_', expand=False)
        tm.assert_series_equal(result, exp)

        # regex split is not supported by rsplit
        values = Series([u('a,b_c'), u('c_d,e'), NA, u('f,g,h')])
        result = values.str.rsplit('[,_]')
        exp = Series([[u('a,b_c')],
                      [u('c_d,e')],
                      NA,
                      [u('f,g,h')]])
        tm.assert_series_equal(result, exp)

        # setting max number of splits, make sure it's from reverse
        values = Series(['a_b_c', 'c_d_e', NA, 'f_g_h'])
        result = values.str.rsplit('_', n=1)
        exp = Series([['a_b', 'c'], ['c_d', 'e'], NA, ['f_g', 'h']])
        tm.assert_series_equal(result, exp)

    def test_split_noargs(self):
        # #1859
        s = Series(['Wes McKinney', 'Travis  Oliphant'])
        result = s.str.split()
        expected = ['Travis', 'Oliphant']
        self.assertEqual(result[1], expected)
        result = s.str.rsplit()
        self.assertEqual(result[1], expected)

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

        with tm.assert_produces_warning(FutureWarning):
            result = s.str.split('_', return_type='frame')

        exp = DataFrame({0: Series(['nosplit', 'alsonosplit'])})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_equal_splits', 'with_no_nans'])
        with tm.assert_produces_warning(FutureWarning):
            result = s.str.split('_', return_type='frame')
        exp = DataFrame({0: ['some', 'with'], 1: ['equal', 'no'],
                         2: ['splits', 'nans']})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_unequal_splits', 'one_of_these_things_is_not'])
        with tm.assert_produces_warning(FutureWarning):
            result = s.str.split('_', return_type='frame')
        exp = DataFrame({0: ['some', 'one'], 1: ['unequal', 'of'],
                         2: ['splits', 'these'], 3: [NA, 'things'],
                         4: [NA, 'is'], 5: [NA, 'not']})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_splits', 'with_index'], index=['preserve', 'me'])
        with tm.assert_produces_warning(FutureWarning):
            result = s.str.split('_', return_type='frame')
        exp = DataFrame({0: ['some', 'with'], 1: ['splits', 'index']},
                        index=['preserve', 'me'])
        tm.assert_frame_equal(result, exp)

        with tm.assertRaisesRegexp(ValueError, "expand must be"):
            with tm.assert_produces_warning(FutureWarning):
                s.str.split('_', return_type="some_invalid_type")

    def test_split_to_dataframe_expand(self):
        s = Series(['nosplit', 'alsonosplit'])
        result = s.str.split('_', expand=True)
        exp = DataFrame({0: Series(['nosplit', 'alsonosplit'])})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_equal_splits', 'with_no_nans'])
        result = s.str.split('_', expand=True)
        exp = DataFrame({0: ['some', 'with'], 1: ['equal', 'no'],
                         2: ['splits', 'nans']})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_unequal_splits', 'one_of_these_things_is_not'])
        result = s.str.split('_', expand=True)
        exp = DataFrame({0: ['some', 'one'], 1: ['unequal', 'of'],
                         2: ['splits', 'these'], 3: [NA, 'things'],
                         4: [NA, 'is'], 5: [NA, 'not']})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_splits', 'with_index'], index=['preserve', 'me'])
        result = s.str.split('_', expand=True)
        exp = DataFrame({0: ['some', 'with'], 1: ['splits', 'index']},
                        index=['preserve', 'me'])
        tm.assert_frame_equal(result, exp)

        with tm.assertRaisesRegexp(ValueError, "expand must be"):
            with tm.assert_produces_warning(FutureWarning):
                s.str.split('_', return_type="some_invalid_type")

    def test_split_to_multiindex_expand(self):
        idx = Index(['nosplit', 'alsonosplit'])
        result = idx.str.split('_', expand=True)
        exp = Index([np.array(['nosplit']), np.array(['alsonosplit'])])
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.nlevels, 1)

        idx = Index(['some_equal_splits', 'with_no_nans'])
        result = idx.str.split('_', expand=True)
        exp = MultiIndex.from_tuples([('some', 'equal', 'splits'),
                                      ('with', 'no', 'nans')])
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.nlevels, 3)

        idx = Index(['some_unequal_splits', 'one_of_these_things_is_not'])
        result = idx.str.split('_', expand=True)
        exp = MultiIndex.from_tuples([('some', 'unequal', 'splits', NA, NA, NA),
                                      ('one', 'of', 'these', 'things', 'is', 'not')])
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.nlevels, 6)

        with tm.assertRaisesRegexp(ValueError, "expand must be"):
            with tm.assert_produces_warning(FutureWarning):
                idx.str.split('_', return_type="some_invalid_type")

    def test_rsplit_to_dataframe_expand(self):
        s = Series(['nosplit', 'alsonosplit'])
        result = s.str.rsplit('_', expand=True)
        exp = DataFrame({0: Series(['nosplit', 'alsonosplit'])})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_equal_splits', 'with_no_nans'])
        result = s.str.rsplit('_', expand=True)
        exp = DataFrame({0: ['some', 'with'], 1: ['equal', 'no'],
                         2: ['splits', 'nans']})
        tm.assert_frame_equal(result, exp)

        result = s.str.rsplit('_', expand=True, n=2)
        exp = DataFrame({0: ['some', 'with'], 1: ['equal', 'no'],
                         2: ['splits', 'nans']})
        tm.assert_frame_equal(result, exp)

        result = s.str.rsplit('_', expand=True, n=1)
        exp = DataFrame({0: ['some_equal', 'with_no'],
                         1: ['splits', 'nans']})
        tm.assert_frame_equal(result, exp)

        s = Series(['some_splits', 'with_index'], index=['preserve', 'me'])
        result = s.str.rsplit('_', expand=True)
        exp = DataFrame({0: ['some', 'with'], 1: ['splits', 'index']},
                        index=['preserve', 'me'])
        tm.assert_frame_equal(result, exp)

    def test_rsplit_to_multiindex_expand(self):
        idx = Index(['nosplit', 'alsonosplit'])
        result = idx.str.rsplit('_', expand=True)
        exp = Index([np.array(['nosplit']), np.array(['alsonosplit'])])
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.nlevels, 1)

        idx = Index(['some_equal_splits', 'with_no_nans'])
        result = idx.str.rsplit('_', expand=True)
        exp = MultiIndex.from_tuples([('some', 'equal', 'splits'),
                                      ('with', 'no', 'nans')])
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.nlevels, 3)

        idx = Index(['some_equal_splits', 'with_no_nans'])
        result = idx.str.rsplit('_', expand=True, n=1)
        exp = MultiIndex.from_tuples([('some_equal', 'splits'),
                                      ('with_no', 'nans')])
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.nlevels, 2)

    def test_partition_series(self):
        values = Series(['a_b_c', 'c_d_e', NA, 'f_g_h'])

        result = values.str.partition('_', expand=False)
        exp = Series([['a', '_', 'b_c'], ['c', '_', 'd_e'], NA, ['f', '_', 'g_h']])
        tm.assert_series_equal(result, exp)

        result = values.str.rpartition('_', expand=False)
        exp = Series([['a_b', '_', 'c'], ['c_d', '_', 'e'], NA, ['f_g', '_', 'h']])
        tm.assert_series_equal(result, exp)

        # more than one char
        values = Series(['a__b__c', 'c__d__e', NA, 'f__g__h'])
        result = values.str.partition('__', expand=False)
        exp = Series([['a', '__', 'b__c'], ['c', '__', 'd__e'], NA, ['f', '__', 'g__h']])
        tm.assert_series_equal(result, exp)

        result = values.str.rpartition('__', expand=False)
        exp = Series([['a__b', '__', 'c'], ['c__d', '__', 'e'], NA, ['f__g', '__', 'h']])
        tm.assert_series_equal(result, exp)

        # None
        values = Series(['a b c', 'c d e', NA, 'f g h'])
        result = values.str.partition(expand=False)
        exp = Series([['a', ' ', 'b c'], ['c', ' ', 'd e'], NA, ['f', ' ', 'g h']])
        tm.assert_series_equal(result, exp)

        result = values.str.rpartition(expand=False)
        exp = Series([['a b', ' ', 'c'], ['c d', ' ', 'e'], NA, ['f g', ' ', 'h']])
        tm.assert_series_equal(result, exp)

        # Not splited
        values = Series(['abc', 'cde', NA, 'fgh'])
        result = values.str.partition('_', expand=False)
        exp = Series([['abc', '', ''], ['cde', '', ''], NA, ['fgh', '', '']])
        tm.assert_series_equal(result, exp)

        result = values.str.rpartition('_', expand=False)
        exp = Series([['', '', 'abc'], ['', '', 'cde'], NA, ['', '', 'fgh']])
        tm.assert_series_equal(result, exp)

        # unicode
        values = Series([u('a_b_c'), u('c_d_e'), NA, u('f_g_h')])

        result = values.str.partition('_', expand=False)
        exp = Series([[u('a'), u('_'), u('b_c')], [u('c'), u('_'), u('d_e')],
                      NA, [u('f'), u('_'), u('g_h')]])
        tm.assert_series_equal(result, exp)

        result = values.str.rpartition('_', expand=False)
        exp = Series([[u('a_b'), u('_'), u('c')], [u('c_d'), u('_'), u('e')],
                      NA, [u('f_g'), u('_'), u('h')]])
        tm.assert_series_equal(result, exp)

        # compare to standard lib
        values = Series(['A_B_C', 'B_C_D', 'E_F_G', 'EFGHEF'])
        result = values.str.partition('_', expand=False).tolist()
        self.assertEqual(result, [v.partition('_') for v in values])
        result = values.str.rpartition('_', expand=False).tolist()
        self.assertEqual(result, [v.rpartition('_') for v in values])

    def test_partition_index(self):
        values = Index(['a_b_c', 'c_d_e', 'f_g_h'])

        result = values.str.partition('_', expand=False)
        exp = Index(np.array([('a', '_', 'b_c'), ('c', '_', 'd_e'), ('f', '_', 'g_h')]))
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.nlevels, 1)

        result = values.str.rpartition('_', expand=False)
        exp = Index(np.array([('a_b', '_', 'c'), ('c_d', '_', 'e'), ('f_g', '_', 'h')]))
        tm.assert_index_equal(result, exp)
        self.assertEqual(result.nlevels, 1)

        result = values.str.partition('_')
        exp = Index([('a', '_', 'b_c'), ('c', '_', 'd_e'), ('f', '_', 'g_h')])
        tm.assert_index_equal(result, exp)
        self.assertTrue(isinstance(result, MultiIndex))
        self.assertEqual(result.nlevels, 3)

        result = values.str.rpartition('_')
        exp = Index([('a_b', '_', 'c'), ('c_d', '_', 'e'), ('f_g', '_', 'h')])
        tm.assert_index_equal(result, exp)
        self.assertTrue(isinstance(result, MultiIndex))
        self.assertEqual(result.nlevels, 3)

    def test_partition_to_dataframe(self):
        values = Series(['a_b_c', 'c_d_e', NA, 'f_g_h'])
        result = values.str.partition('_')
        exp = DataFrame({0: ['a', 'c', np.nan, 'f'],
                         1: ['_', '_', np.nan, '_'],
                         2: ['b_c', 'd_e', np.nan, 'g_h']})
        tm.assert_frame_equal(result, exp)

        result = values.str.rpartition('_')
        exp = DataFrame({0: ['a_b', 'c_d', np.nan, 'f_g'],
                         1: ['_', '_', np.nan, '_'],
                         2: ['c', 'e', np.nan, 'h']})
        tm.assert_frame_equal(result, exp)

        values = Series(['a_b_c', 'c_d_e', NA, 'f_g_h'])
        result = values.str.partition('_', expand=True)
        exp = DataFrame({0: ['a', 'c', np.nan, 'f'],
                         1: ['_', '_', np.nan, '_'],
                         2: ['b_c', 'd_e', np.nan, 'g_h']})
        tm.assert_frame_equal(result, exp)

        result = values.str.rpartition('_', expand=True)
        exp = DataFrame({0: ['a_b', 'c_d', np.nan, 'f_g'],
                         1: ['_', '_', np.nan, '_'],
                         2: ['c', 'e', np.nan, 'h']})
        tm.assert_frame_equal(result, exp)

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

        tm.assertIsInstance(rs, Series)
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
        values = Series(['short', 'a bit longer', 'evenlongerthanthat', '', NA])

        exp = Series(['shrt', 'a it longer', 'evnlongerthanthat', '', NA])
        result = values.str.slice_replace(2, 3)
        tm.assert_series_equal(result, exp)

        exp = Series(['shzrt', 'a zit longer', 'evznlongerthanthat', 'z', NA])
        result = values.str.slice_replace(2, 3, 'z')
        tm.assert_series_equal(result, exp)

        exp = Series(['shzort', 'a zbit longer', 'evzenlongerthanthat', 'z', NA])
        result = values.str.slice_replace(2, 2, 'z')
        tm.assert_series_equal(result, exp)

        exp = Series(['shzort', 'a zbit longer', 'evzenlongerthanthat', 'z', NA])
        result = values.str.slice_replace(2, 1, 'z')
        tm.assert_series_equal(result, exp)

        exp = Series(['shorz', 'a bit longez', 'evenlongerthanthaz', 'z', NA])
        result = values.str.slice_replace(-1, None, 'z')
        tm.assert_series_equal(result, exp)

        exp = Series(['zrt', 'zer', 'zat', 'z', NA])
        result = values.str.slice_replace(None, -2, 'z')
        tm.assert_series_equal(result, exp)

        exp = Series(['shortz', 'a bit znger', 'evenlozerthanthat', 'z', NA])
        result = values.str.slice_replace(6, 8, 'z')
        tm.assert_series_equal(result, exp)

        exp = Series(['zrt', 'a zit longer', 'evenlongzerthanthat', 'z', NA])
        result = values.str.slice_replace(-10, 3, 'z')
        tm.assert_series_equal(result, exp)

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

        tm.assertIsInstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.lstrip()
        xp = Series(['aa  ', NA, 'bb \t\n', NA, NA,
                     NA, NA, NA])

        tm.assertIsInstance(rs, Series)
        tm.assert_almost_equal(rs, xp)

        rs = Series(mixed).str.rstrip()
        xp = Series(['  aa', NA, ' bb', NA, NA,
                     NA, NA, NA])

        tm.assertIsInstance(rs, Series)
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

        tm.assertIsInstance(rs, Series)
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

        with tm.assert_produces_warning(FutureWarning):
            result = data.str.match(pat, flags=re.IGNORECASE)
        self.assertEqual(result[0], ('dave', 'google', 'com'))

        result = data.str.findall(pat, flags=re.IGNORECASE)
        self.assertEqual(result[0][0], ('dave', 'google', 'com'))

        result = data.str.count(pat, flags=re.IGNORECASE)
        self.assertEqual(result[0], 1)

        with tm.assert_produces_warning(UserWarning):
            result = data.str.contains(pat, flags=re.IGNORECASE)
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

    def test_normalize(self):
        values = ['ABC', u'ＡＢＣ', u'１２３', np.nan, u'ｱｲｴ']
        s = Series(values, index=['a', 'b', 'c', 'd', 'e'])

        normed = [u'ABC', u'ABC', u'123', np.nan, u'アイエ']
        expected = Series(normed, index=['a', 'b', 'c', 'd', 'e'])

        result = s.str.normalize('NFKC')
        tm.assert_series_equal(result, expected)

        expected = Series([u'ABC', u'ＡＢＣ', u'１２３', np.nan, u'ｱｲｴ'],
                          index=['a', 'b', 'c', 'd', 'e'])

        result = s.str.normalize('NFC')
        tm.assert_series_equal(result, expected)

        with tm.assertRaisesRegexp(ValueError, "invalid normalization form"):
            s.str.normalize('xxx')

        s = Index([u'ＡＢＣ', u'１２３', u'ｱｲｴ'])
        expected = Index([u'ABC', u'123', u'アイエ'])
        result = s.str.normalize('NFKC')
        tm.assert_index_equal(result, expected)

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


    def test_index_str_accessor_visibility(self):
        from pandas.core.strings import StringMethods

        if not compat.PY3:
            cases = [(['a', 'b'], 'string'),
                     (['a', u('b')], 'mixed'),
                     ([u('a'), u('b')], 'unicode'),
                     (['a', 'b', 1], 'mixed-integer'),
                     (['a', 'b', 1.3], 'mixed'),
                     (['a', 'b', 1.3, 1], 'mixed-integer'),
                     (['aa', datetime(2011, 1, 1)], 'mixed')]
        else:
            cases = [(['a', 'b'], 'string'),
                     (['a', u('b')], 'string'),
                     ([u('a'), u('b')], 'string'),
                     (['a', 'b', 1], 'mixed-integer'),
                     (['a', 'b', 1.3], 'mixed'),
                     (['a', 'b', 1.3, 1], 'mixed-integer'),
                     (['aa', datetime(2011, 1, 1)], 'mixed')]
        for values, tp in cases:
            idx = Index(values)
            self.assertTrue(isinstance(Series(values).str, StringMethods))
            self.assertTrue(isinstance(idx.str, StringMethods))
            self.assertEqual(idx.inferred_type, tp)

        for values, tp in cases:
            idx = Index(values)
            self.assertTrue(isinstance(Series(values).str, StringMethods))
            self.assertTrue(isinstance(idx.str, StringMethods))
            self.assertEqual(idx.inferred_type, tp)

        cases = [([1, np.nan], 'floating'),
                 ([datetime(2011, 1, 1)], 'datetime64'),
                 ([timedelta(1)], 'timedelta64')]
        for values, tp in cases:
            idx = Index(values)
            message = 'Can only use .str accessor with string values'
            with self.assertRaisesRegexp(AttributeError, message):
                Series(values).str
            with self.assertRaisesRegexp(AttributeError, message):
                idx.str
            self.assertEqual(idx.inferred_type, tp)

        # MultiIndex has mixed dtype, but not allow to use accessor
        idx = MultiIndex.from_tuples([('a', 'b'), ('a', 'b')])
        self.assertEqual(idx.inferred_type, 'mixed')
        message = 'Can only use .str accessor with Index, not MultiIndex'
        with self.assertRaisesRegexp(AttributeError, message):
            idx.str

    def test_str_accessor_no_new_attributes(self):
        # https://github.com/pydata/pandas/issues/10673
        s = Series(list('aabbcde'))
        with tm.assertRaisesRegexp(AttributeError, "You cannot add any new attribute"):
            s.str.xlabel = "a"

    def test_method_on_bytes(self):
        lhs = Series(np.array(list('abc'), 'S1').astype(object))
        rhs = Series(np.array(list('def'), 'S1').astype(object))
        if compat.PY3:
            self.assertRaises(TypeError, lhs.str.cat, rhs)
        else:
            result = lhs.str.cat(rhs)
            expected = Series(np.array(['ad', 'be', 'cf'],
                                       'S2').astype(object))
            tm.assert_series_equal(result, expected)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
