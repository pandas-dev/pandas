# pylint: disable-msg=E1101,W0612

from datetime import datetime, timedelta, date
import os
import operator
import unittest

import nose

from numpy import nan as NA
import numpy as np

from pandas import (Index, Series, TimeSeries, DataFrame, isnull, notnull,
                    bdate_range, date_range)
import pandas.core.common as com

from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm

import pandas.core.strings as strings

class TestStringMethods(unittest.TestCase):

    def test_cat(self):
        one = ['a', 'a', 'b', 'b', 'c', NA]
        two = ['a', NA, 'b', 'd', 'foo', NA]

        # single array
        result = strings.str_cat(one)
        self.assert_(isnull(result))

        result = strings.str_cat(one, na_rep='NA')
        exp = 'aabbcNA'
        self.assertEquals(result, exp)

        result = strings.str_cat(one, na_rep='-')
        exp = 'aabbc-'
        self.assertEquals(result, exp)

        result = strings.str_cat(one, sep='_', na_rep='NA')
        exp = 'a_a_b_b_c_NA'
        self.assertEquals(result, exp)

        # Multiple arrays
        result = strings.str_cat(one, [two], na_rep='NA')
        exp = ['aa', 'aNA', 'bb', 'bd', 'cfoo', 'NANA']
        self.assert_(np.array_equal(result, exp))

        result = strings.str_cat(one, two)
        exp = ['aa', NA, 'bb', 'bd', 'cfoo', NA]
        tm.assert_almost_equal(result, exp)

    def test_count(self):
        values = ['foo', 'foofoo', NA, 'foooofooofommmfoo']

        result = strings.str_count(values, 'f[o]+')
        exp = [1, 2, NA, 4]
        tm.assert_almost_equal(result, exp)

        result = Series(values).str.count('f[o]+')
        self.assert_(isinstance(result, Series))
        tm.assert_almost_equal(result, exp)

    def test_contains(self):
        values = ['foo', NA, 'fooommm__foo', 'mmm_']
        pat = 'mmm[_]+'

        result = strings.str_contains(values, pat)
        expected = [False, np.nan, True, True]
        tm.assert_almost_equal(result, expected)

        values = ['foo', 'xyz', 'fooommm__foo', 'mmm_']
        result = strings.str_contains(values, pat)
        expected = [False, False, True, True]
        self.assert_(result.dtype == np.bool_)
        tm.assert_almost_equal(result, expected)

    def test_startswith(self):
        values = Series(['om', NA, 'foo_nom', 'nom', 'bar_foo', NA, 'foo'])

        result = values.str.startswith('foo')
        exp = Series([False, NA, True, False, False, NA, True])
        tm.assert_series_equal(result, exp)

    def test_endswith(self):
        values = Series(['om', NA, 'foo_nom', 'nom', 'bar_foo', NA, 'foo'])

        result = values.str.endswith('foo')
        exp = Series([False, NA, False, False, True, NA, True])
        tm.assert_series_equal(result, exp)

    def test_lower_upper(self):
        values = Series(['om', NA, 'nom', 'nom'])

        result = values.str.upper()
        exp = Series(['OM', NA, 'NOM', 'NOM'])
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

    def test_repeat(self):
        values = Series(['a', 'b', NA, 'c', NA, 'd'])

        result = values.str.repeat(3)
        exp = Series(['aaa', 'bbb', NA, 'ccc', NA, 'ddd'])
        tm.assert_series_equal(result, exp)

        result = values.str.repeat([1, 2, 3, 4, 5, 6])
        exp = Series(['a', 'bb', NA, 'cccc', NA, 'dddddd'])
        tm.assert_series_equal(result, exp)

    def test_match(self):
        values = Series(['fooBAD__barBAD', NA, 'foo'])

        result = values.str.match('.*(BAD[_]+).*(BAD)')
        exp = Series([('BAD__', 'BAD'), NA, []])
        tm.assert_series_equal(result, exp)

    def test_join(self):
        values = Series(['a_b_c', 'c_d_e', np.nan, 'f_g_h'])
        result = values.str.split('_').str.join('_')
        tm.assert_series_equal(values, result)

    def test_len(self):
        values = Series(['foo', 'fooo', 'fooooo', np.nan, 'fooooooo'])

        result = values.str.len()
        exp = values.map(lambda x: len(x) if com.notnull(x) else NA)
        tm.assert_series_equal(result, exp)

    def test_findall(self):
        values = Series(['fooBAD__barBAD', NA, 'foo', 'BAD'])

        result = values.str.findall('BAD[_]*')
        exp = Series([['BAD__', 'BAD'], NA, [], ['BAD']])
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

    def test_center(self):
        values = Series(['a', 'b', NA, 'c', NA, 'eeeeee'])

        result = values.str.center(5)
        exp = Series(['  a  ', '  b  ', NA, '  c  ', NA, 'eeeeee'])
        tm.assert_almost_equal(result, exp)

    def test_split(self):
        values = Series(['a_b_c', 'c_d_e', NA, 'f_g_h'])

        result = values.str.split('_')
        exp = Series([['a', 'b', 'c'], ['c', 'd', 'e'], NA, ['f', 'g', 'h']])
        tm.assert_series_equal(result, exp)

    def test_slice(self):
        values = Series(['aafootwo','aabartwo', NA, 'aabazqux'])

        result = values.str.slice(2, 5)
        exp = Series(['foo', 'bar', NA, 'baz'])
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

    def test_wrap(self):
        pass

    def test_get(self):
        values = Series(['a_b_c', 'c_d_e', np.nan, 'f_g_h'])

        result = values.str.split('_').str.get(1)
        expected = Series(['b', 'd', np.nan, 'g'])
        tm.assert_series_equal(result, expected)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)
