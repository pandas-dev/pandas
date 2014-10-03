#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division, absolute_import, print_function

import sys
import numpy as np
from numpy.testing import *
from numpy.compat import sixu

class TestArrayRepr(object):
    def test_nan_inf(self):
        x = np.array([np.nan, np.inf])
        assert_equal(repr(x), 'array([ nan,  inf])')

class TestComplexArray(TestCase):
    def test_str(self):
        rvals = [0, 1, -1, np.inf, -np.inf, np.nan]
        cvals = [complex(rp, ip) for rp in rvals for ip in rvals]
        dtypes = [np.complex64, np.cdouble, np.clongdouble]
        actual = [str(np.array([c], dt)) for c in cvals for dt in dtypes]
        wanted = [
            '[ 0.+0.j]',    '[ 0.+0.j]',    '[ 0.0+0.0j]',
            '[ 0.+1.j]',    '[ 0.+1.j]',    '[ 0.0+1.0j]',
            '[ 0.-1.j]',    '[ 0.-1.j]',    '[ 0.0-1.0j]',
            '[ 0.+infj]',   '[ 0.+infj]',   '[ 0.0+infj]',
            '[ 0.-infj]',   '[ 0.-infj]',   '[ 0.0-infj]',
            '[ 0.+nanj]',   '[ 0.+nanj]',   '[ 0.0+nanj]',
            '[ 1.+0.j]',    '[ 1.+0.j]',    '[ 1.0+0.0j]',
            '[ 1.+1.j]',    '[ 1.+1.j]',    '[ 1.0+1.0j]',
            '[ 1.-1.j]',    '[ 1.-1.j]',    '[ 1.0-1.0j]',
            '[ 1.+infj]',   '[ 1.+infj]',   '[ 1.0+infj]',
            '[ 1.-infj]',   '[ 1.-infj]',   '[ 1.0-infj]',
            '[ 1.+nanj]',   '[ 1.+nanj]',   '[ 1.0+nanj]',
            '[-1.+0.j]',    '[-1.+0.j]',    '[-1.0+0.0j]',
            '[-1.+1.j]',    '[-1.+1.j]',    '[-1.0+1.0j]',
            '[-1.-1.j]',    '[-1.-1.j]',    '[-1.0-1.0j]',
            '[-1.+infj]',   '[-1.+infj]',   '[-1.0+infj]',
            '[-1.-infj]',   '[-1.-infj]',   '[-1.0-infj]',
            '[-1.+nanj]',   '[-1.+nanj]',   '[-1.0+nanj]',
            '[ inf+0.j]',   '[ inf+0.j]',   '[ inf+0.0j]',
            '[ inf+1.j]',   '[ inf+1.j]',   '[ inf+1.0j]',
            '[ inf-1.j]',   '[ inf-1.j]',   '[ inf-1.0j]',
            '[ inf+infj]',  '[ inf+infj]',  '[ inf+infj]',
            '[ inf-infj]',  '[ inf-infj]',  '[ inf-infj]',
            '[ inf+nanj]',  '[ inf+nanj]',  '[ inf+nanj]',
            '[-inf+0.j]',   '[-inf+0.j]',   '[-inf+0.0j]',
            '[-inf+1.j]',   '[-inf+1.j]',   '[-inf+1.0j]',
            '[-inf-1.j]',   '[-inf-1.j]',   '[-inf-1.0j]',
            '[-inf+infj]',  '[-inf+infj]',  '[-inf+infj]',
            '[-inf-infj]',  '[-inf-infj]',  '[-inf-infj]',
            '[-inf+nanj]',  '[-inf+nanj]',  '[-inf+nanj]',
            '[ nan+0.j]',   '[ nan+0.j]',   '[ nan+0.0j]',
            '[ nan+1.j]',   '[ nan+1.j]',   '[ nan+1.0j]',
            '[ nan-1.j]',   '[ nan-1.j]',   '[ nan-1.0j]',
            '[ nan+infj]',  '[ nan+infj]',  '[ nan+infj]',
            '[ nan-infj]',  '[ nan-infj]',  '[ nan-infj]',
            '[ nan+nanj]',  '[ nan+nanj]',  '[ nan+nanj]']

        for res, val in zip(actual, wanted):
            assert_(res == val)

class TestArray2String(TestCase):
    def test_basic(self):
        """Basic test of array2string."""
        a = np.arange(3)
        assert_(np.array2string(a) == '[0 1 2]')
        assert_(np.array2string(a, max_line_width=4) == '[0 1\n 2]')

    def test_style_keyword(self):
        """This should only apply to 0-D arrays. See #1218."""
        stylestr = np.array2string(np.array(1.5),
                                   style=lambda x: "Value in 0-D array: " + str(x))
        assert_(stylestr == 'Value in 0-D array: 1.5')

    def test_format_function(self):
        """Test custom format function for each element in array."""
        def _format_function(x):
            if np.abs(x) < 1:
                return '.'
            elif np.abs(x) < 2:
                return 'o'
            else:
                return 'O'
        x = np.arange(3)
        if sys.version_info[0] >= 3:
            x_hex = "[0x0 0x1 0x2]"
            x_oct = "[0o0 0o1 0o2]"
        else:
            x_hex = "[0x0L 0x1L 0x2L]"
            x_oct = "[0L 01L 02L]"
        assert_(np.array2string(x, formatter={'all':_format_function}) == \
                "[. o O]")
        assert_(np.array2string(x, formatter={'int_kind':_format_function}) ==\
                "[. o O]")
        assert_(np.array2string(x, formatter={'all':lambda x: "%.4f" % x}) == \
                "[0.0000 1.0000 2.0000]")
        assert_equal(np.array2string(x, formatter={'int':lambda x: hex(x)}), \
                x_hex)
        assert_equal(np.array2string(x, formatter={'int':lambda x: oct(x)}), \
                x_oct)

        x = np.arange(3.)
        assert_(np.array2string(x, formatter={'float_kind':lambda x: "%.2f" % x}) == \
                "[0.00 1.00 2.00]")
        assert_(np.array2string(x, formatter={'float':lambda x: "%.2f" % x}) == \
                "[0.00 1.00 2.00]")

        s = np.array(['abc', 'def'])
        assert_(np.array2string(s, formatter={'numpystr':lambda s: s*2}) == \
            '[abcabc defdef]')


class TestPrintOptions:
    """Test getting and setting global print options."""
    def setUp(self):
        self.oldopts = np.get_printoptions()

    def tearDown(self):
        np.set_printoptions(**self.oldopts)

    def test_basic(self):
        x = np.array([1.5, 0, 1.234567890])
        assert_equal(repr(x), "array([ 1.5       ,  0.        ,  1.23456789])")
        np.set_printoptions(precision=4)
        assert_equal(repr(x), "array([ 1.5   ,  0.    ,  1.2346])")

    def test_formatter(self):
        x = np.arange(3)
        np.set_printoptions(formatter={'all':lambda x: str(x-1)})
        assert_equal(repr(x), "array([-1, 0, 1])")

    def test_formatter_reset(self):
        x = np.arange(3)
        np.set_printoptions(formatter={'all':lambda x: str(x-1)})
        assert_equal(repr(x), "array([-1, 0, 1])")
        np.set_printoptions(formatter={'int':None})
        assert_equal(repr(x), "array([0, 1, 2])")

        np.set_printoptions(formatter={'all':lambda x: str(x-1)})
        assert_equal(repr(x), "array([-1, 0, 1])")
        np.set_printoptions(formatter={'all':None})
        assert_equal(repr(x), "array([0, 1, 2])")

        np.set_printoptions(formatter={'int':lambda x: str(x-1)})
        assert_equal(repr(x), "array([-1, 0, 1])")
        np.set_printoptions(formatter={'int_kind':None})
        assert_equal(repr(x), "array([0, 1, 2])")

        x = np.arange(3.)
        np.set_printoptions(formatter={'float':lambda x: str(x-1)})
        assert_equal(repr(x), "array([-1.0, 0.0, 1.0])")
        np.set_printoptions(formatter={'float_kind':None})
        assert_equal(repr(x), "array([ 0.,  1.,  2.])")

def test_unicode_object_array():
    import sys
    if sys.version_info[0] >= 3:
        expected = "array(['é'], dtype=object)"
    else:
        expected = "array([u'\\xe9'], dtype=object)"
    x = np.array([sixu('\xe9')], dtype=object)
    assert_equal(repr(x), expected)



if __name__ == "__main__":
    run_module_suite()
