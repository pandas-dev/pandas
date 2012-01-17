from datetime import datetime

import unittest

from pandas import Series, DataFrame
from pandas.core.common import notnull, isnull
import pandas.core.common as common
import pandas.util.testing as tm

import numpy as np

def test_notnull():
    assert notnull(1.)
    assert not notnull(None)
    assert not notnull(np.NaN)
    assert not notnull(np.inf)
    assert not notnull(-np.inf)

    float_series = Series(np.random.randn(5))
    obj_series = Series(np.random.randn(5), dtype=object)
    assert(isinstance(notnull(float_series), Series))
    assert(isinstance(notnull(obj_series), Series))

def test_isnull():
    assert not isnull(1.)
    assert isnull(None)
    assert isnull(np.NaN)
    assert isnull(np.inf)
    assert isnull(-np.inf)

    float_series = Series(np.random.randn(5))
    obj_series = Series(np.random.randn(5), dtype=object)
    assert(isinstance(isnull(float_series), Series))
    assert(isinstance(isnull(obj_series), Series))

    # call on DataFrame
    df = DataFrame(np.random.randn(10, 5))
    df['foo'] = 'bar'
    result = isnull(df)
    expected = result.apply(isnull)
    tm.assert_frame_equal(result, expected)

def test_isnull_datetime():
    assert (not isnull(datetime.now()))
    assert notnull(datetime.now())

def test_any_none():
    assert(common._any_none(1, 2, 3, None))
    assert(not common._any_none(1, 2, 3, 4))

def test_all_not_none():
    assert(common._all_not_none(1, 2, 3, 4))
    assert(not common._all_not_none(1, 2, 3, None))
    assert(not common._all_not_none(None, None, None, None))

def test_rands():
    r = common.rands(10)
    assert(len(r) == 10)

def test_adjoin():
    data = [['a', 'b', 'c'],
            ['dd', 'ee', 'ff'],
            ['ggg', 'hhh', 'iii']]
    expected = 'a  dd  ggg\nb  ee  hhh\nc  ff  iii'

    adjoined = common.adjoin(2, *data)

    assert(adjoined == expected)

def test_iterpairs():
    data = [1, 2, 3, 4]
    expected = [(1, 2),
                (2, 3),
                (3, 4)]

    result = list(common.iterpairs(data))

    assert(result == expected)

def test_indent():
    s = 'a b c\nd e f'
    result = common.indent(s, spaces=6)

    assert(result == '      a b c\n      d e f')

def test_banner():
    ban = common.banner('hi')
    assert(ban == ('%s\nhi\n%s' % ('=' * 80, '=' * 80)))

def test_map_indices_py():
    data = [4, 3, 2, 1]
    expected = {4 : 0, 3 : 1, 2 : 2, 1 : 3}

    result = common.map_indices_py(data)

    assert(result == expected)

def test_union():
    a = [1, 2, 3]
    b = [4, 5, 6]

    union = sorted(common.union(a, b))

    assert((a + b) == union)

def test_difference():
    a = [1, 2, 3]
    b = [1, 2, 3, 4, 5, 6]

    inter = sorted(common.difference(b, a))

    assert([4, 5, 6] == inter)

def test_intersection():
    a = [1, 2, 3]
    b = [1, 2, 3, 4, 5, 6]

    inter = sorted(common.intersection(a, b))

    assert(a == inter)

def test_groupby():
    values = ['foo', 'bar', 'baz', 'baz2', 'qux', 'foo3']
    expected = {'f' : ['foo', 'foo3'],
                'b' : ['bar', 'baz', 'baz2'],
                'q' : ['qux']}

    grouped = common.groupby(values, lambda x: x[0])

    for k, v in grouped:
        assert v == expected[k]

class TestTake(unittest.TestCase):

    def test_1d_with_out(self):
        def _test_dtype(dtype):
            out = np.empty(5, dtype=dtype)
            arr = np.random.randn(10).astype(dtype)
            indexer = [0, 2, 4, 7, 1]

            arr.take(indexer, out=out)
            expected = arr.take(indexer)
            tm.assert_almost_equal(out, expected)

        _test_dtype(np.float64)
        _test_dtype(np.float32)
        _test_dtype(np.int32)
        _test_dtype(np.int64)
        _test_dtype(np.object_)
        _test_dtype(np.bool)

    def test_1d_upcast_with_out(self):
        def _test_dtype(dtype):
            out = np.empty(4, dtype=dtype)
            data = np.random.randint(0, 2, 5).astype(dtype)

            indexer = [2, 1, 0, -1]
            self.assertRaises(Exception, common.take_1d, data,
                              indexer, out=out)

        _test_dtype(np.int64)
        _test_dtype(np.int32)
        _test_dtype(np.int16)
        _test_dtype(np.int8)
        _test_dtype(np.bool)

    def test_2d_upcast_with_out(self):
        def _test_dtype(dtype):
            out0 = np.empty((4, 3), dtype=dtype)
            out1 = np.empty((5, 4), dtype=dtype)

            data = np.random.randint(0, 2, (5, 3)).astype(dtype)

            indexer = [2, 1, 0, -1]
            self.assertRaises(Exception, common.take_2d, data,
                              indexer, out=out0, axis=0)
            self.assertRaises(Exception, common.take_2d, data,
                              indexer, out=out1, axis=1)

            # no exception o/w
            data.take(indexer, out=out0, axis=0)
            data.take(indexer, out=out1, axis=1)

        _test_dtype(np.int64)
        _test_dtype(np.int32)
        _test_dtype(np.int16)
        _test_dtype(np.int8)
        _test_dtype(np.bool)

    def test_1d_other_dtypes(self):
        arr = np.random.randn(10).astype(np.float32)

        indexer = [1, 2, 3, -1]
        result = common.take_1d(arr, indexer)
        expected = arr.take(indexer)
        expected[-1] = np.nan
        tm.assert_almost_equal(result, expected)

    def test_2d_other_dtypes(self):
        arr = np.random.randn(10, 5).astype(np.float32)

        indexer = [1, 2, 3, -1]

        # axis=0
        result = common.take_2d(arr, indexer, axis=0)
        expected = arr.take(indexer, axis=0)
        expected[-1] = np.nan
        tm.assert_almost_equal(result, expected)

        # axis=1
        result = common.take_2d(arr, indexer, axis=1)
        expected = arr.take(indexer, axis=1)
        expected[:, -1] = np.nan
        tm.assert_almost_equal(result, expected)

    def test_1d_bool(self):
        arr = np.array([0, 1, 0], dtype=bool)

        result = common.take_1d(arr, [0, 2, 2, 1])
        expected = arr.take([0, 2, 2, 1])
        self.assert_(np.array_equal(result, expected))

        result = common.take_1d(arr, [0, 2, -1])
        self.assert_(result.dtype == np.object_)

    def test_2d_bool(self):
        arr = np.array([[0, 1, 0],
                        [1, 0, 1],
                        [0, 1, 1]], dtype=bool)

        result = common.take_2d(arr, [0, 2, 2, 1])
        expected = arr.take([0, 2, 2, 1], axis=0)
        self.assert_(np.array_equal(result, expected))

        result = common.take_2d(arr, [0, 2, 2, 1], axis=1)
        expected = arr.take([0, 2, 2, 1], axis=1)
        self.assert_(np.array_equal(result, expected))

        result = common.take_2d(arr, [0, 2, -1])
        self.assert_(result.dtype == np.object_)

    def test_2d_float32(self):
        arr = np.random.randn(4, 3).astype(np.float32)
        indexer = [0, 2, -1, 1, -1]

        # axis=0
        result = common.take_2d(arr, indexer)
        result2 = np.empty_like(result)
        common.take_2d(arr, indexer, out=result2)
        tm.assert_almost_equal(result, result)

        expected = arr.take(indexer, axis=0)
        expected[[2, 4]] = np.nan
        tm.assert_almost_equal(result, expected)

        # test with float64 out buffer
        out = np.empty((len(indexer), arr.shape[1]), dtype='f8')
        common.take_2d(arr, indexer, out=out) # it works!

        # axis=1
        result = common.take_2d(arr, indexer, axis=1)
        result2 = np.empty_like(result)
        common.take_2d(arr, indexer, axis=1, out=result2)
        tm.assert_almost_equal(result, result)

        expected = arr.take(indexer, axis=1)
        expected[:, [2, 4]] = np.nan
        tm.assert_almost_equal(result, expected)

class TestEngFormatter(unittest.TestCase):
     def compare(self, formatter, input, output):
         formatted_input = formatter(input)
         msg = "formatting of %s results in '%s', expected '%s'" % (str(input),
                 formatted_input, output)
         self.assertEqual(formatted_input, output, msg)

     def compare_all(self, formatter, in_out):
         """
         Parameters:
         -----------
         formatter: EngFormatter under test
         in_out: list of tuples. Each tuple = (number, expected_formatting)

         It is tested if 'formatter(number) == expected_formatting'.
         *number* should be >= 0 because formatter(-number) == fmt is also
         tested. *fmt* is derived from *expected_formatting*
         """
         for input, output in in_out:
             self.compare(formatter, input, output)
             self.compare(formatter, -input, "-" + output[1:])

     def test_exponents_with_eng_prefix(self):
         formatter = common.EngFormatter(accuracy=3, use_eng_prefix=True)
         f = np.sqrt(2)
         in_out = [(f * 10 ** -24, " 1.414y"),
                   (f * 10 ** -23, " 14.142y"),
                   (f * 10 ** -22, " 141.421y"),
                   (f * 10 ** -21, " 1.414z"),
                   (f * 10 ** -20, " 14.142z"),
                   (f * 10 ** -19, " 141.421z"),
                   (f * 10 ** -18, " 1.414a"),
                   (f * 10 ** -17, " 14.142a"),
                   (f * 10 ** -16, " 141.421a"),
                   (f * 10 ** -15, " 1.414f"),
                   (f * 10 ** -14, " 14.142f"),
                   (f * 10 ** -13, " 141.421f"),
                   (f * 10 ** -12, " 1.414p"),
                   (f * 10 ** -11, " 14.142p"),
                   (f * 10 ** -10, " 141.421p"),
                   (f * 10 ** -9, " 1.414n"),
                   (f * 10 ** -8, " 14.142n"),
                   (f * 10 ** -7, " 141.421n"),
                   (f * 10 ** -6, " 1.414u"),
                   (f * 10 ** -5, " 14.142u"),
                   (f * 10 ** -4, " 141.421u"),
                   (f * 10 ** -3, " 1.414m"),
                   (f * 10 ** -2, " 14.142m"),
                   (f * 10 ** -1, " 141.421m"),
                   (f * 10 ** 0, " 1.414"),
                   (f * 10 ** 1, " 14.142"),
                   (f * 10 ** 2, " 141.421"),
                   (f * 10 ** 3, " 1.414k"),
                   (f * 10 ** 4, " 14.142k"),
                   (f * 10 ** 5, " 141.421k"),
                   (f * 10 ** 6, " 1.414M"),
                   (f * 10 ** 7, " 14.142M"),
                   (f * 10 ** 8, " 141.421M"),
                   (f * 10 ** 9, " 1.414G"),
                   (f * 10 ** 10, " 14.142G"),
                   (f * 10 ** 11, " 141.421G"),
                   (f * 10 ** 12, " 1.414T"),
                   (f * 10 ** 13, " 14.142T"),
                   (f * 10 ** 14, " 141.421T"),
                   (f * 10 ** 15, " 1.414P"),
                   (f * 10 ** 16, " 14.142P"),
                   (f * 10 ** 17, " 141.421P"),
                   (f * 10 ** 18, " 1.414E"),
                   (f * 10 ** 19, " 14.142E"),
                   (f * 10 ** 20, " 141.421E"),
                   (f * 10 ** 21, " 1.414Z"),
                   (f * 10 ** 22, " 14.142Z"),
                   (f * 10 ** 23, " 141.421Z"),
                   (f * 10 ** 24, " 1.414Y"),
                   (f * 10 ** 25, " 14.142Y"),
                   (f * 10 ** 26, " 141.421Y")]
         self.compare_all(formatter, in_out)

     def test_exponents_without_eng_prefix(self):
         formatter = common.EngFormatter(accuracy=4, use_eng_prefix=False)
         f = np.pi
         in_out = [(f * 10 ** -24, " 3.1416E-24"),
                   (f * 10 ** -23, " 31.4159E-24"),
                   (f * 10 ** -22, " 314.1593E-24"),
                   (f * 10 ** -21, " 3.1416E-21"),
                   (f * 10 ** -20, " 31.4159E-21"),
                   (f * 10 ** -19, " 314.1593E-21"),
                   (f * 10 ** -18, " 3.1416E-18"),
                   (f * 10 ** -17, " 31.4159E-18"),
                   (f * 10 ** -16, " 314.1593E-18"),
                   (f * 10 ** -15, " 3.1416E-15"),
                   (f * 10 ** -14, " 31.4159E-15"),
                   (f * 10 ** -13, " 314.1593E-15"),
                   (f * 10 ** -12, " 3.1416E-12"),
                   (f * 10 ** -11, " 31.4159E-12"),
                   (f * 10 ** -10, " 314.1593E-12"),
                   (f * 10 ** -9, " 3.1416E-09"),
                   (f * 10 ** -8, " 31.4159E-09"),
                   (f * 10 ** -7, " 314.1593E-09"),
                   (f * 10 ** -6, " 3.1416E-06"),
                   (f * 10 ** -5, " 31.4159E-06"),
                   (f * 10 ** -4, " 314.1593E-06"),
                   (f * 10 ** -3, " 3.1416E-03"),
                   (f * 10 ** -2, " 31.4159E-03"),
                   (f * 10 ** -1, " 314.1593E-03"),
                   (f * 10 ** 0, " 3.1416E+00"),
                   (f * 10 ** 1, " 31.4159E+00"),
                   (f * 10 ** 2, " 314.1593E+00"),
                   (f * 10 ** 3, " 3.1416E+03"),
                   (f * 10 ** 4, " 31.4159E+03"),
                   (f * 10 ** 5, " 314.1593E+03"),
                   (f * 10 ** 6, " 3.1416E+06"),
                   (f * 10 ** 7, " 31.4159E+06"),
                   (f * 10 ** 8, " 314.1593E+06"),
                   (f * 10 ** 9, " 3.1416E+09"),
                   (f * 10 ** 10, " 31.4159E+09"),
                   (f * 10 ** 11, " 314.1593E+09"),
                   (f * 10 ** 12, " 3.1416E+12"),
                   (f * 10 ** 13, " 31.4159E+12"),
                   (f * 10 ** 14, " 314.1593E+12"),
                   (f * 10 ** 15, " 3.1416E+15"),
                   (f * 10 ** 16, " 31.4159E+15"),
                   (f * 10 ** 17, " 314.1593E+15"),
                   (f * 10 ** 18, " 3.1416E+18"),
                   (f * 10 ** 19, " 31.4159E+18"),
                   (f * 10 ** 20, " 314.1593E+18"),
                   (f * 10 ** 21, " 3.1416E+21"),
                   (f * 10 ** 22, " 31.4159E+21"),
                   (f * 10 ** 23, " 314.1593E+21"),
                   (f * 10 ** 24, " 3.1416E+24"),
                   (f * 10 ** 25, " 31.4159E+24"),
                   (f * 10 ** 26, " 314.1593E+24")]
         self.compare_all(formatter, in_out)

     def test_rounding(self):
         formatter = common.EngFormatter(accuracy=3, use_eng_prefix=True)
         in_out = [(5.55555, ' 5.556'),
                   (55.5555, ' 55.556'),
                   (555.555, ' 555.555'),
                   (5555.55, ' 5.556k'),
                   (55555.5, ' 55.556k'),
                   (555555, ' 555.555k')]
         self.compare_all(formatter, in_out)

         formatter = common.EngFormatter(accuracy=1, use_eng_prefix=True)
         in_out = [(5.55555, ' 5.6'),
                   (55.5555, ' 55.6'),
                   (555.555, ' 555.6'),
                   (5555.55, ' 5.6k'),
                   (55555.5, ' 55.6k'),
                   (555555, ' 555.6k')]
         self.compare_all(formatter, in_out)

         formatter = common.EngFormatter(accuracy=0, use_eng_prefix=True)
         in_out = [(5.55555, ' 6'),
                   (55.5555, ' 56'),
                   (555.555, ' 556'),
                   (5555.55, ' 6k'),
                   (55555.5, ' 56k'),
                   (555555, ' 556k')]
         self.compare_all(formatter, in_out)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

