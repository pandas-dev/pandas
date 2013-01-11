# pylint: disable-msg=W0612,E1101
from copy import deepcopy
from datetime import datetime, timedelta
from StringIO import StringIO
import cPickle as pickle
import operator
import os
import unittest

import nose

from pandas import DataFrame

from numpy import nan
import numpy as np

from pandas.core.reshape import melt, convert_dummies, lreshape
import pandas.util.testing as tm

_multiprocess_can_split_ = True


def test_melt():
    df = tm.makeTimeDataFrame()[:10]
    df['id1'] = (df['A'] > 0).astype(int)
    df['id2'] = (df['B'] > 0).astype(int)

    molten1 = melt(df)
    molten2 = melt(df, id_vars=['id1'])
    molten3 = melt(df, id_vars=['id1', 'id2'])
    molten4 = melt(df, id_vars=['id1', 'id2'],
                   value_vars='A')
    molten5 = melt(df, id_vars=['id1', 'id2'],
                   value_vars=['A', 'B'])


def test_convert_dummies():
    df = DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                          'foo', 'bar', 'foo', 'foo'],
                    'B': ['one', 'one', 'two', 'three',
                          'two', 'two', 'one', 'three'],
                    'C': np.random.randn(8),
                    'D': np.random.randn(8)})

    result = convert_dummies(df, ['A', 'B'])
    result2 = convert_dummies(df, ['A', 'B'], prefix_sep='.')

    expected = DataFrame({'A_foo': [1, 0, 1, 0, 1, 0, 1, 1],
                          'A_bar': [0, 1, 0, 1, 0, 1, 0, 0],
                          'B_one': [1, 1, 0, 0, 0, 0, 1, 0],
                          'B_two': [0, 0, 1, 0, 1, 1, 0, 0],
                          'B_three': [0, 0, 0, 1, 0, 0, 0, 1],
                          'C': df['C'].values,
                          'D': df['D'].values},
                         columns=result.columns, dtype=float)
    expected2 = expected.rename(columns=lambda x: x.replace('_', '.'))

    tm.assert_frame_equal(result, expected)
    tm.assert_frame_equal(result2, expected2)


class Test_lreshape(unittest.TestCase):

    def test_pairs(self):
        data = {'birthdt': ['08jan2009', '20dec2008', '30dec2008',
                            '21dec2008', '11jan2009'],
                'birthwt': [1766, 3301, 1454, 3139, 4133],
                'id': [101, 102, 103, 104, 105],
                'sex': ['Male', 'Female', 'Female', 'Female', 'Female'],
                'visitdt1': ['11jan2009', '22dec2008', '04jan2009',
                             '29dec2008', '20jan2009'],
                'visitdt2': ['21jan2009', nan, '22jan2009', '31dec2008', '03feb2009'],
                'visitdt3': ['05feb2009', nan, nan, '02jan2009', '15feb2009'],
                'wt1': [1823, 3338, 1549, 3298, 4306],
                'wt2': [2011.0, nan, 1892.0, 3338.0, 4575.0],
                'wt3': [2293.0, nan, nan, 3377.0, 4805.0]}

        df = DataFrame(data)

        spec = {'visitdt': ['visitdt%d' % i for i in range(1, 4)],
                'wt': ['wt%d' % i for i in range(1, 4)]}
        result = lreshape(df, spec)

        exp_data = {'birthdt': ['08jan2009', '20dec2008', '30dec2008',
                                '21dec2008', '11jan2009', '08jan2009',
                                '30dec2008', '21dec2008', '11jan2009',
                                '08jan2009', '21dec2008', '11jan2009'],
                    'birthwt': [1766, 3301, 1454, 3139, 4133, 1766,
                                1454, 3139, 4133, 1766, 3139, 4133],
                    'id': [101, 102, 103, 104, 105, 101,
                           103, 104, 105, 101, 104, 105],
                    'sex': ['Male', 'Female', 'Female', 'Female', 'Female',
                            'Male', 'Female', 'Female', 'Female', 'Male',
                            'Female', 'Female'],
                    'visitdt': ['11jan2009', '22dec2008', '04jan2009', '29dec2008',
                                '20jan2009', '21jan2009', '22jan2009', '31dec2008',
                                '03feb2009', '05feb2009', '02jan2009', '15feb2009'],
                    'wt': [1823.0, 3338.0, 1549.0, 3298.0, 4306.0, 2011.0,
                           1892.0, 3338.0, 4575.0, 2293.0, 3377.0, 4805.0]}
        exp = DataFrame(exp_data, columns=result.columns)
        tm.assert_frame_equal(result, exp)

        result = lreshape(df, spec, dropna=False)
        exp_data = {'birthdt': ['08jan2009', '20dec2008', '30dec2008',
                                '21dec2008', '11jan2009',
                                '08jan2009', '20dec2008', '30dec2008',
                                '21dec2008', '11jan2009',
                                '08jan2009', '20dec2008', '30dec2008',
                                '21dec2008', '11jan2009'],
                    'birthwt': [1766, 3301, 1454, 3139, 4133,
                                1766, 3301, 1454, 3139, 4133,
                                1766, 3301, 1454, 3139, 4133],
                    'id': [101, 102, 103, 104, 105,
                           101, 102, 103, 104, 105,
                           101, 102, 103, 104, 105],
                    'sex': ['Male', 'Female', 'Female', 'Female', 'Female',
                            'Male', 'Female', 'Female', 'Female', 'Female',
                            'Male', 'Female', 'Female', 'Female', 'Female'],
                    'visitdt': ['11jan2009', '22dec2008', '04jan2009',
                                '29dec2008', '20jan2009',
                                '21jan2009', nan, '22jan2009',
                                '31dec2008', '03feb2009',
                                '05feb2009', nan, nan, '02jan2009', '15feb2009'],
                    'wt': [1823.0, 3338.0, 1549.0, 3298.0, 4306.0, 2011.0,
                           nan, 1892.0, 3338.0, 4575.0, 2293.0, nan, nan,
                           3377.0, 4805.0]}
        exp = DataFrame(exp_data, columns=result.columns)
        tm.assert_frame_equal(result, exp)

        spec = {'visitdt': ['visitdt%d' % i for i in range(1, 3)],
                'wt': ['wt%d' % i for i in range(1, 4)]}
        self.assertRaises(ValueError, lreshape, df, spec)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
