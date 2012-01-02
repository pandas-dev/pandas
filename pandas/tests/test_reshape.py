from pandas import DataFrame

import numpy as np

from pandas.core.reshape import melt, convert_dummies
import pandas.util.testing as tm

def test_melt():
    df = tm.makeTimeDataFrame()[:10]
    df['id1'] = (df['A'] > 0).astype(int)
    df['id2'] = (df['B'] > 0).astype(int)

    molten1 = melt(df)
    molten2 = melt(df, id_vars=['id1'])
    molten3 = melt(df, id_vars=['id1', 'id2'])

def test_convert_dummies():
    df = DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
                           'foo', 'bar', 'foo', 'foo'],
                    'B' : ['one', 'one', 'two', 'three',
                           'two', 'two', 'one', 'three'],
                    'C' : np.random.randn(8),
                    'D' : np.random.randn(8)})

    result = convert_dummies(df, ['A', 'B'])
    result2 = convert_dummies(df, ['A', 'B'], prefix_sep='.')

    expected = DataFrame({'A_foo' : [1, 0, 1, 0, 1, 0, 1, 1],
                          'A_bar' : [0, 1, 0, 1, 0, 1, 0, 0],
                          'B_one' : [1, 1, 0, 0, 0, 0, 1, 0],
                          'B_two' : [0, 0, 1, 0, 1, 1, 0, 0],
                          'B_three' : [0, 0, 0, 1, 0, 0, 0, 1],
                          'C' : df['C'].values,
                          'D' : df['D'].values},
                         columns=result.columns, dtype=float)
    expected2 = expected.rename(columns=lambda x: x.replace('_', '.'))

    tm.assert_frame_equal(result, expected)
    tm.assert_frame_equal(result2, expected2)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

