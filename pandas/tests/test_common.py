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

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

