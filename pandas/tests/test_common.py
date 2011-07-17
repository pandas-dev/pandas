from pandas.core.common import notnull, isnull
import pandas.core.common as common

import numpy as np

def test_notnull():
    assert notnull(1.)
    assert not notnull(None)
    assert not notnull(np.NaN)
    assert not notnull(np.inf)
    assert not notnull(-np.inf)

def test_isnull():
    assert not isnull(1.)
    assert isnull(None)
    assert isnull(np.NaN)
    assert isnull(np.inf)
    assert isnull(-np.inf)

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

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

