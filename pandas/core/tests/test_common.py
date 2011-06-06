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

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

