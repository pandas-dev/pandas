'''
The tests in this package are to ensure the proper resultant dtypes of
set operations.
'''
import itertools as it
import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas.core.dtypes.dtypes import PeriodDtype, CategoricalDtype, \
                                      IntervalDtype


def makeEmptyIndex(_=None):
    return pd.Index([])


INDEXES = dict(
    unicodeIndex=tm.makeUnicodeIndex,
    strIndex=tm.makeStringIndex,
    dateIndex=tm.makeDateIndex,
    periodIndex=tm.makePeriodIndex,
    tdIndex=tm.makeTimedeltaIndex,
    intIndex=tm.makeIntIndex,
    uintIndex=tm.makeUIntIndex,
    rangeIndex=tm.makeRangeIndex,
    floatIndex=tm.makeFloatIndex,
    catIndex=tm.makeCategoricalIndex,
    emptyIndex=makeEmptyIndex,
    intervalIndex=tm.makeIntervalIndex,
)   


COMPATIBLE_INCONSISTENT_PAIRS = {
    ('intIndex', 'rangeIndex')
}


@pytest.mark.parametrize('idxType', 
    INDEXES.keys()
)
def test_union_same_types(idxType):
    idx1 = INDEXES[idxType](10)
    idx2 = INDEXES[idxType](20)
    assert idx1.union(idx2).dtype == idx1.dtype

    # Note: catIndex reflects only left dtype, should it reflect both?


@pytest.mark.parametrize('idxType1,idxType2', 
    list(it.combinations(INDEXES, 2))
)
def test_union_different_types(idxType1, idxType2):
    if tuple(sorted([idxType1, idxType2])) in COMPATIBLE_INCONSISTENT_PAIRS:
        return    

    idx1 = INDEXES[idxType1](10)
    idx2 = INDEXES[idxType2](20)

    # A union with a CategoricalIndex (even as dtype('O')) and a 
    # non-CategoricalIndex can only be made if both indices are monotonic.
    # This is true before this PR as well.
    if idxType1 == 'catIndex' or idxType2 == 'catIndex':
        idx1 = idx1.sort_values()
        idx2 = idx2.sort_values()

    assert idx1.union(idx2).dtype == np.dtype('O')
    assert idx2.union(idx1).dtype == np.dtype('O')


@pytest.mark.parametrize('idxType1,idxType2',
    COMPATIBLE_INCONSISTENT_PAIRS
)
def test_compatible_inconsistent_pairs(idxType1, idxType2):
    idx1 = INDEXES[idxType1](10)
    idx2 = INDEXES[idxType2](20)

    res1 = idx1.union(idx2)
    res2 = idx2.union(idx1)

    assert res1.dtype in (idx1.dtype, idx2.dtype)
    assert res2.dtype in (idx1.dtype, idx2.dtype)
