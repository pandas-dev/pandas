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
    unicodeIndex=(tm.makeUnicodeIndex, np.dtype('O')),
    strIndex=(tm.makeStringIndex, np.dtype('O')),
    dateIndex=(tm.makeDateIndex, np.dtype('<M8[ns]')),
    periodIndex=(tm.makePeriodIndex, PeriodDtype('B')),
    tdIndex=(tm.makeTimedeltaIndex, np.dtype('<m8[ns]')),
    intIndex=(tm.makeIntIndex, np.int64),
    uintIndex=(tm.makeUIntIndex, np.uint64),
    rangeIndex=(tm.makeRangeIndex, np.int64),
    floatIndex=(tm.makeFloatIndex, np.float64),
    catIndex=(tm.makeCategoricalIndex, ),
    emptyIndex=(makeEmptyIndex, np.dtype('O')),
    intervalIndex=(tm.makeIntervalIndex, ),
)   


COMPATIBLE_INCONSISTENT_PAIRS = {
    ('intIndex', 'rangeIndex')
}


@pytest.mark.parametrize('idxType', 
    INDEXES.keys()
)
def test_union_same_types(idxType):
    idx1 = INDEXES[idxType][0](10)
    idx2 = INDEXES[idxType][0](20)
    assert idx1.union(idx2).dtype == idx1.dtype

    # Note: catIndex reflects only left dtype, should it reflect both?


@pytest.mark.parametrize('idxType1,idxType2', 
    list(it.combinations(INDEXES, 2))
)
def test_union_different_types(idxType1, idxType2):
    if idxType1 == 'catIndex' or idxType2 == 'catIndex':
        return

    if tuple(sorted([idxType1, idxType2])) in COMPATIBLE_INCONSISTENT_PAIRS:
        return    

    idx1 = INDEXES[idxType1][0](10)
    idx2 = INDEXES[idxType2][0](20)
    assert idx1.union(idx2).dtype == np.dtype('O')

    assert idx2.union(idx1).dtype == np.dtype('O')
