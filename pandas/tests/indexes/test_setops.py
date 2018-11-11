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
    list(it.combinations([x for x in INDEXES if x != 'catIndex'], 2))
)
def test_union_different_types(idxType1, idxType2):
    if tuple(sorted([idxType1, idxType2])) in COMPATIBLE_INCONSISTENT_PAIRS:
        return    

    idx1 = INDEXES[idxType1](10)
    idx2 = INDEXES[idxType2](20)
    assert idx1.union(idx2).dtype == np.dtype('O')
    assert idx2.union(idx1).dtype == np.dtype('O')
