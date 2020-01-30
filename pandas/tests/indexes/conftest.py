import pytest

import pandas._testing as tm
from pandas.core.indexes.api import Index, MultiIndex

indices_dict = {
    "unicode": tm.makeUnicodeIndex(100),
    "string": tm.makeStringIndex(100),
    "datetime": tm.makeDateIndex(100),
    "datetime-tz": tm.makeDateIndex(100, tz="US/Pacific"),
    "period": tm.makePeriodIndex(100),
    "timedelta": tm.makeTimedeltaIndex(100),
    "int": tm.makeIntIndex(100),
    "uint": tm.makeUIntIndex(100),
    "range": tm.makeRangeIndex(100),
    "float": tm.makeFloatIndex(100),
    "bool": Index([True, False]),
    "categorical": tm.makeCategoricalIndex(100),
    "interval": tm.makeIntervalIndex(100),
    "empty": Index([]),
    "tuples": MultiIndex.from_tuples(zip(["foo", "bar", "baz"], [1, 2, 3])),
    "repeats": Index([0, 0, 1, 1, 2, 2]),
}


@pytest.fixture(params=indices_dict.keys())
def indices(request):
    # copy to avoid mutation, e.g. setting .name
    return indices_dict[request.param].copy()
