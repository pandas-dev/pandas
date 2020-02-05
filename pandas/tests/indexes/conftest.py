import pytest
import numpy as np

import pandas._testing as tm
from pandas.core.indexes.api import Index, MultiIndex


def _gen_mi():
    # a MultiIndex used to test the general functionality of the
    # general functionality of this object

    # See Also: tests.multi.conftest.idx
    major_axis = Index(["foo", "bar", "baz", "qux"])
    minor_axis = Index(["one", "two"])

    major_codes = np.array([0, 0, 1, 2, 3, 3])
    minor_codes = np.array([0, 1, 0, 1, 0, 1])
    index_names = ["first", "second"]
    mi = MultiIndex(
        levels=[major_axis, minor_axis],
        codes=[major_codes, minor_codes],
        names=index_names,
        verify_integrity=False,
    )
    return mi


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
    "bool": tm.makeBoolIndex(2),
    "categorical": tm.makeCategoricalIndex(100),
    "interval": tm.makeIntervalIndex(100),
    "empty": Index([]),
    "tuples": MultiIndex.from_tuples(zip(["foo", "bar", "baz"], [1, 2, 3])),
    "multi": _gen_mi(),
    "repeats": Index([0, 0, 1, 1, 2, 2]),
}


@pytest.fixture(params=indices_dict.keys())
def indices(request):
    # copy to avoid mutation, e.g. setting .name
    return indices_dict[request.param].copy()
