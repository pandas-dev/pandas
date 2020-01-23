import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.indexes.api import Index, MultiIndex

indices_dict = {
    "unicode": tm.makeUnicodeIndex(100),
    "string": tm.makeStringIndex(100),
    "datetime": tm.makeDateIndex(100),
    "localized-datetime": tm.makeDateIndex(100, tz="US/Eastern"),
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
    "repeats": Index([0, 0, 1, 1, 2, 2]),
}


@pytest.fixture(params=indices_dict.keys())
def indices(request):
    # copy to avoid mutation, e.g. setting .name
    return indices_dict[request.param].copy()


def _create_series(index):
    """ Helper for the _series dict """
    data = np.random.randn(len(index))
    return pd.Series(data, index=index, name=index.name)


_series = {
    f"series-with-{i_id}-index": _create_series(i) for i_id, i in indices_dict.items()
}


def _create_narrow_series(data_dtype):
    """ Helper for the _narrow_series dict """
    index = indices_dict["int"].copy()
    size = len(index)
    if np.issubdtype(data_dtype, np.float):
        data = np.random.choice(size, size=size, replace=False)
    elif np.issubdtype(data_dtype, np.integer):
        data = np.random.randn(size)
    else:
        raise ValueError(f"Received an unexpected data_dtype: {data_dtype}")

    return pd.Series(data.astype(data_dtype), index=index, name="a")


_narrow_series = {
    "float32-series": _create_narrow_series(np.float32),
    "int8-series": _create_narrow_series(np.int8),
    "int16-series": _create_narrow_series(np.int16),
    "int32-series": _create_narrow_series(np.int32),
    "uint8-series": _create_narrow_series(np.uint8),
    "uint16-series": _create_narrow_series(np.uint16),
    "uint32-series": _create_narrow_series(np.uint32),
}

_all_objs = {**indices_dict, **_series, **_narrow_series}


@pytest.fixture(params=_all_objs.keys())
def index_or_series_obj(request):
    """
    Fixture for tests on indexes, series and series with a narrow dtype

    copy to avoid mutation, e.g. setting .name
    """
    return _all_objs[request.param].copy(deep=True)
