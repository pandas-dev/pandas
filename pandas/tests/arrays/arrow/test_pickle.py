from __future__ import annotations

import pickle

import pytest

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")

from pandas.core.arrays.arrow.array import ArrowExtensionArray


def test_pickle_roundtrip(data):
    # GH 42600
    expected = pd.Series(data)
    expected_sliced = expected.head(2)
    full_pickled = pickle.dumps(expected)
    sliced_pickled = pickle.dumps(expected_sliced)

    assert len(full_pickled) > len(sliced_pickled)

    result = pickle.loads(full_pickled)
    tm.assert_series_equal(result, expected)

    result_sliced = pickle.loads(sliced_pickled)
    tm.assert_series_equal(result_sliced, expected_sliced)


# _data was renamed to _pa_data
class OldArrowExtensionArray(ArrowExtensionArray):
    def __getstate__(self):
        state = super().__getstate__()
        state["_data"] = state.pop("_pa_array")
        return state


def test_pickle_old_arrowextensionarray():
    data = pa.array([1])
    expected = OldArrowExtensionArray(data)
    result = pickle.loads(pickle.dumps(expected))
    tm.assert_extension_array_equal(result, expected)
    assert result._pa_array == pa.chunked_array(data)
    assert not hasattr(result, "_data")
