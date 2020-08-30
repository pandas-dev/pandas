import numpy as np
import pytest

import pandas._testing as tm


def test_transpose(index_or_series_obj):
    obj = index_or_series_obj
    tm.assert_equal(obj.transpose(), obj)


def test_transpose_non_default_axes(index_or_series_obj):
    msg = "the 'axes' parameter is not supported"
    obj = index_or_series_obj
    with pytest.raises(ValueError, match=msg):
        obj.transpose(1)
    with pytest.raises(ValueError, match=msg):
        obj.transpose(axes=1)


def test_numpy_transpose(index_or_series_obj):
    msg = "the 'axes' parameter is not supported"
    obj = index_or_series_obj
    tm.assert_equal(np.transpose(obj), obj)

    with pytest.raises(ValueError, match=msg):
        np.transpose(obj, axes=1)
