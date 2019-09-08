import numpy as np
import pytest

from pandas.core.dtypes.cast import maybe_convert_objects


@pytest.mark.parametrize("data", [[1, 2], ["apply", "banana"]])
def test_maybe_convert_objects_copy(data):
    arr = np.array(data)
    out = maybe_convert_objects(arr)

    assert arr is not out
