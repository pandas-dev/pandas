# -*- coding: utf-8 -*-

import numpy as np
import pytest

from pandas.core.dtypes.cast import maybe_convert_objects


@pytest.mark.parametrize("data", [[1, 2], ["apply", "banana"]])
@pytest.mark.parametrize("copy", [True, False])
def test_maybe_convert_objects_copy(data, copy):
    arr = np.array(data)
    out = maybe_convert_objects(arr, copy=copy)

    assert (arr is out) is (not copy)
