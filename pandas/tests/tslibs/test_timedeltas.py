# -*- coding: utf-8 -*-
import numpy as np
import pytest

from pandas._libs.tslibs.timedeltas import delta_to_nanoseconds

import pandas as pd


def test_delta_to_nanoseconds():
    obj = np.timedelta64(14, 'D')
    result = delta_to_nanoseconds(obj)
    assert result == 14 * 24 * 3600 * 1e9

    obj = pd.Timedelta(minutes=-7)
    result = delta_to_nanoseconds(obj)
    assert result == -7 * 60 * 1e9

    obj = pd.Timedelta(minutes=-7).to_pytimedelta()
    result = delta_to_nanoseconds(obj)
    assert result == -7 * 60 * 1e9

    obj = pd.offsets.Nano(125)
    result = delta_to_nanoseconds(obj)
    assert result == 125

    obj = 1
    result = delta_to_nanoseconds(obj)
    assert obj == 1

    obj = np.int64(2)
    result = delta_to_nanoseconds(obj)
    assert obj == 2

    obj = np.int32(3)
    result = delta_to_nanoseconds(obj)
    assert result == 3

    obj = np.array([123456789], dtype='m8[ns]')
    with pytest.raises(TypeError):
        delta_to_nanoseconds(obj)
