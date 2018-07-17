# -*- coding: utf-8 -*-

import numpy as np
import pandas.util.testing as tm
from pandas.util.testing import assert_copy


def test_astype(idx):
    expected = idx.copy()
    actual = idx.astype('O')
    assert_copy(actual.levels, expected.levels)
    assert_copy(actual.labels, expected.labels)
    assert [level.name for level in actual.levels] == list(expected.names)

    with tm.assert_raises_regex(TypeError, "^Setting.*dtype.*object"):
        idx.astype(np.dtype(int))
