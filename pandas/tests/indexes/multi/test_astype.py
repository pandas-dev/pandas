# -*- coding: utf-8 -*-

import numpy as np
import pytest

import pandas.util.testing as tm
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.util.testing import assert_copy


def test_astype(idx):
    expected = idx.copy()
    actual = idx.astype('O')
    assert_copy(actual.levels, expected.levels)
    assert_copy(actual.labels, expected.labels)
    assert [level.name for level in actual.levels] == list(expected.names)

    with tm.assert_raises_regex(TypeError, "^Setting.*dtype.*object"):
        idx.astype(np.dtype(int))


@pytest.mark.parametrize('ordered', [True, False])
def test_astype_category(idx, ordered):
    # GH 18630
    msg = '> 1 ndim Categorical are not supported at this time'
    with tm.assert_raises_regex(NotImplementedError, msg):
        idx.astype(CategoricalDtype(ordered=ordered))

    if ordered is False:
        # dtype='category' defaults to ordered=False, so only test once
        with tm.assert_raises_regex(NotImplementedError, msg):
            idx.astype('category')
