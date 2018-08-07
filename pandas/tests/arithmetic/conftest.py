# -*- coding: utf-8 -*-
import pytest

import numpy as np
import pandas as pd


@pytest.fixture(params=[1, np.array(1, dtype=np.int64)])
def one(request):
    # zero-dim integer array behaves like an integer
    return request.param


# ------------------------------------------------------------------

@pytest.fixture(params=[pd.Index, pd.Series, pd.DataFrame],
                ids=lambda x: x.__name__)
def box(request):
    """
    Several array-like containers that should have effectively identical
    behavior with respect to arithmetic operations.
    """
    return request.param


@pytest.fixture(params=[
    pd.Index,
    pd.Series,
    pytest.param(pd.DataFrame,
                 marks=pytest.mark.xfail(reason="Tries to broadcast "
                                                "incorrectly",
                                         strict=True, raises=ValueError))
], ids=lambda x: x.__name__)
def box_df_broadcast_failure(request):
    """
    Fixture equivalent to `box` but with the common failing case where
    the DataFrame operation tries to broadcast incorrectly.
    """
    return request.param
