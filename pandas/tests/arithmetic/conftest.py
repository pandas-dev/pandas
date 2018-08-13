# -*- coding: utf-8 -*-
import pytest

import numpy as np
import pandas as pd

from pandas.compat import long


@pytest.fixture(params=[1, np.array(1, dtype=np.int64)])
def one(request):
    # zero-dim integer array behaves like an integer
    return request.param


zeros = [box_cls([0] * 5, dtype=dtype)
         for box_cls in [pd.Index, np.array]
         for dtype in [np.int64, np.uint64, np.float64]]
zeros.extend([np.array(0, dtype=dtype)
              for dtype in [np.int64, np.uint64, np.float64]])
zeros.extend([0, 0.0, long(0)])


@pytest.fixture(params=zeros)
def zero(request):
    # For testing division by (or of) zero for Index with length 5, this
    # gives several scalar-zeros and length-5 vector-zeros
    return request.param


@pytest.fixture(params=[pd.Float64Index(np.arange(5, dtype='float64')),
                        pd.Int64Index(np.arange(5, dtype='int64')),
                        pd.UInt64Index(np.arange(5, dtype='uint64'))],
                ids=lambda x: type(x).__name__)
def idx(request):
    return request.param


@pytest.fixture(params=[pd.Timedelta('5m4s').to_pytimedelta(),
                        pd.Timedelta('5m4s'),
                        pd.Timedelta('5m4s').to_timedelta64()],
                ids=lambda x: type(x).__name__)
def scalar_td(request):
    """
    Several variants of Timedelta scalars representing 5 minutes and 4 seconds
    """
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
