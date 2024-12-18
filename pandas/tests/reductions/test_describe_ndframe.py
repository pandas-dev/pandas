# -*- encoding: utf-8 -*-

"""
We test the describe_ndframe function.
"""

import numpy as np
import pytest

from pandas.core.methods.describe import _refine_percentiles

@pytest.mark.parametrize(
    "percentiles_, expected", [
        (None, np.array([0.25, 0.5, 0.75])),
        ([], np.array([0.5])),
        ([0.3, 0.6], np.array([0.3, 0.6])),
    ]
)
def test_refine_percentiles(percentiles_, expected):
    """
    Check the performance of the _refine_percentiles when multiple
    values are passed.
    """

    assert np.array_equal(_refine_percentiles(percentiles_), expected)
