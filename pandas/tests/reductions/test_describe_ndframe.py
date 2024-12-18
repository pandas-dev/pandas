# -*- encoding: utf-8 -*-

"""
We test the describe_ndframe function.
"""

import pytest
import numpy as np

from pandas.core.methods.describe import _refine_percentiles

def test_refine_percentiles():
    """
    Check the performance of the _refine_percentiles when multiple
    values are passed.
    """

    # by default 0.25, 0.50, 0.75 is returned
    # or, when None is passed return behavior is the same
    assert _refine_percentiles() == np.array([0.25, 0.5, 0.75])
    assert _refine_percentiles(percentiles = None) == np.array([0.25, 0.5, 0.75])

    # when any value is passed, then the function should return
    percentiles_ = [0.3, 0.6]
    assert _refine_percentiles(percentiles_) == np.array(percentiles_)

    # when a blank list is passed, then should return only 0.5
    assert _refine_percentiles(percentiles = []) == np.array([0.5])
