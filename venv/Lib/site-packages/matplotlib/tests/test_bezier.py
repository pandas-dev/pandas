"""
Tests specific to the bezier module.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from matplotlib.bezier import (
    _real_roots_in_01, inside_circle, split_bezier_intersecting_with_closedpath
)


@pytest.mark.parametrize("roots, expected_in_01", [
    ([0.5], [0.5]),
    ([0.25, 0.75], [0.25, 0.75]),
    ([0.2, 0.5, 0.8], [0.2, 0.5, 0.8]),
    ([0.1, 0.2, 0.3, 0.4], [0.1, 0.2, 0.3, 0.4]),
    ([0.0, 0.5], [0.0, 0.5]),
    ([0.5, 1.0], [0.5, 1.0]),
    ([2.0], []),                 # outside [0, 1]
    ([0.5, 2.0], [0.5]),         # one in, one out
    ([-1j, 1j], []),             # complex roots
    ([0.5, -1j, 1j], [0.5]),     # mix of real and complex
    ([0.3, 0.3], [0.3, 0.3]),    # repeated root
])
def test_real_roots_in_01(roots, expected_in_01):
    roots = np.array(roots)
    coeffs = np.poly(roots)[::-1]  # np.poly gives descending, we need ascending
    result = _real_roots_in_01(coeffs.real)
    assert_allclose(result, expected_in_01, atol=1e-10)


@pytest.mark.parametrize("coeffs", [[5], [0, 0, 0]])
def test_real_roots_in_01_no_roots(coeffs):
    assert len(_real_roots_in_01(coeffs)) == 0


def test_split_bezier_with_large_values():
    # These numbers come from gh-27753
    arrow_path = [(96950809781500.0, 804.7503795623779),
                  (96950809781500.0, 859.6242585800646),
                  (96950809781500.0, 914.4981375977513)]
    in_f = inside_circle(96950809781500.0, 804.7503795623779, 0.06)
    split_bezier_intersecting_with_closedpath(arrow_path, in_f)
    # All we are testing is that this completes
    # The failure case is an infinite loop resulting from floating point precision
    # pytest will timeout if that occurs
