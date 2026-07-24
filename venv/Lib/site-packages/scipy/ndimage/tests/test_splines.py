"""Tests for spline filtering."""
import pytest

import numpy as np
from scipy._lib._array_api import (
    assert_almost_equal, xp_assert_close, make_xp_test_case,
)
from scipy import ndimage

xfail_xp_backends = pytest.mark.xfail_xp_backends


def get_spline_knot_values(order):
    """Knot values to the right of a B-spline's center."""
    knot_values = {0: [1],
                   1: [1],
                   2: [6, 1],
                   3: [4, 1],
                   4: [230, 76, 1],
                   5: [66, 26, 1]}

    return knot_values[order]


def make_spline_knot_matrix(xp, n, order, mode='mirror'):
    """Matrix to invert to find the spline coefficients."""
    knot_values = get_spline_knot_values(order)

    # NB: do computations with numpy, convert to xp as the last step only

    matrix = np.zeros((n, n))
    for diag, knot_value in enumerate(knot_values):
        indices = np.arange(diag, n)
        if diag == 0:
            matrix[indices, indices] = knot_value
        else:
            matrix[indices, indices - diag] = knot_value
            matrix[indices - diag, indices] = knot_value

    knot_values_sum = knot_values[0] + 2 * sum(knot_values[1:])

    if mode == 'mirror':
        start, step = 1, 1
    elif mode == 'reflect':
        start, step = 0, 1
    elif mode == 'grid-wrap':
        start, step = -1, -1
    else:
        raise ValueError(f'unsupported mode {mode}')

    for row in range(len(knot_values) - 1):
        for idx, knot_value in enumerate(knot_values[row + 1:]):
            matrix[row, start + step*idx] += knot_value
            matrix[-row - 1, -start - 1 - step*idx] += knot_value

    return xp.asarray(matrix / knot_values_sum)


@make_xp_test_case(ndimage.spline_filter1d)  # type:ignore[attr-defined]
@pytest.mark.parametrize('order', [0, 1, 2, 3, 4, 5])
@pytest.mark.parametrize('mode', ['mirror', 'grid-wrap', 'reflect'])
def test_spline_filter_vs_matrix_solution(order, mode, xp):
    n = 100
    eye = xp.eye(n, dtype=xp.float64)
    spline_filter_axis_0 = ndimage.spline_filter1d(eye, axis=0, order=order,
                                                   mode=mode)
    spline_filter_axis_1 = ndimage.spline_filter1d(eye, axis=1, order=order,
                                                   mode=mode)
    matrix = make_spline_knot_matrix(xp, n, order, mode=mode)
    assert_almost_equal(eye, spline_filter_axis_0 @ matrix)
    assert_almost_equal(eye, spline_filter_axis_1 @ matrix.T)


@make_xp_test_case(ndimage.spline_filter1d)  # type:ignore[attr-defined]
@xfail_xp_backends("cupy", reason="CuPy spline_filter1d has the same aliasing bug")
@pytest.mark.parametrize('order', [2, 3, 4, 5])
@pytest.mark.parametrize('n', [2, 3, 4, 5])
def test_spline_filter_reflect_small_n(order, n, xp):
    # Regression test for gh-24550: the causal reflect initialization had an
    # aliasing bug where c[0] was read back after mutation via c[n-1-i].
    # For large n the error is negligible, but for small n it is significant.
    eye = xp.eye(n, dtype=xp.float64)
    filtered = ndimage.spline_filter1d(eye, axis=0, order=order, mode='reflect')
    matrix = make_spline_knot_matrix(xp, n, order, mode='reflect')
    xp_assert_close(filtered @ matrix, eye, atol=1e-12)
