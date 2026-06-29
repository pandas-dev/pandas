import numpy as np

import pandas._testing as tm

from pandas.core.array_algos.putmask import putmask_without_repeat


def test_putmask_without_repeat_scalar_value():
    values = np.array([[1.0, np.nan], [2.0, 3.0]])
    mask = np.array([[True, False], [False, True]])

    putmask_without_repeat(values, mask, np.float64(0))

    expected = np.array([[0.0, np.nan], [2.0, 0.0]])
    tm.assert_numpy_array_equal(values, expected)


def test_putmask_without_repeat_listlike_exact_length():
    values = np.array([1.0, 2.0, 3.0, 4.0])
    mask = np.array([True, False, True, False])

    putmask_without_repeat(values, mask, np.array([10.0, 20.0]))

    expected = np.array([10.0, 2.0, 20.0, 4.0])
    tm.assert_numpy_array_equal(values, expected)
