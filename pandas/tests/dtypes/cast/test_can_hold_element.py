import numpy as np

from pandas.core.dtypes.cast import can_hold_element


def test_can_hold_element_range(any_int_numpy_dtype):
    dtype = np.dtype(any_int_numpy_dtype)
    arr = np.array([], dtype=dtype)

    rng = range(2, 127)
    assert can_hold_element(arr, rng)

    # negatives -> can't be held by uint dtypes
    rng = range(-2, 127)
    if dtype.kind == "i":
        assert can_hold_element(arr, rng)
    else:
        assert not can_hold_element(arr, rng)

    rng = range(2, 255)
    if dtype == "int8":
        assert not can_hold_element(arr, rng)
    else:
        assert can_hold_element(arr, rng)

    rng = range(-255, 65537)
    if dtype.kind == "u":
        assert not can_hold_element(arr, rng)
    elif dtype.itemsize < 4:
        assert not can_hold_element(arr, rng)
    else:
        assert can_hold_element(arr, rng)

    # empty
    rng = range(-(10 ** 10), -(10 ** 10))
    assert len(rng) == 0
    # assert can_hold_element(arr, rng)

    rng = range(10 ** 10, 10 ** 10)
    assert len(rng) == 0
    assert can_hold_element(arr, rng)
