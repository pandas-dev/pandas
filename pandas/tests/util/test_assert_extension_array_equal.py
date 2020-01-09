import numpy as np
import pytest

import pandas._testing as tm
from pandas.core.arrays.sparse import SparseArray


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(),  # Default is check_exact=False
        dict(check_exact=False),
        dict(check_exact=True),
    ],
)
def test_assert_extension_array_equal_not_exact(kwargs):
    # see gh-23709
    arr1 = SparseArray([-0.17387645482451206, 0.3414148016424936])
    arr2 = SparseArray([-0.17387645482451206, 0.3414148016424937])

    if kwargs.get("check_exact", False):
        msg = """\
ExtensionArray are different

ExtensionArray values are different \\(50\\.0 %\\)
\\[left\\]:  \\[-0\\.17387645482.*, 0\\.341414801642.*\\]
\\[right\\]: \\[-0\\.17387645482.*, 0\\.341414801642.*\\]"""

        with pytest.raises(AssertionError, match=msg):
            tm.assert_extension_array_equal(arr1, arr2, **kwargs)
    else:
        tm.assert_extension_array_equal(arr1, arr2, **kwargs)


@pytest.mark.parametrize(
    "check_less_precise", [True, False, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
)
def test_assert_extension_array_equal_less_precise(check_less_precise):
    arr1 = SparseArray([0.5, 0.123456])
    arr2 = SparseArray([0.5, 0.123457])

    kwargs = dict(check_less_precise=check_less_precise)

    if check_less_precise is False or check_less_precise >= 5:
        msg = """\
ExtensionArray are different

ExtensionArray values are different \\(50\\.0 %\\)
\\[left\\]:  \\[0\\.5, 0\\.123456\\]
\\[right\\]: \\[0\\.5, 0\\.123457\\]"""

        with pytest.raises(AssertionError, match=msg):
            tm.assert_extension_array_equal(arr1, arr2, **kwargs)
    else:
        tm.assert_extension_array_equal(arr1, arr2, **kwargs)


def test_assert_extension_array_equal_dtype_mismatch(check_dtype):
    end = 5
    kwargs = dict(check_dtype=check_dtype)

    arr1 = SparseArray(np.arange(end, dtype="int64"))
    arr2 = SparseArray(np.arange(end, dtype="int32"))

    if check_dtype:
        msg = """\
ExtensionArray are different

Attribute "dtype" are different
\\[left\\]:  Sparse\\[int64, 0\\]
\\[right\\]: Sparse\\[int32, 0\\]"""

        with pytest.raises(AssertionError, match=msg):
            tm.assert_extension_array_equal(arr1, arr2, **kwargs)
    else:
        tm.assert_extension_array_equal(arr1, arr2, **kwargs)


def test_assert_extension_array_equal_missing_values():
    arr1 = SparseArray([np.nan, 1, 2, np.nan])
    arr2 = SparseArray([np.nan, 1, 2, 3])

    msg = """\
ExtensionArray NA mask are different

ExtensionArray NA mask values are different \\(25\\.0 %\\)
\\[left\\]:  \\[True, False, False, True\\]
\\[right\\]: \\[True, False, False, False\\]"""

    with pytest.raises(AssertionError, match=msg):
        tm.assert_extension_array_equal(arr1, arr2)


@pytest.mark.parametrize("side", ["left", "right"])
def test_assert_extension_array_equal_non_extension_array(side):
    numpy_array = np.arange(5)
    extension_array = SparseArray(numpy_array)

    msg = f"{side} is not an ExtensionArray"
    args = (
        (numpy_array, extension_array)
        if side == "left"
        else (extension_array, numpy_array)
    )

    with pytest.raises(AssertionError, match=msg):
        tm.assert_extension_array_equal(*args)
