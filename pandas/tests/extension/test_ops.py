import pytest

from pandas.tests.extension.decimal.array import to_decimal
import pandas.util.testing as tm


@pytest.mark.parametrize("reverse, expected_div, expected_mod", [
    (False, [0, 1, 1, 2], [1, 0, 1, 0]),
    (True, [2, 1, 0, 0], [0, 0, 2, 2]),
])
def test_divmod(reverse, expected_div, expected_mod):
    # https://github.com/pandas-dev/pandas/issues/22930
    arr = to_decimal([1, 2, 3, 4])
    if reverse:
        div, mod = divmod(2, arr)
    else:
        div, mod = divmod(arr, 2)
    expected_div = to_decimal(expected_div)
    expected_mod = to_decimal(expected_mod)

    tm.assert_extension_array_equal(div, expected_div)
    tm.assert_extension_array_equal(mod, expected_mod)
