import pytest

from pandas import DataFrame
import pandas._testing as tm


def test_methods_iloc_warn():
    df = DataFrame({"a": [1, 2, 3], "b": 1})
    with tm.assert_cow_warning(match="A value"):
        df.iloc[:, 0].replace(1, 5, inplace=True)

    with tm.assert_cow_warning(match="A value"):
        df.iloc[:, 0].fillna(1, inplace=True)

    with tm.assert_cow_warning(match="A value"):
        df.iloc[:, 0].interpolate(inplace=True)

    with tm.assert_cow_warning(match="A value"):
        df.iloc[:, 0].ffill(inplace=True)

    with tm.assert_cow_warning(match="A value"):
        df.iloc[:, 0].bfill(inplace=True)


@pytest.mark.parametrize(
    "func, args",
    [
        ("replace", (1, 5)),
        ("fillna", (1,)),
        ("interpolate", ()),
        ("bfill", ()),
        ("ffill", ()),
    ],
)
def test_methods_iloc_getitem_dont_warn(func, args):
    df = DataFrame({"a": [1, 2, 3], "b": 1})
    ser = df.iloc[:, 0]
    getattr(ser, func)(*args, inplace=True)

    # parent that holds item_cache is dead, so don't increase ref count
    ser = df.copy()["a"]
    getattr(ser, func)(*args, inplace=True)
