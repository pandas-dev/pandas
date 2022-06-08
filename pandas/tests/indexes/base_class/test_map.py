import pytest

from pandas import (
    Index,
    MultiIndex,
)
import pandas._testing as tm


@pytest.fixture
def mi():
    return MultiIndex.from_tuples(
        [("a", "1/2"), ("b", "3/4"), ("d", "5/6")], names=["1", "2"]
    )


@pytest.mark.parametrize(
    "f,exp",
    [
        (lambda x: "/".join(x), Index(["a/1/2", "b/3/4", "d/5/6"])),
        (
            lambda x: tuple(y.replace("/", "") for y in x),
            MultiIndex.from_tuples([("a", "12"), ("b", "34"), ("d", "56")]),
        ),
        (
            lambda x: (x[0], *x[1].split("/")),
            MultiIndex.from_tuples([("a", "1", "2"), ("b", "3", "4"), ("d", "5", "6")]),
        ),
    ],
)
def test_map_multiindex(mi, f, exp):
    # GH 47173

    result = mi.map(f)

    tm.assert_index_equal(result, exp)
