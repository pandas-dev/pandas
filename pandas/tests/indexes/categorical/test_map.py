import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "data, categories",
    [
        (list("abcbca"), list("cab")),
        (pd.interval_range(0, 3).repeat(3), pd.interval_range(0, 3)),
    ],
    ids=["string", "interval"],
)
def test_map_str(data, categories, ordered):
    # GH 31202 - override base class since we want to maintain categorical/ordered
    index = pd.CategoricalIndex(data, categories=categories, ordered=ordered)
    result = index.map(str)
    expected = pd.CategoricalIndex(
        map(str, data), categories=map(str, categories), ordered=ordered
    )
    tm.assert_index_equal(result, expected)


def test_map():
    ci = pd.CategoricalIndex(list("ABABC"), categories=list("CBA"), ordered=True)
    result = ci.map(lambda x: x.lower())
    exp = pd.CategoricalIndex(list("ababc"), categories=list("cba"), ordered=True)
    tm.assert_index_equal(result, exp)

    ci = pd.CategoricalIndex(
        list("ABABC"), categories=list("BAC"), ordered=False, name="XXX"
    )
    result = ci.map(lambda x: x.lower())
    exp = pd.CategoricalIndex(
        list("ababc"), categories=list("bac"), ordered=False, name="XXX"
    )
    tm.assert_index_equal(result, exp)

    # GH 12766: Return an index not an array
    tm.assert_index_equal(
        ci.map(lambda x: 1), pd.Index(np.array([1] * 5, dtype=np.int64), name="XXX")
    )

    # change categories dtype
    ci = pd.CategoricalIndex(list("ABABC"), categories=list("BAC"), ordered=False)

    def f(x):
        return {"A": 10, "B": 20, "C": 30}.get(x)

    result = ci.map(f)
    exp = pd.CategoricalIndex(
        [10, 20, 10, 20, 30], categories=[20, 10, 30], ordered=False
    )
    tm.assert_index_equal(result, exp)

    result = ci.map(pd.Series([10, 20, 30], index=["A", "B", "C"]))
    tm.assert_index_equal(result, exp)

    result = ci.map({"A": 10, "B": 20, "C": 30})
    tm.assert_index_equal(result, exp)


def test_map_with_categorical_series():
    # GH 12756
    a = pd.Index([1, 2, 3, 4])
    b = pd.Series(["even", "odd", "even", "odd"], dtype="category")
    c = pd.Series(["even", "odd", "even", "odd"])

    exp = pd.CategoricalIndex(["odd", "even", "odd", np.nan])
    tm.assert_index_equal(a.map(b), exp)
    exp = pd.Index(["odd", "even", "odd", np.nan])
    tm.assert_index_equal(a.map(c), exp)


@pytest.mark.parametrize(
    ("data", "f", "expected"),
    (
        ([1, 1, np.nan], pd.isna, pd.CategoricalIndex([False, False, np.nan])),
        ([1, 2, np.nan], pd.isna, pd.Index([False, False, np.nan])),
        ([1, 1, np.nan], {1: False}, pd.CategoricalIndex([False, False, np.nan])),
        ([1, 2, np.nan], {1: False, 2: False}, pd.Index([False, False, np.nan])),
        (
            [1, 1, np.nan],
            pd.Series([False, False]),
            pd.CategoricalIndex([False, False, np.nan]),
        ),
        (
            [1, 2, np.nan],
            pd.Series([False, False, False]),
            pd.Index([False, False, np.nan]),
        ),
    ),
)
def test_map_with_nan_ignore(data, f, expected):  # GH 24241
    values = pd.CategoricalIndex(data)
    result = values.map(f, na_action="ignore")
    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize(
    ("data", "f", "expected"),
    (
        ([1, 1, np.nan], pd.isna, pd.Index([False, False, True])),
        ([1, 2, np.nan], pd.isna, pd.Index([False, False, True])),
        ([1, 1, np.nan], {1: False}, pd.CategoricalIndex([False, False, np.nan])),
        ([1, 2, np.nan], {1: False, 2: False}, pd.Index([False, False, np.nan])),
        (
            [1, 1, np.nan],
            pd.Series([False, False]),
            pd.CategoricalIndex([False, False, np.nan]),
        ),
        (
            [1, 2, np.nan],
            pd.Series([False, False, False]),
            pd.Index([False, False, np.nan]),
        ),
    ),
)
def test_map_with_nan_none(data, f, expected):  # GH 24241
    values = pd.CategoricalIndex(data)
    result = values.map(f, na_action=None)
    tm.assert_index_equal(result, expected)


def test_map_with_dict_or_series():
    orig_values = ["a", "B", 1, "a"]
    new_values = ["one", 2, 3.0, "one"]
    cur_index = pd.CategoricalIndex(orig_values, name="XXX")
    expected = pd.CategoricalIndex(new_values, name="XXX", categories=[3.0, 2, "one"])

    mapper = pd.Series(new_values[:-1], index=orig_values[:-1])
    result = cur_index.map(mapper)
    # Order of categories in result can be different
    tm.assert_index_equal(result, expected)

    mapper = dict(zip(orig_values[:-1], new_values[:-1], strict=True))
    result = cur_index.map(mapper)
    # Order of categories in result can be different
    tm.assert_index_equal(result, expected)
