import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Index
import pandas._testing as tm


@pytest.mark.parametrize(
    "interpolation", ["linear", "lower", "higher", "nearest", "midpoint"]
)
@pytest.mark.parametrize(
    "a_vals,b_vals",
    [
        # Ints
        ([1, 2, 3, 4, 5], [5, 4, 3, 2, 1]),
        ([1, 2, 3, 4], [4, 3, 2, 1]),
        ([1, 2, 3, 4, 5], [4, 3, 2, 1]),
        # Floats
        ([1.0, 2.0, 3.0, 4.0, 5.0], [5.0, 4.0, 3.0, 2.0, 1.0]),
        # Missing data
        ([1.0, np.nan, 3.0, np.nan, 5.0], [5.0, np.nan, 3.0, np.nan, 1.0]),
        ([np.nan, 4.0, np.nan, 2.0, np.nan], [np.nan, 4.0, np.nan, 2.0, np.nan]),
        # Timestamps
        (
            list(pd.date_range("1/1/18", freq="D", periods=5)),
            list(pd.date_range("1/1/18", freq="D", periods=5))[::-1],
        ),
        # All NA
        ([np.nan] * 5, [np.nan] * 5),
    ],
)
@pytest.mark.parametrize("q", [0, 0.25, 0.5, 0.75, 1])
def test_quantile(interpolation, a_vals, b_vals, q):
    if interpolation == "nearest" and q == 0.5 and b_vals == [4, 3, 2, 1]:
        pytest.skip(
            "Unclear numpy expectation for nearest result with equidistant data"
        )

    a_expected = pd.Series(a_vals).quantile(q, interpolation=interpolation)
    b_expected = pd.Series(b_vals).quantile(q, interpolation=interpolation)

    df = DataFrame(
        {"key": ["a"] * len(a_vals) + ["b"] * len(b_vals), "val": a_vals + b_vals}
    )

    expected = DataFrame(
        [a_expected, b_expected], columns=["val"], index=Index(["a", "b"], name="key")
    )
    result = df.groupby("key").quantile(q, interpolation=interpolation)

    tm.assert_frame_equal(result, expected)


def test_quantile_array():
    # https://github.com/pandas-dev/pandas/issues/27526
    df = pd.DataFrame({"A": [0, 1, 2, 3, 4]})
    result = df.groupby([0, 0, 1, 1, 1]).quantile([0.25])

    index = pd.MultiIndex.from_product([[0, 1], [0.25]])
    expected = pd.DataFrame({"A": [0.25, 2.50]}, index=index)
    tm.assert_frame_equal(result, expected)

    df = pd.DataFrame({"A": [0, 1, 2, 3], "B": [4, 5, 6, 7]})
    index = pd.MultiIndex.from_product([[0, 1], [0.25, 0.75]])

    result = df.groupby([0, 0, 1, 1]).quantile([0.25, 0.75])
    expected = pd.DataFrame(
        {"A": [0.25, 0.75, 2.25, 2.75], "B": [4.25, 4.75, 6.25, 6.75]}, index=index
    )
    tm.assert_frame_equal(result, expected)


def test_quantile_array2():
    # https://github.com/pandas-dev/pandas/pull/28085#issuecomment-524066959
    df = pd.DataFrame(
        np.random.RandomState(0).randint(0, 5, size=(10, 3)), columns=list("ABC")
    )
    result = df.groupby("A").quantile([0.3, 0.7])
    expected = pd.DataFrame(
        {
            "B": [0.9, 2.1, 2.2, 3.4, 1.6, 2.4, 2.3, 2.7, 0.0, 0.0],
            "C": [1.2, 2.8, 1.8, 3.0, 0.0, 0.0, 1.9, 3.1, 3.0, 3.0],
        },
        index=pd.MultiIndex.from_product(
            [[0, 1, 2, 3, 4], [0.3, 0.7]], names=["A", None]
        ),
    )
    tm.assert_frame_equal(result, expected)


def test_quantile_array_no_sort():
    df = pd.DataFrame({"A": [0, 1, 2], "B": [3, 4, 5]})
    result = df.groupby([1, 0, 1], sort=False).quantile([0.25, 0.5, 0.75])
    expected = pd.DataFrame(
        {"A": [0.5, 1.0, 1.5, 1.0, 1.0, 1.0], "B": [3.5, 4.0, 4.5, 4.0, 4.0, 4.0]},
        index=pd.MultiIndex.from_product([[1, 0], [0.25, 0.5, 0.75]]),
    )
    tm.assert_frame_equal(result, expected)

    result = df.groupby([1, 0, 1], sort=False).quantile([0.75, 0.25])
    expected = pd.DataFrame(
        {"A": [1.5, 0.5, 1.0, 1.0], "B": [4.5, 3.5, 4.0, 4.0]},
        index=pd.MultiIndex.from_product([[1, 0], [0.75, 0.25]]),
    )
    tm.assert_frame_equal(result, expected)


def test_quantile_array_multiple_levels():
    df = pd.DataFrame(
        {"A": [0, 1, 2], "B": [3, 4, 5], "c": ["a", "a", "a"], "d": ["a", "a", "b"]}
    )
    result = df.groupby(["c", "d"]).quantile([0.25, 0.75])
    index = pd.MultiIndex.from_tuples(
        [("a", "a", 0.25), ("a", "a", 0.75), ("a", "b", 0.25), ("a", "b", 0.75)],
        names=["c", "d", None],
    )
    expected = pd.DataFrame(
        {"A": [0.25, 0.75, 2.0, 2.0], "B": [3.25, 3.75, 5.0, 5.0]}, index=index
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("frame_size", [(2, 3), (100, 10)])
@pytest.mark.parametrize("groupby", [[0], [0, 1]])
@pytest.mark.parametrize("q", [[0.5, 0.6]])
def test_groupby_quantile_with_arraylike_q_and_int_columns(frame_size, groupby, q):
    # GH30289
    nrow, ncol = frame_size
    df = pd.DataFrame(
        np.array([ncol * [_ % 4] for _ in range(nrow)]), columns=range(ncol)
    )

    idx_levels = [list(range(min(nrow, 4)))] * len(groupby) + [q]
    idx_codes = [[x for x in range(min(nrow, 4)) for _ in q]] * len(groupby) + [
        list(range(len(q))) * min(nrow, 4)
    ]
    expected_index = pd.MultiIndex(
        levels=idx_levels, codes=idx_codes, names=groupby + [None]
    )
    expected_values = [
        [float(x)] * (ncol - len(groupby)) for x in range(min(nrow, 4)) for _ in q
    ]
    expected_columns = [x for x in range(ncol) if x not in groupby]
    expected = pd.DataFrame(
        expected_values, index=expected_index, columns=expected_columns
    )
    result = df.groupby(groupby).quantile(q)

    tm.assert_frame_equal(result, expected)


def test_quantile_raises():
    df = pd.DataFrame(
        [["foo", "a"], ["foo", "b"], ["foo", "c"]], columns=["key", "val"]
    )

    with pytest.raises(TypeError, match="cannot be performed against 'object' dtypes"):
        df.groupby("key").quantile()


def test_quantile_out_of_bounds_q_raises():
    # https://github.com/pandas-dev/pandas/issues/27470
    df = pd.DataFrame(dict(a=[0, 0, 0, 1, 1, 1], b=range(6)))
    g = df.groupby([0, 0, 0, 1, 1, 1])
    with pytest.raises(ValueError, match="Got '50.0' instead"):
        g.quantile(50)

    with pytest.raises(ValueError, match="Got '-1.0' instead"):
        g.quantile(-1)


def test_quantile_missing_group_values_no_segfaults():
    # GH 28662
    data = np.array([1.0, np.nan, 1.0])
    df = pd.DataFrame(dict(key=data, val=range(3)))

    # Random segfaults; would have been guaranteed in loop
    grp = df.groupby("key")
    for _ in range(100):
        grp.quantile()


@pytest.mark.parametrize(
    "key, val, expected_key, expected_val",
    [
        ([1.0, np.nan, 3.0, np.nan], range(4), [1.0, 3.0], [0.0, 2.0]),
        ([1.0, np.nan, 2.0, 2.0], range(4), [1.0, 2.0], [0.0, 2.5]),
        (["a", "b", "b", np.nan], range(4), ["a", "b"], [0, 1.5]),
        ([0], [42], [0], [42.0]),
        ([], [], np.array([], dtype="float64"), np.array([], dtype="float64")),
    ],
)
def test_quantile_missing_group_values_correct_results(
    key, val, expected_key, expected_val
):
    # GH 28662, GH 33200, GH 33569
    df = pd.DataFrame({"key": key, "val": val})

    expected = pd.DataFrame(
        expected_val, index=pd.Index(expected_key, name="key"), columns=["val"]
    )

    grp = df.groupby("key")

    result = grp.quantile(0.5)
    tm.assert_frame_equal(result, expected)

    result = grp.quantile()
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "values",
    [
        pd.array([1, 0, None] * 2, dtype="Int64"),
        pd.array([True, False, None] * 2, dtype="boolean"),
    ],
)
@pytest.mark.parametrize("q", [0.5, [0.0, 0.5, 1.0]])
def test_groupby_quantile_nullable_array(values, q):
    # https://github.com/pandas-dev/pandas/issues/33136
    df = pd.DataFrame({"a": ["x"] * 3 + ["y"] * 3, "b": values})
    result = df.groupby("a")["b"].quantile(q)

    if isinstance(q, list):
        idx = pd.MultiIndex.from_product((["x", "y"], q), names=["a", None])
        true_quantiles = [0.0, 0.5, 1.0]
    else:
        idx = pd.Index(["x", "y"], name="a")
        true_quantiles = [0.5]

    expected = pd.Series(true_quantiles * 2, index=idx, name="b")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("q", [0.5, [0.0, 0.5, 1.0]])
def test_groupby_quantile_skips_invalid_dtype(q):
    df = pd.DataFrame({"a": [1], "b": [2.0], "c": ["x"]})
    result = df.groupby("a").quantile(q)
    expected = df.groupby("a")[["b"]].quantile(q)
    tm.assert_frame_equal(result, expected)
