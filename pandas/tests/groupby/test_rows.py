""" Test GroupBy.rows positional grouped indexing GH#42864"""

import random

import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture()
def small_df():
    data = [
        [0, "a", "a0_at_0"],
        [1, "b", "b0_at_1"],
        [2, "a", "a1_at_2"],
        [3, "b", "b1_at_3"],
        [4, "c", "c0_at_4"],
        [5, "a", "a2_at_5"],
        [6, "a", "a3_at_6"],
        [7, "a", "a4_at_7"],
    ]
    df = pd.DataFrame(data, columns=["Index", "Category", "Value"])
    return df.set_index("Index")


@pytest.mark.parametrize(
    "arg, expected_rows",
    [
        [0, [0, 1, 4]],
        [2, [5]],
        [5, []],
        [-1, [3, 4, 7]],
        [-2, [1, 6]],
        [-6, []],
    ],
)
def test_int(small_df, arg, expected_rows):
    """Test single integer"""

    result = small_df.groupby("Category").rows[arg]
    expected = small_df.iloc[expected_rows]

    tm.assert_frame_equal(result, expected)


def test_slice(small_df):
    """Test single slice"""

    result = small_df.groupby("Category").rows[0:3:2]
    expected = small_df.iloc[[0, 1, 4, 5]]

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "arg, expected_rows",
    [
        [[0, 2], [0, 1, 4, 5]],
        [[0, 2, -1], [0, 1, 3, 4, 5, 7]],
        [range(0, 3, 2), [0, 1, 4, 5]],
        [{0, 2}, [0, 1, 4, 5]],
    ],
    ids=[
        "list",
        "negative",
        "range",
        "set",
    ],
)
def test_list(small_df, arg, expected_rows):
    """Test lists of integers and integer valued iterables"""

    result = small_df.groupby("Category").rows[arg]
    expected = small_df.iloc[expected_rows]

    tm.assert_frame_equal(result, expected)


def test_ints(small_df):
    """Test tuple of ints"""

    result = small_df.groupby("Category").rows[0, 2, -1]
    expected = small_df.iloc[[0, 1, 3, 4, 5, 7]]

    tm.assert_frame_equal(result, expected)


def test_slices(small_df):
    """Test tuple of slices"""

    result = small_df.groupby("Category").rows[:2, -2:]
    expected = small_df.iloc[[0, 1, 2, 3, 4, 6, 7]]

    tm.assert_frame_equal(result, expected)


def test_mix(small_df):
    """Test mixed tuple of ints and slices"""

    result = small_df.groupby("Category").rows[0, 1, -2:]
    expected = small_df.iloc[[0, 1, 2, 3, 4, 6, 7]]

    tm.assert_frame_equal(result, expected)


def test_doc_examples():
    """Test the examples in the documentation"""

    df = pd.DataFrame(
        [["a", 1], ["a", 2], ["a", 3], ["b", 4], ["b", 5]], columns=["A", "B"]
    )

    grouped = df.groupby("A")

    result = grouped.rows[1:2]
    expected = pd.DataFrame([["a", 2], ["b", 5]], columns=["A", "B"], index=[1, 4])

    tm.assert_frame_equal(result, expected)

    result = grouped.rows[1, -1]
    expected = pd.DataFrame(
        [["a", 2], ["a", 3], ["b", 5]], columns=["A", "B"], index=[1, 2, 4]
    )

    tm.assert_frame_equal(result, expected)


def test_multiindex():
    """Test the multiindex mentioned as the use-case in the documentation"""

    def make_df_from_data(data):
        rows = {}
        for date in dates:
            for level in data[date]:
                rows[(date, level[0])] = {"A": level[1], "B": level[2]}

        df = pd.DataFrame.from_dict(rows, orient="index")
        df.index.names = ("Date", "Item")
        return df

    ndates = 100
    nitems = 20
    dates = pd.date_range("20130101", periods=ndates, freq="D")
    items = [f"item {i}" for i in range(nitems)]

    data = {}
    for date in dates:
        nitems_for_date = nitems - random.randint(0, 12)
        levels = [
            (item, random.randint(0, 10000) / 100, random.randint(0, 10000) / 100)
            for item in items[:nitems_for_date]
        ]
        levels.sort(key=lambda x: x[1])
        data[date] = levels

    df = make_df_from_data(data)
    result = df.groupby("Date").rows[3:-3]

    sliced = {date: data[date][3:-3] for date in dates}
    expected = make_df_from_data(sliced)

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("arg", [1, 5, 30, 1000])
@pytest.mark.parametrize("method", ["head", "tail"])
@pytest.mark.parametrize("simulated", [True, False])
def test_against_head_and_tail(arg, method, simulated):
    """Test gives the same results as grouped head and tail"""

    n_groups = 100
    n_rows_per_group = 30

    data = {
        "group": [
            f"group {g}" for j in range(n_rows_per_group) for g in range(n_groups)
        ],
        "value": [
            f"group {g} row {j}"
            for j in range(n_rows_per_group)
            for g in range(n_groups)
        ],
    }
    df = pd.DataFrame(data)
    grouped = df.groupby("group")

    if method == "head":
        result = grouped.rows[:arg]

        if simulated:
            indices = []
            for j in range(arg):
                for i in range(n_groups):
                    if j * n_groups + i < n_groups * n_rows_per_group:
                        indices.append(j * n_groups + i)

            expected = df.iloc[indices]

        else:
            expected = grouped.head(arg)

    else:
        result = grouped.rows[-arg:]

        if simulated:
            indices = []
            for j in range(arg):
                for i in range(n_groups):
                    if (n_rows_per_group + j - arg) * n_groups + i >= 0:
                        indices.append((n_rows_per_group + j - arg) * n_groups + i)

            expected = df.iloc[indices]

        else:
            expected = grouped.tail(arg)

    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("start", [None, 0, 1, 10, -1, -10])
@pytest.mark.parametrize("stop", [None, 0, 1, 10, -1, -10])
@pytest.mark.parametrize("step", [None, 1, 5])
def test_against_df_iloc(start, stop, step):
    """Test that a single group gives the same results as DataFame.iloc"""

    n_rows = 30

    data = {
        "group": ["group 0"] * n_rows,
        "value": list(range(n_rows)),
    }
    df = pd.DataFrame(data)
    grouped = df.groupby("group")

    result = grouped.rows[start:stop:step]
    expected = df.iloc[start:stop:step]

    tm.assert_frame_equal(result, expected)


def test_series():
    """Test grouped Series"""

    ser = pd.Series([1, 2, 3, 4, 5], index=["a", "a", "a", "b", "b"])
    grouped = ser.groupby(level=0)
    result = grouped.rows[1:2]
    expected = pd.Series([2, 5], index=["a", "b"])

    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("step", [1, 2, 3, 4, 5])
def test_step(step):
    """Test slice with various step values"""

    data = [["x", f"x{i}"] for i in range(5)]
    data += [["y", f"y{i}"] for i in range(4)]
    data += [["z", f"z{i}"] for i in range(3)]
    df = pd.DataFrame(data, columns=["A", "B"])

    grouped = df.groupby("A")

    result = grouped.rows[::step]

    data = [["x", f"x{i}"] for i in range(0, 5, step)]
    data += [["y", f"y{i}"] for i in range(0, 4, step)]
    data += [["z", f"z{i}"] for i in range(0, 3, step)]

    index = [0 + i for i in range(0, 5, step)]
    index += [5 + i for i in range(0, 4, step)]
    index += [9 + i for i in range(0, 3, step)]

    expected = pd.DataFrame(data, columns=["A", "B"], index=index)

    tm.assert_frame_equal(result, expected)
