import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
from pandas import DataFrame
import pandas._testing as tm

msg = (
    "Passing labels, `like`, or `regex` to DataFrame.filter is deprecated and "
    "will be removed in a future version; filter will only accept a boolean mask"
)


class TestDataFrameFilter:
    def test_filter(self, float_frame, float_string_frame):
        # Items
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            filtered = float_frame.filter(["A", "B", "E"])
        assert len(filtered.columns) == 2
        assert "E" not in filtered

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            filtered = float_frame.filter(["A", "B", "E"], axis="columns")
        assert len(filtered.columns) == 2
        assert "E" not in filtered

        # Other axis
        idx = float_frame.index[0:4]
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            filtered = float_frame.filter(idx, axis="index")
        expected = float_frame.reindex(index=idx)
        tm.assert_frame_equal(filtered, expected)

        # like
        fcopy = float_frame.copy()
        fcopy["AA"] = 1

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            filtered = fcopy.filter(like="A")
        assert len(filtered.columns) == 2
        assert "AA" in filtered

        # like with ints in column names
        df = DataFrame(0.0, index=[0, 1, 2], columns=[0, 1, "_A", "_B"])
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            filtered = df.filter(like="_")
        assert len(filtered.columns) == 2

        # regex with ints in column names
        # from PR #10384
        df = DataFrame(0.0, index=[0, 1, 2], columns=["A1", 1, "B", 2, "C"])
        expected = DataFrame(
            0.0, index=[0, 1, 2], columns=pd.Index([1, 2], dtype=object)
        )
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            filtered = df.filter(regex="^[0-9]+$")
        tm.assert_frame_equal(filtered, expected)

        expected = DataFrame(0.0, index=[0, 1, 2], columns=[0, "0", 1, "1"])
        # shouldn't remove anything
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            filtered = expected.filter(regex="^[0-9]+$")
        tm.assert_frame_equal(filtered, expected)

        # pass in None
        no_arg_msg = (
            "DataFrame.filter requires a boolean mask \\(or the deprecated "
            "`items`, `like`, or `regex`\\)"
        )
        with pytest.raises(TypeError, match=no_arg_msg):
            float_frame.filter()
        with pytest.raises(TypeError, match=no_arg_msg):
            float_frame.filter(items=None)
        with pytest.raises(TypeError, match=no_arg_msg):
            float_frame.filter(axis=1)

        # test mutually exclusive arguments
        with pytest.raises(TypeError, match="mutually exclusive"):
            float_frame.filter(items=["one", "three"], regex="e$", like="bbi")
        with pytest.raises(TypeError, match="mutually exclusive"):
            float_frame.filter(items=["one", "three"], regex="e$", axis=1)
        with pytest.raises(TypeError, match="mutually exclusive"):
            float_frame.filter(items=["one", "three"], regex="e$")
        with pytest.raises(TypeError, match="mutually exclusive"):
            float_frame.filter(items=["one", "three"], like="bbi", axis=0)
        with pytest.raises(TypeError, match="mutually exclusive"):
            float_frame.filter(items=["one", "three"], like="bbi")

        # objects
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            filtered = float_string_frame.filter(like="foo")
        assert "foo" in filtered

        # unicode columns, won't ascii-encode
        df = float_frame.rename(columns={"B": "\u2202"})
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            filtered = df.filter(like="C")
        assert "C" in filtered

    def test_filter_regex_search(self, float_frame):
        fcopy = float_frame.copy()
        fcopy["AA"] = 1

        # regex
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            filtered = fcopy.filter(regex="[A]+")
        assert len(filtered.columns) == 2
        assert "AA" in filtered

        # doesn't have to be at beginning
        df = DataFrame(
            {"aBBa": [1, 2], "BBaBB": [1, 2], "aCCa": [1, 2], "aCCaBB": [1, 2]}
        )

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = df.filter(regex="BB")
        exp = df[[x for x in df.columns if "BB" in x]]
        tm.assert_frame_equal(result, exp)

    @pytest.mark.parametrize(
        "name,expected_data",
        [
            ("a", {"a": [1, 2]}),
            ("あ", {"あ": [3, 4]}),
        ],
    )
    def test_filter_unicode(self, name, expected_data):
        # GH13101
        df = DataFrame({"a": [1, 2], "あ": [3, 4]})
        expected = DataFrame(expected_data)

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            tm.assert_frame_equal(df.filter(like=name), expected)
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            tm.assert_frame_equal(df.filter(regex=name), expected)

    def test_filter_bytestring(self):
        # GH13101
        name = "a"
        df = DataFrame({b"a": [1, 2], b"b": [3, 4]})
        expected = DataFrame({b"a": [1, 2]})

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            tm.assert_frame_equal(df.filter(like=name), expected)
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            tm.assert_frame_equal(df.filter(regex=name), expected)

    def test_filter_corner(self):
        empty = DataFrame()

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = empty.filter([])
        tm.assert_frame_equal(result, empty)

        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = empty.filter(like="foo")
        tm.assert_frame_equal(result, empty)

    def test_filter_regex_non_string(self):
        # GH#5798 trying to filter on non-string columns should drop,
        #  not raise
        df = DataFrame(np.random.default_rng(2).random((3, 2)), columns=["STRING", 123])
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = df.filter(regex="STRING")
        expected = df[["STRING"]]
        tm.assert_frame_equal(result, expected)

    def test_filter_keep_order(self):
        # GH#54980
        df = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = df.filter(items=["B", "A"])
        expected = df[["B", "A"]]
        tm.assert_frame_equal(result, expected)

    def test_filter_different_dtype(self):
        # GH#54980
        df = DataFrame({1: [1, 2, 3], 2: [4, 5, 6]})
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = df.filter(items=["B", "A"])
        expected = df[[]]
        tm.assert_frame_equal(result, expected)


@pytest.fixture
def df():
    return DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}, index=["x", "y", "z"])


@pytest.mark.parametrize(
    "mask",
    [
        [False, True, True],
        np.array([False, True, True]),
        pd.array([False, True, True], dtype="boolean"),
        pd.Index([False, True, True]),
        lambda df: df["a"] > 1,
        pd.col("a") > 1,
    ],
)
def test_filter_mask_rows(df, mask):
    # GH#61317
    result = df.filter(mask)
    expected = df.iloc[1:]
    tm.assert_frame_equal(result, expected)


def test_filter_mask_series_aligns(df):
    # GH#61317
    mask = pd.Series([True, True, False], index=["z", "y", "x"])
    result = df.filter(mask)
    expected = df.iloc[1:]
    tm.assert_frame_equal(result, expected)


def test_filter_mask_series_unalignable(df):
    # GH#61317
    mask = pd.Series([True, True, False], index=["z", "y", "w"])
    msg = (
        "Unalignable boolean Series provided as indexer \\(index of the boolean "
        "Series and of the indexed object do not match\\)."
    )
    with pytest.raises(pd.errors.IndexingError, match=msg):
        df.filter(mask)


def test_filter_mask_wrong_length(df):
    # GH#61317
    msg = "Boolean index has wrong length: 2 instead of 3"
    with pytest.raises(IndexError, match=msg):
        df.filter(np.array([True, False]))


@pytest.mark.parametrize("axis", [1, "columns"])
def test_filter_mask_columns(df, axis):
    # GH#61317
    result = df.filter([False, True], axis=axis)
    expected = df[["b"]]
    tm.assert_frame_equal(result, expected)


def test_filter_mask_2d_raises(df):
    # GH#61317
    msg = "The mask passed to DataFrame.filter must be one-dimensional"
    with pytest.raises(ValueError, match=msg):
        df.filter(df > 1)


def test_filter_callable_must_return_mask(df):
    # GH#61317
    msg = "The callable passed to DataFrame.filter must evaluate to a boolean mask"
    with pytest.raises(TypeError, match=msg):
        df.filter(lambda df: ["a"])


def test_filter_expression_must_return_mask(df):
    # GH#61317
    msg = "The expression passed to DataFrame.filter must evaluate to a boolean mask"
    with pytest.raises(TypeError, match=msg):
        df.filter(pd.col("a"))


@pytest.mark.parametrize(
    "mask",
    [
        [True, None, False],
        np.array([True, None, False], dtype=object),
        pd.array([True, None, False], dtype="boolean"),
        pd.Series([True, None, False], dtype="boolean", index=["x", "y", "z"]),
    ],
)
def test_filter_mask_na(df, mask):
    # GH#61317
    # missing values are treated as False by default, matching df[mask]
    result = df.filter(mask)
    expected = df.iloc[[0]]
    tm.assert_frame_equal(result, expected)

    result = df.filter(mask, na=False)
    tm.assert_frame_equal(result, expected)

    msg = "The mask contains missing values"
    with pytest.raises(ValueError, match=msg):
        df.filter(mask, na="raise")

    result = df.filter(mask, na=True)
    expected = df.iloc[[0, 1]]
    tm.assert_frame_equal(result, expected)


def test_filter_mask_pyarrow_na(df):
    # GH#61317
    pytest.importorskip("pyarrow")
    mask = pd.array([True, None, False], dtype="bool[pyarrow]")
    msg = "The mask contains missing values"
    with pytest.raises(ValueError, match=msg):
        df.filter(mask, na="raise")
    result = df.filter(mask)
    expected = df.iloc[[0]]
    tm.assert_frame_equal(result, expected)
    result = df.filter(mask, na=True)
    expected = df.iloc[[0, 1]]
    tm.assert_frame_equal(result, expected)


def test_filter_invalid_na(df):
    # GH#61317
    msg = "na must be 'raise', True, or False, got 'ignore'"
    with pytest.raises(ValueError, match=msg):
        df.filter([True, False, True], na="ignore")


@pytest.mark.parametrize(
    "mask",
    [
        [True, False],
        np.array([True, False]),
        pd.Series([True, False]),
        pd.Index([True, False]),
    ],
)
def test_filter_bool_labels_deprecated(mask):
    # GH#61317
    # Boolean values on an axis with boolean labels keep selecting
    # labels until the deprecation is enforced.
    df = DataFrame({True: [1], False: [2], "c": [3]})
    msg = (
        "DataFrame.filter with boolean values currently selects the labels True "
        "and False when the axis contains boolean labels. In a future version the "
        "values will be treated as a boolean mask. Use obj.loc with "
        "obj.columns.isin\\(items\\) \\(or obj.index.isin\\(items\\)\\) to select "
        "these labels instead."
    )
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = df.filter(mask)
    expected = df.iloc[:, :2]
    tm.assert_frame_equal(result, expected)

    # the legacy default axis (columns) has no boolean labels, so this is
    # a mask on the index even though the index has boolean labels
    df = DataFrame({"a": [1, 2]}, index=[True, False])
    if isinstance(mask, pd.Series):
        mask = mask.set_axis(df.index)
    result = df.filter(mask)
    expected = df.iloc[[0]]
    tm.assert_frame_equal(result, expected)


def test_filter_bool_list_int_labels_is_mask():
    # GH#61317
    # int labels are not boolean labels, even though True == 1
    df = DataFrame({0: [1], 1: [2]})
    result = df.filter([True, False], axis=1)
    expected = df[[0]]
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "kwargs, hint",
    [
        (
            {"items": ["a"]},
            "obj.loc\\(axis='columns'\\)\\[obj.columns.intersection\\(items\\)\\]",
        ),
        (
            {"items": ["x"], "axis": 0},
            "obj.loc\\(axis='index'\\)\\[obj.index.intersection\\(items\\)\\]",
        ),
        (
            {"like": "a"},
            "obj.filter\\(obj.columns.str.contains\\(like, regex=False\\), "
            "axis='columns'\\)",
        ),
        (
            {"like": "x", "axis": "index"},
            "obj.filter\\(obj.index.str.contains\\(like, regex=False\\),...*",
        ),
        (
            {"regex": "a"},
            "obj.filter\\(obj.columns.str.contains\\(regex\\), axis='columns'\\)",
        ),
        (
            {"regex": "x", "axis": 0},
            "obj.filter\\(obj.index.str.contains\\(regex\\), axis='index'\\)",
        ),
    ],
)
def test_filter_labels_deprecated(df, kwargs, hint):
    # GH#61317
    msg = (
        "Passing labels, `like`, or `regex` to DataFrame.filter is deprecated and "
        "will be removed in a future version; filter will only accept a boolean "
        f"mask. Use {hint} instead."
    )
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = df.filter(**kwargs)
    if kwargs.get("axis") in (0, "index"):
        expected = df.iloc[[0]]
    else:
        expected = df[["a"]]
    tm.assert_frame_equal(result, expected)


def test_filter_mask_default_axis_is_index(df):
    # GH#61317
    # unlike the deprecated label-based usage, a mask defaults to the index
    result = df.filter([True, False, True])
    expected = df.iloc[[0, 2]]
    tm.assert_frame_equal(result, expected)
