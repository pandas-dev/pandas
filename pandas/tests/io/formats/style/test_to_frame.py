import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm

pytest.importorskip("jinja2")


@pytest.fixture
def mi_df():
    index = pd.MultiIndex.from_product([["a", "b"], [1, 2]], names=["l0", "l1"])
    columns = pd.MultiIndex.from_product([["X", "Y"], ["p", "q"]], names=["c0", "c1"])
    return pd.DataFrame(np.arange(16.0).reshape(4, 4), index=index, columns=columns)


def test_to_frame_formatted_values():
    # GH#36680
    df = pd.DataFrame([[1.9999, 1.9999, 1.9999, 1.9999]], columns=list("ABCD"))
    styler = df.style.format(
        {
            "A": lambda x: float(f"{x:.0f}"),
            "B": lambda x: int(f"{x:.0f}"),
            "C": "{:.1f}",
        },
        precision=3,
    )
    result = styler.to_frame()
    expected = pd.DataFrame(
        [[2.0, 2, "2.0", "2.000"]],
        index=pd.Index(["0"], dtype=object),
        columns=pd.Index(list("ABCD"), dtype=object),
        dtype=object,
    )
    tm.assert_frame_equal(result, expected)
    assert (
        result.to_json()
        == '{"A":{"0":2.0},"B":{"0":2},"C":{"0":"2.0"},"D":{"0":"2.000"}}'
    )


def test_to_frame_na_rep():
    df = pd.DataFrame({"A": [1.5, np.nan]})
    result = df.style.format(na_rep="-", precision=1).to_frame()
    expected = pd.DataFrame(
        {"A": ["1.5", "-"]},
        index=pd.Index(["0", "1"], dtype=object),
        columns=pd.Index(["A"], dtype=object),
        dtype=object,
    )
    tm.assert_frame_equal(result, expected)


def test_to_frame_formatters_get_python_scalars():
    # as when rendering, formatters receive Python rather than NumPy scalars
    df = pd.DataFrame({"a": [1000, 2]}, index=[10, 20])
    styler = df.style.format(lambda x: f"{x:,}" if isinstance(x, int) else "?")
    styler = styler.format_index(lambda x: type(x).__name__)
    result = styler.to_frame()
    expected = pd.DataFrame(
        {"a": ["1,000", "2"]},
        index=pd.Index(["int", "int"], dtype=object),
        dtype=object,
    )
    expected.columns = pd.Index(["a"], dtype=object)
    tm.assert_frame_equal(result, expected)


def test_to_frame_non_numeric_unchanged():
    ts = pd.Timestamp("2020-01-01")
    df = pd.DataFrame({"a": [True, None], "b": [ts, pd.NaT]})
    result = df.style.to_frame()
    assert result["a"].tolist() == [True, None]
    assert result["b"].tolist()[0] == ts
    assert result["b"].tolist()[1] is pd.NaT


def test_to_frame_hidden_and_formatted_labels(mi_df):
    styler = (
        mi_df.style.format("{:.1f}")
        .format_index(str.upper, axis=0, level=0)
        .format_index_names(lambda x: x + "!", axis=1)
        .hide([("a", 2)], axis=0)
        .hide(subset=[("X", "p")], axis=1)
    )
    result = styler.to_frame()
    expected = pd.DataFrame(
        [["1.0", "2.0", "3.0"], ["9.0", "10.0", "11.0"], ["13.0", "14.0", "15.0"]],
        index=pd.MultiIndex.from_arrays(
            [
                pd.Index(["A", "B", "B"], dtype=object),
                pd.Index(["1", "1", "2"], dtype=object),
            ],
            names=["l0", "l1"],
        ),
        columns=pd.MultiIndex.from_arrays(
            [
                pd.Index(["X", "Y", "Y"], dtype=object),
                pd.Index(["q", "p", "q"], dtype=object),
            ],
            names=["c0!", "c1!"],
        ),
        dtype=object,
    )
    tm.assert_frame_equal(result, expected)


def test_to_frame_hidden_levels(mi_df):
    result = mi_df.style.format("{:.0f}").hide(level=0, axis=0).hide(axis=1, level="c1")
    result = result.to_frame()
    expected = pd.DataFrame(
        [[str(4 * r + c) for c in range(4)] for r in range(4)],
        index=pd.Index(["1", "2", "1", "2"], name="l1", dtype=object),
        columns=pd.Index(["X", "X", "Y", "Y"], name="c0", dtype=object),
        dtype=object,
    )
    tm.assert_frame_equal(result, expected)


def test_to_frame_hide_index(mi_df):
    result = mi_df.style.hide(axis=0).to_frame()
    assert isinstance(result.index, pd.RangeIndex)
    tm.assert_index_equal(result.index, pd.RangeIndex(4))


@pytest.mark.parametrize("axis", [0, 1])
def test_to_frame_hide_names(mi_df, axis):
    result = mi_df.style.hide(names=True, axis=axis).to_frame()
    assert list(result.axes[axis].names) == [None, None]
    assert list(result.axes[1 - axis].names) != [None, None]


def test_to_frame_concat(mi_df):
    other = mi_df.iloc[:1].style.format("{:.2f}")
    styler = mi_df.style.format("{:.0f}").hide(subset=[("X", "p")], axis=1)
    result = styler.concat(other).to_frame()
    assert result.shape == (5, 3)
    assert list(result.iloc[-1]) == ["1.00", "2.00", "3.00"]
    # the concatenated Styler itself is left unchanged
    assert other.hidden_columns == []


def test_to_frame_concat_uses_outer_headers():
    df = pd.DataFrame({"a": [1, 2], "b": [3, 4]})
    styler = df.style.format_index(str.upper, axis=1).hide(axis=0)
    result = styler.concat(df.style.format("{:.1f}")).to_frame()
    expected = pd.DataFrame(
        [["1", "3"], ["2", "4"], ["1.0", "3.0"], ["2.0", "4.0"]],
        columns=pd.Index(["A", "B"], dtype=object),
        dtype=object,
    )
    tm.assert_frame_equal(result, expected)


def test_to_frame_concat_keeps_index_names():
    df = pd.DataFrame({"a": [1.0, 2.0]}, index=pd.Index(["x", "y"], name="row"))
    result = df.style.concat(df.agg(["sum"]).style).to_frame()
    assert result.index.names == ["row"]


def test_to_frame_tuple_labels():
    df = pd.DataFrame({"a": [1]}, index=pd.Index([("x", "y")], tupleize_cols=False))
    result = df.style.format_index(lambda x: (x, x), axis=1).to_frame()
    assert result.index.nlevels == 1
    assert result.index[0] == ("x", "y")
    assert result.columns.nlevels == 1
    assert result.columns[0] == ("a", "a")


def test_to_frame_concat_hidden_level(mi_df):
    styler = mi_df.style.format("{:.0f}").hide(level=0, axis=0)
    result = styler.concat(mi_df.iloc[:1].style).to_frame()
    assert list(result.index) == ["1", "2", "1", "2", "1"]
    assert result.index.names == ["l1"]


@pytest.mark.parametrize("axis", [0, 1])
def test_to_frame_single_level_multiindex(axis):
    df = pd.DataFrame({"v": [1]}, index=pd.MultiIndex.from_arrays([["x"]], names=["k"]))
    df = df if axis == 0 else df.T
    result = df.style.format_index(str.upper, axis=axis).to_frame()
    tm.assert_index_equal(result.axes[axis], pd.Index(["X"], name="k", dtype=object))


def test_to_frame_partially_named_levels():
    index = pd.MultiIndex.from_product([["a"], [1]], names=[None, "l1"])
    df = pd.DataFrame({"v": [1]}, index=index)
    result = df.style.format_index_names(lambda x: x + "!", axis=0).to_frame()
    assert result.index.names == [None, "l1!"]


def test_to_frame_empty():
    result = pd.DataFrame().style.to_frame()
    assert result.shape == (0, 0)


def test_to_frame_duplicate_columns():
    df = pd.DataFrame([[1, 2]], columns=["a", "a"])
    result = df.style.to_frame()
    expected = pd.DataFrame(
        [["1", "2"]],
        index=pd.Index(["0"], dtype=object),
        columns=pd.Index(["a", "a"], dtype=object),
        dtype=object,
    )
    tm.assert_frame_equal(result, expected)
