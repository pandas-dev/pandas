import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_agg_relabel():
    # GH 26513
    df = pd.DataFrame({"A": [1, 2, 1, 2], "B": [1, 2, 3, 4], "C": [3, 4, 5, 6]})

    # simplest case with one column, one func
    result = df.agg(foo=("B", "sum"))
    expected = pd.DataFrame({"B": [10]}, index=pd.Index(["foo"]))
    tm.assert_frame_equal(result, expected)

    # test on same column with different methods
    result = df.agg(foo=("B", "sum"), bar=("B", "min"))
    expected = pd.DataFrame({"B": [10, 1]}, index=pd.Index(["foo", "bar"]))

    tm.assert_frame_equal(result, expected)


def test_agg_relabel_multi_columns_multi_methods():
    # GH 26513, test on multiple columns with multiple methods
    df = pd.DataFrame({"A": [1, 2, 1, 2], "B": [1, 2, 3, 4], "C": [3, 4, 5, 6]})
    result = df.agg(
        foo=("A", "sum"),
        bar=("B", "mean"),
        cat=("A", "min"),
        dat=("B", "max"),
        f=("A", "max"),
        g=("C", "min"),
    )
    expected = pd.DataFrame(
        {
            "A": [6.0, np.nan, 1.0, np.nan, 2.0, np.nan],
            "B": [np.nan, 2.5, np.nan, 4.0, np.nan, np.nan],
            "C": [np.nan, np.nan, np.nan, np.nan, np.nan, 3.0],
        },
        index=pd.Index(["foo", "bar", "cat", "dat", "f", "g"]),
    )
    tm.assert_frame_equal(result, expected)


def test_agg_relabel_partial_functions():
    # GH 26513, test on partial, functools or more complex cases
    df = pd.DataFrame({"A": [1, 2, 1, 2], "B": [1, 2, 3, 4], "C": [3, 4, 5, 6]})
    result = df.agg(foo=("A", np.mean), bar=("A", "mean"), cat=("A", min))
    expected = pd.DataFrame(
        {"A": [1.5, 1.5, 1.0]}, index=pd.Index(["foo", "bar", "cat"])
    )
    tm.assert_frame_equal(result, expected)

    result = df.agg(
        foo=("A", min),
        bar=("B", np.min),
        cat=("B", max),
        dat=("C", "min"),
        f=("B", np.sum),
        kk=("B", lambda x: min(x)),
    )
    expected = pd.DataFrame(
        {
            "A": [1.0, np.nan, np.nan, np.nan, np.nan, np.nan],
            "B": [np.nan, 1.0, 4.0, np.nan, 10.0, 1.0],
            "C": [np.nan, np.nan, np.nan, 3.0, np.nan, np.nan],
        },
        index=pd.Index(["foo", "bar", "cat", "dat", "f", "kk"]),
    )
    tm.assert_frame_equal(result, expected)


def test_agg_namedtuple():
    # GH 26513
    df = pd.DataFrame({"A": [0, 1], "B": [1, 2]})
    result = df.agg(
        foo=pd.NamedAgg("B", "sum"),
        bar=pd.NamedAgg("B", "min"),
        cat=pd.NamedAgg(column="B", aggfunc="count"),
        fft=pd.NamedAgg("B", aggfunc="max"),
    )

    expected = pd.DataFrame(
        {"B": [3, 1, 2, 2]}, index=pd.Index(["foo", "bar", "cat", "fft"])
    )
    tm.assert_frame_equal(result, expected)

    result = df.agg(
        foo=pd.NamedAgg("A", "min"),
        bar=pd.NamedAgg(column="B", aggfunc="max"),
        cat=pd.NamedAgg(column="A", aggfunc="max"),
    )
    expected = pd.DataFrame(
        {"A": [0.0, np.nan, 1.0], "B": [np.nan, 2.0, np.nan]},
        index=pd.Index(["foo", "bar", "cat"]),
    )
    tm.assert_frame_equal(result, expected)


def test_reconstruct_func():
    # GH 28472, test to ensure reconstruct_func isn't moved;
    # This method is used by other libraries (e.g. dask)
    result = pd.core.apply.reconstruct_func("min")
    expected = (False, "min", None, None)
    tm.assert_equal(result, expected)


def test_reconstruct_func_allow_skip_normalization():
    # GH#63743 the normalization fast path is opt-in: only callers that pass
    #  allow_skip_normalization can receive the un-normalized {column: aggfunc}
    #  form, and only when every output name equals its source column name
    kwargs = {"B": ("B", "sum")}

    result = pd.core.apply.reconstruct_func(None, True, **kwargs)
    assert result == (False, {"B": "sum"}, None, None)

    # default (every caller other than DataFrameGroupBy.aggregate) normalizes
    relabeling, func, columns, order = pd.core.apply.reconstruct_func(None, **kwargs)
    assert relabeling
    assert func == {"B": ["sum"]}
    assert columns == ("B",)
    tm.assert_numpy_array_equal(order, np.array([0], dtype=np.intp))

    # a renamed output still normalizes even when the fast path is allowed
    relabeling, _, columns, _ = pd.core.apply.reconstruct_func(
        None, True, x=("B", "sum")
    )
    assert relabeling
    assert columns == ("x",)


def test_agg_relabel_output_name_matches_column():
    # GH#63743 named aggregation returns a DataFrame indexed by the output
    #  names even when those names match the source column names
    df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})

    result = df.agg(A=("A", "min"))
    expected = pd.DataFrame({"A": [1]}, index=pd.Index(["A"]))
    tm.assert_frame_equal(result, expected)

    result = df.agg(A=("A", "min"), B=("B", "max"))
    expected = pd.DataFrame(
        {"A": [1.0, np.nan], "B": [np.nan, 6.0]}, index=pd.Index(["A", "B"])
    )
    tm.assert_frame_equal(result, expected)


def test_agg_relabel_output_name_matches_column_axis_1():
    # GH#63743 named aggregation with axis=1 raises regardless of whether the
    #  output name matches the column name
    df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, index=["A", "B", "C"])
    msg = "Named aggregation is not supported when axis=1."

    with pytest.raises(NotImplementedError, match=msg):
        df.agg(A=("A", "min"), axis=1)

    with pytest.raises(NotImplementedError, match=msg):
        df.agg(x=("A", "min"), axis=1)
