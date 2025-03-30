from __future__ import annotations

import numpy as np
import pytest

from dask.dataframe import from_array
from dask.dataframe.dask_expr import (
    DataFrame,
    FrameBase,
    Len,
    Series,
    concat,
    from_dict,
    from_pandas,
)
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

# Set DataFrame backend for this module
pd = _backend_library()


@pytest.fixture
def pdf():
    pdf = pd.DataFrame({"x": range(100)})
    pdf["y"] = pdf.x * 10.0
    yield pdf


@pytest.fixture
def df(pdf):
    yield from_pandas(pdf, npartitions=10)


def test_concat_str(df):
    result = str(concat([df, df], join="inner"))
    expected = "Expr=Concat(frames=[df, df], join='inner')"
    assert expected in result, result


def test_concat(pdf, df):
    result = concat([df, df])
    expected = pd.concat([pdf, pdf])
    assert_eq(result, expected)
    assert all(div is None for div in result.divisions)


def test_concat_pdf(pdf, df):
    result = concat([df, pdf])
    expected = pd.concat([pdf, pdf])
    assert_eq(result, expected)
    assert all(div is None for div in result.divisions)


def test_concat_divisions(pdf, df):
    pdf2 = pdf.set_index(np.arange(200, 300))
    df2 = from_pandas(pdf2, npartitions=10)
    result = concat([df, df2])
    expected = pd.concat([pdf, pdf2])
    assert_eq(result, expected)
    assert not any(div is None for div in result.divisions)


@pytest.mark.parametrize("join", ["right", "left"])
def test_invalid_joins(join):
    with pytest.raises(ValueError, match="'join' must be"):
        concat([df, df], join=join)


def test_concat_invalid():
    with pytest.raises(TypeError, match="dfs must"):
        concat(df)
    with pytest.raises(ValueError, match="No objects to"):
        concat([])


def test_concat_one_object(df, pdf):
    result = concat([df])
    expected = pd.concat([pdf])
    assert_eq(result, expected)
    assert not any(div is None for div in result.divisions)


def test_concat_one_no_columns(df, pdf):
    result = concat([df, df[[]]])
    expected = pd.concat([pdf, pdf[[]]])
    assert_eq(result, expected)


def test_concat_simplify(pdf, df):
    pdf2 = pdf.copy()
    pdf2["z"] = 1
    df2 = from_pandas(pdf2)
    q = concat([df, df2])[["z", "x"]]
    result = q.simplify()
    expected = concat([df[["x"]], df2[["x", "z"]]]).simplify()[["z", "x"]]
    assert result._name == expected._name

    assert_eq(q, pd.concat([pdf, pdf2])[["z", "x"]])


def test_concat_simplify_projection_not_added(pdf, df):
    pdf2 = pdf.copy()
    pdf2["z"] = 1
    df2 = from_pandas(pdf2)
    q = concat([df, df2])[["y", "x"]]
    result = q.simplify()
    expected = concat([df, df2[["x", "y"]]]).simplify()[["y", "x"]]
    assert result._name == expected._name

    assert_eq(q, pd.concat([pdf, pdf2])[["y", "x"]])


def test_concat_axis_one_co_aligned(pdf, df):
    df2 = df.add_suffix("_2")
    pdf2 = pdf.add_suffix("_2")
    assert_eq(concat([df, df2], axis=1), pd.concat([pdf, pdf2], axis=1))


def test_concat_axis_one_all_divisions_unknown(pdf):
    pdf = pdf.sort_values(by="x", ascending=False, ignore_index=True)
    pdf.index = list(reversed(pdf.index))
    df = from_pandas(pdf, npartitions=2, sort=False)
    pdf2 = pdf.add_suffix("_2")
    df2 = from_pandas(pdf2, npartitions=2, sort=False)
    with pytest.warns(UserWarning):
        assert_eq(concat([df, df2], axis=1), pd.concat([pdf, pdf2], axis=1))
    assert_eq(
        concat([df, df2], axis=1, ignore_unknown_divisions=True),
        pd.concat([pdf, pdf2], axis=1),
    )


def test_concat_axis_one_drop_dfs_not_selected(pdf, df):
    df2 = df.add_suffix("_2")
    pdf2 = pdf.add_suffix("_2")
    df3 = df.add_suffix("_3")
    pdf3 = pdf.add_suffix("_3")
    result = concat([df, df2, df3], axis=1)[["x", "y", "x_2"]].simplify()
    expected = concat([df, df2[["x_2"]]], axis=1).simplify()
    assert result._name == expected._name
    assert_eq(result, pd.concat([pdf, pdf2, pdf3], axis=1)[["x", "y", "x_2"]])


def test_concat_ignore_order():
    pdf1 = pd.DataFrame(
        {
            "x": pd.Categorical(
                ["a", "b", "c", "a"], categories=["a", "b", "c"], ordered=True
            )
        }
    )
    ddf1 = from_pandas(pdf1, 2)
    pdf2 = pd.DataFrame(
        {"x": pd.Categorical(["c", "b", "a"], categories=["c", "b", "a"], ordered=True)}
    )
    ddf2 = from_pandas(pdf2, 2)
    expected = pd.concat([pdf1, pdf2])
    expected["x"] = expected["x"].astype("category")
    result = concat([ddf1, ddf2], ignore_order=True)
    assert_eq(result, expected)


def test_concat_index(df, pdf):
    df2 = from_pandas(pdf, npartitions=3)
    result = concat([df, df2])
    expected = pd.concat([pdf, pdf])
    assert_eq(len(result), len(expected))

    query = Len(result.expr).optimize(fuse=False)
    expected = (0 + Len(df.expr) + Len(df2.expr)).optimize(fuse=False)
    assert query._name == expected._name


def test_concat_one_series(df):
    c = concat([df.x], axis=0)
    assert isinstance(c, Series)

    c = concat([df.x], axis=1)
    assert isinstance(c, DataFrame)


def test_concat_dataframe_empty():
    df = pd.DataFrame({"a": [100, 200, 300]}, dtype="int64")
    empty_df = pd.DataFrame([], dtype="int64")
    df_concat = pd.concat([df, empty_df])

    ddf = from_pandas(df, npartitions=1)
    empty_ddf = from_pandas(empty_df, npartitions=1)
    ddf_concat = concat([ddf, empty_ddf])
    assert_eq(df_concat, ddf_concat)


def test_concat_after_merge():
    pdf1 = pd.DataFrame(
        {"x": range(10), "y": [1, 2, 3, 4, 5] * 2, "z": ["cat", "dog"] * 5}
    )
    pdf2 = pd.DataFrame(
        {"i": range(10), "j": [2, 3, 4, 5, 6] * 2, "k": ["bird", "dog"] * 5}
    )
    _pdf1 = pdf1[pdf1["z"] == "cat"].merge(pdf2, left_on="y", right_on="j")
    _pdf2 = pdf1[pdf1["z"] == "dog"].merge(pdf2, left_on="y", right_on="j")
    ptotal = pd.concat([_pdf1, _pdf2])

    df1 = from_pandas(pdf1, npartitions=2)
    df2 = from_pandas(pdf2, npartitions=2)
    _df1 = df1[df1["z"] == "cat"].merge(df2, left_on="y", right_on="j")
    _df2 = df1[df1["z"] == "dog"].merge(df2, left_on="y", right_on="j")
    total = concat([_df1, _df2])

    assert_eq(total, ptotal, check_index=False)


def test_concat4_interleave_partitions():
    pdf1 = pd.DataFrame(
        np.random.randn(10, 5), columns=list("ABCDE"), index=list("abcdefghij")
    )
    pdf2 = pd.DataFrame(
        np.random.randn(13, 5), columns=list("ABCDE"), index=list("fghijklmnopqr")
    )
    pdf3 = pd.DataFrame(
        np.random.randn(13, 6), columns=list("CDEXYZ"), index=list("fghijklmnopqr")
    )

    ddf1 = from_pandas(pdf1, 2)
    ddf2 = from_pandas(pdf2, 3)
    ddf3 = from_pandas(pdf3, 2)

    cases = [
        [ddf1, ddf1],
        [ddf1, ddf2],
        [ddf1, ddf3],
        [ddf2, ddf1],
        [ddf2, ddf3],
        [ddf3, ddf1],
        [ddf3, ddf2],
    ]
    for case in cases:
        pdcase = [c.compute() for c in case]

        assert_eq(
            concat(case, interleave_partitions=True), pd.concat(pdcase, sort=False)
        )
        assert_eq(
            concat(case, join="inner", interleave_partitions=True),
            pd.concat(pdcase, join="inner", sort=False),
        )

    with pytest.raises(ValueError, match="'join' must be 'inner' or 'outer'"):
        concat([ddf1, ddf1], join="invalid", interleave_partitions=True)


def test_concat5():
    pdf1 = pd.DataFrame(
        np.random.randn(7, 5), columns=list("ABCDE"), index=list("abcdefg")
    )
    pdf2 = pd.DataFrame(
        np.random.randn(7, 6), columns=list("FGHIJK"), index=list("abcdefg")
    )
    pdf3 = pd.DataFrame(
        np.random.randn(7, 6), columns=list("FGHIJK"), index=list("cdefghi")
    )
    pdf4 = pd.DataFrame(
        np.random.randn(7, 5), columns=list("FGHAB"), index=list("cdefghi")
    )
    pdf5 = pd.DataFrame(
        np.random.randn(7, 5), columns=list("FGHAB"), index=list("fklmnop")
    )

    ddf1 = from_pandas(pdf1, 2)
    ddf2 = from_pandas(pdf2, 3)
    ddf3 = from_pandas(pdf3, 2)
    ddf4 = from_pandas(pdf4, 2)
    ddf5 = from_pandas(pdf5, 3)

    cases = [
        [ddf1, ddf2],
        [ddf1, ddf3],
        [ddf1, ddf4],
        [ddf1, ddf5],
        [ddf3, ddf4],
        [ddf3, ddf5],
        [ddf5, ddf1, ddf4],
        [ddf5, ddf3],
        [ddf1.A, ddf4.A],
        [ddf2.F, ddf3.F],
        [ddf4.A, ddf5.A],
        [ddf1.A, ddf4.F],
        [ddf2.F, ddf3.H],
        [ddf4.A, ddf5.B],
        [ddf1, ddf4.A],
        [ddf3.F, ddf2],
        [ddf5, ddf1.A, ddf2],
    ]

    for case in cases:
        pdcase = [c.compute() for c in case]

        assert_eq(
            concat(case, interleave_partitions=True),
            pd.concat(pdcase, sort=False),
        )

        assert_eq(
            concat(case, join="inner", interleave_partitions=True),
            pd.concat(pdcase, join="inner"),
        )

        assert_eq(concat(case, axis=1), pd.concat(pdcase, axis=1))

        assert_eq(
            concat(case, axis=1, join="inner"),
            pd.concat(pdcase, axis=1, join="inner"),
        )

    # Dask + pandas
    cases = [
        [ddf1, pdf2],
        [ddf1, pdf3],
        [pdf1, ddf4],
        [pdf1.A, ddf4.A],
        [ddf2.F, pdf3.F],
        [ddf1, pdf4.A],
        [ddf3.F, pdf2],
        [ddf2, pdf1, ddf3.F],
    ]

    for case in cases:
        pdcase = [c.compute() if isinstance(c, FrameBase) else c for c in case]

        assert_eq(concat(case, interleave_partitions=True), pd.concat(pdcase))

        assert_eq(
            concat(case, join="inner", interleave_partitions=True),
            pd.concat(pdcase, join="inner"),
        )

        assert_eq(concat(case, axis=1), pd.concat(pdcase, axis=1))

        assert_eq(
            concat(case, axis=1, join="inner"),
            pd.concat(pdcase, axis=1, join="inner"),
        )


def test_concat_series(pdf):
    pdf["z"] = 1
    df = from_pandas(pdf, npartitions=5)
    q = concat([df.y, df.x, df.z], axis=1)[["x", "y"]]
    df2 = df[["x", "y"]]
    expected = concat([df2.y, df2.x], axis=1)[["x", "y"]]
    assert q.optimize(fuse=False)._name == expected.optimize(fuse=False)._name
    assert_eq(q, pd.concat([pdf.y, pdf.x, pdf.z], axis=1)[["x", "y"]])


def test_concat_series_and_projection(df, pdf):
    result = concat([df.x, df.y], axis=1)["x"]
    expected = pd.concat([pdf.x, pdf.y], axis=1)["x"]
    assert_eq(result, expected)

    result = concat([df.x, df.y], axis=1)[["x"]]
    expected = pd.concat([pdf.x, pdf.y], axis=1)[["x"]]
    assert_eq(result, expected)


@pytest.mark.parametrize("npartitions", [1, 2])
@pytest.mark.parametrize("join", ["inner", "outer"])
def test_concat_single_partition_known_divisions(join, npartitions):
    df1 = from_dict({"a": [1, 2, 3], "b": [1, 2, 3]}, npartitions=npartitions)
    df2 = from_dict({"c": [1, 2]}, npartitions=npartitions)

    result = concat([df1, df2], axis=1, join=join)
    expected = pd.concat([df1.compute(), df2.compute()], axis=1, join=join)
    assert_eq(result, expected)


def test_concat_mixed_dtype_columns():
    arr = np.random.random(size=(500, 4))
    df1 = from_array(arr, chunksize=200, columns=[0, 1, 2, "_y"])
    df2 = from_array(arr, chunksize=200, columns=[0, 1, 2, "_y"])
    result = concat([df1, df2]).drop(["_y"], axis=1)
    expected = pd.concat([df1.compute(), df2.compute()]).drop(columns="_y")
    assert_eq(result, expected)


@pytest.mark.parametrize("axis", [0, 1])
@pytest.mark.parametrize("interleave_partitions", [True, False])
def test_concat_optimize_project(axis, interleave_partitions):
    df1 = from_dict({"a": range(10), "b": range(10)}, npartitions=2)
    df2 = from_dict({"c": [5, 2, 3] + list(range(7))}, npartitions=3)
    df2.clear_divisions()
    concatenated = concat(
        [df1, df2], axis=axis, interleave_partitions=interleave_partitions
    )
    optimized = concatenated.optimize()
    optimized_project = optimized[["a", "c"]]
    assert_eq(optimized_project, concatenated[["a", "c"]])
    assert (
        optimized_project.optimize()._name == concatenated[["a", "c"]].optimize()._name
    )
