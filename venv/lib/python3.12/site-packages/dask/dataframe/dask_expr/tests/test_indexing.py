from __future__ import annotations

import numpy as np
import pytest

from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr.io import FromPandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

pd = _backend_library()


@pytest.fixture
def pdf():
    pdf = pd.DataFrame({"x": range(20)})
    pdf["y"] = pdf.x * 10.0
    yield pdf


@pytest.fixture
def df(pdf):
    yield from_pandas(pdf, npartitions=4)


def test_iloc(df, pdf):
    assert_eq(df.iloc[:, 1], pdf.iloc[:, 1])
    assert_eq(df.iloc[:, [1]], pdf.iloc[:, [1]])
    assert_eq(df.iloc[:, [0, 1]], pdf.iloc[:, [0, 1]])
    assert_eq(df.iloc[:, []], pdf.iloc[:, []])


def test_iloc_errors(df):
    with pytest.raises(NotImplementedError):
        df.iloc[1]
    with pytest.raises(NotImplementedError):
        df.iloc[1, 1]
    with pytest.raises(ValueError, match="Too many"):
        df.iloc[(1, 2, 3)]


def test_loc_slice(pdf, df):
    pdf.columns = [10, 20]
    df.columns = [10, 20]
    assert_eq(df.loc[:, :15], pdf.loc[:, :15])
    assert_eq(df.loc[:, 15:], pdf.loc[:, 15:])
    assert_eq(df.loc[:, 25:], pdf.loc[:, 25:])  # no columns
    assert_eq(df.loc[:, ::-1], pdf.loc[:, ::-1])


def test_iloc_slice(df, pdf):
    assert_eq(df.iloc[:, :1], pdf.iloc[:, :1])
    assert_eq(df.iloc[:, 1:], pdf.iloc[:, 1:])
    assert_eq(df.iloc[:, 99:], pdf.iloc[:, 99:])  # no columns
    assert_eq(df.iloc[:, ::-1], pdf.iloc[:, ::-1])


@pytest.mark.parametrize("loc", [False, True])
@pytest.mark.parametrize("update", [False, True])
def test_columns_dtype_on_empty_slice(df, pdf, loc, update):
    pdf.columns = [10, 20]
    if update:
        df.columns = [10, 20]
    else:
        df = from_pandas(pdf, npartitions=10)

    assert df.columns.dtype == pdf.columns.dtype
    assert df.compute().columns.dtype == pdf.columns.dtype
    assert_eq(df, pdf)

    if loc:
        df = df.loc[:, []]
        pdf = pdf.loc[:, []]
    else:
        df = df[[]]
        pdf = pdf[[]]

    assert df.columns.dtype == pdf.columns.dtype
    assert df.compute().columns.dtype == pdf.columns.dtype
    assert_eq(df, pdf)


def test_loc(df, pdf):
    assert_eq(df.loc[:, "x"], pdf.loc[:, "x"])
    assert_eq(df.loc[:, ["x"]], pdf.loc[:, ["x"]])
    assert_eq(df.loc[:, []], pdf.loc[:, []])

    assert_eq(df.loc[df.y == 20, "x"], pdf.loc[pdf.y == 20, "x"])
    assert_eq(df.loc[df.y == 20, ["x"]], pdf.loc[pdf.y == 20, ["x"]])
    assert df.loc[3:8].divisions[0] == 3
    assert df.loc[3:8].divisions[-1] == 8

    assert df.loc[5].divisions == (5, 5)

    assert_eq(df.loc[5], pdf.loc[5:5])
    assert_eq(df.loc[3:8], pdf.loc[3:8])
    assert_eq(df.loc[:8], pdf.loc[:8])
    assert_eq(df.loc[3:], pdf.loc[3:])
    assert_eq(df.loc[[5]], pdf.loc[[5]])

    assert_eq(df.x.loc[5], pdf.x.loc[5:5])
    assert_eq(df.x.loc[3:8], pdf.x.loc[3:8])
    assert_eq(df.x.loc[:8], pdf.x.loc[:8])
    assert_eq(df.x.loc[3:], pdf.x.loc[3:])
    assert_eq(df.x.loc[[5]], pdf.x.loc[[5]])
    assert_eq(df.x.loc[[]], pdf.x.loc[[]])
    assert_eq(df.x.loc[np.array([])], pdf.x.loc[np.array([])])

    pytest.raises(KeyError, lambda: df.loc[1000])
    assert_eq(df.loc[1000:], pdf.loc[1000:])
    assert_eq(df.loc[1000:2000], pdf.loc[1000:2000])
    assert_eq(df.loc[:-1000], pdf.loc[:-1000])
    assert_eq(df.loc[-2000:-1000], pdf.loc[-2000:-1000])


def test_loc_non_informative_index():
    df = pd.DataFrame({"x": [1, 2, 3, 4]}, index=[10, 20, 30, 40])
    ddf = from_pandas(df, npartitions=2, sort=True).clear_divisions()
    assert not ddf.known_divisions

    ddf.loc[20:30].compute(scheduler="sync")

    assert_eq(ddf.loc[20:30], df.loc[20:30])

    df = pd.DataFrame({"x": [1, 2, 3, 4]}, index=[10, 20, 20, 40])
    ddf = from_pandas(df, npartitions=2, sort=True)
    assert_eq(ddf.loc[20], df.loc[20:20])


def test_loc_with_series(df, pdf):
    assert_eq(df.loc[df.x % 2 == 0], pdf.loc[pdf.x % 2 == 0])


def test_loc_with_array(df, pdf):
    assert_eq(df.loc[(df.x % 2 == 0).values], pdf.loc[(pdf.x % 2 == 0).values])


def test_loc_with_function(df, pdf):
    assert_eq(df.loc[lambda df: df["x"] > 3, :], pdf.loc[lambda df: df["x"] > 3, :])

    def _col_loc_fun(_df):
        return _df.columns.str.contains("y")

    assert_eq(df.loc[:, _col_loc_fun], pdf.loc[:, _col_loc_fun])


def test_getitem_align():
    df = pd.DataFrame(
        {
            "A": [1, 2, 3, 4, 5, 6, 7, 8, 9],
            "B": [9, 8, 7, 6, 5, 4, 3, 2, 1],
            "C": [True, False, True] * 3,
        },
        columns=list("ABC"),
    )
    ddf = from_pandas(df, 2)
    assert_eq(ddf[ddf.C.repartition([0, 2, 5, 8])], df[df.C])


def test_loc_bool_cindex():
    # https://github.com/dask/dask/issues/11015
    pdf = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    ddf = from_pandas(pdf, npartitions=1)
    indexer = [True, False]
    assert_eq(pdf.loc[:, indexer], ddf.loc[:, indexer])


def test_loc_slicing():
    npartitions = 10
    pdf = pd.DataFrame(
        {
            "A": np.random.randn(npartitions * 10),
        },
        index=pd.date_range("2024-01-01", "2024-12-31", npartitions * 10),
    )
    df = from_pandas(pdf, npartitions=npartitions)
    result = df["2024-03-01":"2024-09-30"]["A"]
    assert_eq(result, pdf["2024-03-01":"2024-09-30"]["A"])


def test_reverse_indexing(df, pdf):
    assert_eq(df.loc[::-1], pdf.loc[::-1])
    result = df.loc[::-1][["x"]]
    expected = pdf.loc[::-1][["x"]]
    assert_eq(result, expected)
    expr = result.expr.optimize(fuse=False)
    assert expr.frame.frame.columns == ["x"] and isinstance(
        expr.frame.frame, FromPandas
    )
    with pytest.raises(
        ValueError, match="Can not use loc on DataFrame without known divisions"
    ):
        df.loc[1::-1].loc[5]
    with pytest.raises(
        ValueError, match="Can not use loc on DataFrame without known divisions"
    ):
        df.loc[:5:-1].loc[5]


def test_indexing_element_index():
    pdf = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    result = from_pandas(pdf, 2).loc[2].index
    pd.testing.assert_index_equal(result.compute(), pdf.loc[[2]].index)

    result = from_pandas(pdf, 2).loc[[2]].index
    pd.testing.assert_index_equal(result.compute(), pdf.loc[[2]].index)
