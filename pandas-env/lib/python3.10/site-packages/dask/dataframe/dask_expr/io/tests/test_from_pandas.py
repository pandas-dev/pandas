from __future__ import annotations

import copy

import pytest

import dask
from dask.dataframe.dask_expr import from_pandas, repartition
from dask.dataframe.dask_expr.tests._util import _backend_library
from dask.dataframe.utils import assert_eq, pyarrow_strings_enabled

pd = _backend_library()


@pytest.fixture(params=["Series", "DataFrame"])
def pdf(request):
    out = pd.DataFrame({"x": [1, 4, 3, 2, 0, 5]})
    return out["x"] if request.param == "Series" else out


@pytest.mark.parametrize("sort", [True, False])
def test_from_pandas(pdf, sort):
    df = from_pandas(pdf, npartitions=2, sort=sort)

    assert df.npartitions == 2
    assert df.divisions == (0, 3, 5)
    assert_eq(df, pdf, sort_results=sort)
    assert all(tsk.data_producer for tsk in df.dask.values())


def test_from_pandas_noargs(pdf):
    df = from_pandas(pdf)

    assert df.npartitions == 1
    assert df.divisions == (0, 5)
    assert_eq(df, pdf)


def test_from_pandas_empty(pdf):
    pdf = pdf.iloc[:0]
    df = from_pandas(pdf, npartitions=2)
    assert_eq(pdf, df)


def test_from_pandas_immutable(pdf):
    expected = pdf.copy()
    df = from_pandas(pdf)
    pdf.iloc[0] = 100
    assert_eq(df, expected)


def test_from_pandas_sort_and_different_partitions():
    pdf = pd.DataFrame({"a": [1, 2, 3] * 3, "b": 1}).set_index("a")
    df = from_pandas(pdf, npartitions=4, sort=True)
    assert_eq(pdf.sort_index(), df, sort_results=False)

    pdf = pd.DataFrame({"a": [1, 2, 3] * 3, "b": 1}).set_index("a")
    df = from_pandas(pdf, npartitions=4, sort=False)
    assert_eq(pdf, df, sort_results=False)


def test_from_pandas_sort():
    pdf = pd.DataFrame({"a": [1, 2, 3, 1, 2, 2]}, index=[6, 5, 4, 3, 2, 1])
    df = from_pandas(pdf, npartitions=2)
    assert_eq(df, pdf.sort_index(), sort_results=False)


def test_from_pandas_divisions():
    pdf = pd.DataFrame({"a": [1, 2, 3, 1, 2, 2]}, index=[7, 6, 4, 3, 2, 1])
    df = repartition(pdf, (1, 5, 8))
    assert_eq(df, pdf.sort_index())

    pdf = pd.DataFrame({"a": [1, 2, 3, 1, 2, 2]}, index=[7, 6, 4, 3, 2, 1])
    df = repartition(pdf, (1, 4, 8))
    assert_eq(df.partitions[1], pd.DataFrame({"a": [3, 2, 1]}, index=[4, 6, 7]))

    df = repartition(df, divisions=(1, 3, 8), force=True)
    assert_eq(df, pdf.sort_index())


def test_from_pandas_empty_projection():
    pdf = pd.DataFrame({"a": [1, 2, 3], "b": 1})
    df = from_pandas(pdf)
    assert_eq(df[[]], pdf[[]])


def test_from_pandas_divisions_duplicated():
    pdf = pd.DataFrame({"a": 1}, index=[1, 2, 3, 4, 5, 5, 5, 6, 8])
    df = repartition(pdf, (1, 5, 7, 10))
    assert_eq(df, pdf)
    assert_eq(df.partitions[0], pdf.loc[1:4])
    assert_eq(df.partitions[1], pdf.loc[5:6])
    assert_eq(df.partitions[2], pdf.loc[8:])


@pytest.mark.parametrize("npartitions", [1, 3, 6, 7])
@pytest.mark.parametrize("sort", [True, False])
def test_from_pandas_npartitions(pdf, npartitions, sort):
    df = from_pandas(pdf, sort=sort, npartitions=npartitions)
    assert df.npartitions == min(pdf.shape[0], npartitions)
    assert "pandas" in df._name
    assert_eq(df, pdf, sort_results=sort)


@pytest.mark.parametrize("chunksize,npartitions", [(1, 6), (2, 3), (6, 1), (7, 1)])
@pytest.mark.parametrize("sort", [True, False])
def test_from_pandas_chunksize(pdf, chunksize, npartitions, sort):
    df = from_pandas(pdf, sort=sort, chunksize=chunksize)
    assert df.npartitions == npartitions
    assert "pandas" in df._name
    assert_eq(df, pdf, sort_results=sort)


def test_from_pandas_npartitions_and_chunksize(pdf):
    with pytest.raises(ValueError, match="npartitions and chunksize"):
        from_pandas(pdf, npartitions=2, chunksize=3)


def test_from_pandas_string_option():
    pdf = pd.DataFrame({"x": [1, 2, 3], "y": "a"}, index=["a", "b", "c"])
    df = from_pandas(pdf, npartitions=2)
    dtype = "string" if pyarrow_strings_enabled() else "object"
    assert df.dtypes["y"] == dtype
    assert df.index.dtype == dtype
    assert df.compute().dtypes["y"] == dtype
    assert df.compute().index.dtype == dtype
    assert_eq(df, pdf)

    with dask.config.set({"dataframe.convert-string": False}):
        df = from_pandas(pdf, npartitions=2)
        assert df.dtypes["y"] == "object"
        assert df.index.dtype == "object"
        assert df.compute().dtypes["y"] == "object"
        assert df.compute().index.dtype == "object"
        assert_eq(df, pdf)


def test_from_pandas_deepcopy():
    pdf = pd.DataFrame({"col1": [1, 2, 3, 4, 5, 6]})
    df = from_pandas(pdf, npartitions=3)
    df_dict = {"dataset": df}
    result = copy.deepcopy(df_dict)
    assert_eq(result["dataset"], pdf)


def test_from_pandas_empty_chunksize():
    pdf = pd.DataFrame()
    df = from_pandas(pdf, chunksize=10_000)
    assert_eq(pdf, df)
