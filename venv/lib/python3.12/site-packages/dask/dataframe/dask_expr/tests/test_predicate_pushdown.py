from __future__ import annotations

import pytest

from dask.dataframe.dask_expr import from_pandas
from dask.dataframe.dask_expr._expr import rewrite_filters
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

pd = _backend_library()


@pytest.fixture
def pdf():
    pdf = pd.DataFrame({"x": range(100), "a": 1, "b": 2})
    pdf["y"] = pdf.x // 7  # Not unique; duplicates span different partitions
    yield pdf


@pytest.fixture
def df(pdf):
    yield from_pandas(pdf, npartitions=10)


def test_rewrite_filters(df):
    predicate = (df.x == 1) | (df.x == 1)
    expected = df.x == 1
    assert rewrite_filters(predicate.expr)._name == expected._name

    predicate = (df.x == 1) | ((df.x == 1) & (df.y == 2))
    expected = df.x == 1
    assert rewrite_filters(predicate.expr)._name == expected._name

    predicate = ((df.x == 1) & (df.y == 3)) | ((df.x == 1) & (df.y == 2))
    expected = (df.x == 1) & ((df.y == 3) | (df.y == 2))
    assert rewrite_filters(predicate.expr)._name == expected._name

    predicate = ((df.x == 1) & (df.y == 3) & (df.a == 1)) | (
        (df.x == 1) & (df.y == 2) & (df.a == 1)
    )
    expected = ((df.x == 1) & (df.a == 1)) & ((df.y == 3) | (df.y == 2))
    assert rewrite_filters(predicate.expr)._name == expected._name

    predicate = (df.x == 1) | (df.y == 1)
    assert rewrite_filters(predicate.expr)._name == predicate._name

    predicate = df.x == 1
    assert rewrite_filters(predicate.expr)._name == predicate._name

    predicate = (df.x.isin([1, 2, 3]) & (df.y == 3)) | (
        df.x.isin([1, 2, 3]) & (df.y == 2)
    )
    expected = df.x.isin([1, 2, 3]) & ((df.y == 3) | (df.y == 2))
    assert rewrite_filters(predicate.expr)._name == expected._name


def test_rewrite_filters_query(df, pdf):
    result = df[((df.x == 1) & (df.y > 1)) | ((df.x == 1) & (df.y > 2))]
    result = result[["x"]]
    expected = pdf[((pdf.x == 1) & (pdf.y > 1)) | ((pdf.x == 1) & (pdf.y > 2))]
    expected = expected[["x"]]
    assert_eq(result, expected)
