from __future__ import annotations

import numpy as np
import pytest

import dask.array as da
from dask.dataframe.dask_expr import Index, from_pandas
from dask.dataframe.dask_expr.tests._util import _backend_library, assert_eq

pd = _backend_library()


@pytest.fixture
def pdf():
    pdf = pd.DataFrame({"x": range(100)})
    pdf["y"] = pdf.x // 7  # Not unique; duplicates span different partitions
    yield pdf


@pytest.fixture
def df(pdf):
    yield from_pandas(pdf, npartitions=10)


def test_ufunc(df, pdf):
    ufunc = "conj"
    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    pandas_type = pdf.__class__
    dask_type = df.__class__

    assert isinstance(dafunc(df), dask_type)
    assert_eq(dafunc(df), npfunc(pdf))

    if isinstance(npfunc, np.ufunc):
        assert isinstance(npfunc(df), dask_type)
    else:
        assert isinstance(npfunc(df), pandas_type)
    assert_eq(npfunc(df), npfunc(pdf))

    # applying Dask ufunc to normal Series triggers computation
    assert isinstance(dafunc(pdf), pandas_type)
    assert_eq(dafunc(df), npfunc(pdf))

    assert isinstance(dafunc(df.index), Index)
    assert_eq(dafunc(df.index), npfunc(pdf.index))

    if isinstance(npfunc, np.ufunc):
        assert isinstance(npfunc(df.index), Index)
    else:
        assert isinstance(npfunc(df.index), pd.Index)

    assert_eq(npfunc(df.index), npfunc(pdf.index))


@pytest.mark.parametrize(
    "ufunc",
    [
        np.mean,
        np.std,
        np.sum,
        np.cumsum,
        np.cumprod,
        np.var,
        np.min,
        np.max,
        np.all,
        np.any,
        np.prod,
    ],
)
def test_reducers(pdf, df, ufunc):
    assert_eq(ufunc(pdf, axis=0), ufunc(df, axis=0))


def test_clip(pdf, df):
    assert_eq(np.clip(pdf, 10, 20), np.clip(df, 10, 20))
