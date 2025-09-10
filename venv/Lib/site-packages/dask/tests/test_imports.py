from __future__ import annotations

import importlib

import pytest


@pytest.fixture(autouse=True)
def invalidate_caches():
    importlib.invalidate_caches()


@pytest.mark.parametrize(
    "module",
    [
        "dask",
        "dask.bag",
        "dask.base",
        "dask.delayed",
        "dask.graph_manipulation",
        "dask.layers",
        "dask.multiprocessing",
        "dask.optimization",
        "dask.threaded",
    ],
)
def test_defaults(module):
    __import__(module)


def test_array():
    pytest.importorskip("numpy")
    import dask.array  # noqa: F401


def test_pandas_pyarrow():
    pytest.importorskip("pandas")
    pytest.importorskip("pyarrow")
    import dask.dataframe  # noqa: F401


def test_bokeh():
    pytest.importorskip("bokeh")
    import dask.diagnostics  # noqa: F401


def test_distributed():
    pytest.importorskip("distributed")
    import dask.distributed  # noqa: F401
