from __future__ import annotations

from operator import add
from time import sleep

import pytest

from dask.cache import Cache
from dask.callbacks import Callback
from dask.local import get_sync
from dask.threaded import get

cachey = pytest.importorskip("cachey")


flag = []


def inc(x):
    flag.append(x)
    return x + 1


def test_cache():
    c = cachey.Cache(10000)
    cc = Cache(c)

    with cc:
        assert get({"x": (inc, 1)}, "x") == 2

    assert flag == [1]
    assert c.data["x"] == 2

    assert not cc.starttimes
    assert not cc.durations

    while flag:
        flag.pop()
    dsk = {"x": (inc, 1), "y": (inc, 2), "z": (add, "x", "y")}
    with cc:
        assert get(dsk, "z") == 5

    assert flag == [2]  # no x present

    assert not Callback.active


def test_cache_with_number():
    c = Cache(10000, limit=1)
    assert isinstance(c.cache, cachey.Cache)
    assert c.cache.available_bytes == 10000
    assert c.cache.limit == 1


def test_cache_correctness():
    # https://github.com/dask/dask/issues/3631
    c = Cache(10000)
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    from numpy import ones, zeros

    z = da.from_array(zeros(1), chunks=10)
    o = da.from_array(ones(1), chunks=10)
    with c:
        assert (z.compute() == 0).all()
        assert (o.compute() == 1).all()


def f(duration, size, *args):
    sleep(duration)
    return [0] * size


def test_prefer_cheap_dependent():
    dsk = {"x": (f, 0.01, 10), "y": (f, 0.000001, 1, "x")}
    c = Cache(10000)
    with c:
        get_sync(dsk, "y")

    assert c.cache.scorer.cost["x"] < c.cache.scorer.cost["y"]
