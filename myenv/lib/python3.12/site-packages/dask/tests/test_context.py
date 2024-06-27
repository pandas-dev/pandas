from __future__ import annotations

import pytest

import dask
from dask.context import globalmethod


def test_with_get():
    pytest.importorskip("numpy")
    da = pytest.importorskip("dask.array")
    var = [0]

    def myget(dsk, keys, **kwargs):
        var[0] = var[0] + 1
        return dask.get(dsk, keys, **kwargs)

    x = da.ones(10, chunks=(5,))

    assert x.sum().compute() == 10
    assert var[0] == 0

    with dask.config.set(scheduler=myget):
        assert x.sum().compute() == 10
    assert var[0] == 1

    # Make sure we've cleaned up
    assert x.sum().compute() == 10
    assert var[0] == 1


def foo():
    return "foo"


def bar():
    return "bar"


class Foo:
    @globalmethod(key="f")
    def f():  # type: ignore
        return 1

    g = globalmethod(foo, key="g", falsey=bar)


def test_globalmethod():
    x = Foo()

    assert x.f() == 1

    with dask.config.set(f=lambda: 2):
        assert x.f() == 2

    with dask.config.set(f=foo):
        assert x.f is foo
        assert x.f() == "foo"

    assert x.g is foo
    assert x.g() == "foo"

    with dask.config.set(g=False):
        assert x.g is bar
        assert x.g() == "bar"
