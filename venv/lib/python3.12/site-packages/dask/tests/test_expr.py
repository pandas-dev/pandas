from __future__ import annotations

import functools
import pickle

import pytest

from dask._expr import Expr, ProhibitReuse, SingletonExpr, _ExprSequence
from dask._task_spec import DataNode


class MyExpr(Expr):
    _parameters = ["foo", "bar"]


class MyExpr2(MyExpr):
    # A subclass that inherits parameters
    pass


def test_setattr():
    e = MyExpr(foo=1, bar=2)
    e.bar = 3
    assert e.bar == 3
    with pytest.raises(AttributeError):
        e.baz = 4


def test_setattr2():
    e = MyExpr2(foo=1, bar=2)
    e.bar = 3
    assert e.bar == 3
    with pytest.raises(AttributeError):
        e.baz = 4


class MyExprCachedProperty(Expr):
    called_cached_property = False
    _parameters = ["foo", "bar"]

    @property
    def baz(self):
        return self.foo + self.bar

    @functools.cached_property
    def cached_property(self):
        if MyExprCachedProperty.called_cached_property:
            raise RuntimeError("No!")
        MyExprCachedProperty.called_cached_property = True
        return self.foo + self.bar


@pytest.mark.slow()
def test_pickle_cached_properties():
    pytest.importorskip("distributed")
    from distributed import Nanny
    from distributed.utils_test import gen_cluster

    @gen_cluster(client=True, Worker=Nanny, nthreads=[("", 1)])
    async def test(c, s, a):

        expr = MyExprCachedProperty(foo=1, bar=2)
        for _ in range(10):
            assert expr.baz == 3
            assert expr.cached_property == 3

        assert MyExprCachedProperty.called_cached_property is True

        rt = pickle.loads(pickle.dumps(expr))
        assert rt.cached_property == 3
        assert MyExprCachedProperty.called_cached_property is True

        # But this does
        expr3 = MyExprCachedProperty(foo=1, bar=3)
        with pytest.raises(RuntimeError):
            expr3.cached_property

        def f(expr):
            # We want the cache to be part of the pickle, i.e. this is a
            # different process such that the type is reset and the property can
            # be accessed without side effects
            assert MyExprCachedProperty.called_cached_property is False
            assert expr.cached_property == 3
            assert MyExprCachedProperty.called_cached_property is False

        await c.submit(f, expr)

    test()


class MySingleton(SingletonExpr): ...


class MySingletonWithCustomInit(SingletonExpr):
    def __init__(self, *args, **kwargs): ...


class MySingletonInheritsCustomInit(MySingletonWithCustomInit): ...


class Mixin:
    def __init__(self, *args, **kwargs): ...


class MySingletonInheritsCustomInitAsMixin(MySingleton, Mixin): ...


def test_singleton_expr():
    assert MySingleton(1, 2) is MySingleton(1, 2)
    # We don't want to deduplicate if there is an __init__ since that may
    # mutatate our singleton reference and we have no way to know
    assert MySingletonWithCustomInit(1, 2) is not MySingletonWithCustomInit(1, 2)
    assert MySingletonInheritsCustomInit(1, 2) is not MySingletonInheritsCustomInit(
        1, 2
    )
    assert MySingletonInheritsCustomInitAsMixin(
        1, 2
    ) is not MySingletonInheritsCustomInitAsMixin(1, 2)


@pytest.mark.slow()
def test_refcounting_futures():
    pd = pytest.importorskip("pandas")
    dd = pytest.importorskip("dask.dataframe")
    distributed = pytest.importorskip("distributed")

    # See https://github.com/dask/distributed/issues/9041
    # Didn't reproduce with any of our fixtures
    with distributed.Client(
        n_workers=2, worker_class=distributed.Worker, dashboard_address=":0"
    ) as client:

        def gen(i):
            return pd.DataFrame({"A": [i]}, index=[i])

        futures = [client.submit(gen, i) for i in range(3)]

        meta = gen(0)[:0]
        df = dd.from_delayed(futures, meta)
        df.compute()

        del futures

        df.compute()


class FooExpr(Expr):
    def _layer(self) -> dict:
        return {"foo": DataNode("foo", 42)}


def test_prohibit_reuse():
    once = FooExpr()
    ProhibitReuse._ALLOWED_TYPES.append(FooExpr)
    try:
        dsk = _ExprSequence(once, ProhibitReuse(once)).optimize().__dask_graph__()

        assert len(dsk) == 2
        first = dsk.pop("foo")()
        key, val = dsk.popitem()
        assert key.startswith("foo") and key != "foo"
        # We don't want to chain anything but actually _hide_ the task
        assert not val.dependencies
        # Task is wrapped
        assert val() is first
    finally:
        ProhibitReuse._ALLOWED_TYPES.remove(FooExpr)
