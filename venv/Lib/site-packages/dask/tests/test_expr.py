from __future__ import annotations

import pytest

from dask._expr import Expr


def test_setattr():
    class MyExpr(Expr):
        _parameters = ["foo", "bar"]

    e = MyExpr(foo=1, bar=2)
    e.bar = 3
    assert e.bar == 3
    with pytest.raises(AttributeError):
        e.baz = 4
