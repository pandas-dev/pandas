from __future__ import annotations

import pytest

from dask._expr import Expr


class ExprB(Expr):
    def _simplify_down(self):
        return ExprA()

    def _operands_for_repr(self):
        return []


class ExprA(Expr):
    def _simplify_down(self):
        return ExprB()

    def _operands_for_repr(self):
        return []


def test_endless_simplify():
    expr = ExprA()
    with pytest.raises(RuntimeError, match="converge"):
        expr.simplify()
