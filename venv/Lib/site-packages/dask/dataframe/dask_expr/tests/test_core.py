from __future__ import annotations

import pickle
import uuid

import pytest

from dask._expr import Expr
from dask.tokenize import tokenize


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


class NondeterminisitcToken:
    def __dask_tokenize__(self):
        return uuid.uuid4()


def test_expr_nondeterministic_token_pickle_roundtrip():

    foo = NondeterminisitcToken()
    assert tokenize(foo) != tokenize(foo)
    assert tokenize(Expr(foo)) != tokenize(Expr(foo))

    inst = Expr(foo)
    tok = tokenize(inst)
    assert tokenize(inst) == tok
    rt = pickle.loads(pickle.dumps(inst))
    assert tokenize(rt) == tok
