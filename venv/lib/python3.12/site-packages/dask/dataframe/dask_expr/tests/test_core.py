from __future__ import annotations

import pickle
import uuid

import pandas as pd
import pytest

import dask
import dask.dataframe as dd
from dask._expr import Expr
from dask.dataframe.utils import assert_eq
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


def test_combine_multiple_dataframes():
    df = pd.DataFrame({"A": [1, 2, 3, 4], "B": [0, 0, 1, 1]})
    ddf = dd.from_pandas(df, npartitions=2)

    X = ddf[["A", "B"]]
    y = ddf[["B"]]

    x_per, y_per = dask.persist(X, y)
    assert_eq(x_per, X)
    assert_eq(y_per, y)
