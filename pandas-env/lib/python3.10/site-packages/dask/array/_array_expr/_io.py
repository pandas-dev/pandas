from __future__ import annotations

import functools

from dask import istask
from dask.array._array_expr._expr import ArrayExpr
from dask.tokenize import _tokenize_deterministic


class IO(ArrayExpr):
    pass


class FromGraph(IO):
    _parameters = ["layer", "_meta", "chunks", "keys", "name_prefix"]

    @functools.cached_property
    def _meta(self):
        return self.operand("_meta")

    @functools.cached_property
    def chunks(self):
        return self.operand("chunks")

    @functools.cached_property
    def _name(self):
        return (
            self.operand("name_prefix") + "-" + _tokenize_deterministic(*self.operands)
        )

    def _layer(self):
        dsk = dict(self.operand("layer"))
        # The name may not actually match the layers name therefore rewrite this
        # using an alias
        for k in self.operand("keys"):
            if not isinstance(k, tuple):
                raise TypeError(f"Expected tuple, got {type(k)}")
            orig = dsk[k]
            if not istask(orig):
                del dsk[k]
                dsk[(self._name, *k[1:])] = orig
            else:
                dsk[(self._name, *k[1:])] = k
        return dsk
