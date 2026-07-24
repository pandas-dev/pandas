from __future__ import annotations

import functools

from dask.dataframe.dask_expr._accessor import Accessor, FunctionMap
from dask.dataframe.dask_expr._expr import Blockwise
from dask.dataframe.dask_expr._reductions import Reduction
from dask.dataframe.dispatch import make_meta, meta_nonempty


class StringAccessor(Accessor):
    """Accessor object for string properties of the Series values.

    Examples
    --------

    >>> s.str.lower()  # doctest: +SKIP
    """

    _accessor_name = "str"

    _accessor_methods = (
        "capitalize",
        "casefold",
        "center",
        "contains",
        "count",
        "decode",
        "encode",
        "endswith",
        "extract",
        "extractall",
        "find",
        "findall",
        "fullmatch",
        "get",
        "index",
        "isalnum",
        "isalpha",
        "isdecimal",
        "isdigit",
        "islower",
        "isnumeric",
        "isspace",
        "istitle",
        "isupper",
        "join",
        "len",
        "ljust",
        "lower",
        "lstrip",
        "match",
        "normalize",
        "pad",
        "partition",
        "removeprefix",
        "removesuffix",
        "repeat",
        "replace",
        "rfind",
        "rindex",
        "rjust",
        "rpartition",
        "rstrip",
        "slice",
        "slice_replace",
        "startswith",
        "strip",
        "swapcase",
        "title",
        "translate",
        "upper",
        "wrap",
        "zfill",
    )
    _accessor_properties = ()

    def _split(self, method, pat=None, n=-1, expand=False):
        from dask.dataframe.dask_expr import new_collection

        if expand:
            if n == -1:
                raise NotImplementedError(
                    "To use the expand parameter you must specify the number of "
                    "expected splits with the n= parameter. Usually n splits "
                    "result in n+1 output columns."
                )
            return new_collection(
                SplitMap(
                    self._series,
                    self._accessor_name,
                    method,
                    (),
                    {"pat": pat, "n": n, "expand": expand},
                )
            )
        return self._function_map(method, pat=pat, n=n, expand=expand)

    def split(self, pat=None, n=-1, expand=False):
        """Known inconsistencies: ``expand=True`` with unknown ``n`` will raise a ``NotImplementedError``."""
        return self._split("split", pat=pat, n=n, expand=expand)

    def rsplit(self, pat=None, n=-1, expand=False):
        return self._split("rsplit", pat=pat, n=n, expand=expand)

    def cat(self, others=None, sep=None, na_rep=None):
        import pandas as pd

        from dask.dataframe.dask_expr._collection import Index, Series, new_collection

        if others is None:
            return new_collection(Cat(self._series, sep, na_rep))

        valid_types = (Series, Index, pd.Series, pd.Index)
        if isinstance(others, valid_types):
            others = [others]
        elif not all(isinstance(a, valid_types) for a in others):
            raise TypeError("others must be Series/Index")

        return new_collection(CatBlockwise(self._series, sep, na_rep, *others))

    def __getitem__(self, index):
        return self._function_map("__getitem__", index)


class CatBlockwise(Blockwise):
    _parameters = ["frame", "sep", "na_rep"]
    _keyword_only = ["sep", "na_rep"]

    @property
    def _args(self) -> list:
        return [self.frame] + self.operands[len(self._parameters) :]

    @staticmethod
    def operation(ser, *args, **kwargs):
        return ser.str.cat(list(args), **kwargs)


class Cat(Reduction):
    _parameters = ["frame", "sep", "na_rep"]

    @property
    def chunk_kwargs(self):
        return {"sep": self.sep, "na_rep": self.na_rep}

    @property
    def combine_kwargs(self):
        return self.chunk_kwargs

    @property
    def aggregate_kwargs(self):
        return self.chunk_kwargs

    @staticmethod
    def reduction_chunk(ser, *args, **kwargs):
        return ser.str.cat(*args, **kwargs)

    @staticmethod
    def reduction_combine(ser, *args, **kwargs):
        return Cat.reduction_chunk(ser, *args, **kwargs)

    @staticmethod
    def reduction_aggregate(ser, *args, **kwargs):
        return Cat.reduction_chunk(ser, *args, **kwargs)


class SplitMap(FunctionMap):
    _parameters = ["frame", "accessor", "attr", "args", "kwargs"]

    @property
    def n(self):
        return self.kwargs["n"]

    @property
    def pat(self):
        return self.kwargs["pat"]

    @property
    def expand(self):
        return self.kwargs["expand"]

    @functools.cached_property
    def _meta(self):
        delimiter = " " if self.pat is None else self.pat
        meta = meta_nonempty(self.frame._meta)
        meta = self.frame._meta._constructor(
            [delimiter.join(["a"] * (self.n + 1))],
            index=meta.iloc[:1].index,
        )
        return make_meta(
            getattr(meta.str, self.attr)(n=self.n, expand=self.expand, pat=self.pat)
        )
