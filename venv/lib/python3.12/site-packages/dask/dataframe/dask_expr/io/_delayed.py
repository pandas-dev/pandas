from __future__ import annotations

import functools
from collections.abc import Collection
from typing import TYPE_CHECKING

import pandas as pd

from dask._task_spec import Alias, Task, TaskRef
from dask.dataframe.dask_expr._expr import (
    ArrowStringConversion,
    DelayedsExpr,
    PartitionsFiltered,
)
from dask.dataframe.dask_expr.io import BlockwiseIO
from dask.dataframe.dispatch import make_meta
from dask.dataframe.utils import check_meta, pyarrow_strings_enabled
from dask.delayed import Delayed, delayed
from dask.typing import Key

if TYPE_CHECKING:
    import distributed


class FromDelayed(PartitionsFiltered, BlockwiseIO):
    _parameters = [
        "delayed_container",
        "meta",
        "user_divisions",
        "verify_meta",
        "_partitions",
        "prefix",
    ]
    _defaults = {
        "meta": None,
        "_partitions": None,
        "user_divisions": None,
        "verify_meta": True,
        "prefix": None,
    }

    @functools.cached_property
    def _name(self):
        if self.prefix is None:
            return super()._name
        return self.prefix + "-" + self.deterministic_token

    @functools.cached_property
    def _meta(self):
        if self.operand("meta") is not None:
            return self.operand("meta")

        return delayed(make_meta)(self.delayed_container.operands[0]).compute()

    def _divisions(self):
        if self.operand("user_divisions") is not None:
            return self.operand("user_divisions")
        else:
            return self.delayed_container.divisions

    def _filtered_task(self, name: Key, index: int) -> Task:
        if self.verify_meta:
            return Task(
                name,
                functools.partial(check_meta, meta=self._meta, funcname="from_delayed"),
                TaskRef((self.delayed_container._name, index)),
                _data_producer=True,
            )
        else:
            return Alias((self.delayed_container._name, index))  # type: ignore


def identity(x):
    return x


def from_delayed(
    dfs: Delayed | distributed.Future | Collection[Delayed | distributed.Future],
    meta=None,
    divisions: tuple | None = None,
    prefix: str | None = None,
    verify_meta: bool = True,
):
    """Create Dask DataFrame from many Dask Delayed objects

    .. warning::
        ``from_delayed`` should only be used if the objects that create
        the data are complex and cannot be easily represented as a single
        function in an embarrassingly parallel fashion.

        ``from_map`` is recommended if the query can be expressed as a single
        function like:

        def read_xml(path):
            return pd.read_xml(path)

        ddf = dd.from_map(read_xml, paths)

        ``from_delayed`` might be deprecated in the future.

    Parameters
    ----------
    dfs :
        A ``dask.delayed.Delayed``, a ``distributed.Future``, or an iterable of either
        of these objects, e.g. returned by ``client.submit``. These comprise the
        individual partitions of the resulting dataframe.
        If a single object is provided (not an iterable), then the resulting dataframe
        will have only one partition.
    $META
    divisions :
        Partition boundaries along the index.
        For tuple, see https://docs.dask.org/en/latest/dataframe-design.html#partitions
        If None, then won't use index information
    prefix :
        Prefix to prepend to the keys.
    verify_meta :
        If True check that the partitions have consistent metadata, defaults to True.
    """
    if isinstance(dfs, Delayed) or hasattr(dfs, "key"):
        dfs = [dfs]

    if len(dfs) == 0:
        raise TypeError("Must supply at least one delayed object")

    if meta is None:
        meta = delayed(make_meta)(dfs[0]).compute()  # type: ignore

    if divisions == "sorted":
        raise NotImplementedError(
            "divisions='sorted' not supported, please calculate the divisions "
            "yourself."
        )
    elif divisions is not None:
        divs = list(divisions)
        if len(divs) != len(dfs) + 1:
            raise ValueError("divisions should be a tuple of len(dfs) + 1")

    futures = [v for v in dfs if isinstance(v, TaskRef)]
    if len(futures) == len(dfs):
        # All futures. Fast path
        dfs = futures
    else:
        # Every Delayed generates a Layer, i.e. this path is much more expensive
        # if there are many input values.
        dfs = [
            delayed(v) if not isinstance(v, (Delayed,)) and hasattr(v, "key") else v
            for v in dfs
        ]

    for item in dfs:
        if not (isinstance(item, (Delayed, TaskRef))):
            raise TypeError("Expected Delayed object, got %s" % type(item).__name__)

    from dask.dataframe.dask_expr._collection import new_collection

    result = FromDelayed(
        DelayedsExpr(*dfs), make_meta(meta), divisions, verify_meta, None, prefix
    )
    if pyarrow_strings_enabled() and any(
        pd.api.types.is_object_dtype(dtype)
        for dtype in (result.dtypes.values if result.ndim == 2 else [result.dtypes])
    ):
        return new_collection(ArrowStringConversion(result))
    return new_collection(result)
