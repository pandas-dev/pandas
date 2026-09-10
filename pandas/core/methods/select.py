"""
Implementation of DataFrame.select.
"""

from __future__ import annotations

from typing import (
    TYPE_CHECKING,
)

from pandas.core.dtypes.common import is_list_like
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCSeries,
)

from pandas.core import common as com
from pandas.core.col import Expression

if TYPE_CHECKING:
    from collections.abc import Hashable

    from pandas import DataFrame


def _parse_select_args(args: tuple) -> list:
    """
    Parse the positional arguments of DataFrame.select.

    A single list-like argument (other than a tuple, which is a single
    column label) is treated as the sequence of items to select; otherwise
    each argument is one item. Mixing a sequence with other positional
    arguments raises, as do sets, whose iteration order is undefined.
    """

    def is_sequence(arg: object) -> bool:
        if isinstance(arg, (set, frozenset)):
            raise TypeError(
                "`select` does not support sets, as the resulting column "
                "order would be undefined; pass a list instead"
            )
        return not isinstance(arg, (tuple, Expression)) and is_list_like(arg)

    if len(args) == 1 and is_sequence(args[0]):
        return list(args[0])
    if any(is_sequence(arg) for arg in args):
        raise TypeError(
            "`select` supports individual columns "
            "`df.select('col1', 'col2', ...)` or a single list "
            "`df.select(['col1', 'col2', ...])`, but not both. "
            "You can unpack the list if you have a mix: "
            "`df.select(*['col1', 'col2'], 'col3')`."
        )
    return list(args)


def select(
    df: DataFrame, select_args: tuple, select_kwargs: dict[str, object]
) -> DataFrame:
    """
    Implementation of DataFrame.select.
    """
    from pandas.core.reshape.concat import concat

    items = _parse_select_args(select_args)
    nlevels = df.columns.nlevels

    # Items are resolved in order against ``data``. Computed columns are
    # written into ``data`` by name, as in DataFrame.assign, so later
    # items can refer to them. The result is built from ``chunks`` so
    # that a name requested more than once is returned more than once.
    data = df.copy(deep=False)
    chunks: list[DataFrame] = []
    labels: list[Hashable] = []

    def flush_labels() -> None:
        if labels:
            indexer = data.columns._get_indexer_strict(labels, "columns")[1]
            chunks.append(data.take(indexer, axis=1))
            labels.clear()

    def set_column(name: Hashable, value: object) -> None:
        data[name] = value
        loc = data.columns.get_loc(name)
        if isinstance(loc, slice):
            # a duplicated label; every occurrence now holds ``value``,
            # and a computed column contributes one column to the result
            loc = loc.start
        elif not isinstance(loc, int):
            # boolean mask for a non-contiguous duplicated label
            loc = int(loc.argmax())
        chunks.append(data.iloc[:, loc : loc + 1])

    for item in items:
        if isinstance(item, Expression):
            flush_labels()
            result = item._eval_expression(data)
            if isinstance(result, ABCSeries):
                if result.name is None:
                    raise TypeError(
                        f"expression {item!r} evaluated to an unnamed "
                        "Series; use .rename(...) to name it or pass it "
                        "as a keyword argument"
                    )
                if nlevels > 1 and (
                    not isinstance(result.name, tuple) or len(result.name) != nlevels
                ):
                    raise TypeError(
                        f"expression {item!r} evaluated to a Series "
                        f"named {result.name!r}, but the DataFrame's "
                        f"columns have {nlevels} levels; use "
                        f".rename(...) with a tuple of length {nlevels}"
                    )
                set_column(result.name, result)
            elif isinstance(result, ABCDataFrame):
                if result.columns.nlevels != nlevels:
                    raise TypeError(
                        f"expression {item!r} evaluated to a DataFrame "
                        f"whose columns have {result.columns.nlevels} "
                        f"levels, but the DataFrame's columns have "
                        f"{nlevels} levels"
                    )
                for i, name in enumerate(result.columns):
                    set_column(name, result.iloc[:, i])
            else:
                raise TypeError(
                    f"expression {item!r} evaluated to an object of type "
                    f"{type(result).__name__}, expected a Series or "
                    "DataFrame; pass it as a keyword argument instead"
                )
        else:
            if (
                nlevels > 1
                and labels
                and isinstance(item, tuple) != isinstance(labels[-1], tuple)
            ):
                # on a MultiIndex, single-level keys and full tuples
                # cannot be resolved in one indexer batch
                flush_labels()
            labels.append(item)
    flush_labels()

    for key, value in select_kwargs.items():
        if nlevels > 1:
            raise TypeError(
                f"cannot create column {key!r} from a keyword argument: "
                f"the DataFrame's columns have {nlevels} levels, so new "
                f"column names must be tuples of length {nlevels}; pass "
                "a positional expression with .rename(...) instead"
            )
        set_column(key, com.apply_if_callable(value, data))

    if not chunks:
        return df.iloc[:, :0]
    return concat(chunks, axis=1)
