"""
The ``ax`` accessor: apply a length-preserving Index method to an axis of a
Series or DataFrame and return a new object with that axis replaced.
"""

from __future__ import annotations

import functools
from typing import (
    TYPE_CHECKING,
    Any,
)

from pandas.core.indexes.api import Index

if TYPE_CHECKING:
    from collections.abc import Callable

    from pandas.core.generic import NDFrame


# may permute the labels, which the length check in _apply can't catch
_REORDERING_METHODS = frozenset(
    {
        "sort_values",
        "take",
        "union",
        "intersection",
        "difference",
        "symmetric_difference",
        "join",
    }
)


class AxisAccessor:
    """
    Pass-through to the methods of a Series or DataFrame axis.

    ``obj.ax.method(*args, axis=0, **kwargs)`` is equivalent to
    ``obj.set_axis(obj.index.method(*args, **kwargs), axis=0)``, with
    ``axis=1`` using ``obj.columns`` instead; ``axis`` is keyword-only. Methods
    that do not return an Index of the same length as the axis raise TypeError.
    Methods that can reorder the labels without moving the data, such as
    ``sort_values``, ``take`` and the set operations, are not supported.

    Parameters
    ----------
    data : Series or DataFrame
        The object whose axes the methods are applied to.

    See Also
    --------
    DataFrame.set_axis : Assign desired index to given axis.
    DataFrame.rename_axis : Set the name of the axis for the index or columns.

    Examples
    --------
    >>> df = pd.DataFrame(
    ...     {"a": [1, 2]}, index=pd.MultiIndex.from_tuples([("x", 1), ("y", 2)])
    ... )
    >>> df.ax.to_flat_index()
            a
    (x, 1)  1
    (y, 2)  2

    >>> df.ax.str.upper(axis=1)
         A
    x 1  1
    y 2  2
    """

    # non-empty for a nested Index accessor, e.g. ("str",) for ``obj.ax.str``
    _path: tuple[str, ...] = ()

    def __init__(self, data: NDFrame) -> None:
        self._obj = data

    def _sources(self) -> list[Any]:
        axes = self._obj.axes
        if not self._path:
            return [type(labels) for labels in axes]
        name = self._path[0]
        return [
            getattr(type(labels), name)
            for labels in axes
            if name in labels._accessors and hasattr(labels, name)
        ]

    def __getattr__(self, name: str) -> Callable[..., NDFrame] | AxisAccessor:
        if name.startswith("_"):
            raise AttributeError(name)
        if not self._path and name in _REORDERING_METHODS:
            raise AttributeError(
                f"ax does not support '{name}', which can reorder the labels "
                "without moving the data"
            )
        axes = self._obj.axes
        if not self._path and any(name in labels._accessors for labels in axes):
            if not any(hasattr(labels, name) for labels in axes):
                # e.g. ``.str`` on a datetime axis
                raise AttributeError(
                    f"no {type(self._obj).__name__} axis supports the '{name}' accessor"
                )
            nested = AxisAccessor(self._obj)
            nested._path = (name,)
            return nested

        path = (*self._path, name)
        source = next(
            (
                getattr(src, name)
                for src in self._sources()
                if callable(getattr(src, name, None))
            ),
            None,
        )
        if source is None:
            raise AttributeError(
                f"no {type(self._obj).__name__} axis has an Index method "
                f"'{'.'.join(path)}'"
            )

        def method(*args: Any, axis: Any = 0, **kwargs: Any) -> NDFrame:
            return _apply(self._obj, axis, path, args, kwargs)

        method.__name__ = name
        method.__doc__ = source.__doc__
        return method

    def __dir__(self) -> list[str]:
        names = set().union(*(_public_methods(src) for src in self._sources()))
        if not self._path:
            names |= {
                name
                for labels in self._obj.axes
                for name in labels._accessors
                if hasattr(labels, name)
            }
            names -= _REORDERING_METHODS
        return sorted(names)


@functools.cache
def _public_methods(cls: type) -> frozenset[str]:
    # accessor classes are callable too; __dir__ adds the valid ones back
    accessors = getattr(cls, "_accessors", ())
    return frozenset(
        name
        for name in dir(cls)
        if not name.startswith("_")
        and name not in accessors
        and callable(getattr(cls, name))
    )


def _apply(
    obj: NDFrame,
    axis: Any,
    path: tuple[str, ...],
    args: tuple[Any, ...],
    kwargs: dict[str, Any],
) -> NDFrame:
    if kwargs.get("inplace"):
        # would mutate the original axis before the result check could raise
        raise ValueError("ax does not support inplace=True")
    axis_number = obj._get_axis_number(axis)
    labels = obj._get_axis(axis_number)

    func: Any = labels
    for name in path:
        func = getattr(func, name)
    qualname = ".".join(path)

    result = func(*args, **kwargs)
    if not isinstance(result, Index) or len(result) != len(labels):
        raise TypeError(
            f"{type(labels).__name__}.{qualname} did not return an Index of the "
            "same length as the axis; ax only supports length-preserving Index "
            "methods"
        )
    return obj.set_axis(result, axis=axis_number)
