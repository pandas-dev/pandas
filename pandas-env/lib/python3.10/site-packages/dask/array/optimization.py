from __future__ import annotations

from collections.abc import Callable
from itertools import zip_longest
from numbers import Integral
from typing import Any

import numpy as np

from dask import config
from dask._task_spec import (
    Alias,
    DataNode,
    GraphNode,
    Task,
    _execute_subgraph,
    convert_legacy_graph,
    fuse_linear_task_spec,
)
from dask.array.chunk import getitem
from dask.array.core import getter, getter_inline, getter_nofancy
from dask.blockwise import fuse_roots, optimize_blockwise
from dask.core import flatten
from dask.highlevelgraph import HighLevelGraph
from dask.utils import ensure_dict

# All get* functions the optimizations know about


GETTERS = (getter, getter_nofancy, getter_inline, getitem)
# These get* functions aren't ever completely removed from the graph,
# even if the index should be a no-op by numpy semantics. Some array-like's
# don't completely follow semantics, making indexing always necessary.
GETNOREMOVE = (getter, getter_nofancy)


def optimize(dsk, keys, **kwargs):
    """Optimize dask for array computation

    1.  Cull tasks not necessary to evaluate keys
    2.  Perform linear fusion
    """
    if not isinstance(keys, (list, set)):
        keys = [keys]
    keys = list(flatten(keys))

    if not isinstance(dsk, HighLevelGraph):
        dsk = HighLevelGraph.from_collections(id(dsk), dsk, dependencies=())

    dsk = optimize_blockwise(dsk, keys=keys)
    dsk = fuse_roots(dsk, keys=keys)
    dsk = dsk.cull(set(keys))

    # Perform low-level fusion unless the user has
    # specified False explicitly.
    if config.get("optimization.fuse.active") is False:
        return dsk

    dsk = ensure_dict(dsk)

    dsk = convert_legacy_graph(dsk)
    dsk = fuse_linear_task_spec(dsk, keys=keys)
    dsk = _optimize_slices(dsk)

    return dsk


def _is_getter_task(
    value: GraphNode,
) -> tuple[Callable, Any, Any, bool, bool | None] | None:
    """Check if a value in a Dask graph looks like a getter.
    1. Is it a tuple with the first element a known getter.
    If a getter is found, it returns a tuple with (getter, array, index, asarray, lock).
    Otherwise it returns ``None``.
    """
    if not isinstance(value, Task):
        return None
    func = value.func
    if func in GETTERS:
        get = func
    else:
        return None

    length = len(value.args)
    if length == 2:
        # getter defaults to asarray=True, getitem is semantically False
        return get, value.args[0], value.args[1], get is not getitem, None
    elif length == 4:
        return get, *list(value.args)

    return None


def _optimize_slices(dsk):
    """Optimize slices
    1.  Fuse repeated slices, like x[5:][2:6] -> x[7:11]

    This is generally not very important since we are fusing those tasks anyway. There
    is one specific exception to how xarray implements opening netcdf files and subsequent
    slices. Not merging them together can cause reading the whole netcdf file before
    we drop the unnecessary data. Fusing slices avoids that pattern.

    See https://github.com/pydata/xarray/issues/9926
    """
    fancy_ind_types = (list, np.ndarray)
    dsk = dsk.copy()
    for _, v in dsk.items():
        if not (isinstance(v, Task) and v.func is _execute_subgraph):
            continue

        inner_graph: dict = v.args[0]  # type: ignore

        seen = set()

        for inner_key, inner_value in inner_graph.items():
            if inner_key in seen:
                continue
            else:
                seen.add(inner_key)

            if a_task := _is_getter_task(inner_value):
                temp_key = inner_key
                get, a, a_index, a_asarray, a_lock = a_task
                fused = False

                while a.key in inner_graph and (
                    b_task := _is_getter_task(inner_graph[a.key])
                ):
                    seen.add(a.key)
                    f2, b, b_index, b_asarray, b_lock = b_task

                    if isinstance(a_index, DataNode):
                        a_index = a_index.value
                    if isinstance(b_index, DataNode):
                        b_index = b_index.value

                    if a_lock and a_lock is not b_lock:
                        break
                    if (type(a_index) is tuple) != (type(b_index) is tuple):
                        break
                    if type(a_index) is tuple:
                        indices = b_index + a_index
                        if len(a_index) != len(b_index) and any(
                            i is None for i in indices
                        ):
                            break
                        if f2 is getter_nofancy and not any(
                            isinstance(i, fancy_ind_types) for i in indices
                        ):
                            break
                    elif f2 is getter_nofancy and (
                        type(a_index) in fancy_ind_types
                        or type(b_index) in fancy_ind_types
                    ):
                        break
                    try:
                        c_index = fuse_slice(b_index, a_index)
                        # rely on fact that nested gets never decrease in
                        # strictness e.g. `(getter_nofancy, (getter, ...))` never
                        # happens
                        get = getter if f2 is getter_inline else f2

                        inner_graph[temp_key] = Alias(temp_key, a.key)
                        temp_key = a.key
                        fused = True
                    except NotImplementedError:
                        break
                    a, a_index, a_lock = b, c_index, b_lock
                    a_asarray |= b_asarray

                    if not isinstance(a, Task):
                        break

                if not fused:
                    pass
                elif get is getitem or (a_asarray and not a_lock):
                    # default settings are fine, drop the extra parameters Since we
                    # always fallback to inner `get` functions, `get is getitem`
                    # can only occur if all gets are getitem, meaning all
                    # parameters must be getitem defaults.
                    inner_graph[temp_key] = Task(temp_key, get, a, a_index)
                else:
                    inner_graph[temp_key] = Task(
                        temp_key, get, a, a_index, a_asarray, a_lock
                    )

    return dsk


def normalize_slice(s):
    """Replace Nones in slices with integers

    >>> normalize_slice(slice(None, None, None))
    slice(0, None, 1)
    """
    start, stop, step = s.start, s.stop, s.step
    if start is None:
        start = 0
    if step is None:
        step = 1
    if start < 0 or step < 0 or stop is not None and stop < 0:
        raise NotImplementedError()
    return slice(start, stop, step)


def check_for_nonfusible_fancy_indexing(fancy, normal):
    # Check for fancy indexing and normal indexing, where the fancy
    # indexed dimensions != normal indexed dimensions with integers. E.g.:
    # disallow things like:
    # x[:, [1, 2], :][0, :, :] -> x[0, [1, 2], :] or
    # x[0, :, :][:, [1, 2], :] -> x[0, [1, 2], :]
    for f, n in zip_longest(fancy, normal, fillvalue=slice(None)):
        if type(f) is not list and isinstance(n, Integral):
            raise NotImplementedError(
                "Can't handle normal indexing with "
                "integers and fancy indexing if the "
                "integers and fancy indices don't "
                "align with the same dimensions."
            )


def fuse_slice(a, b):
    """Fuse stacked slices together

    Fuse a pair of repeated slices into a single slice:

    >>> fuse_slice(slice(1000, 2000), slice(10, 15))
    slice(1010, 1015, None)

    This also works for tuples of slices

    >>> fuse_slice((slice(100, 200), slice(100, 200, 10)),
    ...            (slice(10, 15), [5, 2]))
    (slice(110, 115, None), [150, 120])

    And a variety of other interesting cases

    >>> fuse_slice(slice(1000, 2000), 10)  # integers
    1010

    >>> fuse_slice(slice(1000, 2000, 5), slice(10, 20, 2))
    slice(1050, 1100, 10)

    >>> fuse_slice(slice(1000, 2000, 5), [1, 2, 3])  # lists
    [1005, 1010, 1015]

    >>> fuse_slice(None, slice(None, None))  # doctest: +SKIP
    None
    """
    # None only works if the second side is a full slice
    if a is None and isinstance(b, slice) and b == slice(None, None):
        return None

    # Replace None with 0 and one in start and step
    if isinstance(a, slice):
        a = normalize_slice(a)
    if isinstance(b, slice):
        b = normalize_slice(b)

    if isinstance(a, slice) and isinstance(b, Integral):
        if b < 0:
            raise NotImplementedError()
        return a.start + b * a.step

    if isinstance(a, slice) and isinstance(b, slice):
        start = a.start + a.step * b.start
        if b.stop is not None:
            stop = a.start + a.step * b.stop
        else:
            stop = None
        if a.stop is not None:
            if stop is not None:
                stop = min(a.stop, stop)
            else:
                stop = a.stop
        step = a.step * b.step
        if step == 1:
            step = None
        return slice(start, stop, step)

    if isinstance(b, list):
        return [fuse_slice(a, bb) for bb in b]
    if isinstance(a, list) and isinstance(b, (Integral, slice)):
        return a[b]

    if isinstance(a, tuple) and not isinstance(b, tuple):
        b = (b,)

    # If given two tuples walk through both, being mindful of uneven sizes
    # and newaxes
    if isinstance(a, tuple) and isinstance(b, tuple):
        # Check for non-fusible cases with fancy-indexing
        a_has_lists = any(isinstance(item, list) for item in a)
        b_has_lists = any(isinstance(item, list) for item in b)
        if a_has_lists and b_has_lists:
            raise NotImplementedError("Can't handle multiple list indexing")
        elif a_has_lists:
            check_for_nonfusible_fancy_indexing(a, b)
        elif b_has_lists:
            check_for_nonfusible_fancy_indexing(b, a)

        j = 0
        result = list()
        for i in range(len(a)):
            #  axis ceased to exist  or we're out of b
            if isinstance(a[i], Integral) or j == len(b):
                result.append(a[i])
                continue
            while b[j] is None:  # insert any Nones on the rhs
                result.append(None)
                j += 1
            result.append(fuse_slice(a[i], b[j]))  # Common case
            j += 1
        while j < len(b):  # anything leftover on the right?
            result.append(b[j])
            j += 1
        return tuple(result)
    raise NotImplementedError()
