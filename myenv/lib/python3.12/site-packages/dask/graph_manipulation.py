"""Tools to modify already existing dask graphs. Unlike in :mod:`dask.optimization`, the
output collections produced by this module are typically not functionally equivalent to
their inputs.
"""
from __future__ import annotations

import uuid
from collections.abc import Callable, Hashable
from typing import Literal, TypeVar

from dask.base import (
    clone_key,
    get_collection_names,
    get_name_from_key,
    replace_name_in_key,
    tokenize,
    unpack_collections,
)
from dask.blockwise import blockwise
from dask.core import flatten
from dask.delayed import Delayed, delayed
from dask.highlevelgraph import HighLevelGraph, Layer, MaterializedLayer
from dask.typing import Graph, Key

__all__ = ("bind", "checkpoint", "clone", "wait_on")

T = TypeVar("T")


def checkpoint(
    *collections,
    split_every: float | Literal[False] | None = None,
) -> Delayed:
    """Build a :doc:`delayed` which waits until all chunks of the input collection(s)
    have been computed before returning None.

    Parameters
    ----------
    collections
        Zero or more Dask collections or nested data structures containing zero or more
        collections
    split_every: int >= 2 or False, optional
        Determines the depth of the recursive aggregation. If greater than the number of
        input keys, the aggregation will be performed in multiple steps; the depth of
        the aggregation graph will be :math:`log_{split_every}(input keys)`. Setting to
        a low value can reduce cache size and network transfers, at the cost of more CPU
        and a larger dask graph.

        Set to False to disable. Defaults to 8.

    Returns
    -------
    :doc:`delayed` yielding None
    """
    if split_every is None:
        split_every = 8
    elif split_every is not False:
        split_every = int(split_every)
        if split_every < 2:
            raise ValueError("split_every must be False, None, or >= 2")

    collections, _ = unpack_collections(*collections)
    if len(collections) == 1:
        return _checkpoint_one(collections[0], split_every)
    else:
        return delayed(chunks.checkpoint)(
            *(_checkpoint_one(c, split_every) for c in collections)
        )


def _checkpoint_one(collection, split_every) -> Delayed:
    tok = tokenize(collection)
    name = "checkpoint-" + tok

    keys_iter = flatten(collection.__dask_keys__())
    try:
        next(keys_iter)
        next(keys_iter)
    except StopIteration:
        # Collection has 0 or 1 keys; no need for a map step
        layer: Graph = {name: (chunks.checkpoint, collection.__dask_keys__())}
        dsk = HighLevelGraph.from_collections(name, layer, dependencies=(collection,))
        return Delayed(name, dsk)

    # Collection has 2+ keys; apply a two-step map->reduce algorithm so that we
    # transfer over the network and store in RAM only a handful of None's instead of
    # the full computed collection's contents
    dsks = []
    map_names = set()
    map_keys = []

    for prev_name in get_collection_names(collection):
        map_name = "checkpoint_map-" + tokenize(prev_name, tok)
        map_names.add(map_name)
        map_layer = _build_map_layer(chunks.checkpoint, prev_name, map_name, collection)
        map_keys += list(map_layer.get_output_keys())
        dsks.append(
            HighLevelGraph.from_collections(
                map_name, map_layer, dependencies=(collection,)
            )
        )

    # recursive aggregation
    reduce_layer: dict = {}
    while split_every and len(map_keys) > split_every:
        k = (name, len(reduce_layer))
        reduce_layer[k] = (chunks.checkpoint, map_keys[:split_every])
        map_keys = map_keys[split_every:] + [k]
    reduce_layer[name] = (chunks.checkpoint, map_keys)

    dsks.append(HighLevelGraph({name: reduce_layer}, dependencies={name: map_names}))
    dsk = HighLevelGraph.merge(*dsks)

    return Delayed(name, dsk)


def _can_apply_blockwise(collection) -> bool:
    """Return True if _map_blocks can be sped up via blockwise operations; False
    otherwise.

    FIXME this returns False for collections that wrap around around da.Array, such as
          pint.Quantity, xarray DataArray, Dataset, and Variable.
    """
    try:
        from dask.bag import Bag

        if isinstance(collection, Bag):
            return True
    except ImportError:
        pass
    try:
        from dask.array import Array

        if isinstance(collection, Array):
            return True
    except ImportError:
        pass
    try:
        from dask.dataframe import DataFrame, Series

        return isinstance(collection, (DataFrame, Series))
    except ImportError:
        return False


def _build_map_layer(
    func: Callable,
    prev_name: str,
    new_name: str,
    collection,
    dependencies: tuple[Delayed, ...] = (),
) -> Layer:
    """Apply func to all keys of collection. Create a Blockwise layer whenever possible;
    fall back to MaterializedLayer otherwise.

    Parameters
    ----------
    func
        Callable to be invoked on the graph node
    prev_name : str
        name of the layer to map from; in case of dask base collections, this is the
        collection name. Note how third-party collections, e.g. xarray.Dataset, can
        have multiple names.
    new_name : str
        name of the layer to map to
    collection
        Arbitrary dask collection
    dependencies
        Zero or more Delayed objects, which will be passed as arbitrary variadic args to
        func after the collection's chunk
    """
    if _can_apply_blockwise(collection):
        # Use a Blockwise layer
        try:
            numblocks = collection.numblocks
        except AttributeError:
            numblocks = (collection.npartitions,)
        indices = tuple(i for i, _ in enumerate(numblocks))
        kwargs = {"_deps": [d.key for d in dependencies]} if dependencies else {}

        return blockwise(
            func,
            new_name,
            indices,
            prev_name,
            indices,
            numblocks={prev_name: numblocks},
            dependencies=dependencies,
            **kwargs,
        )
    else:
        # Delayed, bag.Item, dataframe.core.Scalar, or third-party collection;
        # fall back to MaterializedLayer
        dep_keys = tuple(d.key for d in dependencies)
        return MaterializedLayer(
            {
                replace_name_in_key(k, {prev_name: new_name}): (func, k) + dep_keys
                for k in flatten(collection.__dask_keys__())
                if get_name_from_key(k) == prev_name
            }
        )


def bind(
    children: T,
    parents,
    *,
    omit=None,
    seed: Hashable | None = None,
    assume_layers: bool = True,
    split_every: float | Literal[False] | None = None,
) -> T:
    """
    Make ``children`` collection(s), optionally omitting sub-collections, dependent on
    ``parents`` collection(s). Two examples follow.

    The first example creates an array ``b2`` whose computation first computes an array
    ``a`` completely and then computes ``b`` completely, recomputing ``a`` in the
    process:

    >>> import dask
    >>> import dask.array as da
    >>> a = da.ones(4, chunks=2)
    >>> b = a + 1
    >>> b2 = bind(b, a)
    >>> len(b2.dask)
    9
    >>> b2.compute()
    array([2., 2., 2., 2.])

    The second example creates arrays ``b3`` and ``c3``, whose computation first
    computes an array ``a`` and then computes the additions, this time not
    recomputing ``a`` in the process:

    >>> c = a + 2
    >>> b3, c3 = bind((b, c), a, omit=a)
    >>> len(b3.dask), len(c3.dask)
    (7, 7)
    >>> dask.compute(b3, c3)
    (array([2., 2., 2., 2.]), array([3., 3., 3., 3.]))

    Parameters
    ----------
    children
        Dask collection or nested structure of Dask collections
    parents
        Dask collection or nested structure of Dask collections
    omit
        Dask collection or nested structure of Dask collections
    seed
        Hashable used to seed the key regeneration. Omit to default to a random number
        that will produce different keys at every call.
    assume_layers
        True
            Use a fast algorithm that works at layer level, which assumes that all
            collections in ``children`` and ``omit``

            #. use :class:`~dask.highlevelgraph.HighLevelGraph`,
            #. define the ``__dask_layers__()`` method, and
            #. never had their graphs squashed and rebuilt between the creation of the
               ``omit`` collections and the ``children`` collections; in other words if
               the keys of the ``omit`` collections can be found among the keys of the
               ``children`` collections, then the same must also hold true for the
               layers.
        False
            Use a slower algorithm that works at keys level, which makes none of the
            above assumptions.
    split_every
        See :func:`checkpoint`

    Returns
    -------
    Same as ``children``
        Dask collection or structure of dask collection equivalent to ``children``,
        which compute to the same values. All nodes of ``children`` will be regenerated,
        up to and excluding the nodes of ``omit``. Nodes immediately above ``omit``, or
        the leaf nodes if the collections in ``omit`` are not found, are prevented from
        computing until all collections in ``parents`` have been fully computed.
        The keys of the regenerated nodes will be different from the original ones, so
        that they can be used within the same graph.
    """
    if seed is None:
        seed = uuid.uuid4().bytes

    # parents=None is a special case invoked by the one-liner wrapper clone() below
    blocker = (
        checkpoint(parents, split_every=split_every) if parents is not None else None
    )

    omit, _ = unpack_collections(omit)
    if assume_layers:
        # Set of all the top-level layers of the collections in omit
        omit_layers = {layer for coll in omit for layer in coll.__dask_layers__()}
        omit_keys = set()
    else:
        omit_layers = set()
        # Set of *all* the keys, not just the top-level ones, of the collections in omit
        omit_keys = {key for coll in omit for key in coll.__dask_graph__()}

    unpacked_children, repack = unpack_collections(children)
    return repack(
        [
            _bind_one(child, blocker, omit_layers, omit_keys, seed)
            for child in unpacked_children
        ]
    )[0]


def _bind_one(
    child: T,
    blocker: Delayed | None,
    omit_layers: set[str],
    omit_keys: set[Key],
    seed: Hashable,
) -> T:
    prev_coll_names = get_collection_names(child)
    if not prev_coll_names:
        # Collection with no keys; this is a legitimate use case but, at the moment of
        # writing, can only happen with third-party collections
        return child

    dsk = child.__dask_graph__()  # type: ignore
    new_layers: dict[str, Layer] = {}
    new_deps: dict[str, set[str]] = {}

    if isinstance(dsk, HighLevelGraph):
        try:
            layers_to_clone = set(child.__dask_layers__())  # type: ignore
        except AttributeError:
            layers_to_clone = prev_coll_names.copy()
    else:
        if len(prev_coll_names) == 1:
            hlg_name = next(iter(prev_coll_names))
        else:
            hlg_name = tokenize(*prev_coll_names)
        dsk = HighLevelGraph.from_collections(hlg_name, dsk)
        layers_to_clone = {hlg_name}

    clone_keys = dsk.get_all_external_keys() - omit_keys
    for layer_name in omit_layers:
        try:
            layer = dsk.layers[layer_name]
        except KeyError:
            continue
        clone_keys -= layer.get_output_keys()
    # Note: when assume_layers=True, clone_keys can contain keys of the omit collections
    # that are not top-level. This is OK, as they will never be encountered inside the
    # values of their dependent layers.

    if blocker is not None:
        blocker_key = blocker.key
        blocker_dsk = blocker.__dask_graph__()
        assert isinstance(blocker_dsk, HighLevelGraph)
        new_layers.update(blocker_dsk.layers)
        new_deps.update(blocker_dsk.dependencies)
    else:
        blocker_key = None

    layers_to_copy_verbatim = set()

    while layers_to_clone:
        prev_layer_name = layers_to_clone.pop()
        new_layer_name = clone_key(prev_layer_name, seed=seed)
        if new_layer_name in new_layers:
            continue

        layer = dsk.layers[prev_layer_name]
        layer_deps = dsk.dependencies[prev_layer_name]
        layer_deps_to_clone = layer_deps - omit_layers
        layer_deps_to_omit = layer_deps & omit_layers
        layers_to_clone |= layer_deps_to_clone
        layers_to_copy_verbatim |= layer_deps_to_omit

        new_layers[new_layer_name], is_bound = layer.clone(
            keys=clone_keys, seed=seed, bind_to=blocker_key
        )
        new_dep = {
            clone_key(dep, seed=seed) for dep in layer_deps_to_clone
        } | layer_deps_to_omit
        if is_bound:
            new_dep.add(blocker_key)
        new_deps[new_layer_name] = new_dep

    # Add the layers of the collections from omit from child.dsk. Note that, when
    # assume_layers=False, it would be unsafe to simply do HighLevelGraph.merge(dsk,
    # omit[i].dsk). Also, collections in omit may or may not be parents of this specific
    # child, or of any children at all.
    while layers_to_copy_verbatim:
        layer_name = layers_to_copy_verbatim.pop()
        if layer_name in new_layers:
            continue
        layer_deps = dsk.dependencies[layer_name]
        layers_to_copy_verbatim |= layer_deps
        new_deps[layer_name] = layer_deps
        new_layers[layer_name] = dsk.layers[layer_name]

    rebuild, args = child.__dask_postpersist__()  # type: ignore
    return rebuild(
        HighLevelGraph(new_layers, new_deps),
        *args,
        rename={prev_name: clone_key(prev_name, seed) for prev_name in prev_coll_names},
    )


def clone(*collections, omit=None, seed: Hashable = None, assume_layers: bool = True):
    """Clone dask collections, returning equivalent collections that are generated from
    independent calculations.

    Examples
    --------
    (tokens have been simplified for the sake of brevity)

    >>> import dask.array as da
    >>> x_i = da.asarray([1, 1, 1, 1], chunks=2)
    >>> y_i = x_i + 1
    >>> z_i = y_i + 2
    >>> dict(z_i.dask)  # doctest: +SKIP
    {('array-1', 0): array([1, 1]),
     ('array-1', 1): array([1, 1]),
     ('add-2', 0): (<function operator.add>, ('array-1', 0), 1),
     ('add-2', 1): (<function operator.add>, ('array-1', 1), 1),
     ('add-3', 0): (<function operator.add>, ('add-2', 0), 1),
     ('add-3', 1): (<function operator.add>, ('add-2', 1), 1)}
    >>> w_i = clone(z_i, omit=x_i)
    >>> w_i.compute()
    array([4, 4, 4, 4])
    >>> dict(w_i.dask)  # doctest: +SKIP
    {('array-1', 0): array([1, 1]),
     ('array-1', 1): array([1, 1]),
     ('add-4', 0): (<function operator.add>, ('array-1', 0), 1),
     ('add-4', 1): (<function operator.add>, ('array-1', 1), 1),
     ('add-5', 0): (<function operator.add>, ('add-4', 0), 1),
     ('add-5', 1): (<function operator.add>, ('add-4', 1), 1)}

    The typical usage pattern for clone() is the following:

    >>> x = cheap_computation_with_large_output()  # doctest: +SKIP
    >>> y = expensive_and_long_computation(x)  # doctest: +SKIP
    >>> z = wrap_up(clone(x), y)  # doctest: +SKIP

    In the above code, the chunks of x will be forgotten as soon as they are consumed by
    the chunks of y, and then they'll be regenerated from scratch at the very end of the
    computation. Without clone(), x would only be computed once and then kept in memory
    throughout the whole computation of y, needlessly consuming memory.

    Parameters
    ----------
    collections
        Zero or more Dask collections or nested structures of Dask collections
    omit
        Dask collection or nested structure of Dask collections which will not be cloned
    seed
        See :func:`bind`
    assume_layers
        See :func:`bind`

    Returns
    -------
    Same as ``collections``
        Dask collections of the same type as the inputs, which compute to the same
        value, or nested structures equivalent to the inputs, where the original
        collections have been replaced.
        The keys of the regenerated nodes in the new collections will be different from
        the original ones, so that they can be used within the same graph.
    """
    out = bind(
        collections, parents=None, omit=omit, seed=seed, assume_layers=assume_layers
    )
    return out[0] if len(collections) == 1 else out


def wait_on(
    *collections,
    split_every: float | Literal[False] | None = None,
):
    """Ensure that all chunks of all input collections have been computed before
    computing the dependents of any of the chunks.

    The following example creates a dask array ``u`` that, when used in a computation,
    will only proceed when all chunks of the array ``x`` have been computed, but
    otherwise matches ``x``:

    >>> import dask.array as da
    >>> x = da.ones(10, chunks=5)
    >>> u = wait_on(x)

    The following example will create two arrays ``u`` and ``v`` that, when used in a
    computation, will only proceed when all chunks of the arrays ``x`` and ``y`` have
    been computed but otherwise match ``x`` and ``y``:

    >>> x = da.ones(10, chunks=5)
    >>> y = da.zeros(10, chunks=5)
    >>> u, v = wait_on(x, y)

    Parameters
    ----------
    collections
        Zero or more Dask collections or nested structures of Dask collections
    split_every
        See :func:`checkpoint`

    Returns
    -------
    Same as ``collections``
        Dask collection of the same type as the input, which computes to the same value,
        or a nested structure equivalent to the input where the original collections
        have been replaced.
        The keys of the regenerated nodes of the new collections will be different from
        the original ones, so that they can be used within the same graph.
    """
    blocker = checkpoint(*collections, split_every=split_every)

    def block_one(coll):
        tok = tokenize(coll, blocker)
        dsks = []
        rename = {}
        for prev_name in get_collection_names(coll):
            new_name = "wait_on-" + tokenize(prev_name, tok)
            rename[prev_name] = new_name
            layer = _build_map_layer(
                chunks.bind, prev_name, new_name, coll, dependencies=(blocker,)
            )
            dsks.append(
                HighLevelGraph.from_collections(
                    new_name, layer, dependencies=(coll, blocker)
                )
            )
        dsk = HighLevelGraph.merge(*dsks)
        rebuild, args = coll.__dask_postpersist__()
        return rebuild(dsk, *args, rename=rename)

    unpacked, repack = unpack_collections(*collections)
    out = repack([block_one(coll) for coll in unpacked])
    return out[0] if len(collections) == 1 else out


class chunks:
    """Callables to be inserted in the Dask graph"""

    @staticmethod
    def bind(node: T, *args, **kwargs) -> T:
        """Dummy graph node of :func:`bind` and :func:`wait_on`.
        Wait for both node and all variadic args to complete; then return node.
        """
        return node

    @staticmethod
    def checkpoint(*args, **kwargs) -> None:
        """Dummy graph node of :func:`checkpoint`.
        Wait for all variadic args to complete; then return None.
        """
        pass
