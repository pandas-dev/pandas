from __future__ import annotations

import itertools
import os
from collections.abc import Hashable, Iterable, Mapping, Sequence
from dataclasses import fields, is_dataclass, replace
from itertools import product
from math import prod
from typing import Any

import tlz as toolz

import dask
from dask import base, utils
from dask._task_spec import Dict, GraphNode, List, Task, TaskRef, parse_input
from dask.base import clone_key, get_name_from_key, tokenize
from dask.core import flatten, ishashable, keys_in_tasks, reverse_dict
from dask.delayed import Delayed, finalize
from dask.highlevelgraph import HighLevelGraph, Layer
from dask.optimization import fuse
from dask.typing import Key
from dask.utils import _deprecated, ensure_dict, homogeneous_deepmap


class BlockwiseDep:
    """Blockwise-IO argument

    This is the base class for indexable Blockwise-IO arguments.
    When constructing a ``Blockwise`` Layer, one or more of the
    collection tuples passed in with ``indices`` may contain a
    ``BlockwiseDep`` instance (in place of a "real" collection name).
    This allows a new collection to be created (via IO) within a
    ``Blockwise`` layer.

    Parameters
    ----------
    numblocks: tuple[int, ...]
        The number of blocks/partitions the object can support
        along each dimension.
    produces_tasks: bool
        Whether any nested tasks will be passed to the Blockwise
        function.

    See Also
    --------
    dask.blockwise.Blockwise
    dask.blockwise.BlockwiseDepDict
    """

    numblocks: tuple[int, ...]
    produces_tasks: bool

    def __getitem__(self, idx: tuple[int, ...]) -> Any:
        """Return Blockwise-function arguments for a specific index"""
        raise NotImplementedError(
            "Must define `__getitem__` for `BlockwiseDep` subclass."
        )

    def get(self, idx: tuple[int, ...], default) -> Any:
        """BlockwiseDep ``__getitem__`` Wrapper"""
        try:
            return self.__getitem__(idx)
        except KeyError:
            return default

    @property
    def produces_keys(self) -> bool:
        """Whether this object will produce external key dependencies.

        An external key corresponds to a task key or ``Delayed``-object
        key that does not originate from within the ``Blockwise`` layer
        that is including this ``BlockwiseDep`` object in its ``indices``.
        A ``BlockwiseDep`` object should only return external-key
        dependencies when those dependencies do not correspond to a
        blockwise-compatible Dask collection (otherwise the collection
        name should just be included in ``indices`` list instead).
        """
        return False

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.numblocks}>"


class BlockwiseDepDict(BlockwiseDep):
    """Dictionary-based Blockwise-IO argument

    This is a dictionary-backed instance of ``BlockwiseDep``.
    The purpose of this class is to simplify the construction
    of IO-based Blockwise Layers with block/partition-dependent
    function arguments that are difficult to calculate at
    graph-materialization time.

    Examples
    --------

    Specify an IO-based function for the Blockwise Layer. Note
    that the function will be passed a single input object when
    the task is executed (e.g. a single ``tuple`` or ``dict``):

    >>> import pandas as pd
    >>> func = lambda x: pd.read_csv(**x)

    Use ``BlockwiseDepDict`` to define the input argument to
    ``func`` for each block/partition:

    >>> dep = BlockwiseDepDict(
    ...     mapping={
    ...         (0,) : {
    ...             "filepath_or_buffer": "data.csv",
    ...             "skiprows": 1,
    ...             "nrows": 2,
    ...             "names": ["a", "b"],
    ...         },
    ...         (1,) : {
    ...             "filepath_or_buffer": "data.csv",
    ...             "skiprows": 3,
    ...             "nrows": 2,
    ...             "names": ["a", "b"],
    ...         },
    ...     }
    ... )

    Construct a Blockwise Layer with ``dep`` specified
    in the ``indices`` list:

    >>> layer = Blockwise(
    ...     output="collection-name",
    ...     output_indices="i",
    ...     task=Task("collection-name", func, TaskRef("_0")),
    ...     indices=[(dep, "i")],
    ...     numblocks={},
    ... )

    See Also
    --------
    dask.blockwise.Blockwise
    dask.blockwise.BlockwiseDep
    """

    def __init__(
        self,
        mapping: dict,
        numblocks: tuple[int, ...] | None = None,
        produces_tasks: bool = False,
        produces_keys: bool = False,
    ):
        self.mapping = mapping
        self.produces_tasks = produces_tasks

        # By default, assume 1D shape
        self.numblocks = numblocks or (len(mapping),)

        # Whether `mapping` values are real task keys
        # (e.g. Delayed objects)
        self._produces_keys = produces_keys

    @property
    def produces_keys(self) -> bool:
        return self._produces_keys

    def __getitem__(self, idx: tuple[int, ...]) -> Any:
        try:
            return self.mapping[idx]
        except KeyError as err:
            # If a DataFrame collection was converted
            # to an Array collection, the dimension of
            # `idx` may not agree with the keys in
            # `self.mapping`. In this case, we can
            # use `self.numblocks` to check for a key
            # match in the leading elements of `idx`
            flat_idx = idx[: len(self.numblocks)]
            if flat_idx in self.mapping:
                return self.mapping[flat_idx]
            raise err

    def __len__(self) -> int:
        return len(self.mapping)


class BlockIndex(BlockwiseDep):
    """Index BlockwiseDep argument

    The purpose of this class is to provide each
    block of a ``Blockwise``-based operation with
    the current block index.
    """

    produces_tasks: bool = False

    def __init__(self, numblocks: tuple[int, ...]):
        # NOTE: Unused - Just needs to be set to
        # follow the `BlockwiseDep` interface
        self.numblocks = numblocks

    def __getitem__(self, idx: tuple[int, ...]) -> tuple[int, ...]:
        return idx


def subs(task, substitution):
    """Create a new task with the values substituted

    This is like dask.core.subs, but takes a dict of many substitutions to
    perform simultaneously.  It is not as concerned with micro performance.
    """
    if isinstance(task, dict):
        return {k: subs(v, substitution) for k, v in task.items()}
    if type(task) in (tuple, list, set):
        return type(task)([subs(x, substitution) for x in task])
    try:
        return substitution[task]
    except (KeyError, TypeError):
        return task


def index_subs(ind, substitution):
    """A simple subs function that works both on tuples and strings"""
    if ind is None:
        return ind
    else:
        return tuple(substitution.get(c, c) for c in ind)


_BLOCKWISE_DEFAULT_PREFIX = "__dask_blockwise__"


def blockwise_token(i, prefix=_BLOCKWISE_DEFAULT_PREFIX):
    return prefix + "%d" % i


def blockwise(
    func,
    output,
    output_indices,
    *arrind_pairs,
    numblocks=None,
    concatenate=None,
    new_axes=None,
    dependencies=(),
    _data_producer=False,
    **kwargs,
):
    """Create a Blockwise symbolic mutable mapping

    Applies a function, ``func``, across blocks from many different input
    collections.  We arrange the pattern with which those blocks interact with
    sets of matching indices.  E.g.::

        blockwise(func, 'z', 'i', 'x', 'i', 'y', 'i')

    yield an embarrassingly parallel communication pattern and is read as

        $$ z_i = func(x_i, y_i) $$

    More complex patterns may emerge, including multiple indices::

        blockwise(func, 'z', 'ij', 'x', 'ij', 'y', 'ji')

        $$ z_{ij} = func(x_{ij}, y_{ji}) $$

    Indices missing in the output but present in the inputs results in many
    inputs being sent to one function (see examples).

    Examples
    --------
    Simple embarrassing map operation

    >>> def inc(x):
    ...    return x + 1
    >>> dict(
    ...     blockwise(
    ...         inc,
    ...         'z', 'ij',
    ...         'x', 'ij',
    ...         numblocks={'x': (2, 2)}
    ...     )
    ... )  # doctest: +NORMALIZE_WHITESPACE
    {('z', 0, 0): <Task ('z', 0, 0) inc(Alias(('x', 0, 0)))>,
    ('z', 0, 1): <Task ('z', 0, 1) inc(Alias(('x', 0, 1)))>,
    ('z', 1, 0): <Task ('z', 1, 0) inc(Alias(('x', 1, 0)))>,
    ('z', 1, 1): <Task ('z', 1, 1) inc(Alias(('x', 1, 1)))>}

    Simple operation on two datasets x and y

    >>> def add(x, y):
    ...     return x + y
    >>> dict(
    ...     blockwise(
    ...         add,
    ...         'z', 'ij',
    ...         'x', 'ij',
    ...         'y', 'ij',
    ...         numblocks={'x': (2, 2), 'y': (2, 2)}
    ...     )
    ... )  # doctest: +NORMALIZE_WHITESPACE
    {('z', 0, 0): <Task ('z', 0, 0) add(Alias(('x', 0, 0)), Alias(('y', 0, 0)))>,
    ('z', 0, 1): <Task ('z', 0, 1) add(Alias(('x', 0, 1)), Alias(('y', 0, 1)))>,
    ('z', 1, 0): <Task ('z', 1, 0) add(Alias(('x', 1, 0)), Alias(('y', 1, 0)))>,
    ('z', 1, 1): <Task ('z', 1, 1) add(Alias(('x', 1, 1)), Alias(('y', 1, 1)))>}

    Operation that flips one of the datasets

    >>> def addT(x, y):
    ...     return x + y.T
    >>> dict(
    ...     blockwise(
    ...         addT,
    ...         'z', 'ij',
    ...         'x', 'ij',
    ...         'y', 'ji',
    ...         numblocks={'x': (2, 2), 'y': (2, 2)}
    ...     )
    ... )  # doctest: +NORMALIZE_WHITESPACE
    {('z', 0, 0): <Task ('z', 0, 0) addT(Alias(('x', 0, 0)), Alias(('y', 0, 0)))>,
    ('z', 0, 1): <Task ('z', 0, 1) addT(Alias(('x', 0, 1)), Alias(('y', 1, 0)))>,
    ('z', 1, 0): <Task ('z', 1, 0) addT(Alias(('x', 1, 0)), Alias(('y', 0, 1)))>,
    ('z', 1, 1): <Task ('z', 1, 1) addT(Alias(('x', 1, 1)), Alias(('y', 1, 1)))>}

    Dot product with contraction over ``j`` index.  Yields list arguments

    >>> from dask.array.core import dotmany
    >>> dict(
    ...     blockwise(
    ...         dotmany,
    ...         'z', 'ik',
    ...         'x', 'ij',
    ...         'y', 'jk',
    ...         numblocks={'x': (2, 2), 'y': (2, 2)}
    ...     )
    ... )  # doctest: +SKIP
    {('z', 0,  0):
        <Task ('z', 0, 0) dotmany(
        List((Alias(('x', 0, 0)), Alias(('x', 0, 1)))),
        List((Alias(('y', 0, 0)), Alias(('y', 1, 0)))))
        >,
    ('z', 0,  1):
        <Task ('z', 0, 1) dotmany(
        List((Alias(('x', 0, 0)), Alias(('x', 0, 1)))),
        List((Alias(('y', 0, 1)), Alias(('y', 1, 1)))))
        >,
    ('z', 1,  0):
        <Task ('z', 1, 0) dotmany(
        List((Alias(('x', 1, 0)), Alias(('x', 1, 1)))),
        List((Alias(('y', 0, 0)), Alias(('y', 1, 0)))))
        >,
    ('z', 1,  1):
        <Task ('z', 1, 1) dotmany(
        List((Alias(('x', 1, 0)), Alias(('x', 1, 1)))),
        List((Alias(('y', 0, 1)), Alias(('y', 1, 1)))))
        >
    }

    Pass ``concatenate=True`` to concatenate arrays ahead of time.
    ``concatenate_axes`` is then called on the lists before the arguments are
    passed to the function.

    Supports Broadcasting rules

    >>> dict(
    ...     blockwise(
    ...         add,
    ...         'z', 'ij',
    ...         'x', 'ij',
    ...         'y', 'ij',
    ...         concatenate=True,
    ...         numblocks={'x': (1, 2), 'y': (2, 2)}
    ...     )
    ... )  # doctest: +SKIP
    {
    ('z', 0, 0): <Task ('z', 0, 0) add(Alias(('x', 0, 0)), Alias(('y', 0, 0)))>,
    ('z', 0, 1): <Task ('z', 0, 1) add(Alias(('x', 0, 1)), Alias(('y', 0, 1)))>,
    ('z', 1, 0): <Task ('z', 1, 0) add(Alias(('x', 0, 0)), Alias(('y', 1, 0)))>,
    ('z', 1, 1): <Task ('z', 1, 1) add(Alias(('x', 0, 1)), Alias(('y', 1, 1)))>
    }

    Include literals by indexing with ``None``

    >>> dict(
    ...     blockwise(
    ...         add,
    ...         'z', 'i',
    ...         'x', 'i',
    ...         100, None,
    ...         numblocks={'x': (2, )}
    ...     )
    ... )  # doctest: +NORMALIZE_WHITESPACE
    {('z', 0): <Task ('z', 0) add(Alias(('x', 0)), DataNode(100))>,
    ('z', 1): <Task ('z', 1) add(Alias(('x', 1)), DataNode(100))>}

    If the broadcasted value is a delayed object or other dask collection, the
    key has to be wrapped appropriately as an Alias object.

    >>> from dask import delayed
    >>> from dask._task_spec import Alias
    >>> delayed_inc = delayed(inc)(5, dask_key_name='dinc')
    >>> dict(
    ...     blockwise(
    ...         add,
    ...         'z', 'i',
    ...         'x', 'i',
    ...         Alias(delayed_inc.key), None,
    ...         numblocks={'x': (2, )}
    ...     )
    ... )  # doctest: +NORMALIZE_WHITESPACE
    {('z', 0): <Task ('z', 0) add(Alias(('x', 0)), Alias('dinc'))>,
    ('z', 1): <Task ('z', 1) add(Alias(('x', 1)), Alias('dinc'))>}

    See Also
    --------
    dask.array.blockwise
    Blockwise
    """
    new_axes = new_axes or {}

    arrind_pairs = list(arrind_pairs)

    # Transform indices to canonical elements
    # We use terms like _0, and _1 rather than provided index elements
    unique_indices = {
        i for ii in arrind_pairs[1::2] if ii is not None for i in ii
    } | set(output_indices)
    sub = {k: blockwise_token(i, ".") for i, k in enumerate(sorted(unique_indices))}
    output_indices = index_subs(tuple(output_indices), sub)
    a_pairs_list = []
    for a in arrind_pairs[1::2]:
        if a is not None:
            val = tuple(a)
        else:
            val = a
        a_pairs_list.append(index_subs(val, sub))

    arrind_pairs[1::2] = a_pairs_list
    new_axes = {index_subs((k,), sub)[0]: v for k, v in new_axes.items()}

    # Unpack dask values in non-array arguments
    inputs = []
    inputs_indices = []
    for name, index in toolz.partition(2, arrind_pairs):
        inputs.append(name)
        inputs_indices.append(index)

    task_args = list(map(TaskRef, map(blockwise_token, range(len(inputs)))))
    # Unpack delayed objects in kwargs
    new_keys = {n for c in dependencies for n in c.__dask_layers__()}
    if kwargs:
        # replace keys in kwargs with _0 tokens
        new_tokens = tuple(
            blockwise_token(i) for i in range(len(inputs), len(inputs) + len(new_keys))
        )
        sub = dict(zip(new_keys, new_tokens))
        inputs.extend(map(TaskRef, new_keys))
        inputs_indices.extend((None,) * len(new_keys))

    indices = [(k, v) for k, v in zip(inputs, inputs_indices)]
    task = Task(output, func, *task_args, _data_producer=_data_producer, **kwargs)
    subgraph = Blockwise(
        output,
        output_indices,
        task,
        indices,
        numblocks=numblocks,
        concatenate=concatenate,
        new_axes=new_axes,
    )
    return subgraph


class Blockwise(Layer):
    """Tensor Operation

    This is a lazily constructed mapping for tensor operation graphs.
    This defines a dictionary using an operation and an indexing pattern.
    It is built for many operations like elementwise, transpose, tensordot, and
    so on.  We choose to keep these as symbolic mappings rather than raw
    dictionaries because we are able to fuse them during optimization,
    sometimes resulting in much lower overhead.

    Parameters
    ----------
    output: str
        The name of the output collection.  Used in keynames
    output_indices: tuple
        The output indices, like ``('i', 'j', 'k')`` used to determine the
        structure of the block computations
    dsk: dict
        A small graph to apply per-output-block.  May include keys from the
        input indices.
    indices: tuple[tuple[str, tuple[str, ...] | None], ...]
        An ordered mapping from input key name, like ``'x'``
        to input indices, like ``('i', 'j')``
        Or includes literals, which have ``None`` for an index value.
        In place of input-key names, the first tuple element may also be a
        ``BlockwiseDep`` object.
    numblocks: Mapping[key, Sequence[int]]
        Number of blocks along each dimension for each input
    concatenate: bool
        Whether or not to pass contracted dimensions as a list of inputs or a
        single input to the block function
    new_axes: Mapping
        New index dimensions that may have been created and their size,
        e.g. ``{'j': 2, 'k': 3}``
    output_blocks: set[tuple[int, ...]]
        Specify a specific set of required output blocks. Since the graph
        will only contain the necessary tasks to generate these outputs,
        this kwarg can be used to "cull" the abstract layer (without needing
        to materialize the low-level graph).
    annotations: dict (optional)
        Layer annotations
    io_deps: dict[str, BlockwiseDep] (optional)
        Dictionary containing the mapping between "place-holder" collection
        keys and ``BlockwiseDep``-based objects.
        **WARNING**: This argument should only be used internally (for culling,
        fusion and cloning of existing Blockwise layers). Explicit use of this
        argument will be deprecated in the future.

    See Also
    --------
    dask.blockwise.blockwise
    dask.array.blockwise
    """

    output: str
    output_indices: tuple[str, ...]
    task: Task
    indices: tuple[tuple[str | TaskRef, tuple[str, ...] | None], ...]
    numblocks: Mapping[str, Sequence[int]]
    concatenate: bool | None
    new_axes: Mapping[str, int]
    output_blocks: set[tuple[int, ...]] | None
    io_deps: Mapping[str, BlockwiseDep]

    def __init__(
        self,
        output: str,
        output_indices: Iterable[str],
        task: Task,
        indices: Iterable[tuple[str | TaskRef | BlockwiseDep, Iterable[str] | None]],
        numblocks: Mapping[str, Sequence[int]],
        concatenate: bool | None = None,
        new_axes: Mapping[str, int] | None = None,
        output_blocks: set[tuple[int, ...]] | None = None,
        annotations: Mapping[str, Any] | None = None,
        io_deps: Mapping[str, BlockwiseDep] | None = None,
    ):
        super().__init__(annotations=annotations)
        self.output = output
        self.output_indices = tuple(output_indices)
        self.output_blocks = output_blocks
        self.task = task
        assert isinstance(task, Task)

        # Remove `BlockwiseDep` arguments from input indices
        # and add them to `self.io_deps`.
        # TODO: Remove `io_deps` and handle indexable objects
        # in `self.indices` throughout `Blockwise`.
        _tmp_indices = []
        numblocks = dict(numblocks)
        io_deps = dict(io_deps or {})
        if indices:
            for dep, ind in indices:
                if ind is not None:
                    # FIXME: The Blockwise API is a little weird this way
                    assert not isinstance(
                        dep, TaskRef
                    ), "TaskRef objects are only allowed for broadcasted inputs with None as index."
                if isinstance(dep, BlockwiseDep):
                    name = tokenize(dep)
                    io_deps[name] = dep
                    numblocks[name] = dep.numblocks
                else:
                    name = dep  # type: ignore
                _tmp_indices.append((name, tuple(ind) if ind is not None else ind))
        self.numblocks = numblocks
        self.io_deps = io_deps or {}
        self.indices = tuple(_tmp_indices)

        # optimize_blockwise won't merge where `concatenate` doesn't match, so
        # enforce a canonical value if there are no axes for reduction.
        output_indices_set = set(self.output_indices)
        if concatenate is not None and all(
            i in output_indices_set
            for name, ind in self.indices
            if ind is not None
            for i in ind
        ):
            concatenate = None
        self.concatenate = concatenate
        self.new_axes = new_axes or {}

    @property
    def dims(self):
        """Returns a dictionary mapping between each index specified in
        `self.indices` and the number of output blocks for that indice.
        """
        if not hasattr(self, "_dims"):
            self._dims = _make_dims(self.indices, self.numblocks, self.new_axes)
        return self._dims

    def __repr__(self):
        return f"Blockwise<{self.indices} -> {self.output}>"

    @property
    def _dict(self):
        if hasattr(self, "_cached_dict"):
            return self._cached_dict["dsk"]
        else:
            keys = tuple(map(blockwise_token, range(len(self.indices))))

            dsk = _make_blockwise_graph(
                self.task,
                self.output,
                self.output_indices,
                *list(toolz.concat(self.indices)),
                new_axes=self.new_axes,
                numblocks=self.numblocks,
                concatenate=self.concatenate,
                output_blocks=self.output_blocks,
                dims=self.dims,
                io_deps=self.io_deps,
                keys=keys,
            )

            self._cached_dict = {"dsk": dsk}
        return self._cached_dict["dsk"]

    def get_output_keys(self):
        if self.output_blocks:
            # Culling has already generated a list of output blocks
            return {(self.output, *p) for p in self.output_blocks}

        # Return all possible output keys (no culling)
        return {
            (self.output, *p)
            for p in itertools.product(
                *[range(self.dims[i]) for i in self.output_indices]
            )
        }

    def __getitem__(self, key):
        return self._dict[key]

    def __iter__(self):
        return iter(self._dict)

    def __len__(self) -> int:
        # same method as `get_output_keys`, without manifesting the keys themselves
        return (
            len(self.output_blocks)
            if self.output_blocks
            else prod(self.dims[i] for i in self.output_indices)
        )

    def is_materialized(self):
        return hasattr(self, "_cached_dict")

    def _cull_dependencies(self, all_hlg_keys, output_blocks):
        """Determine the necessary dependencies to produce `output_blocks`.

        This method does not require graph materialization.
        """

        # Check `concatenate` option
        concatenate = None
        if self.concatenate is True:
            from dask.array.core import concatenate_axes as concatenate

        # Generate coordinate map
        (coord_maps, concat_axes, dummies) = _get_coord_mapping(
            self.dims,
            self.output_indices,
            self.numblocks,
            self.indices,
            concatenate,
        )

        # Gather constant dependencies (for all output keys)
        const_deps = set()
        for arg, ind in self.indices:
            if isinstance(arg, TaskRef):
                arg = arg.key
            if ind is None:
                try:
                    if arg in all_hlg_keys:
                        const_deps.add(arg)
                except TypeError:
                    pass  # unhashable

        # Get dependencies for each output block
        key_deps = {}
        for out_coords in output_blocks:
            deps = set()
            coords = out_coords + dummies
            for cmap, axes, (arg, ind) in zip(coord_maps, concat_axes, self.indices):
                if ind is not None and arg not in self.io_deps:
                    arg_coords = tuple(coords[c] for c in cmap)
                    if axes:
                        tups = _lol_product((arg,), arg_coords)
                        deps.update(flatten(tups))
                        if concatenate:
                            tups = (concatenate, tups, axes)
                    else:
                        tups = (arg,) + arg_coords
                        deps.add(tups)
            key_deps[(self.output,) + out_coords] = deps | const_deps

        # Add valid-key dependencies from io_deps
        for key, io_dep in self.io_deps.items():
            if io_dep.produces_keys:
                for out_coords in output_blocks:
                    key = (self.output,) + out_coords
                    valid_key_dep = io_dep[out_coords]
                    if isinstance(valid_key_dep, TaskRef):
                        valid_key_dep = valid_key_dep.key
                    key_deps[key] |= {valid_key_dep}

        return key_deps

    def _cull(self, output_blocks):
        return Blockwise(
            self.output,
            self.output_indices,
            self.task,
            self.indices,
            self.numblocks,
            concatenate=self.concatenate,
            new_axes=self.new_axes,
            output_blocks=output_blocks,
            annotations=self.annotations,
            io_deps=self.io_deps,
        )

    def cull(
        self, keys: set, all_hlg_keys: Iterable
    ) -> tuple[Layer, Mapping[Key, set[Key]]]:
        # Culling is simple for Blockwise layers.  We can just
        # collect a set of required output blocks (tuples), and
        # only construct graph for these blocks in `make_blockwise_graph`

        output_blocks: set[tuple[int, ...]] = set()
        for key in keys:
            if key[0] == self.output:
                output_blocks.add(key[1:])
        culled_deps = self._cull_dependencies(all_hlg_keys, output_blocks)
        out_size_iter = (self.dims[i] for i in self.output_indices)
        if prod(out_size_iter) != len(culled_deps):
            culled_layer = self._cull(output_blocks)
            return culled_layer, culled_deps
        else:
            return self, culled_deps

    def clone(
        self,
        keys: set[Key],
        seed: Hashable,
        bind_to: Key | None = None,
    ) -> tuple[Layer, bool]:
        names = {get_name_from_key(k) for k in keys}
        # We assume that 'keys' will contain either all or none of the output keys of
        # each of the layers, because clone/bind are always invoked at collection level.
        # Asserting this is very expensive, so we only check it during unit tests.
        if "PYTEST_CURRENT_TEST" in os.environ:
            assert not self.get_output_keys() - keys
            for name, nb in self.numblocks.items():
                if name in names:
                    for block in product(*(list(range(nbi)) for nbi in nb)):
                        assert (name, *block) in keys

        is_leaf = True

        indices = []
        k: Key | TaskRef
        for k, idxv in self.indices:
            # Note: k may not be a key and thus not be hashable in the case where
            # one or more args of blockwise() are sequences of literals;
            # e.g. k = (list, [0, 1, 2])
            # See https://github.com/dask/dask/issues/8978

            if ishashable(k) and k in names:
                is_leaf = False
                k = clone_key(k, seed)  # type: ignore
            elif isinstance(k, TaskRef) and k.key in names:
                is_leaf = False
                k = TaskRef(clone_key(k.key, seed))

            indices.append((k, idxv))

        numblocks: dict[str, Sequence[int]] = {}
        for k, nbv in self.numblocks.items():
            if k in names:
                is_leaf = False
                k = clone_key(k, seed)
            numblocks[k] = nbv

        if bind_to is not None and is_leaf:
            from dask.graph_manipulation import chunks

            # It's always a Delayed generated by dask.graph_manipulation.checkpoint;
            # the layer name always matches the key
            assert isinstance(bind_to, str)
            newtask = Task(
                clone_key(self.task.key, seed),
                chunks.bind,
                self.task,
                TaskRef(blockwise_token(len(indices))),
                _data_producer=self.task.data_producer,
            )
            indices.append((TaskRef(bind_to), None))
        else:
            newtask = self.task.substitute({}, key=clone_key(self.task.key, seed))

        return (
            Blockwise(
                output=clone_key(self.output, seed),
                output_indices=self.output_indices,
                task=newtask,
                indices=indices,
                numblocks=numblocks,
                concatenate=self.concatenate,
                new_axes=self.new_axes,
                output_blocks=self.output_blocks,
                annotations=self.annotations,
                io_deps=self.io_deps,
            ),
            (bind_to is not None and is_leaf),
        )


def _get_coord_mapping(
    dims,
    out_indices,
    numblocks,
    argpairs,
    concatenate,
):
    """Calculate coordinate mapping for graph construction.

    This function handles the high-level logic behind Blockwise graph
    construction. The output is a tuple containing: The mapping between
    input and output block coordinates (`coord_maps`), the axes along
    which to concatenate for each input (`concat_axes`), and the dummy
    indices needed for broadcasting (`dummies`).

    Used by `make_blockwise_graph` and `Blockwise._cull_dependencies`.

    Parameters
    ----------
    dims : dict
        Mapping between each index specified in `argpairs` and
        the number of output blocks for that index. Corresponds
        to the Blockwise `dims` attribute.
    out_indices : tuple
        Corresponds to the Blockwise `output_indices` attribute.
    numblocks : dict
        Corresponds to the Blockwise `numblocks` attribute.
    argpairs : tuple
        Corresponds to the Blockwise `indices` attribute.
    concatenate : bool
        Corresponds to the Blockwise `concatenate` attribute.
    """

    block_names = set()
    all_indices = set()
    for name, ind in argpairs:
        if ind is not None:
            block_names.add(name)
            for x in ind:
                all_indices.add(x)
    assert set(numblocks) == block_names, (numblocks, block_names)

    dummy_indices = all_indices - set(out_indices)

    # For each position in the output space, we'll construct a
    # "coordinate set" that consists of
    # - the output indices
    # - the dummy indices
    # - the dummy indices, with indices replaced by zeros (for broadcasting), we
    #   are careful to only emit a single dummy zero when concatenate=True to not
    #   concatenate the same array with itself several times.
    # - a 0 to assist with broadcasting.

    index_pos, zero_pos = {}, {}
    for i, ind in enumerate(out_indices):
        index_pos[ind] = i
        zero_pos[ind] = -1

    _dummies_list = []
    for i, ind in enumerate(dummy_indices):
        index_pos[ind] = 2 * i + len(out_indices)
        zero_pos[ind] = 2 * i + 1 + len(out_indices)
        reps = 1 if concatenate else dims[ind]
        _dummies_list.append([list(range(dims[ind])), [0] * reps])

    # ([0, 1, 2], [0, 0, 0], ...)  For a dummy index of dimension 3
    dummies = tuple(itertools.chain.from_iterable(_dummies_list))
    dummies += (0,)

    # For each coordinate position in each input, gives the position in
    # the coordinate set.
    coord_maps = []

    # Axes along which to concatenate, for each input
    concat_axes = []
    for arg, ind in argpairs:
        if ind is not None:
            coord_maps.append(
                [
                    zero_pos[i] if nb == 1 else index_pos[i]
                    for i, nb in zip(ind, numblocks[arg])
                ]
            )
            concat_axes.append([n for n, i in enumerate(ind) if i in dummy_indices])
        else:
            coord_maps.append(None)
            concat_axes.append(None)

    return coord_maps, concat_axes, dummies


def _make_blockwise_graph(
    task,
    output,
    out_indices,
    *arrind_pairs,
    numblocks=None,
    concatenate=None,
    new_axes=None,
    output_blocks=None,
    dims=None,
    io_deps=None,
    keys=None,
):
    if numblocks is None:
        raise ValueError("Missing required numblocks argument.")
    new_axes = new_axes or {}
    io_deps = io_deps or {}
    argpairs = list(toolz.partition(2, arrind_pairs))

    if concatenate is True:
        from dask.array.core import concatenate_axes as concatenate

    # Dictionary mapping {i: 3, j: 4, ...} for i, j, ... the dimensions
    dims = dims or _make_dims(argpairs, numblocks, new_axes)

    # Generate the abstract "plan" before constructing
    # the actual graph
    (coord_maps, concat_axes, dummies) = _get_coord_mapping(
        dims,
        out_indices,
        numblocks,
        argpairs,
        concatenate,
    )

    # Apply Culling.
    # Only need to construct the specified set of output blocks.
    # Note that we must convert itertools.product to list,
    # because we may need to loop through output_blocks more than
    # once below (itertools.product already uses an internal list,
    # so this is not a memory regression)
    output_blocks = output_blocks or list(
        itertools.product(*[range(dims[i]) for i in out_indices])
    )
    from dask._task_spec import DataNode

    dsk = {}
    for out_coords in output_blocks:
        this_task = task
        coords = out_coords + dummies
        args = []
        for cmap, axes, (arg, ind), key in zip(
            coord_maps, concat_axes, argpairs, keys, strict=True
        ):
            if key not in task.dependencies:
                # FIXME: This feels like a bug
                continue
            if ind is None:
                if not isinstance(arg, (GraphNode, TaskRef)):
                    this_task = this_task.substitute({key: DataNode(None, arg)})
                else:
                    this_task = this_task.substitute({key: arg})
                continue
            arg_coords = tuple(coords[c] for c in cmap)
            if arg in io_deps:
                val = parse_input(io_deps[arg].get(arg_coords, arg_coords))
                if not isinstance(val, GraphNode):
                    val = DataNode(None, val)
                this_task = this_task.substitute({key: val})
            else:
                subs = {}
                if axes:
                    tups = _lol_product((arg,), arg_coords, as_taskref=True)
                    if concatenate:
                        tups = Task(key, concatenate, tups, axes)
                    subs[key] = tups
                else:
                    subs[key] = (arg, *arg_coords)
                this_task = this_task.substitute(subs)
        new_key = (output,) + out_coords
        assert isinstance(this_task, Task)
        dsk[new_key] = Task.fuse(this_task, *args, key=new_key)
    return dsk


def _lol_product(head, values, as_taskref=False):
    """List of list of tuple keys, similar to `itertools.product`.

    Parameters
    ----------
    head : tuple
        Prefix prepended to all results.
    values : sequence
        Mix of singletons and lists. Each list is substituted with every
        possible value and introduces another level of list in the output.

    Examples
    --------
    >>> _lol_product(('x',), (1, 2, 3))
    ('x', 1, 2, 3)
    >>> _lol_product(('x',), (1, [2, 3], 4, [5, 6]))  # doctest: +NORMALIZE_WHITESPACE
    [[('x', 1, 2, 4, 5), ('x', 1, 2, 4, 6)],
     [('x', 1, 3, 4, 5), ('x', 1, 3, 4, 6)]]
    """
    if not values:
        if as_taskref:
            return TaskRef(head)
        return head
    elif isinstance(values[0], list):
        # FIXME: Constructor of List is odd
        if as_taskref:
            return List(
                *(
                    _lol_product(head + (x,), values[1:], as_taskref=as_taskref)
                    for x in values[0]
                )
            )
        else:
            return list(
                _lol_product(head + (x,), values[1:], as_taskref=as_taskref)
                for x in values[0]
            )
    else:
        return _lol_product(head + (values[0],), values[1:], as_taskref=as_taskref)


def lol_tuples(head, ind, values, dummies):
    """List of list of tuple keys

    Parameters
    ----------
    head : tuple
        The known tuple so far
    ind : Iterable
        An iterable of indices not yet covered
    values : dict
        Known values for non-dummy indices
    dummies : dict
        Ranges of values for dummy indices

    Examples
    --------
    >>> lol_tuples(('x',), 'ij', {'i': 1, 'j': 0}, {})
    ('x', 1, 0)

    >>> lol_tuples(('x',), 'ij', {'i': 1}, {'j': range(3)})
    [('x', 1, 0), ('x', 1, 1), ('x', 1, 2)]

    >>> lol_tuples(('x',), 'ijk', {'i': 1}, {'j': [0, 1, 2], 'k': [0, 1]}) # doctest: +NORMALIZE_WHITESPACE
    [[('x', 1, 0, 0), ('x', 1, 0, 1)],
     [('x', 1, 1, 0), ('x', 1, 1, 1)],
     [('x', 1, 2, 0), ('x', 1, 2, 1)]]
    """
    if not ind:
        return head
    if ind[0] not in dummies:
        return lol_tuples(head + (values[ind[0]],), ind[1:], values, dummies)
    else:
        return [
            lol_tuples(head + (v,), ind[1:], values, dummies) for v in dummies[ind[0]]
        ]


def optimize_blockwise(graph, keys=()):
    """High level optimization of stacked Blockwise layers

    For operations that have multiple Blockwise operations one after the other, like
    ``x.T + 123`` we can fuse these into a single Blockwise operation.  This happens
    before any actual tasks are generated, and so can reduce overhead.

    This finds groups of Blockwise operations that can be safely fused, and then
    passes them to ``rewrite_blockwise`` for rewriting.

    Parameters
    ----------
    graph : HighLevelGraph
    keys : Iterable
        The keys of all outputs of all collections.
        Used to make sure that we don't fuse a layer needed by an output

    Returns
    -------
    HighLevelGraph

    See Also
    --------
    rewrite_blockwise
    """
    out = _optimize_blockwise(graph, keys=keys)
    while out.dependencies != graph.dependencies:
        graph = out
        out = _optimize_blockwise(graph, keys=keys)
    return out


def _optimize_blockwise(full_graph, keys=()):
    keep = {k[0] if type(k) is tuple else k for k in keys}
    layers = full_graph.layers
    dependents = reverse_dict(full_graph.dependencies)
    roots = {k for k in full_graph.layers if not dependents.get(k)}
    stack = list(roots)

    out = {}
    dependencies = {}
    seen = set()
    io_names = set()

    while stack:
        layer = stack.pop()
        if layer in seen or layer not in layers:
            continue
        seen.add(layer)

        # Outer loop walks through possible output Blockwise layers
        if isinstance(layers[layer], Blockwise):
            blockwise_layers = {layer}
            deps = set(blockwise_layers)
            io_names |= layers[layer].io_deps.keys()
            while deps:  # we gather as many sub-layers as we can
                dep = deps.pop()

                if dep not in layers:
                    stack.append(dep)
                    continue
                if not isinstance(layers[dep], Blockwise):
                    stack.append(dep)
                    continue
                if dep != layer and dep in keep:
                    stack.append(dep)
                    continue
                if layers[dep].concatenate != layers[layer].concatenate:
                    stack.append(dep)
                    continue
                if (
                    sum(k == dep for k, ind in layers[layer].indices if ind is not None)
                    > 1
                ):
                    stack.append(dep)
                    continue
                if blockwise_layers and not _can_fuse_annotations(
                    layers[next(iter(blockwise_layers))].annotations,
                    layers[dep].annotations,
                ):
                    stack.append(dep)
                    continue

                # passed everything, proceed
                blockwise_layers.add(dep)

                # traverse further to this child's children
                output_indices = set(layers[dep].output_indices)
                input_indices = {
                    i for _, ind in layers[dep].indices if ind for i in ind
                }
                is_io_superset = output_indices.issuperset(input_indices)
                for d in full_graph.dependencies.get(dep, ()):
                    # Don't allow reductions to proceed
                    if is_io_superset and len(dependents[d]) <= 1:
                        deps.add(d)
                    else:
                        stack.append(d)

            # Merge these Blockwise layers into one
            new_layer = rewrite_blockwise([layers[l] for l in blockwise_layers])
            out[layer] = new_layer

            # Get the new (external) dependencies for this layer.
            # This corresponds to the dependencies defined in
            # full_graph.dependencies and are not in blockwise_layers
            new_deps = set()
            for l in blockwise_layers:
                new_deps |= set(
                    {
                        d
                        for d in full_graph.dependencies[l]
                        if d not in blockwise_layers and d in full_graph.dependencies
                    }
                )
            for k, v in new_layer.indices:
                if v is None:
                    new_deps |= keys_in_tasks(full_graph.dependencies, [k])
                elif k not in io_names:
                    new_deps.add(k)
            dependencies[layer] = new_deps
        else:
            out[layer] = layers[layer]
            dependencies[layer] = full_graph.dependencies.get(layer, set())
            stack.extend(full_graph.dependencies.get(layer, ()))

    return HighLevelGraph(out, dependencies)


def _unique_dep(dep, ind):
    # Append blockwise index information to dependency name
    return dep + "_" + "_".join(str(i) for i in list(ind))


def _can_fuse_annotations(a: dict | None, b: dict | None) -> bool:
    """
    Treat the special annotation keys, as fusable since we can apply simple
    rules to capture their intent in a fused layer.
    """
    if a == b:
        return True

    if dask.config.get("optimization.annotations.fuse") is False:
        return False

    fusable = {"retries", "priority", "resources", "workers", "allow_other_workers"}
    if (not a or all(k in fusable for k in a)) and (
        not b or all(k in fusable for k in b)
    ):
        return True

    return False


def _fuse_annotations(*args: dict) -> dict:
    """
    Given an iterable of annotations dictionaries, fuse them according
    to some simple rules.
    """
    # First, do a basic dict merge -- we are presuming that these have already
    # been gated by `_can_fuse_annotations`.
    annotations = toolz.merge(*args)
    # Max of layer retries
    retries = [a["retries"] for a in args if "retries" in a]
    if retries:
        annotations["retries"] = max(retries)
    # Max of layer priorities
    priorities = [a["priority"] for a in args if "priority" in a]
    if priorities:
        annotations["priority"] = max(priorities)
    # Max of all the layer resources
    resources = [a["resources"] for a in args if "resources" in a]
    if resources:
        annotations["resources"] = toolz.merge_with(max, *resources)
    # Intersection of all the worker restrictions
    workers = [a["workers"] for a in args if "workers" in a]
    if workers:
        annotations["workers"] = list(set.intersection(*[set(w) for w in workers]))
    # More restrictive of allow_other_workers
    allow_other_workers = [
        a["allow_other_workers"] for a in args if "allow_other_workers" in a
    ]
    if allow_other_workers:
        annotations["allow_other_workers"] = all(allow_other_workers)

    return annotations


def rewrite_blockwise(inputs):
    """Rewrite a stack of Blockwise expressions into a single blockwise expression

    Given a set of Blockwise layers, combine them into a single layer.  The provided
    layers are expected to fit well together.  That job is handled by
    ``optimize_blockwise``

    Parameters
    ----------
    inputs : list[Blockwise]

    Returns
    -------
    blockwise: Blockwise

    See Also
    --------
    optimize_blockwise
    """
    if len(inputs) == 1:
        # Fast path: if there's only one input we can just use it as-is.
        return inputs[0]

    fused_annotations = _fuse_annotations(
        *[i.annotations for i in inputs if i.annotations]
    )
    inputs = {inp.output: inp for inp in inputs}
    dependencies = {
        inp.output: {d for d, v in inp.indices if v is not None and d in inputs}
        for inp in inputs.values()
    }
    dependents = reverse_dict(dependencies)

    new_index_iter = (
        c + (str(d) if d else "")  # A, B, ... A1, B1, ...
        for d in itertools.count()
        for c in "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    )

    [root] = [k for k, v in dependents.items() if not v]

    # Our final results.  These will change during fusion below
    indices = list(inputs[root].indices)
    new_axes = inputs[root].new_axes
    concatenate = inputs[root].concatenate
    task = inputs[root].task
    dsk = {task.key: task}

    changed = True
    while changed:
        changed = False
        for i, (dep, current_dep_indices) in enumerate(indices):
            if current_dep_indices is None:
                continue
            if dep not in inputs:
                continue

            changed = True

            # Change dep name to avoid fusing the same dep
            # (in different iteration orders) into a single
            # subgraph key/dependency
            # (see: https://github.com/dask/dask/issues/8535)
            local_dep = dep if dep == root else _unique_dep(dep, current_dep_indices)

            # Replace _n with dep name in existing tasks
            # (inc, _0) -> (inc, 'b')
            dsk = {
                k: v.substitute({blockwise_token(i): local_dep}) for k, v in dsk.items()
            }

            # Remove current input from input indices
            # [('a', 'i'), ('b', 'i')] -> [('a', 'i')]
            indices.pop(i)

            sub = {
                blockwise_token(i): blockwise_token(i - 1)
                for i in range(i + 1, len(indices) + 1)
            }
            dsk = {k: v.substitute(sub) for k, v in dsk.items()}

            # Change new input_indices to match give index from current computation
            # [('c', j')] -> [('c', 'i')]
            new_indices = inputs[dep].indices
            sub = dict(zip(inputs[dep].output_indices, current_dep_indices))
            contracted = {
                x
                for _, j in new_indices
                if j is not None
                for x in j
                if x not in inputs[dep].output_indices
            }
            extra = dict(zip(contracted, new_index_iter))
            sub.update(extra)
            new_indices = [(x, index_subs(j, sub)) for x, j in new_indices]

            # Update new_axes
            for k, v in inputs[dep].new_axes.items():
                new_axes[sub[k]] = v

            # Bump new inputs up in list
            sub = {}
            # Map from (id(key), inds or None) -> index in indices. Used to deduplicate indices.
            index_map = {(id(k), inds): n for n, (k, inds) in enumerate(indices)}
            for ii, index in enumerate(new_indices):
                id_key = (id(index[0]), index[1])
                if id_key in index_map:  # use old inputs if available
                    sub[blockwise_token(ii)] = blockwise_token(index_map[id_key])
                else:
                    index_map[id_key] = len(indices)
                    sub[blockwise_token(ii)] = blockwise_token(len(indices))
                    indices.append(index)
            if dep != local_dep:
                key = local_dep
            else:
                key = dep
            new_dsk = {key: inputs[dep].task.substitute(sub, key=key)}

            # indices.extend(new_indices)
            dsk.update(new_dsk)

    # De-duplicate indices like [(a, ij), (b, i), (a, ij)] -> [(a, ij), (b, i)]
    # Make sure that we map everything else appropriately as we remove inputs
    new_indices = []
    seen = {}
    sub = {}  # like {_0: _0, _1: _0, _2: _1}
    for i, x in enumerate(indices):
        if x[1] is not None and x in seen:
            sub[i] = seen[x]
        else:
            if x[1] is not None:
                seen[x] = len(new_indices)
            sub[i] = len(new_indices)
            new_indices.append(x)

    sub = {blockwise_token(k): blockwise_token(v) for k, v in sub.items()}
    dsk = {k: v.substitute(sub) for k, v in dsk.items() if k not in sub.keys()}
    task = Task.fuse(*dsk.values(), key=root)
    indices_check = {k for k, v in indices if v is not None}
    numblocks = toolz.merge([inp.numblocks for inp in inputs.values()])
    numblocks = {k: v for k, v in numblocks.items() if v is None or k in indices_check}

    # Update IO-dependency information
    io_deps = {}
    for v in inputs.values():
        io_deps.update(v.io_deps)

    return Blockwise(
        root,
        inputs[root].output_indices,
        task,
        new_indices,
        numblocks=numblocks,
        new_axes=new_axes,
        concatenate=concatenate,
        annotations=fused_annotations,
        io_deps=io_deps,
    )


@_deprecated()
def zero_broadcast_dimensions(lol, nblocks):
    """
    >>> lol = [('x', 1, 0), ('x', 1, 1), ('x', 1, 2)]
    >>> nblocks = (4, 1, 2)  # note singleton dimension in second place
    >>> lol = [[('x', 1, 0, 0), ('x', 1, 0, 1)],
    ...        [('x', 1, 1, 0), ('x', 1, 1, 1)],
    ...        [('x', 1, 2, 0), ('x', 1, 2, 1)]]

    >>> zero_broadcast_dimensions(lol, nblocks)  # doctest: +SKIP
    [[('x', 1, 0, 0), ('x', 1, 0, 1)],
     [('x', 1, 0, 0), ('x', 1, 0, 1)],
     [('x', 1, 0, 0), ('x', 1, 0, 1)]]

    See Also
    --------
    lol_tuples
    """
    f = lambda t: (t[0],) + tuple(0 if d == 1 else i for i, d in zip(t[1:], nblocks))
    return homogeneous_deepmap(f, lol)


def broadcast_dimensions(argpairs, numblocks, sentinels=(1, (1,)), consolidate=None):
    """Find block dimensions from arguments

    Parameters
    ----------
    argpairs : iterable
        name, ijk index pairs
    numblocks : dict
        maps {name: number of blocks}
    sentinels : iterable (optional)
        values for singleton dimensions
    consolidate : func (optional)
        use this to reduce each set of common blocks into a smaller set

    Examples
    --------
    >>> argpairs = [('x', 'ij'), ('y', 'ji')]
    >>> numblocks = {'x': (2, 3), 'y': (3, 2)}
    >>> broadcast_dimensions(argpairs, numblocks)
    {'i': 2, 'j': 3}

    Supports numpy broadcasting rules

    >>> argpairs = [('x', 'ij'), ('y', 'ij')]
    >>> numblocks = {'x': (2, 1), 'y': (1, 3)}
    >>> broadcast_dimensions(argpairs, numblocks)
    {'i': 2, 'j': 3}

    Works in other contexts too

    >>> argpairs = [('x', 'ij'), ('y', 'ij')]
    >>> d = {'x': ('Hello', 1), 'y': (1, (2, 3))}
    >>> broadcast_dimensions(argpairs, d)
    {'i': 'Hello', 'j': (2, 3)}
    """
    # List like [('i', 2), ('j', 1), ('i', 1), ('j', 2)]
    argpairs2 = [(a, ind) for a, ind in argpairs if ind is not None]
    L = toolz.concat(
        [
            zip(inds, dims)
            for (x, inds), (x, dims) in toolz.join(
                toolz.first, argpairs2, toolz.first, numblocks.items()
            )
        ]
    )

    g = toolz.groupby(0, L)
    g = {k: {d for i, d in v} for k, v in g.items()}

    g2 = {k: v - set(sentinels) if len(v) > 1 else v for k, v in g.items()}

    if consolidate:
        return toolz.valmap(consolidate, g2)

    if g2 and not set(map(len, g2.values())) == {1}:
        raise ValueError("Shapes do not align %s" % g)

    return toolz.valmap(toolz.first, g2)


def _make_dims(indices, numblocks, new_axes):
    """Returns a dictionary mapping between each index specified in
    `indices` and the number of output blocks for that indice.
    """
    dims = broadcast_dimensions(indices, numblocks)
    for k, v in new_axes.items():
        dims[k] = len(v) if isinstance(v, tuple) else 1
    return dims


def fuse_roots(graph: HighLevelGraph, keys: list):
    """
    Fuse nearby layers if they don't have dependencies

    Often Blockwise sections of the graph fill out all of the computation
    except for the initial data access or data loading layers::

      Large Blockwise Layer
        |       |       |
        X       Y       Z

    This can be troublesome because X, Y, and Z tasks may be executed on
    different machines, and then require communication to move around.

    This optimization identifies this situation, lowers all of the graphs to
    concrete dicts, and then calls ``fuse`` on them, with a width equal to the
    number of layers like X, Y, and Z.

    This is currently used within array and dataframe optimizations.

    Parameters
    ----------
    graph : HighLevelGraph
        The full graph of the computation
    keys : list
        The output keys of the computation, to be passed on to fuse

    See Also
    --------
    Blockwise
    fuse
    """
    layers = ensure_dict(graph.layers, copy=True)
    dependencies = ensure_dict(graph.dependencies, copy=True)
    dependents = reverse_dict(dependencies)

    for name, layer in graph.layers.items():
        deps = graph.dependencies[name]
        if (
            isinstance(layer, Blockwise)
            and len(deps) > 1
            and not any(dependencies[dep] for dep in deps)  # no need to fuse if 0 or 1
            and all(len(dependents[dep]) == 1 for dep in deps)
            and all(layer.annotations == graph.layers[dep].annotations for dep in deps)
        ):
            new = toolz.merge(layer, *[layers[dep] for dep in deps])
            new, _ = fuse(new, keys, ave_width=len(deps))

            for dep in deps:
                del layers[dep]
                del dependencies[dep]

            layers[name] = new
            dependencies[name] = set()

    return HighLevelGraph(layers, dependencies)


def _blockwise_unpack_collections_task_spec(expr):
    # FIXME This is a copy of the delayed.unpack_collections function with the
    # addition of the TaskSpec class. Eventually this should all be consolidated
    # but to reduce the number of changes we'll vendor this here
    # FIXME: There is also a dask.base version of unpack_collections that looks
    # similar but is different. At the very least the names should be fixed
    if isinstance(expr, Delayed):
        return TaskRef(expr._key), (expr,)

    if base.is_dask_collection(expr):
        if hasattr(expr, "optimize"):
            # Optimize dask-expr collections
            expr = expr.optimize()

        finalized = finalize(expr)
        return TaskRef(finalized._key), (finalized,)

    if type(expr) is type(iter(list())):
        expr = list(expr)
    elif type(expr) is type(iter(tuple())):
        expr = tuple(expr)
    elif type(expr) is type(iter(set())):
        expr = set(expr)

    typ = type(expr)

    if typ in (list, tuple, set):
        args, collections = utils.unzip(
            (_blockwise_unpack_collections_task_spec(e) for e in expr), 2
        )
        collections = tuple(toolz.unique(toolz.concat(collections), key=id))
        if not collections:
            return expr, ()
        args = List(*args)
        # Ensure output type matches input type
        if typ is not list:
            args = Task(None, typ, args)
        return args, collections

    if typ is dict:
        args, collections = _blockwise_unpack_collections_task_spec(
            [[k, v] for k, v in expr.items()]
        )
        if not collections:
            return expr, ()
        return Dict(args), collections

    if typ is slice:
        args, collections = _blockwise_unpack_collections_task_spec(
            [expr.start, expr.stop, expr.step]
        )
        if not collections:
            return expr, ()
        return Task(None, slice, *args), collections

    if is_dataclass(expr):
        args, collections = _blockwise_unpack_collections_task_spec(
            [
                [f.name, getattr(expr, f.name)]
                for f in fields(expr)
                if hasattr(expr, f.name)  # if init=False, field might not exist
            ]
        )
        if not collections:
            return expr, ()
        try:
            _fields = {
                f.name: getattr(expr, f.name)
                for f in fields(expr)
                if hasattr(expr, f.name)
            }
            replace(expr, **_fields)
        except (TypeError, ValueError) as e:
            if isinstance(e, ValueError) or "is declared with init=False" in str(e):
                raise ValueError(
                    f"Failed to unpack {typ} instance. "
                    "Note that using fields with `init=False` are not supported."
                ) from e
            else:
                raise TypeError(
                    f"Failed to unpack {typ} instance. "
                    "Note that using a custom __init__ is not supported."
                ) from e
        return Task(None, typ, **dict(args)), collections

    if utils.is_namedtuple_instance(expr):
        args, collections = _blockwise_unpack_collections_task_spec([v for v in expr])
        if not collections:
            return expr, ()
        return Task(None, typ, *args), collections

    return expr, ()
