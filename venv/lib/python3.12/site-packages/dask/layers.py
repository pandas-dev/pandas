from __future__ import annotations

import functools
import operator
from itertools import product
from typing import TYPE_CHECKING

import tlz as toolz
from tlz.curried import map

from dask._task_spec import Task, TaskRef
from dask.base import tokenize
from dask.blockwise import Blockwise, BlockwiseDep, BlockwiseDepDict, blockwise_token
from dask.core import flatten
from dask.highlevelgraph import Layer
from dask.tokenize import normalize_token
from dask.utils import cached_cumsum

if TYPE_CHECKING:
    import numpy as np

#
##
###  General Utilities
##
#


class CallableLazyImport:
    """Function Wrapper for Lazy Importing.

    This Class should only be used when materializing a graph
    on a distributed scheduler.
    """

    def __init__(self, function_path):
        self.function_path = function_path

    def __call__(self, *args, **kwargs):
        from distributed.utils import import_term

        return import_term(self.function_path)(*args, **kwargs)


#
##
###  Array Layers & Utilities
##
#


class ArrayBlockwiseDep(BlockwiseDep):
    """
    Blockwise dep for array-likes, which only needs chunking
    information to compute its data.
    """

    chunks: tuple[tuple[int, ...], ...]
    numblocks: tuple[int, ...]
    produces_tasks: bool = False

    def __init__(self, chunks: tuple[tuple[int, ...], ...]):
        self.chunks = chunks
        self.numblocks = tuple(len(chunk) for chunk in chunks)
        self.produces_tasks = False

    def __getitem__(self, idx: tuple[int, ...]):
        raise NotImplementedError("Subclasses must implement __getitem__")


class ArrayChunkShapeDep(ArrayBlockwiseDep):
    """Produce chunk shapes given a chunk index"""

    def __getitem__(self, idx: tuple[int, ...]):
        return tuple(chunk[i] for i, chunk in zip(idx, self.chunks))


class ArraySliceDep(ArrayBlockwiseDep):
    """Produce slice(s) into the full-sized array given a chunk index"""

    starts: tuple[tuple[int, ...], ...]

    def __init__(self, chunks: tuple[tuple[int, ...], ...]):
        super().__init__(chunks)
        self.starts = tuple(cached_cumsum(c, initial_zero=True) for c in chunks)

    def __getitem__(self, idx: tuple):
        loc = tuple((start[i], start[i + 1]) for i, start in zip(idx, self.starts))
        return tuple(slice(*s, None) for s in loc)


class ArrayBlockIdDep(ArrayBlockwiseDep):
    def __getitem__(self, idx: tuple):
        if not isinstance(idx, tuple):
            # This happens if something else except graph generation calls us
            raise NotImplementedError("ArrayBlockIdDep requires a tuple index")
        return idx


class ArrayValuesDep(ArrayBlockwiseDep):
    def __init__(self, chunks: tuple[tuple[int, ...], ...], values: np.ndarray | dict):
        super().__init__(chunks)
        self.values = values

    def __getitem__(self, idx: tuple):
        return self.values[idx]


@normalize_token.register(ArraySliceDep)
def normalize_array_slice_dep(dep):
    return "ArraySliceDep", dep.chunks


@normalize_token.register(ArrayBlockIdDep)
def normalize_array_block_id_dep(dep):
    return "ArrayBlockIdDep", dep.chunks


@normalize_token.register(ArrayValuesDep)
def normalize_array_values_dep(dep):
    return "ArrayValuesDep", dep.chunks, dep.values


class ArrayOverlapLayer(Layer):
    """Simple HighLevelGraph array overlap layer.

    Lazily computed High-level graph layer for a array overlap operations.

    Parameters
    ----------
    name : str
        Name of new output overlap array.
    array : Dask array
    axes: Mapping
        Axes dictionary indicating overlap in each dimension,
        e.g. ``{'0': 1, '1': 1}``
    """

    def __init__(
        self,
        name,
        axes,
        chunks,
        numblocks,
        token,
    ):
        super().__init__()
        self.name = name
        self.axes = axes
        self.chunks = chunks
        self.numblocks = numblocks
        self.token = token
        self._cached_keys = None

    def __repr__(self):
        return f"ArrayOverlapLayer<name='{self.name}'"

    @property
    def _dict(self):
        """Materialize full dict representation"""
        if hasattr(self, "_cached_dict"):
            return self._cached_dict
        else:
            dsk = self._construct_graph()
            self._cached_dict = dsk
        return self._cached_dict

    def __getitem__(self, key):
        return self._dict[key]

    def __iter__(self):
        return iter(self._dict)

    def __len__(self):
        return len(self._dict)

    def is_materialized(self):
        return hasattr(self, "_cached_dict")

    def get_output_keys(self):
        return self.keys()  # FIXME! this implementation materializes the graph

    def _dask_keys(self):
        if self._cached_keys is not None:
            return self._cached_keys

        name, chunks, numblocks = self.name, self.chunks, self.numblocks

        def keys(*args):
            if not chunks:
                return [(name,)]
            ind = len(args)
            if ind + 1 == len(numblocks):
                result = [(name,) + args + (i,) for i in range(numblocks[ind])]
            else:
                result = [keys(*(args + (i,))) for i in range(numblocks[ind])]
            return result

        self._cached_keys = result = keys()
        return result

    def _construct_graph(self, deserializing=False):
        """Construct graph for a simple overlap operation."""
        axes = self.axes
        chunks = self.chunks
        name = self.name
        dask_keys = self._dask_keys()

        getitem_name = "getitem-" + self.token
        overlap_name = "overlap-" + self.token

        if deserializing:
            # Use CallableLazyImport objects to avoid importing dataframe
            # module on the scheduler
            concatenate_shaped = CallableLazyImport(
                "dask.array.core.concatenate_shaped"
            )
        else:
            # Not running on distributed scheduler - Use explicit functions
            from dask.array.core import concatenate_shaped

        dims = list(map(len, chunks))
        expand_key2 = functools.partial(
            _expand_keys_around_center, dims=dims, axes=axes
        )

        # Make keys for each of the surrounding sub-arrays
        interior_keys = toolz.pipe(
            dask_keys,
            flatten,
            map(expand_key2),
            map(lambda a: a[0]),
            map(flatten),
            toolz.concat,
            list,
        )
        interior_slices = {}
        overlap_blocks = {}
        for k in interior_keys:
            frac_slice = fractional_slice((name,) + k, axes)
            if frac_slice is False:
                continue
            if (name,) + k != frac_slice:
                interior_slices[(getitem_name,) + k] = frac_slice
            else:
                interior_slices[(getitem_name,) + k] = (name,) + k
                overlap_blocks[(overlap_name,) + k] = (
                    concatenate_shaped,
                    *(expand_key2((None,) + k, name=getitem_name)),
                )

        dsk = toolz.merge(interior_slices, overlap_blocks)
        return dsk


def _expand_keys_around_center(k, dims, name=None, axes=None):
    """Get all neighboring keys around center

    Parameters
    ----------
    k: Key
        The key around which to generate new keys
    dims: Sequence[int]
        The number of chunks in each dimension
    name: Option[str]
        The name to include in the output keys, or none to include no name
    axes: Dict[int, int]
        The axes active in the expansion.  We don't expand on non-active axes

    Examples
    --------
    >>> _expand_keys_around_center(('x', 2, 3), dims=[5, 5], name='y', axes={0: 1, 1: 1})  # noqa: E501 # doctest: +NORMALIZE_WHITESPACE
    ([('y', 1.1, 2.1), ('y', 1.1, 3), ('y', 1.1, 3.9), ('y',   2, 2.1), ('y',   2, 3), ('y',   2, 3.9), ('y', 2.9, 2.1), ('y', 2.9, 3), ('y', 2.9, 3.9)], (3, 3))

    >>> _expand_keys_around_center(('x', 0, 4), dims=[5, 5], name='y', axes={0: 1, 1: 1})  # noqa: E501 # doctest: +NORMALIZE_WHITESPACE
    ([('y',   0, 3.1), ('y',   0,   4), ('y', 0.9, 3.1), ('y', 0.9,   4)], (2, 2))
    """

    def convert_depth(depth):
        if not isinstance(depth, tuple):
            depth = (depth, depth)
        return depth

    def inds(i, ind, depth):
        depth = convert_depth(depth)
        rv = []
        if ind - 0.9 > 0 and depth[0] != 0:
            rv.append(ind - 0.9)
        rv.append(ind)
        if ind + 0.9 < dims[i] - 1 and depth[1] != 0:
            rv.append(ind + 0.9)
        return rv

    shape = []
    for i, ind in enumerate(k[1:]):
        depth = convert_depth(axes.get(i, 0))
        num = 1
        if ind > 0 and depth[0] != 0:
            num += 1
        if ind < dims[i] - 1 and depth[1] != 0:
            num += 1
        shape.append(num)

    def _valid_depth(depth):
        if isinstance(depth, tuple):
            return any(x != 0 for x in depth)
        else:
            return depth != 0

    args = [
        inds(i, ind, axes.get(i, 0)) if _valid_depth(axes.get(i, 0)) else [ind]
        for i, ind in enumerate(k[1:])
    ]
    if name is not None:
        args = [[name]] + args
    seq = list(product(*args))
    shape2 = tuple(
        d if _valid_depth(axes.get(i, 0)) else 1 for i, d in enumerate(shape)
    )
    return seq, shape2


def fractional_slice(task, axes):
    """

    >>> fractional_slice(('x', 5.1), {0: 2})
    (<built-in function getitem>, ('x', 5), (slice(-2, None, None),))

    >>> fractional_slice(('x', 3, 5.1), {0: 2, 1: 3})
    (<built-in function getitem>, ('x', 3, 5), (slice(None, None, None), slice(-3, None, None)))

    >>> fractional_slice(('x', 2.9, 5.1), {0: 2, 1: 3})
    (<built-in function getitem>, ('x', 3, 5), (slice(0, 2, None), slice(-3, None, None)))
    """
    rounded = (task[0],) + tuple(int(round(i)) for i in task[1:])

    index = []
    for i, (t, r) in enumerate(zip(task[1:], rounded[1:])):
        depth = axes.get(i, 0)
        if isinstance(depth, tuple):
            left_depth = depth[0]
            right_depth = depth[1]
        else:
            left_depth = depth
            right_depth = depth

        if t == r:
            index.append(slice(None, None, None))
        elif t < r and right_depth:
            index.append(slice(0, right_depth))
        elif t > r and left_depth:
            index.append(slice(-left_depth, None))
        else:
            return False
    index = tuple(index)

    if all(ind == slice(None, None, None) for ind in index):
        return task
    else:
        return (operator.getitem, rounded, index)


#
##
###  DataFrame Layers & Utilities
##
#


class DataFrameIOLayer(Blockwise):
    """DataFrame-based Blockwise Layer with IO

    Parameters
    ----------
    name : str
        Name to use for the constructed layer.
    columns : str, list or None
        Field name(s) to read in as columns in the output.
    inputs : list or BlockwiseDep
        List of arguments to be passed to ``io_func`` so
        that the materialized task to produce partition ``i``
        will be: ``(<io_func>, inputs[i])``.  Note that each
        element of ``inputs`` is typically a tuple of arguments.
    io_func : callable
        A callable function that takes in a single tuple
        of arguments, and outputs a DataFrame partition.
        Column projection will be supported for functions
        that satisfy the ``DataFrameIOFunction`` protocol.
    label : str (optional)
        String to use as a prefix in the place-holder collection
        name. If nothing is specified (default), "subset-" will
        be used.
    produces_tasks : bool (optional)
        Whether one or more elements of `inputs` is expected to
        contain a nested task. This argument in only used for
        serialization purposes, and will be deprecated in the
        future. Default is False.
    creation_info: dict (optional)
        Dictionary containing the callable function ('func'),
        positional arguments ('args'), and key-word arguments
        ('kwargs') used to produce the dask collection with
        this underlying ``DataFrameIOLayer``.
    annotations: dict (optional)
        Layer annotations to pass through to Blockwise.
    """

    def __init__(
        self,
        name,
        columns,
        inputs,
        io_func,
        label=None,
        produces_tasks=False,
        creation_info=None,
        annotations=None,
    ):
        self.name = name
        self._columns = columns
        self.inputs = inputs
        self.io_func = io_func
        self.label = label
        self.produces_tasks = produces_tasks
        self.annotations = annotations
        self.creation_info = creation_info

        if not isinstance(inputs, BlockwiseDep):
            # Define mapping between key index and "part"
            io_arg_map = BlockwiseDepDict(
                {(i,): inp for i, inp in enumerate(self.inputs)},
                produces_tasks=self.produces_tasks,
            )
        else:
            io_arg_map = inputs

        # Use Blockwise initializer
        task = Task(self.name, io_func, TaskRef(blockwise_token(0)))
        super().__init__(
            output=self.name,
            output_indices="i",
            task=task,
            indices=[(io_arg_map, "i")],
            numblocks={},
            annotations=annotations,
        )

    @property
    def columns(self):
        """Current column projection for this layer"""
        return self._columns

    def project_columns(self, columns):
        """Produce a column projection for this IO layer.
        Given a list of required output columns, this method
        returns the projected layer.
        """
        from dask.dataframe.io.utils import DataFrameIOFunction

        columns = list(columns)

        if self.columns is None or set(self.columns).issuperset(columns):
            # Apply column projection in IO function.
            # Must satisfy `DataFrameIOFunction` protocol
            if isinstance(self.io_func, DataFrameIOFunction):
                io_func = self.io_func.project_columns(columns)
            else:
                io_func = self.io_func

            layer = DataFrameIOLayer(
                (self.label or "subset") + "-" + tokenize(self.name, columns),
                columns,
                self.inputs,
                io_func,
                label=self.label,
                produces_tasks=self.produces_tasks,
                annotations=self.annotations,
            )
            return layer
        else:
            # Default behavior
            return self

    def __repr__(self):
        return "DataFrameIOLayer<name='{}', n_parts={}, columns={}>".format(
            self.name, len(self.inputs), self.columns
        )
