"""
The code in this module is an experiment in going from N=1 to N=2 parallel computing frameworks in xarray.
It could later be used as the basis for a public interface allowing any N frameworks to interoperate with xarray,
but for now it is just a private experiment.
"""

from __future__ import annotations

import functools
import sys
from abc import ABC, abstractmethod
from collections.abc import Iterable, Sequence
from importlib.metadata import EntryPoint, entry_points
from typing import TYPE_CHECKING, Any, Callable, Generic, Protocol, TypeVar

import numpy as np

from xarray.core.utils import emit_user_level_warning
from xarray.namedarray.pycompat import is_chunked_array

if TYPE_CHECKING:
    from xarray.namedarray._typing import (
        _Chunks,
        _DType,
        _DType_co,
        _NormalizedChunks,
        _ShapeType,
        duckarray,
    )


class ChunkedArrayMixinProtocol(Protocol):
    def rechunk(self, chunks: Any, **kwargs: Any) -> Any: ...

    @property
    def dtype(self) -> np.dtype[Any]: ...

    @property
    def chunks(self) -> _NormalizedChunks: ...

    def compute(
        self, *data: Any, **kwargs: Any
    ) -> tuple[np.ndarray[Any, _DType_co], ...]: ...


T_ChunkedArray = TypeVar("T_ChunkedArray", bound=ChunkedArrayMixinProtocol)


@functools.lru_cache(maxsize=1)
def list_chunkmanagers() -> dict[str, ChunkManagerEntrypoint[Any]]:
    """
    Return a dictionary of available chunk managers and their ChunkManagerEntrypoint subclass objects.

    Returns
    -------
    chunkmanagers : dict
        Dictionary whose values are registered ChunkManagerEntrypoint subclass instances, and whose values
        are the strings under which they are registered.

    Notes
    -----
    # New selection mechanism introduced with Python 3.10. See GH6514.
    """
    if sys.version_info >= (3, 10):
        entrypoints = entry_points(group="xarray.chunkmanagers")
    else:
        entrypoints = entry_points().get("xarray.chunkmanagers", ())

    return load_chunkmanagers(entrypoints)


def load_chunkmanagers(
    entrypoints: Sequence[EntryPoint],
) -> dict[str, ChunkManagerEntrypoint[Any]]:
    """Load entrypoints and instantiate chunkmanagers only once."""

    loaded_entrypoints = {}
    for entrypoint in entrypoints:
        try:
            loaded_entrypoints[entrypoint.name] = entrypoint.load()
        except ModuleNotFoundError as e:
            emit_user_level_warning(
                f"Failed to load chunk manager entrypoint {entrypoint.name} due to {e}. Skipping.",
            )
            pass

    available_chunkmanagers = {
        name: chunkmanager()
        for name, chunkmanager in loaded_entrypoints.items()
        if chunkmanager.available
    }
    return available_chunkmanagers


def guess_chunkmanager(
    manager: str | ChunkManagerEntrypoint[Any] | None,
) -> ChunkManagerEntrypoint[Any]:
    """
    Get namespace of chunk-handling methods, guessing from what's available.

    If the name of a specific ChunkManager is given (e.g. "dask"), then use that.
    Else use whatever is installed, defaulting to dask if there are multiple options.
    """

    chunkmanagers = list_chunkmanagers()

    if manager is None:
        if len(chunkmanagers) == 1:
            # use the only option available
            manager = next(iter(chunkmanagers.keys()))
        else:
            # default to trying to use dask
            manager = "dask"

    if isinstance(manager, str):
        if manager not in chunkmanagers:
            raise ValueError(
                f"unrecognized chunk manager {manager} - must be one of: {list(chunkmanagers)}"
            )

        return chunkmanagers[manager]
    elif isinstance(manager, ChunkManagerEntrypoint):
        # already a valid ChunkManager so just pass through
        return manager
    else:
        raise TypeError(
            f"manager must be a string or instance of ChunkManagerEntrypoint, but received type {type(manager)}"
        )


def get_chunked_array_type(*args: Any) -> ChunkManagerEntrypoint[Any]:
    """
    Detects which parallel backend should be used for given set of arrays.

    Also checks that all arrays are of same chunking type (i.e. not a mix of cubed and dask).
    """

    # TODO this list is probably redundant with something inside xarray.apply_ufunc
    ALLOWED_NON_CHUNKED_TYPES = {int, float, np.ndarray}

    chunked_arrays = [
        a
        for a in args
        if is_chunked_array(a) and type(a) not in ALLOWED_NON_CHUNKED_TYPES
    ]

    # Asserts all arrays are the same type (or numpy etc.)
    chunked_array_types = {type(a) for a in chunked_arrays}
    if len(chunked_array_types) > 1:
        raise TypeError(
            f"Mixing chunked array types is not supported, but received multiple types: {chunked_array_types}"
        )
    elif len(chunked_array_types) == 0:
        raise TypeError("Expected a chunked array but none were found")

    # iterate over defined chunk managers, seeing if each recognises this array type
    chunked_arr = chunked_arrays[0]
    chunkmanagers = list_chunkmanagers()
    selected = [
        chunkmanager
        for chunkmanager in chunkmanagers.values()
        if chunkmanager.is_chunked_array(chunked_arr)
    ]
    if not selected:
        raise TypeError(
            f"Could not find a Chunk Manager which recognises type {type(chunked_arr)}"
        )
    elif len(selected) >= 2:
        raise TypeError(f"Multiple ChunkManagers recognise type {type(chunked_arr)}")
    else:
        return selected[0]


class ChunkManagerEntrypoint(ABC, Generic[T_ChunkedArray]):
    """
    Interface between a particular parallel computing framework and xarray.

    This abstract base class must be subclassed by libraries implementing chunked array types, and
    registered via the ``chunkmanagers`` entrypoint.

    Abstract methods on this class must be implemented, whereas non-abstract methods are only required in order to
    enable a subset of xarray functionality, and by default will raise a ``NotImplementedError`` if called.

    Attributes
    ----------
    array_cls
        Type of the array class this parallel computing framework provides.

        Parallel frameworks need to provide an array class that supports the array API standard.
        This attribute is used for array instance type checking at runtime.
    """

    array_cls: type[T_ChunkedArray]
    available: bool = True

    @abstractmethod
    def __init__(self) -> None:
        """Used to set the array_cls attribute at import time."""
        raise NotImplementedError()

    def is_chunked_array(self, data: duckarray[Any, Any]) -> bool:
        """
        Check if the given object is an instance of this type of chunked array.

        Compares against the type stored in the array_cls attribute by default.

        Parameters
        ----------
        data : Any

        Returns
        -------
        is_chunked : bool

        See Also
        --------
        dask.is_dask_collection
        """
        return isinstance(data, self.array_cls)

    @abstractmethod
    def chunks(self, data: T_ChunkedArray) -> _NormalizedChunks:
        """
        Return the current chunks of the given array.

        Returns chunks explicitly as a tuple of tuple of ints.

        Used internally by xarray objects' .chunks and .chunksizes properties.

        Parameters
        ----------
        data : chunked array

        Returns
        -------
        chunks : tuple[tuple[int, ...], ...]

        See Also
        --------
        dask.array.Array.chunks
        cubed.Array.chunks
        """
        raise NotImplementedError()

    @abstractmethod
    def normalize_chunks(
        self,
        chunks: _Chunks | _NormalizedChunks,
        shape: _ShapeType | None = None,
        limit: int | None = None,
        dtype: _DType | None = None,
        previous_chunks: _NormalizedChunks | None = None,
    ) -> _NormalizedChunks:
        """
        Normalize given chunking pattern into an explicit tuple of tuples representation.

        Exposed primarily because different chunking backends may want to make different decisions about how to
        automatically chunk along dimensions not given explicitly in the input chunks.

        Called internally by xarray.open_dataset.

        Parameters
        ----------
        chunks : tuple, int, dict, or string
            The chunks to be normalized.
        shape : Tuple[int]
            The shape of the array
        limit : int (optional)
            The maximum block size to target in bytes,
            if freedom is given to choose
        dtype : np.dtype
        previous_chunks : Tuple[Tuple[int]], optional
            Chunks from a previous array that we should use for inspiration when
            rechunking dimensions automatically.

        See Also
        --------
        dask.array.core.normalize_chunks
        """
        raise NotImplementedError()

    @abstractmethod
    def from_array(
        self, data: duckarray[Any, Any], chunks: _Chunks, **kwargs: Any
    ) -> T_ChunkedArray:
        """
        Create a chunked array from a non-chunked numpy-like array.

        Generally input should have a ``.shape``, ``.ndim``, ``.dtype`` and support numpy-style slicing.

        Called when the .chunk method is called on an xarray object that is not already chunked.
        Also called within open_dataset (when chunks is not None) to create a chunked array from
        an xarray lazily indexed array.

        Parameters
        ----------
        data : array_like
        chunks : int, tuple
            How to chunk the array.

        See Also
        --------
        dask.array.from_array
        cubed.from_array
        """
        raise NotImplementedError()

    def rechunk(
        self,
        data: T_ChunkedArray,
        chunks: _NormalizedChunks | tuple[int, ...] | _Chunks,
        **kwargs: Any,
    ) -> Any:
        """
        Changes the chunking pattern of the given array.

        Called when the .chunk method is called on an xarray object that is already chunked.

        Parameters
        ----------
        data : dask array
            Array to be rechunked.
        chunks :  int, tuple, dict or str, optional
            The new block dimensions to create. -1 indicates the full size of the
            corresponding dimension. Default is "auto" which automatically
            determines chunk sizes.

        Returns
        -------
        chunked array

        See Also
        --------
        dask.array.Array.rechunk
        cubed.Array.rechunk
        """
        return data.rechunk(chunks, **kwargs)

    @abstractmethod
    def compute(
        self, *data: T_ChunkedArray | Any, **kwargs: Any
    ) -> tuple[np.ndarray[Any, _DType_co], ...]:
        """
        Computes one or more chunked arrays, returning them as eager numpy arrays.

        Called anytime something needs to computed, including multiple arrays at once.
        Used by `.compute`, `.persist`, `.values`.

        Parameters
        ----------
        *data : object
            Any number of objects. If an object is an instance of the chunked array type, it is computed
            and the in-memory result returned as a numpy array. All other types should be passed through unchanged.

        Returns
        -------
        objs
            The input, but with all chunked arrays now computed.

        See Also
        --------
        dask.compute
        cubed.compute
        """
        raise NotImplementedError()

    @property
    def array_api(self) -> Any:
        """
        Return the array_api namespace following the python array API standard.

        See https://data-apis.org/array-api/latest/ . Currently used to access the array API function
        ``full_like``, which is called within the xarray constructors ``xarray.full_like``, ``xarray.ones_like``,
        ``xarray.zeros_like``, etc.

        See Also
        --------
        dask.array
        cubed.array_api
        """
        raise NotImplementedError()

    def reduction(
        self,
        arr: T_ChunkedArray,
        func: Callable[..., Any],
        combine_func: Callable[..., Any] | None = None,
        aggregate_func: Callable[..., Any] | None = None,
        axis: int | Sequence[int] | None = None,
        dtype: _DType_co | None = None,
        keepdims: bool = False,
    ) -> T_ChunkedArray:
        """
        A general version of array reductions along one or more axes.

        Used inside some reductions like nanfirst, which is used by ``groupby.first``.

        Parameters
        ----------
        arr : chunked array
            Data to be reduced along one or more axes.
        func : Callable(x_chunk, axis, keepdims)
            First function to be executed when resolving the dask graph.
            This function is applied in parallel to all original chunks of x.
            See below for function parameters.
        combine_func : Callable(x_chunk, axis, keepdims), optional
            Function used for intermediate recursive aggregation (see
            split_every below). If omitted, it defaults to aggregate_func.
        aggregate_func : Callable(x_chunk, axis, keepdims)
            Last function to be executed, producing the final output. It is always invoked, even when the reduced
            Array counts a single chunk along the reduced axes.
        axis : int or sequence of ints, optional
            Axis or axes to aggregate upon. If omitted, aggregate along all axes.
        dtype : np.dtype
            data type of output. This argument was previously optional, but
            leaving as ``None`` will now raise an exception.
        keepdims : boolean, optional
            Whether the reduction function should preserve the reduced axes,
            leaving them at size ``output_size``, or remove them.

        Returns
        -------
        chunked array

        See Also
        --------
        dask.array.reduction
        cubed.core.reduction
        """
        raise NotImplementedError()

    def scan(
        self,
        func: Callable[..., Any],
        binop: Callable[..., Any],
        ident: float,
        arr: T_ChunkedArray,
        axis: int | None = None,
        dtype: _DType_co | None = None,
        **kwargs: Any,
    ) -> T_ChunkedArray:
        """
        General version of a 1D scan, also known as a cumulative array reduction.

        Used in ``ffill`` and ``bfill`` in xarray.

        Parameters
        ----------
        func: callable
            Cumulative function like np.cumsum or np.cumprod
        binop: callable
            Associated binary operator like ``np.cumsum->add`` or ``np.cumprod->mul``
        ident: Number
            Associated identity like ``np.cumsum->0`` or ``np.cumprod->1``
        arr: dask Array
        axis: int, optional
        dtype: dtype

        Returns
        -------
        Chunked array

        See also
        --------
        dask.array.cumreduction
        """
        raise NotImplementedError()

    @abstractmethod
    def apply_gufunc(
        self,
        func: Callable[..., Any],
        signature: str,
        *args: Any,
        axes: Sequence[tuple[int, ...]] | None = None,
        keepdims: bool = False,
        output_dtypes: Sequence[_DType_co] | None = None,
        vectorize: bool | None = None,
        **kwargs: Any,
    ) -> Any:
        """
        Apply a generalized ufunc or similar python function to arrays.

        ``signature`` determines if the function consumes or produces core
        dimensions. The remaining dimensions in given input arrays (``*args``)
        are considered loop dimensions and are required to broadcast
        naturally against each other.

        In other terms, this function is like ``np.vectorize``, but for
        the blocks of chunked arrays. If the function itself shall also
        be vectorized use ``vectorize=True`` for convenience.

        Called inside ``xarray.apply_ufunc``, which is called internally for most xarray operations.
        Therefore this method must be implemented for the vast majority of xarray computations to be supported.

        Parameters
        ----------
        func : callable
            Function to call like ``func(*args, **kwargs)`` on input arrays
            (``*args``) that returns an array or tuple of arrays. If multiple
            arguments with non-matching dimensions are supplied, this function is
            expected to vectorize (broadcast) over axes of positional arguments in
            the style of NumPy universal functions [1]_ (if this is not the case,
            set ``vectorize=True``). If this function returns multiple outputs,
            ``output_core_dims`` has to be set as well.
        signature: string
            Specifies what core dimensions are consumed and produced by ``func``.
            According to the specification of numpy.gufunc signature [2]_
        *args : numeric
            Input arrays or scalars to the callable function.
        axes: List of tuples, optional, keyword only
            A list of tuples with indices of axes a generalized ufunc should operate on.
            For instance, for a signature of ``"(i,j),(j,k)->(i,k)"`` appropriate for
            matrix multiplication, the base elements are two-dimensional matrices
            and these are taken to be stored in the two last axes of each argument. The
            corresponding axes keyword would be ``[(-2, -1), (-2, -1), (-2, -1)]``.
            For simplicity, for generalized ufuncs that operate on 1-dimensional arrays
            (vectors), a single integer is accepted instead of a single-element tuple,
            and for generalized ufuncs for which all outputs are scalars, the output
            tuples can be omitted.
        keepdims: bool, optional, keyword only
            If this is set to True, axes which are reduced over will be left in the result as
            a dimension with size one, so that the result will broadcast correctly against the
            inputs. This option can only be used for generalized ufuncs that operate on inputs
            that all have the same number of core dimensions and with outputs that have no core
            dimensions , i.e., with signatures like ``"(i),(i)->()"`` or ``"(m,m)->()"``.
            If used, the location of the dimensions in the output can be controlled with axes
            and axis.
        output_dtypes : Optional, dtype or list of dtypes, keyword only
            Valid numpy dtype specification or list thereof.
            If not given, a call of ``func`` with a small set of data
            is performed in order to try to automatically determine the
            output dtypes.
        vectorize: bool, keyword only
            If set to ``True``, ``np.vectorize`` is applied to ``func`` for
            convenience. Defaults to ``False``.
        **kwargs : dict
            Extra keyword arguments to pass to `func`

        Returns
        -------
        Single chunked array or tuple of chunked arrays

        See Also
        --------
        dask.array.gufunc.apply_gufunc
        cubed.apply_gufunc

        References
        ----------
        .. [1] https://docs.scipy.org/doc/numpy/reference/ufuncs.html
        .. [2] https://docs.scipy.org/doc/numpy/reference/c-api/generalized-ufuncs.html
        """
        raise NotImplementedError()

    def map_blocks(
        self,
        func: Callable[..., Any],
        *args: Any,
        dtype: _DType_co | None = None,
        chunks: tuple[int, ...] | None = None,
        drop_axis: int | Sequence[int] | None = None,
        new_axis: int | Sequence[int] | None = None,
        **kwargs: Any,
    ) -> Any:
        """
        Map a function across all blocks of a chunked array.

        Called in elementwise operations, but notably not (currently) called within xarray.map_blocks.

        Parameters
        ----------
        func : callable
            Function to apply to every block in the array.
            If ``func`` accepts ``block_info=`` or ``block_id=``
            as keyword arguments, these will be passed dictionaries
            containing information about input and output chunks/arrays
            during computation. See examples for details.
        args : dask arrays or other objects
        dtype : np.dtype, optional
            The ``dtype`` of the output array. It is recommended to provide this.
            If not provided, will be inferred by applying the function to a small
            set of fake data.
        chunks : tuple, optional
            Chunk shape of resulting blocks if the function does not preserve
            shape. If not provided, the resulting array is assumed to have the same
            block structure as the first input array.
        drop_axis : number or iterable, optional
            Dimensions lost by the function.
        new_axis : number or iterable, optional
            New dimensions created by the function. Note that these are applied
            after ``drop_axis`` (if present).
        **kwargs :
            Other keyword arguments to pass to function. Values must be constants
            (not dask.arrays)

        See Also
        --------
        dask.array.map_blocks
        cubed.map_blocks
        """
        raise NotImplementedError()

    def blockwise(
        self,
        func: Callable[..., Any],
        out_ind: Iterable[Any],
        *args: Any,  # can't type this as mypy assumes args are all same type, but dask blockwise args alternate types
        adjust_chunks: dict[Any, Callable[..., Any]] | None = None,
        new_axes: dict[Any, int] | None = None,
        align_arrays: bool = True,
        **kwargs: Any,
    ) -> Any:
        """
        Tensor operation: Generalized inner and outer products.

        A broad class of blocked algorithms and patterns can be specified with a
        concise multi-index notation.  The ``blockwise`` function applies an in-memory
        function across multiple blocks of multiple inputs in a variety of ways.
        Many chunked array operations are special cases of blockwise including
        elementwise, broadcasting, reductions, tensordot, and transpose.

        Currently only called explicitly in xarray when performing multidimensional interpolation.

        Parameters
        ----------
        func : callable
            Function to apply to individual tuples of blocks
        out_ind : iterable
            Block pattern of the output, something like 'ijk' or (1, 2, 3)
        *args : sequence of Array, index pairs
            You may also pass literal arguments, accompanied by None index
            e.g. (x, 'ij', y, 'jk', z, 'i', some_literal, None)
        **kwargs : dict
            Extra keyword arguments to pass to function
        adjust_chunks : dict
            Dictionary mapping index to function to be applied to chunk sizes
        new_axes : dict, keyword only
            New indexes and their dimension lengths
        align_arrays: bool
            Whether or not to align chunks along equally sized dimensions when
            multiple arrays are provided.  This allows for larger chunks in some
            arrays to be broken into smaller ones that match chunk sizes in other
            arrays such that they are compatible for block function mapping. If
            this is false, then an error will be thrown if arrays do not already
            have the same number of blocks in each dimension.

        See Also
        --------
        dask.array.blockwise
        cubed.core.blockwise
        """
        raise NotImplementedError()

    def unify_chunks(
        self,
        *args: Any,  # can't type this as mypy assumes args are all same type, but dask unify_chunks args alternate types
        **kwargs: Any,
    ) -> tuple[dict[str, _NormalizedChunks], list[T_ChunkedArray]]:
        """
        Unify chunks across a sequence of arrays.

        Called by xarray.unify_chunks.

        Parameters
        ----------
        *args: sequence of Array, index pairs
            Sequence like (x, 'ij', y, 'jk', z, 'i')

        See Also
        --------
        dask.array.core.unify_chunks
        cubed.core.unify_chunks
        """
        raise NotImplementedError()

    def store(
        self,
        sources: T_ChunkedArray | Sequence[T_ChunkedArray],
        targets: Any,
        **kwargs: dict[str, Any],
    ) -> Any:
        """
        Store chunked arrays in array-like objects, overwriting data in target.

        This stores chunked arrays into object that supports numpy-style setitem
        indexing (e.g. a Zarr Store). Allows storing values chunk by chunk so that it does not have to
        fill up memory. For best performance you likely want to align the block size of
        the storage target with the block size of your array.

        Used when writing to any registered xarray I/O backend.

        Parameters
        ----------
        sources: Array or collection of Arrays
        targets: array-like or collection of array-likes
            These should support setitem syntax ``target[10:20] = ...``.
            If sources is a single item, targets must be a single item; if sources is a
            collection of arrays, targets must be a matching collection.
        kwargs:
            Parameters passed to compute/persist (only used if compute=True)

        See Also
        --------
        dask.array.store
        cubed.store
        """
        raise NotImplementedError()
