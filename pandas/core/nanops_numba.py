# Module for implementing nanops in numba
from collections.abc import Callable
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    TypeVar,
)

import numba as nb
from numba.core.extending import overload
from numba.typed import List as NumbaList
import numpy as np

if TYPE_CHECKING:
    from pandas._typing import (
        AxisInt,
        npt,
    )

T = TypeVar("T")
R = TypeVar("R")
F = TypeVar("F", bound=Callable[..., Any])


MIN_INT = np.iinfo(np.int64).min


def is_null(x: Any) -> bool:
    """
    Check if a value is considered null/NA.

    Parameters
    ----------
    x : scalar
        Value to check

    Returns
    -------
    bool
        True if value is null, False otherwise

    Notes
    -----
    This function is overloaded with specialized implementations for
    various numeric types via Numba's overload mechanism.
    """
    dtype = np.asarray(x).dtype
    if np.issubdtype(dtype, np.float64):
        return np.isnan(x)

    elif np.issubdtype(dtype, np.int64):
        return x == MIN_INT

    else:
        return False


@overload(is_null)  # type: ignore[misc]
def jit_is_null(x: Any) -> Any:
    if isinstance(x, nb.types.Integer):

        def is_null(x: Any) -> bool:
            return x == MIN_INT

        return is_null

    elif isinstance(x, nb.types.Float) or isinstance(x, float):

        def is_null(x: Any) -> bool:
            return np.isnan(x)

        return is_null

    elif isinstance(x, nb.types.Boolean):

        def is_null(x: Any) -> bool:
            return False

        return is_null

    else:
        return np.isnat


def null_for_np_type(dtype: np.dtype) -> np.ndarray:
    """
    Return the appropriate null value for a given numpy dtype.

    Parameters
    ----------
    dtype : np.dtype
        NumPy data type to get null value for

    Returns
    -------
    scalar
        NaT for datetime/timedelta types, np.nan for other types

    Notes
    -----
    For datetime64 and timedelta64 dtypes (kind 'm' or 'M'), returns
    the appropriate NaT (Not a Time) value. For all other dtypes,
    returns np.nan.
    """
    if dtype.kind in "mM":
        return np.array(["NaT"], dtype=dtype)[0]
    else:
        return np.array(np.nan)


@nb.njit(nogil=True)
def _get_initial_value(
    arr: np.ndarray, skipna: bool = True, mask: np.ndarray | None = None
) -> tuple[int, Any]:
    """
    Find the first value in an array for use in a reduction, and its location.
    If skipna is True then we find the first non-null value,
    otherwise just the first value.
    If the array is empty, return a null value.

    Parameters
    ----------
    arr : array-like
        Array to search for non-null values
    mask: np.ndarray
        Boolean mask where True indicates null values in an extension array

    Returns
    -------
    tuple
        (index, value) of first non-null value, or (-1, np.nan) if all values are null

    Notes
    -----
    This function is JIT-compiled with Numba for performance.
    """
    if not skipna:
        if mask is not None and mask[0]:
            return 0, np.nan
        elif is_null(arr[0]):
            return 0, np.nan
        else:
            return 0, arr[0]
    elif mask is not None:
        for i, x in enumerate(arr):
            if not mask[i]:
                return i, x
    else:
        for i, x in enumerate(arr):
            if not is_null(x):
                return i, x
    return -1, np.nan


_SCALAR_SIGNATURES = [
    "float64(float64, float64)",
    "uint64(uint64, uint64)",
    "int64(int64, int64)",
]


def _njit_scalar_reduce(func: Callable) -> Any:
    """
    Decorator to create numba-compiled scalar reduction functions.

    Parameters
    ----------
    func : callable
        Scalar reduction function taking two arguments

    Returns
    -------
    staticmethod
        Numba-compiled version of the function with standard signatures

    Notes
    -----
    This decorator compiles the function with predefined signatures for
    common numeric types (float64, uint64, int64) and enables nogil mode
    for better performance in multithreaded environments.
    """
    return staticmethod(nb.njit(_SCALAR_SIGNATURES, nogil=True)(func))


class NumbaReductionOps:
    """
    Collection of numba implementations of scalar reduction ops"""

    @_njit_scalar_reduce
    def count(x: Any, y: Any) -> Any:
        return x + 1

    @_njit_scalar_reduce
    def min(x: Any, y: Any) -> Any:
        return x if x <= y else y

    @_njit_scalar_reduce
    def max(x: Any, y: Any) -> Any:
        return x if x >= y else y

    @_njit_scalar_reduce
    def sum(x: Any, y: Any) -> Any:
        return x + y

    @_njit_scalar_reduce
    def prod(x: Any, y: Any) -> Any:
        return x + y

    @_njit_scalar_reduce
    def sum_square(x: Any, y: Any) -> Any:
        return x + float(y) ** 2

    @_njit_scalar_reduce
    def any(x: Any, y: Any) -> Any:
        return x or y

    @_njit_scalar_reduce
    def all(x: Any, y: Any) -> Any:
        return x and y


@nb.njit(nogil=True)
def _nb_reduce_single_arr(
    reduce_func: Callable,
    arr: np.ndarray,
    skipna: bool = True,
    find_initial_value: bool = True,
    mask: np.ndarray | None = None,
) -> tuple[float | int, int]:
    """
    Apply a reduction function to a numpy array, with NA/null handling.
    Returns the count of non-nulls as well as the reduction.

    Parameters
    ----------
    reduce_func : callable
        Function that combines two values (e.g., min, max, sum)
    arr : array-like
        Array to reduce
    skipna : bool, default True
        Whether to skip NA/null values
    initial_value:
        Initial_value for each reduction. Should be 0 or None.
        If None, we find the first_non_null value before commencing the reduction

    Returns
    -------
    scalar
        Result of the reduction operation

    Notes
    -----
    This function is JIT-compiled with Numba for performance.
    """
    if not find_initial_value:
        initial_value = 0.0
        initial_loc = -1
        count = 0

    else:
        # find the initial non-null value to pass through the reduction
        # If the array is empty then this returns the type-appropriate null
        initial_loc, initial_value = _get_initial_value(arr, skipna=skipna, mask=mask)
        if is_null(
            initial_value
        ):  # skipna is False and initial value is null, or all values are null
            return np.nan, 0
        else:
            count = 1

    result = initial_value
    if mask is not None:
        for i in range(initial_loc + 1, len(arr)):
            if mask[i]:
                if skipna:
                    continue
                else:
                    return np.nan, count

            result = reduce_func(result, arr[i])
            count += 1

    else:
        for x in arr[initial_loc + 1 :]:
            if is_null(x):
                if skipna:
                    continue
                else:
                    return np.nan, count

            result = reduce_func(result, x)
            count += 1

    return result, count


@nb.njit(nogil=True, parallel=True)
def _nb_reduce_arr_list_in_parallel(
    reduce_func: Callable,
    arr_list: NumbaList[np.ndarray] | np.ndarray,
    target: np.ndarray,
    mask_list: np.ndarray | None,
    skipna: bool = True,
    find_initial_value: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    counts = np.zeros(len(arr_list), dtype=np.int64)
    for i in nb.prange(len(arr_list)):
        arr = arr_list[i]
        if mask_list is None:
            mask = None
        else:
            mask = mask_list[i]

        target[i], counts[i] = _nb_reduce_single_arr(
            reduce_func,
            arr,
            skipna=skipna,
            find_initial_value=find_initial_value,
            mask=mask,
        )

    return target, counts


def reduction_return_type_and_empty_result_for_op_and_type(
    dtype: np.dtype, op: Literal["count", "min", "max", "sum", "sum_square", "mean"]
) -> tuple[np.dtype, Any]:
    """
    Determine the return dtype and empty result value for a reduction operation.

    Parameters
    ----------
    dtype : np.dtype
        Input array dtype
    op : {"count", "min", "max", "sum", "sum_square", "mean"}
        Reduction operation to perform

    Returns
    -------
    tuple
        (return_dtype, empty_result_value) for the given operation and input dtype

    Notes
    -----
    This function defines the type promotion rules and empty result values
    for various reduction operations on different input dtypes.
    """
    if op == "count":
        return np.dtype(np.int64), 0
    elif op in ("min", "max"):
        return dtype, null_for_np_type(dtype)
    elif op == "sum":
        match dtype.kind:
            case "f":
                return np.dtype("float64"), 0.0
            case "u":
                return np.dtype("uint64"), 0
            case "m":
                return dtype, np.timedelta64(0)
            case "M":
                return dtype, np.datetime64(0, "ns")
            case _:
                return np.dtype("int64"), 0
    elif op == "mean":
        # always use floats for mean/var/std calculation to avoid overflow
        if dtype.kind in "mM":
            return dtype, null_for_np_type(dtype)
        else:
            return np.dtype("float64"), np.nan
    elif op == "sum_square":
        return np.dtype("float64"), np.nan
    else:
        raise ValueError(
            'op must be one of ["count", "min", "max", "sum", "sum_square"]'
        )


def _nullify_below_mincount(
    result: np.ndarray, count: np.ndarray, min_count: int
) -> np.ndarray:
    """
    Set result elements to null where count is below minimum threshold.

    Parameters
    ----------
    result : np.ndarray
        Result array to modify
    count : np.ndarray
        Count of valid values for each result element
    min_count : int
        Minimum number of non-null values required

    Returns
    -------
    np.ndarray
        Modified result array with nullified values

    Notes
    -----
    For unsigned integer dtypes, uses MIN_INT as null value.
    For all other dtypes, uses np.nan as null value.
    """
    if result.dtype.kind in "ui":
        null = MIN_INT
    else:
        null = np.nan  # type: ignore[assignment]

    result[count < min_count] = null

    return result


def _reduce_empty_array(
    op: Literal["count", "min", "max", "sum", "sum_square", "mean"],
    values: np.ndarray,
    axis: int | None,
    min_count: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    return_type, empty_result = reduction_return_type_and_empty_result_for_op_and_type(
        values.dtype, op
    )
    if min_count > 0:
        empty_result = null_for_np_type(return_type)
    if values.ndim == 2 and axis is not None:
        n = values.shape[1 - axis]
        return np.full(n, empty_result), np.zeros(n, dtype=int)
    else:
        return empty_result, np.array(0)


def _arr_to_iterable_for_multi_threading(
    arr: np.ndarray,
    multi_threading: bool,
    axis: int | None,
    mask: np.ndarray | None = None,
) -> tuple[NumbaList | np.ndarray, NumbaList | np.ndarray, int]:
    """
    Split arrays into chunks for potential parallel processing in reduction operations.
    In cases where the array is 2-D and there is sufficient length in the axis
    of reduction we simply return the array (or its transpose) which serves
    as the iterable. This is an important optimization as converting to a very
    long list is expensive, e.g. if reducing along the column axis of a DataFrame
    with many rows. If the reduction axis is shorter than the number of threads
    we chunk the other axis as well.

    Parameters
    ----------
    arr : np.ndarray
        Input array to be chunked. Must be 1D or 2D.
    multi_threading : bool
        If True, split array into multiple chunks for parallel processing.
        If False, return single chunk (no parallelization).
    axis : int or None
        Reduction axis. For 2D arrays:
        - axis=0: transpose array so reduction operates along columns
        - axis=1: keep array as-is, reduction operates along rows
        - axis=None:  keep array as-is, reduction operates along rows, then columns
    mask : np.ndarray, optional
        Boolean mask indicating null values. If provided, will be split
        consistently with values array.

    Returns
    -------
    tuple
        - arr_list : NumbaList
            List of array chunks ready for parallel processing
        - mask_list : NumbaList
            List of corresponding mask chunks (empty if mask=None)
        - final_length : int
            Length of the final reduction dimension. 0 for 1D arrays,
            number of columns/rows for 2D arrays.

    Notes
    -----
    Thread count is determined automatically based on array size when
    multi_threading=True, with a maximum of 6 threads and minimum of 1.
    Arrays smaller than 1 million elements use single threading.

    For 1D arrays, the array is split into n_threads chunks along axis 0.
    For 2D arrays, the array is either transposed (axis=0) or used as-is
    (axis=1) to prepare for row-wise or column-wise reductions.

    Raises
    ------
    ValueError
        If input array has more than 2 dimensions.

    Examples
    --------
    >>> arr = np.array([[1, 2, 3], [4, 5, 6]])
    >>> arr_list, mask_list, final_length = _arr_to_iterable_for_multi_threading(
    ...     arr, multi_threading=False, axis=0
    ... )
    >>> final_length
    3
    >>> len(arr_list)
    3
    """
    ndim = arr.ndim
    if multi_threading:
        # TODO: be smarter about this choice. numba is handling the distribution
        # of the compute so don't need to worry about setting it too high
        max_n_threads = min(6, int(arr.size // 1e6))
        n_threads = max(1, max_n_threads)
    else:
        n_threads = 1

    if mask is None:
        mask = np.array([])

    if ndim == 1:
        axis = 0
        arr_list = NumbaList(np.array_split(arr, n_threads))
        mask_list = NumbaList(np.array_split(mask, n_threads))
        final_length = 0

    elif ndim == 2:
        if axis == 0:
            arr_list = arr.T
            mask_list = mask.T
        else:
            arr_list = arr
            mask_list = mask
        final_length = 0 if axis is None else len(arr_list)
    else:
        raise ValueError("Only arrays of 1 or 2 dimensions are supported")

    if len(arr_list) < n_threads:
        arr_list = NumbaList(
            [chunk for row in arr_list for chunk in np.array_split(row, n_threads)]
        )
        mask_list = NumbaList(
            [chunk for row in mask_list for chunk in np.array_split(row, n_threads)]
        )

    return arr_list, mask_list, final_length


def _reduce_chunked_results(
    op: str,
    chunk_results: np.ndarray,
    counts: np.ndarray,
    final_length: int,
    **kwargs: Any,
) -> tuple[np.ndarray, np.ndarray]:
    chunk_reducer = (
        NumbaReductionOps.sum if op == "sum_square" else getattr(NumbaReductionOps, op)
    )
    result: np.ndarray

    if final_length == 0:
        # we chunked and want to reduce both axes
        result, _ = _nb_reduce_single_arr(chunk_reducer, arr=chunk_results, **kwargs)  # type: ignore[assignment]
        result = np.array(result)
        count = counts.sum()
    elif len(chunk_results) > final_length:
        # We chunked along both axes and want to reduce a single axis
        arr_list = np.array_split(chunk_results, final_length)
        target = np.zeros(final_length, dtype=np.float64)
        result, _ = _nb_reduce_arr_list_in_parallel(
            chunk_reducer, arr_list=arr_list, mask_list=None, target=target, **kwargs
        )
        count = [c.sum() for c in np.array_split(counts, final_length)]
    else:
        result, count = chunk_results, counts

    result, count = map(np.atleast_1d, (result, count))

    return result, count


def nb_reduce(
    op: Literal["count", "min", "max", "sum", "sum_square", "mean"],
    values: np.ndarray,
    axis: int | None = None,
    skipna: bool = True,
    min_count: int = 0,
    mask: np.ndarray | None = None,
    multi_threading: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Apply a reduction operation to a numpy array using Numba-accelerated functions.

    Parameters
    ----------
    op : {"count", "min", "max", "sum", "sum_square"}
        The reduction operation to perform
    arr : np.ndarray
        Input array to reduce (1D or 2D)
    axis : int, optional
        Axis along which to perform the reduction. If None, reduces over all elements.
        For 2D arrays, axis=0 reduces along rows, axis=1 reduces along columns.
    skipna : bool, default True
        Whether to skip NA/null values during reduction
    multi_threading : bool, default True
        Whether to use parallel processing by splitting array into chunks (1D only)

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Two-element tuple containing:
        - Reduction result: scalar for 1D or axis=None, array for 2D with specified axis
        - Count of non-null values processed (scalar or array matching result shape)

    Notes
    -----
    This function provides high-performance reduction operations by leveraging
    Numba's JIT compilation and optional parallel processing. For 1D arrays with
    multi_threading=True, the array is split into chunks processed in parallel.

    Supports arrays up to 2 dimensions:
    - 1D arrays: reduces to scalar
    - 2D arrays: reduces along specified axis or to scalar if axis=None

    The function handles null values according to the skipna parameter:
    - If skipna=True: null values are ignored in the reduction
    - If skipna=False: any null value causes early termination

    For integer arrays, MIN_INT is used as the null sentinel value.
    For float arrays, NaN is used as the null value.

    Examples
    --------
    >>> import numpy as np
    >>> # 1D array reduction
    >>> arr = np.array([1.0, 2.0, np.nan, 4.0])
    >>> result, count = nb_reduce("sum", arr, skipna=True)
    >>> result, count
    (np.float64(7.0), 3)

    >>> # 2D array reduction along axis
    >>> arr_2d = np.array([[1.0, 2.0], [3.0, np.nan]])
    >>> result, count = nb_reduce("sum", arr_2d, axis=0, skipna=True)
    >>> result, count
    (array([4., 2.]), array([2, 1]))
    """
    values = np.asarray(values)
    if values.dtype.kind in ("b", "u", "i") and mask is None:
        skipna = False

    elif values.dtype.kind == "c":
        kwargs = locals().copy()
        real_piece, count = nb_reduce(**(kwargs | {"values": values.real}))  # type: ignore[attr-defined]
        imaginary_piece, count = nb_reduce(**(kwargs | {"values": values.imag}))  # type: ignore[attr-defined]
        return real_piece + 1j * imaginary_piece, count

    if values.size == 0:
        return _reduce_empty_array(op, values, axis=axis, min_count=min_count)

    return_dtype, _ = reduction_return_type_and_empty_result_for_op_and_type(
        values.dtype, op
    )

    is_timelike = values.dtype.kind in "mM"
    if is_timelike:
        values = values.view(int)

    ndim = np.ndim(values)
    if not (axis is None or axis < ndim):
        raise ValueError(f"axis {axis} out-of-bounds for array of dimension {ndim}")

    return_scalar = ndim == 1 or axis is None

    reduce_op = "sum" if op == "mean" else op
    reduce_func = getattr(NumbaReductionOps, reduce_op)

    arr_list, mask_list, final_length = _arr_to_iterable_for_multi_threading(
        values, multi_threading=multi_threading, axis=axis, mask=mask
    )

    kwargs = {
        "skipna": skipna,
        "find_initial_value": "sum" not in reduce_op,
    }
    target = np.zeros(len(arr_list), dtype=np.float64)

    result, count = _nb_reduce_arr_list_in_parallel(
        reduce_func=reduce_func,
        arr_list=arr_list,
        target=target,
        mask_list=None if mask is None else mask_list,
        **kwargs,
    )
    result, count = _reduce_chunked_results(
        reduce_op,
        result,
        count,
        final_length=final_length,
        **kwargs,
    )

    if op in ["mean", "sum_square"]:
        if op == "mean":
            with np.errstate(invalid="ignore", divide="ignore"):
                result = result / count
        if not skipna:
            # null integers need to be nullified here as dividing by the count
            # causes MIN_INT results to increase
            null = count < values.shape[axis or 0]
            result[null] = np.nan

    if min_count > 0:
        result = _nullify_below_mincount(result, count, min_count)

    if return_dtype.kind in "mM":
        result = _cast_to_timelike(result, return_dtype)

    elif return_dtype.kind == "f" or not np.isnan(result).any():
        result = result.astype(return_dtype, copy=False)

    if return_scalar:
        result, count = result[0], int(count[0])  # type: ignore[assignment]
        result = result.dtype.type(result)

    return result, count  # type: ignore[return-value]


def _cast_to_timelike(arr: np.ndarray, to_dtype: np.dtype) -> np.ndarray:
    """
    Convert a float array to timelike (datetime/timedelta) dtype.

    Parameters
    ----------
    arr : np.ndarray
        Float array to convert
    to_dtype : np.dtype
        Target datetime or timedelta dtype

    Returns
    -------
    np.ndarray
        Array converted to timelike dtype with NaN values replaced by MIN_INT

    Notes
    -----
    This function is used to convert float arrays back to timelike dtypes
    after reduction operations. NaN values are replaced with MIN_INT before
    conversion to preserve null representation in integer-based time types.
    """
    isnan = np.isnan(arr)
    if isnan.any():
        arr[isnan] = MIN_INT
    arr = arr.astype(int, copy=False).astype(to_dtype, copy=False)

    return arr


def nanmax(
    values: "np.ndarray",
    *,
    axis: "AxisInt | None" = None,
    skipna: "bool" = True,
    min_count: "int" = 0,
    mask: "npt.NDArray[np.bool_] | None" = None,
    multi_threading: bool = True,
) -> np.ndarray:
    return nb_reduce("max", **locals())[0]


def nanmin(
    values: "np.ndarray",
    *,
    axis: "AxisInt | None" = None,
    skipna: "bool" = True,
    min_count: "int" = 0,
    mask: "npt.NDArray[np.bool_] | None" = None,
    multi_threading: bool = True,
) -> np.ndarray:
    return nb_reduce("min", **locals())[0]


def nansum(
    values: "np.ndarray",
    *,
    axis: "AxisInt | None" = None,
    skipna: "bool" = True,
    min_count: "int" = 0,
    mask: "npt.NDArray[np.bool_] | None" = None,
    multi_threading: bool = True,
) -> np.ndarray:
    return nb_reduce("sum", **locals())[0]


def nanmean(
    values: "np.ndarray",
    *,
    axis: "AxisInt | None" = None,
    skipna: "bool" = True,
    min_count: "int" = 0,
    mask: "npt.NDArray[np.bool_] | None" = None,
    multi_threading: bool = True,
) -> np.ndarray:
    return nb_reduce("mean", **locals())[0]


def _nanvar_std_sem(
    values: "np.ndarray",
    *,
    axis: "AxisInt | None" = None,
    skipna: "bool" = True,
    ddof: int = 1,
    mask: "npt.NDArray[np.bool_] | None" = None,
    multi_threading: bool = True,
    std: bool = False,
    sem: bool = False,
) -> np.ndarray:
    kwargs = locals().copy()
    if values.dtype.kind == "c":
        return sum(  # type: ignore[return-value]
            _nanvar_std_sem(**(kwargs | {"values": x}))
            for x in [values.real, values.imag]  # type: ignore[attr-defined]
        )

    dtype = values.dtype
    is_timelike = dtype.kind in "mM"

    del kwargs["ddof"], kwargs["std"], kwargs["sem"]
    mean, count = nb_reduce("mean", **kwargs)

    if np.ndim(mean) == 1 and len(mean) == len(values):
        kwargs["values"] = (values.T - mean).T
    else:
        kwargs["values"] = values - mean

    sum_of_squares, count = nb_reduce("sum_square", **kwargs)

    if np.ndim(mean) == 0:
        if is_null(mean) or count <= ddof:
            return np.timedelta64(MIN_INT) if is_timelike else np.nan  # type: ignore[return-value]

    result = sum_of_squares / (count - ddof)

    if std or sem:
        result = result.astype(float, copy=False) ** 0.5
        if sem:
            result = result / np.sqrt(count)
        if is_timelike:
            result = np.array(result).astype(dtype.str.replace("M", "m"))
            if np.ndim(result) == 0:
                result = np.atleast_1d(result)[0]

    return result


def nanvar(
    values: "np.ndarray",
    *,
    axis: "AxisInt | None" = None,
    skipna: "bool" = True,
    ddof: int = 1,
    mask: "npt.NDArray[np.bool_] | None" = None,
    multi_threading: bool = True,
) -> np.ndarray:
    return _nanvar_std_sem(**locals())


def nanstd(
    values: "np.ndarray",
    *,
    axis: "AxisInt | None" = None,
    skipna: "bool" = True,
    ddof: int = 1,
    mask: "npt.NDArray[np.bool_] | None" = None,
    multi_threading: bool = True,
) -> np.ndarray:
    return _nanvar_std_sem(**locals(), std=True)


def nansem(
    values: "np.ndarray",
    *,
    axis: "AxisInt | None" = None,
    skipna: "bool" = True,
    ddof: int = 1,
    mask: "npt.NDArray[np.bool_] | None" = None,
    multi_threading: bool = True,
) -> np.ndarray:
    return _nanvar_std_sem(**locals(), sem=True)
