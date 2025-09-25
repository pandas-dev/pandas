# Module for implementing nanops in numba
from typing import Any, Callable, Tuple, TypeVar, Literal, Optional

import numpy as np
import numba as nb

from numba.core.extending import overload
from numba.typed import List as NumbaList


T = TypeVar("T")
R = TypeVar("R")
F = TypeVar("F", bound=Callable[..., Any])


MIN_INT = np.iinfo(np.int64).min


def is_null(x):
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


@overload(is_null)
def jit_is_null(x):
    if isinstance(x, nb.types.Integer):

        def is_null(x):
            return x == MIN_INT

        return is_null

    elif isinstance(x, nb.types.Float) or isinstance(x, float):

        def is_null(x):
            return np.isnan(x)

        return is_null

    elif isinstance(x, nb.types.Boolean):

        def is_null(x):
            return False

        return is_null

    else:

        return np.isnat


def null_for_np_type(dtype):
    if dtype.kind in "mM":
        return np.array(["NaT"], dtype=dtype)[0]
    else:
        return np.nan


@nb.njit(nogil=True)
def _get_initial_value(
    arr, skipna: bool = True, mask: Optional[np.ndarray] = None
) -> Tuple[int, T]:
    """
    Find the first value in an array for use in a reduction, and its location.
    If skipna is True then we find the first non-null value, other wise just the first value.
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


def _njit_scalar_reduce(func):
    return staticmethod(nb.njit(_SCALAR_SIGNATURES, nogil=True)(func))


class NumbaReductionOps:
    """
    Collection of numba implementations of scalar reduction ops"""

    @_njit_scalar_reduce
    def count(x, y):
        return x + 1

    @_njit_scalar_reduce
    def min(x, y):
        return x if x <= y else y

    @_njit_scalar_reduce
    def max(x, y):
        return x if x >= y else y

    @_njit_scalar_reduce
    def sum(x, y):
        return x + y

    @_njit_scalar_reduce
    def prod(x, y):
        return x + y

    @_njit_scalar_reduce
    def sum_square(x, y):
        return x + float(y) ** 2

    @_njit_scalar_reduce
    def any(x, y):
        return x or y

    @_njit_scalar_reduce
    def all(x, y):
        return x and y


@nb.njit(nogil=True)
def _nb_reduce_single_arr(
    reduce_func: Callable,
    arr: np.ndarray,
    skipna: bool = True,
    find_initial_value: bool = True,
    mask: Optional[np.ndarray] = None,
) -> Tuple[float | int, int]:
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

    arr = arr[initial_loc + 1 :]
    result = initial_value
    if mask is not None:
        mask = mask[initial_loc + 1 :]
        for x, mask_i in zip(arr, mask):
            if mask_i:
                if skipna:
                    continue
                else:
                    return np.nan, count

            result = reduce_func(result, x)
            count += 1

    else:
        for x in arr:
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
    arr_list: NumbaList[np.ndarray] | np.ndarray,  # type: ignore
    target: np.ndarray,
    mask_list: Optional[np.ndarray] | np.ndarray,
    skipna: bool = True,
    find_initial_value: bool = True,
) -> Tuple[float | int, int]:
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
    dtype, op: Literal["count", "min", "max", "sum", "sum_square", "mean"]
):
    if op == "count":
        return np.int64, 0
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


def _nullify_below_mincount(result, count, min_count):
    if result.dtype.kind in "ui":
        null = MIN_INT
    else:
        null = np.nan

    result[count < min_count] = null

    return result


def _reduce_empty_array(op, values: np.ndarray, axis: int, min_count: int = 0):
    return_type, empty_result = reduction_return_type_and_empty_result_for_op_and_type(
        values.dtype, op
    )
    if min_count > 0:
        empty_result = null_for_np_type(return_type)
    if values.ndim == 2 and axis is not None:
        n = values.shape[1 - axis]
        return np.full(n, empty_result), np.zeros(n)
    else:
        return empty_result, 0


def _chunk_arr_into_arr_list(
    values: np.ndarray,
    multi_threading: bool,
    axis: Optional[int],
    mask: Optional[np.ndarray] = None,
) -> NumbaList:
    ndim = values.ndim
    if multi_threading:
        # TODO: be smarter about this choice. numba is handling the distribution of the compute
        # so don't need to worry about setting it too high
        max_n_threads = min(6, values.size // 1e6)
        n_threads = max(1, max_n_threads)
    else:
        n_threads = 1

    if mask is None:
        mask = np.array([])

    if ndim == 1:
        axis = 0
        arr_list = np.array_split(values, n_threads)  # type: ignore
        mask_list = np.array_split(mask, n_threads)  # type: ignore
        final_length = 0

    elif ndim == 2:
        if axis == 0:
            arr_list = values.T
            mask_list = mask.T
        else:
            arr_list = values
            mask_list = mask
        final_length = 0 if axis is None else len(arr_list)
    else:
        raise ValueError("Only arrays of 1 or 2 dimensions are supported")

    if len(arr_list) < n_threads:
        arr_list = [
            chunk for row in arr_list for chunk in np.array_split(row, n_threads)
        ]
        mask_list = [
            chunk for row in mask_list for chunk in np.array_split(row, n_threads)
        ]

    return arr_list, mask_list, final_length


def _reduce_chunked_results(
    op,
    chunk_results: np.ndarray,
    counts: np.ndarray,
    final_length: int,
    return_dtype: np.dtype,
    **kwargs,
):
    chunk_reducer = (
        NumbaReductionOps.sum if op == "sum_square" else getattr(NumbaReductionOps, op)
    )

    if final_length == 0:
        # we chunked and want to reduce both axes
        result, _ = _nb_reduce_single_arr(chunk_reducer, arr=chunk_results, **kwargs)
        count = counts.sum()
    elif len(chunk_results) > final_length:
        # We chunked along both axes and want to reduce a single axis
        arr_list = np.array_split(chunk_results, final_length)
        target = np.zeros(final_length, dtype=np.float64)  # type: ignore
        result, _ = _nb_reduce_arr_list_in_parallel(
            chunk_reducer, arr_list=arr_list, mask_list=None, target=target, **kwargs
        )
        count = [c.sum() for c in np.array_split(counts, final_length)]
    else:
        result, count = chunk_results, counts

    result, count = map(np.atleast_1d, (result, count))

    return result, count


def nb_reduce(
    op: Literal["count", "min", "max", "sum", "sum_square"],
    values: np.ndarray,
    axis: Optional[int] = None,
    skipna: bool = True,
    min_count: int = 0,
    mask: Optional[np.ndarray] = None,
    multi_threading: bool = True,
) -> Tuple[np.ndarray, np.ndarray]:
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
    tuple[float | int | np.ndarray, int | np.ndarray]
        Two-element tuple containing:
        - Reduction result (scalar for 1D or axis=None, array for 2D with specified axis)
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
    (7.0, 3)

    >>> # 2D array reduction along axis
    >>> arr_2d = np.array([[1.0, 2.0], [3.0, np.nan]])
    >>> result, count = nb_reduce("sum", arr_2d, axis=0, skipna=True)
    >>> result, count
    (array([4.0, 2.0]), array([2, 1]))
    """
    values = np.asarray(values)
    if values.dtype.kind in "bui" and mask is None:
        skipna = False

    elif values.dtype.kind == "c":
        kwargs = locals().copy()
        real_piece, count = nb_reduce(**(kwargs | {"values": values.real}))
        imaginary_piece, count = nb_reduce(**(kwargs | {"values": values.imag}))
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

    arr_list, mask_list, final_length = _chunk_arr_into_arr_list(
        values, multi_threading=multi_threading, axis=axis, mask=mask
    )

    kwargs = {
        "skipna": skipna,
        "find_initial_value": "sum" not in reduce_op,
    }
    target = np.zeros(len(arr_list), dtype=np.float64)  # type: ignore

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
        return_dtype=return_dtype,
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
        result, count = result[0], int(count[0])
        result = result.dtype.type(result)

    return result, count


def _cast_to_timelike(arr, to_dtype):
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
        return sum(
            _nanvar_std_sem(**(kwargs | {"values": x}))
            for x in [values.real, values.imag]
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
            return np.timedelta64(MIN_INT) if is_timelike else np.nan

    result = sum_of_squares / (count - ddof)

    if std or sem:
        result = result.astype(float, copy=False) ** 0.5
        if sem:
            result = result / np.sqrt(count)
        if is_timelike:
            result = np.array(result).astype(dtype.str.replace("M", "m"))
            if np.ndim(result) == 0:
                result = np.timedelta64(result)

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
