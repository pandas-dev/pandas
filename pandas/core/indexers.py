"""
Low-dependency indexing utilities.
"""
import numpy as np

from pandas._typing import AnyArrayLike

from pandas.core.dtypes.common import is_list_like
from pandas.core.dtypes.generic import ABCIndexClass, ABCSeries

# -----------------------------------------------------------
# Indexer Identification


def is_list_like_indexer(key) -> bool:
    """
    Check if we have a list-like indexer that is *not* a NamedTuple.

    Parameters
    ----------
    key : object

    Returns
    -------
    bool
    """
    # allow a list_like, but exclude NamedTuples which can be indexers
    return is_list_like(key) and not (isinstance(key, tuple) and type(key) is not tuple)


def is_scalar_indexer(indexer, arr_value) -> bool:
    """
    Return True if we are all scalar indexers.

    Returns
    -------
    bool
    """
    if arr_value.ndim == 1:
        if not isinstance(indexer, tuple):
            indexer = tuple([indexer])
            return any(isinstance(idx, np.ndarray) and len(idx) == 0 for idx in indexer)
    return False


def is_empty_indexer(indexer, arr_value: np.ndarray) -> bool:
    """
    Check if we have an empty indexer.

    Parameters
    ----------
    indexer : object
    arr_value : np.ndarray

    Returns
    -------
    bool
    """
    if is_list_like(indexer) and not len(indexer):
        return True
    if arr_value.ndim == 1:
        if not isinstance(indexer, tuple):
            indexer = tuple([indexer])
        return any(isinstance(idx, np.ndarray) and len(idx) == 0 for idx in indexer)
    return False


# -----------------------------------------------------------
# Indexer Validation


def check_setitem_lengths(indexer, value, values) -> None:
    """
    Validate that value and indexer are the same length.

    An special-case is allowed for when the indexer is a boolean array
    and the number of true values equals the length of ``value``. In
    this case, no exception is raised.

    Parameters
    ----------
    indexer : sequence
        Key for the setitem.
    value : array-like
        Value for the setitem.
    values : array-like
        Values being set into.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        When the indexer is an ndarray or list and the lengths don't match.
    """
    # boolean with truth values == len of the value is ok too
    if isinstance(indexer, (np.ndarray, list)):
        if is_list_like(value) and len(indexer) != len(value):
            if not (
                isinstance(indexer, np.ndarray)
                and indexer.dtype == np.bool_
                and len(indexer[indexer]) == len(value)
            ):
                raise ValueError(
                    "cannot set using a list-like indexer "
                    "with a different length than the value"
                )

    elif isinstance(indexer, slice):
        # slice
        if is_list_like(value) and len(values):
            if len(value) != length_of_indexer(indexer, values):
                raise ValueError(
                    "cannot set using a slice indexer with a "
                    "different length than the value"
                )


def validate_indices(indices: np.ndarray, n: int) -> None:
    """
    Perform bounds-checking for an indexer.

    -1 is allowed for indicating missing values.

    Parameters
    ----------
    indices : ndarray
    n : int
        Length of the array being indexed.

    Raises
    ------
    ValueError

    Examples
    --------
    >>> validate_indices([1, 2], 3)
    # OK
    >>> validate_indices([1, -2], 3)
    ValueError
    >>> validate_indices([1, 2, 3], 3)
    IndexError
    >>> validate_indices([-1, -1], 0)
    # OK
    >>> validate_indices([0, 1], 0)
    IndexError
    """
    if len(indices):
        min_idx = indices.min()
        if min_idx < -1:
            msg = f"'indices' contains values less than allowed ({min_idx} < -1)"
            raise ValueError(msg)

        max_idx = indices.max()
        if max_idx >= n:
            raise IndexError("indices are out-of-bounds")


# -----------------------------------------------------------
# Indexer Conversion


def maybe_convert_indices(indices, n: int):
    """
    Attempt to convert indices into valid, positive indices.

    If we have negative indices, translate to positive here.
    If we have indices that are out-of-bounds, raise an IndexError.

    Parameters
    ----------
    indices : array-like
        Array of indices that we are to convert.
    n : int
        Number of elements in the array that we are indexing.

    Returns
    -------
    array-like
        An array-like of positive indices that correspond to the ones
        that were passed in initially to this function.

    Raises
    ------
    IndexError
        One of the converted indices either exceeded the number of,
        elements (specified by `n`), or was still negative.
    """
    if isinstance(indices, list):
        indices = np.array(indices)
        if len(indices) == 0:
            # If `indices` is empty, np.array will return a float,
            # and will cause indexing errors.
            return np.empty(0, dtype=np.intp)

    mask = indices < 0
    if mask.any():
        indices = indices.copy()
        indices[mask] += n

    mask = (indices >= n) | (indices < 0)
    if mask.any():
        raise IndexError("indices are out-of-bounds")
    return indices


# -----------------------------------------------------------
# Unsorted


def length_of_indexer(indexer, target=None) -> int:
    """
    Return the length of a single non-tuple indexer which could be a slice.

    Returns
    -------
    int
    """
    if target is not None and isinstance(indexer, slice):
        target_len = len(target)
        start = indexer.start
        stop = indexer.stop
        step = indexer.step
        if start is None:
            start = 0
        elif start < 0:
            start += target_len
        if stop is None or stop > target_len:
            stop = target_len
        elif stop < 0:
            stop += target_len
        if step is None:
            step = 1
        elif step < 0:
            start, stop = stop + 1, start + 1
            step = -step
        return (stop - start + step - 1) // step
    elif isinstance(indexer, (ABCSeries, ABCIndexClass, np.ndarray, list)):
        return len(indexer)
    elif not is_list_like_indexer(indexer):
        return 1
    raise AssertionError("cannot find the length of the indexer")


def check_bool_array_indexer(array: AnyArrayLike, mask: AnyArrayLike) -> np.ndarray:
    """
    Check if `mask` is a valid boolean indexer for `array`.

    `array` and `mask` are checked to have the same length, and the
    dtype is validated.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    array : array
        The array that's being masked.
    mask : array
        The boolean array that's masking.

    Returns
    -------
    numpy.ndarray
        The validated boolean mask.

    Raises
    ------
    IndexError
        When the lengths don't match.
    ValueError
        When `mask` cannot be converted to a bool-dtype ndarray.

    See Also
    --------
    api.types.is_bool_dtype : Check if `key` is of boolean dtype.

    Examples
    --------
    A boolean ndarray is returned when the arguments are all valid.

    >>> mask = pd.array([True, False])
    >>> arr = pd.array([1, 2])
    >>> pd.api.extensions.check_bool_array_indexer(arr, mask)
    array([ True, False])

    An IndexError is raised when the lengths don't match.

    >>> mask = pd.array([True, False, True])
    >>> pd.api.extensions.check_bool_array_indexer(arr, mask)
    Traceback (most recent call last):
    ...
    IndexError: Item wrong length 3 instead of 2.

    A ValueError is raised when the mask cannot be converted to
    a bool-dtype ndarray.

    >>> mask = pd.array([True, pd.NA])
    >>> pd.api.extensions.check_bool_array_indexer(arr, mask)
    Traceback (most recent call last):
    ...
    ValueError: cannot convert to bool numpy array in presence of missing values
    """
    result = np.asarray(mask, dtype=bool)
    # GH26658
    if len(result) != len(array):
        raise IndexError(f"Item wrong length {len(result)} instead of {len(array)}.")
    return result
