"""
Low-dependency indexing utilities.
"""
import warnings

import numpy as np

from pandas._typing import Any, AnyArrayLike

from pandas.core.dtypes.common import (
    is_array_like,
    is_bool_dtype,
    is_integer_dtype,
    is_list_like,
)
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


def deprecate_ndim_indexing(result):
    """
    Helper function to raise the deprecation warning for multi-dimensional
    indexing on 1D Series/Index.

    GH#27125 indexer like idx[:, None] expands dim, but we cannot do that
    and keep an index, so we currently return ndarray, which is deprecated
    (Deprecation GH#30588).
    """
    if np.ndim(result) > 1:
        warnings.warn(
            "Support for multi-dimensional indexing (e.g. `index[:, None]`) "
            "on an Index is deprecated and will be removed in a future "
            "version.  Convert to a numpy array before indexing instead.",
            DeprecationWarning,
            stacklevel=3,
        )


# -----------------------------------------------------------
# Public indexer validation


def check_array_indexer(array: AnyArrayLike, indexer: Any) -> Any:
    """
    Check if `indexer` is a valid array indexer for `array`.

    For a boolean mask, `array` and `indexer` are checked to have the same
    length. The dtype is validated, and if it is an integer or boolean
    ExtensionArray, it is checked if there are missing values present, and
    it is converted to the appropriate numpy array. Other dtypes will raise
    an error.

    Non-array indexers (integer, slice, Ellipsis, tuples, ..) are passed
    through as is.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    array : array-like
        The array that is being indexed (only used for the length).
    indexer : array-like or list-like
        The array-like that's used to index. List-like input that is not yet
        a numpy array or an ExtensionArray is converted to one. Other input
        types are passed through as is.

    Returns
    -------
    numpy.ndarray
        The validated indexer as a numpy array that can be used to index.

    Raises
    ------
    IndexError
        When the lengths don't match.
    ValueError
        When `indexer` cannot be converted to a numpy ndarray to index
        (e.g. presence of missing values).

    See Also
    --------
    api.types.is_bool_dtype : Check if `key` is of boolean dtype.

    Examples
    --------
    When checking a boolean mask, a boolean ndarray is returned when the
    arguments are all valid.

    >>> mask = pd.array([True, False])
    >>> arr = pd.array([1, 2])
    >>> pd.api.indexers.check_array_indexer(arr, mask)
    array([ True, False])

    An IndexError is raised when the lengths don't match.

    >>> mask = pd.array([True, False, True])
    >>> pd.api.indexers.check_array_indexer(arr, mask)
    Traceback (most recent call last):
    ...
    IndexError: Boolean index has wrong length: 3 instead of 2.

    A ValueError is raised when the mask cannot be converted to
    a bool-dtype ndarray.

    >>> mask = pd.array([True, pd.NA])
    >>> pd.api.indexers.check_array_indexer(arr, mask)
    Traceback (most recent call last):
    ...
    ValueError: Cannot mask with a boolean indexer containing NA values

    A numpy boolean mask will get passed through (if the length is correct):

    >>> mask = np.array([True, False])
    >>> pd.api.indexers.check_array_indexer(arr, mask)
    array([ True, False])

    Similarly for integer indexers, an integer ndarray is returned when it is
    a valid indexer, otherwise an error is  (for integer indexers, a matching
    length is not required):

    >>> indexer = pd.array([0, 2], dtype="Int64")
    >>> arr = pd.array([1, 2, 3])
    >>> pd.api.indexers.check_array_indexer(arr, indexer)
    array([0, 2])

    >>> indexer = pd.array([0, pd.NA], dtype="Int64")
    >>> pd.api.indexers.check_array_indexer(arr, indexer)
    Traceback (most recent call last):
    ...
    ValueError: Cannot index with an integer indexer containing NA values

    For non-integer/boolean dtypes, an appropriate error is raised:

    >>> indexer = np.array([0., 2.], dtype="float64")
    >>> pd.api.indexers.check_array_indexer(arr, indexer)
    Traceback (most recent call last):
    ...
    IndexError: arrays used as indices must be of integer or boolean type
    """
    from pandas.core.construction import array as pd_array

    # whathever is not an array-like is returned as-is (possible valid array
    # indexers that are not array-like: integer, slice, Ellipsis, None)
    # In this context, tuples are not considered as array-like, as they have
    # a specific meaning in indexing (multi-dimensional indexing)
    if is_list_like(indexer):
        if isinstance(indexer, tuple):
            return indexer
    else:
        return indexer

    # convert list-likes to array
    if not is_array_like(indexer):
        indexer = pd_array(indexer)
        if len(indexer) == 0:
            # empty list is converted to float array by pd.array
            indexer = np.array([], dtype=np.intp)

    dtype = indexer.dtype
    if is_bool_dtype(dtype):
        try:
            indexer = np.asarray(indexer, dtype=bool)
        except ValueError:
            raise ValueError("Cannot mask with a boolean indexer containing NA values")

        # GH26658
        if len(indexer) != len(array):
            raise IndexError(
                f"Boolean index has wrong length: "
                f"{len(indexer)} instead of {len(array)}"
            )
    elif is_integer_dtype(dtype):
        try:
            indexer = np.asarray(indexer, dtype=np.intp)
        except ValueError:
            raise ValueError(
                "Cannot index with an integer indexer containing NA values"
            )
    else:
        raise IndexError("arrays used as indices must be of integer or boolean type")

    return indexer
