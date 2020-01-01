"""
Misc tools for implementing data structures

Note: pandas.core.common is *not* part of the public API.
"""

import collections
from collections import abc
from datetime import datetime, timedelta
from functools import partial
import inspect
from typing import Any, Collection, Iterable, Union

import numpy as np

from pandas._libs import lib, tslibs
from pandas._typing import T

from pandas.core.dtypes.cast import construct_1d_object_array_from_listlike
from pandas.core.dtypes.common import (
    is_array_like,
    is_bool_dtype,
    is_extension_array_dtype,
    is_integer,
)
from pandas.core.dtypes.generic import ABCIndex, ABCIndexClass, ABCSeries
from pandas.core.dtypes.inference import _iterable_not_string
from pandas.core.dtypes.missing import isna, isnull, notnull  # noqa


class SettingWithCopyError(ValueError):
    pass


class SettingWithCopyWarning(Warning):
    pass


def flatten(l):
    """
    Flatten an arbitrarily nested sequence.

    Parameters
    ----------
    l : sequence
        The non string sequence to flatten

    Notes
    -----
    This doesn't consider strings sequences.

    Returns
    -------
    flattened : generator
    """
    for el in l:
        if _iterable_not_string(el):
            for s in flatten(el):
                yield s
        else:
            yield el


def consensus_name_attr(objs):
    name = objs[0].name
    for obj in objs[1:]:
        try:
            if obj.name != name:
                name = None
        except ValueError:
            name = None
    return name


def maybe_box(indexer, values, obj, key):

    # if we have multiples coming back, box em
    if isinstance(values, np.ndarray):
        return obj[indexer.get_loc(key)]

    # return the value
    return values


def maybe_box_datetimelike(value):
    # turn a datetime like into a Timestamp/timedelta as needed

    if isinstance(value, (np.datetime64, datetime)):
        value = tslibs.Timestamp(value)
    elif isinstance(value, (np.timedelta64, timedelta)):
        value = tslibs.Timedelta(value)

    return value


values_from_object = lib.values_from_object


def is_bool_indexer(key: Any) -> bool:
    """
    Check whether `key` is a valid boolean indexer.

    Parameters
    ----------
    key : Any
        Only list-likes may be considered boolean indexers.
        All other types are not considered a boolean indexer.
        For array-like input, boolean ndarrays or ExtensionArrays
        with ``_is_boolean`` set are considered boolean indexers.

    Returns
    -------
    bool

    Raises
    ------
    ValueError
        When the array is an object-dtype ndarray or ExtensionArray
        and contains missing values.
    """
    na_msg = "cannot index with vector containing NA / NaN values"
    if isinstance(key, (ABCSeries, np.ndarray, ABCIndex)) or (
        is_array_like(key) and is_extension_array_dtype(key.dtype)
    ):
        if key.dtype == np.object_:
            key = np.asarray(values_from_object(key))

            if not lib.is_bool_array(key):
                if isna(key).any():
                    raise ValueError(na_msg)
                return False
            return True
        elif is_bool_dtype(key.dtype):
            # an ndarray with bool-dtype by definition has no missing values.
            # So we only need to check for NAs in ExtensionArrays
            if is_extension_array_dtype(key.dtype):
                if np.any(key.isna()):
                    raise ValueError(na_msg)
            return True
    elif isinstance(key, list):
        try:
            arr = np.asarray(key)
            return arr.dtype == np.bool_ and len(arr) == len(key)
        except TypeError:  # pragma: no cover
            return False

    return False


def cast_scalar_indexer(val):
    """
    To avoid numpy DeprecationWarnings, cast float to integer where valid.

    Parameters
    ----------
    val : scalar

    Returns
    -------
    outval : scalar
    """
    # assumes lib.is_scalar(val)
    if lib.is_float(val) and val == int(val):
        return int(val)
    return val


def not_none(*args):
    """
    Returns a generator consisting of the arguments that are not None.
    """
    return (arg for arg in args if arg is not None)


def any_none(*args):
    """
    Returns a boolean indicating if any argument is None.
    """
    return any(arg is None for arg in args)


def all_none(*args):
    """
    Returns a boolean indicating if all arguments are None.
    """
    return all(arg is None for arg in args)


def any_not_none(*args):
    """
    Returns a boolean indicating if any argument is not None.
    """
    return any(arg is not None for arg in args)


def all_not_none(*args):
    """
    Returns a boolean indicating if all arguments are not None.
    """
    return all(arg is not None for arg in args)


def count_not_none(*args):
    """
    Returns the count of arguments that are not None.
    """
    return sum(x is not None for x in args)


def try_sort(iterable):
    listed = list(iterable)
    try:
        return sorted(listed)
    except TypeError:
        return listed


def asarray_tuplesafe(values, dtype=None):

    if not (isinstance(values, (list, tuple)) or hasattr(values, "__array__")):
        values = list(values)
    elif isinstance(values, ABCIndexClass):
        return values.values

    if isinstance(values, list) and dtype in [np.object_, object]:
        return construct_1d_object_array_from_listlike(values)

    result = np.asarray(values, dtype=dtype)

    if issubclass(result.dtype.type, str):
        result = np.asarray(values, dtype=object)

    if result.ndim == 2:
        # Avoid building an array of arrays:
        values = [tuple(x) for x in values]
        result = construct_1d_object_array_from_listlike(values)

    return result


def index_labels_to_array(labels, dtype=None):
    """
    Transform label or iterable of labels to array, for use in Index.

    Parameters
    ----------
    dtype : dtype
        If specified, use as dtype of the resulting array, otherwise infer.

    Returns
    -------
    array
    """
    if isinstance(labels, (str, tuple)):
        labels = [labels]

    if not isinstance(labels, (list, np.ndarray)):
        try:
            labels = list(labels)
        except TypeError:  # non-iterable
            labels = [labels]

    labels = asarray_tuplesafe(labels, dtype=dtype)

    return labels


def maybe_make_list(obj):
    if obj is not None and not isinstance(obj, (tuple, list)):
        return [obj]
    return obj


def maybe_iterable_to_list(obj: Union[Iterable[T], T]) -> Union[Collection[T], T]:
    """
    If obj is Iterable but not list-like, consume into list.
    """
    if isinstance(obj, abc.Iterable) and not isinstance(obj, abc.Sized):
        return list(obj)
    return obj


def is_null_slice(obj):
    """
    We have a null slice.
    """
    return (
        isinstance(obj, slice)
        and obj.start is None
        and obj.stop is None
        and obj.step is None
    )


def is_true_slices(l):
    """
    Find non-trivial slices in "l": return a list of booleans with same length.
    """
    return [isinstance(k, slice) and not is_null_slice(k) for k in l]


# TODO: used only once in indexing; belongs elsewhere?
def is_full_slice(obj, l):
    """
    We have a full length slice.
    """
    return (
        isinstance(obj, slice) and obj.start == 0 and obj.stop == l and obj.step is None
    )


def get_callable_name(obj):
    # typical case has name
    if hasattr(obj, "__name__"):
        return getattr(obj, "__name__")
    # some objects don't; could recurse
    if isinstance(obj, partial):
        return get_callable_name(obj.func)
    # fall back to class name
    if hasattr(obj, "__call__"):
        return type(obj).__name__
    # everything failed (probably because the argument
    # wasn't actually callable); we return None
    # instead of the empty string in this case to allow
    # distinguishing between no name and a name of ''
    return None


def apply_if_callable(maybe_callable, obj, **kwargs):
    """
    Evaluate possibly callable input using obj and kwargs if it is callable,
    otherwise return as it is.

    Parameters
    ----------
    maybe_callable : possibly a callable
    obj : NDFrame
    **kwargs
    """

    if callable(maybe_callable):
        return maybe_callable(obj, **kwargs)

    return maybe_callable


def dict_compat(d):
    """
    Helper function to convert datetimelike-keyed dicts
    to Timestamp-keyed dict.

    Parameters
    ----------
    d: dict like object

    Returns
    -------
    dict

    """
    return {maybe_box_datetimelike(key): value for key, value in d.items()}


def standardize_mapping(into):
    """
    Helper function to standardize a supplied mapping.

    .. versionadded:: 0.21.0

    Parameters
    ----------
    into : instance or subclass of collections.abc.Mapping
        Must be a class, an initialized collections.defaultdict,
        or an instance of a collections.abc.Mapping subclass.

    Returns
    -------
    mapping : a collections.abc.Mapping subclass or other constructor
        a callable object that can accept an iterator to create
        the desired Mapping.

    See Also
    --------
    DataFrame.to_dict
    Series.to_dict
    """
    if not inspect.isclass(into):
        if isinstance(into, collections.defaultdict):
            return partial(collections.defaultdict, into.default_factory)
        into = type(into)
    if not issubclass(into, abc.Mapping):
        raise TypeError(f"unsupported type: {into}")
    elif into == collections.defaultdict:
        raise TypeError("to_dict() only accepts initialized defaultdicts")
    return into


def random_state(state=None):
    """
    Helper function for processing random_state arguments.

    Parameters
    ----------
    state : int, np.random.RandomState, None.
        If receives an int, passes to np.random.RandomState() as seed.
        If receives an np.random.RandomState object, just returns object.
        If receives `None`, returns np.random.
        If receives anything else, raises an informative ValueError.
        Default None.

    Returns
    -------
    np.random.RandomState
    """

    if is_integer(state):
        return np.random.RandomState(state)
    elif isinstance(state, np.random.RandomState):
        return state
    elif state is None:
        return np.random
    else:
        raise ValueError(
            "random_state must be an integer, a numpy RandomState, or None"
        )


def pipe(obj, func, *args, **kwargs):
    """
    Apply a function ``func`` to object ``obj`` either by passing obj as the
    first argument to the function or, in the case that the func is a tuple,
    interpret the first element of the tuple as a function and pass the obj to
    that function as a keyword argument whose key is the value of the second
    element of the tuple.

    Parameters
    ----------
    func : callable or tuple of (callable, str)
        Function to apply to this object or, alternatively, a
        ``(callable, data_keyword)`` tuple where ``data_keyword`` is a
        string indicating the keyword of `callable`` that expects the
        object.
    *args : iterable, optional
        Positional arguments passed into ``func``.
    **kwargs : dict, optional
        A dictionary of keyword arguments passed into ``func``.

    Returns
    -------
    object : the return type of ``func``.
    """
    if isinstance(func, tuple):
        func, target = func
        if target in kwargs:
            msg = f"{target} is both the pipe target and a keyword argument"
            raise ValueError(msg)
        kwargs[target] = obj
        return func(*args, **kwargs)
    else:
        return func(obj, *args, **kwargs)


def get_rename_function(mapper):
    """
    Returns a function that will map names/labels, dependent if mapper
    is a dict, Series or just a function.
    """
    if isinstance(mapper, (abc.Mapping, ABCSeries)):

        def f(x):
            if x in mapper:
                return mapper[x]
            else:
                return x

    else:
        f = mapper

    return f
