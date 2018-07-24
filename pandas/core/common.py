"""
Misc tools for implementing data structures

Note: pandas.core.common is *not* part of the public API.
"""

from datetime import datetime, timedelta
from functools import partial
import inspect
import collections

import numpy as np
from pandas._libs import lib, tslibs

from pandas import compat
from pandas.compat import iteritems, PY36, OrderedDict
from pandas.core.dtypes.generic import ABCSeries, ABCIndex, ABCIndexClass
from pandas.core.dtypes.common import is_integer
from pandas.core.dtypes.inference import _iterable_not_string
from pandas.core.dtypes.missing import isna, isnull, notnull  # noqa
from pandas.core.dtypes.cast import construct_1d_object_array_from_listlike


class SettingWithCopyError(ValueError):
    pass


class SettingWithCopyWarning(Warning):
    pass


def flatten(l):
    """Flatten an arbitrarily nested sequence.

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


# TODO: only used once in frame.py; belongs elsewhere?
def get_info_slice(obj, indexer):
    """Slice the info axis of `obj` with `indexer`."""
    if not hasattr(obj, '_info_axis_number'):
        msg = 'object of type {typ!r} has no info axis'
        raise TypeError(msg.format(typ=type(obj).__name__))
    slices = [slice(None)] * obj.ndim
    slices[obj._info_axis_number] = indexer
    return tuple(slices)


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


def is_bool_indexer(key):
    if isinstance(key, (ABCSeries, np.ndarray, ABCIndex)):
        if key.dtype == np.object_:
            key = np.asarray(values_from_object(key))

            if not lib.is_bool_array(key):
                if isna(key).any():
                    raise ValueError('cannot index with vector containing '
                                     'NA / NaN values')
                return False
            return True
        elif key.dtype == np.bool_:
            return True
    elif isinstance(key, list):
        try:
            arr = np.asarray(key)
            return arr.dtype == np.bool_ and len(arr) == len(key)
        except TypeError:  # pragma: no cover
            return False

    return False


def _not_none(*args):
    """Returns a generator consisting of the arguments that are not None"""
    return (arg for arg in args if arg is not None)


def _any_none(*args):
    """Returns a boolean indicating if any argument is None"""
    for arg in args:
        if arg is None:
            return True
    return False


def _all_none(*args):
    """Returns a boolean indicating if all arguments are None"""
    for arg in args:
        if arg is not None:
            return False
    return True


def _any_not_none(*args):
    """Returns a boolean indicating if any argument is not None"""
    for arg in args:
        if arg is not None:
            return True
    return False


def _all_not_none(*args):
    """Returns a boolean indicating if all arguments are not None"""
    for arg in args:
        if arg is None:
            return False
    return True


def count_not_none(*args):
    """Returns the count of arguments that are not None"""
    return sum(x is not None for x in args)


def try_sort(iterable):
    listed = list(iterable)
    try:
        return sorted(listed)
    except Exception:
        return listed


def dict_keys_to_ordered_list(mapping):
    # when pandas drops support for Python < 3.6, this function
    # can be replaced by a simple list(mapping.keys())
    if PY36 or isinstance(mapping, OrderedDict):
        keys = list(mapping.keys())
    else:
        keys = try_sort(mapping)
    return keys


def asarray_tuplesafe(values, dtype=None):

    if not (isinstance(values, (list, tuple)) or hasattr(values, '__array__')):
        values = list(values)
    elif isinstance(values, ABCIndexClass):
        return values.values

    if isinstance(values, list) and dtype in [np.object_, object]:
        return construct_1d_object_array_from_listlike(values)

    result = np.asarray(values, dtype=dtype)

    if issubclass(result.dtype.type, compat.string_types):
        result = np.asarray(values, dtype=object)

    if result.ndim == 2:
        # Avoid building an array of arrays:
        # TODO: verify whether any path hits this except #18819 (invalid)
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
    if isinstance(labels, (compat.string_types, tuple)):
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


def is_null_slice(obj):
    """ we have a null slice """
    return (isinstance(obj, slice) and obj.start is None and
            obj.stop is None and obj.step is None)


def is_true_slices(l):
    """
    Find non-trivial slices in "l": return a list of booleans with same length.
    """
    return [isinstance(k, slice) and not is_null_slice(k) for k in l]


# TODO: used only once in indexing; belongs elsewhere?
def is_full_slice(obj, l):
    """ we have a full length slice """
    return (isinstance(obj, slice) and obj.start == 0 and obj.stop == l and
            obj.step is None)


def get_callable_name(obj):
    # typical case has name
    if hasattr(obj, '__name__'):
        return getattr(obj, '__name__')
    # some objects don't; could recurse
    if isinstance(obj, partial):
        return get_callable_name(obj.func)
    # fall back to class name
    if hasattr(obj, '__call__'):
        return obj.__class__.__name__
    # everything failed (probably because the argument
    # wasn't actually callable); we return None
    # instead of the empty string in this case to allow
    # distinguishing between no name and a name of ''
    return None


def apply_if_callable(maybe_callable, obj, **kwargs):
    """
    Evaluate possibly callable input using obj and kwargs if it is callable,
    otherwise return as it is

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
    Helper function to convert datetimelike-keyed dicts to Timestamp-keyed dict

    Parameters
    ----------
    d: dict like object

    Returns
    -------
    dict

    """
    return dict((maybe_box_datetimelike(key), value)
                for key, value in iteritems(d))


def standardize_mapping(into):
    """
    Helper function to standardize a supplied mapping.

    .. versionadded:: 0.21.0

    Parameters
    ----------
    into : instance or subclass of collections.Mapping
        Must be a class, an initialized collections.defaultdict,
        or an instance of a collections.Mapping subclass.

    Returns
    -------
    mapping : a collections.Mapping subclass or other constructor
        a callable object that can accept an iterator to create
        the desired Mapping.

    See Also
    --------
    DataFrame.to_dict
    Series.to_dict
    """
    if not inspect.isclass(into):
        if isinstance(into, collections.defaultdict):
            return partial(
                collections.defaultdict, into.default_factory)
        into = type(into)
    if not issubclass(into, collections.Mapping):
        raise TypeError('unsupported type: {into}'.format(into=into))
    elif into == collections.defaultdict:
        raise TypeError(
            'to_dict() only accepts initialized defaultdicts')
    return into


def sentinel_factory():
    class Sentinel(object):
        pass

    return Sentinel()


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
        raise ValueError("random_state must be an integer, a numpy "
                         "RandomState, or None")


# TODO: only used once in indexes.api; belongs elsewhere?
def get_distinct_objs(objs):
    """
    Return a list with distinct elements of "objs" (different ids).
    Preserves order.
    """
    ids = set()
    res = []
    for obj in objs:
        if not id(obj) in ids:
            ids.add(id(obj))
            res.append(obj)
    return res


def _pipe(obj, func, *args, **kwargs):
    """
    Apply a function ``func`` to object ``obj`` either by passing obj as the
    first argument to the function or, in the case that the func is a tuple,
    interpret the first element of the tuple as a function and pass the obj to
    that function as a keyword argument whose key is the value of the second
    element of the tuple.

    Parameters
    ----------
    func : callable or tuple of (callable, string)
        Function to apply to this object or, alternatively, a
        ``(callable, data_keyword)`` tuple where ``data_keyword`` is a
        string indicating the keyword of `callable`` that expects the
        object.
    args : iterable, optional
        positional arguments passed into ``func``.
    kwargs : dict, optional
        a dictionary of keyword arguments passed into ``func``.

    Returns
    -------
    object : the return type of ``func``.
    """
    if isinstance(func, tuple):
        func, target = func
        if target in kwargs:
            msg = '%s is both the pipe target and a keyword argument' % target
            raise ValueError(msg)
        kwargs[target] = obj
        return func(*args, **kwargs)
    else:
        return func(obj, *args, **kwargs)
