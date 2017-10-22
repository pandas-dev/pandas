"""
Misc tools for implementing data structures
"""

import sys
import warnings
from datetime import datetime, timedelta
from functools import partial
import inspect
import collections

import numpy as np
from pandas._libs import lib, tslib

from pandas import compat
from pandas.compat import long, zip, iteritems
from pandas.core.config import get_option
from pandas.core.dtypes.generic import ABCSeries, ABCIndex
from pandas.core.dtypes.common import _NS_DTYPE
from pandas.core.dtypes.inference import _iterable_not_string
from pandas.core.dtypes.missing import isna, isnull, notnull  # noqa
from pandas.api import types
from pandas.core.dtypes import common

# compat
from pandas.errors import (  # noqa
    PerformanceWarning, UnsupportedFunctionCall, UnsortedIndexError)

# back-compat of public API
# deprecate these functions
m = sys.modules['pandas.core.common']
for t in [t for t in dir(types) if not t.startswith('_')]:

    def outer(t=t):

        def wrapper(*args, **kwargs):
            warnings.warn("pandas.core.common.{t} is deprecated. "
                          "import from the public API: "
                          "pandas.api.types.{t} instead".format(t=t),
                          DeprecationWarning, stacklevel=3)
            return getattr(types, t)(*args, **kwargs)
        return wrapper

    setattr(m, t, outer(t))

# back-compat for non-public functions
# deprecate these functions
for t in ['is_datetime_arraylike',
          'is_datetime_or_timedelta_dtype',
          'is_datetimelike',
          'is_datetimelike_v_numeric',
          'is_datetimelike_v_object',
          'is_datetimetz',
          'is_int_or_datetime_dtype',
          'is_period_arraylike',
          'is_string_like',
          'is_string_like_dtype']:

    def outer(t=t):

        def wrapper(*args, **kwargs):
            warnings.warn("pandas.core.common.{t} is deprecated. "
                          "These are not longer public API functions, "
                          "but can be imported from "
                          "pandas.api.types.{t} instead".format(t=t),
                          DeprecationWarning, stacklevel=3)
            return getattr(common, t)(*args, **kwargs)
        return wrapper

    setattr(m, t, outer(t))


# deprecate array_equivalent

def array_equivalent(*args, **kwargs):
    warnings.warn("'pandas.core.common.array_equivalent' is deprecated and "
                  "is no longer public API", DeprecationWarning, stacklevel=2)
    from pandas.core.dtypes import missing
    return missing.array_equivalent(*args, **kwargs)


class SettingWithCopyError(ValueError):
    pass


class SettingWithCopyWarning(Warning):
    pass


class AbstractMethodError(NotImplementedError):
    """Raise this error instead of NotImplementedError for abstract methods
    while keeping compatibility with Python 2 and Python 3.
    """

    def __init__(self, class_instance):
        self.class_instance = class_instance

    def __str__(self):
        msg = "This method must be defined in the concrete class of {name}"
        return (msg.format(name=self.class_instance.__class__.__name__))


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


def _consensus_name_attr(objs):
    name = objs[0].name
    for obj in objs[1:]:
        if obj.name != name:
            return None
    return name


def _maybe_match_name(a, b):
    a_has = hasattr(a, 'name')
    b_has = hasattr(b, 'name')
    if a_has and b_has:
        if a.name == b.name:
            return a.name
        else:
            return None
    elif a_has:
        return a.name
    elif b_has:
        return b.name
    return None


def _get_info_slice(obj, indexer):
    """Slice the info axis of `obj` with `indexer`."""
    if not hasattr(obj, '_info_axis_number'):
        msg = 'object of type {typ!r} has no info axis'
        raise TypeError(msg.format(typ=type(obj).__name__))
    slices = [slice(None)] * obj.ndim
    slices[obj._info_axis_number] = indexer
    return tuple(slices)


def _maybe_box(indexer, values, obj, key):

    # if we have multiples coming back, box em
    if isinstance(values, np.ndarray):
        return obj[indexer.get_loc(key)]

    # return the value
    return values


def _maybe_box_datetimelike(value):
    # turn a datetime like into a Timestamp/timedelta as needed

    if isinstance(value, (np.datetime64, datetime)):
        value = tslib.Timestamp(value)
    elif isinstance(value, (np.timedelta64, timedelta)):
        value = tslib.Timedelta(value)

    return value


_values_from_object = lib.values_from_object


def is_bool_indexer(key):
    if isinstance(key, (ABCSeries, np.ndarray, ABCIndex)):
        if key.dtype == np.object_:
            key = np.asarray(_values_from_object(key))

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


def _default_index(n):
    from pandas.core.index import RangeIndex
    return RangeIndex(0, n, name=None)


def _mut_exclusive(**kwargs):
    item1, item2 = kwargs.items()
    label1, val1 = item1
    label2, val2 = item2
    if val1 is not None and val2 is not None:
        msg = 'mutually exclusive arguments: {label1!r} and {label2!r}'
        raise TypeError(msg.format(label1=label1, label2=label2))
    elif val1 is not None:
        return val1
    else:
        return val2


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


def _count_not_none(*args):
    """Returns the count of arguments that are not None"""
    return sum(x is not None for x in args)


def _try_sort(iterable):
    listed = list(iterable)
    try:
        return sorted(listed)
    except Exception:
        return listed


def iterpairs(seq):
    """
    Parameters
    ----------
    seq : sequence

    Returns
    -------
    iterator returning overlapping pairs of elements

    Examples
    --------
    >>> list(iterpairs([1, 2, 3, 4]))
    [(1, 2), (2, 3), (3, 4)]
    """
    # input may not be sliceable
    seq_it = iter(seq)
    seq_it_next = iter(seq)
    next(seq_it_next)

    return zip(seq_it, seq_it_next)


def split_ranges(mask):
    """ Generates tuples of ranges which cover all True value in mask

    >>> list(split_ranges([1,0,0,1,0]))
    [(0, 1), (3, 4)]
    """
    ranges = [(0, len(mask))]

    for pos, val in enumerate(mask):
        if not val:  # this pos should be ommited, split off the prefix range
            r = ranges.pop()
            if pos > r[0]:  # yield non-zero range
                yield (r[0], pos)
            if pos + 1 < len(mask):  # save the rest for processing
                ranges.append((pos + 1, len(mask)))
    if ranges:
        yield ranges[-1]


def _long_prod(vals):
    result = long(1)
    for x in vals:
        result *= x
    return result


class groupby(dict):
    """
    A simple groupby different from the one in itertools.

    Does not require the sequence elements to be sorted by keys,
    however it is slower.
    """

    def __init__(self, seq, key=lambda x: x):
        for value in seq:
            k = key(value)
            self.setdefault(k, []).append(value)

    try:
        __iter__ = dict.iteritems
    except AttributeError:  # pragma: no cover
        # Python 3
        def __iter__(self):
            return iter(dict.items(self))


def map_indices_py(arr):
    """
    Returns a dictionary with (element, index) pairs for each element in the
    given array/list
    """
    return dict([(x, i) for i, x in enumerate(arr)])


def union(*seqs):
    result = set([])
    for seq in seqs:
        if not isinstance(seq, set):
            seq = set(seq)
        result |= seq
    return type(seqs[0])(list(result))


def difference(a, b):
    return type(a)(list(set(a) - set(b)))


def intersection(*seqs):
    result = set(seqs[0])
    for seq in seqs:
        if not isinstance(seq, set):
            seq = set(seq)
        result &= seq
    return type(seqs[0])(list(result))


def _asarray_tuplesafe(values, dtype=None):
    from pandas.core.index import Index

    if not (isinstance(values, (list, tuple)) or hasattr(values, '__array__')):
        values = list(values)
    elif isinstance(values, Index):
        return values.values

    if isinstance(values, list) and dtype in [np.object_, object]:
        return lib.list_to_object_array(values)

    result = np.asarray(values, dtype=dtype)

    if issubclass(result.dtype.type, compat.string_types):
        result = np.asarray(values, dtype=object)

    if result.ndim == 2:
        if isinstance(values, list):
            return lib.list_to_object_array(values)
        else:
            # Making a 1D array that safely contains tuples is a bit tricky
            # in numpy, leading to the following
            try:
                result = np.empty(len(values), dtype=object)
                result[:] = values
            except ValueError:
                # we have a list-of-list
                result[:] = [tuple(x) for x in values]

    return result


def _index_labels_to_array(labels):
    if isinstance(labels, (compat.string_types, tuple)):
        labels = [labels]

    if not isinstance(labels, (list, np.ndarray)):
        try:
            labels = list(labels)
        except TypeError:  # non-iterable
            labels = [labels]

    labels = _asarray_tuplesafe(labels)

    return labels


def _maybe_make_list(obj):
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


def is_full_slice(obj, l):
    """ we have a full length slice """
    return (isinstance(obj, slice) and obj.start == 0 and obj.stop == l and
            obj.step is None)


def _get_callable_name(obj):
    # typical case has name
    if hasattr(obj, '__name__'):
        return getattr(obj, '__name__')
    # some objects don't; could recurse
    if isinstance(obj, partial):
        return _get_callable_name(obj.func)
    # fall back to class name
    if hasattr(obj, '__call__'):
        return obj.__class__.__name__
    # everything failed (probably because the argument
    # wasn't actually callable); we return None
    # instead of the empty string in this case to allow
    # distinguishing between no name and a name of ''
    return None


def _apply_if_callable(maybe_callable, obj, **kwargs):
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


def _where_compat(mask, arr1, arr2):
    if arr1.dtype == _NS_DTYPE and arr2.dtype == _NS_DTYPE:
        new_vals = np.where(mask, arr1.view('i8'), arr2.view('i8'))
        return new_vals.view(_NS_DTYPE)

    if arr1.dtype == _NS_DTYPE:
        arr1 = tslib.ints_to_pydatetime(arr1.view('i8'))
    if arr2.dtype == _NS_DTYPE:
        arr2 = tslib.ints_to_pydatetime(arr2.view('i8'))

    return np.where(mask, arr1, arr2)


def _dict_compat(d):
    """
    Helper function to convert datetimelike-keyed dicts to Timestamp-keyed dict

    Parameters
    ----------
    d: dict like object

    Returns
    -------
    dict

    """
    return dict((_maybe_box_datetimelike(key), value)
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


# ----------------------------------------------------------------------
# Detect our environment

def in_interactive_session():
    """ check if we're running in an interactive shell

    returns True if running under python/ipython interactive shell
    """

    def check_main():
        import __main__ as main
        return (not hasattr(main, '__file__') or
                get_option('mode.sim_interactive'))

    try:
        return __IPYTHON__ or check_main()  # noqa
    except:
        return check_main()


def in_qtconsole():
    """
    check if we're inside an IPython qtconsole

    .. deprecated:: 0.14.1
       This is no longer needed, or working, in IPython 3 and above.
    """
    try:
        ip = get_ipython()  # noqa
        front_end = (
            ip.config.get('KernelApp', {}).get('parent_appname', "") or
            ip.config.get('IPKernelApp', {}).get('parent_appname', ""))
        if 'qtconsole' in front_end.lower():
            return True
    except:
        return False
    return False


def in_ipnb():
    """
    check if we're inside an IPython Notebook

    .. deprecated:: 0.14.1
       This is no longer needed, or working, in IPython 3 and above.
    """
    try:
        ip = get_ipython()  # noqa
        front_end = (
            ip.config.get('KernelApp', {}).get('parent_appname', "") or
            ip.config.get('IPKernelApp', {}).get('parent_appname', ""))
        if 'notebook' in front_end.lower():
            return True
    except:
        return False
    return False


def in_ipython_frontend():
    """
    check if we're inside an an IPython zmq frontend
    """
    try:
        ip = get_ipython()  # noqa
        return 'zmq' in str(type(ip)).lower()
    except:
        pass

    return False


def _random_state(state=None):
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

    if types.is_integer(state):
        return np.random.RandomState(state)
    elif isinstance(state, np.random.RandomState):
        return state
    elif state is None:
        return np.random
    else:
        raise ValueError("random_state must be an integer, a numpy "
                         "RandomState, or None")


def _get_distinct_objs(objs):
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
