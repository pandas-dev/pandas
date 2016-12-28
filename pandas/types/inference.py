""" basic inference routines """

import collections
import re
import numpy as np
from numbers import Number
from pandas.compat import (string_types, text_type,
                           string_and_binary_types)
from pandas import lib

is_bool = lib.is_bool

is_integer = lib.is_integer

is_float = lib.is_float

is_complex = lib.is_complex

is_scalar = lib.isscalar

is_decimal = lib.is_decimal


def is_number(obj):
    return isinstance(obj, (Number, np.number))


def is_string_like(obj):
    return isinstance(obj, (text_type, string_types))


def _iterable_not_string(x):
    return (isinstance(x, collections.Iterable) and
            not isinstance(x, string_types))


def is_iterator(obj):
    # python 3 generators have __next__ instead of next
    return hasattr(obj, 'next') or hasattr(obj, '__next__')


def is_re(obj):
    return isinstance(obj, re._pattern_type)


def is_re_compilable(obj):
    try:
        re.compile(obj)
    except TypeError:
        return False
    else:
        return True


def is_list_like(arg):
    return (hasattr(arg, '__iter__') and
            not isinstance(arg, string_and_binary_types))


def is_dict_like(arg):
    return hasattr(arg, '__getitem__') and hasattr(arg, 'keys')


def is_named_tuple(arg):
    return isinstance(arg, tuple) and hasattr(arg, '_fields')


def is_hashable(arg):
    """Return True if hash(arg) will succeed, False otherwise.

    Some types will pass a test against collections.Hashable but fail when they
    are actually hashed with hash().

    Distinguish between these and other types by trying the call to hash() and
    seeing if they raise TypeError.

    Examples
    --------
    >>> a = ([],)
    >>> isinstance(a, collections.Hashable)
    True
    >>> is_hashable(a)
    False
    """
    # unfortunately, we can't use isinstance(arg, collections.Hashable), which
    # can be faster than calling hash, because numpy scalars on Python 3 fail
    # this test

    # reconsider this decision once this numpy bug is fixed:
    # https://github.com/numpy/numpy/issues/5562

    try:
        hash(arg)
    except TypeError:
        return False
    else:
        return True


def is_sequence(x):
    try:
        iter(x)
        len(x)  # it has a length
        return not isinstance(x, string_and_binary_types)
    except (TypeError, AttributeError):
        return False
