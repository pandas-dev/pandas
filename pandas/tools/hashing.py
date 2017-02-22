"""
data hash pandas / numpy objects
"""
import itertools

import numpy as np
from pandas import _hash, Series, factorize, Categorical, Index, MultiIndex
from pandas.lib import is_bool_array
from pandas.types.generic import ABCIndexClass, ABCSeries, ABCDataFrame
from pandas.types.common import (is_categorical_dtype, is_numeric_dtype,
                                 is_datetime64_dtype, is_timedelta64_dtype,
                                 is_list_like)

# 16 byte long hashing key
_default_hash_key = '0123456789123456'


def _combine_hash_arrays(arrays, num_items):
    """
    Parameters
    ----------
    arrays : generator
    num_items : int

    Should be the same as CPython's tupleobject.c
    """
    try:
        first = next(arrays)
    except StopIteration:
        return np.array([], dtype=np.uint64)

    arrays = itertools.chain([first], arrays)

    mult = np.uint64(1000003)
    out = np.zeros_like(first) + np.uint64(0x345678)
    for i, a in enumerate(arrays):
        inverse_i = num_items - i
        out ^= a
        out *= mult
        mult += np.uint64(82520 + inverse_i + inverse_i)
    assert i + 1 == num_items, 'Fed in wrong num_items'
    out += np.uint64(97531)
    return out


def hash_pandas_object(obj, index=True, encoding='utf8', hash_key=None,
                       categorize=True):
    """
    Return a data hash of the Index/Series/DataFrame

    .. versionadded:: 0.19.2

    Parameters
    ----------
    index : boolean, default True
        include the index in the hash (if Series/DataFrame)
    encoding : string, default 'utf8'
        encoding for data & key when strings
    hash_key : string key to encode, default to _default_hash_key
    categorize : bool, default True
        Whether to first categorize object arrays before hashing. This is more
        efficient when the array contains duplicate values.

        .. versionadded:: 0.20.0

    Returns
    -------
    Series of uint64, same length as the object

    """
    if hash_key is None:
        hash_key = _default_hash_key

    if isinstance(obj, MultiIndex):
        return Series(hash_tuples(obj, encoding, hash_key),
                      dtype='uint64', copy=False)

    if isinstance(obj, ABCIndexClass):
        h = hash_array(obj.values, encoding, hash_key,
                       categorize).astype('uint64', copy=False)
        h = Series(h, index=obj, dtype='uint64', copy=False)
    elif isinstance(obj, ABCSeries):
        h = hash_array(obj.values, encoding, hash_key,
                       categorize).astype('uint64', copy=False)
        if index:
            index_iter = (hash_pandas_object(obj.index,
                                             index=False,
                                             encoding=encoding,
                                             hash_key=hash_key,
                                             categorize=categorize).values
                          for _ in [None])
            arrays = itertools.chain([h], index_iter)
            h = _combine_hash_arrays(arrays, 2)

        h = Series(h, index=obj.index, dtype='uint64', copy=False)

    elif isinstance(obj, ABCDataFrame):
        hashes = (hash_array(series.values) for _, series in obj.iteritems())
        num_items = len(obj.columns)
        if index:
            index_hash_generator = (hash_pandas_object(obj.index,
                                                       index=False,
                                                       encoding=encoding,
                                                       hash_key=hash_key,
                                                       categorize=categorize).values  # noqa
                                    for _ in [None])
            num_items += 1
            hashes = itertools.chain(hashes, index_hash_generator)
        h = _combine_hash_arrays(hashes, num_items)

        h = Series(h, index=obj.index, dtype='uint64', copy=False)
    else:
        raise TypeError("Unexpected type for hashing %s" % type(obj))
    return h


def hash_tuples(vals, encoding='utf8', hash_key=None):
    """
    Hash an MultiIndex / list-of-tuples efficiently

    .. versionadded:: 0.20.0

    Parameters
    ----------
    vals : MultiIndex, list-of-tuples, or single tuple
    encoding : string, default 'utf8'
    hash_key : string key to encode, default to _default_hash_key

    Returns
    -------
    ndarray of hashed values array
    """

    is_tuple = False
    if isinstance(vals, tuple):
        vals = [vals]
        is_tuple = True
    elif not is_list_like(vals):
        raise TypeError("must be convertible to a list-of-tuples")

    if not isinstance(vals, MultiIndex):
        vals = MultiIndex.from_tuples(vals)

    # create a list-of-Categoricals
    vals = [Categorical(vals.labels[level],
                        vals.levels[level],
                        ordered=False,
                        fastpath=True)
            for level in range(vals.nlevels)]

    # hash the list-of-ndarrays
    hashes = (_hash_categorical(cat,
                                encoding=encoding,
                                hash_key=hash_key)
              for cat in vals)
    h = _combine_hash_arrays(hashes, len(vals))
    if is_tuple:
        h = h[0]

    return h


def _hash_categorical(c, encoding, hash_key):
    """
    Hash a Categorical by hashing its categories, and then mapping the codes
    to the hashes

    Parameters
    ----------
    c : Categorical
    encoding : string, default 'utf8'
    hash_key : string key to encode, default to _default_hash_key

    Returns
    -------
    ndarray of hashed values array, same size as len(c)
    """
    hashed = hash_array(c.categories.values, encoding, hash_key,
                        categorize=False)

    # we have uint64, as we don't directly support missing values
    # we don't want to use take_nd which will coerce to float
    # instead, directly construt the result with a
    # max(np.uint64) as the missing value indicator
    #
    # TODO: GH 15362

    mask = c.isnull()
    if len(hashed):
        result = hashed.take(c.codes)
    else:
        result = np.zeros(len(mask), dtype='uint64')

    if mask.any():
        result[mask] = np.iinfo(np.uint64).max

    return result


def hash_array(vals, encoding='utf8', hash_key=None, categorize=True):
    """
    Given a 1d array, return an array of deterministic integers.

    .. versionadded:: 0.19.2

    Parameters
    ----------
    vals : ndarray, Categorical
    encoding : string, default 'utf8'
        encoding for data & key when strings
    hash_key : string key to encode, default to _default_hash_key
    categorize : bool, default True
        Whether to first categorize object arrays before hashing. This is more
        efficient when the array contains duplicate values.

        .. versionadded:: 0.20.0

    Returns
    -------
    1d uint64 numpy array of hash values, same length as the vals

    """

    if not hasattr(vals, 'dtype'):
        raise TypeError("must pass a ndarray-like")

    if hash_key is None:
        hash_key = _default_hash_key

    # For categoricals, we hash the categories, then remap the codes to the
    # hash values. (This check is above the complex check so that we don't ask
    # numpy if categorical is a subdtype of complex, as it will choke.
    if is_categorical_dtype(vals.dtype):
        return _hash_categorical(vals, encoding, hash_key)

    # we'll be working with everything as 64-bit values, so handle this
    # 128-bit value early
    if np.issubdtype(vals.dtype, np.complex128):
        return hash_array(vals.real) + 23 * hash_array(vals.imag)

    # First, turn whatever array this is into unsigned 64-bit ints, if we can
    # manage it.
    if is_bool_array(vals):
        vals = vals.astype('u8')
    elif (is_datetime64_dtype(vals) or
          is_timedelta64_dtype(vals)):
        vals = vals.view('i8').astype('u8', copy=False)
    elif (is_numeric_dtype(vals) and vals.dtype.itemsize <= 8):
        vals = vals.view('u{}'.format(vals.dtype.itemsize)).astype('u8')
    else:
        # With repeated values, its MUCH faster to categorize object dtypes,
        # then hash and rename categories. We allow skipping the categorization
        # when the values are known/likely to be unique.
        if categorize:
            codes, categories = factorize(vals, sort=False)
            cat = Categorical(codes, Index(categories),
                              ordered=False, fastpath=True)
            return _hash_categorical(cat, encoding, hash_key)

        try:
            vals = _hash.hash_object_array(vals, hash_key, encoding)
        except TypeError:
            # we have mixed types
            vals = _hash.hash_object_array(vals.astype(str).astype(object),
                                           hash_key, encoding)

    # Then, redistribute these 64-bit ints within the space of 64-bit ints
    vals ^= vals >> 30
    vals *= np.uint64(0xbf58476d1ce4e5b9)
    vals ^= vals >> 27
    vals *= np.uint64(0x94d049bb133111eb)
    vals ^= vals >> 31
    return vals
