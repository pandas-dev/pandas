"""
data hash pandas / numpy objects
"""

import numpy as np
from pandas import _hash, Series, factorize, Categorical, Index, MultiIndex
from pandas.lib import is_bool_array
from pandas.types.generic import ABCIndexClass, ABCSeries, ABCDataFrame
from pandas.types.common import (is_categorical_dtype, is_numeric_dtype,
                                 is_datetime64_dtype, is_timedelta64_dtype,
                                 is_object_dtype)

# 16 byte long hashing key
_default_hash_key = '0123456789123456'


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

    def adder(h, hashed_to_add):
        h = np.multiply(h, np.uint(3), h)
        return np.add(h, hashed_to_add, h)

    if isinstance(obj, MultiIndex):
        return _hash_tuples(obj, encoding, hash_key)

    if isinstance(obj, ABCIndexClass):
        h = hash_array(obj.values, encoding, hash_key,
                       categorize).astype('uint64')
        h = Series(h, index=obj, dtype='uint64')
    elif isinstance(obj, ABCSeries):
        h = hash_array(obj.values, encoding, hash_key,
                       categorize).astype('uint64')
        if index:
            h = adder(h, hash_pandas_object(obj.index,
                                            index=False,
                                            encoding=encoding,
                                            hash_key=hash_key,
                                            categorize=categorize).values)
        h = Series(h, index=obj.index, dtype='uint64')
    elif isinstance(obj, ABCDataFrame):
        cols = obj.iteritems()
        first_series = next(cols)[1]
        h = hash_array(first_series.values, encoding,
                       hash_key, categorize).astype('uint64')
        for _, col in cols:
            h = adder(h, hash_array(col.values, encoding, hash_key,
                                    categorize))
        if index:
            h = adder(h, hash_pandas_object(obj.index,
                                            index=False,
                                            encoding=encoding,
                                            hash_key=hash_key,
                                            categorize=categorize).values)

        h = Series(h, index=obj.index, dtype='uint64')
    else:
        raise TypeError("Unexpected type for hashing %s" % type(obj))
    return h


def _hash_tuples(vals, encoding, hash_key):
    """
    Hash an MultiIndex / array_of_tuples efficiently

    Parameters
    ----------
    vals : MultiIndex or ndarray of tuples
    encoding : string, default 'utf8'
    hash_key : string key to encode, default to _default_hash_key

    Returns
    -------
    ndarray of hashed values array, same size as len(c)
    """

    if not isinstance(vals, MultiIndex):
        vals = MultiIndex.from_tuples(vals)

    # efficiently turn us into a DataFrame and hash
    return hash_pandas_object(vals.to_dataframe(index=False),
                              index=False, encoding=encoding,
                              hash_key=hash_key, categorize=False)


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
    cat_hashed = hash_array(c.categories.values, encoding, hash_key,
                            categorize=False).astype(np.uint64, copy=False)
    return c.rename_categories(cat_hashed).astype(np.uint64)


def hash_array(vals, encoding='utf8', hash_key=None, categorize=True):
    """
    Given a 1d array, return an array of deterministic integers.

    .. versionadded:: 0.19.2

    Parameters
    ----------
    vals : ndarray
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

    if hash_key is None:
        hash_key = _default_hash_key

    if isinstance(vals, list) and len(vals) and isinstance(vals[0], tuple):
        # we hash an list of tuples similar to a MultiIndex
        return _hash_tuples(vals, encoding, hash_key).values

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
    elif ((is_datetime64_dtype(vals) or
           is_timedelta64_dtype(vals) or
           is_numeric_dtype(vals)) and vals.dtype.itemsize <= 8):
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

        vals = _hash.hash_object_array(vals, hash_key, encoding)

    # Then, redistribute these 64-bit ints within the space of 64-bit ints
    vals ^= vals >> 30
    vals *= np.uint64(0xbf58476d1ce4e5b9)
    vals ^= vals >> 27
    vals *= np.uint64(0x94d049bb133111eb)
    vals ^= vals >> 31
    return vals
