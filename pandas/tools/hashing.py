"""
data hash pandas / numpy objects
"""

import numpy as np
from pandas import _hash, Series, factorize, Categorical, Index
from pandas.lib import infer_dtype
from pandas.types.generic import ABCIndexClass, ABCSeries, ABCDataFrame
from pandas.types.common import is_categorical_dtype

# 16 byte long hashing key
_default_hash_key = '0123456789123456'


def hash_pandas_object(obj, index=True, encoding='utf8', hash_key=None):
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

    Returns
    -------
    Series of uint64, same length as the object

    """
    if hash_key is None:
        hash_key = _default_hash_key

    def adder(h, hashed_to_add):
        h = np.multiply(h, np.uint(3), h)
        return np.add(h, hashed_to_add, h)

    if isinstance(obj, ABCIndexClass):
        h = hash_array(obj.values, encoding, hash_key).astype('uint64')
        h = Series(h, index=obj, dtype='uint64')
    elif isinstance(obj, ABCSeries):
        h = hash_array(obj.values, encoding, hash_key).astype('uint64')
        if index:
            h = adder(h, hash_pandas_object(obj.index,
                                            index=False,
                                            encoding=encoding,
                                            hash_key=hash_key).values)
        h = Series(h, index=obj.index, dtype='uint64')
    elif isinstance(obj, ABCDataFrame):
        cols = obj.iteritems()
        first_series = next(cols)[1]
        h = hash_array(first_series.values, encoding,
                       hash_key).astype('uint64')
        for _, col in cols:
            h = adder(h, hash_array(col.values, encoding, hash_key))
        if index:
            h = adder(h, hash_pandas_object(obj.index,
                                            index=False,
                                            encoding=encoding,
                                            hash_key=hash_key).values)

        h = Series(h, index=obj.index, dtype='uint64')
    else:
        raise TypeError("Unexpected type for hashing %s" % type(obj))
    return h


def hash_array(vals, encoding='utf8', hash_key=None):
    """
    Given a 1d array, return an array of deterministic integers.

    .. versionadded:: 0.19.2

    Parameters
    ----------
    vals : ndarray
    encoding : string, default 'utf8'
        encoding for data & key when strings
    hash_key : string key to encode, default to _default_hash_key

    Returns
    -------
    1d uint64 numpy array of hash values, same length as the vals

    """

    # work with cagegoricals as ints. (This check is above the complex
    # check so that we don't ask numpy if categorical is a subdtype of
    # complex, as it will choke.
    if hash_key is None:
        hash_key = _default_hash_key

    if is_categorical_dtype(vals.dtype):
        vals = vals.codes

    # we'll be working with everything as 64-bit values, so handle this
    # 128-bit value early
    if np.issubdtype(vals.dtype, np.complex128):
        return hash_array(vals.real) + 23 * hash_array(vals.imag)

    # MAIN LOGIC:
    inferred = infer_dtype(vals)

    # First, turn whatever array this is into unsigned 64-bit ints, if we can
    # manage it.
    if inferred == 'boolean':
        vals = vals.astype('u8')

    if (np.issubdtype(vals.dtype, np.datetime64) or
       np.issubdtype(vals.dtype, np.timedelta64) or
       np.issubdtype(vals.dtype, np.number)) and vals.dtype.itemsize <= 8:

        vals = vals.view('u{}'.format(vals.dtype.itemsize)).astype('u8')
    else:

        # its MUCH faster to categorize object dtypes, then hash and rename
        codes, categories = factorize(vals, sort=False)
        categories = Index(categories)
        c = Series(Categorical(codes, categories,
                               ordered=False, fastpath=True))
        vals = _hash.hash_object_array(categories.values,
                                       hash_key,
                                       encoding)

        # rename & extract
        vals = c.cat.rename_categories(Index(vals)).astype(np.uint64).values

    # Then, redistribute these 64-bit ints within the space of 64-bit ints
    vals ^= vals >> 30
    vals *= np.uint64(0xbf58476d1ce4e5b9)
    vals ^= vals >> 27
    vals *= np.uint64(0x94d049bb133111eb)
    vals ^= vals >> 31
    return vals
