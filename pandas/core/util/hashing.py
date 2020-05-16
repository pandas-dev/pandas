"""
data hash pandas / numpy objects
"""
import itertools
from typing import Optional

import numpy as np

from pandas._libs import Timestamp
import pandas._libs.hashing as hashing

from pandas.core.dtypes.cast import infer_dtype_from_scalar
from pandas.core.dtypes.common import (
    is_categorical_dtype,
    is_extension_array_dtype,
    is_list_like,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCIndexClass,
    ABCMultiIndex,
    ABCSeries,
)
from pandas.core.dtypes.missing import isna

# 16 byte long hashing key
_default_hash_key = "0123456789123456"


def _combine_hash_arrays(arrays, num_items: int):
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
    assert i + 1 == num_items, "Fed in wrong num_items"
    out += np.uint64(97531)
    return out


def hash_pandas_object(
    obj,
    index: bool = True,
    encoding: str = "utf8",
    hash_key: Optional[str] = _default_hash_key,
    categorize: bool = True,
):
    """
    Return a data hash of the Index/Series/DataFrame.

    Parameters
    ----------
    index : bool, default True
        Include the index in the hash (if Series/DataFrame).
    encoding : str, default 'utf8'
        Encoding for data & key when strings.
    hash_key : str, default _default_hash_key
        Hash_key for string key to encode.
    categorize : bool, default True
        Whether to first categorize object arrays before hashing. This is more
        efficient when the array contains duplicate values.

    Returns
    -------
    Series of uint64, same length as the object
    """
    from pandas import Series

    if hash_key is None:
        hash_key = _default_hash_key

    if isinstance(obj, ABCMultiIndex):
        return Series(hash_tuples(obj, encoding, hash_key), dtype="uint64", copy=False)

    elif isinstance(obj, ABCIndexClass):
        h = hash_array(obj.values, encoding, hash_key, categorize).astype(
            "uint64", copy=False
        )
        h = Series(h, index=obj, dtype="uint64", copy=False)

    elif isinstance(obj, ABCSeries):
        h = hash_array(obj.values, encoding, hash_key, categorize).astype(
            "uint64", copy=False
        )
        if index:
            index_iter = (
                hash_pandas_object(
                    obj.index,
                    index=False,
                    encoding=encoding,
                    hash_key=hash_key,
                    categorize=categorize,
                ).values
                for _ in [None]
            )
            arrays = itertools.chain([h], index_iter)
            h = _combine_hash_arrays(arrays, 2)

        h = Series(h, index=obj.index, dtype="uint64", copy=False)

    elif isinstance(obj, ABCDataFrame):
        hashes = (hash_array(series.values) for _, series in obj.items())
        num_items = len(obj.columns)
        if index:
            index_hash_generator = (
                hash_pandas_object(
                    obj.index,
                    index=False,
                    encoding=encoding,
                    hash_key=hash_key,
                    categorize=categorize,
                ).values  # noqa
                for _ in [None]
            )
            num_items += 1

            # keep `hashes` specifically a generator to keep mypy happy
            _hashes = itertools.chain(hashes, index_hash_generator)
            hashes = (x for x in _hashes)
        h = _combine_hash_arrays(hashes, num_items)

        h = Series(h, index=obj.index, dtype="uint64", copy=False)
    else:
        raise TypeError(f"Unexpected type for hashing {type(obj)}")
    return h


def hash_tuples(vals, encoding="utf8", hash_key: str = _default_hash_key):
    """
    Hash an MultiIndex / list-of-tuples efficiently

    Parameters
    ----------
    vals : MultiIndex, list-of-tuples, or single tuple
    encoding : str, default 'utf8'
    hash_key : str, default _default_hash_key

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

    from pandas import Categorical, MultiIndex

    if not isinstance(vals, ABCMultiIndex):
        vals = MultiIndex.from_tuples(vals)

    # create a list-of-Categoricals
    vals = [
        Categorical(vals.codes[level], vals.levels[level], ordered=False, fastpath=True)
        for level in range(vals.nlevels)
    ]

    # hash the list-of-ndarrays
    hashes = (
        _hash_categorical(cat, encoding=encoding, hash_key=hash_key) for cat in vals
    )
    h = _combine_hash_arrays(hashes, len(vals))
    if is_tuple:
        h = h[0]

    return h


def hash_tuple(val, encoding: str = "utf8", hash_key: str = _default_hash_key):
    """
    Hash a single tuple efficiently

    Parameters
    ----------
    val : single tuple
    encoding : str, default 'utf8'
    hash_key : str, default _default_hash_key

    Returns
    -------
    hash

    """
    hashes = (_hash_scalar(v, encoding=encoding, hash_key=hash_key) for v in val)

    h = _combine_hash_arrays(hashes, len(val))[0]

    return h


def _hash_categorical(c, encoding: str, hash_key: str):
    """
    Hash a Categorical by hashing its categories, and then mapping the codes
    to the hashes

    Parameters
    ----------
    c : Categorical
    encoding : str
    hash_key : str

    Returns
    -------
    ndarray of hashed values array, same size as len(c)
    """
    # Convert ExtensionArrays to ndarrays
    values = np.asarray(c.categories.values)
    hashed = hash_array(values, encoding, hash_key, categorize=False)

    # we have uint64, as we don't directly support missing values
    # we don't want to use take_nd which will coerce to float
    # instead, directly construct the result with a
    # max(np.uint64) as the missing value indicator
    #
    # TODO: GH 15362

    mask = c.isna()
    if len(hashed):
        result = hashed.take(c.codes)
    else:
        result = np.zeros(len(mask), dtype="uint64")

    if mask.any():
        result[mask] = np.iinfo(np.uint64).max

    return result


def hash_array(
    vals,
    encoding: str = "utf8",
    hash_key: str = _default_hash_key,
    categorize: bool = True,
):
    """
    Given a 1d array, return an array of deterministic integers.

    Parameters
    ----------
    vals : ndarray, Categorical
    encoding : str, default 'utf8'
        Encoding for data & key when strings.
    hash_key : str, default _default_hash_key
        Hash_key for string key to encode.
    categorize : bool, default True
        Whether to first categorize object arrays before hashing. This is more
        efficient when the array contains duplicate values.

    Returns
    -------
    1d uint64 numpy array of hash values, same length as the vals
    """

    if not hasattr(vals, "dtype"):
        raise TypeError("must pass a ndarray-like")
    dtype = vals.dtype

    # For categoricals, we hash the categories, then remap the codes to the
    # hash values. (This check is above the complex check so that we don't ask
    # numpy if categorical is a subdtype of complex, as it will choke).
    if is_categorical_dtype(dtype):
        return _hash_categorical(vals, encoding, hash_key)
    elif is_extension_array_dtype(dtype):
        vals, _ = vals._values_for_factorize()
        dtype = vals.dtype

    # we'll be working with everything as 64-bit values, so handle this
    # 128-bit value early
    if np.issubdtype(dtype, np.complex128):
        return hash_array(np.real(vals)) + 23 * hash_array(np.imag(vals))

    # First, turn whatever array this is into unsigned 64-bit ints, if we can
    # manage it.
    elif isinstance(dtype, np.bool):
        vals = vals.astype("u8")
    elif issubclass(dtype.type, (np.datetime64, np.timedelta64)):
        vals = vals.view("i8").astype("u8", copy=False)
    elif issubclass(dtype.type, np.number) and dtype.itemsize <= 8:
        vals = vals.view("u{}".format(vals.dtype.itemsize)).astype("u8")
    else:
        # With repeated values, its MUCH faster to categorize object dtypes,
        # then hash and rename categories. We allow skipping the categorization
        # when the values are known/likely to be unique.
        if categorize:
            from pandas import factorize, Categorical, Index

            codes, categories = factorize(vals, sort=False)
            cat = Categorical(codes, Index(categories), ordered=False, fastpath=True)
            return _hash_categorical(cat, encoding, hash_key)

        try:
            vals = hashing.hash_object_array(vals, hash_key, encoding)
        except TypeError:
            # we have mixed types
            vals = hashing.hash_object_array(
                vals.astype(str).astype(object), hash_key, encoding
            )

    # Then, redistribute these 64-bit ints within the space of 64-bit ints
    vals ^= vals >> 30
    vals *= np.uint64(0xBF58476D1CE4E5B9)
    vals ^= vals >> 27
    vals *= np.uint64(0x94D049BB133111EB)
    vals ^= vals >> 31
    return vals


def _hash_scalar(
    val, encoding: str = "utf8", hash_key: str = _default_hash_key
) -> np.ndarray:
    """
    Hash scalar value.

    Parameters
    ----------
    val : scalar
    encoding : str, default "utf8"
    hash_key : str, default _default_hash_key

    Returns
    -------
    1d uint64 numpy array of hash value, of length 1
    """

    if isna(val):
        # this is to be consistent with the _hash_categorical implementation
        return np.array([np.iinfo(np.uint64).max], dtype="u8")

    if getattr(val, "tzinfo", None) is not None:
        # for tz-aware datetimes, we need the underlying naive UTC value and
        # not the tz aware object or pd extension type (as
        # infer_dtype_from_scalar would do)
        if not isinstance(val, Timestamp):
            val = Timestamp(val)
        val = val.tz_convert(None)

    dtype, val = infer_dtype_from_scalar(val)
    vals = np.array([val], dtype=dtype)

    return hash_array(vals, hash_key=hash_key, encoding=encoding, categorize=False)
