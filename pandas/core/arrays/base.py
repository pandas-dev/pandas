"""An interface for extending pandas with custom arrays.

.. warning::

   This is an experimental API and subject to breaking changes
   without warning.
"""
import textwrap

import numpy as np

from pandas.errors import AbstractMethodError
from pandas.compat import _default_fill_value
from pandas.compat.numpy import function as nv
from pandas.util._decorators import Appender, Substitution

_not_implemented_message = "{} does not implement {}."


_take_docstring = textwrap.dedent("""\
Take elements from an array.

Parameters
----------
%(arr)s\
indexer : sequence of integers
    Indices to be taken. See Notes for how negative indicies
    are handled.
fill_value : any, optional
    Fill value to use for NA-indicies. This has a few behaviors.

    * fill_value is not specified : triggers NumPy's semantics
      where negative values in `indexer` mean slices from the end.
    * fill_value is NA : Fill positions where `indexer` is ``-1``
      with ``self.dtype.na_value``. Anything considered NA by
      :func:`pandas.isna` will result in ``self.dtype.na_value``
      being used to fill.
    * fill_value is not NA : Fill positions where `indexer` is ``-1``
      with `fill_value`.

Returns
-------
ExtensionArray

Raises
------
IndexError
    When the indexer is out of bounds for the array.
ValueError
    When the indexer contains negative values other than ``-1``
    and `fill_value` is specified.

Notes
-----
The meaning of negative values in `indexer` depends on the
`fill_value` argument. By default, we follow the behavior
:meth:`numpy.take` of where negative indices indicate slices
from the end.

When `fill_value` is specified, we follow pandas semantics of ``-1``
indicating a missing value. In this case, positions where `indexer`
is ``-1`` will be filled with `fill_value` or the default NA value
for this type.

ExtensionArray.take is called by ``Series.__getitem__``, ``.loc``,
``iloc``, when the indexer is a sequence of values. Additionally,
it's called by :meth:`Series.reindex` with a `fill_value`.

See Also
--------
numpy.take""")


class ExtensionArray(object):
    """Abstract base class for custom 1-D array types.

    pandas will recognize instances of this class as proper arrays
    with a custom type and will not attempt to coerce them to objects. They
    may be stored directly inside a :class:`DataFrame` or :class:`Series`.

    .. versionadded:: 0.23.0

    Notes
    -----
    The interface includes the following abstract methods that must be
    implemented by subclasses:

    * _from_sequence
    * _from_factorized
    * __getitem__
    * __len__
    * dtype
    * nbytes
    * isna
    * take
    * copy
    * _concat_same_type

    Some additional methods are available to satisfy pandas' internal, private
    block API:

    * _can_hold_na
    * _formatting_values

    Some methods require casting the ExtensionArray to an ndarray of Python
    objects with ``self.astype(object)``, which may be expensive. When
    performance is a concern, we highly recommend overriding the following
    methods:

    * fillna
    * unique
    * factorize / _values_for_factorize
    * argsort / _values_for_argsort
    * take / _values_for_take

    This class does not inherit from 'abc.ABCMeta' for performance reasons.
    Methods and properties required by the interface raise
    ``pandas.errors.AbstractMethodError`` and no ``register`` method is
    provided for registering virtual subclasses.

    ExtensionArrays are limited to 1 dimension.

    They may be backed by none, one, or many NumPy arrays. For example,
    ``pandas.Categorical`` is an extension array backed by two arrays,
    one for codes and one for categories. An array of IPv6 address may
    be backed by a NumPy structured array with two fields, one for the
    lower 64 bits and one for the upper 64 bits. Or they may be backed
    by some other storage type, like Python lists. Pandas makes no
    assumptions on how the data are stored, just that it can be converted
    to a NumPy array.
    The ExtensionArray interface does not impose any rules on how this data
    is stored. However, currently, the backing data cannot be stored in
    attributes called ``.values`` or ``._values`` to ensure full compatibility
    with pandas internals. But other names as ``.data``, ``._data``,
    ``._items``, ... can be freely used.
    """
    # '_typ' is for pandas.core.dtypes.generic.ABCExtensionArray.
    # Don't override this.
    _typ = 'extension'

    # ------------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------------
    @classmethod
    def _from_sequence(cls, scalars):
        """Construct a new ExtensionArray from a sequence of scalars.

        Parameters
        ----------
        scalars : Sequence
            Each element will be an instance of the scalar type for this
            array, ``cls.dtype.type``.
        Returns
        -------
        ExtensionArray
        """
        raise AbstractMethodError(cls)

    @classmethod
    def _from_factorized(cls, values, original):
        """Reconstruct an ExtensionArray after factorization.

        Parameters
        ----------
        values : ndarray
            An integer ndarray with the factorized values.
        original : ExtensionArray
            The original ExtensionArray that factorize was called on.

        See Also
        --------
        pandas.factorize
        ExtensionArray.factorize
        """
        raise AbstractMethodError(cls)

    # ------------------------------------------------------------------------
    # Must be a Sequence
    # ------------------------------------------------------------------------

    def __getitem__(self, item):
        # type (Any) -> Any
        """Select a subset of self.

        Parameters
        ----------
        item : int, slice, or ndarray
            * int: The position in 'self' to get.

            * slice: A slice object, where 'start', 'stop', and 'step' are
              integers or None

            * ndarray: A 1-d boolean NumPy ndarray the same length as 'self'

        Returns
        -------
        item : scalar or ExtensionArray

        Notes
        -----
        For scalar ``item``, return a scalar value suitable for the array's
        type. This should be an instance of ``self.dtype.type``.

        For slice ``key``, return an instance of ``ExtensionArray``, even
        if the slice is length 0 or 1.

        For a boolean mask, return an instance of ``ExtensionArray``, filtered
        to the values where ``item`` is True.
        """
        raise AbstractMethodError(self)

    def __setitem__(self, key, value):
        # type: (Union[int, np.ndarray], Any) -> None
        """Set one or more values inplace.

        This method is not required to satisfy the pandas extension array
        interface.

        Parameters
        ----------
        key : int, ndarray, or slice
            When called from, e.g. ``Series.__setitem__``, ``key`` will be
            one of

            * scalar int
            * ndarray of integers.
            * boolean ndarray
            * slice object

        value : ExtensionDtype.type, Sequence[ExtensionDtype.type], or object
            value or values to be set of ``key``.

        Returns
        -------
        None
        """
        # Some notes to the ExtensionArray implementor who may have ended up
        # here. While this method is not required for the interface, if you
        # *do* choose to implement __setitem__, then some semantics should be
        # observed:
        #
        # * Setting multiple values : ExtensionArrays should support setting
        #   multiple values at once, 'key' will be a sequence of integers and
        #  'value' will be a same-length sequence.
        #
        # * Broadcasting : For a sequence 'key' and a scalar 'value',
        #   each position in 'key' should be set to 'value'.
        #
        # * Coercion : Most users will expect basic coercion to work. For
        #   example, a string like '2018-01-01' is coerced to a datetime
        #   when setting on a datetime64ns array. In general, if the
        #   __init__ method coerces that value, then so should __setitem__
        raise NotImplementedError(_not_implemented_message.format(
            type(self), '__setitem__')
        )

    def __len__(self):
        """Length of this array

        Returns
        -------
        length : int
        """
        # type: () -> int
        raise AbstractMethodError(self)

    def __iter__(self):
        """Iterate over elements of the array.

        """
        # This needs to be implemented so that pandas recognizes extension
        # arrays as list-like. The default implementation makes successive
        # calls to ``__getitem__``, which may be slower than necessary.
        for i in range(len(self)):
            yield self[i]

    # ------------------------------------------------------------------------
    # Required attributes
    # ------------------------------------------------------------------------
    @property
    def dtype(self):
        # type: () -> ExtensionDtype
        """An instance of 'ExtensionDtype'."""
        raise AbstractMethodError(self)

    @property
    def shape(self):
        # type: () -> Tuple[int, ...]
        """Return a tuple of the array dimensions."""
        return (len(self),)

    @property
    def ndim(self):
        # type: () -> int
        """Extension Arrays are only allowed to be 1-dimensional."""
        return 1

    @property
    def nbytes(self):
        # type: () -> int
        """The number of bytes needed to store this object in memory.

        """
        # If this is expensive to compute, return an approximate lower bound
        # on the number of bytes needed.
        raise AbstractMethodError(self)

    # ------------------------------------------------------------------------
    # Additional Methods
    # ------------------------------------------------------------------------
    def astype(self, dtype, copy=True):
        """Cast to a NumPy array with 'dtype'.

        Parameters
        ----------
        dtype : str or dtype
            Typecode or data-type to which the array is cast.
        copy : bool, default True
            Whether to copy the data, even if not necessary. If False,
            a copy is made only if the old dtype does not match the
            new dtype.

        Returns
        -------
        array : ndarray
            NumPy ndarray with 'dtype' for its dtype.
        """
        return np.array(self, dtype=dtype, copy=copy)

    def isna(self):
        # type: () -> np.ndarray
        """Boolean NumPy array indicating if each value is missing.

        This should return a 1-D array the same length as 'self'.
        """
        raise AbstractMethodError(self)

    def _values_for_argsort(self):
        # type: () -> ndarray
        """Return values for sorting.

        Returns
        -------
        ndarray
            The transformed values should maintain the ordering between values
            within the array.

        See Also
        --------
        ExtensionArray.argsort
        """
        # Note: this is used in `ExtensionArray.argsort`.
        return np.array(self)

    def argsort(self, ascending=True, kind='quicksort', *args, **kwargs):
        """
        Return the indices that would sort this array.

        Parameters
        ----------
        ascending : bool, default True
            Whether the indices should result in an ascending
            or descending sort.
        kind : {'quicksort', 'mergesort', 'heapsort'}, optional
            Sorting algorithm.
        *args, **kwargs:
            passed through to :func:`numpy.argsort`.

        Returns
        -------
        index_array : ndarray
            Array of indices that sort ``self``.

        See Also
        --------
        numpy.argsort : Sorting implementation used internally.
        """
        # Implementor note: You have two places to override the behavior of
        # argsort.
        # 1. _values_for_argsort : construct the values passed to np.argsort
        # 2. argsort : total control over sorting.
        ascending = nv.validate_argsort_with_ascending(ascending, args, kwargs)
        values = self._values_for_argsort()
        result = np.argsort(values, kind=kind, **kwargs)
        if not ascending:
            result = result[::-1]
        return result

    def fillna(self, value=None, method=None, limit=None):
        """ Fill NA/NaN values using the specified method.

        Parameters
        ----------
        value : scalar, array-like
            If a scalar value is passed it is used to fill all missing values.
            Alternatively, an array-like 'value' can be given. It's expected
            that the array-like have the same length as 'self'.
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        limit : int, default None
            If method is specified, this is the maximum number of consecutive
            NaN values to forward/backward fill. In other words, if there is
            a gap with more than this number of consecutive NaNs, it will only
            be partially filled. If method is not specified, this is the
            maximum number of entries along the entire axis where NaNs will be
            filled.

        Returns
        -------
        filled : ExtensionArray with NA/NaN filled
        """
        from pandas.api.types import is_array_like
        from pandas.util._validators import validate_fillna_kwargs
        from pandas.core.missing import pad_1d, backfill_1d

        value, method = validate_fillna_kwargs(value, method)

        mask = self.isna()

        if is_array_like(value):
            if len(value) != len(self):
                raise ValueError("Length of 'value' does not match. Got ({}) "
                                 " expected {}".format(len(value), len(self)))
            value = value[mask]

        if mask.any():
            if method is not None:
                func = pad_1d if method == 'pad' else backfill_1d
                new_values = func(self.astype(object), limit=limit,
                                  mask=mask)
                new_values = self._from_sequence(new_values)
            else:
                # fill with value
                new_values = self.copy()
                new_values[mask] = value
        else:
            new_values = self.copy()
        return new_values

    def unique(self):
        """Compute the ExtensionArray of unique values.

        Returns
        -------
        uniques : ExtensionArray
        """
        from pandas import unique

        uniques = unique(self.astype(object))
        return self._from_sequence(uniques)

    def _values_for_factorize(self):
        # type: () -> Tuple[ndarray, Any]
        """Return an array and missing value suitable for factorization.

        Returns
        -------
        values : ndarray
            An array suitable for factoraization. This should maintain order
            and be a supported dtype (Float64, Int64, UInt64, String, Object).
            By default, the extension array is cast to object dtype.
        na_value : object
            The value in `values` to consider missing. This will be treated
            as NA in the factorization routines, so it will be coded as
            `na_sentinal` and not included in `uniques`. By default,
            ``np.nan`` is used.
        """
        return self.astype(object), np.nan

    def factorize(self, na_sentinel=-1):
        # type: (int) -> Tuple[ndarray, ExtensionArray]
        """Encode the extension array as an enumerated type.

        Parameters
        ----------
        na_sentinel : int, default -1
            Value to use in the `labels` array to indicate missing values.

        Returns
        -------
        labels : ndarray
            An interger NumPy array that's an indexer into the original
            ExtensionArray.
        uniques : ExtensionArray
            An ExtensionArray containing the unique values of `self`.

            .. note::

               uniques will *not* contain an entry for the NA value of
               the ExtensionArray if there are any missing values present
               in `self`.

        See Also
        --------
        pandas.factorize : Top-level factorize method that dispatches here.

        Notes
        -----
        :meth:`pandas.factorize` offers a `sort` keyword as well.
        """
        # Impelmentor note: There are two ways to override the behavior of
        # pandas.factorize
        # 1. _values_for_factorize and _from_factorize.
        #    Specify the values passed to pandas' internal factorization
        #    routines, and how to convert from those values back to the
        #    original ExtensionArray.
        # 2. ExtensionArray.factorize.
        #    Complete control over factorization.
        from pandas.core.algorithms import _factorize_array

        arr, na_value = self._values_for_factorize()

        labels, uniques = _factorize_array(arr, na_sentinel=na_sentinel,
                                           na_value=na_value)

        uniques = self._from_factorized(uniques, self)
        return labels, uniques

    # ------------------------------------------------------------------------
    # Indexing methods
    # ------------------------------------------------------------------------
    def _values_for_take(self):
        """Values to use for `take`.

        Coerces to object dtype by default.

        Returns
        -------
        array-like
            Must satisify NumPy's `take` semantics.
        """
        return self.astype(object)

    @Substitution(arr='')
    @Appender(_take_docstring)
    def take(self, indexer, fill_value=_default_fill_value):
        # type: (Sequence[int], Optional[Any]) -> ExtensionArray
        if fill_value is np.nan:
            import pdb; pdb.set_trace()
        from pandas.core.algorithms import take

        data = self._values_for_take()
        result = take(data, indexer, fill_value=fill_value)
        return self._from_sequence(result)

    def copy(self, deep=False):
        # type: (bool) -> ExtensionArray
        """Return a copy of the array.

        Parameters
        ----------
        deep : bool, default False
            Also copy the underlying data backing this array.

        Returns
        -------
        ExtensionArray
        """
        raise AbstractMethodError(self)

    # ------------------------------------------------------------------------
    # Block-related methods
    # ------------------------------------------------------------------------

    def _formatting_values(self):
        # type: () -> np.ndarray
        # At the moment, this has to be an array since we use result.dtype
        """An array of values to be printed in, e.g. the Series repr"""
        return np.array(self)

    @classmethod
    def _concat_same_type(cls, to_concat):
        # type: (Sequence[ExtensionArray]) -> ExtensionArray
        """Concatenate multiple array

        Parameters
        ----------
        to_concat : sequence of this type

        Returns
        -------
        ExtensionArray
        """
        raise AbstractMethodError(cls)

    @property
    def _can_hold_na(self):
        # type: () -> bool
        """Whether your array can hold missing values. True by default.

        Notes
        -----
        Setting this to false will optimize some operations like fillna.
        """
        return True

    @property
    def _ndarray_values(self):
        # type: () -> np.ndarray
        """Internal pandas method for lossy conversion to a NumPy ndarray.

        This method is not part of the pandas interface.

        The expectation is that this is cheap to compute, and is primarily
        used for interacting with our indexers.
        """
        return np.array(self)

