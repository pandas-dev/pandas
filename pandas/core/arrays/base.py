"""An interface for extending pandas with custom arrays."""
import numpy as np

from pandas.errors import AbstractMethodError
from pandas.compat.numpy import function as nv

_not_implemented_message = "{} does not implement {}."


class ExtensionArray(object):
    """Abstract base class for custom 1-D array types.

    pandas will recognize instances of this class as proper arrays
    with a custom type and will not attempt to coerce them to objects. They
    may be stored directly inside a :class:`DataFrame` or :class:`Series`.

    Notes
    -----
    The interface includes the following abstract methods that must be
    implemented by subclasses:

    * _constructor_from_sequence
    * __getitem__
    * __len__
    * dtype
    * nbytes
    * isna
    * take
    * copy
    * _concat_same_type

    Some additional methods are available to satisfy pandas' internal, private
    block API.

    * _can_hold_na
    * _formatting_values

    This class does not inherit from 'abc.ABCMeta' for performance reasons.
    Methods and properties required by the interface raise
    ``pandas.errors.AbstractMethodError`` and no ``register`` method is
    provided for registering virtual subclasses.

    ExtensionArrays are limited to 1 dimension.

    They may be backed by none, one, or many NumPy ararys. For example,
    ``pandas.Categorical`` is an extension array backed by two arrays,
    one for codes and one for categories. An array of IPv6 address may
    be backed by a NumPy structured array with two fields, one for the
    lower 64 bits and one for the upper 64 bits. Or they may be backed
    by some other storage type, like Python lists. Pandas makes no
    assumptions on how the data are stored, just that it can be converted
    to a NumPy array.

    Extension arrays should be able to be constructed with instances of
    the class, i.e. ``ExtensionArray(extension_array)`` should return
    an instance, not error.
    """
    # '_typ' is for pandas.core.dtypes.generic.ABCExtensionArray.
    # Don't override this.
    _typ = 'extension'

    # ------------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------------
    @classmethod
    def _constructor_from_sequence(cls, scalars):
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

    def _simple_ndarray(self):
        # type: () -> ndarray
        """Convert the array to a simple ndarray representaiton.

        Many methods can operate indirectly on a cheap-to-compute array that
        is somehow representative of the extension array. For example, rather
        than sorting an ExtensionArray directly, which might be expensive,
        we could convert the ExtensionArray to a representative ndarray of
        integers, sort the integers, and perform a ``take``.

        The coversion between ExtensionArray and the simple ndarray should be
        strictly monotonic https://en.wikipedia.org/wiki/Monotonic_function,
        and as cheap to compute as possible.

        Returns
        -------
        values : ndarray

        See Also
        --------
        ExtensionArray.argsort
        """
        # Implemnetor note: This method is currently used in
        # - ExtensionArray.argsort

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
        # 1. _simple_ndarray : construct the values passed to np.argsort
        # 2. argsort : total control over sorting.

        ascending = nv.validate_argsort_with_ascending(ascending, args, kwargs)
        values = self._simple_ndarray()
        result = np.argsort(values, kind=kind, **kwargs)
        if not ascending:
            result = result[::-1]
        return result

    def unique(self):
        """Compute the ExtensionArray of unique values.

        Returns
        -------
        uniques : ExtensionArray
        """
        from pandas import unique

        uniques = unique(self.astype(object))
        return self._constructor_from_sequence(uniques)

    # ------------------------------------------------------------------------
    # Indexing methods
    # ------------------------------------------------------------------------
    def take(self, indexer, allow_fill=True, fill_value=None):
        # type: (Sequence[int], bool, Optional[Any]) -> ExtensionArray
        """Take elements from an array.

        Parameters
        ----------
        indexer : sequence of integers
            indices to be taken. -1 is used to indicate values
            that are missing.
        allow_fill : bool, default True
            If False, indexer is assumed to contain no -1 values so no filling
            will be done. This short-circuits computation of a mask. Result is
            undefined if allow_fill == False and -1 is present in indexer.
        fill_value : any, default None
            Fill value to replace -1 values with. If applicable, this should
            use the sentinel missing value for this type.

        Notes
        -----
        This should follow pandas' semantics where -1 indicates missing values.
        Positions where indexer is ``-1`` should be filled with the missing
        value for this type.

        This is called by ``Series.__getitem__``, ``.loc``, ``iloc``, when the
        indexer is a sequence of values.

        Examples
        --------
        Suppose the extension array is backed by a NumPy array stored as
        ``self.data``. Then ``take`` may be written as

        .. code-block:: python

           def take(self, indexer, allow_fill=True, fill_value=None):
               mask = indexer == -1
               result = self.data.take(indexer)
               result[mask] = np.nan  # NA for this type
               return type(self)(result)

        See Also
        --------
        numpy.take
        """
        raise AbstractMethodError(self)

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
