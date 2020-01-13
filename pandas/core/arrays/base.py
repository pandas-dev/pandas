"""An interface for extending pandas with custom arrays.

.. warning::

   This is an experimental API and subject to breaking changes
   without warning.
"""
import operator
from typing import Any, Callable, Dict, Optional, Sequence, Tuple, Union

import numpy as np

from pandas._libs import lib
from pandas._typing import ArrayLike
from pandas.compat import set_function_name
from pandas.compat.numpy import function as nv
from pandas.errors import AbstractMethodError
from pandas.util._decorators import Appender, Substitution
from pandas.util._validators import validate_fillna_kwargs

from pandas.core.dtypes.common import is_array_like, is_list_like
from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.generic import ABCExtensionArray, ABCIndexClass, ABCSeries
from pandas.core.dtypes.missing import isna

from pandas.core import ops
from pandas.core.algorithms import _factorize_array, unique
from pandas.core.missing import backfill_1d, pad_1d
from pandas.core.sorting import nargsort

_extension_array_shared_docs: Dict[str, str] = dict()


def try_cast_to_ea(cls_or_instance, obj, dtype=None):
    """
    Call to `_from_sequence` that returns the object unchanged on Exception.

    Parameters
    ----------
    cls_or_instance : ExtensionArray subclass or instance
    obj : arraylike
        Values to pass to cls._from_sequence
    dtype : ExtensionDtype, optional

    Returns
    -------
    ExtensionArray or obj
    """
    try:
        result = cls_or_instance._from_sequence(obj, dtype=dtype)
    except Exception:
        # We can't predict what downstream EA constructors may raise
        result = obj
    return result


class ExtensionArray:
    """
    Abstract base class for custom 1-D array types.

    pandas will recognize instances of this class as proper arrays
    with a custom type and will not attempt to coerce them to objects. They
    may be stored directly inside a :class:`DataFrame` or :class:`Series`.

    .. versionadded:: 0.23.0

    Attributes
    ----------
    dtype
    nbytes
    ndim
    shape

    Methods
    -------
    argsort
    astype
    copy
    dropna
    factorize
    fillna
    isna
    ravel
    repeat
    searchsorted
    shift
    take
    unique
    view
    _concat_same_type
    _formatter
    _from_factorized
    _from_sequence
    _from_sequence_of_strings
    _ndarray_values
    _reduce
    _values_for_argsort
    _values_for_factorize

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

    A default repr displaying the type, (truncated) data, length,
    and dtype is provided. It can be customized or replaced by
    by overriding:

    * __repr__ : A default repr for the ExtensionArray.
    * _formatter : Print scalars inside a Series or DataFrame.

    Some methods require casting the ExtensionArray to an ndarray of Python
    objects with ``self.astype(object)``, which may be expensive. When
    performance is a concern, we highly recommend overriding the following
    methods:

    * fillna
    * dropna
    * unique
    * factorize / _values_for_factorize
    * argsort / _values_for_argsort
    * searchsorted

    The remaining methods implemented on this class should be performant,
    as they only compose abstract methods. Still, a more efficient
    implementation may be available, and these methods can be overridden.

    One can implement methods to handle array reductions.

    * _reduce

    One can implement methods to handle parsing from strings that will be used
    in methods such as ``pandas.io.parsers.read_csv``.

    * _from_sequence_of_strings

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

    If implementing NumPy's ``__array_ufunc__`` interface, pandas expects
    that

    1. You defer by returning ``NotImplemented`` when any Series are present
       in `inputs`. Pandas will extract the arrays and call the ufunc again.
    2. You define a ``_HANDLED_TYPES`` tuple as an attribute on the class.
       Pandas inspect this to determine whether the ufunc is valid for the
       types present.

    See :ref:`extending.extension.ufunc` for more.
    """

    # '_typ' is for pandas.core.dtypes.generic.ABCExtensionArray.
    # Don't override this.
    _typ = "extension"

    # ------------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------------

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        """
        Construct a new ExtensionArray from a sequence of scalars.

        Parameters
        ----------
        scalars : Sequence
            Each element will be an instance of the scalar type for this
            array, ``cls.dtype.type``.
        dtype : dtype, optional
            Construct for this particular dtype. This should be a Dtype
            compatible with the ExtensionArray.
        copy : bool, default False
            If True, copy the underlying data.

        Returns
        -------
        ExtensionArray
        """
        raise AbstractMethodError(cls)

    @classmethod
    def _from_sequence_of_strings(cls, strings, dtype=None, copy=False):
        """Construct a new ExtensionArray from a sequence of strings.

        .. versionadded:: 0.24.0

        Parameters
        ----------
        strings : Sequence
            Each element will be an instance of the scalar type for this
            array, ``cls.dtype.type``.
        dtype : dtype, optional
            Construct for this particular dtype. This should be a Dtype
            compatible with the ExtensionArray.
        copy : bool, default False
            If True, copy the underlying data.

        Returns
        -------
        ExtensionArray
        """
        raise AbstractMethodError(cls)

    @classmethod
    def _from_factorized(cls, values, original):
        """
        Reconstruct an ExtensionArray after factorization.

        Parameters
        ----------
        values : ndarray
            An integer ndarray with the factorized values.
        original : ExtensionArray
            The original ExtensionArray that factorize was called on.

        See Also
        --------
        factorize
        ExtensionArray.factorize
        """
        raise AbstractMethodError(cls)

    # ------------------------------------------------------------------------
    # Must be a Sequence
    # ------------------------------------------------------------------------

    def __getitem__(self, item):
        # type (Any) -> Any
        """
        Select a subset of self.

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

    def __setitem__(self, key: Union[int, np.ndarray], value: Any) -> None:
        """
        Set one or more values inplace.

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
        # Note, also, that Series/DataFrame.where internally use __setitem__
        # on a copy of the data.
        raise NotImplementedError(f"{type(self)} does not implement __setitem__.")

    def __len__(self) -> int:
        """
        Length of this array

        Returns
        -------
        length : int
        """
        raise AbstractMethodError(self)

    def __iter__(self):
        """
        Iterate over elements of the array.
        """
        # This needs to be implemented so that pandas recognizes extension
        # arrays as list-like. The default implementation makes successive
        # calls to ``__getitem__``, which may be slower than necessary.
        for i in range(len(self)):
            yield self[i]

    def to_numpy(self, dtype=None, copy=False, na_value=lib.no_default):
        """
        Convert to a NumPy ndarray.

        .. versionadded:: 1.0.0

        This is similar to :meth:`numpy.asarray`, but may provide additional control
        over how the conversion is done.

        Parameters
        ----------
        dtype : str or numpy.dtype, optional
            The dtype to pass to :meth:`numpy.asarray`.
        copy : bool, default False
            Whether to ensure that the returned value is a not a view on
            another array. Note that ``copy=False`` does not *ensure* that
            ``to_numpy()`` is no-copy. Rather, ``copy=True`` ensure that
            a copy is made, even if not strictly necessary.
        na_value : Any, optional
            The value to use for missing values. The default value depends
            on `dtype` and the type of the array.

        Returns
        -------
        numpy.ndarray
        """
        result = np.asarray(self, dtype=dtype)
        if copy or na_value is not lib.no_default:
            result = result.copy()
        if na_value is not lib.no_default:
            result[self.isna()] = na_value
        return result

    # ------------------------------------------------------------------------
    # Required attributes
    # ------------------------------------------------------------------------

    @property
    def dtype(self) -> ExtensionDtype:
        """
        An instance of 'ExtensionDtype'.
        """
        raise AbstractMethodError(self)

    @property
    def shape(self) -> Tuple[int, ...]:
        """
        Return a tuple of the array dimensions.
        """
        return (len(self),)

    @property
    def ndim(self) -> int:
        """
        Extension Arrays are only allowed to be 1-dimensional.
        """
        return 1

    @property
    def nbytes(self) -> int:
        """
        The number of bytes needed to store this object in memory.
        """
        # If this is expensive to compute, return an approximate lower bound
        # on the number of bytes needed.
        raise AbstractMethodError(self)

    # ------------------------------------------------------------------------
    # Additional Methods
    # ------------------------------------------------------------------------

    def astype(self, dtype, copy=True):
        """
        Cast to a NumPy array with 'dtype'.

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

    def isna(self) -> ArrayLike:
        """
        A 1-D array indicating if each value is missing.

        Returns
        -------
        na_values : Union[np.ndarray, ExtensionArray]
            In most cases, this should return a NumPy ndarray. For
            exceptional cases like ``SparseArray``, where returning
            an ndarray would be expensive, an ExtensionArray may be
            returned.

        Notes
        -----
        If returning an ExtensionArray, then

        * ``na_values._is_boolean`` should be True
        * `na_values` should implement :func:`ExtensionArray._reduce`
        * ``na_values.any`` and ``na_values.all`` should be implemented
        """
        raise AbstractMethodError(self)

    def _values_for_argsort(self) -> np.ndarray:
        """
        Return values for sorting.

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

    def argsort(
        self, ascending: bool = True, kind: str = "quicksort", *args, **kwargs
    ) -> np.ndarray:
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
        ndarray
            Array of indices that sort ``self``. If NaN values are contained,
            NaN values are placed at the end.

        See Also
        --------
        numpy.argsort : Sorting implementation used internally.
        """
        # Implementor note: You have two places to override the behavior of
        # argsort.
        # 1. _values_for_argsort : construct the values passed to np.argsort
        # 2. argsort : total control over sorting.
        ascending = nv.validate_argsort_with_ascending(ascending, args, kwargs)

        result = nargsort(self, kind=kind, ascending=ascending, na_position="last")
        return result

    def fillna(self, value=None, method=None, limit=None):
        """
        Fill NA/NaN values using the specified method.

        Parameters
        ----------
        value : scalar, array-like
            If a scalar value is passed it is used to fill all missing values.
            Alternatively, an array-like 'value' can be given. It's expected
            that the array-like have the same length as 'self'.
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap.
        limit : int, default None
            If method is specified, this is the maximum number of consecutive
            NaN values to forward/backward fill. In other words, if there is
            a gap with more than this number of consecutive NaNs, it will only
            be partially filled. If method is not specified, this is the
            maximum number of entries along the entire axis where NaNs will be
            filled.

        Returns
        -------
        ExtensionArray
            With NA/NaN filled.
        """
        value, method = validate_fillna_kwargs(value, method)

        mask = self.isna()

        if is_array_like(value):
            if len(value) != len(self):
                raise ValueError(
                    f"Length of 'value' does not match. Got ({len(value)}) "
                    f"expected {len(self)}"
                )
            value = value[mask]

        if mask.any():
            if method is not None:
                func = pad_1d if method == "pad" else backfill_1d
                new_values = func(self.astype(object), limit=limit, mask=mask)
                new_values = self._from_sequence(new_values, dtype=self.dtype)
            else:
                # fill with value
                new_values = self.copy()
                new_values[mask] = value
        else:
            new_values = self.copy()
        return new_values

    def dropna(self):
        """
        Return ExtensionArray without NA values.

        Returns
        -------
        valid : ExtensionArray
        """
        return self[~self.isna()]

    def shift(self, periods: int = 1, fill_value: object = None) -> ABCExtensionArray:
        """
        Shift values by desired number.

        Newly introduced missing values are filled with
        ``self.dtype.na_value``.

        .. versionadded:: 0.24.0

        Parameters
        ----------
        periods : int, default 1
            The number of periods to shift. Negative values are allowed
            for shifting backwards.

        fill_value : object, optional
            The scalar value to use for newly introduced missing values.
            The default is ``self.dtype.na_value``.

            .. versionadded:: 0.24.0

        Returns
        -------
        ExtensionArray
            Shifted.

        Notes
        -----
        If ``self`` is empty or ``periods`` is 0, a copy of ``self`` is
        returned.

        If ``periods > len(self)``, then an array of size
        len(self) is returned, with all values filled with
        ``self.dtype.na_value``.
        """
        # Note: this implementation assumes that `self.dtype.na_value` can be
        # stored in an instance of your ExtensionArray with `self.dtype`.
        if not len(self) or periods == 0:
            return self.copy()

        if isna(fill_value):
            fill_value = self.dtype.na_value

        empty = self._from_sequence(
            [fill_value] * min(abs(periods), len(self)), dtype=self.dtype
        )
        if periods > 0:
            a = empty
            b = self[:-periods]
        else:
            a = self[abs(periods) :]
            b = empty
        return self._concat_same_type([a, b])

    def unique(self):
        """
        Compute the ExtensionArray of unique values.

        Returns
        -------
        uniques : ExtensionArray
        """
        uniques = unique(self.astype(object))
        return self._from_sequence(uniques, dtype=self.dtype)

    def searchsorted(self, value, side="left", sorter=None):
        """
        Find indices where elements should be inserted to maintain order.

        .. versionadded:: 0.24.0

        Find the indices into a sorted array `self` (a) such that, if the
        corresponding elements in `value` were inserted before the indices,
        the order of `self` would be preserved.

        Assuming that `self` is sorted:

        ======  ================================
        `side`  returned index `i` satisfies
        ======  ================================
        left    ``self[i-1] < value <= self[i]``
        right   ``self[i-1] <= value < self[i]``
        ======  ================================

        Parameters
        ----------
        value : array_like
            Values to insert into `self`.
        side : {'left', 'right'}, optional
            If 'left', the index of the first suitable location found is given.
            If 'right', return the last such index.  If there is no suitable
            index, return either 0 or N (where N is the length of `self`).
        sorter : 1-D array_like, optional
            Optional array of integer indices that sort array a into ascending
            order. They are typically the result of argsort.

        Returns
        -------
        array of ints
            Array of insertion points with the same shape as `value`.

        See Also
        --------
        numpy.searchsorted : Similar method from NumPy.
        """
        # Note: the base tests provided by pandas only test the basics.
        # We do not test
        # 1. Values outside the range of the `data_for_sorting` fixture
        # 2. Values between the values in the `data_for_sorting` fixture
        # 3. Missing values.
        arr = self.astype(object)
        return arr.searchsorted(value, side=side, sorter=sorter)

    def _values_for_factorize(self) -> Tuple[np.ndarray, Any]:
        """
        Return an array and missing value suitable for factorization.

        Returns
        -------
        values : ndarray

            An array suitable for factorization. This should maintain order
            and be a supported dtype (Float64, Int64, UInt64, String, Object).
            By default, the extension array is cast to object dtype.
        na_value : object
            The value in `values` to consider missing. This will be treated
            as NA in the factorization routines, so it will be coded as
            `na_sentinal` and not included in `uniques`. By default,
            ``np.nan`` is used.

        Notes
        -----
        The values returned by this method are also used in
        :func:`pandas.util.hash_pandas_object`.
        """
        return self.astype(object), np.nan

    def factorize(self, na_sentinel: int = -1) -> Tuple[np.ndarray, ABCExtensionArray]:
        """
        Encode the extension array as an enumerated type.

        Parameters
        ----------
        na_sentinel : int, default -1
            Value to use in the `codes` array to indicate missing values.

        Returns
        -------
        codes : ndarray
            An integer NumPy array that's an indexer into the original
            ExtensionArray.
        uniques : ExtensionArray
            An ExtensionArray containing the unique values of `self`.

            .. note::

               uniques will *not* contain an entry for the NA value of
               the ExtensionArray if there are any missing values present
               in `self`.

        See Also
        --------
        factorize : Top-level factorize method that dispatches here.

        Notes
        -----
        :meth:`pandas.factorize` offers a `sort` keyword as well.
        """
        # Implementer note: There are two ways to override the behavior of
        # pandas.factorize
        # 1. _values_for_factorize and _from_factorize.
        #    Specify the values passed to pandas' internal factorization
        #    routines, and how to convert from those values back to the
        #    original ExtensionArray.
        # 2. ExtensionArray.factorize.
        #    Complete control over factorization.
        arr, na_value = self._values_for_factorize()

        codes, uniques = _factorize_array(
            arr, na_sentinel=na_sentinel, na_value=na_value
        )

        uniques = self._from_factorized(uniques, self)
        return codes, uniques

    _extension_array_shared_docs[
        "repeat"
    ] = """
        Repeat elements of a %(klass)s.

        Returns a new %(klass)s where each element of the current %(klass)s
        is repeated consecutively a given number of times.

        Parameters
        ----------
        repeats : int or array of ints
            The number of repetitions for each element. This should be a
            non-negative integer. Repeating 0 times will return an empty
            %(klass)s.
        axis : None
            Must be ``None``. Has no effect but is accepted for compatibility
            with numpy.

        Returns
        -------
        repeated_array : %(klass)s
            Newly created %(klass)s with repeated elements.

        See Also
        --------
        Series.repeat : Equivalent function for Series.
        Index.repeat : Equivalent function for Index.
        numpy.repeat : Similar method for :class:`numpy.ndarray`.
        ExtensionArray.take : Take arbitrary positions.

        Examples
        --------
        >>> cat = pd.Categorical(['a', 'b', 'c'])
        >>> cat
        [a, b, c]
        Categories (3, object): [a, b, c]
        >>> cat.repeat(2)
        [a, a, b, b, c, c]
        Categories (3, object): [a, b, c]
        >>> cat.repeat([1, 2, 3])
        [a, b, b, c, c, c]
        Categories (3, object): [a, b, c]
        """

    @Substitution(klass="ExtensionArray")
    @Appender(_extension_array_shared_docs["repeat"])
    def repeat(self, repeats, axis=None):
        nv.validate_repeat(tuple(), dict(axis=axis))
        ind = np.arange(len(self)).repeat(repeats)
        return self.take(ind)

    # ------------------------------------------------------------------------
    # Indexing methods
    # ------------------------------------------------------------------------

    def take(
        self, indices: Sequence[int], allow_fill: bool = False, fill_value: Any = None
    ) -> ABCExtensionArray:
        """
        Take elements from an array.

        Parameters
        ----------
        indices : sequence of int
            Indices to be taken.
        allow_fill : bool, default False
            How to handle negative values in `indices`.

            * False: negative values in `indices` indicate positional indices
              from the right (the default). This is similar to
              :func:`numpy.take`.

            * True: negative values in `indices` indicate
              missing values. These values are set to `fill_value`. Any other
              other negative values raise a ``ValueError``.

        fill_value : any, optional
            Fill value to use for NA-indices when `allow_fill` is True.
            This may be ``None``, in which case the default NA value for
            the type, ``self.dtype.na_value``, is used.

            For many ExtensionArrays, there will be two representations of
            `fill_value`: a user-facing "boxed" scalar, and a low-level
            physical NA value. `fill_value` should be the user-facing version,
            and the implementation should handle translating that to the
            physical version for processing the take if necessary.

        Returns
        -------
        ExtensionArray

        Raises
        ------
        IndexError
            When the indices are out of bounds for the array.
        ValueError
            When `indices` contains negative values other than ``-1``
            and `allow_fill` is True.

        See Also
        --------
        numpy.take
        api.extensions.take

        Notes
        -----
        ExtensionArray.take is called by ``Series.__getitem__``, ``.loc``,
        ``iloc``, when `indices` is a sequence of values. Additionally,
        it's called by :meth:`Series.reindex`, or any other method
        that causes realignment, with a `fill_value`.

        Examples
        --------
        Here's an example implementation, which relies on casting the
        extension array to object dtype. This uses the helper method
        :func:`pandas.api.extensions.take`.

        .. code-block:: python

           def take(self, indices, allow_fill=False, fill_value=None):
               from pandas.core.algorithms import take

               # If the ExtensionArray is backed by an ndarray, then
               # just pass that here instead of coercing to object.
               data = self.astype(object)

               if allow_fill and fill_value is None:
                   fill_value = self.dtype.na_value

               # fill value should always be translated from the scalar
               # type for the array, to the physical storage type for
               # the data, before passing to take.

               result = take(data, indices, fill_value=fill_value,
                             allow_fill=allow_fill)
               return self._from_sequence(result, dtype=self.dtype)
        """
        # Implementer note: The `fill_value` parameter should be a user-facing
        # value, an instance of self.dtype.type. When passed `fill_value=None`,
        # the default of `self.dtype.na_value` should be used.
        # This may differ from the physical storage type your ExtensionArray
        # uses. In this case, your implementation is responsible for casting
        # the user-facing type to the storage type, before using
        # pandas.api.extensions.take
        raise AbstractMethodError(self)

    def copy(self) -> ABCExtensionArray:
        """
        Return a copy of the array.

        Returns
        -------
        ExtensionArray
        """
        raise AbstractMethodError(self)

    def view(self, dtype=None) -> Union[ABCExtensionArray, np.ndarray]:
        """
        Return a view on the array.

        Parameters
        ----------
        dtype : str, np.dtype, or ExtensionDtype, optional
            Default None.

        Returns
        -------
        ExtensionArray
            A view of the :class:`ExtensionArray`.
        """
        # NB:
        # - This must return a *new* object referencing the same data, not self.
        # - The only case that *must* be implemented is with dtype=None,
        #   giving a view with the same dtype as self.
        if dtype is not None:
            raise NotImplementedError(dtype)
        return self[:]

    # ------------------------------------------------------------------------
    # Printing
    # ------------------------------------------------------------------------

    def __repr__(self) -> str:
        from pandas.io.formats.printing import format_object_summary

        # the short repr has no trailing newline, while the truncated
        # repr does. So we include a newline in our template, and strip
        # any trailing newlines from format_object_summary
        data = format_object_summary(
            self, self._formatter(), indent_for_name=False
        ).rstrip(", \n")
        class_name = f"<{type(self).__name__}>\n"
        return f"{class_name}{data}\nLength: {len(self)}, dtype: {self.dtype}"

    def _formatter(self, boxed: bool = False) -> Callable[[Any], Optional[str]]:
        """Formatting function for scalar values.

        This is used in the default '__repr__'. The returned formatting
        function receives instances of your scalar type.

        Parameters
        ----------
        boxed : bool, default False
            An indicated for whether or not your array is being printed
            within a Series, DataFrame, or Index (True), or just by
            itself (False). This may be useful if you want scalar values
            to appear differently within a Series versus on its own (e.g.
            quoted or not).

        Returns
        -------
        Callable[[Any], str]
            A callable that gets instances of the scalar type and
            returns a string. By default, :func:`repr` is used
            when ``boxed=False`` and :func:`str` is used when
            ``boxed=True``.
        """
        if boxed:
            return str
        return repr

    # ------------------------------------------------------------------------
    # Reshaping
    # ------------------------------------------------------------------------

    def ravel(self, order="C") -> ABCExtensionArray:
        """
        Return a flattened view on this array.

        Parameters
        ----------
        order : {None, 'C', 'F', 'A', 'K'}, default 'C'

        Returns
        -------
        ExtensionArray

        Notes
        -----
        - Because ExtensionArrays are 1D-only, this is a no-op.
        - The "order" argument is ignored, is for compatibility with NumPy.
        """
        return self

    @classmethod
    def _concat_same_type(
        cls, to_concat: Sequence[ABCExtensionArray]
    ) -> ABCExtensionArray:
        """
        Concatenate multiple array.

        Parameters
        ----------
        to_concat : sequence of this type

        Returns
        -------
        ExtensionArray
        """
        raise AbstractMethodError(cls)

    # The _can_hold_na attribute is set to True so that pandas internals
    # will use the ExtensionDtype.na_value as the NA value in operations
    # such as take(), reindex(), shift(), etc.  In addition, those results
    # will then be of the ExtensionArray subclass rather than an array
    # of objects
    _can_hold_na = True

    @property
    def _ndarray_values(self) -> np.ndarray:
        """
        Internal pandas method for lossy conversion to a NumPy ndarray.

        This method is not part of the pandas interface.

        The expectation is that this is cheap to compute, and is primarily
        used for interacting with our indexers.

        Returns
        -------
        array : ndarray
        """
        return np.array(self)

    def _reduce(self, name, skipna=True, **kwargs):
        """
        Return a scalar result of performing the reduction operation.

        Parameters
        ----------
        name : str
            Name of the function, supported values are:
            { any, all, min, max, sum, mean, median, prod,
            std, var, sem, kurt, skew }.
        skipna : bool, default True
            If True, skip NaN values.
        **kwargs
            Additional keyword arguments passed to the reduction function.
            Currently, `ddof` is the only supported kwarg.

        Returns
        -------
        scalar

        Raises
        ------
        TypeError : subclass does not define reductions
        """
        raise TypeError(f"cannot perform {name} with type {self.dtype}")


class ExtensionOpsMixin:
    """
    A base class for linking the operators to their dunder names.

    .. note::

       You may want to set ``__array_priority__`` if you want your
       implementation to be called when involved in binary operations
       with NumPy arrays.
    """

    @classmethod
    def _add_arithmetic_ops(cls):
        cls.__add__ = cls._create_arithmetic_method(operator.add)
        cls.__radd__ = cls._create_arithmetic_method(ops.radd)
        cls.__sub__ = cls._create_arithmetic_method(operator.sub)
        cls.__rsub__ = cls._create_arithmetic_method(ops.rsub)
        cls.__mul__ = cls._create_arithmetic_method(operator.mul)
        cls.__rmul__ = cls._create_arithmetic_method(ops.rmul)
        cls.__pow__ = cls._create_arithmetic_method(operator.pow)
        cls.__rpow__ = cls._create_arithmetic_method(ops.rpow)
        cls.__mod__ = cls._create_arithmetic_method(operator.mod)
        cls.__rmod__ = cls._create_arithmetic_method(ops.rmod)
        cls.__floordiv__ = cls._create_arithmetic_method(operator.floordiv)
        cls.__rfloordiv__ = cls._create_arithmetic_method(ops.rfloordiv)
        cls.__truediv__ = cls._create_arithmetic_method(operator.truediv)
        cls.__rtruediv__ = cls._create_arithmetic_method(ops.rtruediv)
        cls.__divmod__ = cls._create_arithmetic_method(divmod)
        cls.__rdivmod__ = cls._create_arithmetic_method(ops.rdivmod)

    @classmethod
    def _add_comparison_ops(cls):
        cls.__eq__ = cls._create_comparison_method(operator.eq)
        cls.__ne__ = cls._create_comparison_method(operator.ne)
        cls.__lt__ = cls._create_comparison_method(operator.lt)
        cls.__gt__ = cls._create_comparison_method(operator.gt)
        cls.__le__ = cls._create_comparison_method(operator.le)
        cls.__ge__ = cls._create_comparison_method(operator.ge)

    @classmethod
    def _add_logical_ops(cls):
        cls.__and__ = cls._create_logical_method(operator.and_)
        cls.__rand__ = cls._create_logical_method(ops.rand_)
        cls.__or__ = cls._create_logical_method(operator.or_)
        cls.__ror__ = cls._create_logical_method(ops.ror_)
        cls.__xor__ = cls._create_logical_method(operator.xor)
        cls.__rxor__ = cls._create_logical_method(ops.rxor)


class ExtensionScalarOpsMixin(ExtensionOpsMixin):
    """
    A mixin for defining  ops on an ExtensionArray.

    It is assumed that the underlying scalar objects have the operators
    already defined.

    Notes
    -----
    If you have defined a subclass MyExtensionArray(ExtensionArray), then
    use MyExtensionArray(ExtensionArray, ExtensionScalarOpsMixin) to
    get the arithmetic operators.  After the definition of MyExtensionArray,
    insert the lines

    MyExtensionArray._add_arithmetic_ops()
    MyExtensionArray._add_comparison_ops()

    to link the operators to your class.

    .. note::

       You may want to set ``__array_priority__`` if you want your
       implementation to be called when involved in binary operations
       with NumPy arrays.
    """

    @classmethod
    def _create_method(cls, op, coerce_to_dtype=True):
        """
        A class method that returns a method that will correspond to an
        operator for an ExtensionArray subclass, by dispatching to the
        relevant operator defined on the individual elements of the
        ExtensionArray.

        Parameters
        ----------
        op : function
            An operator that takes arguments op(a, b)
        coerce_to_dtype : bool, default True
            boolean indicating whether to attempt to convert
            the result to the underlying ExtensionArray dtype.
            If it's not possible to create a new ExtensionArray with the
            values, an ndarray is returned instead.

        Returns
        -------
        Callable[[Any, Any], Union[ndarray, ExtensionArray]]
            A method that can be bound to a class. When used, the method
            receives the two arguments, one of which is the instance of
            this class, and should return an ExtensionArray or an ndarray.

            Returning an ndarray may be necessary when the result of the
            `op` cannot be stored in the ExtensionArray. The dtype of the
            ndarray uses NumPy's normal inference rules.

        Examples
        --------
        Given an ExtensionArray subclass called MyExtensionArray, use

        >>> __add__ = cls._create_method(operator.add)

        in the class definition of MyExtensionArray to create the operator
        for addition, that will be based on the operator implementation
        of the underlying elements of the ExtensionArray
        """

        def _binop(self, other):
            def convert_values(param):
                if isinstance(param, ExtensionArray) or is_list_like(param):
                    ovalues = param
                else:  # Assume its an object
                    ovalues = [param] * len(self)
                return ovalues

            if isinstance(other, (ABCSeries, ABCIndexClass)):
                # rely on pandas to unbox and dispatch to us
                return NotImplemented

            lvalues = self
            rvalues = convert_values(other)

            # If the operator is not defined for the underlying objects,
            # a TypeError should be raised
            res = [op(a, b) for (a, b) in zip(lvalues, rvalues)]

            def _maybe_convert(arr):
                if coerce_to_dtype:
                    # https://github.com/pandas-dev/pandas/issues/22850
                    # We catch all regular exceptions here, and fall back
                    # to an ndarray.
                    res = try_cast_to_ea(self, arr)
                    if not isinstance(res, type(self)):
                        # exception raised in _from_sequence; ensure we have ndarray
                        res = np.asarray(arr)
                else:
                    res = np.asarray(arr)
                return res

            if op.__name__ in {"divmod", "rdivmod"}:
                a, b = zip(*res)
                return _maybe_convert(a), _maybe_convert(b)

            return _maybe_convert(res)

        op_name = ops._get_op_name(op, True)
        return set_function_name(_binop, op_name, cls)

    @classmethod
    def _create_arithmetic_method(cls, op):
        return cls._create_method(op)

    @classmethod
    def _create_comparison_method(cls, op):
        return cls._create_method(op, coerce_to_dtype=False)
