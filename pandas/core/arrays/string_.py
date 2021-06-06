from __future__ import annotations

from typing import (
    TYPE_CHECKING,
    Any,
    Sequence,
    TypeVar,
    cast,
)

import numpy as np

from pandas._config import get_option

from pandas._libs import (
    lib,
    missing as libmissing,
)
from pandas._libs.arrays import NDArrayBacked
from pandas._typing import (
    Dtype,
    NpDtype,
    PositionalIndexer,
    Scalar,
    type_t,
)
from pandas.compat import pa_version_under1p0
from pandas.compat.numpy import function as nv

from pandas.core.dtypes.base import (
    ExtensionDtype,
    register_extension_dtype,
)
from pandas.core.dtypes.common import (
    is_array_like,
    is_bool_dtype,
    is_dtype_equal,
    is_integer_dtype,
    is_object_dtype,
    is_string_dtype,
    pandas_dtype,
)

from pandas.core import ops
from pandas.core.array_algos import masked_reductions
from pandas.core.arraylike import OpsMixin
from pandas.core.arrays import (
    FloatingArray,
    IntegerArray,
    PandasArray,
)
from pandas.core.arrays.base import ExtensionArray
from pandas.core.arrays.floating import FloatingDtype
from pandas.core.arrays.integer import _IntegerDtype
from pandas.core.construction import extract_array
from pandas.core.indexers import check_array_indexer
from pandas.core.missing import isna
from pandas.core.strings.object_array import ObjectStringArrayMixin

if TYPE_CHECKING:
    from typing import Literal

    import pyarrow

    from pandas.core.arrays.string_arrow import ArrowStringArray

    StringStorage = Literal["python", "pyarrow"]


def _validate_string_storage(storage: StringStorage) -> None:
    if storage not in {"python", "pyarrow"}:
        raise ValueError(
            f"Storage must be 'python' or 'pyarrow'. Got {storage} instead."
        )
    if storage == "pyarrow" and pa_version_under1p0:
        raise ImportError("pyarrow>=1.0.0 is required for PyArrow backed StringArray.")


def _get_string_storage(storage: StringStorage | None) -> StringStorage:
    if storage is None:
        storage = get_option("mode.string_storage")
    _validate_string_storage(storage)
    return storage


@register_extension_dtype
class StringDtype(ExtensionDtype):
    """
    Extension dtype for string data.

    .. versionadded:: 1.0.0

    .. warning::

       StringDtype is considered experimental. The implementation and
       parts of the API may change without warning.

       In particular, StringDtype.na_value may change to no longer be
       ``numpy.nan``.

    Attributes
    ----------
    None

    Methods
    -------
    None

    Examples
    --------
    >>> pd.StringDtype()
    StringDtype
    """

    name = "string"

    #: StringDtype.na_value uses pandas.NA
    na_value = libmissing.NA

    @property
    def type(self) -> type[str]:
        return str

    @classmethod
    def construct_array_type(cls) -> type_t[StringArray]:
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return StringArray

    def __repr__(self) -> str:
        return "StringDtype"

    def __from_arrow__(
        self, array: pyarrow.Array | pyarrow.ChunkedArray
    ) -> StringArray:
        """
        Construct StringArray from pyarrow Array/ChunkedArray.
        """
        import pyarrow

        if isinstance(array, pyarrow.Array):
            chunks = [array]
        else:
            # pyarrow.ChunkedArray
            chunks = array.chunks

        results = []
        for arr in chunks:
            # using _from_sequence to ensure None is converted to NA
            str_arr = StringArray._from_sequence(np.array(arr))
            results.append(str_arr)

        if results:
            return StringArray._concat_same_type(results)
        else:
            return StringArray(np.array([], dtype="object"))


StringArrayT = TypeVar("StringArrayT", bound="StringArray")


class StringArray(OpsMixin, ExtensionArray, ObjectStringArrayMixin):
    """
    Extension array for string data.

    .. versionadded:: 1.0.0

    .. warning::

       StringArray is considered experimental. The implementation and
       parts of the API may change without warning.

    Parameters
    ----------
    values : array-like
        The array of data.

        .. warning::

           Currently, this expects an object-dtype ndarray
           where the elements are Python strings or :attr:`pandas.NA`.
           This may change without warning in the future. Use
           :meth:`pandas.array` with ``dtype="string"`` for a stable way of
           creating a `StringArray` from any sequence.

    copy : bool, default False
        Whether to copy the array of data.

    storage : {"python", "pyarrow"}, optional
        If not given, the value of ``pd.options.mode.string_storage``.

    Attributes
    ----------
    None

    Methods
    -------
    None

    See Also
    --------
    array
        The recommended function for creating a StringArray.
    Series.str
        The string methods are available on Series backed by
        a StringArray.

    Notes
    -----
    StringArray returns a BooleanArray for comparison methods.

    Examples
    --------
    >>> pd.array(['This is', 'some text', None, 'data.'], dtype="string")
    <StringArray>
    ['This is', 'some text', <NA>, 'data.']
    Length: 4, dtype: string

    Unlike arrays instantiated with ``dtype="object"``, ``StringArray``
    will convert the values to strings.

    >>> pd.array(['1', 1], dtype="object")
    <PandasArray>
    ['1', 1]
    Length: 2, dtype: object
    >>> pd.array(['1', 1], dtype="string")
    <StringArray>
    ['1', '1']
    Length: 2, dtype: string

    However, instantiating StringArrays directly with non-strings will raise an error.

    For comparison methods, `StringArray` returns a :class:`pandas.BooleanArray`:

    >>> pd.array(["a", None, "c"], dtype="string") == "a"
    <BooleanArray>
    [True, <NA>, False]
    Length: 3, dtype: boolean
    """

    _dtype = StringDtype()
    _array: ObjectStringArray | ArrowStringArray
    _storage: StringStorage

    @property
    def storage(self) -> str:
        return self._storage

    # ------------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------------

    def __init__(self, values, copy=False, *, storage: StringStorage | None = None):
        from pandas.core.arrays.string_arrow import ArrowStringArray

        storage = _get_string_storage(storage)
        self._storage = storage
        klass = ObjectStringArray if storage == "python" else ArrowStringArray
        # error: Incompatible types in assignment (expression has type
        # "ObjectStringArrayMixin", variable has type "Union[ObjectStringArray,
        # ArrowStringArray]")
        self._array = klass(values, copy=copy)  # type: ignore[assignment]

    def _from_array(self, array):
        klass = type(self)
        new_string_array = klass.__new__(klass)
        new_string_array._storage = self._storage
        new_string_array._array = array
        return new_string_array

    def _maybe_wrap_result(self, result):
        if isinstance(result, type(self._array)):
            return self._from_array(result)
        return result

    @classmethod
    def _from_sequence(
        cls,
        scalars,
        *,
        dtype: Dtype | None = None,
        copy=False,
        storage: StringStorage | None = None,
    ):
        from pandas.core.arrays.string_arrow import ArrowStringArray

        if dtype:
            assert dtype == "string"

        new_string_array = cls.__new__(cls)
        storage = _get_string_storage(storage)
        new_string_array._storage = storage
        klass = ObjectStringArray if storage == "python" else ArrowStringArray
        # error: "Type[ObjectStringArrayMixin]" has no attribute "_from_sequence"
        new_string_array._array = klass._from_sequence(  # type: ignore[attr-defined]
            scalars, dtype=dtype, copy=copy
        )
        return new_string_array

    @classmethod
    def _from_sequence_of_strings(
        cls,
        strings,
        *,
        dtype: Dtype | None = None,
        copy=False,
        storage: StringStorage | None = None,
    ):
        from pandas.core.arrays.string_arrow import ArrowStringArray

        if dtype:
            assert dtype == "string"

        new_string_array = cls.__new__(cls)
        storage = _get_string_storage(storage)
        new_string_array._storage = storage
        klass = ObjectStringArray if storage == "python" else ArrowStringArray
        # error: "Type[ObjectStringArrayMixin]" has no attribute
        # "_from_sequence_of_strings"
        tmp = klass._from_sequence_of_strings  # type: ignore[attr-defined]
        new_string_array._array = tmp(strings, dtype=dtype, copy=copy)
        return new_string_array

    # ------------------------------------------------------------------------
    # Must be a Sequence
    # ------------------------------------------------------------------------

    def __getitem__(self, item: PositionalIndexer) -> Any:
        result = self._array[item]
        return self._maybe_wrap_result(result)

    def __setitem__(self, key, value) -> None:
        if isinstance(value, type(self)):
            value = value._array
        self._array[key] = value

    def __len__(self) -> int:
        return len(self._array)

    def to_numpy(
        self,
        dtype=None,
        copy: bool = False,
        na_value=lib.no_default,
    ) -> np.ndarray:
        return self._array.to_numpy(dtype=dtype, copy=copy, na_value=na_value)

    # ------------------------------------------------------------------------
    # Required attributes
    # ------------------------------------------------------------------------

    @property
    def dtype(self) -> StringDtype:
        return self._dtype

    @property
    def nbytes(self) -> int:
        return self._array.nbytes

    # ------------------------------------------------------------------------
    # Additional Methods
    # ------------------------------------------------------------------------

    def astype(self, dtype, copy=True):
        dtype = pandas_dtype(dtype)

        if is_dtype_equal(dtype, self.dtype):
            if copy:
                return self.copy()
            return self

        return self._array.astype(dtype, copy=copy)

    def isna(self) -> np.ndarray:
        return self._array.isna()

    # def _values_for_argsort(self) -> np.ndarray:
    #     return self._array._values_for_argsort()

    # def argsort(
    #     self,
    #     ascending: bool = True,
    #     kind: str = "quicksort",
    #     na_position: str = "last",
    #     *args,
    #     **kwargs,
    # ) -> np.ndarray:

    # def argmin(self, skipna: bool = True) -> int:

    # def argmax(self, skipna: bool = True) -> int:

    # def fillna(
    #     self,
    #     value: object | ArrayLike | None = None,
    #     method: FillnaOptions | None = None,
    #     limit: int | None = None,
    # ):

    # def dropna(self):
    #     """
    #     Return ExtensionArray without NA values.

    #     Returns
    #     -------
    #     valid : ExtensionArray
    #     """
    #     # error: Unsupported operand type for ~ ("ExtensionArray")
    #     return self[~self.isna()]  # type: ignore[operator]

    # def shift(self, periods: int = 1, fill_value: object = None) -> ExtensionArray:
    #     """
    #     Shift values by desired number.

    #     Newly introduced missing values are filled with
    #     ``self.dtype.na_value``.

    #     .. versionadded:: 0.24.0

    #     Parameters
    #     ----------
    #     periods : int, default 1
    #         The number of periods to shift. Negative values are allowed
    #         for shifting backwards.

    #     fill_value : object, optional
    #         The scalar value to use for newly introduced missing values.
    #         The default is ``self.dtype.na_value``.

    #         .. versionadded:: 0.24.0

    #     Returns
    #     -------
    #     ExtensionArray
    #         Shifted.

    #     Notes
    #     -----
    #     If ``self`` is empty or ``periods`` is 0, a copy of ``self`` is
    #     returned.

    #     If ``periods > len(self)``, then an array of size
    #     len(self) is returned, with all values filled with
    #     ``self.dtype.na_value``.
    #     """
    #     # Note: this implementation assumes that `self.dtype.na_value` can be
    #     # stored in an instance of your ExtensionArray with `self.dtype`.
    #     if not len(self) or periods == 0:
    #         return self.copy()

    #     if isna(fill_value):
    #         fill_value = self.dtype.na_value

    #     empty = self._from_sequence(
    #         [fill_value] * min(abs(periods), len(self)), dtype=self.dtype
    #     )
    #     if periods > 0:
    #         a = empty
    #         b = self[:-periods]
    #     else:
    #         a = self[abs(periods) :]
    #         b = empty
    #     return self._concat_same_type([a, b])

    # def unique(self: ExtensionArrayT) -> ExtensionArrayT:
    #     """
    #     Compute the ExtensionArray of unique values.

    #     Returns
    #     -------
    #     uniques : ExtensionArray
    #     """
    #     uniques = unique(self.astype(object))
    #     return self._from_sequence(uniques, dtype=self.dtype)

    # def searchsorted(self, value, side="left", sorter=None):
    #     """
    #     Find indices where elements should be inserted to maintain order.

    #     .. versionadded:: 0.24.0

    #     Find the indices into a sorted array `self` (a) such that, if the
    #     corresponding elements in `value` were inserted before the indices,
    #     the order of `self` would be preserved.

    #     Assuming that `self` is sorted:

    #     ======  ================================
    #     `side`  returned index `i` satisfies
    #     ======  ================================
    #     left    ``self[i-1] < value <= self[i]``
    #     right   ``self[i-1] <= value < self[i]``
    #     ======  ================================

    #     Parameters
    #     ----------
    #     value : array_like
    #         Values to insert into `self`.
    #     side : {'left', 'right'}, optional
    #         If 'left', the index of the first suitable location found is given.
    #         If 'right', return the last such index.  If there is no suitable
    #         index, return either 0 or N (where N is the length of `self`).
    #     sorter : 1-D array_like, optional
    #         Optional array of integer indices that sort array a into ascending
    #         order. They are typically the result of argsort.

    #     Returns
    #     -------
    #     array of ints
    #         Array of insertion points with the same shape as `value`.

    #     See Also
    #     --------
    #     numpy.searchsorted : Similar method from NumPy.
    #     """
    #     # Note: the base tests provided by pandas only test the basics.
    #     # We do not test
    #     # 1. Values outside the range of the `data_for_sorting` fixture
    #     # 2. Values between the values in the `data_for_sorting` fixture
    #     # 3. Missing values.
    #     arr = self.astype(object)
    #     return arr.searchsorted(value, side=side, sorter=sorter)

    def equals(self, other: object) -> bool:
        # TODO: allow ObjectStringArray and ArrowStringArray to compare equal
        if isinstance(other, type(self)):
            other = other._array
        return self._array.equals(other)

    # def isin(self, values) -> np.ndarray:
    #     """
    #     Pointwise comparison for set containment in the given values.

    #     Roughly equivalent to `np.array([x in values for x in self])`

    #     Parameters
    #     ----------
    #     values : Sequence

    #     Returns
    #     -------
    #     np.ndarray[bool]
    #     """
    #     return isin(np.asarray(self), values)

    # def _values_for_factorize(self) -> tuple[np.ndarray, Any]:
    #     """
    #     Return an array and missing value suitable for factorization.

    #     Returns
    #     -------
    #     values : ndarray

    #         An array suitable for factorization. This should maintain order
    #         and be a supported dtype (Float64, Int64, UInt64, String, Object).
    #         By default, the extension array is cast to object dtype.
    #     na_value : object
    #         The value in `values` to consider missing. This will be treated
    #         as NA in the factorization routines, so it will be coded as
    #         `na_sentinel` and not included in `uniques`. By default,
    #         ``np.nan`` is used.

    #     Notes
    #     -----
    #     The values returned by this method are also used in
    #     :func:`pandas.util.hash_pandas_object`.
    #     """
    #     return self.astype(object), np.nan

    def factorize(self, na_sentinel: int = -1) -> tuple[np.ndarray, ExtensionArray]:
        return self._array.factorize(na_sentinel=na_sentinel)

    # def factorize(self, na_sentinel: int = -1) -> tuple[np.ndarray, ExtensionArray]:
    #     """
    #     Encode the extension array as an enumerated type.

    #     Parameters
    #     ----------
    #     na_sentinel : int, default -1
    #         Value to use in the `codes` array to indicate missing values.

    #     Returns
    #     -------
    #     codes : ndarray
    #         An integer NumPy array that's an indexer into the original
    #         ExtensionArray.
    #     uniques : ExtensionArray
    #         An ExtensionArray containing the unique values of `self`.

    #         .. note::

    #            uniques will *not* contain an entry for the NA value of
    #            the ExtensionArray if there are any missing values present
    #            in `self`.

    #     See Also
    #     --------
    #     factorize : Top-level factorize method that dispatches here.

    #     Notes
    #     -----
    #     :meth:`pandas.factorize` offers a `sort` keyword as well.
    #     """
    #     # Implementer note: There are two ways to override the behavior of
    #     # pandas.factorize
    #     # 1. _values_for_factorize and _from_factorize.
    #     #    Specify the values passed to pandas' internal factorization
    #     #    routines, and how to convert from those values back to the
    #     #    original ExtensionArray.
    #     # 2. ExtensionArray.factorize.
    #     #    Complete control over factorization.
    #     arr, na_value = self._values_for_factorize()

    #     codes, uniques = factorize_array(
    #         arr, na_sentinel=na_sentinel, na_value=na_value
    #     )

    #     uniques = self._from_factorized(uniques, self)
    #     # error: Incompatible return value type (got "Tuple[ndarray, ndarray]",
    #     # expected "Tuple[ndarray, ExtensionArray]")
    #     return codes, uniques  # type: ignore[return-value]

    # _extension_array_shared_docs[
    #     "repeat"
    # ] = """
    #     Repeat elements of a %(klass)s.

    #     Returns a new %(klass)s where each element of the current %(klass)s
    #     is repeated consecutively a given number of times.

    #     Parameters
    #     ----------
    #     repeats : int or array of ints
    #         The number of repetitions for each element. This should be a
    #         non-negative integer. Repeating 0 times will return an empty
    #         %(klass)s.
    #     axis : None
    #         Must be ``None``. Has no effect but is accepted for compatibility
    #         with numpy.

    #     Returns
    #     -------
    #     repeated_array : %(klass)s
    #         Newly created %(klass)s with repeated elements.

    #     See Also
    #     --------
    #     Series.repeat : Equivalent function for Series.
    #     Index.repeat : Equivalent function for Index.
    #     numpy.repeat : Similar method for :class:`numpy.ndarray`.
    #     ExtensionArray.take : Take arbitrary positions.

    #     Examples
    #     --------
    #     >>> cat = pd.Categorical(['a', 'b', 'c'])
    #     >>> cat
    #     ['a', 'b', 'c']
    #     Categories (3, object): ['a', 'b', 'c']
    #     >>> cat.repeat(2)
    #     ['a', 'a', 'b', 'b', 'c', 'c']
    #     Categories (3, object): ['a', 'b', 'c']
    #     >>> cat.repeat([1, 2, 3])
    #     ['a', 'b', 'b', 'c', 'c', 'c']
    #     Categories (3, object): ['a', 'b', 'c']
    #     """

    # @Substitution(klass="ExtensionArray")
    # @Appender(_extension_array_shared_docs["repeat"])
    # def repeat(self, repeats: int | Sequence[int], axis: int | None = None):
    #     nv.validate_repeat((), {"axis": axis})
    #     ind = np.arange(len(self)).repeat(repeats)
    #     return self.take(ind)

    # ------------------------------------------------------------------------
    # Indexing methods
    # ------------------------------------------------------------------------

    def take(
        self, indices: Sequence[int], allow_fill: bool = False, fill_value: Any = None
    ):
        result = self._array.take(indices, allow_fill=allow_fill, fill_value=fill_value)
        return self._from_array(result)

    def copy(self) -> StringArray:
        result = self._array.copy()
        return self._from_array(result)

    # def view(self, dtype: Dtype | None = None) -> ArrayLike:
    #     """
    #     Return a view on the array.

    #     Parameters
    #     ----------
    #     dtype : str, np.dtype, or ExtensionDtype, optional
    #         Default None.

    #     Returns
    #     -------
    #     ExtensionArray or np.ndarray
    #         A view on the :class:`ExtensionArray`'s data.
    #     """
    #     # NB:
    #     # - This must return a *new* object referencing the same data, not self.
    #     # - The only case that *must* be implemented is with dtype=None,
    #     #   giving a view with the same dtype as self.
    #     if dtype is not None:
    #         raise NotImplementedError(dtype)
    #     return self[:]

    # ------------------------------------------------------------------------
    # Printing
    # ------------------------------------------------------------------------

    # def __repr__(self) -> str:
    #     from pandas.io.formats.printing import format_object_summary

    #     # the short repr has no trailing newline, while the truncated
    #     # repr does. So we include a newline in our template, and strip
    #     # any trailing newlines from format_object_summary
    #     data = format_object_summary(
    #         self, self._formatter(), indent_for_name=False
    #     ).rstrip(", \n")
    #     class_name = f"<{type(self).__name__}>\n"
    #     return f"{class_name}{data}\nLength: {len(self)}, dtype: {self.dtype}"

    # def _formatter(self, boxed: bool = False) -> Callable[[Any], str | None]:
    #     """
    #     Formatting function for scalar values.

    #     This is used in the default '__repr__'. The returned formatting
    #     function receives instances of your scalar type.

    #     Parameters
    #     ----------
    #     boxed : bool, default False
    #         An indicated for whether or not your array is being printed
    #         within a Series, DataFrame, or Index (True), or just by
    #         itself (False). This may be useful if you want scalar values
    #         to appear differently within a Series versus on its own (e.g.
    #         quoted or not).

    #     Returns
    #     -------
    #     Callable[[Any], str]
    #         A callable that gets instances of the scalar type and
    #         returns a string. By default, :func:`repr` is used
    #         when ``boxed=False`` and :func:`str` is used when
    #         ``boxed=True``.
    #     """
    #     if boxed:
    #         return str
    #     return repr

    # ------------------------------------------------------------------------
    # Reshaping
    # ------------------------------------------------------------------------

    @classmethod
    def _concat_same_type(
        cls: type[StringArrayT], to_concat: Sequence[StringArrayT]
    ) -> StringArrayT:
        from pandas.core.arrays.string_arrow import ArrowStringArray

        result: ObjectStringArray | ArrowStringArray
        if all(arr.storage == "python" for arr in to_concat):
            to_concat_object = cast(
                Sequence[ObjectStringArray], [arr._array for arr in to_concat]
            )
            result = ObjectStringArray._concat_same_type(to_concat_object)
            storage = "python"
        elif all(arr.storage == "pyarrow" for arr in to_concat):
            to_concat_arrow = [arr._array for arr in to_concat]
            result = ArrowStringArray._concat_same_type(to_concat_arrow)
            storage = "pyarrow"
        else:
            raise NotImplementedError

        new_string_array = cls.__new__(cls)
        new_string_array._storage = storage
        new_string_array._array = result
        return new_string_array

    # The _can_hold_na attribute is set to True so that pandas internals
    # will use the ExtensionDtype.na_value as the NA value in operations
    # such as take(), reindex(), shift(), etc.  In addition, those results
    # will then be of the ExtensionArray subclass rather than an array
    # of objects
    _can_hold_na = True

    def _reduce(self, name: str, *, skipna: bool = True, **kwargs):
        return self._array._reduce(name, skipna=skipna, **kwargs)

    # def __hash__(self) -> int:
    #     raise TypeError(f"unhashable type: {repr(type(self).__name__)}")

    # ------------------------------------------------------------------------
    # Other
    # ------------------------------------------------------------------------

    # @classmethod
    # def _empty(cls, shape, dtype) -> StringArray:
    #     values = np.empty(shape, dtype=object)
    #     values[:] = libmissing.NA
    #     return cls(values).astype(dtype, copy=False)

    # def _values_for_factorize(self):
    #     arr = self._ndarray.copy()
    #     mask = self.isna()
    #     arr[mask] = -1
    #     return arr, -1

    # ------------------------------------------------------------------------
    # Additional array methods
    #  These are not part of the EA API, but we implement them because
    #  pandas assumes they're there.
    # ------------------------------------------------------------------------

    def __array__(self, dtype: NpDtype | None = None) -> np.ndarray:
        return np.asarray(self._array, dtype=dtype)

    def __arrow_array__(self, type=None):
        return self._array.__arrow_array__(type)

    def min(self, axis=None, skipna: bool = True, **kwargs) -> Scalar:
        return self._array.min(axis=axis, skipna=skipna, **kwargs)

    def max(self, axis=None, skipna: bool = True, **kwargs) -> Scalar:
        return self._array.max(axis=axis, skipna=skipna, **kwargs)

    def value_counts(self, dropna: bool = True):
        return self._array.value_counts(dropna=dropna)

    def memory_usage(self, deep: bool = False) -> int:
        return self._array.memory_usage(deep=deep)

    # ------------------------------------------------------------------------
    # OpsMixin interface
    # ------------------------------------------------------------------------

    def _cmp_method(self, other, op):
        return self._array._cmp_method(other, op)

    def _logical_method(self, other, op):
        return self._array._logical_method(other, op)

    def _arith_method(self, other, op):
        return self._array._arith_method(other, op)

    # ------------------------------------------------------------------------
    # String methods interface
    # ------------------------------------------------------------------------

    _str_na_value = StringDtype.na_value

    def _str_map(
        self, f, na_value=None, dtype: Dtype | None = None, convert: bool = True
    ):
        result = self._array._str_map(
            f, na_value=na_value, dtype=dtype, convert=convert
        )
        return self._maybe_wrap_result(result)

    # TODO: dispatch all str accessor methods to array instead of wrapping result of
    # object fallback (_str_map)


class ObjectStringArray(PandasArray):
    # undo the PandasArray hack
    _typ = "extension"

    def __init__(self, values, copy=False):
        values = extract_array(values)

        super().__init__(values, copy=copy)
        # error: Incompatible types in assignment (expression has type "StringDtype",
        # variable has type "PandasDtype")
        NDArrayBacked.__init__(self, self._ndarray, StringDtype())
        if not isinstance(values, type(self)):
            self._validate()

    def _validate(self):
        """Validate that we only store NA or strings."""
        if len(self._ndarray) and not lib.is_string_array(self._ndarray, skipna=True):
            raise ValueError("StringArray requires a sequence of strings or pandas.NA")
        if self._ndarray.dtype != "object":
            raise ValueError(
                "StringArray requires a sequence of strings or pandas.NA. Got "
                f"'{self._ndarray.dtype}' dtype instead."
            )

    @classmethod
    def _from_sequence(cls, scalars, *, dtype: Dtype | None = None, copy=False):
        if dtype:
            assert dtype == "string"

        from pandas.core.arrays.masked import BaseMaskedArray

        if isinstance(scalars, BaseMaskedArray):
            # avoid costly conversion to object dtype
            na_values = scalars._mask
            result = scalars._data
            result = lib.ensure_string_array(result, copy=copy, convert_na_value=False)
            result[na_values] = StringDtype.na_value

        else:
            # convert non-na-likes to str, and nan-likes to StringDtype.na_value
            result = lib.ensure_string_array(
                scalars, na_value=StringDtype.na_value, copy=copy
            )

        # Manually creating new array avoids the validation step in the __init__, so is
        # faster. Refactor need for validation?
        new_string_array = cls.__new__(cls)
        NDArrayBacked.__init__(new_string_array, result, StringDtype())

        return new_string_array

    @classmethod
    def _from_sequence_of_strings(
        cls, strings, *, dtype: Dtype | None = None, copy=False
    ):
        return cls._from_sequence(strings, dtype=dtype, copy=copy)

    @classmethod
    def _empty(cls, shape, dtype) -> ObjectStringArray:
        values = np.empty(shape, dtype=object)
        values[:] = libmissing.NA
        return cls(values).astype(dtype, copy=False)

    def __arrow_array__(self, type=None):
        """
        Convert myself into a pyarrow Array.
        """
        import pyarrow as pa

        if type is None:
            type = pa.string()

        values = self._ndarray.copy()
        values[self.isna()] = None
        return pa.array(values, type=type, from_pandas=True)

    def _values_for_factorize(self):
        arr = self._ndarray.copy()
        mask = self.isna()
        arr[mask] = -1
        return arr, -1

    def __setitem__(self, key, value):
        value = extract_array(value, extract_numpy=True)
        if isinstance(value, type(self)):
            # extract_array doesn't extract PandasArray subclasses
            value = value._ndarray

        key = check_array_indexer(self, key)
        scalar_key = lib.is_scalar(key)
        scalar_value = lib.is_scalar(value)
        if scalar_key and not scalar_value:
            raise ValueError("setting an array element with a sequence.")

        # validate new items
        if scalar_value:
            if isna(value):
                value = StringDtype.na_value
            elif not isinstance(value, str):
                raise ValueError(
                    f"Cannot set non-string value '{value}' into a StringArray."
                )
        else:
            if not is_array_like(value):
                value = np.asarray(value, dtype=object)
            if len(value) and not lib.is_string_array(value, skipna=True):
                raise ValueError("Must provide strings.")

        super().__setitem__(key, value)

    def astype(self, dtype, copy=True):
        dtype = pandas_dtype(dtype)

        if is_dtype_equal(dtype, self.dtype):
            if copy:
                return self.copy()
            return self

        elif isinstance(dtype, _IntegerDtype):
            arr = self._ndarray.copy()
            mask = self.isna()
            arr[mask] = 0
            values = arr.astype(dtype.numpy_dtype)
            return IntegerArray(values, mask, copy=False)
        elif isinstance(dtype, FloatingDtype):
            arr = self.copy()
            mask = self.isna()
            arr[mask] = "0"
            values = arr.astype(dtype.numpy_dtype)
            return FloatingArray(values, mask, copy=False)
        elif isinstance(dtype, ExtensionDtype):
            cls = dtype.construct_array_type()
            return cls._from_sequence(self, dtype=dtype, copy=copy)
        elif np.issubdtype(dtype, np.floating):
            arr = self._ndarray.copy()
            mask = self.isna()
            arr[mask] = 0
            values = arr.astype(dtype)
            values[mask] = np.nan
            return values

        return super().astype(dtype, copy)

    def _reduce(self, name: str, *, skipna: bool = True, **kwargs):
        if name in ["min", "max"]:
            return getattr(self, name)(skipna=skipna)

        raise TypeError(f"Cannot perform reduction '{name}' with string dtype")

    def min(self, axis=None, skipna: bool = True, **kwargs) -> Scalar:
        nv.validate_min((), kwargs)
        result = masked_reductions.min(
            values=self.to_numpy(), mask=self.isna(), skipna=skipna
        )
        return self._wrap_reduction_result(axis, result)

    def max(self, axis=None, skipna: bool = True, **kwargs) -> Scalar:
        nv.validate_max((), kwargs)
        result = masked_reductions.max(
            values=self.to_numpy(), mask=self.isna(), skipna=skipna
        )
        return self._wrap_reduction_result(axis, result)

    def value_counts(self, dropna: bool = True):
        from pandas import value_counts

        return value_counts(self._ndarray, dropna=dropna).astype("Int64")

    def memory_usage(self, deep: bool = False) -> int:
        result = self._ndarray.nbytes
        if deep:
            return result + lib.memory_usage_of_objects(self._ndarray)
        return result

    def _cmp_method(self, other, op):
        from pandas.arrays import BooleanArray

        if isinstance(other, ObjectStringArray):
            other = other._ndarray

        mask = isna(self) | isna(other)
        valid = ~mask

        if not lib.is_scalar(other):
            if len(other) != len(self):
                # prevent improper broadcasting when other is 2D
                raise ValueError(
                    f"Lengths of operands do not match: {len(self)} != {len(other)}"
                )

            other = np.asarray(other)
            other = other[valid]

        if op.__name__ in ops.ARITHMETIC_BINOPS:
            result = np.empty_like(self._ndarray, dtype="object")
            result[mask] = StringDtype.na_value
            result[valid] = op(self._ndarray[valid], other)
            return type(self)(result)
        else:
            # logical
            result = np.zeros(len(self._ndarray), dtype="bool")
            result[valid] = op(self._ndarray[valid], other)
            return BooleanArray(result, mask)

    _arith_method = _cmp_method

    # ------------------------------------------------------------------------
    # String methods interface
    _str_na_value = StringDtype.na_value

    def _str_map(
        self, f, na_value=None, dtype: Dtype | None = None, convert: bool = True
    ):
        from pandas.arrays import BooleanArray

        if dtype is None:
            dtype = StringDtype()
        if na_value is None:
            na_value = self.dtype.na_value

        mask = isna(self)
        arr = np.asarray(self)

        if is_integer_dtype(dtype) or is_bool_dtype(dtype):
            constructor: type[IntegerArray] | type[BooleanArray]
            if is_integer_dtype(dtype):
                constructor = IntegerArray
            else:
                constructor = BooleanArray

            na_value_is_na = isna(na_value)
            if na_value_is_na:
                na_value = 1
            result = lib.map_infer_mask(
                arr,
                f,
                mask.view("uint8"),
                convert=False,
                na_value=na_value,
                # error: Value of type variable "_DTypeScalar" of "dtype" cannot be
                # "object"
                # error: Argument 1 to "dtype" has incompatible type
                # "Union[ExtensionDtype, str, dtype[Any], Type[object]]"; expected
                # "Type[object]"
                dtype=np.dtype(dtype),  # type: ignore[type-var,arg-type]
            )

            if not na_value_is_na:
                mask[:] = False

            return constructor(result, mask)

        elif is_string_dtype(dtype) and not is_object_dtype(dtype):
            # i.e. StringDtype
            result = lib.map_infer_mask(
                arr, f, mask.view("uint8"), convert=False, na_value=na_value
            )
            return type(self)(result)
        else:
            # This is when the result type is object. We reach this when
            # -> We know the result type is truly object (e.g. .encode returns bytes
            #    or .findall returns a list).
            # -> We don't know the result type. E.g. `.get` can return anything.
            return lib.map_infer_mask(arr, f, mask.view("uint8"))
