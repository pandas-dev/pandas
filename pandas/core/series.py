"""
Data structure for 1-dimensional cross-sectional and time series data
"""
from io import StringIO
from shutil import get_terminal_size
from textwrap import dedent
from typing import (
    IO,
    TYPE_CHECKING,
    Any,
    Callable,
    Hashable,
    Iterable,
    List,
    Optional,
    Tuple,
    Type,
)
import warnings

import numpy as np

from pandas._config import get_option

from pandas._libs import lib, properties, reshape, tslibs
from pandas._libs.index import validate_numeric_casting
from pandas._typing import Label
from pandas.compat.numpy import function as nv
from pandas.util._decorators import Appender, Substitution
from pandas.util._validators import validate_bool_kwarg, validate_percentile

from pandas.core.dtypes.cast import convert_dtypes
from pandas.core.dtypes.common import (
    _is_unorderable_exception,
    ensure_platform_int,
    is_bool,
    is_categorical_dtype,
    is_datetime64_dtype,
    is_dict_like,
    is_extension_array_dtype,
    is_integer,
    is_iterator,
    is_list_like,
    is_object_dtype,
    is_scalar,
    is_timedelta64_dtype,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCDatetimeIndex,
    ABCMultiIndex,
    ABCPeriodIndex,
    ABCSeries,
    ABCSparseArray,
)
from pandas.core.dtypes.inference import is_hashable
from pandas.core.dtypes.missing import (
    isna,
    na_value_for_dtype,
    notna,
    remove_na_arraylike,
)

import pandas as pd
from pandas.core import algorithms, base, generic, nanops, ops
from pandas.core.accessor import CachedAccessor
from pandas.core.arrays import ExtensionArray, try_cast_to_ea
from pandas.core.arrays.categorical import Categorical, CategoricalAccessor
from pandas.core.arrays.sparse import SparseAccessor
import pandas.core.common as com
from pandas.core.construction import (
    create_series_with_explicit_dtype,
    extract_array,
    is_empty_data,
    sanitize_array,
)
from pandas.core.indexers import maybe_convert_indices
from pandas.core.indexes.accessors import CombinedDatetimelikeProperties
from pandas.core.indexes.api import (
    Float64Index,
    Index,
    IntervalIndex,
    InvalidIndexError,
    MultiIndex,
    ensure_index,
)
import pandas.core.indexes.base as ibase
from pandas.core.indexes.datetimes import DatetimeIndex
from pandas.core.indexes.period import PeriodIndex
from pandas.core.indexes.timedeltas import TimedeltaIndex
from pandas.core.indexing import check_bool_indexer
from pandas.core.internals import SingleBlockManager
from pandas.core.strings import StringMethods
from pandas.core.tools.datetimes import to_datetime

import pandas.io.formats.format as fmt
import pandas.plotting

if TYPE_CHECKING:
    from pandas.core.frame import DataFrame
    from pandas.core.groupby.generic import SeriesGroupBy

__all__ = ["Series"]

_shared_doc_kwargs = dict(
    axes="index",
    klass="Series",
    axes_single_arg="{0 or 'index'}",
    axis="""axis : {0 or 'index'}
        Parameter needed for compatibility with DataFrame.""",
    inplace="""inplace : boolean, default False
        If True, performs operation inplace and returns None.""",
    unique="np.ndarray",
    duplicated="Series",
    optional_by="",
    optional_mapper="",
    optional_labels="",
    optional_axis="",
    versionadded_to_excel="\n    .. versionadded:: 0.20.0\n",
)


def _coerce_method(converter):
    """
    Install the scalar coercion methods.
    """

    def wrapper(self):
        if len(self) == 1:
            return converter(self.iloc[0])
        raise TypeError(f"cannot convert the series to {converter}")

    wrapper.__name__ = f"__{converter.__name__}__"
    return wrapper


# ----------------------------------------------------------------------
# Series class


class Series(base.IndexOpsMixin, generic.NDFrame):
    """
    One-dimensional ndarray with axis labels (including time series).

    Labels need not be unique but must be a hashable type. The object
    supports both integer- and label-based indexing and provides a host of
    methods for performing operations involving the index. Statistical
    methods from ndarray have been overridden to automatically exclude
    missing data (currently represented as NaN).

    Operations between Series (+, -, /, *, **) align values based on their
    associated index values-- they need not be the same length. The result
    index will be the sorted union of the two indexes.

    Parameters
    ----------
    data : array-like, Iterable, dict, or scalar value
        Contains data stored in Series.

        .. versionchanged:: 0.23.0
           If data is a dict, argument order is maintained for Python 3.6
           and later.

    index : array-like or Index (1d)
        Values must be hashable and have the same length as `data`.
        Non-unique index values are allowed. Will default to
        RangeIndex (0, 1, 2, ..., n) if not provided. If both a dict and index
        sequence are used, the index will override the keys found in the
        dict.
    dtype : str, numpy.dtype, or ExtensionDtype, optional
        Data type for the output Series. If not specified, this will be
        inferred from `data`.
        See the :ref:`user guide <basics.dtypes>` for more usages.
    name : str, optional
        The name to give to the Series.
    copy : bool, default False
        Copy input data.
    """

    _typ = "series"

    _name: Optional[Hashable]
    _metadata: List[str] = ["name"]
    _internal_names_set = {"index"} | generic.NDFrame._internal_names_set
    _accessors = {"dt", "cat", "str", "sparse"}
    _deprecations = (
        base.IndexOpsMixin._deprecations
        | generic.NDFrame._deprecations
        | frozenset(["compress", "ptp"])
    )

    # Override cache_readonly bc Series is mutable
    hasnans = property(
        base.IndexOpsMixin.hasnans.func, doc=base.IndexOpsMixin.hasnans.__doc__
    )
    _data: SingleBlockManager
    div: Callable[["Series", Any], "Series"]
    rdiv: Callable[["Series", Any], "Series"]

    # ----------------------------------------------------------------------
    # Constructors

    def __init__(
        self, data=None, index=None, dtype=None, name=None, copy=False, fastpath=False
    ):

        # we are called internally, so short-circuit
        if fastpath:

            # data is an ndarray, index is defined
            if not isinstance(data, SingleBlockManager):
                data = SingleBlockManager(data, index, fastpath=True)
            if copy:
                data = data.copy()
            if index is None:
                index = data.index

        else:

            name = ibase.maybe_extract_name(name, data, type(self))

            if is_empty_data(data) and dtype is None:
                # gh-17261
                warnings.warn(
                    "The default dtype for empty Series will be 'object' instead "
                    "of 'float64' in a future version. Specify a dtype explicitly "
                    "to silence this warning.",
                    DeprecationWarning,
                    stacklevel=2,
                )
                # uncomment the line below when removing the DeprecationWarning
                # dtype = np.dtype(object)

            if index is not None:
                index = ensure_index(index)

            if data is None:
                data = {}
            if dtype is not None:
                dtype = self._validate_dtype(dtype)

            if isinstance(data, MultiIndex):
                raise NotImplementedError(
                    "initializing a Series from a MultiIndex is not supported"
                )
            elif isinstance(data, Index):

                if dtype is not None:
                    # astype copies
                    data = data.astype(dtype)
                else:
                    # need to copy to avoid aliasing issues
                    data = data._values.copy()
                    if isinstance(data, ABCDatetimeIndex) and data.tz is not None:
                        # GH#24096 need copy to be deep for datetime64tz case
                        # TODO: See if we can avoid these copies
                        data = data._values.copy(deep=True)
                copy = False

            elif isinstance(data, np.ndarray):
                if len(data.dtype):
                    # GH#13296 we are dealing with a compound dtype, which
                    #  should be treated as 2D
                    raise ValueError(
                        "Cannot construct a Series from an ndarray with "
                        "compound dtype.  Use DataFrame instead."
                    )
                pass
            elif isinstance(data, ABCSeries):
                if index is None:
                    index = data.index
                else:
                    data = data.reindex(index, copy=copy)
                data = data._data
            elif is_dict_like(data):
                data, index = self._init_dict(data, index, dtype)
                dtype = None
                copy = False
            elif isinstance(data, SingleBlockManager):
                if index is None:
                    index = data.index
                elif not data.index.equals(index) or copy:
                    # GH#19275 SingleBlockManager input should only be called
                    # internally
                    raise AssertionError(
                        "Cannot pass both SingleBlockManager "
                        "`data` argument and a different "
                        "`index` argument. `copy` must be False."
                    )

            elif is_extension_array_dtype(data):
                pass
            elif isinstance(data, (set, frozenset)):
                raise TypeError(f"'{type(data).__name__}' type is unordered")
            elif isinstance(data, ABCSparseArray):
                # handle sparse passed here (and force conversion)
                data = data.to_dense()
            else:
                data = com.maybe_iterable_to_list(data)

            if index is None:
                if not is_list_like(data):
                    data = [data]
                index = ibase.default_index(len(data))
            elif is_list_like(data):

                # a scalar numpy array is list-like but doesn't
                # have a proper length
                try:
                    if len(index) != len(data):
                        raise ValueError(
                            f"Length of passed values is {len(data)}, "
                            f"index implies {len(index)}."
                        )
                except TypeError:
                    pass

            # create/copy the manager
            if isinstance(data, SingleBlockManager):
                if dtype is not None:
                    data = data.astype(dtype=dtype, errors="ignore", copy=copy)
                elif copy:
                    data = data.copy()
            else:
                data = sanitize_array(data, index, dtype, copy, raise_cast_failure=True)

                data = SingleBlockManager(data, index, fastpath=True)

        generic.NDFrame.__init__(self, data, fastpath=True)
        self.name = name
        self._set_axis(0, index, fastpath=True)

    def _init_dict(self, data, index=None, dtype=None):
        """
        Derive the "_data" and "index" attributes of a new Series from a
        dictionary input.

        Parameters
        ----------
        data : dict or dict-like
            Data used to populate the new Series.
        index : Index or index-like, default None
            Index for the new Series: if None, use dict keys.
        dtype : dtype, default None
            The dtype for the new Series: if None, infer from data.

        Returns
        -------
        _data : BlockManager for the new Series
        index : index for the new Series
        """
        # Looking for NaN in dict doesn't work ({np.nan : 1}[float('nan')]
        # raises KeyError), so we iterate the entire dict, and align
        if data:
            keys, values = zip(*data.items())
            values = list(values)
        elif index is not None:
            # fastpath for Series(data=None). Just use broadcasting a scalar
            # instead of reindexing.
            values = na_value_for_dtype(dtype)
            keys = index
        else:
            keys, values = [], []

        # Input is now list-like, so rely on "standard" construction:

        # TODO: passing np.float64 to not break anything yet. See GH-17261
        s = create_series_with_explicit_dtype(
            values, index=keys, dtype=dtype, dtype_if_empty=np.float64
        )

        # Now we just make sure the order is respected, if any
        if data and index is not None:
            s = s.reindex(index, copy=False)
        return s._data, s.index

    # ----------------------------------------------------------------------

    @property
    def _constructor(self) -> Type["Series"]:
        return Series

    @property
    def _constructor_expanddim(self) -> Type["DataFrame"]:
        from pandas.core.frame import DataFrame

        return DataFrame

    # types
    @property
    def _can_hold_na(self):
        return self._data._can_hold_na

    _index = None

    def _set_axis(self, axis, labels, fastpath: bool = False) -> None:
        """
        Override generic, we want to set the _typ here.
        """

        if not fastpath:
            labels = ensure_index(labels)

        is_all_dates = labels.is_all_dates
        if is_all_dates:
            if not isinstance(labels, (DatetimeIndex, PeriodIndex, TimedeltaIndex)):
                try:
                    labels = DatetimeIndex(labels)
                    # need to set here because we changed the index
                    if fastpath:
                        self._data.set_axis(axis, labels)
                except (tslibs.OutOfBoundsDatetime, ValueError):
                    # labels may exceeds datetime bounds,
                    # or not be a DatetimeIndex
                    pass

        object.__setattr__(self, "_index", labels)
        if not fastpath:
            self._data.set_axis(axis, labels)

    def _update_inplace(self, result, **kwargs):
        # we want to call the generic version and not the IndexOpsMixin
        return generic.NDFrame._update_inplace(self, result, **kwargs)

    # ndarray compatibility
    @property
    def dtype(self):
        """
        Return the dtype object of the underlying data.
        """
        return self._data.dtype

    @property
    def dtypes(self):
        """
        Return the dtype object of the underlying data.
        """
        return self._data.dtype

    @property
    def name(self) -> Optional[Hashable]:
        return self._name

    @name.setter
    def name(self, value: Optional[Hashable]) -> None:
        if not is_hashable(value):
            raise TypeError("Series.name must be a hashable type")
        object.__setattr__(self, "_name", value)

    @property
    def values(self):
        """
        Return Series as ndarray or ndarray-like depending on the dtype.

        .. warning::

           We recommend using :attr:`Series.array` or
           :meth:`Series.to_numpy`, depending on whether you need
           a reference to the underlying data or a NumPy array.

        Returns
        -------
        numpy.ndarray or ndarray-like

        See Also
        --------
        Series.array : Reference to the underlying data.
        Series.to_numpy : A NumPy array representing the underlying data.

        Examples
        --------
        >>> pd.Series([1, 2, 3]).values
        array([1, 2, 3])

        >>> pd.Series(list('aabc')).values
        array(['a', 'a', 'b', 'c'], dtype=object)

        >>> pd.Series(list('aabc')).astype('category').values
        [a, a, b, c]
        Categories (3, object): [a, b, c]

        Timezone aware datetime data is converted to UTC:

        >>> pd.Series(pd.date_range('20130101', periods=3,
        ...                         tz='US/Eastern')).values
        array(['2013-01-01T05:00:00.000000000',
               '2013-01-02T05:00:00.000000000',
               '2013-01-03T05:00:00.000000000'], dtype='datetime64[ns]')
        """
        return self._data.external_values()

    @property
    def _values(self):
        """
        Return the internal repr of this data (defined by Block.interval_values).
        This are the values as stored in the Block (ndarray or ExtensionArray
        depending on the Block class), with datetime64[ns] and timedelta64[ns]
        wrapped in ExtensionArrays to match Index._values behavior.

        Differs from the public ``.values`` for certain data types, because of
        historical backwards compatibility of the public attribute (e.g. period
        returns object ndarray and datetimetz a datetime64[ns] ndarray for
        ``.values`` while it returns an ExtensionArray for ``._values`` in those
        cases).

        Differs from ``.array`` in that this still returns the numpy array if
        the Block is backed by a numpy array (except for datetime64 and
        timedelta64 dtypes), while ``.array`` ensures to always return an
        ExtensionArray.

        Differs from ``._ndarray_values``, as that ensures to always return a
        numpy array (it will call ``_ndarray_values`` on the ExtensionArray, if
        the Series was backed by an ExtensionArray).

        Overview:

        dtype       | values        | _values       | array         | _ndarray_values |
        ----------- | ------------- | ------------- | ------------- | --------------- |
        Numeric     | ndarray       | ndarray       | PandasArray   | ndarray         |
        Category    | Categorical   | Categorical   | Categorical   | ndarray[int]    |
        dt64[ns]    | ndarray[M8ns] | DatetimeArray | DatetimeArray | ndarray[M8ns]   |
        dt64[ns tz] | ndarray[M8ns] | DatetimeArray | DatetimeArray | ndarray[M8ns]   |
        td64[ns]    | ndarray[m8ns] | TimedeltaArray| ndarray[m8ns] | ndarray[m8ns]   |
        Period      | ndarray[obj]  | PeriodArray   | PeriodArray   | ndarray[int]    |
        Nullable    | EA            | EA            | EA            | ndarray         |

        """
        return self._data.internal_values()

    @Appender(base.IndexOpsMixin.array.__doc__)  # type: ignore
    @property
    def array(self) -> ExtensionArray:
        return self._data._block.array_values()

    def _internal_get_values(self):
        """
        Same as values (but handles sparseness conversions); is a view.

        Returns
        -------
        numpy.ndarray
            Data of the Series.
        """

        return self._data.get_values()

    # ops
    def ravel(self, order="C"):
        """
        Return the flattened underlying data as an ndarray.

        Returns
        -------
        numpy.ndarray or ndarray-like
            Flattened data of the Series.

        See Also
        --------
        numpy.ndarray.ravel
        """
        return self._values.ravel(order=order)

    def __len__(self) -> int:
        """
        Return the length of the Series.
        """
        return len(self._data)

    def view(self, dtype=None) -> "Series":
        """
        Create a new view of the Series.

        This function will return a new Series with a view of the same
        underlying values in memory, optionally reinterpreted with a new data
        type. The new data type must preserve the same size in bytes as to not
        cause index misalignment.

        Parameters
        ----------
        dtype : data type
            Data type object or one of their string representations.

        Returns
        -------
        Series
            A new Series object as a view of the same data in memory.

        See Also
        --------
        numpy.ndarray.view : Equivalent numpy function to create a new view of
            the same data in memory.

        Notes
        -----
        Series are instantiated with ``dtype=float64`` by default. While
        ``numpy.ndarray.view()`` will return a view with the same data type as
        the original array, ``Series.view()`` (without specified dtype)
        will try using ``float64`` and may fail if the original data type size
        in bytes is not the same.

        Examples
        --------
        >>> s = pd.Series([-2, -1, 0, 1, 2], dtype='int8')
        >>> s
        0   -2
        1   -1
        2    0
        3    1
        4    2
        dtype: int8

        The 8 bit signed integer representation of `-1` is `0b11111111`, but
        the same bytes represent 255 if read as an 8 bit unsigned integer:

        >>> us = s.view('uint8')
        >>> us
        0    254
        1    255
        2      0
        3      1
        4      2
        dtype: uint8

        The views share the same underlying values:

        >>> us[0] = 128
        >>> s
        0   -128
        1     -1
        2      0
        3      1
        4      2
        dtype: int8
        """
        return self._constructor(
            self._values.view(dtype), index=self.index
        ).__finalize__(self)

    # ----------------------------------------------------------------------
    # NDArray Compat
    _HANDLED_TYPES = (Index, ExtensionArray, np.ndarray)

    def __array_ufunc__(
        self, ufunc: Callable, method: str, *inputs: Any, **kwargs: Any
    ):
        # TODO: handle DataFrame
        cls = type(self)

        # for binary ops, use our custom dunder methods
        result = ops.maybe_dispatch_ufunc_to_dunder_op(
            self, ufunc, method, *inputs, **kwargs
        )
        if result is not NotImplemented:
            return result

        # Determine if we should defer.
        no_defer = (np.ndarray.__array_ufunc__, cls.__array_ufunc__)

        for item in inputs:
            higher_priority = (
                hasattr(item, "__array_priority__")
                and item.__array_priority__ > self.__array_priority__
            )
            has_array_ufunc = (
                hasattr(item, "__array_ufunc__")
                and type(item).__array_ufunc__ not in no_defer
                and not isinstance(item, self._HANDLED_TYPES)
            )
            if higher_priority or has_array_ufunc:
                return NotImplemented

        # align all the inputs.
        names = [getattr(x, "name") for x in inputs if hasattr(x, "name")]
        types = tuple(type(x) for x in inputs)
        # TODO: dataframe
        alignable = [x for x, t in zip(inputs, types) if issubclass(t, Series)]

        if len(alignable) > 1:
            # This triggers alignment.
            # At the moment, there aren't any ufuncs with more than two inputs
            # so this ends up just being x1.index | x2.index, but we write
            # it to handle *args.
            index = alignable[0].index
            for s in alignable[1:]:
                index |= s.index
            inputs = tuple(
                x.reindex(index) if issubclass(t, Series) else x
                for x, t in zip(inputs, types)
            )
        else:
            index = self.index

        inputs = tuple(extract_array(x, extract_numpy=True) for x in inputs)
        result = getattr(ufunc, method)(*inputs, **kwargs)

        name: Optional[Hashable]
        if len(set(names)) == 1:
            name = names[0]
        else:
            name = None

        def construct_return(result):
            if lib.is_scalar(result):
                return result
            elif result.ndim > 1:
                # e.g. np.subtract.outer
                if method == "outer":
                    # GH#27198
                    raise NotImplementedError
                return result
            return self._constructor(result, index=index, name=name, copy=False)

        if type(result) is tuple:
            # multiple return values
            return tuple(construct_return(x) for x in result)
        elif method == "at":
            # no return value
            return None
        else:
            return construct_return(result)

    def __array__(self, dtype=None) -> np.ndarray:
        """
        Return the values as a NumPy array.

        Users should not call this directly. Rather, it is invoked by
        :func:`numpy.array` and :func:`numpy.asarray`.

        Parameters
        ----------
        dtype : str or numpy.dtype, optional
            The dtype to use for the resulting NumPy array. By default,
            the dtype is inferred from the data.

        Returns
        -------
        numpy.ndarray
            The values in the series converted to a :class:`numpy.ndarray`
            with the specified `dtype`.

        See Also
        --------
        array : Create a new array from data.
        Series.array : Zero-copy view to the array backing the Series.
        Series.to_numpy : Series method for similar behavior.

        Examples
        --------
        >>> ser = pd.Series([1, 2, 3])
        >>> np.asarray(ser)
        array([1, 2, 3])

        For timezone-aware data, the timezones may be retained with
        ``dtype='object'``

        >>> tzser = pd.Series(pd.date_range('2000', periods=2, tz="CET"))
        >>> np.asarray(tzser, dtype="object")
        array([Timestamp('2000-01-01 00:00:00+0100', tz='CET', freq='D'),
               Timestamp('2000-01-02 00:00:00+0100', tz='CET', freq='D')],
              dtype=object)

        Or the values may be localized to UTC and the tzinfo discarded with
        ``dtype='datetime64[ns]'``

        >>> np.asarray(tzser, dtype="datetime64[ns]")  # doctest: +ELLIPSIS
        array(['1999-12-31T23:00:00.000000000', ...],
              dtype='datetime64[ns]')
        """
        return np.asarray(self.array, dtype)

    # ----------------------------------------------------------------------
    # Unary Methods

    # coercion
    __float__ = _coerce_method(float)
    __long__ = _coerce_method(int)
    __int__ = _coerce_method(int)

    # ----------------------------------------------------------------------

    # indexers
    @property
    def axes(self) -> List[Index]:
        """
        Return a list of the row axis labels.
        """
        return [self.index]

    # ----------------------------------------------------------------------
    # Indexing Methods

    @Appender(generic.NDFrame.take.__doc__)
    def take(self, indices, axis=0, is_copy=None, **kwargs) -> "Series":
        if is_copy is not None:
            warnings.warn(
                "is_copy is deprecated and will be removed in a future version. "
                "'take' always returns a copy, so there is no need to specify this.",
                FutureWarning,
                stacklevel=2,
            )
        nv.validate_take(tuple(), kwargs)

        indices = ensure_platform_int(indices)
        new_index = self.index.take(indices)

        if is_categorical_dtype(self):
            # https://github.com/pandas-dev/pandas/issues/20664
            # TODO: remove when the default Categorical.take behavior changes
            indices = maybe_convert_indices(indices, len(self._get_axis(axis)))
            kwargs = {"allow_fill": False}
        else:
            kwargs = {}
        new_values = self._values.take(indices, **kwargs)

        return self._constructor(
            new_values, index=new_index, fastpath=True
        ).__finalize__(self)

    def _take_with_is_copy(self, indices, axis=0, **kwargs):
        """
        Internal version of the `take` method that sets the `_is_copy`
        attribute to keep track of the parent dataframe (using in indexing
        for the SettingWithCopyWarning). For Series this does the same
        as the public take (it never sets `_is_copy`).

        See the docstring of `take` for full explanation of the parameters.
        """
        return self.take(indices=indices, axis=axis, **kwargs)

    def _ixs(self, i: int, axis: int = 0):
        """
        Return the i-th value or values in the Series by location.

        Parameters
        ----------
        i : int

        Returns
        -------
        scalar (int) or Series (slice, sequence)
        """
        return self._values[i]

    def _slice(self, slobj: slice, axis: int = 0, kind=None) -> "Series":
        slobj = self.index._convert_slice_indexer(slobj, kind=kind or "getitem")
        return self._get_values(slobj)

    def __getitem__(self, key):
        key = com.apply_if_callable(key, self)

        if key is Ellipsis:
            return self

        key_is_scalar = is_scalar(key)

        if key_is_scalar or isinstance(self.index, MultiIndex):
            # Otherwise index.get_value will raise InvalidIndexError
            try:
                result = self.index.get_value(self, key)

                return result
            except InvalidIndexError:
                pass
            except (KeyError, ValueError):
                if isinstance(key, tuple) and isinstance(self.index, MultiIndex):
                    # kludge
                    pass
                else:

                    # we can try to coerce the indexer (or this will raise)
                    new_key = self.index._convert_scalar_indexer(key, kind="getitem")
                    if type(new_key) != type(key):
                        return self.__getitem__(new_key)
                    raise

        if not key_is_scalar:
            # avoid expensive checks if we know we have a scalar
            if is_iterator(key):
                key = list(key)

            if com.is_bool_indexer(key):
                key = check_bool_indexer(self.index, key)
                return self._get_values(key)

        return self._get_with(key)

    def _get_with(self, key):
        # other: fancy integer or otherwise
        if isinstance(key, slice):
            return self._slice(key)
        elif isinstance(key, ABCDataFrame):
            raise TypeError(
                "Indexing a Series with DataFrame is not "
                "supported, use the appropriate DataFrame column"
            )
        elif isinstance(key, tuple):
            try:
                return self._get_values_tuple(key)
            except ValueError:
                # if we don't have a MultiIndex, we may still be able to handle
                #  a 1-tuple.  see test_1tuple_without_multiindex
                if len(key) == 1:
                    key = key[0]
                    if isinstance(key, slice):
                        return self._get_values(key)
                raise

        if not isinstance(key, (list, np.ndarray, Series, Index)):
            key = list(key)

        if isinstance(key, Index):
            key_type = key.inferred_type
        else:
            key_type = lib.infer_dtype(key, skipna=False)

        # Note: The key_type == "boolean" case should be caught by the
        #  com.is_bool_indexer check in __getitem__
        if key_type == "integer":
            if self.index.is_integer() or self.index.is_floating():
                return self.loc[key]
            elif isinstance(self.index, IntervalIndex):
                indexer = self.index.get_indexer_for(key)
                return self.iloc[indexer]
            else:
                return self._get_values(key)

        if isinstance(key, (list, tuple)):
            # TODO: de-dup with tuple case handled above?
            # handle the dup indexing case GH#4246
            if len(key) == 1 and isinstance(key[0], slice):
                # [slice(0, 5, None)] will break if you convert to ndarray,
                # e.g. as requested by np.median
                # FIXME: hack
                return self._get_values(key)

            return self.loc[key]

        return self.reindex(key)

    def _get_values_tuple(self, key):
        # mpl hackaround
        if com.any_none(*key):
            # suppress warning from slicing the index with a 2d indexer.
            # eventually we'll want Series itself to warn.
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore", "Support for multi-dim", DeprecationWarning
                )
                return self._get_values(key)

        if not isinstance(self.index, MultiIndex):
            raise ValueError("Can only tuple-index with a MultiIndex")

        # If key is contained, would have returned by now
        indexer, new_index = self.index.get_loc_level(key)
        return self._constructor(self._values[indexer], index=new_index).__finalize__(
            self
        )

    def _get_values(self, indexer):
        try:
            return self._constructor(
                self._data.get_slice(indexer), fastpath=True
            ).__finalize__(self)
        except ValueError:
            # mpl compat if we look up e.g. ser[:, np.newaxis];
            #  see tests.series.timeseries.test_mpl_compat_hack
            return self._values[indexer]

    def _get_value(self, label, takeable: bool = False):
        """
        Quickly retrieve single value at passed index label.

        Parameters
        ----------
        label : object
        takeable : interpret the index as indexers, default False

        Returns
        -------
        scalar value
        """
        if takeable:
            return self._values[label]
        return self.index.get_value(self, label)

    def __setitem__(self, key, value):
        key = com.apply_if_callable(key, self)
        cacher_needs_updating = self._check_is_chained_assignment_possible()

        try:
            self._set_with_engine(key, value)
        except com.SettingWithCopyError:
            raise
        except (KeyError, ValueError):
            values = self._values
            if is_integer(key) and not self.index.inferred_type == "integer":
                values[key] = value
            elif key is Ellipsis:
                self[:] = value
            else:
                self.loc[key] = value
        except InvalidIndexError:
            # e.g. slice
            self._set_with(key, value)

        except TypeError as e:
            if isinstance(key, tuple) and not isinstance(self.index, MultiIndex):
                raise ValueError("Can only tuple-index with a MultiIndex")

            # python 3 type errors should be raised
            if _is_unorderable_exception(e):
                raise IndexError(key)

            if com.is_bool_indexer(key):
                key = check_bool_indexer(self.index, key)
                try:
                    self._where(~key, value, inplace=True)
                    return
                except InvalidIndexError:
                    pass

            self._set_with(key, value)

        if cacher_needs_updating:
            self._maybe_update_cacher()

    def _set_with_engine(self, key, value):
        # fails with AttributeError for IntervalIndex
        loc = self.index._engine.get_loc(key)
        validate_numeric_casting(self.dtype, value)
        self._values[loc] = value

    def _set_with(self, key, value):
        # other: fancy integer or otherwise
        if isinstance(key, slice):
            indexer = self.index._convert_slice_indexer(key, kind="getitem")
            return self._set_values(indexer, value)

        elif is_scalar(key) and not is_integer(key) and key not in self.index:
            # GH#12862 adding an new key to the Series
            # Note: have to exclude integers because that is ambiguously
            #  position-based
            self.loc[key] = value
            return

        else:
            if isinstance(key, tuple):
                try:
                    # TODO: no test cases that get here
                    self._set_values(key, value)
                except Exception:
                    pass

            if is_scalar(key):
                key = [key]

            if isinstance(key, Index):
                key_type = key.inferred_type
                key = key._values
            else:
                key_type = lib.infer_dtype(key, skipna=False)

            if key_type == "integer":
                if self.index.inferred_type == "integer":
                    self._set_labels(key, value)
                else:
                    return self._set_values(key, value)
            elif key_type == "boolean":
                self._set_values(key.astype(np.bool_), value)
            else:
                self._set_labels(key, value)

    def _set_labels(self, key, value):
        key = com.asarray_tuplesafe(key)
        indexer = self.index.get_indexer(key)
        mask = indexer == -1
        if mask.any():
            raise ValueError(f"{key[mask]} not contained in the index")
        self._set_values(indexer, value)

    def _set_values(self, key, value):
        if isinstance(key, Series):
            key = key._values
        self._data = self._data.setitem(indexer=key, value=value)
        self._maybe_update_cacher()

    def _set_value(self, label, value, takeable: bool = False):
        """
        Quickly set single value at passed label.

        If label is not contained, a new object is created with the label
        placed at the end of the result index.

        Parameters
        ----------
        label : object
            Partial indexing with MultiIndex not allowed.
        value : object
            Scalar value.
        takeable : interpret the index as indexers, default False

        Returns
        -------
        Series
            If label is contained, will be reference to calling Series,
            otherwise a new object.
        """
        try:
            if takeable:
                self._values[label] = value
            else:
                loc = self.index.get_loc(label)
                validate_numeric_casting(self.dtype, value)
                self._values[loc] = value
        except KeyError:

            # set using a non-recursive method
            self.loc[label] = value

        return self

    # ----------------------------------------------------------------------
    # Unsorted

    @property
    def _is_mixed_type(self):
        return False

    def repeat(self, repeats, axis=None) -> "Series":
        """
        Repeat elements of a Series.

        Returns a new Series where each element of the current Series
        is repeated consecutively a given number of times.

        Parameters
        ----------
        repeats : int or array of ints
            The number of repetitions for each element. This should be a
            non-negative integer. Repeating 0 times will return an empty
            Series.
        axis : None
            Must be ``None``. Has no effect but is accepted for compatibility
            with numpy.

        Returns
        -------
        Series
            Newly created Series with repeated elements.

        See Also
        --------
        Index.repeat : Equivalent function for Index.
        numpy.repeat : Similar method for :class:`numpy.ndarray`.

        Examples
        --------
        >>> s = pd.Series(['a', 'b', 'c'])
        >>> s
        0    a
        1    b
        2    c
        dtype: object
        >>> s.repeat(2)
        0    a
        0    a
        1    b
        1    b
        2    c
        2    c
        dtype: object
        >>> s.repeat([1, 2, 3])
        0    a
        1    b
        1    b
        2    c
        2    c
        2    c
        dtype: object
        """
        nv.validate_repeat(tuple(), dict(axis=axis))
        new_index = self.index.repeat(repeats)
        new_values = self._values.repeat(repeats)
        return self._constructor(new_values, index=new_index).__finalize__(self)

    def reset_index(self, level=None, drop=False, name=None, inplace=False):
        """
        Generate a new DataFrame or Series with the index reset.

        This is useful when the index needs to be treated as a column, or
        when the index is meaningless and needs to be reset to the default
        before another operation.

        Parameters
        ----------
        level : int, str, tuple, or list, default optional
            For a Series with a MultiIndex, only remove the specified levels
            from the index. Removes all levels by default.
        drop : bool, default False
            Just reset the index, without inserting it as a column in
            the new DataFrame.
        name : object, optional
            The name to use for the column containing the original Series
            values. Uses ``self.name`` by default. This argument is ignored
            when `drop` is True.
        inplace : bool, default False
            Modify the Series in place (do not create a new object).

        Returns
        -------
        Series or DataFrame
            When `drop` is False (the default), a DataFrame is returned.
            The newly created columns will come first in the DataFrame,
            followed by the original Series values.
            When `drop` is True, a `Series` is returned.
            In either case, if ``inplace=True``, no value is returned.

        See Also
        --------
        DataFrame.reset_index: Analogous function for DataFrame.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4], name='foo',
        ...               index=pd.Index(['a', 'b', 'c', 'd'], name='idx'))

        Generate a DataFrame with default index.

        >>> s.reset_index()
          idx  foo
        0   a    1
        1   b    2
        2   c    3
        3   d    4

        To specify the name of the new column use `name`.

        >>> s.reset_index(name='values')
          idx  values
        0   a       1
        1   b       2
        2   c       3
        3   d       4

        To generate a new Series with the default set `drop` to True.

        >>> s.reset_index(drop=True)
        0    1
        1    2
        2    3
        3    4
        Name: foo, dtype: int64

        To update the Series in place, without generating a new one
        set `inplace` to True. Note that it also requires ``drop=True``.

        >>> s.reset_index(inplace=True, drop=True)
        >>> s
        0    1
        1    2
        2    3
        3    4
        Name: foo, dtype: int64

        The `level` parameter is interesting for Series with a multi-level
        index.

        >>> arrays = [np.array(['bar', 'bar', 'baz', 'baz']),
        ...           np.array(['one', 'two', 'one', 'two'])]
        >>> s2 = pd.Series(
        ...     range(4), name='foo',
        ...     index=pd.MultiIndex.from_arrays(arrays,
        ...                                     names=['a', 'b']))

        To remove a specific level from the Index, use `level`.

        >>> s2.reset_index(level='a')
               a  foo
        b
        one  bar    0
        two  bar    1
        one  baz    2
        two  baz    3

        If `level` is not set, all levels are removed from the Index.

        >>> s2.reset_index()
             a    b  foo
        0  bar  one    0
        1  bar  two    1
        2  baz  one    2
        3  baz  two    3
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        if drop:
            new_index = ibase.default_index(len(self))
            if level is not None:
                if not isinstance(level, (tuple, list)):
                    level = [level]
                level = [self.index._get_level_number(lev) for lev in level]
                if len(level) < self.index.nlevels:
                    new_index = self.index.droplevel(level)

            if inplace:
                self.index = new_index
                # set name if it was passed, otherwise, keep the previous name
                self.name = name or self.name
            else:
                return self._constructor(
                    self._values.copy(), index=new_index
                ).__finalize__(self)
        elif inplace:
            raise TypeError(
                "Cannot reset_index inplace on a Series to create a DataFrame"
            )
        else:
            df = self.to_frame(name)
            return df.reset_index(level=level, drop=drop)

    # ----------------------------------------------------------------------
    # Rendering Methods

    def __repr__(self) -> str:
        """
        Return a string representation for a particular Series.
        """
        buf = StringIO("")
        width, height = get_terminal_size()
        max_rows = (
            height
            if get_option("display.max_rows") == 0
            else get_option("display.max_rows")
        )
        min_rows = (
            height
            if get_option("display.max_rows") == 0
            else get_option("display.min_rows")
        )
        show_dimensions = get_option("display.show_dimensions")

        self.to_string(
            buf=buf,
            name=self.name,
            dtype=self.dtype,
            min_rows=min_rows,
            max_rows=max_rows,
            length=show_dimensions,
        )
        result = buf.getvalue()

        return result

    def to_string(
        self,
        buf=None,
        na_rep="NaN",
        float_format=None,
        header=True,
        index=True,
        length=False,
        dtype=False,
        name=False,
        max_rows=None,
        min_rows=None,
    ):
        """
        Render a string representation of the Series.

        Parameters
        ----------
        buf : StringIO-like, optional
            Buffer to write to.
        na_rep : str, optional
            String representation of NaN to use, default 'NaN'.
        float_format : one-parameter function, optional
            Formatter function to apply to columns' elements if they are
            floats, default None.
        header : bool, default True
            Add the Series header (index name).
        index : bool, optional
            Add index (row) labels, default True.
        length : bool, default False
            Add the Series length.
        dtype : bool, default False
            Add the Series dtype.
        name : bool, default False
            Add the Series name if not None.
        max_rows : int, optional
            Maximum number of rows to show before truncating. If None, show
            all.
        min_rows : int, optional
            The number of rows to display in a truncated repr (when number
            of rows is above `max_rows`).

        Returns
        -------
        str or None
            String representation of Series if ``buf=None``, otherwise None.
        """

        formatter = fmt.SeriesFormatter(
            self,
            name=name,
            length=length,
            header=header,
            index=index,
            dtype=dtype,
            na_rep=na_rep,
            float_format=float_format,
            min_rows=min_rows,
            max_rows=max_rows,
        )
        result = formatter.to_string()

        # catch contract violations
        if not isinstance(result, str):
            raise AssertionError(
                "result must be of type str, type "
                f"of result is {repr(type(result).__name__)}"
            )

        if buf is None:
            return result
        else:
            try:
                buf.write(result)
            except AttributeError:
                with open(buf, "w") as f:
                    f.write(result)

    @Appender(
        """
        Examples
        --------
        >>> s = pd.Series(["elk", "pig", "dog", "quetzal"], name="animal")
        >>> print(s.to_markdown())
        |    | animal   |
        |---:|:---------|
        |  0 | elk      |
        |  1 | pig      |
        |  2 | dog      |
        |  3 | quetzal  |
        """
    )
    @Substitution(klass="Series")
    @Appender(generic._shared_docs["to_markdown"])
    def to_markdown(
        self, buf: Optional[IO[str]] = None, mode: Optional[str] = None, **kwargs
    ) -> Optional[str]:
        return self.to_frame().to_markdown(buf, mode, **kwargs)

    # ----------------------------------------------------------------------

    def items(self) -> Iterable[Tuple[Label, Any]]:
        """
        Lazily iterate over (index, value) tuples.

        This method returns an iterable tuple (index, value). This is
        convenient if you want to create a lazy iterator.

        Returns
        -------
        iterable
            Iterable of tuples containing the (index, value) pairs from a
            Series.

        See Also
        --------
        DataFrame.items : Iterate over (column name, Series) pairs.
        DataFrame.iterrows : Iterate over DataFrame rows as (index, Series) pairs.

        Examples
        --------
        >>> s = pd.Series(['A', 'B', 'C'])
        >>> for index, value in s.items():
        ...     print(f"Index : {index}, Value : {value}")
        Index : 0, Value : A
        Index : 1, Value : B
        Index : 2, Value : C
        """
        return zip(iter(self.index), iter(self))

    @Appender(items.__doc__)
    def iteritems(self) -> Iterable[Tuple[Label, Any]]:
        return self.items()

    # ----------------------------------------------------------------------
    # Misc public methods

    def keys(self) -> Index:
        """
        Return alias for index.

        Returns
        -------
        Index
            Index of the Series.
        """
        return self.index

    def to_dict(self, into=dict):
        """
        Convert Series to {label -> value} dict or dict-like object.

        Parameters
        ----------
        into : class, default dict
            The collections.abc.Mapping subclass to use as the return
            object. Can be the actual class or an empty
            instance of the mapping type you want.  If you want a
            collections.defaultdict, you must pass it initialized.

            .. versionadded:: 0.21.0

        Returns
        -------
        collections.abc.Mapping
            Key-value representation of Series.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4])
        >>> s.to_dict()
        {0: 1, 1: 2, 2: 3, 3: 4}
        >>> from collections import OrderedDict, defaultdict
        >>> s.to_dict(OrderedDict)
        OrderedDict([(0, 1), (1, 2), (2, 3), (3, 4)])
        >>> dd = defaultdict(list)
        >>> s.to_dict(dd)
        defaultdict(<class 'list'>, {0: 1, 1: 2, 2: 3, 3: 4})
        """
        # GH16122
        into_c = com.standardize_mapping(into)
        return into_c(self.items())

    def to_frame(self, name=None) -> "DataFrame":
        """
        Convert Series to DataFrame.

        Parameters
        ----------
        name : object, default None
            The passed name should substitute for the series name (if it has
            one).

        Returns
        -------
        DataFrame
            DataFrame representation of Series.

        Examples
        --------
        >>> s = pd.Series(["a", "b", "c"],
        ...               name="vals")
        >>> s.to_frame()
          vals
        0    a
        1    b
        2    c
        """
        if name is None:
            df = self._constructor_expanddim(self)
        else:
            df = self._constructor_expanddim({name: self})

        return df

    def _set_name(self, name, inplace=False) -> "Series":
        """
        Set the Series name.

        Parameters
        ----------
        name : str
        inplace : bool
            Whether to modify `self` directly or return a copy.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        ser = self if inplace else self.copy()
        ser.name = name
        return ser

    @Appender(
        """
Examples
--------
>>> ser = pd.Series([390., 350., 30., 20.],
...                 index=['Falcon', 'Falcon', 'Parrot', 'Parrot'], name="Max Speed")
>>> ser
Falcon    390.0
Falcon    350.0
Parrot     30.0
Parrot     20.0
Name: Max Speed, dtype: float64
>>> ser.groupby(["a", "b", "a", "b"]).mean()
a    210.0
b    185.0
Name: Max Speed, dtype: float64
>>> ser.groupby(level=0).mean()
Falcon    370.0
Parrot     25.0
Name: Max Speed, dtype: float64
>>> ser.groupby(ser > 100).mean()
Max Speed
False     25.0
True     370.0
Name: Max Speed, dtype: float64

**Grouping by Indexes**

We can groupby different levels of a hierarchical index
using the `level` parameter:

>>> arrays = [['Falcon', 'Falcon', 'Parrot', 'Parrot'],
...           ['Captive', 'Wild', 'Captive', 'Wild']]
>>> index = pd.MultiIndex.from_arrays(arrays, names=('Animal', 'Type'))
>>> ser = pd.Series([390., 350., 30., 20.], index=index, name="Max Speed")
>>> ser
Animal  Type
Falcon  Captive    390.0
        Wild       350.0
Parrot  Captive     30.0
        Wild        20.0
Name: Max Speed, dtype: float64
>>> ser.groupby(level=0).mean()
Animal
Falcon    370.0
Parrot     25.0
Name: Max Speed, dtype: float64
>>> ser.groupby(level="Type").mean()
Type
Captive    210.0
Wild       185.0
Name: Max Speed, dtype: float64
"""
    )
    @Appender(generic._shared_docs["groupby"] % _shared_doc_kwargs)
    def groupby(
        self,
        by=None,
        axis=0,
        level=None,
        as_index: bool = True,
        sort: bool = True,
        group_keys: bool = True,
        squeeze: bool = False,
        observed: bool = False,
    ) -> "SeriesGroupBy":
        from pandas.core.groupby.generic import SeriesGroupBy

        if level is None and by is None:
            raise TypeError("You have to supply one of 'by' and 'level'")
        axis = self._get_axis_number(axis)

        return SeriesGroupBy(
            obj=self,
            keys=by,
            axis=axis,
            level=level,
            as_index=as_index,
            sort=sort,
            group_keys=group_keys,
            squeeze=squeeze,
            observed=observed,
        )

    # ----------------------------------------------------------------------
    # Statistics, overridden ndarray methods

    # TODO: integrate bottleneck

    def count(self, level=None):
        """
        Return number of non-NA/null observations in the Series.

        Parameters
        ----------
        level : int or level name, default None
            If the axis is a MultiIndex (hierarchical), count along a
            particular level, collapsing into a smaller Series.

        Returns
        -------
        int or Series (if level specified)
            Number of non-null values in the Series.

        Examples
        --------
        >>> s = pd.Series([0.0, 1.0, np.nan])
        >>> s.count()
        2
        """
        if level is None:
            return notna(self.array).sum()

        if isinstance(level, str):
            level = self.index._get_level_number(level)

        lev = self.index.levels[level]
        level_codes = np.array(self.index.codes[level], subok=False, copy=True)

        mask = level_codes == -1
        if mask.any():
            level_codes[mask] = cnt = len(lev)
            lev = lev.insert(cnt, lev._na_value)

        obs = level_codes[notna(self.values)]
        out = np.bincount(obs, minlength=len(lev) or None)
        return self._constructor(out, index=lev, dtype="int64").__finalize__(self)

    def mode(self, dropna=True) -> "Series":
        """
        Return the mode(s) of the dataset.

        Always returns Series even if only one value is returned.

        Parameters
        ----------
        dropna : bool, default True
            Don't consider counts of NaN/NaT.

            .. versionadded:: 0.24.0

        Returns
        -------
        Series
            Modes of the Series in sorted order.
        """
        # TODO: Add option for bins like value_counts()
        return algorithms.mode(self, dropna=dropna)

    def unique(self):
        """
        Return unique values of Series object.

        Uniques are returned in order of appearance. Hash table-based unique,
        therefore does NOT sort.

        Returns
        -------
        ndarray or ExtensionArray
            The unique values returned as a NumPy array. See Notes.

        See Also
        --------
        unique : Top-level unique method for any 1-d array-like object.
        Index.unique : Return Index with unique values from an Index object.

        Notes
        -----
        Returns the unique values as a NumPy array. In case of an
        extension-array backed Series, a new
        :class:`~api.extensions.ExtensionArray` of that type with just
        the unique values is returned. This includes

            * Categorical
            * Period
            * Datetime with Timezone
            * Interval
            * Sparse
            * IntegerNA

        See Examples section.

        Examples
        --------
        >>> pd.Series([2, 1, 3, 3], name='A').unique()
        array([2, 1, 3])

        >>> pd.Series([pd.Timestamp('2016-01-01') for _ in range(3)]).unique()
        array(['2016-01-01T00:00:00.000000000'], dtype='datetime64[ns]')

        >>> pd.Series([pd.Timestamp('2016-01-01', tz='US/Eastern')
        ...            for _ in range(3)]).unique()
        <DatetimeArray>
        ['2016-01-01 00:00:00-05:00']
        Length: 1, dtype: datetime64[ns, US/Eastern]

        An unordered Categorical will return categories in the order of
        appearance.

        >>> pd.Series(pd.Categorical(list('baabc'))).unique()
        [b, a, c]
        Categories (3, object): [b, a, c]

        An ordered Categorical preserves the category ordering.

        >>> pd.Series(pd.Categorical(list('baabc'), categories=list('abc'),
        ...                          ordered=True)).unique()
        [b, a, c]
        Categories (3, object): [a < b < c]
        """
        result = super().unique()
        return result

    def drop_duplicates(self, keep="first", inplace=False) -> "Series":
        """
        Return Series with duplicate values removed.

        Parameters
        ----------
        keep : {'first', 'last', ``False``}, default 'first'
            Method to handle dropping duplicates:

            - 'first' : Drop duplicates except for the first occurrence.
            - 'last' : Drop duplicates except for the last occurrence.
            - ``False`` : Drop all duplicates.

        inplace : bool, default ``False``
            If ``True``, performs operation inplace and returns None.

        Returns
        -------
        Series
            Series with duplicates dropped.

        See Also
        --------
        Index.drop_duplicates : Equivalent method on Index.
        DataFrame.drop_duplicates : Equivalent method on DataFrame.
        Series.duplicated : Related method on Series, indicating duplicate
            Series values.

        Examples
        --------
        Generate a Series with duplicated entries.

        >>> s = pd.Series(['lama', 'cow', 'lama', 'beetle', 'lama', 'hippo'],
        ...               name='animal')
        >>> s
        0      lama
        1       cow
        2      lama
        3    beetle
        4      lama
        5     hippo
        Name: animal, dtype: object

        With the 'keep' parameter, the selection behaviour of duplicated values
        can be changed. The value 'first' keeps the first occurrence for each
        set of duplicated entries. The default value of keep is 'first'.

        >>> s.drop_duplicates()
        0      lama
        1       cow
        3    beetle
        5     hippo
        Name: animal, dtype: object

        The value 'last' for parameter 'keep' keeps the last occurrence for
        each set of duplicated entries.

        >>> s.drop_duplicates(keep='last')
        1       cow
        3    beetle
        4      lama
        5     hippo
        Name: animal, dtype: object

        The value ``False`` for parameter 'keep' discards all sets of
        duplicated entries. Setting the value of 'inplace' to ``True`` performs
        the operation inplace and returns ``None``.

        >>> s.drop_duplicates(keep=False, inplace=True)
        >>> s
        1       cow
        3    beetle
        5     hippo
        Name: animal, dtype: object
        """
        return super().drop_duplicates(keep=keep, inplace=inplace)

    def duplicated(self, keep="first") -> "Series":
        """
        Indicate duplicate Series values.

        Duplicated values are indicated as ``True`` values in the resulting
        Series. Either all duplicates, all except the first or all except the
        last occurrence of duplicates can be indicated.

        Parameters
        ----------
        keep : {'first', 'last', False}, default 'first'
            Method to handle dropping duplicates:

            - 'first' : Mark duplicates as ``True`` except for the first
              occurrence.
            - 'last' : Mark duplicates as ``True`` except for the last
              occurrence.
            - ``False`` : Mark all duplicates as ``True``.

        Returns
        -------
        Series
            Series indicating whether each value has occurred in the
            preceding values.

        See Also
        --------
        Index.duplicated : Equivalent method on pandas.Index.
        DataFrame.duplicated : Equivalent method on pandas.DataFrame.
        Series.drop_duplicates : Remove duplicate values from Series.

        Examples
        --------
        By default, for each set of duplicated values, the first occurrence is
        set on False and all others on True:

        >>> animals = pd.Series(['lama', 'cow', 'lama', 'beetle', 'lama'])
        >>> animals.duplicated()
        0    False
        1    False
        2     True
        3    False
        4     True
        dtype: bool

        which is equivalent to

        >>> animals.duplicated(keep='first')
        0    False
        1    False
        2     True
        3    False
        4     True
        dtype: bool

        By using 'last', the last occurrence of each set of duplicated values
        is set on False and all others on True:

        >>> animals.duplicated(keep='last')
        0     True
        1    False
        2     True
        3    False
        4    False
        dtype: bool

        By setting keep on ``False``, all duplicates are True:

        >>> animals.duplicated(keep=False)
        0     True
        1    False
        2     True
        3    False
        4     True
        dtype: bool
        """
        return super().duplicated(keep=keep)

    def idxmin(self, axis=0, skipna=True, *args, **kwargs):
        """
        Return the row label of the minimum value.

        If multiple values equal the minimum, the first row label with that
        value is returned.

        Parameters
        ----------
        axis : int, default 0
            For compatibility with DataFrame.idxmin. Redundant for application
            on Series.
        skipna : bool, default True
            Exclude NA/null values. If the entire Series is NA, the result
            will be NA.
        *args, **kwargs
            Additional arguments and keywords have no effect but might be
            accepted for compatibility with NumPy.

        Returns
        -------
        Index
            Label of the minimum value.

        Raises
        ------
        ValueError
            If the Series is empty.

        See Also
        --------
        numpy.argmin : Return indices of the minimum values
            along the given axis.
        DataFrame.idxmin : Return index of first occurrence of minimum
            over requested axis.
        Series.idxmax : Return index *label* of the first occurrence
            of maximum of values.

        Notes
        -----
        This method is the Series version of ``ndarray.argmin``. This method
        returns the label of the minimum, while ``ndarray.argmin`` returns
        the position. To get the position, use ``series.values.argmin()``.

        Examples
        --------
        >>> s = pd.Series(data=[1, None, 4, 1],
        ...               index=['A', 'B', 'C', 'D'])
        >>> s
        A    1.0
        B    NaN
        C    4.0
        D    1.0
        dtype: float64

        >>> s.idxmin()
        'A'

        If `skipna` is False and there is an NA value in the data,
        the function returns ``nan``.

        >>> s.idxmin(skipna=False)
        nan
        """
        skipna = nv.validate_argmin_with_skipna(skipna, args, kwargs)
        i = nanops.nanargmin(com.values_from_object(self), skipna=skipna)
        if i == -1:
            return np.nan
        return self.index[i]

    def idxmax(self, axis=0, skipna=True, *args, **kwargs):
        """
        Return the row label of the maximum value.

        If multiple values equal the maximum, the first row label with that
        value is returned.

        Parameters
        ----------
        axis : int, default 0
            For compatibility with DataFrame.idxmax. Redundant for application
            on Series.
        skipna : bool, default True
            Exclude NA/null values. If the entire Series is NA, the result
            will be NA.
        *args, **kwargs
            Additional arguments and keywords have no effect but might be
            accepted for compatibility with NumPy.

        Returns
        -------
        Index
            Label of the maximum value.

        Raises
        ------
        ValueError
            If the Series is empty.

        See Also
        --------
        numpy.argmax : Return indices of the maximum values
            along the given axis.
        DataFrame.idxmax : Return index of first occurrence of maximum
            over requested axis.
        Series.idxmin : Return index *label* of the first occurrence
            of minimum of values.

        Notes
        -----
        This method is the Series version of ``ndarray.argmax``. This method
        returns the label of the maximum, while ``ndarray.argmax`` returns
        the position. To get the position, use ``series.values.argmax()``.

        Examples
        --------
        >>> s = pd.Series(data=[1, None, 4, 3, 4],
        ...               index=['A', 'B', 'C', 'D', 'E'])
        >>> s
        A    1.0
        B    NaN
        C    4.0
        D    3.0
        E    4.0
        dtype: float64

        >>> s.idxmax()
        'C'

        If `skipna` is False and there is an NA value in the data,
        the function returns ``nan``.

        >>> s.idxmax(skipna=False)
        nan
        """
        skipna = nv.validate_argmax_with_skipna(skipna, args, kwargs)
        i = nanops.nanargmax(com.values_from_object(self), skipna=skipna)
        if i == -1:
            return np.nan
        return self.index[i]

    def round(self, decimals=0, *args, **kwargs) -> "Series":
        """
        Round each value in a Series to the given number of decimals.

        Parameters
        ----------
        decimals : int, default 0
            Number of decimal places to round to. If decimals is negative,
            it specifies the number of positions to the left of the decimal point.

        Returns
        -------
        Series
            Rounded values of the Series.

        See Also
        --------
        numpy.around : Round values of an np.array.
        DataFrame.round : Round values of a DataFrame.

        Examples
        --------
        >>> s = pd.Series([0.1, 1.3, 2.7])
        >>> s.round()
        0    0.0
        1    1.0
        2    3.0
        dtype: float64
        """
        nv.validate_round(args, kwargs)
        result = com.values_from_object(self).round(decimals)
        result = self._constructor(result, index=self.index).__finalize__(self)

        return result

    def quantile(self, q=0.5, interpolation="linear"):
        """
        Return value at the given quantile.

        Parameters
        ----------
        q : float or array-like, default 0.5 (50% quantile)
            The quantile(s) to compute, which can lie in range: 0 <= q <= 1.
        interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
            This optional parameter specifies the interpolation method to use,
            when the desired quantile lies between two data points `i` and `j`:

                * linear: `i + (j - i) * fraction`, where `fraction` is the
                  fractional part of the index surrounded by `i` and `j`.
                * lower: `i`.
                * higher: `j`.
                * nearest: `i` or `j` whichever is nearest.
                * midpoint: (`i` + `j`) / 2.

        Returns
        -------
        float or Series
            If ``q`` is an array, a Series will be returned where the
            index is ``q`` and the values are the quantiles, otherwise
            a float will be returned.

        See Also
        --------
        core.window.Rolling.quantile
        numpy.percentile

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4])
        >>> s.quantile(.5)
        2.5
        >>> s.quantile([.25, .5, .75])
        0.25    1.75
        0.50    2.50
        0.75    3.25
        dtype: float64
        """

        validate_percentile(q)

        # We dispatch to DataFrame so that core.internals only has to worry
        #  about 2D cases.
        df = self.to_frame()

        result = df.quantile(q=q, interpolation=interpolation, numeric_only=False)
        if result.ndim == 2:
            result = result.iloc[:, 0]

        if is_list_like(q):
            result.name = self.name
            return self._constructor(result, index=Float64Index(q), name=self.name)
        else:
            # scalar
            return result.iloc[0]

    def corr(self, other, method="pearson", min_periods=None) -> float:
        """
        Compute correlation with `other` Series, excluding missing values.

        Parameters
        ----------
        other : Series
            Series with which to compute the correlation.
        method : {'pearson', 'kendall', 'spearman'} or callable
            Method used to compute correlation:

            - pearson : Standard correlation coefficient
            - kendall : Kendall Tau correlation coefficient
            - spearman : Spearman rank correlation
            - callable: Callable with input two 1d ndarrays and returning a float.

            .. versionadded:: 0.24.0
                Note that the returned matrix from corr will have 1 along the
                diagonals and will be symmetric regardless of the callable's
                behavior.
        min_periods : int, optional
            Minimum number of observations needed to have a valid result.

        Returns
        -------
        float
            Correlation with other.

        Examples
        --------
        >>> def histogram_intersection(a, b):
        ...     v = np.minimum(a, b).sum().round(decimals=1)
        ...     return v
        >>> s1 = pd.Series([.2, .0, .6, .2])
        >>> s2 = pd.Series([.3, .6, .0, .1])
        >>> s1.corr(s2, method=histogram_intersection)
        0.3
        """
        this, other = self.align(other, join="inner", copy=False)
        if len(this) == 0:
            return np.nan

        if method in ["pearson", "spearman", "kendall"] or callable(method):
            return nanops.nancorr(
                this.values, other.values, method=method, min_periods=min_periods
            )

        raise ValueError(
            "method must be either 'pearson', "
            "'spearman', 'kendall', or a callable, "
            f"'{method}' was supplied"
        )

    def cov(self, other, min_periods=None) -> float:
        """
        Compute covariance with Series, excluding missing values.

        Parameters
        ----------
        other : Series
            Series with which to compute the covariance.
        min_periods : int, optional
            Minimum number of observations needed to have a valid result.

        Returns
        -------
        float
            Covariance between Series and other normalized by N-1
            (unbiased estimator).

        Examples
        --------
        >>> s1 = pd.Series([0.90010907, 0.13484424, 0.62036035])
        >>> s2 = pd.Series([0.12528585, 0.26962463, 0.51111198])
        >>> s1.cov(s2)
        -0.01685762652715874
        """
        this, other = self.align(other, join="inner", copy=False)
        if len(this) == 0:
            return np.nan
        return nanops.nancov(this.values, other.values, min_periods=min_periods)

    def diff(self, periods=1) -> "Series":
        """
        First discrete difference of element.

        Calculates the difference of a Series element compared with another
        element in the Series (default is element in previous row).

        Parameters
        ----------
        periods : int, default 1
            Periods to shift for calculating difference, accepts negative
            values.

        Returns
        -------
        Series
            First differences of the Series.

        See Also
        --------
        Series.pct_change: Percent change over given number of periods.
        Series.shift: Shift index by desired number of periods with an
            optional time freq.
        DataFrame.diff: First discrete difference of object.

        Notes
        -----
        For boolean dtypes, this uses :meth:`operator.xor` rather than
        :meth:`operator.sub`.

        Examples
        --------
        Difference with previous row

        >>> s = pd.Series([1, 1, 2, 3, 5, 8])
        >>> s.diff()
        0    NaN
        1    0.0
        2    1.0
        3    1.0
        4    2.0
        5    3.0
        dtype: float64

        Difference with 3rd previous row

        >>> s.diff(periods=3)
        0    NaN
        1    NaN
        2    NaN
        3    2.0
        4    4.0
        5    6.0
        dtype: float64

        Difference with following row

        >>> s.diff(periods=-1)
        0    0.0
        1   -1.0
        2   -1.0
        3   -2.0
        4   -3.0
        5    NaN
        dtype: float64
        """
        result = algorithms.diff(self.array, periods)
        return self._constructor(result, index=self.index).__finalize__(self)

    def autocorr(self, lag=1) -> float:
        """
        Compute the lag-N autocorrelation.

        This method computes the Pearson correlation between
        the Series and its shifted self.

        Parameters
        ----------
        lag : int, default 1
            Number of lags to apply before performing autocorrelation.

        Returns
        -------
        float
            The Pearson correlation between self and self.shift(lag).

        See Also
        --------
        Series.corr : Compute the correlation between two Series.
        Series.shift : Shift index by desired number of periods.
        DataFrame.corr : Compute pairwise correlation of columns.
        DataFrame.corrwith : Compute pairwise correlation between rows or
            columns of two DataFrame objects.

        Notes
        -----
        If the Pearson correlation is not well defined return 'NaN'.

        Examples
        --------
        >>> s = pd.Series([0.25, 0.5, 0.2, -0.05])
        >>> s.autocorr()  # doctest: +ELLIPSIS
        0.10355...
        >>> s.autocorr(lag=2)  # doctest: +ELLIPSIS
        -0.99999...

        If the Pearson correlation is not well defined, then 'NaN' is returned.

        >>> s = pd.Series([1, 0, 0, 0])
        >>> s.autocorr()
        nan
        """
        return self.corr(self.shift(lag))

    def dot(self, other):
        """
        Compute the dot product between the Series and the columns of other.

        This method computes the dot product between the Series and another
        one, or the Series and each columns of a DataFrame, or the Series and
        each columns of an array.

        It can also be called using `self @ other` in Python >= 3.5.

        Parameters
        ----------
        other : Series, DataFrame or array-like
            The other object to compute the dot product with its columns.

        Returns
        -------
        scalar, Series or numpy.ndarray
            Return the dot product of the Series and other if other is a
            Series, the Series of the dot product of Series and each rows of
            other if other is a DataFrame or a numpy.ndarray between the Series
            and each columns of the numpy array.

        See Also
        --------
        DataFrame.dot: Compute the matrix product with the DataFrame.
        Series.mul: Multiplication of series and other, element-wise.

        Notes
        -----
        The Series and other has to share the same index if other is a Series
        or a DataFrame.

        Examples
        --------
        >>> s = pd.Series([0, 1, 2, 3])
        >>> other = pd.Series([-1, 2, -3, 4])
        >>> s.dot(other)
        8
        >>> s @ other
        8
        >>> df = pd.DataFrame([[0, 1], [-2, 3], [4, -5], [6, 7]])
        >>> s.dot(df)
        0    24
        1    14
        dtype: int64
        >>> arr = np.array([[0, 1], [-2, 3], [4, -5], [6, 7]])
        >>> s.dot(arr)
        array([24, 14])
        """
        if isinstance(other, (Series, ABCDataFrame)):
            common = self.index.union(other.index)
            if len(common) > len(self.index) or len(common) > len(other.index):
                raise ValueError("matrices are not aligned")

            left = self.reindex(index=common, copy=False)
            right = other.reindex(index=common, copy=False)
            lvals = left.values
            rvals = right.values
        else:
            lvals = self.values
            rvals = np.asarray(other)
            if lvals.shape[0] != rvals.shape[0]:
                raise Exception(
                    f"Dot product shape mismatch, {lvals.shape} vs {rvals.shape}"
                )

        if isinstance(other, ABCDataFrame):
            return self._constructor(
                np.dot(lvals, rvals), index=other.columns
            ).__finalize__(self)
        elif isinstance(other, Series):
            return np.dot(lvals, rvals)
        elif isinstance(rvals, np.ndarray):
            return np.dot(lvals, rvals)
        else:  # pragma: no cover
            raise TypeError(f"unsupported type: {type(other)}")

    def __matmul__(self, other):
        """
        Matrix multiplication using binary `@` operator in Python>=3.5.
        """
        return self.dot(other)

    def __rmatmul__(self, other):
        """
        Matrix multiplication using binary `@` operator in Python>=3.5.
        """
        return self.dot(np.transpose(other))

    @Substitution(klass="Series")
    @Appender(base._shared_docs["searchsorted"])
    def searchsorted(self, value, side="left", sorter=None):
        return algorithms.searchsorted(self._values, value, side=side, sorter=sorter)

    # -------------------------------------------------------------------
    # Combination

    def append(self, to_append, ignore_index=False, verify_integrity=False):
        """
        Concatenate two or more Series.

        Parameters
        ----------
        to_append : Series or list/tuple of Series
            Series to append with self.
        ignore_index : bool, default False
            If True, do not use the index labels.
        verify_integrity : bool, default False
            If True, raise Exception on creating index with duplicates.

        Returns
        -------
        Series
            Concatenated Series.

        See Also
        --------
        concat : General function to concatenate DataFrame or Series objects.

        Notes
        -----
        Iteratively appending to a Series can be more computationally intensive
        than a single concatenate. A better solution is to append values to a
        list and then concatenate the list with the original Series all at
        once.

        Examples
        --------
        >>> s1 = pd.Series([1, 2, 3])
        >>> s2 = pd.Series([4, 5, 6])
        >>> s3 = pd.Series([4, 5, 6], index=[3, 4, 5])
        >>> s1.append(s2)
        0    1
        1    2
        2    3
        0    4
        1    5
        2    6
        dtype: int64

        >>> s1.append(s3)
        0    1
        1    2
        2    3
        3    4
        4    5
        5    6
        dtype: int64

        With `ignore_index` set to True:

        >>> s1.append(s2, ignore_index=True)
        0    1
        1    2
        2    3
        3    4
        4    5
        5    6
        dtype: int64

        With `verify_integrity` set to True:

        >>> s1.append(s2, verify_integrity=True)
        Traceback (most recent call last):
        ...
        ValueError: Indexes have overlapping values: [0, 1, 2]
        """
        from pandas.core.reshape.concat import concat

        if isinstance(to_append, (list, tuple)):
            to_concat = [self]
            to_concat.extend(to_append)
        else:
            to_concat = [self, to_append]
        return concat(
            to_concat, ignore_index=ignore_index, verify_integrity=verify_integrity
        )

    def _binop(self, other, func, level=None, fill_value=None):
        """
        Perform generic binary operation with optional fill value.

        Parameters
        ----------
        other : Series
        func : binary operator
        fill_value : float or object
            Value to substitute for NA/null values. If both Series are NA in a
            location, the result will be NA regardless of the passed fill value.
        level : int or level name, default None
            Broadcast across a level, matching Index values on the
            passed MultiIndex level.

        Returns
        -------
        Series
        """

        if not isinstance(other, Series):
            raise AssertionError("Other operand must be Series")

        new_index = self.index
        this = self

        if not self.index.equals(other.index):
            this, other = self.align(other, level=level, join="outer", copy=False)
            new_index = this.index

        this_vals, other_vals = ops.fill_binop(this.values, other.values, fill_value)

        with np.errstate(all="ignore"):
            result = func(this_vals, other_vals)

        name = ops.get_op_result_name(self, other)
        ret = ops._construct_result(self, result, new_index, name)
        return ret

    def combine(self, other, func, fill_value=None) -> "Series":
        """
        Combine the Series with a Series or scalar according to `func`.

        Combine the Series and `other` using `func` to perform elementwise
        selection for combined Series.
        `fill_value` is assumed when value is missing at some index
        from one of the two objects being combined.

        Parameters
        ----------
        other : Series or scalar
            The value(s) to be combined with the `Series`.
        func : function
            Function that takes two scalars as inputs and returns an element.
        fill_value : scalar, optional
            The value to assume when an index is missing from
            one Series or the other. The default specifies to use the
            appropriate NaN value for the underlying dtype of the Series.

        Returns
        -------
        Series
            The result of combining the Series with the other object.

        See Also
        --------
        Series.combine_first : Combine Series values, choosing the calling
            Series' values first.

        Examples
        --------
        Consider 2 Datasets ``s1`` and ``s2`` containing
        highest clocked speeds of different birds.

        >>> s1 = pd.Series({'falcon': 330.0, 'eagle': 160.0})
        >>> s1
        falcon    330.0
        eagle     160.0
        dtype: float64
        >>> s2 = pd.Series({'falcon': 345.0, 'eagle': 200.0, 'duck': 30.0})
        >>> s2
        falcon    345.0
        eagle     200.0
        duck       30.0
        dtype: float64

        Now, to combine the two datasets and view the highest speeds
        of the birds across the two datasets

        >>> s1.combine(s2, max)
        duck        NaN
        eagle     200.0
        falcon    345.0
        dtype: float64

        In the previous example, the resulting value for duck is missing,
        because the maximum of a NaN and a float is a NaN.
        So, in the example, we set ``fill_value=0``,
        so the maximum value returned will be the value from some dataset.

        >>> s1.combine(s2, max, fill_value=0)
        duck       30.0
        eagle     200.0
        falcon    345.0
        dtype: float64
        """
        if fill_value is None:
            fill_value = na_value_for_dtype(self.dtype, compat=False)

        if isinstance(other, Series):
            # If other is a Series, result is based on union of Series,
            # so do this element by element
            new_index = self.index.union(other.index)
            new_name = ops.get_op_result_name(self, other)
            new_values = []
            for idx in new_index:
                lv = self.get(idx, fill_value)
                rv = other.get(idx, fill_value)
                with np.errstate(all="ignore"):
                    new_values.append(func(lv, rv))
        else:
            # Assume that other is a scalar, so apply the function for
            # each element in the Series
            new_index = self.index
            with np.errstate(all="ignore"):
                new_values = [func(lv, other) for lv in self._values]
            new_name = self.name

        if is_categorical_dtype(self.values):
            pass
        elif is_extension_array_dtype(self.values):
            # The function can return something of any type, so check
            # if the type is compatible with the calling EA.
            new_values = try_cast_to_ea(self._values, new_values)
        return self._constructor(new_values, index=new_index, name=new_name)

    def combine_first(self, other) -> "Series":
        """
        Combine Series values, choosing the calling Series's values first.

        Parameters
        ----------
        other : Series
            The value(s) to be combined with the `Series`.

        Returns
        -------
        Series
            The result of combining the Series with the other object.

        See Also
        --------
        Series.combine : Perform elementwise operation on two Series
            using a given function.

        Notes
        -----
        Result index will be the union of the two indexes.

        Examples
        --------
        >>> s1 = pd.Series([1, np.nan])
        >>> s2 = pd.Series([3, 4])
        >>> s1.combine_first(s2)
        0    1.0
        1    4.0
        dtype: float64
        """
        new_index = self.index.union(other.index)
        this = self.reindex(new_index, copy=False)
        other = other.reindex(new_index, copy=False)
        if this.dtype.kind == "M" and other.dtype.kind != "M":
            other = to_datetime(other)

        return this.where(notna(this), other)

    def update(self, other) -> None:
        """
        Modify Series in place using non-NA values from passed
        Series. Aligns on index.

        Parameters
        ----------
        other : Series

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s.update(pd.Series([4, 5, 6]))
        >>> s
        0    4
        1    5
        2    6
        dtype: int64

        >>> s = pd.Series(['a', 'b', 'c'])
        >>> s.update(pd.Series(['d', 'e'], index=[0, 2]))
        >>> s
        0    d
        1    b
        2    e
        dtype: object

        >>> s = pd.Series([1, 2, 3])
        >>> s.update(pd.Series([4, 5, 6, 7, 8]))
        >>> s
        0    4
        1    5
        2    6
        dtype: int64

        If ``other`` contains NaNs the corresponding values are not updated
        in the original Series.

        >>> s = pd.Series([1, 2, 3])
        >>> s.update(pd.Series([4, np.nan, 6]))
        >>> s
        0    4
        1    2
        2    6
        dtype: int64
        """
        other = other.reindex_like(self)
        mask = notna(other)

        self._data = self._data.putmask(mask=mask, new=other, inplace=True)
        self._maybe_update_cacher()

    # ----------------------------------------------------------------------
    # Reindexing, sorting

    def sort_values(
        self,
        axis=0,
        ascending=True,
        inplace: bool = False,
        kind: str = "quicksort",
        na_position: str = "last",
        ignore_index: bool = False,
    ):
        """
        Sort by the values.

        Sort a Series in ascending or descending order by some
        criterion.

        Parameters
        ----------
        axis : {0 or 'index'}, default 0
            Axis to direct sorting. The value 'index' is accepted for
            compatibility with DataFrame.sort_values.
        ascending : bool, default True
            If True, sort values in ascending order, otherwise descending.
        inplace : bool, default False
            If True, perform operation in-place.
        kind : {'quicksort', 'mergesort' or 'heapsort'}, default 'quicksort'
            Choice of sorting algorithm. See also :func:`numpy.sort` for more
            information. 'mergesort' is the only stable  algorithm.
        na_position : {'first' or 'last'}, default 'last'
            Argument 'first' puts NaNs at the beginning, 'last' puts NaNs at
            the end.
        ignore_index : bool, default False
             If True, the resulting axis will be labeled 0, 1, …, n - 1.

             .. versionadded:: 1.0.0

        Returns
        -------
        Series
            Series ordered by values.

        See Also
        --------
        Series.sort_index : Sort by the Series indices.
        DataFrame.sort_values : Sort DataFrame by the values along either axis.
        DataFrame.sort_index : Sort DataFrame by indices.

        Examples
        --------
        >>> s = pd.Series([np.nan, 1, 3, 10, 5])
        >>> s
        0     NaN
        1     1.0
        2     3.0
        3     10.0
        4     5.0
        dtype: float64

        Sort values ascending order (default behaviour)

        >>> s.sort_values(ascending=True)
        1     1.0
        2     3.0
        4     5.0
        3    10.0
        0     NaN
        dtype: float64

        Sort values descending order

        >>> s.sort_values(ascending=False)
        3    10.0
        4     5.0
        2     3.0
        1     1.0
        0     NaN
        dtype: float64

        Sort values inplace

        >>> s.sort_values(ascending=False, inplace=True)
        >>> s
        3    10.0
        4     5.0
        2     3.0
        1     1.0
        0     NaN
        dtype: float64

        Sort values putting NAs first

        >>> s.sort_values(na_position='first')
        0     NaN
        1     1.0
        2     3.0
        4     5.0
        3    10.0
        dtype: float64

        Sort a series of strings

        >>> s = pd.Series(['z', 'b', 'd', 'a', 'c'])
        >>> s
        0    z
        1    b
        2    d
        3    a
        4    c
        dtype: object

        >>> s.sort_values()
        3    a
        1    b
        4    c
        2    d
        0    z
        dtype: object
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        # Validate the axis parameter
        self._get_axis_number(axis)

        # GH 5856/5853
        if inplace and self._is_cached:
            raise ValueError(
                "This Series is a view of some other array, to "
                "sort in-place you must create a copy"
            )

        def _try_kind_sort(arr):
            # easier to ask forgiveness than permission
            try:
                # if kind==mergesort, it can fail for object dtype
                return arr.argsort(kind=kind)
            except TypeError:
                # stable sort not available for object dtype
                # uses the argsort default quicksort
                return arr.argsort(kind="quicksort")

        arr = self._values
        sorted_index = np.empty(len(self), dtype=np.int32)

        bad = isna(arr)

        good = ~bad
        idx = ibase.default_index(len(self))

        argsorted = _try_kind_sort(arr[good])

        if is_list_like(ascending):
            if len(ascending) != 1:
                raise ValueError(
                    f"Length of ascending ({len(ascending)}) must be 1 for Series"
                )
            ascending = ascending[0]

        if not is_bool(ascending):
            raise ValueError("ascending must be boolean")

        if not ascending:
            argsorted = argsorted[::-1]

        if na_position == "last":
            n = good.sum()
            sorted_index[:n] = idx[good][argsorted]
            sorted_index[n:] = idx[bad]
        elif na_position == "first":
            n = bad.sum()
            sorted_index[n:] = idx[good][argsorted]
            sorted_index[:n] = idx[bad]
        else:
            raise ValueError(f"invalid na_position: {na_position}")

        result = self._constructor(arr[sorted_index], index=self.index[sorted_index])

        if ignore_index:
            result.index = ibase.default_index(len(sorted_index))

        if inplace:
            self._update_inplace(result)
        else:
            return result.__finalize__(self)

    def sort_index(
        self,
        axis=0,
        level=None,
        ascending=True,
        inplace=False,
        kind="quicksort",
        na_position="last",
        sort_remaining=True,
        ignore_index: bool = False,
    ):
        """
        Sort Series by index labels.

        Returns a new Series sorted by label if `inplace` argument is
        ``False``, otherwise updates the original series and returns None.

        Parameters
        ----------
        axis : int, default 0
            Axis to direct sorting. This can only be 0 for Series.
        level : int, optional
            If not None, sort on values in specified index level(s).
        ascending : bool, default true
            Sort ascending vs. descending.
        inplace : bool, default False
            If True, perform operation in-place.
        kind : {'quicksort', 'mergesort', 'heapsort'}, default 'quicksort'
            Choice of sorting algorithm. See also :func:`numpy.sort` for more
            information.  'mergesort' is the only stable algorithm. For
            DataFrames, this option is only applied when sorting on a single
            column or label.
        na_position : {'first', 'last'}, default 'last'
            If 'first' puts NaNs at the beginning, 'last' puts NaNs at the end.
            Not implemented for MultiIndex.
        sort_remaining : bool, default True
            If True and sorting by level and index is multilevel, sort by other
            levels too (in order) after sorting by specified level.
        ignore_index : bool, default False
            If True, the resulting axis will be labeled 0, 1, …, n - 1.

            .. versionadded:: 1.0.0

        Returns
        -------
        Series
            The original Series sorted by the labels.

        See Also
        --------
        DataFrame.sort_index: Sort DataFrame by the index.
        DataFrame.sort_values: Sort DataFrame by the value.
        Series.sort_values : Sort Series by the value.

        Examples
        --------
        >>> s = pd.Series(['a', 'b', 'c', 'd'], index=[3, 2, 1, 4])
        >>> s.sort_index()
        1    c
        2    b
        3    a
        4    d
        dtype: object

        Sort Descending

        >>> s.sort_index(ascending=False)
        4    d
        3    a
        2    b
        1    c
        dtype: object

        Sort Inplace

        >>> s.sort_index(inplace=True)
        >>> s
        1    c
        2    b
        3    a
        4    d
        dtype: object

        By default NaNs are put at the end, but use `na_position` to place
        them at the beginning

        >>> s = pd.Series(['a', 'b', 'c', 'd'], index=[3, 2, 1, np.nan])
        >>> s.sort_index(na_position='first')
        NaN     d
         1.0    c
         2.0    b
         3.0    a
        dtype: object

        Specify index level to sort

        >>> arrays = [np.array(['qux', 'qux', 'foo', 'foo',
        ...                     'baz', 'baz', 'bar', 'bar']),
        ...           np.array(['two', 'one', 'two', 'one',
        ...                     'two', 'one', 'two', 'one'])]
        >>> s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8], index=arrays)
        >>> s.sort_index(level=1)
        bar  one    8
        baz  one    6
        foo  one    4
        qux  one    2
        bar  two    7
        baz  two    5
        foo  two    3
        qux  two    1
        dtype: int64

        Does not sort by remaining levels when sorting by levels

        >>> s.sort_index(level=1, sort_remaining=False)
        qux  one    2
        foo  one    4
        baz  one    6
        bar  one    8
        qux  two    1
        foo  two    3
        baz  two    5
        bar  two    7
        dtype: int64
        """
        # TODO: this can be combined with DataFrame.sort_index impl as
        # almost identical
        inplace = validate_bool_kwarg(inplace, "inplace")
        # Validate the axis parameter
        self._get_axis_number(axis)
        index = self.index

        if level is not None:
            new_index, indexer = index.sortlevel(
                level, ascending=ascending, sort_remaining=sort_remaining
            )
        elif isinstance(index, MultiIndex):
            from pandas.core.sorting import lexsort_indexer

            labels = index._sort_levels_monotonic()
            indexer = lexsort_indexer(
                labels._get_codes_for_sorting(),
                orders=ascending,
                na_position=na_position,
            )
        else:
            from pandas.core.sorting import nargsort

            # Check monotonic-ness before sort an index
            # GH11080
            if (ascending and index.is_monotonic_increasing) or (
                not ascending and index.is_monotonic_decreasing
            ):
                if inplace:
                    return
                else:
                    return self.copy()

            indexer = nargsort(
                index, kind=kind, ascending=ascending, na_position=na_position
            )

        indexer = ensure_platform_int(indexer)
        new_index = index.take(indexer)
        new_index = new_index._sort_levels_monotonic()

        new_values = self._values.take(indexer)
        result = self._constructor(new_values, index=new_index)

        if ignore_index:
            result.index = ibase.default_index(len(result))

        if inplace:
            self._update_inplace(result)
        else:
            return result.__finalize__(self)

    def argsort(self, axis=0, kind="quicksort", order=None) -> "Series":
        """
        Override ndarray.argsort. Argsorts the value, omitting NA/null values,
        and places the result in the same locations as the non-NA values.

        Parameters
        ----------
        axis : {0 or "index"}
            Has no effect but is accepted for compatibility with numpy.
        kind : {'mergesort', 'quicksort', 'heapsort'}, default 'quicksort'
            Choice of sorting algorithm. See np.sort for more
            information. 'mergesort' is the only stable algorithm.
        order : None
            Has no effect but is accepted for compatibility with numpy.

        Returns
        -------
        Series
            Positions of values within the sort order with -1 indicating
            nan values.

        See Also
        --------
        numpy.ndarray.argsort
        """
        values = self._values
        mask = isna(values)

        if mask.any():
            result = Series(-1, index=self.index, name=self.name, dtype="int64")
            notmask = ~mask
            result[notmask] = np.argsort(values[notmask], kind=kind)
            return self._constructor(result, index=self.index).__finalize__(self)
        else:
            return self._constructor(
                np.argsort(values, kind=kind), index=self.index, dtype="int64"
            ).__finalize__(self)

    def nlargest(self, n=5, keep="first") -> "Series":
        """
        Return the largest `n` elements.

        Parameters
        ----------
        n : int, default 5
            Return this many descending sorted values.
        keep : {'first', 'last', 'all'}, default 'first'
            When there are duplicate values that cannot all fit in a
            Series of `n` elements:

            - ``first`` : return the first `n` occurrences in order
                of appearance.
            - ``last`` : return the last `n` occurrences in reverse
                order of appearance.
            - ``all`` : keep all occurrences. This can result in a Series of
                size larger than `n`.

        Returns
        -------
        Series
            The `n` largest values in the Series, sorted in decreasing order.

        See Also
        --------
        Series.nsmallest: Get the `n` smallest elements.
        Series.sort_values: Sort Series by values.
        Series.head: Return the first `n` rows.

        Notes
        -----
        Faster than ``.sort_values(ascending=False).head(n)`` for small `n`
        relative to the size of the ``Series`` object.

        Examples
        --------
        >>> countries_population = {"Italy": 59000000, "France": 65000000,
        ...                         "Malta": 434000, "Maldives": 434000,
        ...                         "Brunei": 434000, "Iceland": 337000,
        ...                         "Nauru": 11300, "Tuvalu": 11300,
        ...                         "Anguilla": 11300, "Monserat": 5200}
        >>> s = pd.Series(countries_population)
        >>> s
        Italy       59000000
        France      65000000
        Malta         434000
        Maldives      434000
        Brunei        434000
        Iceland       337000
        Nauru          11300
        Tuvalu         11300
        Anguilla       11300
        Monserat        5200
        dtype: int64

        The `n` largest elements where ``n=5`` by default.

        >>> s.nlargest()
        France      65000000
        Italy       59000000
        Malta         434000
        Maldives      434000
        Brunei        434000
        dtype: int64

        The `n` largest elements where ``n=3``. Default `keep` value is 'first'
        so Malta will be kept.

        >>> s.nlargest(3)
        France    65000000
        Italy     59000000
        Malta       434000
        dtype: int64

        The `n` largest elements where ``n=3`` and keeping the last duplicates.
        Brunei will be kept since it is the last with value 434000 based on
        the index order.

        >>> s.nlargest(3, keep='last')
        France      65000000
        Italy       59000000
        Brunei        434000
        dtype: int64

        The `n` largest elements where ``n=3`` with all duplicates kept. Note
        that the returned Series has five elements due to the three duplicates.

        >>> s.nlargest(3, keep='all')
        France      65000000
        Italy       59000000
        Malta         434000
        Maldives      434000
        Brunei        434000
        dtype: int64
        """
        return algorithms.SelectNSeries(self, n=n, keep=keep).nlargest()

    def nsmallest(self, n=5, keep="first") -> "Series":
        """
        Return the smallest `n` elements.

        Parameters
        ----------
        n : int, default 5
            Return this many ascending sorted values.
        keep : {'first', 'last', 'all'}, default 'first'
            When there are duplicate values that cannot all fit in a
            Series of `n` elements:

            - ``first`` : return the first `n` occurrences in order
                of appearance.
            - ``last`` : return the last `n` occurrences in reverse
                order of appearance.
            - ``all`` : keep all occurrences. This can result in a Series of
                size larger than `n`.

        Returns
        -------
        Series
            The `n` smallest values in the Series, sorted in increasing order.

        See Also
        --------
        Series.nlargest: Get the `n` largest elements.
        Series.sort_values: Sort Series by values.
        Series.head: Return the first `n` rows.

        Notes
        -----
        Faster than ``.sort_values().head(n)`` for small `n` relative to
        the size of the ``Series`` object.

        Examples
        --------
        >>> countries_population = {"Italy": 59000000, "France": 65000000,
        ...                         "Brunei": 434000, "Malta": 434000,
        ...                         "Maldives": 434000, "Iceland": 337000,
        ...                         "Nauru": 11300, "Tuvalu": 11300,
        ...                         "Anguilla": 11300, "Monserat": 5200}
        >>> s = pd.Series(countries_population)
        >>> s
        Italy       59000000
        France      65000000
        Brunei        434000
        Malta         434000
        Maldives      434000
        Iceland       337000
        Nauru          11300
        Tuvalu         11300
        Anguilla       11300
        Monserat        5200
        dtype: int64

        The `n` smallest elements where ``n=5`` by default.

        >>> s.nsmallest()
        Monserat      5200
        Nauru        11300
        Tuvalu       11300
        Anguilla     11300
        Iceland     337000
        dtype: int64

        The `n` smallest elements where ``n=3``. Default `keep` value is
        'first' so Nauru and Tuvalu will be kept.

        >>> s.nsmallest(3)
        Monserat     5200
        Nauru       11300
        Tuvalu      11300
        dtype: int64

        The `n` smallest elements where ``n=3`` and keeping the last
        duplicates. Anguilla and Tuvalu will be kept since they are the last
        with value 11300 based on the index order.

        >>> s.nsmallest(3, keep='last')
        Monserat     5200
        Anguilla    11300
        Tuvalu      11300
        dtype: int64

        The `n` smallest elements where ``n=3`` with all duplicates kept. Note
        that the returned Series has four elements due to the three duplicates.

        >>> s.nsmallest(3, keep='all')
        Monserat     5200
        Nauru       11300
        Tuvalu      11300
        Anguilla    11300
        dtype: int64
        """
        return algorithms.SelectNSeries(self, n=n, keep=keep).nsmallest()

    def swaplevel(self, i=-2, j=-1, copy=True) -> "Series":
        """
        Swap levels i and j in a :class:`MultiIndex`.

        Default is to swap the two innermost levels of the index.

        Parameters
        ----------
        i, j : int, str
            Level of the indices to be swapped. Can pass level name as string.
        copy : bool, default True
            Whether to copy underlying data.

        Returns
        -------
        Series
            Series with levels swapped in MultiIndex.
        """
        assert isinstance(self.index, ABCMultiIndex)
        new_index = self.index.swaplevel(i, j)
        return self._constructor(self._values, index=new_index, copy=copy).__finalize__(
            self
        )

    def reorder_levels(self, order) -> "Series":
        """
        Rearrange index levels using input order.

        May not drop or duplicate levels.

        Parameters
        ----------
        order : list of int representing new level order
            Reference level by number or key.

        Returns
        -------
        type of caller (new object)
        """
        if not isinstance(self.index, MultiIndex):  # pragma: no cover
            raise Exception("Can only reorder levels on a hierarchical axis.")

        result = self.copy()
        assert isinstance(result.index, ABCMultiIndex)
        result.index = result.index.reorder_levels(order)
        return result

    def explode(self) -> "Series":
        """
        Transform each element of a list-like to a row, replicating the
        index values.

        .. versionadded:: 0.25.0

        Returns
        -------
        Series
            Exploded lists to rows; index will be duplicated for these rows.

        See Also
        --------
        Series.str.split : Split string values on specified separator.
        Series.unstack : Unstack, a.k.a. pivot, Series with MultiIndex
            to produce DataFrame.
        DataFrame.melt : Unpivot a DataFrame from wide format to long format.
        DataFrame.explode : Explode a DataFrame from list-like
            columns to long format.

        Notes
        -----
        This routine will explode list-likes including lists, tuples,
        Series, and np.ndarray. The result dtype of the subset rows will
        be object. Scalars will be returned unchanged. Empty list-likes will
        result in a np.nan for that row.

        Examples
        --------
        >>> s = pd.Series([[1, 2, 3], 'foo', [], [3, 4]])
        >>> s
        0    [1, 2, 3]
        1          foo
        2           []
        3       [3, 4]
        dtype: object

        >>> s.explode()
        0      1
        0      2
        0      3
        1    foo
        2    NaN
        3      3
        3      4
        dtype: object
        """
        if not len(self) or not is_object_dtype(self):
            return self.copy()

        values, counts = reshape.explode(np.asarray(self.array))

        result = Series(values, index=self.index.repeat(counts), name=self.name)
        return result

    def unstack(self, level=-1, fill_value=None):
        """
        Unstack, also known as pivot, Series with MultiIndex to produce DataFrame.
        The level involved will automatically get sorted.

        Parameters
        ----------
        level : int, str, or list of these, default last level
            Level(s) to unstack, can pass level name.
        fill_value : scalar value, default None
            Value to use when replacing NaN values.

        Returns
        -------
        DataFrame
            Unstacked Series.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4],
        ...               index=pd.MultiIndex.from_product([['one', 'two'],
        ...                                                 ['a', 'b']]))
        >>> s
        one  a    1
             b    2
        two  a    3
             b    4
        dtype: int64

        >>> s.unstack(level=-1)
             a  b
        one  1  2
        two  3  4

        >>> s.unstack(level=0)
           one  two
        a    1    3
        b    2    4
        """
        from pandas.core.reshape.reshape import unstack

        return unstack(self, level, fill_value)

    # ----------------------------------------------------------------------
    # function application

    def map(self, arg, na_action=None) -> "Series":
        """
        Map values of Series according to input correspondence.

        Used for substituting each value in a Series with another value,
        that may be derived from a function, a ``dict`` or
        a :class:`Series`.

        Parameters
        ----------
        arg : function, collections.abc.Mapping subclass or Series
            Mapping correspondence.
        na_action : {None, 'ignore'}, default None
            If 'ignore', propagate NaN values, without passing them to the
            mapping correspondence.

        Returns
        -------
        Series
            Same index as caller.

        See Also
        --------
        Series.apply : For applying more complex functions on a Series.
        DataFrame.apply : Apply a function row-/column-wise.
        DataFrame.applymap : Apply a function elementwise on a whole DataFrame.

        Notes
        -----
        When ``arg`` is a dictionary, values in Series that are not in the
        dictionary (as keys) are converted to ``NaN``. However, if the
        dictionary is a ``dict`` subclass that defines ``__missing__`` (i.e.
        provides a method for default values), then this default is used
        rather than ``NaN``.

        Examples
        --------
        >>> s = pd.Series(['cat', 'dog', np.nan, 'rabbit'])
        >>> s
        0      cat
        1      dog
        2      NaN
        3   rabbit
        dtype: object

        ``map`` accepts a ``dict`` or a ``Series``. Values that are not found
        in the ``dict`` are converted to ``NaN``, unless the dict has a default
        value (e.g. ``defaultdict``):

        >>> s.map({'cat': 'kitten', 'dog': 'puppy'})
        0   kitten
        1    puppy
        2      NaN
        3      NaN
        dtype: object

        It also accepts a function:

        >>> s.map('I am a {}'.format)
        0       I am a cat
        1       I am a dog
        2       I am a nan
        3    I am a rabbit
        dtype: object

        To avoid applying the function to missing values (and keep them as
        ``NaN``) ``na_action='ignore'`` can be used:

        >>> s.map('I am a {}'.format, na_action='ignore')
        0     I am a cat
        1     I am a dog
        2            NaN
        3  I am a rabbit
        dtype: object
        """
        new_values = super()._map_values(arg, na_action=na_action)
        return self._constructor(new_values, index=self.index).__finalize__(self)

    def _gotitem(self, key, ndim, subset=None) -> "Series":
        """
        Sub-classes to define. Return a sliced object.

        Parameters
        ----------
        key : string / list of selections
        ndim : 1,2
            Requested ndim of result.
        subset : object, default None
            Subset to act on.
        """
        return self

    _agg_see_also_doc = dedent(
        """
    See Also
    --------
    Series.apply : Invoke function on a Series.
    Series.transform : Transform function producing a Series with like indexes.
    """
    )

    _agg_examples_doc = dedent(
        """
    Examples
    --------
    >>> s = pd.Series([1, 2, 3, 4])
    >>> s
    0    1
    1    2
    2    3
    3    4
    dtype: int64

    >>> s.agg('min')
    1

    >>> s.agg(['min', 'max'])
    min   1
    max   4
    dtype: int64
    """
    )

    @Substitution(
        see_also=_agg_see_also_doc,
        examples=_agg_examples_doc,
        versionadded="\n.. versionadded:: 0.20.0\n",
        **_shared_doc_kwargs,
    )
    @Appender(generic._shared_docs["aggregate"])
    def aggregate(self, func, axis=0, *args, **kwargs):
        # Validate the axis parameter
        self._get_axis_number(axis)
        result, how = self._aggregate(func, *args, **kwargs)
        if result is None:

            # we can be called from an inner function which
            # passes this meta-data
            kwargs.pop("_axis", None)
            kwargs.pop("_level", None)

            # try a regular apply, this evaluates lambdas
            # row-by-row; however if the lambda is expected a Series
            # expression, e.g.: lambda x: x-x.quantile(0.25)
            # this will fail, so we can try a vectorized evaluation

            # we cannot FIRST try the vectorized evaluation, because
            # then .agg and .apply would have different semantics if the
            # operation is actually defined on the Series, e.g. str
            try:
                result = self.apply(func, *args, **kwargs)
            except (ValueError, AttributeError, TypeError):
                result = func(self, *args, **kwargs)

        return result

    agg = aggregate

    @Appender(generic._shared_docs["transform"] % _shared_doc_kwargs)
    def transform(self, func, axis=0, *args, **kwargs):
        # Validate the axis parameter
        self._get_axis_number(axis)
        return super().transform(func, *args, **kwargs)

    def apply(self, func, convert_dtype=True, args=(), **kwds):
        """
        Invoke function on values of Series.

        Can be ufunc (a NumPy function that applies to the entire Series)
        or a Python function that only works on single values.

        Parameters
        ----------
        func : function
            Python function or NumPy ufunc to apply.
        convert_dtype : bool, default True
            Try to find better dtype for elementwise function results. If
            False, leave as dtype=object.
        args : tuple
            Positional arguments passed to func after the series value.
        **kwds
            Additional keyword arguments passed to func.

        Returns
        -------
        Series or DataFrame
            If func returns a Series object the result will be a DataFrame.

        See Also
        --------
        Series.map: For element-wise operations.
        Series.agg: Only perform aggregating type operations.
        Series.transform: Only perform transforming type operations.

        Examples
        --------
        Create a series with typical summer temperatures for each city.

        >>> s = pd.Series([20, 21, 12],
        ...               index=['London', 'New York', 'Helsinki'])
        >>> s
        London      20
        New York    21
        Helsinki    12
        dtype: int64

        Square the values by defining a function and passing it as an
        argument to ``apply()``.

        >>> def square(x):
        ...     return x ** 2
        >>> s.apply(square)
        London      400
        New York    441
        Helsinki    144
        dtype: int64

        Square the values by passing an anonymous function as an
        argument to ``apply()``.

        >>> s.apply(lambda x: x ** 2)
        London      400
        New York    441
        Helsinki    144
        dtype: int64

        Define a custom function that needs additional positional
        arguments and pass these additional arguments using the
        ``args`` keyword.

        >>> def subtract_custom_value(x, custom_value):
        ...     return x - custom_value

        >>> s.apply(subtract_custom_value, args=(5,))
        London      15
        New York    16
        Helsinki     7
        dtype: int64

        Define a custom function that takes keyword arguments
        and pass these arguments to ``apply``.

        >>> def add_custom_values(x, **kwargs):
        ...     for month in kwargs:
        ...         x += kwargs[month]
        ...     return x

        >>> s.apply(add_custom_values, june=30, july=20, august=25)
        London      95
        New York    96
        Helsinki    87
        dtype: int64

        Use a function from the Numpy library.

        >>> s.apply(np.log)
        London      2.995732
        New York    3.044522
        Helsinki    2.484907
        dtype: float64
        """
        if len(self) == 0:
            return self._constructor(dtype=self.dtype, index=self.index).__finalize__(
                self
            )

        # dispatch to agg
        if isinstance(func, (list, dict)):
            return self.aggregate(func, *args, **kwds)

        # if we are a string, try to dispatch
        if isinstance(func, str):
            return self._try_aggregate_string_function(func, *args, **kwds)

        # handle ufuncs and lambdas
        if kwds or args and not isinstance(func, np.ufunc):

            def f(x):
                return func(x, *args, **kwds)

        else:
            f = func

        with np.errstate(all="ignore"):
            if isinstance(f, np.ufunc):
                return f(self)

            # row-wise access
            if is_extension_array_dtype(self.dtype) and hasattr(self._values, "map"):
                # GH#23179 some EAs do not have `map`
                mapped = self._values.map(f)
            else:
                values = self.astype(object).values
                mapped = lib.map_infer(values, f, convert=convert_dtype)

        if len(mapped) and isinstance(mapped[0], Series):
            # GH 25959 use pd.array instead of tolist
            # so extension arrays can be used
            return self._constructor_expanddim(pd.array(mapped), index=self.index)
        else:
            return self._constructor(mapped, index=self.index).__finalize__(self)

    def _reduce(
        self, op, name, axis=0, skipna=True, numeric_only=None, filter_type=None, **kwds
    ):
        """
        Perform a reduction operation.

        If we have an ndarray as a value, then simply perform the operation,
        otherwise delegate to the object.
        """
        delegate = self._values

        if axis is not None:
            self._get_axis_number(axis)

        if isinstance(delegate, Categorical):
            return delegate._reduce(name, skipna=skipna, **kwds)
        elif isinstance(delegate, ExtensionArray):
            # dispatch to ExtensionArray interface
            return delegate._reduce(name, skipna=skipna, **kwds)
        elif is_datetime64_dtype(delegate):
            # use DatetimeIndex implementation to handle skipna correctly
            delegate = DatetimeIndex(delegate)
        elif is_timedelta64_dtype(delegate) and hasattr(TimedeltaIndex, name):
            # use TimedeltaIndex to handle skipna correctly
            # TODO: remove hasattr check after TimedeltaIndex has `std` method
            delegate = TimedeltaIndex(delegate)

        # dispatch to numpy arrays
        elif isinstance(delegate, np.ndarray):
            if numeric_only:
                raise NotImplementedError(
                    f"Series.{name} does not implement numeric_only."
                )
            with np.errstate(all="ignore"):
                return op(delegate, skipna=skipna, **kwds)

        # TODO(EA) dispatch to Index
        # remove once all internals extension types are
        # moved to ExtensionArrays
        return delegate._reduce(
            op=op,
            name=name,
            axis=axis,
            skipna=skipna,
            numeric_only=numeric_only,
            filter_type=filter_type,
            **kwds,
        )

    def _reindex_indexer(self, new_index, indexer, copy):
        if indexer is None:
            if copy:
                return self.copy()
            return self

        new_values = algorithms.take_1d(
            self._values, indexer, allow_fill=True, fill_value=None
        )
        return self._constructor(new_values, index=new_index)

    def _needs_reindex_multi(self, axes, method, level):
        """
        Check if we do need a multi reindex; this is for compat with
        higher dims.
        """
        return False

    @Appender(generic._shared_docs["align"] % _shared_doc_kwargs)
    def align(
        self,
        other,
        join="outer",
        axis=None,
        level=None,
        copy=True,
        fill_value=None,
        method=None,
        limit=None,
        fill_axis=0,
        broadcast_axis=None,
    ):
        return super().align(
            other,
            join=join,
            axis=axis,
            level=level,
            copy=copy,
            fill_value=fill_value,
            method=method,
            limit=limit,
            fill_axis=fill_axis,
            broadcast_axis=broadcast_axis,
        )

    def rename(
        self,
        index=None,
        *,
        axis=None,
        copy=True,
        inplace=False,
        level=None,
        errors="ignore",
    ):
        """
        Alter Series index labels or name.

        Function / dict values must be unique (1-to-1). Labels not contained in
        a dict / Series will be left as-is. Extra labels listed don't throw an
        error.

        Alternatively, change ``Series.name`` with a scalar value.

        See the :ref:`user guide <basics.rename>` for more.

        Parameters
        ----------
        axis : {0 or "index"}
            Unused. Accepted for compatability with DataFrame method only.
        index : scalar, hashable sequence, dict-like or function, optional
            Functions or dict-like are transformations to apply to
            the index.
            Scalar or hashable sequence-like will alter the ``Series.name``
            attribute.

        **kwargs
            Additional keyword arguments passed to the function. Only the
            "inplace" keyword is used.

        Returns
        -------
        Series
            Series with index labels or name altered.

        See Also
        --------
        DataFrame.rename : Corresponding DataFrame method.
        Series.rename_axis : Set the name of the axis.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s
        0    1
        1    2
        2    3
        dtype: int64
        >>> s.rename("my_name")  # scalar, changes Series.name
        0    1
        1    2
        2    3
        Name: my_name, dtype: int64
        >>> s.rename(lambda x: x ** 2)  # function, changes labels
        0    1
        1    2
        4    3
        dtype: int64
        >>> s.rename({1: 3, 2: 5})  # mapping, changes labels
        0    1
        3    2
        5    3
        dtype: int64
        """
        if callable(index) or is_dict_like(index):
            return super().rename(
                index, copy=copy, inplace=inplace, level=level, errors=errors
            )
        else:
            return self._set_name(index, inplace=inplace)

    @Appender(
        """
        >>> s = pd.Series([1, 2, 3])
        >>> s
        0    1
        1    2
        2    3
        dtype: int64

        >>> s.set_axis(['a', 'b', 'c'], axis=0)
        a    1
        b    2
        c    3
        dtype: int64
    """
    )
    @Substitution(
        **_shared_doc_kwargs,
        extended_summary_sub="",
        axis_description_sub="",
        see_also_sub="",
    )
    @Appender(generic.NDFrame.set_axis.__doc__)
    def set_axis(self, labels, axis=0, inplace=False):
        return super().set_axis(labels, axis=axis, inplace=inplace)

    @Substitution(**_shared_doc_kwargs)
    @Appender(generic.NDFrame.reindex.__doc__)
    def reindex(self, index=None, **kwargs):
        return super().reindex(index=index, **kwargs)

    def drop(
        self,
        labels=None,
        axis=0,
        index=None,
        columns=None,
        level=None,
        inplace=False,
        errors="raise",
    ) -> "Series":
        """
        Return Series with specified index labels removed.

        Remove elements of a Series based on specifying the index labels.
        When using a multi-index, labels on different levels can be removed
        by specifying the level.

        Parameters
        ----------
        labels : single label or list-like
            Index labels to drop.
        axis : 0, default 0
            Redundant for application on Series.
        index : single label or list-like
            Redundant for application on Series, but 'index' can be used instead
            of 'labels'.

            .. versionadded:: 0.21.0
        columns : single label or list-like
            No change is made to the Series; use 'index' or 'labels' instead.

            .. versionadded:: 0.21.0
        level : int or level name, optional
            For MultiIndex, level for which the labels will be removed.
        inplace : bool, default False
            If True, do operation inplace and return None.
        errors : {'ignore', 'raise'}, default 'raise'
            If 'ignore', suppress error and only existing labels are dropped.

        Returns
        -------
        Series
            Series with specified index labels removed.

        Raises
        ------
        KeyError
            If none of the labels are found in the index.

        See Also
        --------
        Series.reindex : Return only specified index labels of Series.
        Series.dropna : Return series without null values.
        Series.drop_duplicates : Return Series with duplicate values removed.
        DataFrame.drop : Drop specified labels from rows or columns.

        Examples
        --------
        >>> s = pd.Series(data=np.arange(3), index=['A', 'B', 'C'])
        >>> s
        A  0
        B  1
        C  2
        dtype: int64

        Drop labels B en C

        >>> s.drop(labels=['B', 'C'])
        A  0
        dtype: int64

        Drop 2nd level label in MultiIndex Series

        >>> midx = pd.MultiIndex(levels=[['lama', 'cow', 'falcon'],
        ...                              ['speed', 'weight', 'length']],
        ...                      codes=[[0, 0, 0, 1, 1, 1, 2, 2, 2],
        ...                             [0, 1, 2, 0, 1, 2, 0, 1, 2]])
        >>> s = pd.Series([45, 200, 1.2, 30, 250, 1.5, 320, 1, 0.3],
        ...               index=midx)
        >>> s
        lama    speed      45.0
                weight    200.0
                length      1.2
        cow     speed      30.0
                weight    250.0
                length      1.5
        falcon  speed     320.0
                weight      1.0
                length      0.3
        dtype: float64

        >>> s.drop(labels='weight', level=1)
        lama    speed      45.0
                length      1.2
        cow     speed      30.0
                length      1.5
        falcon  speed     320.0
                length      0.3
        dtype: float64
        """
        return super().drop(
            labels=labels,
            axis=axis,
            index=index,
            columns=columns,
            level=level,
            inplace=inplace,
            errors=errors,
        )

    @Substitution(**_shared_doc_kwargs)
    @Appender(generic.NDFrame.fillna.__doc__)
    def fillna(
        self,
        value=None,
        method=None,
        axis=None,
        inplace=False,
        limit=None,
        downcast=None,
    ) -> Optional["Series"]:
        return super().fillna(
            value=value,
            method=method,
            axis=axis,
            inplace=inplace,
            limit=limit,
            downcast=downcast,
        )

    @Appender(generic._shared_docs["replace"] % _shared_doc_kwargs)
    def replace(
        self,
        to_replace=None,
        value=None,
        inplace=False,
        limit=None,
        regex=False,
        method="pad",
    ):
        return super().replace(
            to_replace=to_replace,
            value=value,
            inplace=inplace,
            limit=limit,
            regex=regex,
            method=method,
        )

    @Appender(generic._shared_docs["shift"] % _shared_doc_kwargs)
    def shift(self, periods=1, freq=None, axis=0, fill_value=None) -> "Series":
        return super().shift(
            periods=periods, freq=freq, axis=axis, fill_value=fill_value
        )

    def memory_usage(self, index=True, deep=False):
        """
        Return the memory usage of the Series.

        The memory usage can optionally include the contribution of
        the index and of elements of `object` dtype.

        Parameters
        ----------
        index : bool, default True
            Specifies whether to include the memory usage of the Series index.
        deep : bool, default False
            If True, introspect the data deeply by interrogating
            `object` dtypes for system-level memory consumption, and include
            it in the returned value.

        Returns
        -------
        int
            Bytes of memory consumed.

        See Also
        --------
        numpy.ndarray.nbytes : Total bytes consumed by the elements of the
            array.
        DataFrame.memory_usage : Bytes consumed by a DataFrame.

        Examples
        --------
        >>> s = pd.Series(range(3))
        >>> s.memory_usage()
        152

        Not including the index gives the size of the rest of the data, which
        is necessarily smaller:

        >>> s.memory_usage(index=False)
        24

        The memory footprint of `object` values is ignored by default:

        >>> s = pd.Series(["a", "b"])
        >>> s.values
        array(['a', 'b'], dtype=object)
        >>> s.memory_usage()
        144
        >>> s.memory_usage(deep=True)
        260
        """
        v = super().memory_usage(deep=deep)
        if index:
            v += self.index.memory_usage(deep=deep)
        return v

    def isin(self, values) -> "Series":
        """
        Check whether `values` are contained in Series.

        Return a boolean Series showing whether each element in the Series
        matches an element in the passed sequence of `values` exactly.

        Parameters
        ----------
        values : set or list-like
            The sequence of values to test. Passing in a single string will
            raise a ``TypeError``. Instead, turn a single string into a
            list of one element.

        Returns
        -------
        Series
            Series of booleans indicating if each element is in values.

        Raises
        ------
        TypeError
          * If `values` is a string

        See Also
        --------
        DataFrame.isin : Equivalent method on DataFrame.

        Examples
        --------
        >>> s = pd.Series(['lama', 'cow', 'lama', 'beetle', 'lama',
        ...                'hippo'], name='animal')
        >>> s.isin(['cow', 'lama'])
        0     True
        1     True
        2     True
        3    False
        4     True
        5    False
        Name: animal, dtype: bool

        Passing a single string as ``s.isin('lama')`` will raise an error. Use
        a list of one element instead:

        >>> s.isin(['lama'])
        0     True
        1    False
        2     True
        3    False
        4     True
        5    False
        Name: animal, dtype: bool
        """
        result = algorithms.isin(self, values)
        return self._constructor(result, index=self.index).__finalize__(self)

    def between(self, left, right, inclusive=True) -> "Series":
        """
        Return boolean Series equivalent to left <= series <= right.

        This function returns a boolean vector containing `True` wherever the
        corresponding Series element is between the boundary values `left` and
        `right`. NA values are treated as `False`.

        Parameters
        ----------
        left : scalar or list-like
            Left boundary.
        right : scalar or list-like
            Right boundary.
        inclusive : bool, default True
            Include boundaries.

        Returns
        -------
        Series
            Series representing whether each element is between left and
            right (inclusive).

        See Also
        --------
        Series.gt : Greater than of series and other.
        Series.lt : Less than of series and other.

        Notes
        -----
        This function is equivalent to ``(left <= ser) & (ser <= right)``

        Examples
        --------
        >>> s = pd.Series([2, 0, 4, 8, np.nan])

        Boundary values are included by default:

        >>> s.between(1, 4)
        0     True
        1    False
        2     True
        3    False
        4    False
        dtype: bool

        With `inclusive` set to ``False`` boundary values are excluded:

        >>> s.between(1, 4, inclusive=False)
        0     True
        1    False
        2    False
        3    False
        4    False
        dtype: bool

        `left` and `right` can be any scalar value:

        >>> s = pd.Series(['Alice', 'Bob', 'Carol', 'Eve'])
        >>> s.between('Anna', 'Daniel')
        0    False
        1     True
        2     True
        3    False
        dtype: bool
        """
        if inclusive:
            lmask = self >= left
            rmask = self <= right
        else:
            lmask = self > left
            rmask = self < right

        return lmask & rmask

    # ----------------------------------------------------------------------
    # Convert to types that support pd.NA

    def _convert_dtypes(
        self: ABCSeries,
        infer_objects: bool = True,
        convert_string: bool = True,
        convert_integer: bool = True,
        convert_boolean: bool = True,
    ) -> "Series":
        input_series = self
        if infer_objects:
            input_series = input_series.infer_objects()
            if is_object_dtype(input_series):
                input_series = input_series.copy()

        if convert_string or convert_integer or convert_boolean:
            inferred_dtype = convert_dtypes(
                input_series._values, convert_string, convert_integer, convert_boolean
            )
            try:
                result = input_series.astype(inferred_dtype)
            except TypeError:
                result = input_series.copy()
        else:
            result = input_series.copy()
        return result

    @Appender(generic._shared_docs["isna"] % _shared_doc_kwargs)
    def isna(self) -> "Series":
        return super().isna()

    @Appender(generic._shared_docs["isna"] % _shared_doc_kwargs)
    def isnull(self) -> "Series":
        return super().isnull()

    @Appender(generic._shared_docs["notna"] % _shared_doc_kwargs)
    def notna(self) -> "Series":
        return super().notna()

    @Appender(generic._shared_docs["notna"] % _shared_doc_kwargs)
    def notnull(self) -> "Series":
        return super().notnull()

    def dropna(self, axis=0, inplace=False, how=None):
        """
        Return a new Series with missing values removed.

        See the :ref:`User Guide <missing_data>` for more on which values are
        considered missing, and how to work with missing data.

        Parameters
        ----------
        axis : {0 or 'index'}, default 0
            There is only one axis to drop values from.
        inplace : bool, default False
            If True, do operation inplace and return None.
        how : str, optional
            Not in use. Kept for compatibility.

        Returns
        -------
        Series
            Series with NA entries dropped from it.

        See Also
        --------
        Series.isna: Indicate missing values.
        Series.notna : Indicate existing (non-missing) values.
        Series.fillna : Replace missing values.
        DataFrame.dropna : Drop rows or columns which contain NA values.
        Index.dropna : Drop missing indices.

        Examples
        --------
        >>> ser = pd.Series([1., 2., np.nan])
        >>> ser
        0    1.0
        1    2.0
        2    NaN
        dtype: float64

        Drop NA values from a Series.

        >>> ser.dropna()
        0    1.0
        1    2.0
        dtype: float64

        Keep the Series with valid entries in the same variable.

        >>> ser.dropna(inplace=True)
        >>> ser
        0    1.0
        1    2.0
        dtype: float64

        Empty strings are not considered NA values. ``None`` is considered an
        NA value.

        >>> ser = pd.Series([np.NaN, 2, pd.NaT, '', None, 'I stay'])
        >>> ser
        0       NaN
        1         2
        2       NaT
        3
        4      None
        5    I stay
        dtype: object
        >>> ser.dropna()
        1         2
        3
        5    I stay
        dtype: object
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        # Validate the axis parameter
        self._get_axis_number(axis or 0)

        if self._can_hold_na:
            result = remove_na_arraylike(self)
            if inplace:
                self._update_inplace(result)
            else:
                return result
        else:
            if inplace:
                # do nothing
                pass
            else:
                return self.copy()

    # ----------------------------------------------------------------------
    # Time series-oriented methods

    def to_timestamp(self, freq=None, how="start", copy=True) -> "Series":
        """
        Cast to DatetimeIndex of Timestamps, at *beginning* of period.

        Parameters
        ----------
        freq : str, default frequency of PeriodIndex
            Desired frequency.
        how : {'s', 'e', 'start', 'end'}
            Convention for converting period to timestamp; start of period
            vs. end.
        copy : bool, default True
            Whether or not to return a copy.

        Returns
        -------
        Series with DatetimeIndex
        """
        new_values = self._values
        if copy:
            new_values = new_values.copy()

        assert isinstance(self.index, (ABCDatetimeIndex, ABCPeriodIndex))
        new_index = self.index.to_timestamp(freq=freq, how=how)
        return self._constructor(new_values, index=new_index).__finalize__(self)

    def to_period(self, freq=None, copy=True) -> "Series":
        """
        Convert Series from DatetimeIndex to PeriodIndex with desired
        frequency (inferred from index if not passed).

        Parameters
        ----------
        freq : str, default None
            Frequency associated with the PeriodIndex.
        copy : bool, default True
            Whether or not to return a copy.

        Returns
        -------
        Series
            Series with index converted to PeriodIndex.
        """
        new_values = self._values
        if copy:
            new_values = new_values.copy()

        assert isinstance(self.index, ABCDatetimeIndex)
        new_index = self.index.to_period(freq=freq)
        return self._constructor(new_values, index=new_index).__finalize__(self)

    # ----------------------------------------------------------------------
    # Add index
    index: "Index" = properties.AxisProperty(
        axis=0, doc="The index (axis labels) of the Series."
    )

    # ----------------------------------------------------------------------
    # Accessor Methods
    # ----------------------------------------------------------------------
    str = CachedAccessor("str", StringMethods)
    dt = CachedAccessor("dt", CombinedDatetimelikeProperties)
    cat = CachedAccessor("cat", CategoricalAccessor)
    plot = CachedAccessor("plot", pandas.plotting.PlotAccessor)
    sparse = CachedAccessor("sparse", SparseAccessor)

    # ----------------------------------------------------------------------
    # Add plotting methods to Series
    hist = pandas.plotting.hist_series


Series._setup_axes(["index"], docs={"index": "The index (axis labels) of the Series."})
Series._add_numeric_operations()
Series._add_series_or_dataframe_operations()

# Add arithmetic!
ops.add_flex_arithmetic_methods(Series)
ops.add_special_arithmetic_methods(Series)
