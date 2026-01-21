"""
Data structure for 1-dimensional cross-sectional and time series data
"""

from __future__ import annotations

from collections.abc import (
    Callable,
    Hashable,
    Iterable,
    Mapping,
    Sequence,
)
import functools
import operator
import sys
from textwrap import dedent
from typing import (
    IO,
    TYPE_CHECKING,
    Any,
    Literal,
    Self,
    cast,
    overload,
)
import warnings

import numpy as np

from pandas._libs import (
    lib,
    properties,
    reshape,
)
from pandas._libs.lib import is_range_indexer
from pandas.compat import CHAINED_WARNING_DISABLED
from pandas.compat._constants import (
    REF_COUNT,
    REF_COUNT_METHOD,
)
from pandas.compat._optional import import_optional_dependency
from pandas.compat.numpy import function as nv
from pandas.errors import (
    ChainedAssignmentError,
    InvalidIndexError,
    Pandas4Warning,
)
from pandas.errors.cow import (
    _chained_assignment_method_update_msg,
    _chained_assignment_msg,
)
from pandas.util._decorators import (
    Appender,
    deprecate_nonkeyword_arguments,
    doc,
    set_module,
)
from pandas.util._exceptions import (
    find_stack_level,
)
from pandas.util._validators import (
    validate_ascending,
    validate_bool_kwarg,
    validate_percentile,
)

from pandas.core.dtypes.astype import astype_is_view
from pandas.core.dtypes.cast import (
    LossySetitemError,
    construct_1d_arraylike_from_scalar,
    find_common_type,
    infer_dtype_from,
    maybe_box_native,
    maybe_unbox_numpy_scalar,
)
from pandas.core.dtypes.common import (
    is_dict_like,
    is_float,
    is_integer,
    is_iterator,
    is_list_like,
    is_object_dtype,
    is_scalar,
    pandas_dtype,
    validate_all_hashable,
)
from pandas.core.dtypes.dtypes import (
    ExtensionDtype,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCSeries,
)
from pandas.core.dtypes.inference import is_hashable
from pandas.core.dtypes.missing import (
    isna,
    na_value_for_dtype,
    notna,
    remove_na_arraylike,
)

from pandas.core import (
    algorithms,
    base,
    common as com,
    nanops,
    ops,
    roperator,
)
from pandas.core.accessor import Accessor
from pandas.core.apply import SeriesApply
from pandas.core.arrays import ExtensionArray
from pandas.core.arrays.arrow import (
    ListAccessor,
    StructAccessor,
)
from pandas.core.arrays.categorical import CategoricalAccessor
from pandas.core.arrays.sparse import SparseAccessor
from pandas.core.construction import (
    array as pd_array,
    extract_array,
    sanitize_array,
)
from pandas.core.generic import NDFrame
from pandas.core.indexers import (
    disallow_ndim_indexing,
    unpack_1tuple,
)
from pandas.core.indexes.accessors import CombinedDatetimelikeProperties
from pandas.core.indexes.api import (
    DatetimeIndex,
    Index,
    MultiIndex,
    PeriodIndex,
    default_index,
    ensure_index,
    maybe_sequence_to_range,
)
import pandas.core.indexes.base as ibase
from pandas.core.indexes.multi import maybe_droplevels
from pandas.core.indexing import (
    check_bool_indexer,
    check_dict_or_set_indexers,
)
from pandas.core.internals import SingleBlockManager
from pandas.core.methods import selectn
from pandas.core.shared_docs import _shared_docs
from pandas.core.sorting import (
    ensure_key_mapped,
    nargsort,
)
from pandas.core.strings.accessor import StringMethods
from pandas.core.tools.datetimes import to_datetime

import pandas.io.formats.format as fmt
from pandas.io.formats.info import (
    SeriesInfo,
)
import pandas.plotting

if TYPE_CHECKING:
    from pandas._libs.internals import BlockValuesRefs
    from pandas._typing import (
        AggFuncType,
        AnyAll,
        AnyArrayLike,
        ArrayLike,
        ArrowArrayExportable,
        ArrowStreamExportable,
        Axis,
        AxisInt,
        CorrelationMethod,
        DropKeep,
        Dtype,
        DtypeObj,
        FilePath,
        Frequency,
        IgnoreRaise,
        IndexKeyFunc,
        IndexLabel,
        Level,
        ListLike,
        MutableMappingT,
        NaPosition,
        NumpySorter,
        NumpyValueArrayLike,
        QuantileInterpolation,
        ReindexMethod,
        Renamer,
        Scalar,
        SortKind,
        StorageOptions,
        Suffixes,
        ValueKeyFunc,
        WriteBuffer,
        npt,
    )

    from pandas.core.frame import DataFrame
    from pandas.core.groupby.generic import SeriesGroupBy

__all__ = ["Series"]

_shared_doc_kwargs = {
    "axes": "index",
    "klass": "Series",
    "axes_single_arg": "{0 or 'index'}",
    "axis": """axis : {0 or 'index'}
        Unused. Parameter needed for compatibility with DataFrame.""",
    "inplace": """inplace : bool, default False
        If True, performs operation inplace and returns None.""",
    "unique": "np.ndarray",
    "duplicated": "Series",
    "optional_by": "",
    "optional_reindex": """
index : array-like, optional
    New labels for the index. Preferably an Index object to avoid
    duplicating data.
axis : int or str, optional
    Unused.""",
}

# ----------------------------------------------------------------------
# Series class


# error: Cannot override final attribute "ndim" (previously declared in base
# class "NDFrame")
# error: Cannot override final attribute "size" (previously declared in base
# class "NDFrame")
# definition in base class "NDFrame"
@set_module("pandas")
class Series(base.IndexOpsMixin, NDFrame):  # type: ignore[misc]
    """
    One-dimensional ndarray with axis labels (including time series).

    Labels need not be unique but must be a hashable type. The object
    supports both integer- and label-based indexing and provides a host of
    methods for performing operations involving the index. Statistical
    methods from ndarray have been overridden to automatically exclude
    missing data (currently represented as NaN).

    Operations between Series (+, -, /, \\*, \\*\\*) align values based on their
    associated index values-- they need not be the same length. The result
    index will be the sorted union of the two indexes.

    Parameters
    ----------
    data : array-like, Iterable, dict, or scalar value
        Contains data stored in Series. If data is a dict, argument order is
        maintained. Unordered sets are not supported.
    index : array-like or Index (1d)
        Values must be hashable and have the same length as `data`.
        Non-unique index values are allowed. Will default to
        RangeIndex (0, 1, 2, ..., n) if not provided. If data is dict-like
        and index is None, then the keys in the data are used as the index. If the
        index is not None, the resulting Series is reindexed with the index values.
    dtype : str, numpy.dtype, or ExtensionDtype, optional
        Data type for the output Series. If not specified, this will be
        inferred from `data`.
        See the :ref:`user guide <basics.dtypes>` for more usages.
    name : Hashable, default None
        The name to give to the Series.
    copy : bool, default None
        Whether to copy input data, only relevant for array, Series, and Index
        inputs (for other input, e.g. a list, a new array is created anyway).
        Defaults to True for array input and False for Index/Series.
        Even when False for Index/Series, a shallow copy of the data is made.
        Set to False to avoid copying array input at your own risk (if you
        know the input data won't be modified elsewhere).
        Set to True to force copying Series/Index input up front.

    See Also
    --------
    DataFrame : Two-dimensional, size-mutable, potentially heterogeneous tabular data.
    Index : Immutable sequence used for indexing and alignment.

    Notes
    -----
    Please reference the :ref:`User Guide <basics.series>` for more information.

    Examples
    --------
    Constructing Series from a dictionary with an Index specified

    >>> d = {"a": 1, "b": 2, "c": 3}
    >>> ser = pd.Series(data=d, index=["a", "b", "c"])
    >>> ser
    a   1
    b   2
    c   3
    dtype: int64

    The keys of the dictionary match with the Index values, hence the Index
    values have no effect.

    >>> d = {"a": 1, "b": 2, "c": 3}
    >>> ser = pd.Series(data=d, index=["x", "y", "z"])
    >>> ser
    x   NaN
    y   NaN
    z   NaN
    dtype: float64

    Note that the Index is first built with the keys from the dictionary.
    After this the Series is reindexed with the given Index values, hence we
    get all NaN as a result.

    Constructing Series from a list with `copy=False`.

    >>> r = [1, 2]
    >>> ser = pd.Series(r, copy=False)
    >>> ser.iloc[0] = 999
    >>> r
    [1, 2]
    >>> ser
    0    999
    1      2
    dtype: int64

    Due to input data type the Series has a `copy` of
    the original data even though `copy=False`, so
    the data is unchanged.

    Constructing Series from a 1d ndarray with `copy=False`.

    >>> r = np.array([1, 2])
    >>> ser = pd.Series(r, copy=False)
    >>> ser.iloc[0] = 999
    >>> r
    array([999,   2])
    >>> ser
    0    999
    1      2
    dtype: int64

    Due to input data type the Series has a `view` on
    the original data, so
    the data is changed as well.
    """

    _typ = "series"
    _HANDLED_TYPES = (Index, ExtensionArray, np.ndarray)

    _name: Hashable
    _metadata: list[str] = ["_name"]
    _internal_names_set = {"index", "name"} | NDFrame._internal_names_set
    _accessors = {"dt", "cat", "str", "sparse"}
    _hidden_attrs = (
        base.IndexOpsMixin._hidden_attrs | NDFrame._hidden_attrs | frozenset([])
    )

    # similar to __array_priority__, positions Series after DataFrame
    #  but before Index and ExtensionArray.  Should NOT be overridden by subclasses.
    __pandas_priority__ = 3000

    # Override cache_readonly bc Series is mutable
    hasnans = property(
        # error: "Callable[[IndexOpsMixin], bool]" has no attribute "fget"
        base.IndexOpsMixin.hasnans.fget,  # type: ignore[attr-defined]
        doc=base.IndexOpsMixin.hasnans.__doc__,
    )
    _mgr: SingleBlockManager

    # ----------------------------------------------------------------------
    # Constructors

    def __init__(
        self,
        data=None,
        index=None,
        dtype: Dtype | None = None,
        name=None,
        copy: bool | None = None,
    ) -> None:
        allow_mgr = False
        if (
            isinstance(data, SingleBlockManager)
            and index is None
            and dtype is None
            and (copy is False or copy is None)
        ):
            if not allow_mgr:
                # GH#52419
                warnings.warn(
                    f"Passing a {type(data).__name__} to {type(self).__name__} "
                    "is deprecated and will raise in a future version. "
                    "Use public APIs instead.",
                    Pandas4Warning,
                    stacklevel=2,
                )
            data = data.copy(deep=False)
            # GH#33357 called with just the SingleBlockManager
            NDFrame.__init__(self, data)
            self.name = name
            return

        if isinstance(data, (ExtensionArray, np.ndarray)):
            if copy is not False:
                if dtype is None or astype_is_view(data.dtype, pandas_dtype(dtype)):
                    data = data.copy()
                    copy = False
        if copy is None:
            copy = False

        if isinstance(data, SingleBlockManager) and not copy:
            data = data.copy(deep=False)

            if not allow_mgr:
                warnings.warn(
                    f"Passing a {type(data).__name__} to {type(self).__name__} "
                    "is deprecated and will raise in a future version. "
                    "Use public APIs instead.",
                    Pandas4Warning,
                    stacklevel=2,
                )
                allow_mgr = True

        name = ibase.maybe_extract_name(name, data, type(self))

        if index is not None:
            index = ensure_index(index)

        if dtype is not None:
            dtype = self._validate_dtype(dtype)

        if data is None:
            index = index if index is not None else default_index(0)
            if len(index) or dtype is not None:
                data = na_value_for_dtype(pandas_dtype(dtype), compat=False)
            else:
                data = []

        if isinstance(data, MultiIndex):
            raise NotImplementedError(
                "initializing a Series from a MultiIndex is not supported"
            )

        refs = None
        if isinstance(data, Index):
            if dtype is not None:
                data = data.astype(dtype)
            if not copy:
                refs = data._references

        elif isinstance(data, np.ndarray):
            if len(data.dtype):
                # GH#13296 we are dealing with a compound dtype, which
                #  should be treated as 2D
                raise ValueError(
                    "Cannot construct a Series from an ndarray with "
                    "compound dtype.  Use DataFrame instead."
                )
        elif isinstance(data, Series):
            if index is None:
                index = data.index
                data = data._mgr.copy(deep=False)
            else:
                data = data.reindex(index)
                data = data._mgr
                if data._has_no_reference(0):
                    copy = False
        elif isinstance(data, Mapping):
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

            if not allow_mgr:
                warnings.warn(
                    f"Passing a {type(data).__name__} to {type(self).__name__} "
                    "is deprecated and will raise in a future version. "
                    "Use public APIs instead.",
                    Pandas4Warning,
                    stacklevel=2,
                )
                allow_mgr = True

        elif isinstance(data, ExtensionArray):
            pass
        else:
            data = com.maybe_iterable_to_list(data)
            if is_list_like(data) and not len(data) and dtype is None:
                # GH 29405: Pre-2.0, this defaulted to float.
                dtype = np.dtype(object)

        if index is None:
            if not is_list_like(data):
                data = [data]
            index = default_index(len(data))
        elif is_list_like(data):
            com.require_length_match(data, index)

        # create/copy the manager
        if isinstance(data, SingleBlockManager):
            if dtype is not None:
                if not astype_is_view(data.dtype, pandas_dtype(dtype)):
                    copy = False
                data = data.astype(dtype=dtype)
            if copy:
                data = data.copy(deep=True)
        else:
            data = sanitize_array(data, index, dtype, copy)
            data = SingleBlockManager.from_array(data, index, refs=refs)

        NDFrame.__init__(self, data)
        self.name = name
        self._set_axis(0, index)

    def _init_dict(
        self, data: Mapping, index: Index | None = None, dtype: DtypeObj | None = None
    ):
        """
        Derive the "_mgr" and "index" attributes of a new Series from a
        dictionary input.

        Parameters
        ----------
        data : dict or dict-like
            Data used to populate the new Series.
        index : Index or None, default None
            Index for the new Series: if None, use dict keys.
        dtype : np.dtype, ExtensionDtype, or None, default None
            The dtype for the new Series: if None, infer from data.

        Returns
        -------
        _data : BlockManager for the new Series
        index : index for the new Series
        """
        # Looking for NaN in dict doesn't work ({np.nan : 1}[float('nan')]
        # raises KeyError), so we iterate the entire dict, and align
        if data:
            # GH:34717, issue was using zip to extract key and values from data.
            # using generators in effects the performance.
            # Below is the new way of extracting the keys and values

            keys = maybe_sequence_to_range(tuple(data.keys()))
            values = list(data.values())  # Generating list of values- faster way
        elif index is not None:
            # fastpath for Series(data=None). Just use broadcasting a scalar
            # instead of reindexing.
            if len(index) or dtype is not None:
                values = na_value_for_dtype(pandas_dtype(dtype), compat=False)
            else:
                values = []
            keys = index
        else:
            keys, values = default_index(0), []

        # Input is now list-like, so rely on "standard" construction:
        s = Series(values, index=keys, dtype=dtype)

        # Now we just make sure the order is respected, if any
        if data and index is not None:
            s = s.reindex(index)
        return s._mgr, s.index

    # ----------------------------------------------------------------------

    def __arrow_c_stream__(self, requested_schema=None):
        """
        Export the pandas Series as an Arrow C stream PyCapsule.

        This relies on pyarrow to convert the pandas Series to the Arrow
        format (and follows the default behavior of ``pyarrow.Array.from_pandas``
        in its handling of the index, i.e. to ignore it).
        This conversion is not necessarily zero-copy.

        Parameters
        ----------
        requested_schema : PyCapsule, default None
            The schema to which the dataframe should be casted, passed as a
            PyCapsule containing a C ArrowSchema representation of the
            requested schema.

        Returns
        -------
        PyCapsule
        """
        pa = import_optional_dependency("pyarrow", min_version="16.0.0")
        type = (
            pa.DataType._import_from_c_capsule(requested_schema)
            if requested_schema is not None
            else None
        )
        ca = pa.array(self, type=type)
        if not isinstance(ca, pa.ChunkedArray):
            ca = pa.chunked_array([ca])
        return ca.__arrow_c_stream__()

    # ----------------------------------------------------------------------

    @property
    def _constructor(self) -> type[Series]:
        return Series

    def _constructor_from_mgr(self, mgr, axes):
        ser = Series._from_mgr(mgr, axes=axes)
        ser._name = None  # caller is responsible for setting real name

        if type(self) is Series:
            # This would also work `if self._constructor is Series`, but
            #  this check is slightly faster, benefiting the most-common case.
            return ser

        # We assume that the subclass __init__ knows how to handle a
        #  pd.Series object.
        return self._constructor(ser)

    @property
    def _constructor_expanddim(self) -> Callable[..., DataFrame]:
        """
        Used when a manipulation result has one higher dimension as the
        original, such as Series.to_frame()
        """
        from pandas.core.frame import DataFrame

        return DataFrame

    def _constructor_expanddim_from_mgr(self, mgr, axes):
        from pandas.core.frame import DataFrame

        df = DataFrame._from_mgr(mgr, axes=mgr.axes)

        if type(self) is Series:
            # This would also work `if self._constructor_expanddim is DataFrame`,
            #  but this check is slightly faster, benefiting the most-common case.
            return df

        # We assume that the subclass __init__ knows how to handle a
        #  pd.DataFrame object.
        return self._constructor_expanddim(df)

    # types
    @property
    def _can_hold_na(self) -> bool:
        return self._mgr._can_hold_na

    # ndarray compatibility
    @property
    def dtype(self) -> DtypeObj:
        """
        Return the dtype object of the underlying data.

        See Also
        --------
        Series.dtypes : Return the dtype object of the underlying data.
        Series.astype : Cast a pandas object to a specified dtype dtype.
        Series.convert_dtypes : Convert columns to the best possible dtypes using dtypes
            supporting pd.NA.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s.dtype
        dtype('int64')
        """
        return self._mgr.dtype

    @property
    def dtypes(self) -> DtypeObj:
        """
        Return the dtype object of the underlying data.

        See Also
        --------
        DataFrame.dtypes :  Return the dtypes in the DataFrame.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s.dtypes
        dtype('int64')
        """
        # DataFrame compatibility
        return self.dtype

    @property
    def name(self) -> Hashable:
        """
        Return the name of the Series.

        The name of a Series becomes its index or column name if it is used
        to form a DataFrame. It is also used whenever displaying the Series
        using the interpreter.

        Returns
        -------
        label (hashable object)
            The name of the Series, also the column name if part of a DataFrame.

        See Also
        --------
        Series.rename : Sets the Series name when given a scalar input.
        Index.name : Corresponding Index property.

        Examples
        --------
        The Series name can be set initially when calling the constructor.

        >>> s = pd.Series([1, 2, 3], dtype=np.int64, name="Numbers")
        >>> s
        0    1
        1    2
        2    3
        Name: Numbers, dtype: int64
        >>> s.name = "Integers"
        >>> s
        0    1
        1    2
        2    3
        Name: Integers, dtype: int64

        The name of a Series within a DataFrame is its column name.

        >>> df = pd.DataFrame(
        ...     [[1, 2], [3, 4], [5, 6]], columns=["Odd Numbers", "Even Numbers"]
        ... )
        >>> df
           Odd Numbers  Even Numbers
        0            1             2
        1            3             4
        2            5             6
        >>> df["Even Numbers"].name
        'Even Numbers'
        """
        return self._name

    @name.setter
    def name(self, value: Hashable) -> None:
        validate_all_hashable(value, error_name=f"{type(self).__name__}.name")
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

        >>> pd.Series(list("aabc")).values
        <ArrowStringArray>
        ['a', 'a', 'b', 'c']
        Length: 4, dtype: str

        >>> pd.Series(list("aabc")).astype("category").values
        ['a', 'a', 'b', 'c']
        Categories (3, str): ['a', 'b', 'c']

        Timezone aware datetime data is converted to UTC:

        >>> pd.Series(pd.date_range("20130101", periods=3, tz="US/Eastern")).values
        array(['2013-01-01T05:00:00.000000',
               '2013-01-02T05:00:00.000000',
               '2013-01-03T05:00:00.000000'], dtype='datetime64[us]')
        """
        return self._mgr.external_values()

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

        Overview:

        dtype       | values        | _values       | array                 |
        ----------- | ------------- | ------------- | --------------------- |
        Numeric     | ndarray       | ndarray       | NumpyExtensionArray   |
        Category    | Categorical   | Categorical   | Categorical           |
        dt64[ns]    | ndarray[M8ns] | DatetimeArray | DatetimeArray         |
        dt64[ns tz] | ndarray[M8ns] | DatetimeArray | DatetimeArray         |
        td64[ns]    | ndarray[m8ns] | TimedeltaArray| TimedeltaArray        |
        Period      | ndarray[obj]  | PeriodArray   | PeriodArray           |
        Nullable    | EA            | EA            | EA                    |

        """
        return self._mgr.internal_values()

    @property
    def _references(self) -> BlockValuesRefs:
        return self._mgr._block.refs

    @Appender(base.IndexOpsMixin.array.__doc__)  # type: ignore[prop-decorator]
    @property
    def array(self) -> ExtensionArray:
        arr = self._mgr.array_values()
        # TODO decide on read-only https://github.com/pandas-dev/pandas/issues/63099
        # arr = arr.view()
        # arr._readonly = True
        return arr

    def __len__(self) -> int:
        """
        Return the length of the Series.
        """
        return len(self._mgr)

    # ----------------------------------------------------------------------
    # NDArray Compat
    def __array__(
        self, dtype: npt.DTypeLike | None = None, copy: bool | None = None
    ) -> np.ndarray:
        """
        Return the values as a NumPy array.

        Users should not call this directly. Rather, it is invoked by
        :func:`numpy.array` and :func:`numpy.asarray`.

        Parameters
        ----------
        dtype : str or numpy.dtype, optional
            The dtype to use for the resulting NumPy array. By default,
            the dtype is inferred from the data.

        copy : bool or None, optional
            See :func:`numpy.asarray`.

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

        >>> tzser = pd.Series(pd.date_range("2000", periods=2, tz="CET"))
        >>> np.asarray(tzser, dtype="object")
        array([Timestamp('2000-01-01 00:00:00+0100', tz='CET'),
               Timestamp('2000-01-02 00:00:00+0100', tz='CET')],
              dtype=object)

        Or the values may be localized to UTC and the tzinfo discarded with
        ``dtype='datetime64[ns]'``

        >>> np.asarray(tzser, dtype="datetime64[ns]")  # doctest: +ELLIPSIS
        array(['1999-12-31T23:00:00.000000000', ...],
              dtype='datetime64[ns]')
        """
        values = self._values
        if copy is None:
            # Note: branch avoids `copy=None` for NumPy 1.x support
            arr = np.asarray(values, dtype=dtype)
        else:
            arr = np.array(values, dtype=dtype, copy=copy)

        if copy is True:
            return arr
        if copy is False or astype_is_view(values.dtype, arr.dtype):
            arr = arr.view()
            arr.flags.writeable = False
        return arr

    # ----------------------------------------------------------------------

    # indexers
    @property
    def axes(self) -> list[Index]:
        """
        Return a list of the row axis labels.
        """
        return [self.index]

    # ----------------------------------------------------------------------
    # Indexing Methods

    def _ixs(self, i: int, axis: AxisInt = 0) -> Any:
        """
        Return the i-th value or values in the Series by location.

        Parameters
        ----------
        i : int

        Returns
        -------
        scalar
        """
        return self._values[i]

    def _slice(self, slobj: slice, axis: AxisInt = 0) -> Series:
        # axis kwarg is retained for compat with NDFrame method
        #  _slice is *always* positional
        mgr = self._mgr.get_slice(slobj, axis=axis)
        out = self._constructor_from_mgr(mgr, axes=mgr.axes)
        out._name = self._name
        return out.__finalize__(self)

    def __getitem__(self, key):
        check_dict_or_set_indexers(key)
        key = com.apply_if_callable(key, self)

        if key is Ellipsis:
            return self.copy(deep=False)

        key_is_scalar = is_scalar(key)
        if isinstance(key, (list, tuple)):
            key = unpack_1tuple(key)

        elif key_is_scalar:
            # Note: GH#50617 in 3.0 we changed int key to always be treated as
            #  a label, matching DataFrame behavior.
            return self._get_value(key)

        # Convert generator to list before going through hashable part
        # (We will iterate through the generator there to check for slices)
        if is_iterator(key):
            key = list(key)

        if is_hashable(key, allow_slice=False):
            # Otherwise index.get_value will raise InvalidIndexError
            try:
                # For labels that don't resolve as scalars like tuples and frozensets
                result = self._get_value(key)

                return result

            except (KeyError, TypeError, InvalidIndexError):
                # InvalidIndexError for e.g. generator
                #  see test_series_getitem_corner_generator
                if isinstance(key, tuple) and isinstance(self.index, MultiIndex):
                    # We still have the corner case where a tuple is a key
                    # in the first level of our MultiIndex
                    return self._get_values_tuple(key)

        if isinstance(key, slice):
            # Do slice check before somewhat-costly is_bool_indexer
            return self._getitem_slice(key)

        if com.is_bool_indexer(key):
            key = check_bool_indexer(self.index, key)
            key = np.asarray(key, dtype=bool)
            return self._get_rows_with_mask(key)

        return self._get_with(key)

    def _get_with(self, key):
        # other: fancy integer or otherwise
        if isinstance(key, ABCDataFrame):
            raise TypeError(
                "Indexing a Series with DataFrame is not "
                "supported, use the appropriate DataFrame column"
            )
        elif isinstance(key, tuple):
            return self._get_values_tuple(key)

        return self.loc[key]

    def _get_values_tuple(self, key: tuple):
        # mpl hackaround
        if com.any_none(*key):
            # mpl compat if we look up e.g. ser[:, np.newaxis];
            #  see tests.series.timeseries.test_mpl_compat_hack
            # the asarray is needed to avoid returning a 2D DatetimeArray
            result = np.asarray(self._values[key])
            disallow_ndim_indexing(result)
            return result

        if not isinstance(self.index, MultiIndex):
            raise KeyError("key of type tuple not found and not a MultiIndex")

        # If key is contained, would have returned by now
        indexer, new_index = self.index.get_loc_level(key)
        new_ser = self._constructor(self._values[indexer], index=new_index, copy=False)
        if isinstance(indexer, slice):
            new_ser._mgr.add_references(self._mgr)
        return new_ser.__finalize__(self)

    def _get_rows_with_mask(self, indexer: npt.NDArray[np.bool_]) -> Series:
        new_mgr = self._mgr.get_rows_with_mask(indexer)
        return self._constructor_from_mgr(new_mgr, axes=new_mgr.axes).__finalize__(self)

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

        # Similar to Index.get_value, but we do not fall back to positional
        loc = self.index.get_loc(label)

        if is_integer(loc):
            return self._values[loc]

        if isinstance(self.index, MultiIndex):
            mi = self.index
            new_values = self._values[loc]
            if len(new_values) == 1 and mi.nlevels == 1:
                # If more than one level left, we can not return a scalar
                return new_values[0]

            new_index = mi[loc]
            new_index = maybe_droplevels(new_index, label)
            new_ser = self._constructor(
                new_values, index=new_index, name=self.name, copy=False
            )
            if isinstance(loc, slice):
                new_ser._mgr.add_references(self._mgr)
            return new_ser.__finalize__(self)

        else:
            return self.iloc[loc]

    def __setitem__(self, key, value) -> None:
        if not CHAINED_WARNING_DISABLED:
            if sys.getrefcount(self) <= REF_COUNT and not com.is_local_in_caller_frame(
                self
            ):
                warnings.warn(
                    _chained_assignment_msg, ChainedAssignmentError, stacklevel=2
                )

        check_dict_or_set_indexers(key)
        key = com.apply_if_callable(key, self)

        if key is Ellipsis:
            key = slice(None)

        if isinstance(key, slice):
            indexer = self.index._convert_slice_indexer(key, kind="getitem")
            return self._set_values(indexer, value)

        try:
            self._set_with_engine(key, value)
        except KeyError:
            # We have a scalar (or for MultiIndex or object-dtype, scalar-like)
            #  key that is not present in self.index.
            # GH#12862 adding a new key to the Series
            self.loc[key] = value

        except (TypeError, ValueError, LossySetitemError):
            # The key was OK, but we cannot set the value losslessly
            indexer = self.index.get_loc(key)
            self._set_values(indexer, value)

        except InvalidIndexError as err:
            if isinstance(key, tuple) and not isinstance(self.index, MultiIndex):
                # cases with MultiIndex don't get here bc they raise KeyError
                # e.g. test_basic_getitem_setitem_corner
                raise KeyError(
                    "key of type tuple not found and not a MultiIndex"
                ) from err

            if com.is_bool_indexer(key):
                key = check_bool_indexer(self.index, key)
                key = np.asarray(key, dtype=bool)

                if (
                    is_list_like(value)
                    and len(value) != len(self)
                    and not isinstance(value, Series)
                    and not is_object_dtype(self.dtype)
                ):
                    # Series will be reindexed to have matching length inside
                    #  _where call below
                    # GH#44265
                    indexer = key.nonzero()[0]
                    self._set_values(indexer, value)
                    return

                # otherwise with listlike other we interpret series[mask] = other
                #  as series[mask] = other[mask]
                try:
                    self._where(~key, value, inplace=True)
                except InvalidIndexError:
                    # test_where_dups
                    self.iloc[key] = value
                return

            else:
                self._set_with(key, value)

    def _set_with_engine(self, key, value) -> None:
        loc = self.index.get_loc(key)

        # this is equivalent to self._values[key] = value
        self._mgr.setitem_inplace(loc, value)

    def _set_with(self, key, value) -> None:
        # We got here via exception-handling off of InvalidIndexError, so
        #  key should always be listlike at this point.
        assert not isinstance(key, tuple)

        if is_iterator(key):
            # Without this, the call to infer_dtype will consume the generator
            key = list(key)

        self._set_labels(key, value)

    def _set_labels(self, key, value) -> None:
        key = com.asarray_tuplesafe(key)
        indexer: np.ndarray = self.index.get_indexer(key)
        mask = indexer == -1
        if mask.any():
            raise KeyError(f"{key[mask]} not in index")
        self._set_values(indexer, value)

    def _set_values(self, key, value) -> None:
        if isinstance(key, (Index, Series)):
            key = key._values

        self._mgr = self._mgr.setitem(indexer=key, value=value)

    def _set_value(self, label, value, takeable: bool = False) -> None:
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
        """
        if not takeable:
            try:
                loc = self.index.get_loc(label)
            except KeyError:
                # set using a non-recursive method
                self.loc[label] = value
                return
        else:
            loc = label

        self._set_values(loc, value)

    # ----------------------------------------------------------------------
    # Unsorted

    def repeat(self, repeats: int | Sequence[int], axis: None = None) -> Series:
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
            Unused. Parameter needed for compatibility with DataFrame.

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
        >>> s = pd.Series(["a", "b", "c"])
        >>> s
        0    a
        1    b
        2    c
        dtype: str
        >>> s.repeat(2)
        0    a
        0    a
        1    b
        1    b
        2    c
        2    c
        dtype: str
        >>> s.repeat([1, 2, 3])
        0    a
        1    b
        1    b
        2    c
        2    c
        2    c
        dtype: str
        """
        nv.validate_repeat((), {"axis": axis})
        new_index = self.index.repeat(repeats)
        new_values = self._values.repeat(repeats)
        return self._constructor(new_values, index=new_index, copy=False).__finalize__(
            self, method="repeat"
        )

    @overload
    def reset_index(
        self,
        level: IndexLabel = ...,
        *,
        drop: Literal[False] = ...,
        name: Level = ...,
        inplace: Literal[False] = ...,
        allow_duplicates: bool = ...,
    ) -> DataFrame: ...

    @overload
    def reset_index(
        self,
        level: IndexLabel = ...,
        *,
        drop: Literal[True],
        name: Level = ...,
        inplace: Literal[False] = ...,
        allow_duplicates: bool = ...,
    ) -> Series: ...

    @overload
    def reset_index(
        self,
        level: IndexLabel = ...,
        *,
        drop: bool = ...,
        name: Level = ...,
        inplace: Literal[True],
        allow_duplicates: bool = ...,
    ) -> None: ...

    def reset_index(
        self,
        level: IndexLabel | None = None,
        *,
        drop: bool = False,
        name: Level = lib.no_default,
        inplace: bool = False,
        allow_duplicates: bool = False,
    ) -> DataFrame | Series | None:
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
        allow_duplicates : bool, default False
            Allow duplicate column labels to be created.

        Returns
        -------
        Series or DataFrame or None
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
        >>> s = pd.Series(
        ...     [1, 2, 3, 4],
        ...     name="foo",
        ...     index=pd.Index(["a", "b", "c", "d"], name="idx"),
        ... )

        Generate a DataFrame with default index.

        >>> s.reset_index()
          idx  foo
        0   a    1
        1   b    2
        2   c    3
        3   d    4

        To specify the name of the new column use `name`.

        >>> s.reset_index(name="values")
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

        The `level` parameter is interesting for Series with a multi-level
        index.

        >>> arrays = [
        ...     np.array(["bar", "bar", "baz", "baz"]),
        ...     np.array(["one", "two", "one", "two"]),
        ... ]
        >>> s2 = pd.Series(
        ...     range(4),
        ...     name="foo",
        ...     index=pd.MultiIndex.from_arrays(arrays, names=["a", "b"]),
        ... )

        To remove a specific level from the Index, use `level`.

        >>> s2.reset_index(level="a")
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
            new_index = default_index(len(self))
            if level is not None:
                level_list: Sequence[Hashable]
                if not isinstance(level, (tuple, list)):
                    level_list = [level]
                else:
                    level_list = level
                level_list = [self.index._get_level_number(lev) for lev in level_list]
                if len(level_list) < self.index.nlevels:
                    new_index = self.index.droplevel(level_list)

            if inplace:
                self.index = new_index
            else:
                new_ser = self.copy(deep=False)
                new_ser.index = new_index
                return new_ser.__finalize__(self, method="reset_index")
        elif inplace:
            raise TypeError(
                "Cannot reset_index inplace on a Series to create a DataFrame"
            )
        else:
            if name is lib.no_default:
                # For backwards compatibility, keep columns as [0] instead of
                #  [None] when self.name is None
                if self.name is None:
                    name = 0
                else:
                    name = self.name

            df = self.to_frame(name)
            return df.reset_index(
                level=level, drop=drop, allow_duplicates=allow_duplicates
            )
        return None

    # ----------------------------------------------------------------------
    # Rendering Methods

    def __repr__(self) -> str:
        """
        Return a string representation for a particular Series.
        """
        repr_params = fmt.get_series_repr_params()
        return self.to_string(**repr_params)

    @overload
    def to_string(
        self,
        buf: None = ...,
        *,
        na_rep: str = ...,
        float_format: str | None = ...,
        header: bool = ...,
        index: bool = ...,
        length: bool = ...,
        dtype=...,
        name=...,
        max_rows: int | None = ...,
        min_rows: int | None = ...,
    ) -> str: ...

    @overload
    def to_string(
        self,
        buf: FilePath | WriteBuffer[str],
        *,
        na_rep: str = ...,
        float_format: str | None = ...,
        header: bool = ...,
        index: bool = ...,
        length: bool = ...,
        dtype=...,
        name=...,
        max_rows: int | None = ...,
        min_rows: int | None = ...,
    ) -> None: ...

    @deprecate_nonkeyword_arguments(
        Pandas4Warning, allowed_args=["self", "buf"], name="to_string"
    )
    def to_string(
        self,
        buf: FilePath | WriteBuffer[str] | None = None,
        na_rep: str = "NaN",
        float_format: str | None = None,
        header: bool = True,
        index: bool = True,
        length: bool = False,
        dtype: bool = False,
        name: bool = False,
        max_rows: int | None = None,
        min_rows: int | None = None,
    ) -> str | None:
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

        See Also
        --------
        Series.to_dict : Convert Series to dict object.
        Series.to_frame : Convert Series to DataFrame object.
        Series.to_markdown : Print Series in Markdown-friendly format.
        Series.to_timestamp : Cast to DatetimeIndex of Timestamps.

        Examples
        --------
        >>> ser = pd.Series([1, 2, 3]).to_string()
        >>> ser
        '0    1\\n1    2\\n2    3'
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
                f"of result is {type(result).__name__!r}"
            )

        if buf is None:
            return result
        elif hasattr(buf, "write"):
            buf.write(result)
        else:
            with open(buf, "w", encoding="utf-8") as f:
                f.write(result)
        return None

    @overload
    def to_markdown(
        self,
        buf: None = ...,
        *,
        mode: str = ...,
        index: bool = ...,
        storage_options: StorageOptions | None = ...,
        **kwargs,
    ) -> str: ...

    @overload
    def to_markdown(
        self,
        buf: IO[str],
        *,
        mode: str = ...,
        index: bool = ...,
        storage_options: StorageOptions | None = ...,
        **kwargs,
    ) -> None: ...

    @overload
    def to_markdown(
        self,
        buf: IO[str] | None,
        *,
        mode: str = ...,
        index: bool = ...,
        storage_options: StorageOptions | None = ...,
        **kwargs,
    ) -> str | None: ...

    @deprecate_nonkeyword_arguments(
        Pandas4Warning, allowed_args=["self", "buf"], name="to_markdown"
    )
    def to_markdown(
        self,
        buf: IO[str] | None = None,
        mode: str = "wt",
        index: bool = True,
        storage_options: StorageOptions | None = None,
        **kwargs,
    ) -> str | None:
        """
        Print Series in Markdown-friendly format.

        Parameters
        ----------
        buf : str, Path or StringIO-like, optional, default None
            Buffer to write to. If None, the output is returned as a string.
        mode : str, optional
            Mode in which file is opened, "wt" by default.
        index : bool, optional, default True
            Add index (row) labels.

        storage_options : dict, optional
            Extra options that make sense for a particular storage connection, e.g.
            host, port, username, password, etc. For HTTP(S) URLs the key-value pairs
            are forwarded to ``urllib.request.Request`` as header options. For other
            URLs (e.g. starting with "s3://", and "gcs://") the key-value pairs are
            forwarded to ``fsspec.open``. Please see ``fsspec`` and ``urllib`` for more
            details, and for more examples on storage options refer `here
            <https://pandas.pydata.org/docs/user_guide/io.html?
            highlight=storage_options#reading-writing-remote-files>`_.

        **kwargs
            These parameters will be passed to `tabulate \
                <https://pypi.org/project/tabulate>`_.

        Returns
        -------
        str
            Series in Markdown-friendly format.

        See Also
        --------
        Series.to_frame : Rrite a text representation of object to the system clipboard.
        Series.to_latex : Render Series to LaTeX-formatted table.

        Notes
        -----
        Requires the `tabulate <https://pypi.org/project/tabulate>`_ package.

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

            Output markdown with a tabulate option.

            >>> print(s.to_markdown(tablefmt="grid"))
            +----+----------+
            |    | animal   |
            +====+==========+
            |  0 | elk      |
            +----+----------+
            |  1 | pig      |
            +----+----------+
            |  2 | dog      |
            +----+----------+
            |  3 | quetzal  |
            +----+----------+
        """
        return self.to_frame().to_markdown(
            buf, mode=mode, index=index, storage_options=storage_options, **kwargs
        )

    # ----------------------------------------------------------------------

    def items(self) -> Iterable[tuple[Hashable, Any]]:
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
        >>> s = pd.Series(["A", "B", "C"])
        >>> for index, value in s.items():
        ...     print(f"Index : {index}, Value : {value}")
        Index : 0, Value : A
        Index : 1, Value : B
        Index : 2, Value : C
        """
        return zip(iter(self.index), iter(self), strict=True)

    # ----------------------------------------------------------------------
    # Misc public methods

    def keys(self) -> Index:
        """
        Return alias for index.

        Returns
        -------
        Index
            Index of the Series.

        See Also
        --------
        Series.index : The index (axis labels) of the Series.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3], index=[0, 1, 2])
        >>> s.keys()
        Index([0, 1, 2], dtype='int64')
        """
        return self.index

    @overload
    def to_dict(
        self, *, into: type[MutableMappingT] | MutableMappingT
    ) -> MutableMappingT: ...

    @overload
    def to_dict(self, *, into: type[dict] = ...) -> dict: ...

    # error: Incompatible default for argument "into" (default has type "type[
    # dict[Any, Any]]", argument has type "type[MutableMappingT] | MutableMappingT")
    def to_dict(
        self,
        *,
        into: type[MutableMappingT] | MutableMappingT = dict,  # type: ignore[assignment]
    ) -> MutableMappingT:
        """
        Convert Series to {label -> value} dict or dict-like object.

        Parameters
        ----------
        into : class, default dict
            The collections.abc.MutableMapping subclass to use as the return
            object. Can be the actual class or an empty instance of the mapping
            type you want.  If you want a collections.defaultdict, you must
            pass it initialized.

        Returns
        -------
        collections.abc.MutableMapping
            Key-value representation of Series.

        See Also
        --------
        Series.to_list: Converts Series to a list of the values.
        Series.to_numpy: Converts Series to NumPy ndarray.
        Series.array: ExtensionArray of the data backing this Series.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4])
        >>> s.to_dict()
        {0: 1, 1: 2, 2: 3, 3: 4}
        >>> from collections import OrderedDict, defaultdict
        >>> s.to_dict(into=OrderedDict)
        OrderedDict([(0, 1), (1, 2), (2, 3), (3, 4)])
        >>> dd = defaultdict(list)
        >>> s.to_dict(into=dd)
        defaultdict(<class 'list'>, {0: 1, 1: 2, 2: 3, 3: 4})
        """
        # GH16122
        into_c = com.standardize_mapping(into)

        if is_object_dtype(self.dtype) or isinstance(self.dtype, ExtensionDtype):
            return into_c((k, maybe_box_native(v)) for k, v in self.items())
        else:
            # Not an object dtype => all types will be the same so let the default
            # indexer return native python type
            return into_c(self.items())

    def to_frame(self, name: Hashable = lib.no_default) -> DataFrame:
        """
        Convert Series to DataFrame.

        Parameters
        ----------
        name : object, optional
            The passed name should substitute for the series name (if it has
            one).

        Returns
        -------
        DataFrame
            DataFrame representation of Series.

        See Also
        --------
        Series.to_dict : Convert Series to dict object.

        Examples
        --------
        >>> s = pd.Series(["a", "b", "c"], name="vals")
        >>> s.to_frame()
          vals
        0    a
        1    b
        2    c
        """
        columns: Index
        if name is lib.no_default:
            name = self.name
            if name is None:
                # default to [0], same as we would get with DataFrame(self)
                columns = default_index(1)
            else:
                columns = Index([name])
        else:
            columns = Index([name])

        mgr = self._mgr.to_2d_mgr(columns)
        df = self._constructor_expanddim_from_mgr(mgr, axes=mgr.axes)
        return df.__finalize__(self, method="to_frame")

    @classmethod
    def from_arrow(cls, data: ArrowArrayExportable | ArrowStreamExportable) -> Series:
        """
        Construct a Series from an array-like Arrow object.

        This function accepts any Arrow-compatible array-like object implementing
        the `Arrow PyCapsule Protocol`_ (i.e. having an ``__arrow_c_array__``
        or ``__arrow_c_stream__`` method).

        This function currently relies on ``pyarrow`` to convert the object
        in Arrow format to pandas.

        .. _Arrow PyCapsule Protocol: https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html

        .. versionadded:: 3.0

        Parameters
        ----------
        data : pyarrow.Array or Arrow-compatible object
            Any array-like object implementing the Arrow PyCapsule Protocol
            (i.e. has an ``__arrow_c_array__`` or ``__arrow_c_stream__``
            method).

        Returns
        -------
        Series

        See Also
        --------
        DataFrame.from_arrow : Construct a DataFrame from an Arrow object.

        Examples
        --------
        >>> import pyarrow as pa
        >>> arrow_array = pa.array([1, 2, 3])
        >>> pd.Series.from_arrow(arrow_array)
        0    1
        1    2
        2    3
        dtype: int64
        """
        pa = import_optional_dependency("pyarrow", min_version="14.0.0")
        if not isinstance(data, (pa.Array, pa.ChunkedArray)):
            if not (
                hasattr(data, "__arrow_c_array__")
                or hasattr(data, "__arrow_c_stream__")
            ):
                # explicitly test this, because otherwise we would accept variour other
                # input types through the pa.chunked_array(..) call
                raise TypeError(
                    "Expected an Arrow-compatible array-like object (i.e. having an "
                    "'_arrow_c_array__' or '__arrow_c_stream__' method), got "
                    f"'{type(data).__name__}' instead."
                )
            # using chunked_array() as it works for both arrays and streams
            pa_array = pa.chunked_array(data)
        else:
            pa_array = data

        ser = pa_array.to_pandas()
        return ser

    def _set_name(self, name, inplace: bool = False) -> Series:
        """
        Set the Series name.

        Parameters
        ----------
        name : str
        inplace : bool
            Whether to modify `self` directly or return a copy.
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        ser = self if inplace else self.copy(deep=False)
        ser.name = name
        return ser

    @Appender(
        dedent(
            """
        Examples
        --------
        >>> ser = pd.Series([390., 350., 30., 20.],
        ...                 index=['Falcon', 'Falcon', 'Parrot', 'Parrot'],
        ...                 name="Max Speed")
        >>> ser
        Falcon    390.0
        Falcon    350.0
        Parrot     30.0
        Parrot     20.0
        Name: Max Speed, dtype: float64

        We can pass a list of values to group the Series data by custom labels:

        >>> ser.groupby(["a", "b", "a", "b"]).mean()
        a    210.0
        b    185.0
        Name: Max Speed, dtype: float64

        Grouping by numeric labels yields similar results:

        >>> ser.groupby([0, 1, 0, 1]).mean()
        0    210.0
        1    185.0
        Name: Max Speed, dtype: float64

        We can group by a level of the index:

        >>> ser.groupby(level=0).mean()
        Falcon    370.0
        Parrot     25.0
        Name: Max Speed, dtype: float64

        We can group by a condition applied to the Series values:

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

        We can also group by the 'Type' level of the hierarchical index
        to get the mean speed for each type:

        >>> ser.groupby(level="Type").mean()
        Type
        Captive    210.0
        Wild       185.0
        Name: Max Speed, dtype: float64

        We can also choose to include `NA` in group keys or not by defining
        `dropna` parameter, the default setting is `True`.

        >>> ser = pd.Series([1, 2, 3, 3], index=["a", 'a', 'b', np.nan])
        >>> ser.groupby(level=0).sum()
        a    3
        b    3
        dtype: int64

        To include `NA` values in the group keys, set `dropna=False`:

        >>> ser.groupby(level=0, dropna=False).sum()
        a    3
        b    3
        NaN  3
        dtype: int64

        We can also group by a custom list with NaN values to handle
        missing group labels:

        >>> arrays = ['Falcon', 'Falcon', 'Parrot', 'Parrot']
        >>> ser = pd.Series([390., 350., 30., 20.], index=arrays, name="Max Speed")
        >>> ser.groupby(["a", "b", "a", np.nan]).mean()
        a    210.0
        b    350.0
        Name: Max Speed, dtype: float64

        >>> ser.groupby(["a", "b", "a", np.nan], dropna=False).mean()
        a    210.0
        b    350.0
        NaN   20.0
        Name: Max Speed, dtype: float64
        """
        )
    )
    @Appender(_shared_docs["groupby"] % _shared_doc_kwargs)
    @deprecate_nonkeyword_arguments(
        Pandas4Warning, allowed_args=["self", "by", "level"], name="groupby"
    )
    def groupby(
        self,
        by=None,
        level: IndexLabel | None = None,
        as_index: bool = True,
        sort: bool = True,
        group_keys: bool = True,
        observed: bool = True,
        dropna: bool = True,
    ) -> SeriesGroupBy:
        from pandas.core.groupby.generic import SeriesGroupBy

        if level is None and by is None:
            raise TypeError("You have to supply one of 'by' and 'level'")
        if not as_index:
            raise TypeError("as_index=False only valid with DataFrame")

        return SeriesGroupBy(
            obj=self,
            keys=by,
            level=level,
            as_index=as_index,
            sort=sort,
            group_keys=group_keys,
            observed=observed,
            dropna=dropna,
        )

    # ----------------------------------------------------------------------
    # Statistics, overridden ndarray methods

    # TODO: integrate bottleneck
    def count(self) -> int:
        """
        Return number of non-NA/null observations in the Series.

        Returns
        -------
        int
            Number of non-null values in the Series.

        See Also
        --------
        DataFrame.count : Count non-NA cells for each column or row.

        Examples
        --------
        >>> s = pd.Series([0.0, 1.0, np.nan])
        >>> s.count()
        2
        """
        return maybe_unbox_numpy_scalar(notna(self._values).sum().astype("int64"))

    def mode(self, dropna: bool = True) -> Series:
        """
        Return the mode(s) of the Series.

        The mode is the value that appears most often. There can be multiple modes.

        Always returns Series even if only one value is returned.

        Parameters
        ----------
        dropna : bool, default True
            Don't consider counts of NaN/NaT.

        Returns
        -------
        Series
            Modes of the Series in sorted order.

        See Also
        --------
        numpy.mode : Equivalent numpy function for computing median.
        Series.sum : Sum of the values.
        Series.median : Median of the values.
        Series.std : Standard deviation of the values.
        Series.var : Variance of the values.
        Series.min : Minimum value.
        Series.max : Maximum value.

        Examples
        --------
        >>> s = pd.Series([2, 4, 2, 2, 4, None])
        >>> s.mode()
        0    2.0
        dtype: float64

        More than one mode:

        >>> s = pd.Series([2, 4, 8, 2, 4, None])
        >>> s.mode()
        0    2.0
        1    4.0
        dtype: float64

        With and without considering null value:

        >>> s = pd.Series([2, 4, None, None, 4, None])
        >>> s.mode(dropna=False)
        0   NaN
        dtype: float64
        >>> s = pd.Series([2, 4, None, None, 4, None])
        >>> s.mode()
        0    4.0
        dtype: float64
        """
        # TODO: Add option for bins like value_counts()
        values = self._values
        if isinstance(values, np.ndarray):
            res_values, _ = algorithms.mode(values, dropna=dropna)
        else:
            res_values = values._mode(dropna=dropna)

        # Ensure index is type stable (should always use int index)
        return self._constructor(
            res_values,
            index=range(len(res_values)),
            name=self.name,
            copy=False,
            dtype=self.dtype,
        ).__finalize__(self, method="mode")

    def unique(self) -> ArrayLike:
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
        Series.drop_duplicates : Return Series with duplicate values removed.
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
            * Datetime without Timezone
            * Timedelta
            * Interval
            * Sparse
            * IntegerNA

        See Examples section.

        Examples
        --------
        >>> pd.Series([2, 1, 3, 3], name="A").unique()
        array([2, 1, 3])

        >>> pd.Series([pd.Timestamp("2016-01-01") for _ in range(3)]).unique()
        <DatetimeArray>
        ['2016-01-01 00:00:00']
        Length: 1, dtype: datetime64[us]

        >>> pd.Series(
        ...     [pd.Timestamp("2016-01-01", tz="US/Eastern") for _ in range(3)]
        ... ).unique()
        <DatetimeArray>
        ['2016-01-01 00:00:00-05:00']
        Length: 1, dtype: datetime64[us, US/Eastern]

        A Categorical will return categories in the order of
        appearance and with the same dtype.

        >>> pd.Series(pd.Categorical(list("baabc"))).unique()
        ['b', 'a', 'c']
        Categories (3, str): ['a', 'b', 'c']
        >>> pd.Series(
        ...     pd.Categorical(list("baabc"), categories=list("abc"), ordered=True)
        ... ).unique()
        ['b', 'a', 'c']
        Categories (3, str): ['a' < 'b' < 'c']
        """
        return super().unique()

    @overload
    def drop_duplicates(
        self,
        *,
        keep: DropKeep = ...,
        inplace: Literal[False] = ...,
        ignore_index: bool = ...,
    ) -> Series: ...

    @overload
    def drop_duplicates(
        self, *, keep: DropKeep = ..., inplace: Literal[True], ignore_index: bool = ...
    ) -> None: ...

    @overload
    def drop_duplicates(
        self, *, keep: DropKeep = ..., inplace: bool = ..., ignore_index: bool = ...
    ) -> Series | None: ...

    def drop_duplicates(
        self,
        *,
        keep: DropKeep = "first",
        inplace: bool = False,
        ignore_index: bool = False,
    ) -> Series | None:
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

        ignore_index : bool, default ``False``
            If ``True``, the resulting axis will be labeled 0, 1, …, n - 1.

            .. versionadded:: 2.0.0

        Returns
        -------
        Series or None
            Series with duplicates dropped or None if ``inplace=True``.

        See Also
        --------
        Index.drop_duplicates : Equivalent method on Index.
        DataFrame.drop_duplicates : Equivalent method on DataFrame.
        Series.duplicated : Related method on Series, indicating duplicate
            Series values.
        Series.unique : Return unique values as an array.

        Examples
        --------
        Generate a Series with duplicated entries.

        >>> s = pd.Series(
        ...     ["llama", "cow", "llama", "beetle", "llama", "hippo"], name="animal"
        ... )
        >>> s
        0     llama
        1       cow
        2     llama
        3    beetle
        4     llama
        5     hippo
        Name: animal, dtype: str

        With the 'keep' parameter, the selection behavior of duplicated values
        can be changed. The value 'first' keeps the first occurrence for each
        set of duplicated entries. The default value of keep is 'first'.

        >>> s.drop_duplicates()
        0     llama
        1       cow
        3    beetle
        5     hippo
        Name: animal, dtype: str

        The value 'last' for parameter 'keep' keeps the last occurrence for
        each set of duplicated entries.

        >>> s.drop_duplicates(keep="last")
        1       cow
        3    beetle
        4     llama
        5     hippo
        Name: animal, dtype: str

        The value ``False`` for parameter 'keep' discards all sets of
        duplicated entries.

        >>> s.drop_duplicates(keep=False)
        1       cow
        3    beetle
        5     hippo
        Name: animal, dtype: str
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        result = super().drop_duplicates(keep=keep)

        if ignore_index:
            result.index = default_index(len(result))

        if inplace:
            self._update_inplace(result)
            return None
        else:
            return result

    def duplicated(self, keep: DropKeep = "first") -> Series:
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
        Series[bool]
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

        >>> animals = pd.Series(["llama", "cow", "llama", "beetle", "llama"])
        >>> animals.duplicated()
        0    False
        1    False
        2     True
        3    False
        4     True
        dtype: bool

        which is equivalent to

        >>> animals.duplicated(keep="first")
        0    False
        1    False
        2     True
        3    False
        4     True
        dtype: bool

        By using 'last', the last occurrence of each set of duplicated values
        is set on False and all others on True:

        >>> animals.duplicated(keep="last")
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
        res = self._duplicated(keep=keep)
        result = self._constructor(res, index=self.index, copy=False)
        return result.__finalize__(self, method="duplicated")

    def idxmin(self, axis: Axis = 0, skipna: bool = True, *args, **kwargs) -> Hashable:
        """
        Return the row label of the minimum value.

        If multiple values equal the minimum, the first row label with that
        value is returned.

        Parameters
        ----------
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.
        skipna : bool, default True
            Exclude NA/null values. If the entire Series is NA, or if ``skipna=False``
            and there is an NA value, this method will raise a ``ValueError``.
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
        >>> s = pd.Series(data=[1, None, 4, 1], index=["A", "B", "C", "D"])
        >>> s
        A    1.0
        B    NaN
        C    4.0
        D    1.0
        dtype: float64

        >>> s.idxmin()
        'A'
        """
        axis = self._get_axis_number(axis)
        iloc = self.argmin(axis, skipna, *args, **kwargs)
        return self.index[iloc]

    def idxmax(self, axis: Axis = 0, skipna: bool = True, *args, **kwargs) -> Hashable:
        """
        Return the row label of the maximum value.

        If multiple values equal the maximum, the first row label with that
        value is returned.

        Parameters
        ----------
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.
        skipna : bool, default True
            Exclude NA/null values. If the entire Series is NA, or if ``skipna=False``
            and there is an NA value, this method will raise a ``ValueError``.
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
        >>> s = pd.Series(data=[1, None, 4, 3, 4], index=["A", "B", "C", "D", "E"])
        >>> s
        A    1.0
        B    NaN
        C    4.0
        D    3.0
        E    4.0
        dtype: float64

        >>> s.idxmax()
        'C'
        """
        axis = self._get_axis_number(axis)
        iloc = self.argmax(axis, skipna, *args, **kwargs)
        return self.index[iloc]

    def round(self, decimals: int = 0, *args, **kwargs) -> Series:
        """
        Round each value in a Series to the given number of decimals.

        Parameters
        ----------
        decimals : int, default 0
            Number of decimal places to round to. If decimals is negative,
            it specifies the number of positions to the left of the decimal point.
        *args, **kwargs
            Additional arguments and keywords have no effect but might be
            accepted for compatibility with NumPy.

        Returns
        -------
        Series
            Rounded values of the Series.

        See Also
        --------
        numpy.around : Round values of an np.array.
        DataFrame.round : Round values of a DataFrame.
        Series.dt.round : Round values of data to the specified freq.

        Notes
        -----
        For values exactly halfway between rounded decimal values, pandas rounds
        to the nearest even value (e.g. -0.5 and 0.5 round to 0.0, 1.5 and 2.5
        round to 2.0, etc.).

        Examples
        --------
        >>> s = pd.Series([-0.5, 0.1, 2.5, 1.3, 2.7])
        >>> s.round()
        0   -0.0
        1    0.0
        2    2.0
        3    1.0
        4    3.0
        dtype: float64
        """

        nv.validate_round(args, kwargs)

        if len(self) == 0:
            return self.copy()

        if is_object_dtype(self.dtype):
            values = self._values
            result = lib.map_infer(values, lambda x: round(x, decimals), convert=False)
            return self._constructor(result, index=self.index, copy=False).__finalize__(
                self, method="round"
            )
        new_mgr = self._mgr.round(decimals=decimals)
        return self._constructor_from_mgr(new_mgr, axes=new_mgr.axes).__finalize__(
            self, method="round"
        )

    @overload
    def quantile(
        self, q: float = ..., interpolation: QuantileInterpolation = ...
    ) -> float: ...

    @overload
    def quantile(
        self,
        q: Sequence[float] | AnyArrayLike,
        interpolation: QuantileInterpolation = ...,
    ) -> Series: ...

    @overload
    def quantile(
        self,
        q: float | Sequence[float] | AnyArrayLike = ...,
        interpolation: QuantileInterpolation = ...,
    ) -> float | Series: ...

    def quantile(
        self,
        q: float | Sequence[float] | AnyArrayLike = 0.5,
        interpolation: QuantileInterpolation = "linear",
    ) -> float | Series:
        """
        Return value at the given quantile.

        Parameters
        ----------
        q : float or array-like, default 0.5 (50% quantile)
            The quantile(s) to compute, which can lie in range: 0 <= q <= 1.
        interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
            This optional parameter specifies the interpolation method to use,
            when the desired quantile lies between two data points `i` and `j`:

                * linear: `i + (j - i) * (x-i)/(j-i)`, where `(x-i)/(j-i)` is
                  the fractional part of the index surrounded by `i > j`.
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
        core.window.Rolling.quantile : Calculate the rolling quantile.
        numpy.percentile : Returns the q-th percentile(s) of the array elements.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4])
        >>> s.quantile(0.5)
        2.5
        >>> s.quantile([0.25, 0.5, 0.75])
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
            idx = Index(q, dtype=np.float64)
            return self._constructor(result, index=idx, name=self.name)
        else:
            # scalar
            return maybe_unbox_numpy_scalar(result.iloc[0])

    def corr(
        self,
        other: Series,
        method: CorrelationMethod = "pearson",
        min_periods: int | None = None,
    ) -> float:
        """
        Compute correlation with `other` Series, excluding missing values.

        The two `Series` objects are not required to be the same length and will be
        aligned internally before the correlation function is applied.

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

            .. warning::
                Note that the returned matrix from corr will have 1 along the
                diagonals and will be symmetric regardless of the callable's
                behavior.
        min_periods : int, optional
            Minimum number of observations needed to have a valid result.

        Returns
        -------
        float
            Correlation with other.

        See Also
        --------
        DataFrame.corr : Compute pairwise correlation between columns.
        DataFrame.corrwith : Compute pairwise correlation with another
            DataFrame or Series.

        Notes
        -----
        Pearson, Kendall and Spearman correlation are currently computed using pairwise complete observations.

        * `Pearson correlation coefficient <https://en.wikipedia.org/wiki/Pearson_correlation_coefficient>`_
        * `Kendall rank correlation coefficient <https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient>`_
        * `Spearman's rank correlation coefficient <https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient>`_

        Automatic data alignment: as with all pandas operations, automatic data alignment is performed for this method.
        ``corr()`` automatically considers values with matching indices.

        Examples
        --------
        >>> def histogram_intersection(a, b):
        ...     v = np.minimum(a, b).sum().round(decimals=1)
        ...     return v
        >>> s1 = pd.Series([0.2, 0.0, 0.6, 0.2])
        >>> s2 = pd.Series([0.3, 0.6, 0.0, 0.1])
        >>> s1.corr(s2, method=histogram_intersection)
        0.3

        Pandas auto-aligns the values with matching indices

        >>> s1 = pd.Series([1, 2, 3], index=[0, 1, 2])
        >>> s2 = pd.Series([1, 2, 3], index=[2, 1, 0])
        >>> s1.corr(s2)
        -1.0

        If the input is a constant array, the correlation is not defined in this case,
        and ``np.nan`` is returned.

        >>> s1 = pd.Series([0.45, 0.45])
        >>> s1.corr(s1)
        nan
        """  # noqa: E501
        this, other = self.align(other, join="inner")
        if len(this) == 0:
            return np.nan

        this_values = this.to_numpy(dtype=float, na_value=np.nan, copy=False)
        other_values = other.to_numpy(dtype=float, na_value=np.nan, copy=False)

        if method in ["pearson", "spearman", "kendall"] or callable(method):
            result = nanops.nancorr(
                this_values, other_values, method=method, min_periods=min_periods
            )
            result = maybe_unbox_numpy_scalar(result)
            return result

        raise ValueError(
            "method must be either 'pearson', "
            "'spearman', 'kendall', or a callable, "
            f"'{method}' was supplied"
        )

    def cov(
        self,
        other: Series,
        min_periods: int | None = None,
        ddof: int | None = 1,
    ) -> float:
        """
        Compute covariance with Series, excluding missing values.

        The two `Series` objects are not required to be the same length and
        will be aligned internally before the covariance is calculated.

        Parameters
        ----------
        other : Series
            Series with which to compute the covariance.
        min_periods : int, optional
            Minimum number of observations needed to have a valid result.
        ddof : int, default 1
            Delta degrees of freedom.  The divisor used in calculations
            is ``N - ddof``, where ``N`` represents the number of elements.

        Returns
        -------
        float
            Covariance between Series and other normalized by N-1
            (unbiased estimator).

        See Also
        --------
        DataFrame.cov : Compute pairwise covariance of columns.

        Examples
        --------
        >>> s1 = pd.Series([0.90010907, 0.13484424, 0.62036035])
        >>> s2 = pd.Series([0.12528585, 0.26962463, 0.51111198])
        >>> s1.cov(s2)
        -0.01685762652715874
        """
        this, other = self.align(other, join="inner")
        if len(this) == 0:
            return np.nan
        this_values = this.to_numpy(dtype=float, na_value=np.nan, copy=False)
        other_values = other.to_numpy(dtype=float, na_value=np.nan, copy=False)
        result = nanops.nancov(
            this_values, other_values, min_periods=min_periods, ddof=ddof
        )
        result = maybe_unbox_numpy_scalar(result)
        return result

    def diff(self, periods: int = 1) -> Series:
        """
        First discrete difference of Series elements.

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
        The result is calculated according to current dtype in Series,
        however dtype of the result is always float64.

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

        Overflow in input dtype

        >>> s = pd.Series([1, 0], dtype=np.uint8)
        >>> s.diff()
        0      NaN
        1    255.0
        dtype: float64
        """
        if not lib.is_integer(periods):
            if not (is_float(periods) and periods.is_integer()):
                raise ValueError("periods must be an integer")
        result = algorithms.diff(self._values, periods)
        return self._constructor(
            result, index=self.index.view(), copy=False
        ).__finalize__(self, method="diff")

    def autocorr(self, lag: int = 1) -> float:
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
        return self.corr(cast(Series, self.shift(lag)))

    def dot(self, other: AnyArrayLike | DataFrame) -> Series | np.ndarray:
        """
        Compute the dot product between the Series and the columns of other.

        This method computes the dot product between the Series and another
        one, or the Series and each columns of a DataFrame, or the Series and
        each columns of an array.

        It can also be called using `self @ other`.

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

            left = self.reindex(index=common)
            right = other.reindex(index=common)
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
            common_type = find_common_type([self.dtypes, *list(other.dtypes)])
            return self._constructor(
                np.dot(lvals, rvals), index=other.columns, copy=False, dtype=common_type
            ).__finalize__(self, method="dot")
        elif isinstance(other, Series):
            result = np.dot(lvals, rvals)
        elif isinstance(rvals, np.ndarray):
            result = np.dot(lvals, rvals)
        else:  # pragma: no cover
            raise TypeError(f"unsupported type: {type(other)}")
        return maybe_unbox_numpy_scalar(result)

    def __matmul__(self, other):
        """
        Matrix multiplication using binary `@` operator.
        """
        return self.dot(other)

    def __rmatmul__(self, other):
        """
        Matrix multiplication using binary `@` operator.
        """
        return self.dot(np.transpose(other))

    # Signature of "searchsorted" incompatible with supertype "IndexOpsMixin"
    def searchsorted(  # type: ignore[override]
        self,
        value: NumpyValueArrayLike | ExtensionArray,
        side: Literal["left", "right"] = "left",
        sorter: NumpySorter | None = None,
    ) -> npt.NDArray[np.intp] | np.intp:
        """
        Find indices where elements should be inserted to maintain order.

        Find the indices into a sorted Series `self` such that, if the
        corresponding elements in `value` were inserted before the indices,
        the order of `self` would be preserved.

        .. note::
            The Series *must* be monotonically sorted, otherwise
            wrong locations will likely be returned. Pandas does *not*
            check this for you.

        Parameters
        ----------
        value : array-like or scalar
            Values to insert into `self`.
        side : {'left', 'right'}, optional
            If 'left', the index of the first suitable location found is given.
            If 'right', return the last such index.  If there is no suitable
            index, return either 0 or N (where N is the length of `self`).
        sorter : 1-D array-like, optional
            Optional array of integer indices that sort `self` into ascending
            order. They are typically the result of ``np.argsort``.

        Returns
        -------
        int or array of int
            A scalar or array of insertion points with the
            same shape as `value`.

        See Also
        --------
        sort_values : Sort by the values along either axis.
        numpy.searchsorted : Similar method from NumPy.

        Notes
        -----
        Binary search is used to find the required insertion points.

        Examples
        --------
        >>> ser = pd.Series([1, 2, 3])
        >>> ser
        0    1
        1    2
        2    3
        dtype: int64
        >>> ser.searchsorted(4)
        np.int64(3)
        >>> ser.searchsorted([0, 4])
        array([0, 3])
        >>> ser.searchsorted([1, 3], side="left")
        array([0, 2])
        >>> ser.searchsorted([1, 3], side="right")
        array([1, 3])
        >>> ser = pd.Series(pd.to_datetime(["3/11/2000", "3/12/2000", "3/13/2000"]))
        >>> ser
        0   2000-03-11
        1   2000-03-12
        2   2000-03-13
        dtype: datetime64[us]
        >>> ser.searchsorted("3/14/2000")
        np.int64(3)
        >>> ser = pd.Categorical(
        ...     ["apple", "bread", "bread", "cheese", "milk"], ordered=True
        ... )
        >>> ser
        ['apple', 'bread', 'bread', 'cheese', 'milk']
        Categories (4, str): ['apple' < 'bread' < 'cheese' < 'milk']
        >>> ser.searchsorted("bread")
        np.int64(1)
        >>> ser.searchsorted(["bread"], side="right")
        array([3])

        If the values are not monotonically sorted, wrong locations
        may be returned:

        >>> ser = pd.Series([2, 1, 3])
        >>> ser
        0    2
        1    1
        2    3
        dtype: int64
        >>> ser.searchsorted(1)  # doctest: +SKIP
        0  # wrong result, correct would be 1
        """
        return base.IndexOpsMixin.searchsorted(self, value, side=side, sorter=sorter)

    # -------------------------------------------------------------------
    # Combination

    def _append_internal(self, to_append: Series, ignore_index: bool = False) -> Series:
        from pandas.core.reshape.concat import concat

        return concat([self, to_append], ignore_index=ignore_index)

    def compare(
        self,
        other: Series,
        align_axis: Axis = 1,
        keep_shape: bool = False,
        keep_equal: bool = False,
        result_names: Suffixes = ("self", "other"),
    ) -> DataFrame | Series:
        """
        Compare to another Series and show the differences.

        Parameters
        ----------
        other : Series
            Object to compare with.

        align_axis : {{0 or 'index', 1 or 'columns'}}, default 1
            Determine which axis to align the comparison on.

            * 0, or 'index' : Resulting differences are stacked vertically
              with rows drawn alternately from self and other.
            * 1, or 'columns' : Resulting differences are aligned horizontally
              with columns drawn alternately from self and other.

        keep_shape : bool, default False
            If true, all rows and columns are kept.
            Otherwise, only the ones with different values are kept.

        keep_equal : bool, default False
            If true, the result keeps values that are equal.
            Otherwise, equal values are shown as NaNs.

        result_names : tuple, default ('self', 'other')
            Set the dataframes names in the comparison.

        Returns
        -------
        Series or DataFrame
            If axis is 0 or 'index' the result will be a Series.
            The resulting index will be a MultiIndex with 'self' and 'other'
            stacked alternately at the inner level.

            If axis is 1 or 'columns' the result will be a DataFrame.
            It will have two columns namely 'self' and 'other'.

        See Also
        --------
        DataFrame.compare : Compare with another DataFrame and show differences.

        Notes
        -----
        Matching NaNs will not appear as a difference.

        Examples
        --------
        >>> s1 = pd.Series(["a", "b", "c", "d", "e"])
        >>> s2 = pd.Series(["a", "a", "c", "b", "e"])

        Align the differences on columns

        >>> s1.compare(s2)
          self other
        1    b     a
        3    d     b

        Stack the differences on indices

        >>> s1.compare(s2, align_axis=0)
        1  self     b
           other    a
        3  self     d
           other    b
        dtype: str

        Keep all original rows

        >>> s1.compare(s2, keep_shape=True)
          self other
        0  NaN   NaN
        1    b     a
        2  NaN   NaN
        3    d     b
        4  NaN   NaN

        Keep all original rows and also all original values

        >>> s1.compare(s2, keep_shape=True, keep_equal=True)
          self other
        0    a     a
        1    b     a
        2    c     c
        3    d     b
        4    e     e
        """

        return super().compare(
            other=other,
            align_axis=align_axis,
            keep_shape=keep_shape,
            keep_equal=keep_equal,
            result_names=result_names,
        )

    def combine(
        self,
        other: Series | Hashable,
        func: Callable[[Hashable, Hashable], Hashable],
        fill_value: Hashable | None = None,
    ) -> Series:
        """
        Combine the Series with a Series or scalar according to `func`.

        Combine the Series and `other` using `func` to perform elementwise
        selection for combined Series.
        `fill_value` is assumed when value is not present at some index
        from one of the two Series being combined.

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

        >>> s1 = pd.Series({"falcon": 330.0, "eagle": 160.0})
        >>> s1
        falcon    330.0
        eagle     160.0
        dtype: float64
        >>> s2 = pd.Series({"falcon": 345.0, "eagle": 200.0, "duck": 30.0})
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
            new_values = np.empty(len(new_index), dtype=object)
            with np.errstate(all="ignore"):
                for i, idx in enumerate(new_index):
                    lv = self.get(idx, fill_value)
                    rv = other.get(idx, fill_value)
                    new_values[i] = func(lv, rv)
        else:
            # Assume that other is a scalar, so apply the function for
            # each element in the Series
            new_index = self.index
            new_values = np.empty(len(new_index), dtype=object)
            with np.errstate(all="ignore"):
                new_values[:] = [func(lv, other) for lv in self._values]
            new_name = self.name

        res_values = self.array._cast_pointwise_result(new_values)
        return self._constructor(
            res_values,
            dtype=res_values.dtype,
            index=new_index,
            name=new_name,
            copy=False,
        )

    def combine_first(self, other) -> Series:
        """
        Update null elements with value in the same location in 'other'.

        Combine two Series objects by filling null values in one Series with
        non-null values from the other Series. Result index will be the union
        of the two indexes.

        Parameters
        ----------
        other : Series
            The value(s) to be used for filling null values.

        Returns
        -------
        Series
            The result of combining the provided Series with the other object.

        See Also
        --------
        Series.combine : Perform element-wise operation on two Series
            using a given function.

        Examples
        --------
        >>> s1 = pd.Series([1, np.nan])
        >>> s2 = pd.Series([3, 4, 5])
        >>> s1.combine_first(s2)
        0    1.0
        1    4.0
        2    5.0
        dtype: float64

        Null values still persist if the location of that null value
        does not exist in `other`

        >>> s1 = pd.Series({"falcon": np.nan, "eagle": 160.0})
        >>> s2 = pd.Series({"eagle": 200.0, "duck": 30.0})
        >>> s1.combine_first(s2)
        duck       30.0
        eagle     160.0
        falcon      NaN
        dtype: float64
        """
        from pandas.core.reshape.concat import concat

        if self.dtype == other.dtype:
            if self.index.equals(other.index):
                return self.mask(self.isna(), other)

        new_index = self.index.union(other.index)

        this = self
        # identify the index subset to keep for each series
        keep_other = other.index.difference(this.index[notna(this)])
        keep_this = this.index.difference(keep_other)

        this = this.reindex(keep_this)
        other = other.reindex(keep_other)

        if this.dtype.kind == "M" and other.dtype.kind != "M":
            # TODO: try to match resos?
            other = to_datetime(other)
            warnings.warn(
                # GH#62931
                "Silently casting non-datetime 'other' to datetime in "
                "Series.combine_first is deprecated and will be removed "
                "in a future version. Explicitly cast before calling "
                "combine_first instead.",
                Pandas4Warning,
                stacklevel=find_stack_level(),
            )

        combined = concat([this, other])
        combined = combined.reindex(new_index)
        return combined.__finalize__(self, method="combine_first")

    def update(self, other: Series | Sequence | Mapping) -> None:
        """
        Modify Series in place using values from passed Series.

        Uses non-NA values from passed Series to make updates. Aligns
        on index.

        Parameters
        ----------
        other : Series, or object coercible into Series
            Other Series that provides values to update the current Series.

        See Also
        --------
        Series.combine : Perform element-wise operation on two Series
            using a given function.
        Series.transform: Modify a Series using a function.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s.update(pd.Series([4, 5, 6]))
        >>> s
        0    4
        1    5
        2    6
        dtype: int64

        >>> s = pd.Series(["a", "b", "c"])
        >>> s.update(pd.Series(["d", "e"], index=[0, 2]))
        >>> s
        0    d
        1    b
        2    e
        dtype: str

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

        ``other`` can also be a non-Series object type
        that is coercible into a Series

        >>> s = pd.Series([1, 2, 3])
        >>> s.update([4, np.nan, 6])
        >>> s
        0    4
        1    2
        2    6
        dtype: int64

        >>> s = pd.Series([1, 2, 3])
        >>> s.update({1: 9})
        >>> s
        0    1
        1    9
        2    3
        dtype: int64
        """
        if not CHAINED_WARNING_DISABLED:
            if sys.getrefcount(
                self
            ) <= REF_COUNT_METHOD and not com.is_local_in_caller_frame(self):
                warnings.warn(
                    _chained_assignment_method_update_msg,
                    ChainedAssignmentError,
                    stacklevel=2,
                )

        if not isinstance(other, Series):
            other = Series(other)

        other = other.reindex_like(self)
        mask = notna(other)

        self._mgr = self._mgr.putmask(mask=mask, new=other)

    # ----------------------------------------------------------------------
    # Reindexing, sorting

    @overload
    def sort_values(
        self,
        *,
        axis: Axis = ...,
        ascending: bool | Sequence[bool] = ...,
        inplace: Literal[False] = ...,
        kind: SortKind = ...,
        na_position: NaPosition = ...,
        ignore_index: bool = ...,
        key: ValueKeyFunc = ...,
    ) -> Series: ...

    @overload
    def sort_values(
        self,
        *,
        axis: Axis = ...,
        ascending: bool | Sequence[bool] = ...,
        inplace: Literal[True],
        kind: SortKind = ...,
        na_position: NaPosition = ...,
        ignore_index: bool = ...,
        key: ValueKeyFunc = ...,
    ) -> None: ...

    @overload
    def sort_values(
        self,
        *,
        axis: Axis = ...,
        ascending: bool | Sequence[bool] = ...,
        inplace: bool = ...,
        kind: SortKind = ...,
        na_position: NaPosition = ...,
        ignore_index: bool = ...,
        key: ValueKeyFunc = ...,
    ) -> Series | None: ...

    def sort_values(
        self,
        *,
        axis: Axis = 0,
        ascending: bool | Sequence[bool] = True,
        inplace: bool = False,
        kind: SortKind = "quicksort",
        na_position: NaPosition = "last",
        ignore_index: bool = False,
        key: ValueKeyFunc | None = None,
    ) -> Series | None:
        """
        Sort by the values.

        Sort a Series in ascending or descending order by some
        criterion.

        Parameters
        ----------
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.
        ascending : bool or list of bools, default True
            If True, sort values in ascending order, otherwise descending.
        inplace : bool, default False
            If True, perform operation in-place.
        kind : {'quicksort', 'mergesort', 'heapsort', 'stable'}, default 'quicksort'
            Choice of sorting algorithm. See also :func:`numpy.sort` for more
            information. 'mergesort' and 'stable' are the only stable  algorithms.
        na_position : {'first' or 'last'}, default 'last'
            Argument 'first' puts NaNs at the beginning, 'last' puts NaNs at
            the end.
        ignore_index : bool, default False
            If True, the resulting axis will be labeled 0, 1, …, n - 1.
        key : callable, optional
            If not None, apply the key function to the series values
            before sorting. This is similar to the `key` argument in the
            builtin :meth:`sorted` function, with the notable difference that
            this `key` function should be *vectorized*. It should expect a
            ``Series`` and return an array-like.

        Returns
        -------
        Series or None
            Series ordered by values or None if ``inplace=True``.

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

        Sort values ascending order (default behavior)

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

        Sort values putting NAs first

        >>> s.sort_values(na_position="first")
        0     NaN
        1     1.0
        2     3.0
        4     5.0
        3    10.0
        dtype: float64

        Sort a series of strings

        >>> s = pd.Series(["z", "b", "d", "a", "c"])
        >>> s
        0    z
        1    b
        2    d
        3    a
        4    c
        dtype: str

        >>> s.sort_values()
        3    a
        1    b
        4    c
        2    d
        0    z
        dtype: str

        Sort using a key function. Your `key` function will be
        given the ``Series`` of values and should return an array-like.

        >>> s = pd.Series(["a", "B", "c", "D", "e"])
        >>> s.sort_values()
        1    B
        3    D
        0    a
        2    c
        4    e
        dtype: str
        >>> s.sort_values(key=lambda x: x.str.lower())
        0    a
        1    B
        2    c
        3    D
        4    e
        dtype: str

        NumPy ufuncs work well here. For example, we can
        sort by the ``sin`` of the value

        >>> s = pd.Series([-4, -2, 0, 2, 4])
        >>> s.sort_values(key=np.sin)
        1   -2
        4    4
        2    0
        0   -4
        3    2
        dtype: int64

        More complicated user-defined functions can be used,
        as long as they expect a Series and return an array-like

        >>> s.sort_values(key=lambda x: (np.tan(x.cumsum())))
        0   -4
        3    2
        4    4
        1   -2
        2    0
        dtype: int64
        """
        inplace = validate_bool_kwarg(inplace, "inplace")
        # Validate the axis parameter
        self._get_axis_number(axis)

        if is_list_like(ascending):
            ascending = cast(Sequence[bool], ascending)
            if len(ascending) != 1:
                raise ValueError(
                    f"Length of ascending ({len(ascending)}) must be 1 for Series"
                )
            ascending = ascending[0]

        ascending = validate_ascending(ascending)

        if na_position not in ["first", "last"]:
            raise ValueError(f"invalid na_position: {na_position}")

        # GH 35922. Make sorting stable by leveraging nargsort
        if key:
            values_to_sort = cast(Series, ensure_key_mapped(self, key))._values
        else:
            values_to_sort = self._values
        sorted_index = nargsort(values_to_sort, kind, bool(ascending), na_position)

        if is_range_indexer(sorted_index, len(sorted_index)):
            if inplace:
                return self._update_inplace(self)
            return self.copy(deep=False)

        result = self._constructor(
            self._values[sorted_index], index=self.index[sorted_index], copy=False
        )

        if ignore_index:
            result.index = default_index(len(sorted_index))

        if not inplace:
            return result.__finalize__(self, method="sort_values")
        self._update_inplace(result)
        return None

    @overload
    def sort_index(
        self,
        *,
        axis: Axis = ...,
        level: IndexLabel = ...,
        ascending: bool | Sequence[bool] = ...,
        inplace: Literal[True],
        kind: SortKind = ...,
        na_position: NaPosition = ...,
        sort_remaining: bool = ...,
        ignore_index: bool = ...,
        key: IndexKeyFunc = ...,
    ) -> None: ...

    @overload
    def sort_index(
        self,
        *,
        axis: Axis = ...,
        level: IndexLabel = ...,
        ascending: bool | Sequence[bool] = ...,
        inplace: Literal[False] = ...,
        kind: SortKind = ...,
        na_position: NaPosition = ...,
        sort_remaining: bool = ...,
        ignore_index: bool = ...,
        key: IndexKeyFunc = ...,
    ) -> Series: ...

    @overload
    def sort_index(
        self,
        *,
        axis: Axis = ...,
        level: IndexLabel = ...,
        ascending: bool | Sequence[bool] = ...,
        inplace: bool = ...,
        kind: SortKind = ...,
        na_position: NaPosition = ...,
        sort_remaining: bool = ...,
        ignore_index: bool = ...,
        key: IndexKeyFunc = ...,
    ) -> Series | None: ...

    def sort_index(
        self,
        *,
        axis: Axis = 0,
        level: IndexLabel | None = None,
        ascending: bool | Sequence[bool] = True,
        inplace: bool = False,
        kind: SortKind = "quicksort",
        na_position: NaPosition = "last",
        sort_remaining: bool = True,
        ignore_index: bool = False,
        key: IndexKeyFunc | None = None,
    ) -> Series | None:
        """
        Sort Series by index labels.

        Returns a new Series sorted by label if `inplace` argument is
        ``False``, otherwise updates the original series and returns None.

        Parameters
        ----------
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.
        level : int, optional
            If not None, sort on values in specified index level(s).
        ascending : bool or list-like of bools, default True
            Sort ascending vs. descending. When the index is a MultiIndex the
            sort direction can be controlled for each level individually.
        inplace : bool, default False
            If True, perform operation in-place.
        kind : {'quicksort', 'mergesort', 'heapsort', 'stable'}, default 'quicksort'
            Choice of sorting algorithm. See also :func:`numpy.sort` for more
            information. 'mergesort' and 'stable' are the only stable algorithms. For
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
        key : callable, optional
            If not None, apply the key function to the index values
            before sorting. This is similar to the `key` argument in the
            builtin :meth:`sorted` function, with the notable difference that
            this `key` function should be *vectorized*. It should expect an
            ``Index`` and return an ``Index`` of the same shape.

        Returns
        -------
        Series or None
            The original Series sorted by the labels or None if ``inplace=True``.

        See Also
        --------
        DataFrame.sort_index: Sort DataFrame by the index.
        DataFrame.sort_values: Sort DataFrame by the value.
        Series.sort_values : Sort Series by the value.

        Examples
        --------
        >>> s = pd.Series(["a", "b", "c", "d"], index=[3, 2, 1, 4])
        >>> s.sort_index()
        1    c
        2    b
        3    a
        4    d
        dtype: str

        Sort Descending

        >>> s.sort_index(ascending=False)
        4    d
        3    a
        2    b
        1    c
        dtype: str

        By default NaNs are put at the end, but use `na_position` to place
        them at the beginning

        >>> s = pd.Series(["a", "b", "c", "d"], index=[3, 2, 1, np.nan])
        >>> s.sort_index(na_position="first")
        NaN     d
         1.0    c
         2.0    b
         3.0    a
        dtype: str

        Specify index level to sort

        >>> arrays = [
        ...     np.array(["qux", "qux", "foo", "foo", "baz", "baz", "bar", "bar"]),
        ...     np.array(["two", "one", "two", "one", "two", "one", "two", "one"]),
        ... ]
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

        Apply a key function before sorting

        >>> s = pd.Series([1, 2, 3, 4], index=["A", "b", "C", "d"])
        >>> s.sort_index(key=lambda x: x.str.lower())
        A    1
        b    2
        C    3
        d    4
        dtype: int64
        """

        return super().sort_index(
            axis=axis,
            level=level,
            ascending=ascending,
            inplace=inplace,
            kind=kind,
            na_position=na_position,
            sort_remaining=sort_remaining,
            ignore_index=ignore_index,
            key=key,
        )

    def argsort(
        self,
        axis: Axis = 0,
        kind: SortKind = "quicksort",
        order: None = None,
        stable: None = None,
    ) -> Series:
        """
        Return the integer indices that would sort the Series values.

        Override ndarray.argsort. Argsorts the value, omitting NA/null values,
        and places the result in the same locations as the non-NA values.

        Parameters
        ----------
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.
        kind : {'mergesort', 'quicksort', 'heapsort', 'stable'}, default 'quicksort'
            Choice of sorting algorithm. See :func:`numpy.sort` for more
            information. 'mergesort' and 'stable' are the only stable algorithms.
        order : None
            Has no effect but is accepted for compatibility with numpy.
        stable : None
            Has no effect but is accepted for compatibility with numpy.

        Returns
        -------
        Series[np.intp]
            Positions of values within the sort order with -1 indicating
            nan values.

        See Also
        --------
        numpy.ndarray.argsort : Returns the indices that would sort this array.

        Examples
        --------
        >>> s = pd.Series([3, 2, 1])
        >>> s.argsort()
        0    2
        1    1
        2    0
        dtype: int64
        """
        if axis != -1:
            # GH#54257 We allow -1 here so that np.argsort(series) works
            self._get_axis_number(axis)

        result = self.array.argsort(kind=kind)

        res = self._constructor(
            result, index=self.index, name=self.name, dtype=np.intp, copy=False
        )
        return res.__finalize__(self, method="argsort")

    def nlargest(
        self, n: int = 5, keep: Literal["first", "last", "all"] = "first"
    ) -> Series:
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
        >>> countries_population = {
        ...     "Italy": 59000000,
        ...     "France": 65000000,
        ...     "Malta": 434000,
        ...     "Maldives": 434000,
        ...     "Brunei": 434000,
        ...     "Iceland": 337000,
        ...     "Nauru": 11300,
        ...     "Tuvalu": 11300,
        ...     "Anguilla": 11300,
        ...     "Montserrat": 5200,
        ... }
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
        Montserrat      5200
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

        >>> s.nlargest(3, keep="last")
        France      65000000
        Italy       59000000
        Brunei        434000
        dtype: int64

        The `n` largest elements where ``n=3`` with all duplicates kept. Note
        that the returned Series has five elements due to the three duplicates.

        >>> s.nlargest(3, keep="all")
        France      65000000
        Italy       59000000
        Malta         434000
        Maldives      434000
        Brunei        434000
        dtype: int64
        """
        return selectn.SelectNSeries(self, n=n, keep=keep).nlargest()

    def nsmallest(
        self, n: int = 5, keep: Literal["first", "last", "all"] = "first"
    ) -> Series:
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
        >>> countries_population = {
        ...     "Italy": 59000000,
        ...     "France": 65000000,
        ...     "Brunei": 434000,
        ...     "Malta": 434000,
        ...     "Maldives": 434000,
        ...     "Iceland": 337000,
        ...     "Nauru": 11300,
        ...     "Tuvalu": 11300,
        ...     "Anguilla": 11300,
        ...     "Montserrat": 5200,
        ... }
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
        Montserrat      5200
        dtype: int64

        The `n` smallest elements where ``n=5`` by default.

        >>> s.nsmallest()
        Montserrat    5200
        Nauru        11300
        Tuvalu       11300
        Anguilla     11300
        Iceland     337000
        dtype: int64

        The `n` smallest elements where ``n=3``. Default `keep` value is
        'first' so Nauru and Tuvalu will be kept.

        >>> s.nsmallest(3)
        Montserrat   5200
        Nauru       11300
        Tuvalu      11300
        dtype: int64

        The `n` smallest elements where ``n=3`` and keeping the last
        duplicates. Anguilla and Tuvalu will be kept since they are the last
        with value 11300 based on the index order.

        >>> s.nsmallest(3, keep="last")
        Montserrat   5200
        Anguilla    11300
        Tuvalu      11300
        dtype: int64

        The `n` smallest elements where ``n=3`` with all duplicates kept. Note
        that the returned Series has four elements due to the three duplicates.

        >>> s.nsmallest(3, keep="all")
        Montserrat   5200
        Nauru       11300
        Tuvalu      11300
        Anguilla    11300
        dtype: int64
        """
        return selectn.SelectNSeries(self, n=n, keep=keep).nsmallest()

    def swaplevel(
        self, i: Level = -2, j: Level = -1, copy: bool | lib.NoDefault = lib.no_default
    ) -> Series:
        """
        Swap levels i and j in a :class:`MultiIndex`.

        Default is to swap the two innermost levels of the index.

        Parameters
        ----------
        i, j : int or str
            Levels of the indices to be swapped. Can pass level name as string.
        copy : bool, default False
            This keyword is now ignored; changing its value will have no
            impact on the method.

            .. deprecated:: 3.0.0

                This keyword is ignored and will be removed in pandas 4.0. Since
                pandas 3.0, this method always returns a new object using a lazy
                copy mechanism that defers copies until necessary
                (Copy-on-Write). See the `user guide on Copy-on-Write
                <https://pandas.pydata.org/docs/dev/user_guide/copy_on_write.html>`__
                for more details.

        Returns
        -------
        Series
            Series with levels swapped in MultiIndex.

        See Also
        --------
        DataFrame.swaplevel : Swap levels i and j in a :class:`DataFrame`.
        Series.reorder_levels : Rearrange index levels using input order.
        MultiIndex.swaplevel : Swap levels i and j in a :class:`MultiIndex`.

        Examples
        --------
        >>> s = pd.Series(
        ...     ["A", "B", "A", "C"],
        ...     index=[
        ...         ["Final exam", "Final exam", "Coursework", "Coursework"],
        ...         ["History", "Geography", "History", "Geography"],
        ...         ["January", "February", "March", "April"],
        ...     ],
        ... )
        >>> s
        Final exam  History    January     A
                    Geography  February    B
        Coursework  History    March       A
                    Geography  April       C
        dtype: str

        In the following example, we will swap the levels of the indices.
        Here, we will swap the levels column-wise, but levels can be swapped row-wise
        in a similar manner. Note that column-wise is the default behavior.
        By not supplying any arguments for i and j, we swap the last and second to
        last indices.

        >>> s.swaplevel()
        Final exam  January   History       A
                    February  Geography     B
        Coursework  March     History       A
                    April     Geography     C
        dtype: str

        By supplying one argument, we can choose which index to swap the last
        index with. We can for example swap the first index with the last one as
        follows.

        >>> s.swaplevel(0)
        January     History     Final exam      A
        February    Geography   Final exam      B
        March       History     Coursework      A
        April       Geography   Coursework      C
        dtype: str

        We can also define explicitly which indices we want to swap by supplying values
        for both i and j. Here, we for example swap the first and second indices.

        >>> s.swaplevel(0, 1)
        History     Final exam  January         A
        Geography   Final exam  February        B
        History     Coursework  March           A
        Geography   Coursework  April           C
        dtype: str
        """
        self._check_copy_deprecation(copy)
        assert isinstance(self.index, MultiIndex)
        result = self.copy(deep=False)
        result.index = self.index.swaplevel(i, j)
        return result

    def reorder_levels(self, order: Sequence[Level]) -> Series:
        """
        Rearrange index levels using input order.

        May not drop or duplicate levels.

        Parameters
        ----------
        order : list of int representing new level order
            Reference level by number or key.

        Returns
        -------
        Series
            Type of caller with index as MultiIndex (new object).

        See Also
        --------
        DataFrame.reorder_levels : Rearrange index or column levels using
            input ``order``.

        Examples
        --------
        >>> arrays = [
        ...     np.array(["dog", "dog", "cat", "cat", "bird", "bird"]),
        ...     np.array(["white", "black", "white", "black", "white", "black"]),
        ... ]
        >>> s = pd.Series([1, 2, 3, 3, 5, 2], index=arrays)
        >>> s
        dog   white    1
              black    2
        cat   white    3
              black    3
        bird  white    5
              black    2
        dtype: int64
        >>> s.reorder_levels([1, 0])
        white  dog     1
        black  dog     2
        white  cat     3
        black  cat     3
        white  bird    5
        black  bird    2
        dtype: int64
        """
        if not isinstance(self.index, MultiIndex):  # pragma: no cover
            raise Exception("Can only reorder levels on a hierarchical axis.")

        result = self.copy(deep=False)
        assert isinstance(result.index, MultiIndex)
        result.index = result.index.reorder_levels(order)
        return result

    def explode(self, ignore_index: bool = False) -> Series:
        """
        Transform each element of a list-like to a row.

        Parameters
        ----------
        ignore_index : bool, default False
            If True, the resulting index will be labeled 0, 1, …, n - 1.

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
        This routine will explode list-likes including lists, tuples, sets,
        Series, and np.ndarray. The result dtype of the subset rows will
        be object. Scalars will be returned unchanged, and empty list-likes will
        result in an np.nan for that row. In addition, the ordering of elements in
        the output will be non-deterministic when exploding sets.

        Reference :ref:`the user guide <reshaping.explode>` for more examples.

        Examples
        --------
        >>> s = pd.Series([[1, 2, 3], "foo", [], [3, 4]])
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
        if isinstance(self.dtype, ExtensionDtype):
            values, counts = self._values._explode()
        elif len(self) and is_object_dtype(self.dtype):
            values, counts = reshape.explode(np.asarray(self._values))
        else:
            result = self.copy()
            return result.reset_index(drop=True) if ignore_index else result

        if ignore_index:
            index: Index = default_index(len(values))
        else:
            index = self.index.repeat(counts)

        return self._constructor(values, index=index, name=self.name, copy=False)

    def unstack(
        self,
        level: IndexLabel = -1,
        fill_value: Hashable | None = None,
        sort: bool = True,
    ) -> DataFrame:
        """
        Unstack, also known as pivot, Series with MultiIndex to produce DataFrame.

        Parameters
        ----------
        level : int, str, or list of these, default last level
            Level(s) to unstack, can pass level name.
        fill_value : scalar value, default None
            Value to use when replacing NaN values.
        sort : bool, default True
            Sort the level(s) in the resulting MultiIndex columns.

        Returns
        -------
        DataFrame
            Unstacked Series.

        See Also
        --------
        DataFrame.unstack : Pivot the MultiIndex of a DataFrame.

        Notes
        -----
        Reference :ref:`the user guide <reshaping.stacking>` for more examples.

        Examples
        --------
        >>> s = pd.Series(
        ...     [1, 2, 3, 4],
        ...     index=pd.MultiIndex.from_product([["one", "two"], ["a", "b"]]),
        ... )
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

        return unstack(self, level, fill_value, sort)

    # ----------------------------------------------------------------------
    # function application

    def map(
        self,
        func: Callable | Mapping | Series | None = None,
        na_action: Literal["ignore"] | None = None,
        engine: Callable | None = None,
        **kwargs,
    ) -> Series:
        """
        Map values of Series according to an input mapping or function.

        Used for substituting each value in a Series with another value,
        that may be derived from a function, a ``dict`` or
        a :class:`Series`.

        Parameters
        ----------
        func : function, collections.abc.Mapping subclass or Series
            Function or mapping correspondence.
        na_action : {None, 'ignore'}, default None
            If 'ignore', propagate NaN values, without passing them to the
            mapping correspondence.
        engine : decorator, optional
            Choose the execution engine to use to run the function. Only used for
            functions. If ``map`` is called with a mapping or ``Series``, an
            exception will be raised. If ``engine`` is not provided the function will
            be executed by the regular Python interpreter.

            Options include JIT compilers such as Numba, Bodo or Blosc2, which in some
            cases can speed up the execution. To use an executor you can provide the
            decorators ``numba.jit``, ``numba.njit``, ``bodo.jit`` or ``blosc2.jit``.
            You can also provide the decorator with parameters, like
            ``numba.jit(nogit=True)``.

            Not all functions can be executed with all execution engines. In general,
            JIT compilers will require type stability in the function (no variable
            should change data type during the execution). And not all pandas and
            NumPy APIs are supported. Check the engine documentation for limitations.

            .. versionadded:: 3.0.0

        **kwargs
            Additional keyword arguments to pass as keywords arguments to
            `arg`.

            .. versionadded:: 3.0.0

        Returns
        -------
        Series
            Same index as caller.

        See Also
        --------
        Series.apply : For applying more complex functions on a Series.
        Series.replace: Replace values given in `to_replace` with `value`.
        DataFrame.apply : Apply a function row-/column-wise.
        DataFrame.map : Apply a function elementwise on a whole DataFrame.

        Notes
        -----
        When ``arg`` is a dictionary, values in Series that are not in the
        dictionary (as keys) are converted to ``NaN``. However, if the
        dictionary is a ``dict`` subclass that defines ``__missing__`` (i.e.
        provides a method for default values), then this default is used
        rather than ``NaN``.

        Examples
        --------
        >>> s = pd.Series(["cat", "dog", np.nan, "rabbit"])
        >>> s
        0      cat
        1      dog
        2      NaN
        3   rabbit
        dtype: str

        ``map`` accepts a ``dict`` or a ``Series``. Values that are not found
        in the ``dict`` are converted to ``NaN``, unless the dict has a default
        value (e.g. ``defaultdict``):

        >>> s.map({"cat": "kitten", "dog": "puppy"})
        0   kitten
        1    puppy
        2      NaN
        3      NaN
        dtype: str

        It also accepts a function:

        >>> s.map("I am a {}".format)
        0       I am a cat
        1       I am a dog
        2       I am a nan
        3    I am a rabbit
        dtype: str

        To avoid applying the function to missing values (and keep them as
        ``NaN``) ``na_action='ignore'`` can be used:

        >>> s.map("I am a {}".format, na_action="ignore")
        0     I am a cat
        1     I am a dog
        2            NaN
        3  I am a rabbit
        dtype: str

        For categorical data, the function is only applied to the categories:

        >>> s = pd.Series(list("cabaa"))
        >>> s.map(print)
        c
        a
        b
        a
        a
        0    None
        1    None
        2    None
        3    None
        4    None
        dtype: object

        >>> s_cat = s.astype("category")
        >>> s_cat.map(print)  # function called once per unique category
        a
        b
        c
        0    None
        1    None
        2    None
        3    None
        4    None
        dtype: object
        """
        if func is None:
            if "arg" in kwargs:
                # `.map(arg=my_func)`
                func = kwargs.pop("arg")
                # https://github.com/pandas-dev/pandas/pull/61264
                warnings.warn(
                    "The parameter `arg` has been renamed to `func`, and it "
                    "will stop being supported in a future version of pandas.",
                    Pandas4Warning,
                    stacklevel=find_stack_level(),
                )
            else:
                raise ValueError("The `func` parameter is required")

        if engine is not None:
            if not callable(func):
                raise ValueError(
                    "The engine argument can only be specified when func is a function"
                )
            if not hasattr(engine, "__pandas_udf__"):
                raise ValueError(f"Not a valid engine: {engine!r}")
            result = engine.__pandas_udf__.map(  # type: ignore[attr-defined]
                data=self,
                func=func,
                args=(),
                kwargs=kwargs,
                decorator=engine,
                skip_na=na_action == "ignore",
            )
            if not isinstance(result, Series):
                result = Series(result, index=self.index, name=self.name)
            return result.__finalize__(self, method="map")

        if callable(func):
            func = functools.partial(func, **kwargs)
        new_values = self._map_values(func, na_action=na_action)
        return self._constructor(new_values, index=self.index, copy=False).__finalize__(
            self, method="map"
        )

    def _gotitem(self, key, ndim, subset=None) -> Self:
        """
        Sub-classes to define. Return a sliced object.

        Parameters
        ----------
        key : string / list of selections
        ndim : {1, 2}
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

    def aggregate(self, func=None, axis: Axis = 0, *args, **kwargs):
        """
        Aggregate using one or more operations over the specified axis.

        Parameters
        ----------
        func : function, str, list or dict
            Function to use for aggregating the data. If a function, must either
            work when passed a Series or when passed to Series.apply.

            Accepted combinations are:

            - function
            - string function name
            - list of functions and/or function names, e.g. ``[np.sum, 'mean']``
            - dict of axis labels -> functions, function names or list of such.
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.
        *args
            Positional arguments to pass to `func`.
        **kwargs
            Keyword arguments to pass to `func`.

        Returns
        -------
        scalar, Series or DataFrame
            The return can be:

            * scalar : when Series.agg is called with single function
            * Series : when DataFrame.agg is called with a single function
            * DataFrame : when DataFrame.agg is called with several functions

        See Also
        --------
        Series.apply : Invoke function on a Series.
        Series.transform : Transform function producing a Series with like indexes.

        Notes
        -----
        The aggregation operations are always performed over an axis, either the
        index (default) or the column axis. This behavior is different from
        `numpy` aggregation functions (`mean`, `median`, `prod`, `sum`, `std`,
        `var`), where the default is to compute the aggregation of the flattened
        array, e.g., ``numpy.mean(arr_2d)`` as opposed to
        ``numpy.mean(arr_2d, axis=0)``.

        `agg` is an alias for `aggregate`. Use the alias.

        Functions that mutate the passed object can produce unexpected
        behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
        for more details.

        A passed user-defined-function will be passed a Series for evaluation.

        If ``func`` defines an index relabeling, ``axis`` must be ``0`` or ``index``.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4])
        >>> s
        0    1
        1    2
        2    3
        3    4
        dtype: int64

        >>> s.agg("min")
        1

        >>> s.agg(["min", "max"])
        min   1
        max   4
        dtype: int64
        """

        # Validate the axis parameter
        self._get_axis_number(axis)

        # if func is None, will switch to user-provided "named aggregation" kwargs
        if func is None:
            func = dict(kwargs.items())

        op = SeriesApply(self, func, args=args, kwargs=kwargs)
        result = op.agg()
        return result

    agg = aggregate

    def transform(
        self, func: AggFuncType, axis: Axis = 0, *args, **kwargs
    ) -> DataFrame | Series:
        """
        Call ``func`` on self producing a Series with the same axis shape as self.

        Parameters
        ----------
        func : function, str, list-like or dict-like
            Function to use for transforming the data. If a function, must either
            work when passed a Series or when passed to Series.apply. If func
            is both list-like and dict-like, dict-like behavior takes precedence.

            Accepted combinations are:

            - function
            - string function name
            - list-like of functions and/or function names, e.g. ``[np.exp, 'sqrt']``
            - dict-like of axis labels -> functions, function names or list-like of such

        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.

        *args
            Positional arguments to pass to `func`.
        **kwargs
            Keyword arguments to pass to `func`.

        Returns
        -------
        Series
            A Series that must have the same length as self.

        Raises
        ------
        ValueError : If the returned Series has a different length than self.

        See Also
        --------
        Series.agg : Only perform aggregating type operations.
        Series.apply : Invoke function on a Series.

        Notes
        -----
        Functions that mutate the passed object can produce unexpected
        behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
        for more details.

        Examples
        --------
        >>> df = pd.DataFrame({"A": range(3), "B": range(1, 4)})
        >>> df
        A  B
        0  0  1
        1  1  2
        2  2  3
        >>> df.transform(lambda x: x + 1)
        A  B
        0  1  2
        1  2  3
        2  3  4

        Even though the resulting Series must have the same length as the
        input Series, it is possible to provide several input functions:

        >>> s = pd.Series(range(3))
        >>> s
        0    0
        1    1
        2    2
        dtype: int64
        >>> s.transform([np.sqrt, np.exp])
            sqrt        exp
        0  0.000000   1.000000
        1  1.000000   2.718282
        2  1.414214   7.389056

        You can call transform on a GroupBy object:

        >>> df = pd.DataFrame(
        ...     {
        ...         "Date": [
        ...             "2015-05-08",
        ...             "2015-05-07",
        ...             "2015-05-06",
        ...             "2015-05-05",
        ...             "2015-05-08",
        ...             "2015-05-07",
        ...             "2015-05-06",
        ...             "2015-05-05",
        ...         ],
        ...         "Data": [5, 8, 6, 1, 50, 100, 60, 120],
        ...     }
        ... )
        >>> df
                Date  Data
        0  2015-05-08     5
        1  2015-05-07     8
        2  2015-05-06     6
        3  2015-05-05     1
        4  2015-05-08    50
        5  2015-05-07   100
        6  2015-05-06    60
        7  2015-05-05   120
        >>> df.groupby("Date")["Data"].transform("sum")
        0     55
        1    108
        2     66
        3    121
        4     55
        5    108
        6     66
        7    121
        Name: Data, dtype: int64

        >>> df = pd.DataFrame(
        ...     {
        ...         "c": [1, 1, 1, 2, 2, 2, 2],
        ...         "type": ["m", "n", "o", "m", "m", "n", "n"],
        ...     }
        ... )
        >>> df
        c type
        0  1    m
        1  1    n
        2  1    o
        3  2    m
        4  2    m
        5  2    n
        6  2    n
        >>> df["size"] = df.groupby("c")["type"].transform(len)
        >>> df
        c type size
        0  1    m    3
        1  1    n    3
        2  1    o    3
        3  2    m    4
        4  2    m    4
        5  2    n    4
        6  2    n    4
        """
        # Validate axis argument
        self._get_axis_number(axis)
        ser = self.copy(deep=False)
        result = SeriesApply(ser, func=func, args=args, kwargs=kwargs).transform()
        return result

    def apply(
        self,
        func: AggFuncType,
        args: tuple[Any, ...] = (),
        *,
        by_row: Literal[False, "compat"] = "compat",
        **kwargs,
    ) -> DataFrame | Series:
        """
        Invoke function on values of Series.

        Can be ufunc (a NumPy function that applies to the entire Series)
        or a Python function that only works on single values.

        Parameters
        ----------
        func : function
            Python function or NumPy ufunc to apply.
        args : tuple
            Positional arguments passed to func after the series value.
        by_row : False or "compat", default "compat"
            If ``"compat"`` and func is a callable, func will be passed each element of
            the Series, like ``Series.map``. If func is a list or dict of
            callables, will first try to translate each func into pandas methods. If
            that doesn't work, will try call to apply again with ``by_row="compat"``
            and if that fails, will call apply again with ``by_row=False``
            (backward compatible).
            If False, the func will be passed the whole Series at once.

            ``by_row`` has no effect when ``func`` is a string.

            .. versionadded:: 2.1.0
        **kwargs
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

        Notes
        -----
        Functions that mutate the passed object can produce unexpected
        behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
        for more details.

        Examples
        --------
        Create a series with typical summer temperatures for each city.

        >>> s = pd.Series([20, 21, 12], index=["London", "New York", "Helsinki"])
        >>> s
        London      20
        New York    21
        Helsinki    12
        dtype: int64

        Square the values by defining a function and passing it as an
        argument to ``apply()``.

        >>> def square(x):
        ...     return x**2
        >>> s.apply(square)
        London      400
        New York    441
        Helsinki    144
        dtype: int64

        Square the values by passing an anonymous function as an
        argument to ``apply()``.

        >>> s.apply(lambda x: x**2)
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
        return SeriesApply(
            self,
            func,
            by_row=by_row,
            args=args,
            kwargs=kwargs,
        ).apply()

    def _reindex_indexer(
        self,
        new_index: Index | None,
        indexer: npt.NDArray[np.intp] | None,
    ) -> Series:
        # Note: new_index is None iff indexer is None
        # if not None, indexer is np.intp
        if indexer is None and (
            new_index is None or new_index.names == self.index.names
        ):
            return self.copy(deep=False)

        new_values = algorithms.take_nd(
            self._values, indexer, allow_fill=True, fill_value=None
        )
        return self._constructor(new_values, index=new_index, copy=False)

    def _needs_reindex_multi(self, axes, method, level) -> bool:
        """
        Check if we do need a multi reindex; this is for compat with
        higher dims.
        """
        return False

    @overload
    def rename(
        self,
        index: Renamer | Hashable | None = ...,
        *,
        axis: Axis | None = ...,
        copy: bool | lib.NoDefault = ...,
        inplace: Literal[True],
        level: Level | None = ...,
        errors: IgnoreRaise = ...,
    ) -> Series | None: ...

    @overload
    def rename(
        self,
        index: Renamer | Hashable | None = ...,
        *,
        axis: Axis | None = ...,
        copy: bool | lib.NoDefault = ...,
        inplace: Literal[False] = ...,
        level: Level | None = ...,
        errors: IgnoreRaise = ...,
    ) -> Series: ...

    def rename(
        self,
        index: Renamer | Hashable | None = None,
        *,
        axis: Axis | None = None,
        copy: bool | lib.NoDefault = lib.no_default,
        inplace: bool = False,
        level: Level | None = None,
        errors: IgnoreRaise = "ignore",
    ) -> Series | None:
        """
        Alter Series index labels or name.

        Function / dict values must be unique (1-to-1). Labels not contained in
        a dict / Series will be left as-is. Extra labels listed don't throw an
        error.

        Alternatively, change ``Series.name`` with a scalar value.

        See the :ref:`user guide <basics.rename>` for more.

        Parameters
        ----------
        index : scalar, hashable sequence, dict-like or function optional
            Functions or dict-like are transformations to apply to
            the index.
            Scalar or hashable sequence-like will alter the ``Series.name``
            attribute.
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.
        copy : bool, default False
            This keyword is now ignored; changing its value will have no
            impact on the method.

            .. deprecated:: 3.0.0

                This keyword is ignored and will be removed in pandas 4.0. Since
                pandas 3.0, this method always returns a new object using a lazy
                copy mechanism that defers copies until necessary
                (Copy-on-Write). See the `user guide on Copy-on-Write
                <https://pandas.pydata.org/docs/dev/user_guide/copy_on_write.html>`__
                for more details.

        inplace : bool, default False
            Whether to return a new Series. If True the value of copy is ignored.
        level : int or level name, default None
            In case of MultiIndex, only rename labels in the specified level.
        errors : {'ignore', 'raise'}, default 'ignore'
            If 'raise', raise `KeyError` when a `dict-like mapper` or
            `index` contains labels that are not present in the index being transformed.
            If 'ignore', existing keys will be renamed and extra keys will be ignored.

        Returns
        -------
        Series
            A shallow copy with index labels or name altered, or the same object
            if ``inplace=True`` and index is not a dict or callable else None.

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
        >>> s.rename(lambda x: x**2)  # function, changes labels
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
        self._check_copy_deprecation(copy)
        if axis is not None:
            # Make sure we raise if an invalid 'axis' is passed.
            axis = self._get_axis_number(axis)

        if callable(index) or is_dict_like(index):
            # error: Argument 1 to "_rename" of "NDFrame" has incompatible
            # type "Union[Union[Mapping[Any, Hashable], Callable[[Any],
            # Hashable]], Hashable, None]"; expected "Union[Mapping[Any,
            # Hashable], Callable[[Any], Hashable], None]"
            return super()._rename(
                index,  # type: ignore[arg-type]
                inplace=inplace,
                level=level,
                errors=errors,
            )
        else:
            return self._set_name(index, inplace=inplace)

    def set_axis(
        self,
        labels,
        *,
        axis: Axis = 0,
        copy: bool | lib.NoDefault = lib.no_default,
    ) -> Series:
        """
        Assign desired index to given axis.

        .. deprecated:: 3.0.0
            This keyword is ignored and will be removed in pandas 4.0. Since
            pandas 3.0, this method always returns a new object using a lazy
            copy mechanism that defers copies until necessary
            (Copy-on-Write). See the `user guide on Copy-on-Write
            <https://pandas.pydata.org/docs/dev/user_guide/copy_on_write.html>`__
            for more details.

        Indexes for row labels can be changed by assigning a list-like or Index.

        Parameters
        ----------
        labels : list-like or Index
            The values for the new index.
        axis : {0 or 'index'}, default 0
            The axis to update. The value 0 identifies the rows. For `Series`
            this parameter is unused and defaults to 0.
        copy : bool, default False
            This keyword is now ignored; changing its value will have no
            impact on the method.

        Returns
        -------
        Series
            A shallow copy of the object with axis altered to the given index.

        See Also
        --------
        Series.rename_axis : Alter the name of the index.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s
        0    1
        1    2
        2    3
        dtype: int64
        >>> s.set_axis(["a", "b", "c"], axis=0)
        a    1
        b    2
        c    3
        dtype: int64
        """

        return super().set_axis(labels, axis=axis, copy=copy)

    # error: Cannot determine type of 'reindex'

    def reindex(  # type: ignore[override]
        self,
        index=None,
        *,
        axis: Axis | None = None,
        method: ReindexMethod | None = None,
        copy: bool | lib.NoDefault = lib.no_default,
        level: Level | None = None,
        fill_value: Scalar | None = None,
        limit: int | None = None,
        tolerance=None,
    ) -> Series:
        """
        Conform Series to new index with optional filling logic.

        Places NA/NaN in locations having no value in the previous index. A new object
        is produced unless the new index is equivalent to the current one and
        ``copy=False``.

        Parameters
        ----------
        index : scalar, list-like, dict-like or function, optional
            A scalar, list-like, dict-like or functions transformations to
            apply to that axis' values.
        axis : {0 or 'index'}, default 0
            The axis to rename. For `Series` this parameter is unused and defaults to 0.
        method : {{None, 'backfill'/'bfill', 'pad'/'ffill', 'nearest'}}
            Method to use for filling holes in reindexed DataFrame.
            Please note: this is only applicable to DataFrames/Series with a
            monotonically increasing/decreasing index.

            * None (default): don't fill gaps
            * pad / ffill: Propagate last valid observation forward to next
              valid.
            * backfill / bfill: Use next valid observation to fill gap.
            * nearest: Use nearest valid observations to fill gap.

        copy : bool, default False
            This keyword is now ignored; changing its value will have no
            impact on the method.

            .. deprecated:: 3.0.0

                This keyword is ignored and will be removed in pandas 4.0. Since
                pandas 3.0, this method always returns a new object using a lazy
                copy mechanism that defers copies until necessary
                (Copy-on-Write). See the `user guide on Copy-on-Write
                <https://pandas.pydata.org/docs/dev/user_guide/copy_on_write.html>`__
                for more details.

        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level.
        fill_value : scalar, default np.nan
            Value to use for missing values. Defaults to NaN, but can be any
            "compatible" value.
        limit : int, default None
            Maximum number of consecutive elements to forward or backward fill.
        tolerance : optional
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations most
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.

            Tolerance may be a scalar value, which applies the same tolerance
            to all values, or list-like, which applies variable tolerance per
            element. List-like includes list, tuple, array, Series, and must be
            the same size as the index and its dtype must exactly match the
            index's type.

        Returns
        -------
        Series
            Series with changed index.

        See Also
        --------
        DataFrame.set_index : Set row labels.
        DataFrame.reset_index : Remove row labels or move them to new columns.
        DataFrame.reindex_like : Change to same indices as other DataFrame.

        Examples
        --------
        ``DataFrame.reindex`` supports two calling conventions

        * ``(index=index_labels, columns=column_labels, ...)``
        * ``(labels, axis={{'index', 'columns'}}, ...)``

        We *highly* recommend using keyword arguments to clarify your
        intent.

        Create a DataFrame with some fictional data.

        >>> index = ["Firefox", "Chrome", "Safari", "IE10", "Konqueror"]
        >>> columns = ["http_status", "response_time"]
        >>> df = pd.DataFrame(
        ...     [[200, 0.04], [200, 0.02], [404, 0.07], [404, 0.08], [301, 1.0]],
        ...     columns=columns,
        ...     index=index,
        ... )
        >>> df
                   http_status  response_time
        Firefox            200           0.04
        Chrome             200           0.02
        Safari             404           0.07
        IE10               404           0.08
        Konqueror          301           1.00

        Create a new index and reindex the DataFrame. By default
        values in the new index that do not have corresponding
        records in the DataFrame are assigned ``NaN``.

        >>> new_index = ["Safari", "Iceweasel", "Comodo Dragon", "IE10", "Chrome"]
        >>> df.reindex(new_index)
                       http_status  response_time
        Safari               404.0           0.07
        Iceweasel              NaN            NaN
        Comodo Dragon          NaN            NaN
        IE10                 404.0           0.08
        Chrome               200.0           0.02

        We can fill in the missing values by passing a value to
        the keyword ``fill_value``. Because the index is not monotonically
        increasing or decreasing, we cannot use arguments to the keyword
        ``method`` to fill the ``NaN`` values.

        >>> df.reindex(new_index, fill_value=0)
                       http_status  response_time
        Safari                 404           0.07
        Iceweasel                0           0.00
        Comodo Dragon            0           0.00
        IE10                   404           0.08
        Chrome                 200           0.02

        >>> df.reindex(new_index, fill_value="missing")
                      http_status response_time
        Safari                404          0.07
        Iceweasel         missing       missing
        Comodo Dragon     missing       missing
        IE10                  404          0.08
        Chrome                200          0.02

        We can also reindex the columns.

        >>> df.reindex(columns=["http_status", "user_agent"])
                   http_status  user_agent
        Firefox            200         NaN
        Chrome             200         NaN
        Safari             404         NaN
        IE10               404         NaN
        Konqueror          301         NaN

        Or we can use "axis-style" keyword arguments

        >>> df.reindex(["http_status", "user_agent"], axis="columns")
                   http_status  user_agent
        Firefox            200         NaN
        Chrome             200         NaN
        Safari             404         NaN
        IE10               404         NaN
        Konqueror          301         NaN

        To further illustrate the filling functionality in
        ``reindex``, we will create a DataFrame with a
        monotonically increasing index (for example, a sequence
        of dates).

        >>> date_index = pd.date_range("1/1/2010", periods=6, freq="D")
        >>> df2 = pd.DataFrame(
        ...     {"prices": [100, 101, np.nan, 100, 89, 88]}, index=date_index
        ... )
        >>> df2
                    prices
        2010-01-01   100.0
        2010-01-02   101.0
        2010-01-03     NaN
        2010-01-04   100.0
        2010-01-05    89.0
        2010-01-06    88.0

        Suppose we decide to expand the DataFrame to cover a wider
        date range.

        >>> date_index2 = pd.date_range("12/29/2009", periods=10, freq="D")
        >>> df2.reindex(date_index2)
                    prices
        2009-12-29     NaN
        2009-12-30     NaN
        2009-12-31     NaN
        2010-01-01   100.0
        2010-01-02   101.0
        2010-01-03     NaN
        2010-01-04   100.0
        2010-01-05    89.0
        2010-01-06    88.0
        2010-01-07     NaN

        The index entries that did not have a value in the original data frame
        (for example, '2009-12-29') are by default filled with ``NaN``.
        If desired, we can fill in the missing values using one of several
        options.

        For example, to back-propagate the last valid value to fill the ``NaN``
        values, pass ``bfill`` as an argument to the ``method`` keyword.

        >>> df2.reindex(date_index2, method="bfill")
                    prices
        2009-12-29   100.0
        2009-12-30   100.0
        2009-12-31   100.0
        2010-01-01   100.0
        2010-01-02   101.0
        2010-01-03     NaN
        2010-01-04   100.0
        2010-01-05    89.0
        2010-01-06    88.0
        2010-01-07     NaN

        Please note that the ``NaN`` value present in the original DataFrame
        (at index value 2010-01-03) will not be filled by any of the
        value propagation schemes. This is because filling while reindexing
        does not look at DataFrame values, but only compares the original and
        desired indexes. If you do want to fill in the ``NaN`` values present
        in the original DataFrame, use the ``fillna()`` method.

        See the :ref:`user guide <basics.reindexing>` for more.
        """
        return super().reindex(
            index=index,
            method=method,
            level=level,
            fill_value=fill_value,
            limit=limit,
            tolerance=tolerance,
            copy=copy,
        )

    @overload  # type: ignore[override]
    def rename_axis(
        self,
        mapper: IndexLabel | lib.NoDefault = ...,
        *,
        index=...,
        axis: Axis = ...,
        copy: bool | lib.NoDefault = ...,
        inplace: Literal[True],
    ) -> None: ...

    @overload
    def rename_axis(
        self,
        mapper: IndexLabel | lib.NoDefault = ...,
        *,
        index=...,
        axis: Axis = ...,
        copy: bool | lib.NoDefault = ...,
        inplace: Literal[False] = ...,
    ) -> Self: ...

    @overload
    def rename_axis(
        self,
        mapper: IndexLabel | lib.NoDefault = ...,
        *,
        index=...,
        axis: Axis = ...,
        copy: bool | lib.NoDefault = ...,
        inplace: bool = ...,
    ) -> Self | None: ...

    def rename_axis(
        self,
        mapper: IndexLabel | lib.NoDefault = lib.no_default,
        *,
        index=lib.no_default,
        axis: Axis = 0,
        copy: bool | lib.NoDefault = lib.no_default,
        inplace: bool = False,
    ) -> Self | None:
        """
        Set the name of the axis for the index.

        Parameters
        ----------
        mapper : scalar, list-like, optional
            Value to set the axis name attribute.

            Use either ``mapper`` and ``axis`` to
            specify the axis to target with ``mapper``, or ``index``.

        index : scalar, list-like, dict-like or function, optional
            A scalar, list-like, dict-like or functions transformations to
            apply to that axis' values.
        axis : {0 or 'index'}, default 0
            The axis to rename. For `Series` this parameter is unused and defaults to 0.
        copy : bool, default False
            This keyword is now ignored; changing its value will have no
            impact on the method.

            .. deprecated:: 3.0.0

                This keyword is ignored and will be removed in pandas 4.0. Since
                pandas 3.0, this method always returns a new object using a lazy
                copy mechanism that defers copies until necessary
                (Copy-on-Write). See the `user guide on Copy-on-Write
                <https://pandas.pydata.org/docs/dev/user_guide/copy_on_write.html>`__
                for more details.

        inplace : bool, default False
            Modifies the object directly, instead of creating a new Series
            or DataFrame.

        Returns
        -------
        Series, or None
            The same type as the caller or None if ``inplace=True``.

        See Also
        --------
        Series.rename : Alter Series index labels or name.
        DataFrame.rename : Alter DataFrame index labels or name.
        Index.rename : Set new names on index.

        Examples
        --------

        >>> s = pd.Series(["dog", "cat", "monkey"])
        >>> s
        0       dog
        1       cat
        2    monkey
        dtype: str
        >>> s.rename_axis("animal")
        animal
        0    dog
        1    cat
        2    monkey
        dtype: str
        """
        return super().rename_axis(
            mapper=mapper,
            index=index,
            axis=axis,
            inplace=inplace,
            copy=copy,
        )

    @overload
    def drop(
        self,
        labels: IndexLabel | ListLike = ...,
        *,
        axis: Axis = ...,
        index: IndexLabel | ListLike = ...,
        columns: IndexLabel | ListLike = ...,
        level: Level | None = ...,
        inplace: Literal[True],
        errors: IgnoreRaise = ...,
    ) -> None: ...

    @overload
    def drop(
        self,
        labels: IndexLabel | ListLike = ...,
        *,
        axis: Axis = ...,
        index: IndexLabel | ListLike = ...,
        columns: IndexLabel | ListLike = ...,
        level: Level | None = ...,
        inplace: Literal[False] = ...,
        errors: IgnoreRaise = ...,
    ) -> Series: ...

    @overload
    def drop(
        self,
        labels: IndexLabel | ListLike = ...,
        *,
        axis: Axis = ...,
        index: IndexLabel | ListLike = ...,
        columns: IndexLabel | ListLike = ...,
        level: Level | None = ...,
        inplace: bool = ...,
        errors: IgnoreRaise = ...,
    ) -> Series | None: ...

    def drop(
        self,
        labels: IndexLabel | ListLike = None,
        *,
        axis: Axis = 0,
        index: IndexLabel | ListLike = None,
        columns: IndexLabel | ListLike = None,
        level: Level | None = None,
        inplace: bool = False,
        errors: IgnoreRaise = "raise",
    ) -> Series | None:
        """
        Return Series with specified index labels removed.

        Remove elements of a Series based on specifying the index labels.
        When using a multi-index, labels on different levels can be removed
        by specifying the level.

        Parameters
        ----------
        labels : single label or list-like
            Index labels to drop.
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.
        index : single label or list-like
            Redundant for application on Series, but 'index' can be used instead
            of 'labels'.
        columns : single label or list-like
            No change is made to the Series; use 'index' or 'labels' instead.
        level : int or level name, optional
            For MultiIndex, level for which the labels will be removed.
        inplace : bool, default False
            If True, do operation inplace and return None.
        errors : {'ignore', 'raise'}, default 'raise'
            If 'ignore', suppress error and only existing labels are dropped.

        Returns
        -------
        Series or None
            Series with specified index labels removed or None if ``inplace=True``.

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
        >>> s = pd.Series(data=np.arange(3), index=["A", "B", "C"])
        >>> s
        A  0
        B  1
        C  2
        dtype: int64

        Drop labels B and C

        >>> s.drop(labels=["B", "C"])
        A  0
        dtype: int64

        Drop 2nd level label in MultiIndex Series

        >>> midx = pd.MultiIndex(
        ...     levels=[["llama", "cow", "falcon"], ["speed", "weight", "length"]],
        ...     codes=[[0, 0, 0, 1, 1, 1, 2, 2, 2], [0, 1, 2, 0, 1, 2, 0, 1, 2]],
        ... )
        >>> s = pd.Series([45, 200, 1.2, 30, 250, 1.5, 320, 1, 0.3], index=midx)
        >>> s
        llama   speed      45.0
                weight    200.0
                length      1.2
        cow     speed      30.0
                weight    250.0
                length      1.5
        falcon  speed     320.0
                weight      1.0
                length      0.3
        dtype: float64

        >>> s.drop(labels="weight", level=1)
        llama   speed      45.0
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

    def pop(self, item: Hashable) -> Any:
        """
        Return item and drops from series. Raise KeyError if not found.

        Parameters
        ----------
        item : label
            Index of the element that needs to be removed.

        Returns
        -------
        scalar
            Value that is popped from series.

        See Also
        --------
        Series.drop: Drop specified values from Series.
        Series.drop_duplicates: Return Series with duplicate values removed.

        Examples
        --------
        >>> ser = pd.Series([1, 2, 3])

        >>> ser.pop(0)
        1

        >>> ser
        1    2
        2    3
        dtype: int64
        """
        return maybe_unbox_numpy_scalar(super().pop(item=item))

    def info(
        self,
        verbose: bool | None = None,
        buf: IO[str] | None = None,
        max_cols: int | None = None,
        memory_usage: bool | str | None = None,
        show_counts: bool = True,
    ) -> None:
        """
        Print a concise summary of a Series.

        This method prints information about a Series including
        the index dtype, non-NA values and memory usage.

        Parameters
        ----------
        verbose : bool, optional
            Whether to print the full summary. By default, the setting in
            ``pandas.options.display.max_info_columns`` is followed.
        buf : writable buffer, defaults to sys.stdout
            Where to send the output. By default, the output is printed to
            sys.stdout. Pass a writable buffer if you need to further process
            the output.
        max_cols : int, optional
            Unused, exists only for compatibility with DataFrame.info.
        memory_usage : bool, str, optional
            Specifies whether total memory usage of the Series
            elements (including the index) should be displayed. By default,
            this follows the ``pandas.options.display.memory_usage`` setting.

            True always show memory usage. False never shows memory usage.
            A value of 'deep' is equivalent to "True with deep introspection".
            Memory usage is shown in human-readable units (base-2
            representation). Without deep introspection a memory estimation is
            made based in column dtype and number of rows assuming values
            consume the same memory amount for corresponding dtypes. With deep
            memory introspection, a real memory usage calculation is performed
            at the cost of computational resources. See the
            :ref:`Frequently Asked Questions <df-memory-usage>` for more
            details.
        show_counts : bool, optional
            Whether to show the non-null counts. By default, this is shown
            only if the DataFrame is smaller than
            ``pandas.options.display.max_info_rows`` and
            ``pandas.options.display.max_info_columns``. A value of True always
            shows the counts, and False never shows the counts.

        Returns
        -------
        None
            This method prints a summary of a Series and returns None.

        See Also
        --------
        Series.describe: Generate descriptive statistics of Series.
        Series.memory_usage: Memory usage of Series.

        Examples
        --------
        >>> int_values = [1, 2, 3, 4, 5]
        >>> text_values = ["alpha", "beta", "gamma", "delta", "epsilon"]
        >>> s = pd.Series(text_values, index=int_values)
        >>> s.info()
        <class 'pandas.Series'>
        Index: 5 entries, 1 to 5
        Series name: None
        Non-Null Count  Dtype
        --------------  -----
        5 non-null      str
        dtypes: str(1)
        memory usage: 106.0 bytes

        Prints a summary excluding information about its values:

        >>> s.info(verbose=False)
        <class 'pandas.Series'>
        Index: 5 entries, 1 to 5
        dtypes: str(1)
        memory usage: 106.0 bytes

        Pipe output of Series.info to buffer instead of sys.stdout, get
        buffer content and writes to a text file:

        >>> import io
        >>> buffer = io.StringIO()
        >>> s.info(buf=buffer)
        >>> s = buffer.getvalue()
        >>> with open("df_info.txt", "w", encoding="utf-8") as f:  # doctest: +SKIP
        ...     f.write(s)
        260

        The `memory_usage` parameter allows deep introspection mode, specially
        useful for big Series and fine-tune memory optimization:

        >>> random_strings_array = np.random.choice(["a", "b", "c"], 10**6)
        >>> s = pd.Series(np.random.choice(["a", "b", "c"], 10**6))
        >>> s.info()
        <class 'pandas.Series'>
        RangeIndex: 1000000 entries, 0 to 999999
        Series name: None
        Non-Null Count    Dtype
        --------------    -----
        1000000 non-null  str
        dtypes: str(1)
        memory usage: 8.6 MB

        >>> s.info(memory_usage="deep")
        <class 'pandas.Series'>
        RangeIndex: 1000000 entries, 0 to 999999
        Series name: None
        Non-Null Count    Dtype
        --------------    -----
        1000000 non-null  str
        dtypes: str(1)
        memory usage: 8.6 MB
        """
        return SeriesInfo(self, memory_usage).render(
            buf=buf,
            max_cols=max_cols,
            verbose=verbose,
            show_counts=show_counts,
        )

    def memory_usage(self, index: bool = True, deep: bool = False) -> int:
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
        156

        Not including the index gives the size of the rest of the data, which
        is necessarily smaller:

        >>> s.memory_usage(index=False)
        24

        The memory footprint of `object` values is ignored by default:

        >>> s = pd.Series(["a", "b"])
        >>> s.values
        <ArrowStringArray>
        ['a', 'b']
        Length: 2, dtype: str
        >>> s.memory_usage()
        150
        >>> s.memory_usage(deep=True)
        150
        """
        v = self._memory_usage(deep=deep)
        if index:
            v += self.index.memory_usage(deep=deep)
        return v

    def isin(self, values) -> Series:
        """
        Whether elements in Series are contained in `values`.

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
        >>> s = pd.Series(
        ...     ["llama", "cow", "llama", "beetle", "llama", "hippo"], name="animal"
        ... )
        >>> s.isin(["cow", "llama"])
        0     True
        1     True
        2     True
        3    False
        4     True
        5    False
        Name: animal, dtype: bool

        To invert the boolean values, use the ``~`` operator:

        >>> ~s.isin(["cow", "llama"])
        0    False
        1    False
        2    False
        3     True
        4    False
        5     True
        Name: animal, dtype: bool

        Passing a single string as ``s.isin('llama')`` will raise an error. Use
        a list of one element instead:

        >>> s.isin(["llama"])
        0     True
        1    False
        2     True
        3    False
        4     True
        5    False
        Name: animal, dtype: bool

        Strings and integers are distinct and are therefore not comparable:

        >>> pd.Series([1]).isin(["1"])
        0    False
        dtype: bool
        >>> pd.Series([1.1]).isin(["1.1"])
        0    False
        dtype: bool
        """
        result = algorithms.isin(self._values, values)
        return self._constructor(result, index=self.index, copy=False).__finalize__(
            self, method="isin"
        )

    def between(
        self,
        left,
        right,
        inclusive: Literal["both", "neither", "left", "right"] = "both",
    ) -> Series:
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
        inclusive : {"both", "neither", "left", "right"}
            Include boundaries. Whether to set each bound as closed or open.

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

        With `inclusive` set to ``"neither"`` boundary values are excluded:

        >>> s.between(1, 4, inclusive="neither")
        0     True
        1    False
        2    False
        3    False
        4    False
        dtype: bool

        `left` and `right` can be any scalar value:

        >>> s = pd.Series(["Alice", "Bob", "Carol", "Eve"])
        >>> s.between("Anna", "Daniel")
        0    False
        1     True
        2     True
        3    False
        dtype: bool
        """
        if inclusive == "both":
            lmask = self >= left
            rmask = self <= right
        elif inclusive == "left":
            lmask = self >= left
            rmask = self < right
        elif inclusive == "right":
            lmask = self > left
            rmask = self <= right
        elif inclusive == "neither":
            lmask = self > left
            rmask = self < right
        else:
            raise ValueError(
                "Inclusive has to be either string of 'both',"
                "'left', 'right', or 'neither'."
            )

        return lmask & rmask

    def case_when(
        self,
        caselist: list[
            tuple[
                ArrayLike | Callable[[Series], Series | np.ndarray | Sequence[bool]],
                ArrayLike | Scalar | Callable[[Series], Series | np.ndarray],
            ],
        ],
    ) -> Series:
        """
        Replace values where the conditions are True.

        .. versionadded:: 2.2.0

        Parameters
        ----------
        caselist : A list of tuples of conditions and expected replacements
            Takes the form:  ``(condition0, replacement0)``,
            ``(condition1, replacement1)``, ... .
            ``condition`` should be a 1-D boolean array-like object
            or a callable. If ``condition`` is a callable,
            it is computed on the Series
            and should return a boolean Series or array.
            The callable must not change the input Series
            (though pandas doesn`t check it). ``replacement`` should be a
            1-D array-like object, a scalar or a callable.
            If ``replacement`` is a callable, it is computed on the Series
            and should return a scalar or Series. The callable
            must not change the input Series
            (though pandas doesn`t check it).

        Returns
        -------
        Series
            A new Series with values replaced based on the provided conditions.

        See Also
        --------
        Series.mask : Replace values where the condition is True.

        Examples
        --------
        >>> c = pd.Series([6, 7, 8, 9], name="c")
        >>> a = pd.Series([0, 0, 1, 2])
        >>> b = pd.Series([0, 3, 4, 5])

        >>> c.case_when(
        ...     caselist=[
        ...         (a.gt(0), a),  # condition, replacement
        ...         (b.gt(0), b),
        ...     ]
        ... )
        0    6
        1    3
        2    1
        3    2
        Name: c, dtype: int64
        """
        if not isinstance(caselist, list):
            raise TypeError(
                f"The caselist argument should be a list; instead got {type(caselist)}"
            )

        if not caselist:
            raise ValueError(
                "provide at least one boolean condition, "
                "with a corresponding replacement."
            )

        for num, entry in enumerate(caselist):
            if not isinstance(entry, tuple):
                raise TypeError(
                    f"Argument {num} must be a tuple; instead got {type(entry)}."
                )
            if len(entry) != 2:
                raise ValueError(
                    f"Argument {num} must have length 2; "
                    "a condition and replacement; "
                    f"instead got length {len(entry)}."
                )
        caselist = [
            (
                com.apply_if_callable(condition, self),
                com.apply_if_callable(replacement, self),
            )
            for condition, replacement in caselist
        ]
        default = self.copy(deep=False)
        conditions, replacements = zip(*caselist, strict=True)
        common_dtypes = [infer_dtype_from(arg)[0] for arg in [*replacements, default]]
        if len(set(common_dtypes)) > 1:
            common_dtype = find_common_type(common_dtypes)
            updated_replacements = []
            for condition, replacement in zip(conditions, replacements, strict=True):
                if is_scalar(replacement):
                    replacement = construct_1d_arraylike_from_scalar(
                        value=replacement, length=len(condition), dtype=common_dtype
                    )
                elif isinstance(replacement, ABCSeries):
                    replacement = replacement.astype(common_dtype)
                else:
                    replacement = pd_array(replacement, dtype=common_dtype)
                updated_replacements.append(replacement)
            replacements = updated_replacements
            default = default.astype(common_dtype)

        counter = range(len(conditions) - 1, -1, -1)
        for position, condition, replacement in zip(
            counter, reversed(conditions), reversed(replacements), strict=True
        ):
            try:
                default = default.mask(
                    condition, other=replacement, axis=0, inplace=False, level=None
                )
            except Exception as error:
                raise ValueError(
                    f"Failed to apply condition{position} and replacement{position}."
                ) from error
        return default

    # error: Cannot determine type of 'isna'
    def isna(self) -> Series:
        """
        Detect missing values.

        Return a boolean same-sized Series indicating if the values are NA.
        NA values, such as None or :attr:`numpy.NaN`, get mapped to True
        values.
        Everything else gets mapped to False values. Characters such as empty
        strings ``''`` or :attr:`numpy.inf` are not considered NA values.

        Returns
        -------
        Series
            Mask of bool values for each element in Series that
            indicates whether an element is an NA value.

        See Also
        --------
        DataFrame.isna : Detect missing values.
        DataFrame.isnull : Alias of isna.
        Series.notna : Boolean inverse of isna.
        DataFrame.notna : Boolean inverse of isna.
        Series.notnull : Alias of notna.
        DataFrame.notnull : Alias of notna.
        Series.dropna : Omit axes labels with missing values.
        DataFrame.dropna : Omit axes labels with missing values.
        isna : Top-level isna.

        Examples
        --------
        Show which entries in a Series are NA.

        >>> ser = pd.Series([5, 6, np.nan])
        >>> ser
        0    5.0
        1    6.0
        2    NaN
        dtype: float64
        >>> ser.isna()
        0    False
        1    False
        2     True
        dtype: bool
        """
        return NDFrame.isna(self)

    # error: Cannot determine type of 'isna'
    @doc(NDFrame.isna, klass=_shared_doc_kwargs["klass"])
    def isnull(self) -> Series:
        """
        Series.isnull is an alias for Series.isna.
        """
        return super().isnull()

    # error: Cannot determine type of 'notna'
    def notna(self) -> Series:
        """
        Detect existing (non-missing) values.

        Return a boolean same-sized Series indicating if the values are not NA.
        Non-missing values get mapped to True. Characters such as empty
        strings ``''`` or :attr:`numpy.inf` are not considered NA values.
        NA values, such as None or :attr:`numpy.NaN`, get mapped to False
        values.

        Returns
        -------
        Series
            Mask of bool values for each element in Series that
            indicates whether an element is not an NA value.

        See Also
        --------
        Series.isna : Detect missing values.
        DataFrame.isna : Detect missing values.
        Series.isnull : Alias of isna.
        DataFrame.isnull : Alias of isna.
        DataFrame.notna : Boolean inverse of isna.
        DataFrame.notnull : Alias of notna.
        Series.dropna : Omit axes labels with missing values.
        DataFrame.dropna : Omit axes labels with missing values.
        notna : Top-level notna.

        Examples
        --------
        Show which entries in a Series are not NA.

        >>> ser = pd.Series([5, 6, np.nan])
        >>> ser
        0    5.0
        1    6.0
        2    NaN
        dtype: float64
        >>> ser.notna()
        0     True
        1     True
        2    False
        dtype: bool
        """
        return super().notna()

    # error: Cannot determine type of 'notna'
    @doc(NDFrame.notna, klass=_shared_doc_kwargs["klass"])
    def notnull(self) -> Series:
        """
        Series.notnull is an alias for Series.notna.
        """
        return super().notnull()

    @overload
    def dropna(
        self,
        *,
        axis: Axis = ...,
        inplace: Literal[False] = ...,
        how: AnyAll | None = ...,
        ignore_index: bool = ...,
    ) -> Series: ...

    @overload
    def dropna(
        self,
        *,
        axis: Axis = ...,
        inplace: Literal[True],
        how: AnyAll | None = ...,
        ignore_index: bool = ...,
    ) -> None: ...

    def dropna(
        self,
        *,
        axis: Axis = 0,
        inplace: bool = False,
        how: AnyAll | None = None,
        ignore_index: bool = False,
    ) -> Series | None:
        """
        Return a new Series with missing values removed.

        See the :ref:`User Guide <missing_data>` for more on which values are
        considered missing, and how to work with missing data.

        Parameters
        ----------
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.
        inplace : bool, default False
            If True, do operation inplace and return None.
        how : str, optional
            Not in use. Kept for compatibility.
        ignore_index : bool, default ``False``
            If ``True``, the resulting axis will be labeled 0, 1, …, n - 1.

            .. versionadded:: 2.0.0

        Returns
        -------
        Series or None
            Series with NA entries dropped from it or None if ``inplace=True``.

        See Also
        --------
        Series.isna: Indicate missing values.
        Series.notna : Indicate existing (non-missing) values.
        Series.fillna : Replace missing values.
        DataFrame.dropna : Drop rows or columns which contain NA values.
        Index.dropna : Drop missing indices.

        Examples
        --------
        >>> ser = pd.Series([1.0, 2.0, np.nan])
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

        Empty strings are not considered NA values. ``None`` is considered an
        NA value.

        >>> ser = pd.Series([np.nan, 2, pd.NaT, "", None, "I stay"])
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
        ignore_index = validate_bool_kwarg(ignore_index, "ignore_index")
        # Validate the axis parameter
        self._get_axis_number(axis or 0)

        if self._can_hold_na:
            result = remove_na_arraylike(self)
        elif not inplace:
            result = self.copy(deep=False)
        else:
            result = self

        if ignore_index:
            result.index = default_index(len(result))

        if inplace:
            return self._update_inplace(result)
        else:
            return result

    # ----------------------------------------------------------------------
    # Time series-oriented methods

    def to_timestamp(
        self,
        freq: Frequency | None = None,
        how: Literal["s", "e", "start", "end"] = "start",
        copy: bool | lib.NoDefault = lib.no_default,
    ) -> Series:
        """
        Cast to DatetimeIndex of Timestamps, at *beginning* of period.

        This can be changed to the *end* of the period, by specifying `how="e"`.

        Parameters
        ----------
        freq : str, default frequency of PeriodIndex
            Desired frequency.
        how : {'s', 'e', 'start', 'end'}
            Convention for converting period to timestamp; start of period
            vs. end.
        copy : bool, default False
            This keyword is now ignored; changing its value will have no
            impact on the method.

            .. deprecated:: 3.0.0

                This keyword is ignored and will be removed in pandas 4.0. Since
                pandas 3.0, this method always returns a new object using a lazy
                copy mechanism that defers copies until necessary
                (Copy-on-Write). See the `user guide on Copy-on-Write
                <https://pandas.pydata.org/docs/dev/user_guide/copy_on_write.html>`__
                for more details.

        Returns
        -------
        Series with DatetimeIndex
            Series with the PeriodIndex cast to DatetimeIndex.

        See Also
        --------
        Series.to_period: Inverse method to cast DatetimeIndex to PeriodIndex.
        DataFrame.to_timestamp: Equivalent method for DataFrame.

        Examples
        --------
        >>> idx = pd.PeriodIndex(["2023", "2024", "2025"], freq="Y")
        >>> s1 = pd.Series([1, 2, 3], index=idx)
        >>> s1
        2023    1
        2024    2
        2025    3
        Freq: Y-DEC, dtype: int64

        The resulting frequency of the Timestamps is `YearBegin`

        >>> s1 = s1.to_timestamp()
        >>> s1
        2023-01-01    1
        2024-01-01    2
        2025-01-01    3
        Freq: YS-JAN, dtype: int64

        Using `freq` which is the offset that the Timestamps will have

        >>> s2 = pd.Series([1, 2, 3], index=idx)
        >>> s2 = s2.to_timestamp(freq="M")
        >>> s2
        2023-01-31    1
        2024-01-31    2
        2025-01-31    3
        Freq: YE-JAN, dtype: int64
        """
        self._check_copy_deprecation(copy)
        if not isinstance(self.index, PeriodIndex):
            raise TypeError(f"unsupported Type {type(self.index).__name__}")

        new_obj = self.copy(deep=False)
        new_index = self.index.to_timestamp(freq=freq, how=how)
        setattr(new_obj, "index", new_index)
        return new_obj

    def to_period(
        self,
        freq: str | None = None,
        copy: bool | lib.NoDefault = lib.no_default,
    ) -> Series:
        """
        Convert Series from DatetimeIndex to PeriodIndex.

        Parameters
        ----------
        freq : str, default None
            Frequency associated with the PeriodIndex.
        copy : bool, default False
            This keyword is now ignored; changing its value will have no
            impact on the method.

            .. deprecated:: 3.0.0

                This keyword is ignored and will be removed in pandas 4.0. Since
                pandas 3.0, this method always returns a new object using a lazy
                copy mechanism that defers copies until necessary
                (Copy-on-Write). See the `user guide on Copy-on-Write
                <https://pandas.pydata.org/docs/dev/user_guide/copy_on_write.html>`__
                for more details.

        Returns
        -------
        Series
            Series with index converted to PeriodIndex.

        See Also
        --------
        DataFrame.to_period: Equivalent method for DataFrame.
        Series.dt.to_period: Convert DateTime column values.

        Examples
        --------
        >>> idx = pd.DatetimeIndex(["2023", "2024", "2025"])
        >>> s = pd.Series([1, 2, 3], index=idx)
        >>> s = s.to_period()
        >>> s
        2023    1
        2024    2
        2025    3
        Freq: Y-DEC, dtype: int64

        Viewing the index

        >>> s.index
        PeriodIndex(['2023', '2024', '2025'], dtype='period[Y-DEC]')
        """
        self._check_copy_deprecation(copy)
        if not isinstance(self.index, DatetimeIndex):
            raise TypeError(f"unsupported Type {type(self.index).__name__}")

        new_obj = self.copy(deep=False)
        new_index = self.index.to_period(freq=freq)
        setattr(new_obj, "index", new_index)
        return new_obj

    # ----------------------------------------------------------------------
    # Add index
    _AXIS_ORDERS: list[Literal["index", "columns"]] = ["index"]
    _AXIS_LEN = len(_AXIS_ORDERS)
    _info_axis_number: Literal[0] = 0
    _info_axis_name: Literal["index"] = "index"

    index = properties.AxisProperty(
        axis=0,
        doc="""
        The index (axis labels) of the Series.

        The index of a Series is used to label and identify each element of the
        underlying data. The index can be thought of as an immutable ordered set
        (technically a multi-set, as it may contain duplicate labels), and is
        used to index and align data in pandas.

        Returns
        -------
        Index
            The index labels of the Series.

        See Also
        --------
        Series.reindex : Conform Series to new index.
        Index : The base pandas index type.

        Notes
        -----
        For more information on pandas indexing, see the `indexing user guide
        <https://pandas.pydata.org/docs/user_guide/indexing.html>`__.

        Examples
        --------
        To create a Series with a custom index and view the index labels:

        >>> cities = ['Kolkata', 'Chicago', 'Toronto', 'Lisbon']
        >>> populations = [14.85, 2.71, 2.93, 0.51]
        >>> city_series = pd.Series(populations, index=cities)
        >>> city_series.index
        Index(['Kolkata', 'Chicago', 'Toronto', 'Lisbon'], dtype='object')

        To change the index labels of an existing Series:

        >>> city_series.index = ['KOL', 'CHI', 'TOR', 'LIS']
        >>> city_series.index
        Index(['KOL', 'CHI', 'TOR', 'LIS'], dtype='object')
        """,
    )

    # ----------------------------------------------------------------------
    # Accessor Methods
    # ----------------------------------------------------------------------
    str = Accessor("str", StringMethods)
    dt = Accessor("dt", CombinedDatetimelikeProperties)
    cat = Accessor("cat", CategoricalAccessor)
    plot = Accessor("plot", pandas.plotting.PlotAccessor)
    sparse = Accessor("sparse", SparseAccessor)
    struct = Accessor("struct", StructAccessor)
    list = Accessor("list", ListAccessor)

    # ----------------------------------------------------------------------
    # Add plotting methods to Series
    hist = pandas.plotting.hist_series

    # ----------------------------------------------------------------------
    # Template-Based Arithmetic/Comparison Methods

    def _cmp_method(self, other, op):
        res_name = ops.get_op_result_name(self, other)

        if isinstance(other, Series) and not self._indexed_same(other):
            raise ValueError("Can only compare identically-labeled Series objects")

        lvalues = self._values
        rvalues = extract_array(other, extract_numpy=True, extract_range=True)

        res_values = ops.comparison_op(lvalues, rvalues, op)

        return self._construct_result(res_values, name=res_name, other=other)

    def _logical_method(self, other, op):
        res_name = ops.get_op_result_name(self, other)
        self, other = self._align_for_op(other, align_asobject=True)

        lvalues = self._values
        rvalues = extract_array(other, extract_numpy=True, extract_range=True)

        res_values = ops.logical_op(lvalues, rvalues, op)
        return self._construct_result(res_values, name=res_name, other=other)

    def _arith_method(self, other, op):
        self, other = self._align_for_op(other)
        return base.IndexOpsMixin._arith_method(self, other, op)

    def _align_for_op(self, right, align_asobject: bool = False):
        """align lhs and rhs Series"""
        # TODO: Different from DataFrame._align_for_op, list, tuple and ndarray
        # are not coerced here
        # because Series has inconsistencies described in GH#13637
        left = self

        if isinstance(right, Series):
            # avoid repeated alignment
            if not left.index.equals(right.index):
                if align_asobject:
                    if left.dtype not in (object, np.bool_) or right.dtype not in (
                        object,
                        np.bool_,
                    ):
                        pass
                        # GH#52538 no longer cast in these cases
                    else:
                        # to keep original value's dtype for bool ops
                        left = left.astype(object)
                        right = right.astype(object)

                left, right = left.align(right)

        return left, right

    def _binop(self, other: Series, func, level=None, fill_value=None) -> Series:
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
        this = self

        if not self.index.equals(other.index):
            this, other = self.align(other, level=level, join="outer")

        this_vals, other_vals = ops.fill_binop(this._values, other._values, fill_value)

        with np.errstate(all="ignore"):
            result = func(this_vals, other_vals)

        name = ops.get_op_result_name(self, other)

        out = this._construct_result(result, name, other)
        return cast(Series, out)

    def _construct_result(
        self,
        result: ArrayLike | tuple[ArrayLike, ArrayLike],
        name: Hashable,
        other: AnyArrayLike | DataFrame,
    ) -> Series | tuple[Series, Series]:
        """
        Construct an appropriately-labelled Series from the result of an op.

        Parameters
        ----------
        result : ndarray or ExtensionArray
        name : Label
        other : Series, DataFrame or array-like

        Returns
        -------
        Series
            In the case of __divmod__ or __rdivmod__, a 2-tuple of Series.
        """
        if isinstance(result, tuple):
            # produced by divmod or rdivmod

            res1 = self._construct_result(result[0], name=name, other=other)
            res2 = self._construct_result(result[1], name=name, other=other)

            # GH#33427 assertions to keep mypy happy
            assert isinstance(res1, Series)
            assert isinstance(res2, Series)
            return (res1, res2)

        # TODO: result should always be ArrayLike, but this fails for some
        #  JSONArray tests
        dtype = getattr(result, "dtype", None)
        out = self._constructor(result, index=self.index, dtype=dtype, copy=False)
        out = out.__finalize__(self)
        out = out.__finalize__(other)

        # Set the result's name after __finalize__ is called because __finalize__
        #  would set it back to self.name
        out.name = name
        return out

    def _flex_method(self, other, op, *, level=None, fill_value=None, axis: Axis = 0):
        if axis is not None:
            self._get_axis_number(axis)

        res_name = ops.get_op_result_name(self, other)

        if isinstance(other, Series):
            return self._binop(other, op, level=level, fill_value=fill_value)
        elif isinstance(other, (np.ndarray, list, tuple, ExtensionArray)):
            if len(other) != len(self):
                raise ValueError("Lengths must be equal")
            other = self._constructor(other, self.index, copy=False)
            result = self._binop(other, op, level=level, fill_value=fill_value)
            result._name = res_name
            return result
        elif isinstance(other, ABCDataFrame):
            # GH#46179
            raise TypeError(
                f"Series.{op.__name__.strip('_')} does not support a DataFrame "
                f"`other`. Use df.{op.__name__.strip('_')}(ser) instead."
            )
        else:
            if fill_value is not None:
                if isna(other):
                    return op(self, fill_value)
                self = self.fillna(fill_value)

            return op(self, other)

    def eq(
        self,
        other,
        level: Level | None = None,
        fill_value: float | None = None,
        axis: Axis = 0,
    ) -> Series:
        """
        Return Equal to of series and other, element-wise (binary operator `eq`).

        Equivalent to ``series == other``, but with support to substitute a fill_value
        for missing data in either one of the inputs.

        Parameters
        ----------
        other : object
            When a Series is provided, will align on indexes. For all other types,
            will behave the same as ``==`` but with possibly different results due
            to the other arguments.
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level.
        fill_value : None or float value, default None (NaN)
            Fill existing missing (NaN) values, and any new element needed for
            successful Series alignment, with this value before computation.
            If data in both corresponding Series locations is missing
            the result of filling (at that location) will be missing.
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.

        Returns
        -------
        Series
            The result of the operation.

        See Also
        --------
        Series.ge : Return elementwise Greater than or equal to of series and other.
        Series.le : Return elementwise Less than or equal to of series and other.
        Series.gt : Return elementwise Greater than of series and other.
        Series.lt : Return elementwise Less than of series and other.

        Examples
        --------
        >>> a = pd.Series([1, 1, 1, np.nan], index=["a", "b", "c", "d"])
        >>> a
        a    1.0
        b    1.0
        c    1.0
        d    NaN
        dtype: float64
        >>> b = pd.Series([1, np.nan, 1, np.nan], index=["a", "b", "d", "e"])
        >>> b
        a    1.0
        b    NaN
        d    1.0
        e    NaN
        dtype: float64
        >>> a.eq(b, fill_value=0)
        a     True
        b    False
        c    False
        d    False
        e    False
        dtype: bool
        """
        return self._flex_method(
            other, operator.eq, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("ne", "series"))
    def ne(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, operator.ne, level=level, fill_value=fill_value, axis=axis
        )

    def le(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        """
        Return Less than or equal to of series and other, \
        element-wise (binary operator `le`).

        Equivalent to ``series <= other``, but with support to substitute a
        fill_value for missing data in either one of the inputs.

        Parameters
        ----------
        other : object
            When a Series is provided, will align on indexes. For all other types,
            will behave the same as ``==`` but with possibly different results due
            to the other arguments.
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level.
        fill_value : None or float value, default None (NaN)
            Fill existing missing (NaN) values, and any new element needed for
            successful Series alignment, with this value before computation.
            If data in both corresponding Series locations is missing
            the result of filling (at that location) will be missing.
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.

        Returns
        -------
        Series
            The result of the operation.

        See Also
        --------
        Series.ge : Return elementwise Greater than or equal to of series and other.
        Series.lt : Return elementwise Less than of series and other.
        Series.gt : Return elementwise Greater than of series and other.
        Series.eq : Return elementwise equal to of series and other.

        Examples
        --------
        >>> a = pd.Series([1, 1, 1, np.nan, 1], index=['a', 'b', 'c', 'd', 'e'])
        >>> a
        a    1.0
        b    1.0
        c    1.0
        d    NaN
        e    1.0
        dtype: float64
        >>> b = pd.Series([0, 1, 2, np.nan, 1], index=['a', 'b', 'c', 'd', 'f'])
        >>> b
        a    0.0
        b    1.0
        c    2.0
        d    NaN
        f    1.0
        dtype: float64
        >>> a.le(b, fill_value=0)
        a    False
        b     True
        c     True
        d    False
        e    False
        f     True
        dtype: bool
        """
        return self._flex_method(
            other, operator.le, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("lt", "series"))
    def lt(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, operator.lt, level=level, fill_value=fill_value, axis=axis
        )

    def ge(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        """
        Return Greater than or equal to of series and other, \
        element-wise (binary operator `ge`).

        Equivalent to ``series >= other``, but with support to substitute a
        fill_value for missing data in either one of the inputs.

        Parameters
        ----------
        other : object
            When a Series is provided, will align on indexes. For all other types,
            will behave the same as ``==`` but with possibly different results due
            to the other arguments.
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level.
        fill_value : None or float value, default None (NaN)
            Fill existing missing (NaN) values, and any new element needed for
            successful Series alignment, with this value before computation.
            If data in both corresponding Series locations is missing
            the result of filling (at that location) will be missing.
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.

        Returns
        -------
        Series
            The result of the operation.

        See Also
        --------
        Series.gt : Greater than comparison, element-wise.
        Series.le : Less than or equal to comparison, element-wise.
        Series.lt : Less than comparison, element-wise.
        Series.eq : Equal to comparison, element-wise.
        Series.ne : Not equal to comparison, element-wise.

        Examples
        --------
        >>> a = pd.Series([1, 1, 1, np.nan, 1], index=["a", "b", "c", "d", "e"])
        >>> a
        a    1.0
        b    1.0
        c    1.0
        d    NaN
        e    1.0
        dtype: float64
        >>> b = pd.Series([0, 1, 2, np.nan, 1], index=["a", "b", "c", "d", "f"])
        >>> b
        a    0.0
        b    1.0
        c    2.0
        d    NaN
        f    1.0
        dtype: float64
        >>> a.ge(b, fill_value=0)
        a     True
        b     True
        c    False
        d    False
        e     True
        f    False
        dtype: bool
        """
        return self._flex_method(
            other, operator.ge, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("gt", "series"))
    def gt(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, operator.gt, level=level, fill_value=fill_value, axis=axis
        )

    def add(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        """
        Return Addition of series and other, element-wise (binary operator `add`).

        Equivalent to ``series + other``, but with support to substitute a fill_value
        for missing data in either one of the inputs.

        Parameters
        ----------
        other : Series or scalar value
            With which to compute the addition.
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level.
        fill_value : None or float value, default None (NaN)
            Fill existing missing (NaN) values, and any new element needed for
            successful Series alignment, with this value before computation.
            If data in both corresponding Series locations is missing
            the result of filling (at that location) will be missing.
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.

        Returns
        -------
        Series
            The result of the operation.

        See Also
        --------
        Series.radd : Reverse of the Addition operator, see
            `Python documentation
            <https://docs.python.org/3/reference/datamodel.html#emulating-numeric-types>`_
            for more details.

        Examples
        --------
        >>> a = pd.Series([1, 1, 1, np.nan], index=["a", "b", "c", "d"])
        >>> a
        a    1.0
        b    1.0
        c    1.0
        d    NaN
        dtype: float64
        >>> b = pd.Series([1, np.nan, 1, np.nan], index=["a", "b", "d", "e"])
        >>> b
        a    1.0
        b    NaN
        d    1.0
        e    NaN
        dtype: float64
        >>> a.add(b, fill_value=0)
        a    2.0
        b    1.0
        c    1.0
        d    1.0
        e    NaN
        dtype: float64
        """
        return self._flex_method(
            other, operator.add, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("radd", "series"))
    def radd(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, roperator.radd, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("sub", "series"))
    def sub(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, operator.sub, level=level, fill_value=fill_value, axis=axis
        )

    subtract = sub

    @Appender(ops.make_flex_doc("rsub", "series"))
    def rsub(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, roperator.rsub, level=level, fill_value=fill_value, axis=axis
        )

    def mul(
        self,
        other,
        level: Level | None = None,
        fill_value: float | None = None,
        axis: Axis = 0,
    ) -> Series:
        """
        Return Multiplication of series and other, element-wise (binary operator `mul`).

        Equivalent to ``series * other``, but with support to substitute
        a fill_value for missing data in either one of the inputs.

        Parameters
        ----------
        other : Series or scalar value
            With which to compute the multiplication.
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level.
        fill_value : None or float value, default None (NaN)
            Fill existing missing (NaN) values, and any new element needed for
            successful Series alignment, with this value before computation.
            If data in both corresponding Series locations is missing
            the result of filling (at that location) will be missing.
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.

        Returns
        -------
        Series
            The result of the operation.

        See Also
        --------
        Series.rmul : Reverse of the Multiplication operator, see
            `Python documentation
            <https://docs.python.org/3/reference/datamodel.html#emulating-numeric-types>`_
            for more details.

        Examples
        --------
        >>> a = pd.Series([1, 1, 1, np.nan], index=["a", "b", "c", "d"])
        >>> a
        a    1.0
        b    1.0
        c    1.0
        d    NaN
        dtype: float64
        >>> b = pd.Series([1, np.nan, 1, np.nan], index=["a", "b", "d", "e"])
        >>> b
        a    1.0
        b    NaN
        d    1.0
        e    NaN
        dtype: float64
        >>> a.multiply(b, fill_value=0)
        a    1.0
        b    0.0
        c    0.0
        d    0.0
        e    NaN
        dtype: float64
        >>> a.mul(5, fill_value=0)
        a    5.0
        b    5.0
        c    5.0
        d    0.0
        dtype: float64
        """
        return self._flex_method(
            other, operator.mul, level=level, fill_value=fill_value, axis=axis
        )

    multiply = mul

    @Appender(ops.make_flex_doc("rmul", "series"))
    def rmul(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, roperator.rmul, level=level, fill_value=fill_value, axis=axis
        )

    def truediv(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        """
        Return Floating division of series and other, \
        element-wise (binary operator `truediv`).

        Equivalent to ``series / other``, but with support to substitute a
        fill_value for missing data in either one of the inputs.

        Parameters
        ----------
        other : Series or scalar value
            Series with which to compute division.
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level.
        fill_value : None or float value, default None (NaN)
            Fill existing missing (NaN) values, and any new element needed for
            successful Series alignment, with this value before computation.
            If data in both corresponding Series locations is missing
            the result of filling (at that location) will be missing.
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.

        Returns
        -------
        Series
            The result of the operation.

        See Also
        --------
        Series.rtruediv : Reverse of the Floating division operator, see
            `Python documentation
            <https://docs.python.org/3/reference/datamodel.html#emulating-numeric-types>`_
            for more details.

        Examples
        --------
        >>> a = pd.Series([1, 1, 1, np.nan], index=["a", "b", "c", "d"])
        >>> a
        a    1.0
        b    1.0
        c    1.0
        d    NaN
        dtype: float64
        >>> b = pd.Series([1, np.nan, 1, np.nan], index=["a", "b", "d", "e"])
        >>> b
        a    1.0
        b    NaN
        d    1.0
        e    NaN
        dtype: float64
        >>> a.divide(b, fill_value=0)
        a    1.0
        b    inf
        c    inf
        d    0.0
        e    NaN
        dtype: float64
        """
        return self._flex_method(
            other, operator.truediv, level=level, fill_value=fill_value, axis=axis
        )

    div = truediv
    divide = truediv

    @Appender(ops.make_flex_doc("rtruediv", "series"))
    def rtruediv(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, roperator.rtruediv, level=level, fill_value=fill_value, axis=axis
        )

    rdiv = rtruediv

    @Appender(ops.make_flex_doc("floordiv", "series"))
    def floordiv(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, operator.floordiv, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("rfloordiv", "series"))
    def rfloordiv(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, roperator.rfloordiv, level=level, fill_value=fill_value, axis=axis
        )

    def mod(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        """
        Return Modulo of series and other, element-wise (binary operator `mod`).

        Equivalent to ``series % other``, but with support to substitute a
        fill_value for missing data in either one of the inputs.

        Parameters
        ----------
        other : Series or scalar value
            Series with which to compute modulo.
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level.
        fill_value : None or float value, default None (NaN)
            Fill existing missing (NaN) values, and any new element needed for
            successful Series alignment, with this value before computation.
            If data in both corresponding Series locations is missing
            the result of filling (at that location) will be missing.
        axis : {0 or 'index'}
            Unused. Parameter needed for compatibility with DataFrame.

        Returns
        -------
        Series
            The result of the operation.

        See Also
        --------
        Series.rmod : Reverse of the Modulo operator, see
            `Python documentation
            <https://docs.python.org/3/reference/datamodel.html#emulating-numeric-types>`_
            for more details.

        Examples
        --------
        >>> a = pd.Series([1, 1, 1, np.nan], index=["a", "b", "c", "d"])
        >>> a
        a    1.0
        b    1.0
        c    1.0
        d    NaN
        dtype: float64
        >>> b = pd.Series([1, np.nan, 1, np.nan], index=["a", "b", "d", "e"])
        >>> b
        a    1.0
        b    NaN
        d    1.0
        e    NaN
        dtype: float64
        >>> a.mod(b, fill_value=0)
        a    0.0
        b    NaN
        c    NaN
        d    0.0
        e    NaN
        dtype: float64
        """
        return self._flex_method(
            other, operator.mod, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("rmod", "series"))
    def rmod(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, roperator.rmod, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("pow", "series"))
    def pow(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, operator.pow, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("rpow", "series"))
    def rpow(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, roperator.rpow, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("divmod", "series"))
    def divmod(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, divmod, level=level, fill_value=fill_value, axis=axis
        )

    @Appender(ops.make_flex_doc("rdivmod", "series"))
    def rdivmod(self, other, level=None, fill_value=None, axis: Axis = 0) -> Series:
        return self._flex_method(
            other, roperator.rdivmod, level=level, fill_value=fill_value, axis=axis
        )

    # ----------------------------------------------------------------------
    # Reductions

    def _reduce(
        self,
        op,
        # error: Variable "pandas.core.series.Series.str" is not valid as a type
        name: str,  # type: ignore[valid-type]
        *,
        axis: Axis = 0,
        skipna: bool = True,
        numeric_only: bool = False,
        filter_type=None,
        **kwds,
    ):
        """
        Perform a reduction operation.

        If we have an ndarray as a value, then simply perform the operation,
        otherwise delegate to the object.
        """
        delegate = self._values

        if axis is not None:
            self._get_axis_number(axis)

        if isinstance(delegate, ExtensionArray):
            # dispatch to ExtensionArray interface
            result = delegate._reduce(name, skipna=skipna, **kwds)

        else:
            # dispatch to numpy arrays
            if numeric_only and self.dtype.kind not in "iufcb":
                # i.e. not is_numeric_dtype(self.dtype)
                kwd_name = "numeric_only"
                if name in ["any", "all"]:
                    kwd_name = "bool_only"
                # GH#47500 - change to TypeError to match other methods
                raise TypeError(
                    f"Series.{name} does not allow {kwd_name}={numeric_only} "
                    "with non-numeric dtypes."
                )
            result = op(delegate, skipna=skipna, **kwds)

        result = maybe_unbox_numpy_scalar(result)
        return result

    # error: Signature of "any" incompatible with supertype "NDFrame"
    def any(  # type: ignore[override]
        self,
        *,
        axis: Axis = 0,
        bool_only: bool = False,
        skipna: bool = True,
        **kwargs,
    ) -> bool:
        """
        Return whether any element is True, potentially over an axis.

        Returns False unless there is at least one element within a series or
        along a Dataframe axis that is True or equivalent (e.g. non-zero or
        non-empty).

        Parameters
        ----------
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Indicate which axis or axes should be reduced. For `Series` this parameter
            is unused and defaults to 0.

            * 0 / 'index' : reduce the index, return a Series whose index is the
              original column labels.
            * 1 / 'columns' : reduce the columns, return a Series whose index is the
              original index.
            * None : reduce all axes, return a scalar.

        bool_only : bool, default False
            Include only boolean columns. Not implemented for Series.
        skipna : bool, default True
            Exclude NA/null values. If the entire row/column is NA and skipna is
            True, then the result will be False, as for an empty row/column.
            If skipna is False, then NA are treated as True, because these are not
            equal to zero.
        **kwargs : any, default None
            Additional keywords have no effect but might be accepted for
            compatibility with NumPy.

        Returns
        -------
        Series or scalar
            If axis=None, then a scalar boolean is returned.
            Otherwise a Series is returned with index matching the index argument.

        See Also
        --------
        numpy.any : Numpy version of this method.
        Series.any : Return whether any element is True.
        Series.all : Return whether all elements are True.
        DataFrame.any : Return whether any element is True over requested axis.
        DataFrame.all : Return whether all elements are True over requested axis.

        Examples
        --------
        **Series**

        For Series input, the output is a scalar indicating whether any element
        is True.

        >>> pd.Series([False, False]).any()
        False
        >>> pd.Series([True, False]).any()
        True
        >>> pd.Series([], dtype="float64").any()
        False
        >>> pd.Series([np.nan]).any()
        False
        >>> pd.Series([np.nan]).any(skipna=False)
        True

        **DataFrame**

        Whether each column contains at least one True element (the default).

        >>> df = pd.DataFrame({"A": [1, 2], "B": [0, 2], "C": [0, 0]})
        >>> df
           A  B  C
        0  1  0  0
        1  2  2  0

        >>> df.any()
        A     True
        B     True
        C    False
        dtype: bool

        Aggregating over the columns.

        >>> df = pd.DataFrame({"A": [True, False], "B": [1, 2]})
        >>> df
               A  B
        0   True  1
        1  False  2

        >>> df.any(axis="columns")
        0    True
        1    True
        dtype: bool

        >>> df = pd.DataFrame({"A": [True, False], "B": [1, 0]})
        >>> df
               A  B
        0   True  1
        1  False  0

        >>> df.any(axis="columns")
        0    True
        1    False
        dtype: bool

        Aggregating over the entire DataFrame with ``axis=None``.

        >>> df.any(axis=None)
        True

        `any` for an empty DataFrame is an empty Series.

        >>> pd.DataFrame([]).any()
        Series([], dtype: bool)
        """
        nv.validate_logical_func((), kwargs, fname="any")
        validate_bool_kwarg(skipna, "skipna", none_allowed=False)
        return self._reduce(
            nanops.nanany,
            name="any",
            axis=axis,
            numeric_only=bool_only,
            skipna=skipna,
            filter_type="bool",
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="all")
    def all(
        self,
        axis: Axis = 0,
        bool_only: bool = False,
        skipna: bool = True,
        **kwargs,
    ) -> bool:
        """
        Return whether all elements are True, potentially over an axis.

        Returns True unless there at least one element within a series or
        along a Dataframe axis that is False or equivalent (e.g. zero or
        empty).

        Parameters
        ----------
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            Indicate which axis or axes should be reduced. For `Series` this parameter
            is unused and defaults to 0.

            * 0 / 'index' : reduce the index, return a Series whose index is the
              original column labels.
            * 1 / 'columns' : reduce the columns, return a Series whose index is the
              original index.
            * None : reduce all axes, return a scalar.

        bool_only : bool, default False
            Include only boolean columns. Not implemented for Series.
        skipna : bool, default True
            Exclude NA/null values. If the entire row/column is NA and skipna is
            True, then the result will be True, as for an empty row/column.
            If skipna is False, then NA are treated as True, because these are not
            equal to zero.
        **kwargs : any, default None
            Additional keywords have no effect but might be accepted for
            compatibility with NumPy.

        Returns
        -------
        Series or scalar
            If axis=None, then a scalar boolean is returned.
            Otherwise a Series is returned with index matching the index argument.

        See Also
        --------
        Series.all : Return True if all elements are True.
        DataFrame.any : Return True if one (or more) elements are True.

        Examples
        --------
        **Series**

        >>> pd.Series([True, True]).all()
        True
        >>> pd.Series([True, False]).all()
        False
        >>> pd.Series([], dtype="float64").all()
        True
        >>> pd.Series([np.nan]).all()
        True
        >>> pd.Series([np.nan]).all(skipna=False)
        True

        **DataFrames**

        Create a DataFrame from a dictionary.

        >>> df = pd.DataFrame({"col1": [True, True], "col2": [True, False]})
        >>> df
           col1   col2
        0  True   True
        1  True  False

        Default behaviour checks if values in each column all return True.

        >>> df.all()
        col1     True
        col2    False
        dtype: bool

        Specify ``axis='columns'`` to check if values in each row all return True.

        >>> df.all(axis="columns")
        0     True
        1    False
        dtype: bool

        Or ``axis=None`` for whether every value is True.

        >>> df.all(axis=None)
        False
        """
        nv.validate_logical_func((), kwargs, fname="all")
        validate_bool_kwarg(skipna, "skipna", none_allowed=False)
        return self._reduce(
            nanops.nanall,
            name="all",
            axis=axis,
            numeric_only=bool_only,
            skipna=skipna,
            filter_type="bool",
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="min")
    def min(
        self,
        axis: Axis | None = 0,
        skipna: bool = True,
        numeric_only: bool = False,
        **kwargs,
    ):
        """
        Return the minimum of the values over the requested axis.

        If you want the *index* of the minimum, use ``idxmin``.
        This is the equivalent of the ``numpy.ndarray`` method ``argmin``.

        Parameters
        ----------
        axis : {index (0)}
            Axis for the function to be applied on.
            For `Series` this parameter is unused and defaults to 0.

            For DataFrames, specifying ``axis=None`` will apply the aggregation
            across both axes.

            .. versionadded:: 2.0.0

        skipna : bool, default True
            Exclude NA/null values when computing the result.
        numeric_only : bool, default False
            Include only float, int, boolean columns.
        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        scalar or Series (if level specified)
            The minimum of the values in the Series.

        See Also
        --------
        numpy.min : Equivalent numpy function for arrays.
        Series.min : Return the minimum.
        Series.max : Return the maximum.
        Series.idxmin : Return the index of the minimum.
        Series.idxmax : Return the index of the maximum.
        DataFrame.min : Return the minimum over the requested axis.
        DataFrame.max : Return the maximum over the requested axis.
        DataFrame.idxmin : Return the index of the minimum over the requested axis.
        DataFrame.idxmax : Return the index of the maximum over the requested axis.

        Examples
        --------
        >>> idx = pd.MultiIndex.from_arrays(
        ...     [["warm", "warm", "cold", "cold"], ["dog", "falcon", "fish", "spider"]],
        ...     names=["blooded", "animal"],
        ... )
        >>> s = pd.Series([4, 2, 0, 8], name="legs", index=idx)
        >>> s
        blooded  animal
        warm     dog       4
                 falcon    2
        cold     fish      0
                 spider    8
        Name: legs, dtype: int64

        >>> s.min()
        0
        """
        return NDFrame.min(
            self, axis=axis, skipna=skipna, numeric_only=numeric_only, **kwargs
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="max")
    def max(
        self,
        axis: Axis | None = 0,
        skipna: bool = True,
        numeric_only: bool = False,
        **kwargs,
    ):
        """
        Return the maximum of the values over the requested axis.

        If you want the *index* of the maximum, use ``idxmax``.
        This is the equivalent of the ``numpy.ndarray`` method ``argmax``.

        Parameters
        ----------
        axis : {index (0)}
            Axis for the function to be applied on.
            For `Series` this parameter is unused and defaults to 0.

            For DataFrames, specifying ``axis=None`` will apply the aggregation
            across both axes.

            .. versionadded:: 2.0.0

        skipna : bool, default True
            Exclude NA/null values when computing the result.
        numeric_only : bool, default False
            Include only float, int, boolean columns.
        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        scalar or Series (if level specified)
            The maximum of the values in the Series.

        See Also
        --------
        numpy.max : Equivalent numpy function for arrays.
        Series.min : Return the minimum.
        Series.max : Return the maximum.
        Series.idxmin : Return the index of the minimum.
        Series.idxmax : Return the index of the maximum.
        DataFrame.min : Return the minimum over the requested axis.
        DataFrame.max : Return the maximum over the requested axis.
        DataFrame.idxmin : Return the index of the minimum over the requested axis.
        DataFrame.idxmax : Return the index of the maximum over the requested axis.

        Examples
        --------
        >>> idx = pd.MultiIndex.from_arrays(
        ...     [["warm", "warm", "cold", "cold"], ["dog", "falcon", "fish", "spider"]],
        ...     names=["blooded", "animal"],
        ... )
        >>> s = pd.Series([4, 2, 0, 8], name="legs", index=idx)
        >>> s
        blooded  animal
        warm     dog       4
                 falcon    2
        cold     fish      0
                 spider    8
        Name: legs, dtype: int64

        >>> s.max()
        8
        """
        return NDFrame.max(
            self, axis=axis, skipna=skipna, numeric_only=numeric_only, **kwargs
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="sum")
    def sum(
        self,
        axis: Axis | None = None,
        skipna: bool = True,
        numeric_only: bool = False,
        min_count: int = 0,
        **kwargs,
    ):
        """
        Return the sum of the values over the requested axis.

        This is equivalent to the method ``numpy.sum``.

        Parameters
        ----------
        axis : {index (0)}
            Axis for the function to be applied on.
            For `Series` this parameter is unused and defaults to 0.

            .. warning::

                The behavior of DataFrame.sum with ``axis=None`` is deprecated,
                in a future version this will reduce over both axes and return a scalar
                To retain the old behavior, pass axis=0 (or do not pass axis).

            .. versionadded:: 2.0.0

        skipna : bool, default True
            Exclude NA/null values when computing the result.
        numeric_only : bool, default False
            Include only float, int, boolean columns. Not implemented for Series.

        min_count : int, default 0
            The required number of valid values to perform the operation. If fewer than
            ``min_count`` non-NA values are present the result will be NA.
        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        scalar or Series (if level specified)
            Sum of the values for the requested axis.

        See Also
        --------
        numpy.sum : Equivalent numpy function for computing sum.
        Series.mean : Mean of the values.
        Series.median : Median of the values.
        Series.std : Standard deviation of the values.
        Series.var : Variance of the values.
        Series.min : Minimum value.
        Series.max : Maximum value.

        Examples
        --------
        >>> idx = pd.MultiIndex.from_arrays(
        ...     [["warm", "warm", "cold", "cold"], ["dog", "falcon", "fish", "spider"]],
        ...     names=["blooded", "animal"],
        ... )
        >>> s = pd.Series([4, 2, 0, 8], name="legs", index=idx)
        >>> s
        blooded  animal
        warm     dog       4
                 falcon    2
        cold     fish      0
                 spider    8
        Name: legs, dtype: int64

        >>> s.sum()
        14

        By default, the sum of an empty or all-NA Series is ``0``.

        >>> pd.Series([], dtype="float64").sum()  # min_count=0 is the default
        0.0

        This can be controlled with the ``min_count`` parameter. For example, if
        you'd like the sum of an empty series to be NaN, pass ``min_count=1``.

        >>> pd.Series([], dtype="float64").sum(min_count=1)
        nan

        Thanks to the ``skipna`` parameter, ``min_count`` handles all-NA and
        empty series identically.

        >>> pd.Series([np.nan]).sum()
        0.0

        >>> pd.Series([np.nan]).sum(min_count=1)
        nan
        """
        return NDFrame.sum(
            self,
            axis=axis,
            skipna=skipna,
            numeric_only=numeric_only,
            min_count=min_count,
            **kwargs,
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="prod")
    def prod(
        self,
        axis: Axis | None = None,
        skipna: bool = True,
        numeric_only: bool = False,
        min_count: int = 0,
        **kwargs,
    ):
        """
        Return the product of the values over the requested axis.

        By default, missing values are skipped. To include them in the calculation,
        set ``skipna`` parameter to False.

        Parameters
        ----------
        axis : {index (0)}
            Axis for the function to be applied on.
            For `Series` this parameter is unused and defaults to 0.

            .. warning::
                The behavior of DataFrame.prod with ``axis=None`` is deprecated,
                in a future version this will reduce over both axes and return a scalar
                To retain the old behavior, pass axis=0 (or do not pass axis).

            .. versionadded:: 2.0.0
        skipna : bool, default True
            Exclude NA/null values when computing the result.
        numeric_only : bool, default False
            Include only float, int, boolean columns. Not implemented for Series.
        min_count : int, default 0
            The required number of valid values to perform the operation. If fewer than
            ``min_count`` non-NA values are present the result will be NA.
        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        scalar
            Value containing the calculation referenced in the description.

        See Also
        --------
        Series.sum : Return the sum.
        Series.min : Return the minimum.
        Series.max : Return the maximum.
        Series.idxmin : Return the index of the minimum.
        Series.idxmax : Return the index of the maximum.

        DataFrame.sum : Return the sum over the requested axis.
        DataFrame.min : Return the minimum over the requested axis.
        DataFrame.max : Return the maximum over the requested axis.
        DataFrame.idxmin : Return the index of the minimum over the requested axis.
        DataFrame.idxmax : Return the index of the maximum over the requested axis.

        Examples
        --------
        By default, the product of an empty or all-NA Series is ``1``

        >>> pd.Series([], dtype="float64").prod()
        1.0

        This can be controlled with the ``min_count`` parameter

        >>> pd.Series([], dtype="float64").prod(min_count=1)
        nan

        Thanks to the ``skipna`` parameter, ``min_count`` handles all-NA and
        empty series identically.

        >>> pd.Series([np.nan]).prod()
        1.0
        >>> pd.Series([np.nan]).prod(min_count=1)
        nan
        """
        return NDFrame.prod(
            self,
            axis=axis,
            skipna=skipna,
            numeric_only=numeric_only,
            min_count=min_count,
            **kwargs,
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="mean")
    def mean(
        self,
        axis: Axis | None = 0,
        skipna: bool = True,
        numeric_only: bool = False,
        **kwargs,
    ) -> Any:
        """
        Return the mean of the values over the requested axis.

        Parameters
        ----------
        axis : {index (0)}
            Axis for the function to be applied on.
            For `Series` this parameter is unused and defaults to 0.

            For DataFrames, specifying ``axis=None`` will apply the aggregation
            across both axes.

            .. versionadded:: 2.0.0

        skipna : bool, default True
            Exclude NA/null values when computing the result.
        numeric_only : bool, default False
            Include only float, int, boolean columns.
        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        scalar or Series (if level specified)
            Mean of the values for the requested axis.

        See Also
        --------
        numpy.median : Equivalent numpy function for computing median.
        Series.sum : Sum of the values.
        Series.median : Median of the values.
        Series.std : Standard deviation of the values.
        Series.var : Variance of the values.
        Series.min : Minimum value.
        Series.max : Maximum value.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s.mean()
        2.0
        """
        return NDFrame.mean(
            self, axis=axis, skipna=skipna, numeric_only=numeric_only, **kwargs
        )

    @deprecate_nonkeyword_arguments(
        Pandas4Warning, allowed_args=["self"], name="median"
    )
    def median(
        self,
        axis: Axis | None = 0,
        skipna: bool = True,
        numeric_only: bool = False,
        **kwargs,
    ) -> Any:
        """
        Return the median of the values over the requested axis.

        Parameters
        ----------
        axis : {index (0)}
            Axis for the function to be applied on.
            For `Series` this parameter is unused and defaults to 0.

            For DataFrames, specifying ``axis=None`` will apply the aggregation
            across both axes.

            .. versionadded:: 2.0.0

        skipna : bool, default True
            Exclude NA/null values when computing the result.
        numeric_only : bool, default False
            Include only float, int, boolean columns.
        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        scalar or Series (if level specified)
            Median of the values for the requested axis.

        See Also
        --------
        numpy.median : Equivalent numpy function for computing median.
        Series.sum : Sum of the values.
        Series.median : Median of the values.
        Series.std : Standard deviation of the values.
        Series.var : Variance of the values.
        Series.min : Minimum value.
        Series.max : Maximum value.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s.median()
        2.0

        With a DataFrame

        >>> df = pd.DataFrame({"a": [1, 2], "b": [2, 3]}, index=["tiger", "zebra"])
        >>> df
               a   b
        tiger  1   2
        zebra  2   3
        >>> df.median()
        a   1.5
        b   2.5
        dtype: float64

        Using axis=1

        >>> df.median(axis=1)
        tiger   1.5
        zebra   2.5
        dtype: float64

        In this case, `numeric_only` should be set to `True`
        to avoid getting an error.

        >>> df = pd.DataFrame({"a": [1, 2], "b": ["T", "Z"]}, index=["tiger", "zebra"])
        >>> df.median(numeric_only=True)
        a   1.5
        dtype: float64
        """
        return NDFrame.median(
            self, axis=axis, skipna=skipna, numeric_only=numeric_only, **kwargs
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="sem")
    def sem(
        self,
        axis: Axis | None = None,
        skipna: bool = True,
        ddof: int = 1,
        numeric_only: bool = False,
        **kwargs,
    ):
        """
        Return unbiased standard error of the mean over requested axis.

        Normalized by N-1 by default. This can be changed using the ddof argument

        Parameters
        ----------
        axis : {index (0)}
            This parameter is unused and defaults to 0.
        skipna : bool, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA.
        ddof : int, default 1
            Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
            where N represents the number of elements.
        numeric_only : bool, default False
            Include only float, int, boolean columns. Not implemented for Series.
        **kwargs :
            Additional keywords have no effect but might be accepted
            for compatibility with NumPy.

        Returns
        -------
        scalar or Series (if level specified)
            Unbiased standard error of the mean over requested axis.

        See Also
        --------
        scipy.stats.sem : Compute standard error of the mean.
        Series.std : Return sample standard deviation over requested axis.
        Series.var : Return unbiased variance over requested axis.
        Series.mean : Return the mean of the values over the requested axis.
        Series.median : Return the median of the values over the requested axis.
        Series.mode : Return the mode(s) of the Series.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> round(s.sem(), 6)
        0.57735
        """
        return NDFrame.sem(
            self,
            axis=axis,
            skipna=skipna,
            ddof=ddof,
            numeric_only=numeric_only,
            **kwargs,
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="var")
    def var(
        self,
        axis: Axis | None = None,
        skipna: bool = True,
        ddof: int = 1,
        numeric_only: bool = False,
        **kwargs,
    ):
        """
        Return unbiased variance over requested axis.

        Normalized by N-1 by default. This can be changed using the ddof argument.

        Parameters
        ----------
        axis : {index (0)}
            For `Series` this parameter is unused and defaults to 0.

            .. warning::

                The behavior of DataFrame.var with ``axis=None`` is deprecated,
                in a future version this will reduce over both axes and return a scalar
                To retain the old behavior, pass axis=0 (or do not pass axis).

        skipna : bool, default True
            Exclude NA/null values. If an entire row/column is NA, the result
            will be NA.
        ddof : int, default 1
            Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
            where N represents the number of elements.
        numeric_only : bool, default False
            Include only float, int, boolean columns. Not implemented for Series.
        **kwargs :
            Additional keywords passed.

        Returns
        -------
        scalar or Series (if level specified)
            Unbiased variance over requested axis.

        See Also
        --------
        numpy.var : Equivalent function in NumPy.
        Series.std : Returns the standard deviation of the Series.
        DataFrame.var : Returns the variance of the DataFrame.
        DataFrame.std : Return standard deviation of the values over
            the requested axis.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {
        ...         "person_id": [0, 1, 2, 3],
        ...         "age": [21, 25, 62, 43],
        ...         "height": [1.61, 1.87, 1.49, 2.01],
        ...     }
        ... ).set_index("person_id")
        >>> df
                   age  height
        person_id
        0           21    1.61
        1           25    1.87
        2           62    1.49
        3           43    2.01

        >>> df.var()
        age       352.916667
        height      0.056367
        dtype: float64

        Alternatively, ``ddof=0`` can be set to normalize by N instead of N-1:

        >>> df.var(ddof=0)
        age       264.687500
        height      0.042275
        dtype: float64
        """
        return NDFrame.var(
            self,
            axis=axis,
            skipna=skipna,
            ddof=ddof,
            numeric_only=numeric_only,
            **kwargs,
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="std")
    def std(
        self,
        axis: Axis | None = None,
        skipna: bool = True,
        ddof: int = 1,
        numeric_only: bool = False,
        **kwargs,
    ):
        """
        Return sample standard deviation.

        Normalized by N-1 by default. This can be changed using the ddof argument.

        Parameters
        ----------
        axis : {index (0)}
            This parameter is unused and defaults to 0.
        skipna : bool, default True
            Exclude NA/null values. If Series is NA, the result
            will be NA.
        ddof : int, default 1
            Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
            where N represents the number of elements.
        numeric_only : bool, default False
            Not implemented for Series.
        **kwargs :
            Additional keywords have no effect but might be accepted
            for compatibility with NumPy.

        Returns
        -------
        scalar
            Standard deviation over all values in the Series.

        See Also
        --------
        numpy.std : Compute the standard deviation along the specified axis.
        Series.var : Return unbiased variance over requested axis.
        Series.sem : Return unbiased standard error of the mean over requested axis.
        Series.mean : Return the mean of the values over the requested axis.
        Series.median : Return the median of the values over the requested axis.
        Series.mode : Return the mode(s) of the Series.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s.std()
        1.0

        Alternatively, ``ddof=0`` can be set to normalize by $N$ instead of $N-1$:

        >>> s.std(ddof=0)
        0.816496580927726
        """
        return NDFrame.std(
            self,
            axis=axis,
            skipna=skipna,
            ddof=ddof,
            numeric_only=numeric_only,
            **kwargs,
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="skew")
    def skew(
        self,
        axis: Axis | None = 0,
        skipna: bool = True,
        numeric_only: bool = False,
        **kwargs,
    ):
        """
        Return unbiased skew over requested axis.

        Normalized by N-1.

        Parameters
        ----------
        axis : {index (0)}
            This parameter is unused and defaults to 0.
        skipna : bool, default True
            Exclude NA/null values when computing the result.
        numeric_only : bool, default False
            Unused.
        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        scalar
            Unbiased skew of the Series.

        See Also
        --------

        Series.var : Return unbiased variance over requested axis.
        Series.std : Return unbiased standard deviation over requested axis.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3])
        >>> s.skew()
        0.0
        """
        return NDFrame.skew(
            self, axis=axis, skipna=skipna, numeric_only=numeric_only, **kwargs
        )

    @deprecate_nonkeyword_arguments(Pandas4Warning, allowed_args=["self"], name="kurt")
    def kurt(
        self,
        axis: Axis | None = 0,
        skipna: bool = True,
        numeric_only: bool = False,
        **kwargs,
    ):
        """
        Return unbiased kurtosis over requested axis.

        Kurtosis obtained using Fisher's definition of
        kurtosis (kurtosis of normal == 0.0). Normalized by N-1.

        Parameters
        ----------
        axis : {index (0)}
            Axis for the function to be applied on.
            For `Series` this parameter is unused and defaults to 0.

            For DataFrames, specifying ``axis=None`` will apply the aggregation
            across both axes.

            .. versionadded:: 2.0.0

        skipna : bool, default True
            Exclude NA/null values when computing the result.
        numeric_only : bool, default False
            Include only float, int, boolean columns.

        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        scalar
            Unbiased kurtosis.

        See Also
        --------
        Series.skew : Return unbiased skew over requested axis.
        Series.var : Return unbiased variance over requested axis.
        Series.std : Return unbiased standard deviation over requested axis.

        Examples
        --------
        >>> s = pd.Series([1, 2, 2, 3], index=["cat", "dog", "dog", "mouse"])
        >>> s
        cat    1
        dog    2
        dog    2
        mouse  3
        dtype: int64
        >>> s.kurt()
        1.5
        """
        return NDFrame.kurt(
            self, axis=axis, skipna=skipna, numeric_only=numeric_only, **kwargs
        )

    kurtosis = kurt
    product = prod

    def cummin(self, axis: Axis = 0, skipna: bool = True, *args, **kwargs) -> Self:
        """
        Return cumulative minimum over a Series.

        Returns a Series of the same size containing the cumulative
        minimum.

        Parameters
        ----------
        axis : {0 or 'index'}, default 0
            This parameter is unused and defaults to 0.
        skipna : bool, default True
            If the entire series is NA, the result will be NA.
        *args, **kwargs
            Additional keywords have no effect but might be accepted for
            compatibility with NumPy.

        Returns
        -------
        Series
            Return cumulative minimum of the Series.

        See Also
        --------
        core.window.expanding.Expanding.min : Similar functionality
            but ignores ``NaN`` values.
        Series.min : Return the minimum value of the Series.
        Series.cummax : Return cumulative maximum.
        Series.cumsum : Return cumulative sum.
        Series.cumprod : Return cumulative product.

        Examples
        --------
        >>> s = pd.Series([2, np.nan, 5, -1, 0])
        >>> s
        0    2.0
        1    NaN
        2    5.0
        3   -1.0
        4    0.0
        dtype: float64

        By default, NA values are ignored.

        >>> s.cummin()
        0    2.0
        1    NaN
        2    2.0
        3   -1.0
        4   -1.0
        dtype: float64

        To include NA values in the operation, use ``skipna=False``

        >>> s.cummin(skipna=False)
        0    2.0
        1    NaN
        2    NaN
        3    NaN
        4    NaN
        dtype: float64
        """
        return NDFrame.cummin(self, axis, skipna, *args, **kwargs)

    def cummax(self, axis: Axis = 0, skipna: bool = True, *args, **kwargs) -> Self:
        """
        Return cumulative maximum over a Series.

        Returns a Series of the same size containing the cumulative
        maximum.

        Parameters
        ----------
        axis : {0 or 'index'}, default 0
            This parameter is unused and defaults to 0.
        skipna : bool, default True
            Exclude NA/null values. If the series is NA, the result is NA.
        *args, **kwargs
            Additional keywords have no effect but might be accepted for
            compatibility with NumPy.

        Returns
        -------
        Series
            Return cumulative maximum of Series.

        See Also
        --------
        core.window.expanding.Expanding.max : Similar functionality
            but ignores ``NaN`` values.
        Series.max : Return the maximum over a Series.
        Series.cummin : Return cumulative minimum.
        Series.cumsum : Return cumulative sum.
        Series.cumprod : Return cumulative product.

        Examples
        --------
        >>> s = pd.Series([2, np.nan, 5, -1, 0])
        >>> s
        0    2.0
        1    NaN
        2    5.0
        3   -1.0
        4    0.0
        dtype: float64

        By default, NA values are ignored.

        >>> s.cummax()
        0    2.0
        1    NaN
        2    5.0
        3    5.0
        4    5.0
        dtype: float64

        To include NA values in the operation, use ``skipna=False``

        >>> s.cummax(skipna=False)
        0    2.0
        1    NaN
        2    NaN
        3    NaN
        4    NaN
        dtype: float64
        """
        return NDFrame.cummax(self, axis, skipna, *args, **kwargs)

    def cumsum(self, axis: Axis = 0, skipna: bool = True, *args, **kwargs) -> Self:
        """
        Return cumulative sum over a Series.

        Returns a Series of the same size containing the cumulative sum.

        Parameters
        ----------
        axis : {0 or 'index'}, default 0
            This parameter is unused and defaults to 0.
        skipna : bool, default True
            Exclude NA/null values. If entire series is NA, the result will be NA.
        *args, **kwargs
            Additional keywords have no effect but might be accepted for
            compatibility with NumPy.

        Returns
        -------
        Series
            Return cumulative sum of Series.

        See Also
        --------
        core.window.expanding.Expanding.sum : Similar functionality
            but ignores ``NaN`` values.
        Series.sum : Return the sum over Series.
        Series.cummax : Return cumulative maximum.
        Series.cummin : Return cumulative minimum.
        Series.cumprod : Return cumulative product.

        Examples
        --------
        >>> s = pd.Series([2, np.nan, 5, -1, 0])
        >>> s
        0    2.0
        1    NaN
        2    5.0
        3   -1.0
        4    0.0
        dtype: float64

        By default, NA values are ignored.

        >>> s.cumsum()
        0    2.0
        1    NaN
        2    7.0
        3    6.0
        4    6.0
        dtype: float64

        To include NA values in the operation, use ``skipna=False``

        >>> s.cumsum(skipna=False)
        0    2.0
        1    NaN
        2    NaN
        3    NaN
        4    NaN
        dtype: float64
        """
        return NDFrame.cumsum(self, axis, skipna, *args, **kwargs)

    def cumprod(self, axis: Axis = 0, skipna: bool = True, *args, **kwargs) -> Self:
        """
        Return cumulative product over a Series.

        Returns a Series of the same size containing the cumulative
        product.

        Parameters
        ----------
        axis : {0 or 'index'}, default 0
            This parameter is unused and defaults to 0.
        skipna : bool, default True
            Exclude NA/null values. If entire Series is NA, the result will be NA.
        *args, **kwargs
            Additional keywords have no effect but might be accepted for
            compatibility with NumPy.

        Returns
        -------
        Series
            Return cumulative product of Series.

        See Also
        --------
        core.window.expanding.Expanding.prod : Similar functionality
            but ignores ``NaN`` values.
        Series.prod : Return the product over Series.
        Series.cummax : Return cumulative maximum.
        Series.cummin : Return cumulative minimum.
        Series.cumsum : Return cumulative sum.

        Examples
        --------
        >>> s = pd.Series([2, np.nan, 5, -1, 0])
        >>> s
        0    2.0
        1    NaN
        2    5.0
        3   -1.0
        4    0.0
        dtype: float64

        By default, NA values are ignored.

        >>> s.cumprod()
        0     2.0
        1     NaN
        2    10.0
        3   -10.0
        4    -0.0
        dtype: float64

        To include NA values in the operation, use ``skipna=False``

        >>> s.cumprod(skipna=False)
        0    2.0
        1    NaN
        2    NaN
        3    NaN
        4    NaN
        dtype: float64
        """
        return NDFrame.cumprod(self, axis, skipna, *args, **kwargs)
