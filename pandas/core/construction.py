"""
Constructor functions intended to be shared by pd.array, Series.__init__,
and Index.__new__.

These should not depend on core.internals.
"""
from __future__ import annotations

from collections import abc
from typing import TYPE_CHECKING, Any, Optional, Sequence, Union, cast

import numpy as np
import numpy.ma as ma

from pandas._libs import lib
from pandas._libs.tslibs import IncompatibleFrequency, OutOfBoundsDatetime
from pandas._typing import AnyArrayLike, ArrayLike, Dtype, DtypeObj

from pandas.core.dtypes.base import ExtensionDtype, registry
from pandas.core.dtypes.cast import (
    construct_1d_arraylike_from_scalar,
    construct_1d_ndarray_preserving_na,
    construct_1d_object_array_from_listlike,
    infer_dtype_from_scalar,
    maybe_cast_to_datetime,
    maybe_cast_to_integer_array,
    maybe_castable,
    maybe_convert_platform,
    maybe_upcast,
)
from pandas.core.dtypes.common import (
    is_datetime64_ns_dtype,
    is_extension_array_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_iterator,
    is_list_like,
    is_object_dtype,
    is_sparse,
    is_string_dtype,
    is_timedelta64_ns_dtype,
)
from pandas.core.dtypes.generic import (
    ABCExtensionArray,
    ABCIndexClass,
    ABCPandasArray,
    ABCSeries,
)
from pandas.core.dtypes.missing import isna

import pandas.core.common as com

if TYPE_CHECKING:
    from pandas import ExtensionArray, Index, Series


def array(
    data: Union[Sequence[object], AnyArrayLike],
    dtype: Optional[Dtype] = None,
    copy: bool = True,
) -> ExtensionArray:
    """
    Create an array.

    .. versionadded:: 0.24.0

    Parameters
    ----------
    data : Sequence of objects
        The scalars inside `data` should be instances of the
        scalar type for `dtype`. It's expected that `data`
        represents a 1-dimensional array of data.

        When `data` is an Index or Series, the underlying array
        will be extracted from `data`.

    dtype : str, np.dtype, or ExtensionDtype, optional
        The dtype to use for the array. This may be a NumPy
        dtype or an extension type registered with pandas using
        :meth:`pandas.api.extensions.register_extension_dtype`.

        If not specified, there are two possibilities:

        1. When `data` is a :class:`Series`, :class:`Index`, or
           :class:`ExtensionArray`, the `dtype` will be taken
           from the data.
        2. Otherwise, pandas will attempt to infer the `dtype`
           from the data.

        Note that when `data` is a NumPy array, ``data.dtype`` is
        *not* used for inferring the array type. This is because
        NumPy cannot represent all the types of data that can be
        held in extension arrays.

        Currently, pandas will infer an extension dtype for sequences of

        ============================== =====================================
        Scalar Type                    Array Type
        ============================== =====================================
        :class:`pandas.Interval`       :class:`pandas.arrays.IntervalArray`
        :class:`pandas.Period`         :class:`pandas.arrays.PeriodArray`
        :class:`datetime.datetime`     :class:`pandas.arrays.DatetimeArray`
        :class:`datetime.timedelta`    :class:`pandas.arrays.TimedeltaArray`
        :class:`int`                   :class:`pandas.arrays.IntegerArray`
        :class:`float`                 :class:`pandas.arrays.FloatingArray`
        :class:`str`                   :class:`pandas.arrays.StringArray`
        :class:`bool`                  :class:`pandas.arrays.BooleanArray`
        ============================== =====================================

        For all other cases, NumPy's usual inference rules will be used.

        .. versionchanged:: 1.0.0

           Pandas infers nullable-integer dtype for integer data,
           string dtype for string data, and nullable-boolean dtype
           for boolean data.

        .. versionchanged:: 1.2.0

            Pandas now also infers nullable-floating dtype for float-like
            input data

    copy : bool, default True
        Whether to copy the data, even if not necessary. Depending
        on the type of `data`, creating the new array may require
        copying data, even if ``copy=False``.

    Returns
    -------
    ExtensionArray
        The newly created array.

    Raises
    ------
    ValueError
        When `data` is not 1-dimensional.

    See Also
    --------
    numpy.array : Construct a NumPy array.
    Series : Construct a pandas Series.
    Index : Construct a pandas Index.
    arrays.PandasArray : ExtensionArray wrapping a NumPy array.
    Series.array : Extract the array stored within a Series.

    Notes
    -----
    Omitting the `dtype` argument means pandas will attempt to infer the
    best array type from the values in the data. As new array types are
    added by pandas and 3rd party libraries, the "best" array type may
    change. We recommend specifying `dtype` to ensure that

    1. the correct array type for the data is returned
    2. the returned array type doesn't change as new extension types
       are added by pandas and third-party libraries

    Additionally, if the underlying memory representation of the returned
    array matters, we recommend specifying the `dtype` as a concrete object
    rather than a string alias or allowing it to be inferred. For example,
    a future version of pandas or a 3rd-party library may include a
    dedicated ExtensionArray for string data. In this event, the following
    would no longer return a :class:`arrays.PandasArray` backed by a NumPy
    array.

    >>> pd.array(['a', 'b'], dtype=str)
    <PandasArray>
    ['a', 'b']
    Length: 2, dtype: str32

    This would instead return the new ExtensionArray dedicated for string
    data. If you really need the new array to be backed by a  NumPy array,
    specify that in the dtype.

    >>> pd.array(['a', 'b'], dtype=np.dtype("<U1"))
    <PandasArray>
    ['a', 'b']
    Length: 2, dtype: str32

    Finally, Pandas has arrays that mostly overlap with NumPy

      * :class:`arrays.DatetimeArray`
      * :class:`arrays.TimedeltaArray`

    When data with a ``datetime64[ns]`` or ``timedelta64[ns]`` dtype is
    passed, pandas will always return a ``DatetimeArray`` or ``TimedeltaArray``
    rather than a ``PandasArray``. This is for symmetry with the case of
    timezone-aware data, which NumPy does not natively support.

    >>> pd.array(['2015', '2016'], dtype='datetime64[ns]')
    <DatetimeArray>
    ['2015-01-01 00:00:00', '2016-01-01 00:00:00']
    Length: 2, dtype: datetime64[ns]

    >>> pd.array(["1H", "2H"], dtype='timedelta64[ns]')
    <TimedeltaArray>
    ['0 days 01:00:00', '0 days 02:00:00']
    Length: 2, dtype: timedelta64[ns]

    Examples
    --------
    If a dtype is not specified, pandas will infer the best dtype from the values.
    See the description of `dtype` for the types pandas infers for.

    >>> pd.array([1, 2])
    <IntegerArray>
    [1, 2]
    Length: 2, dtype: Int64

    >>> pd.array([1, 2, np.nan])
    <IntegerArray>
    [1, 2, <NA>]
    Length: 3, dtype: Int64

    >>> pd.array([1.1, 2.2])
    <FloatingArray>
    [1.1, 2.2]
    Length: 2, dtype: Float64

    >>> pd.array(["a", None, "c"])
    <StringArray>
    ['a', <NA>, 'c']
    Length: 3, dtype: string

    >>> pd.array([pd.Period('2000', freq="D"), pd.Period("2000", freq="D")])
    <PeriodArray>
    ['2000-01-01', '2000-01-01']
    Length: 2, dtype: period[D]

    You can use the string alias for `dtype`

    >>> pd.array(['a', 'b', 'a'], dtype='category')
    ['a', 'b', 'a']
    Categories (2, object): ['a', 'b']

    Or specify the actual dtype

    >>> pd.array(['a', 'b', 'a'],
    ...          dtype=pd.CategoricalDtype(['a', 'b', 'c'], ordered=True))
    ['a', 'b', 'a']
    Categories (3, object): ['a' < 'b' < 'c']

    If pandas does not infer a dedicated extension type a
    :class:`arrays.PandasArray` is returned.

    >>> pd.array([1 + 1j, 3 + 2j])
    <PandasArray>
    [(1+1j), (3+2j)]
    Length: 2, dtype: complex128

    As mentioned in the "Notes" section, new extension types may be added
    in the future (by pandas or 3rd party libraries), causing the return
    value to no longer be a :class:`arrays.PandasArray`. Specify the `dtype`
    as a NumPy dtype if you need to ensure there's no future change in
    behavior.

    >>> pd.array([1, 2], dtype=np.dtype("int32"))
    <PandasArray>
    [1, 2]
    Length: 2, dtype: int32

    `data` must be 1-dimensional. A ValueError is raised when the input
    has the wrong dimensionality.

    >>> pd.array(1)
    Traceback (most recent call last):
      ...
    ValueError: Cannot pass scalar '1' to 'pandas.array'.
    """
    from pandas.core.arrays import (
        BooleanArray,
        DatetimeArray,
        FloatingArray,
        IntegerArray,
        IntervalArray,
        PandasArray,
        StringArray,
        TimedeltaArray,
        period_array,
    )

    if lib.is_scalar(data):
        msg = f"Cannot pass scalar '{data}' to 'pandas.array'."
        raise ValueError(msg)

    if dtype is None and isinstance(
        data, (ABCSeries, ABCIndexClass, ABCExtensionArray)
    ):
        dtype = data.dtype

    data = extract_array(data, extract_numpy=True)

    # this returns None for not-found dtypes.
    if isinstance(dtype, str):
        dtype = registry.find(dtype) or dtype

    if is_extension_array_dtype(dtype):
        cls = cast(ExtensionDtype, dtype).construct_array_type()
        return cls._from_sequence(data, dtype=dtype, copy=copy)

    if dtype is None:
        inferred_dtype = lib.infer_dtype(data, skipna=True)
        if inferred_dtype == "period":
            try:
                return period_array(data, copy=copy)
            except IncompatibleFrequency:
                # We may have a mixture of frequencies.
                # We choose to return an ndarray, rather than raising.
                pass
        elif inferred_dtype == "interval":
            try:
                return IntervalArray(data, copy=copy)
            except ValueError:
                # We may have a mixture of `closed` here.
                # We choose to return an ndarray, rather than raising.
                pass

        elif inferred_dtype.startswith("datetime"):
            # datetime, datetime64
            try:
                return DatetimeArray._from_sequence(data, copy=copy)
            except ValueError:
                # Mixture of timezones, fall back to PandasArray
                pass

        elif inferred_dtype.startswith("timedelta"):
            # timedelta, timedelta64
            return TimedeltaArray._from_sequence(data, copy=copy)

        elif inferred_dtype == "string":
            return StringArray._from_sequence(data, copy=copy)

        elif inferred_dtype == "integer":
            return IntegerArray._from_sequence(data, copy=copy)

        elif inferred_dtype in ("floating", "mixed-integer-float"):
            return FloatingArray._from_sequence(data, copy=copy)

        elif inferred_dtype == "boolean":
            return BooleanArray._from_sequence(data, copy=copy)

    # Pandas overrides NumPy for
    #   1. datetime64[ns]
    #   2. timedelta64[ns]
    # so that a DatetimeArray is returned.
    if is_datetime64_ns_dtype(dtype):
        return DatetimeArray._from_sequence(data, dtype=dtype, copy=copy)
    elif is_timedelta64_ns_dtype(dtype):
        return TimedeltaArray._from_sequence(data, dtype=dtype, copy=copy)

    result = PandasArray._from_sequence(data, dtype=dtype, copy=copy)
    return result


def extract_array(obj: object, extract_numpy: bool = False) -> Union[Any, ArrayLike]:
    """
    Extract the ndarray or ExtensionArray from a Series or Index.

    For all other types, `obj` is just returned as is.

    Parameters
    ----------
    obj : object
        For Series / Index, the underlying ExtensionArray is unboxed.
        For Numpy-backed ExtensionArrays, the ndarray is extracted.

    extract_numpy : bool, default False
        Whether to extract the ndarray from a PandasArray

    Returns
    -------
    arr : object

    Examples
    --------
    >>> extract_array(pd.Series(['a', 'b', 'c'], dtype='category'))
    ['a', 'b', 'c']
    Categories (3, object): ['a', 'b', 'c']

    Other objects like lists, arrays, and DataFrames are just passed through.

    >>> extract_array([1, 2, 3])
    [1, 2, 3]

    For an ndarray-backed Series / Index a PandasArray is returned.

    >>> extract_array(pd.Series([1, 2, 3]))
    <PandasArray>
    [1, 2, 3]
    Length: 3, dtype: int64

    To extract all the way down to the ndarray, pass ``extract_numpy=True``.

    >>> extract_array(pd.Series([1, 2, 3]), extract_numpy=True)
    array([1, 2, 3])
    """
    if isinstance(obj, (ABCIndexClass, ABCSeries)):
        obj = obj.array

    if extract_numpy and isinstance(obj, ABCPandasArray):
        obj = obj.to_numpy()

    return obj


def ensure_wrapped_if_datetimelike(arr):
    """
    Wrap datetime64 and timedelta64 ndarrays in DatetimeArray/TimedeltaArray.
    """
    if isinstance(arr, np.ndarray):
        if arr.dtype.kind == "M":
            from pandas.core.arrays import DatetimeArray

            return DatetimeArray._from_sequence(arr)

        elif arr.dtype.kind == "m":
            from pandas.core.arrays import TimedeltaArray

            return TimedeltaArray._from_sequence(arr)

    return arr


def sanitize_array(
    data,
    index: Optional[Index],
    dtype: Optional[DtypeObj] = None,
    copy: bool = False,
    raise_cast_failure: bool = False,
) -> ArrayLike:
    """
    Sanitize input data to an ndarray or ExtensionArray, copy if specified,
    coerce to the dtype if specified.
    """

    if isinstance(data, ma.MaskedArray):
        mask = ma.getmaskarray(data)
        if mask.any():
            data, fill_value = maybe_upcast(data, copy=True)
            data.soften_mask()  # set hardmask False if it was True
            data[mask] = fill_value
        else:
            data = data.copy()

    # extract ndarray or ExtensionArray, ensure we have no PandasArray
    data = extract_array(data, extract_numpy=True)

    # GH#846
    if isinstance(data, np.ndarray):

        if dtype is not None and is_float_dtype(data.dtype) and is_integer_dtype(dtype):
            # possibility of nan -> garbage
            try:
                subarr = _try_cast(data, dtype, copy, True)
            except ValueError:
                if copy:
                    subarr = data.copy()
                else:
                    subarr = np.array(data, copy=False)
        else:
            # we will try to copy be-definition here
            subarr = _try_cast(data, dtype, copy, raise_cast_failure)

    elif isinstance(data, ABCExtensionArray):
        # it is already ensured above this is not a PandasArray
        subarr = data

        if dtype is not None:
            subarr = subarr.astype(dtype, copy=copy)
        elif copy:
            subarr = subarr.copy()
        return subarr

    elif isinstance(data, (list, tuple, abc.Set, abc.ValuesView)) and len(data) > 0:
        if isinstance(data, set):
            # Raise only for unordered sets, e.g., not for dict_keys
            raise TypeError("Set type is unordered")
        data = list(data)

        if dtype is not None:
            subarr = _try_cast(data, dtype, copy, raise_cast_failure)
        else:
            subarr = maybe_convert_platform(data)

        subarr = maybe_cast_to_datetime(subarr, dtype)

    elif isinstance(data, range):
        # GH#16804
        arr = np.arange(data.start, data.stop, data.step, dtype="int64")
        subarr = _try_cast(arr, dtype, copy, raise_cast_failure)
    elif lib.is_scalar(data) and index is not None and dtype is not None:
        data = maybe_cast_to_datetime(data, dtype)
        if not lib.is_scalar(data):
            data = data[0]
        subarr = construct_1d_arraylike_from_scalar(data, len(index), dtype)
    else:
        subarr = _try_cast(data, dtype, copy, raise_cast_failure)

    # scalar like, GH
    if getattr(subarr, "ndim", 0) == 0:
        if isinstance(data, list):  # pragma: no cover
            subarr = np.array(data, dtype=object)
        elif index is not None:
            value = data

            # figure out the dtype from the value (upcast if necessary)
            if dtype is None:
                dtype, value = infer_dtype_from_scalar(value, pandas_dtype=True)
            else:
                # need to possibly convert the value here
                value = maybe_cast_to_datetime(value, dtype)

            subarr = construct_1d_arraylike_from_scalar(value, len(index), dtype)

        else:
            return subarr.item()

    # the result that we want
    elif subarr.ndim == 1:
        if index is not None:

            # a 1-element ndarray
            if len(subarr) != len(index) and len(subarr) == 1:
                subarr = construct_1d_arraylike_from_scalar(
                    subarr[0], len(index), subarr.dtype
                )

    elif subarr.ndim > 1:
        if isinstance(data, np.ndarray):
            raise ValueError("Data must be 1-dimensional")
        else:
            subarr = com.asarray_tuplesafe(data, dtype=dtype)

    if not (is_extension_array_dtype(subarr.dtype) or is_extension_array_dtype(dtype)):
        # This is to prevent mixed-type Series getting all casted to
        # NumPy string type, e.g. NaN --> '-1#IND'.
        if issubclass(subarr.dtype.type, str):
            # GH#16605
            # If not empty convert the data to dtype
            # GH#19853: If data is a scalar, subarr has already the result
            if not lib.is_scalar(data):
                if not np.all(isna(data)):
                    data = np.array(data, dtype=dtype, copy=False)
                subarr = np.array(data, dtype=object, copy=copy)

        is_object_or_str_dtype = is_object_dtype(dtype) or is_string_dtype(dtype)
        if is_object_dtype(subarr.dtype) and not is_object_or_str_dtype:
            inferred = lib.infer_dtype(subarr, skipna=False)
            if inferred in {"interval", "period"}:
                subarr = array(subarr)

    return subarr


def _try_cast(arr, dtype: Optional[DtypeObj], copy: bool, raise_cast_failure: bool):
    """
    Convert input to numpy ndarray and optionally cast to a given dtype.

    Parameters
    ----------
    arr : ndarray, scalar, list, tuple, iterator (catchall)
        Excludes: ExtensionArray, Series, Index.
    dtype : np.dtype, ExtensionDtype or None
    copy : bool
        If False, don't copy the data if not needed.
    raise_cast_failure : bool
        If True, and if a dtype is specified, raise errors during casting.
        Otherwise an object array is returned.
    """
    # perf shortcut as this is the most common case
    if isinstance(arr, np.ndarray):
        if maybe_castable(arr) and not copy and dtype is None:
            return arr

    if isinstance(dtype, ExtensionDtype) and (dtype.kind != "M" or is_sparse(dtype)):
        # create an extension array from its dtype
        # DatetimeTZ case needs to go through maybe_cast_to_datetime but
        # SparseDtype does not
        array_type = dtype.construct_array_type()._from_sequence
        subarr = array_type(arr, dtype=dtype, copy=copy)
        return subarr

    try:
        # GH#15832: Check if we are requesting a numeric dtype and
        # that we can convert the data to the requested dtype.
        if is_integer_dtype(dtype):
            # this will raise if we have e.g. floats
            maybe_cast_to_integer_array(arr, dtype)
            subarr = arr
        else:
            subarr = maybe_cast_to_datetime(arr, dtype)

        # Take care in creating object arrays (but iterators are not
        # supported):
        if is_object_dtype(dtype) and (
            is_list_like(subarr)
            and not (is_iterator(subarr) or isinstance(subarr, np.ndarray))
        ):
            subarr = construct_1d_object_array_from_listlike(subarr)
        elif not is_extension_array_dtype(subarr):
            subarr = construct_1d_ndarray_preserving_na(subarr, dtype, copy=copy)
    except OutOfBoundsDatetime:
        # in case of out of bound datetime64 -> always raise
        raise
    except (ValueError, TypeError):
        if dtype is not None and raise_cast_failure:
            raise
        else:
            subarr = np.array(arr, dtype=object, copy=copy)
    return subarr


def is_empty_data(data: Any) -> bool:
    """
    Utility to check if a Series is instantiated with empty data,
    which does not contain dtype information.

    Parameters
    ----------
    data : array-like, Iterable, dict, or scalar value
        Contains data stored in Series.

    Returns
    -------
    bool
    """
    is_none = data is None
    is_list_like_without_dtype = is_list_like(data) and not hasattr(data, "dtype")
    is_simple_empty = is_list_like_without_dtype and not data
    return is_none or is_simple_empty


def create_series_with_explicit_dtype(
    data: Any = None,
    index: Optional[Union[ArrayLike, Index]] = None,
    dtype: Optional[Dtype] = None,
    name: Optional[str] = None,
    copy: bool = False,
    fastpath: bool = False,
    dtype_if_empty: Dtype = object,
) -> Series:
    """
    Helper to pass an explicit dtype when instantiating an empty Series.

    This silences a DeprecationWarning described in GitHub-17261.

    Parameters
    ----------
    data : Mirrored from Series.__init__
    index : Mirrored from Series.__init__
    dtype : Mirrored from Series.__init__
    name : Mirrored from Series.__init__
    copy : Mirrored from Series.__init__
    fastpath : Mirrored from Series.__init__
    dtype_if_empty : str, numpy.dtype, or ExtensionDtype
        This dtype will be passed explicitly if an empty Series will
        be instantiated.

    Returns
    -------
    Series
    """
    from pandas.core.series import Series

    if is_empty_data(data) and dtype is None:
        dtype = dtype_if_empty
    return Series(
        data=data, index=index, dtype=dtype, name=name, copy=copy, fastpath=fastpath
    )
