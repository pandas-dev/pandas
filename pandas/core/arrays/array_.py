from typing import Optional, Sequence, Union

import numpy as np

from pandas._libs import lib, tslibs

from pandas.core.dtypes.common import (
    is_datetime64_ns_dtype, is_extension_array_dtype, is_timedelta64_ns_dtype)
from pandas.core.dtypes.dtypes import ExtensionDtype, registry


def array(data,         # type: Sequence[object]
          dtype=None,   # type: Optional[Union[str, np.dtype, ExtensionDtype]]
          copy=True,    # type: bool
          ):
    # type: (...) -> ExtensionArray
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
        ============================== =====================================

        For all other cases, NumPy's usual inference rules will be used.

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

    Or use the dedicated constructor for the array you're expecting, and
    wrap that in a PandasArray

    >>> pd.array(np.array(['a', 'b'], dtype='<U1'))
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
    ['01:00:00', '02:00:00']
    Length: 2, dtype: timedelta64[ns]

    Examples
    --------
    If a dtype is not specified, `data` is passed through to
    :meth:`numpy.array`, and a :class:`arrays.PandasArray` is returned.

    >>> pd.array([1, 2])
    <PandasArray>
    [1, 2]
    Length: 2, dtype: int64

    Or the NumPy dtype can be specified

    >>> pd.array([1, 2], dtype=np.dtype("int32"))
    <PandasArray>
    [1, 2]
    Length: 2, dtype: int32

    You can use the string alias for `dtype`

    >>> pd.array(['a', 'b', 'a'], dtype='category')
    [a, b, a]
    Categories (2, object): [a, b]

    Or specify the actual dtype

    >>> pd.array(['a', 'b', 'a'],
    ...          dtype=pd.CategoricalDtype(['a', 'b', 'c'], ordered=True))
    [a, b, a]
    Categories (3, object): [a < b < c]

    Because omitting the `dtype` passes the data through to NumPy,
    a mixture of valid integers and NA will return a floating-point
    NumPy array.

    >>> pd.array([1, 2, np.nan])
    <PandasArray>
    [1.0,  2.0, nan]
    Length: 3, dtype: float64

    To use pandas' nullable :class:`pandas.arrays.IntegerArray`, specify
    the dtype:

    >>> pd.array([1, 2, np.nan], dtype='Int64')
    <IntegerArray>
    [1, 2, NaN]
    Length: 3, dtype: Int64

    Pandas will infer an ExtensionArray for some types of data:

    >>> pd.array([pd.Period('2000', freq="D"), pd.Period("2000", freq="D")])
    <PeriodArray>
    ['2000-01-01', '2000-01-01']
    Length: 2, dtype: period[D]

    `data` must be 1-dimensional. A ValueError is raised when the input
    has the wrong dimensionality.

    >>> pd.array(1)
    Traceback (most recent call last):
      ...
    ValueError: Cannot pass scalar '1' to 'pandas.array'.
    """
    from pandas.core.arrays import (
        period_array, ExtensionArray, IntervalArray, PandasArray,
        DatetimeArray,
        TimedeltaArray,
    )
    from pandas.core.internals.arrays import extract_array

    if lib.is_scalar(data):
        msg = (
            "Cannot pass scalar '{}' to 'pandas.array'."
        )
        raise ValueError(msg.format(data))

    data = extract_array(data, extract_numpy=True)

    if dtype is None and isinstance(data, ExtensionArray):
        dtype = data.dtype

    # this returns None for not-found dtypes.
    if isinstance(dtype, str):
        dtype = registry.find(dtype) or dtype

    if is_extension_array_dtype(dtype):
        cls = dtype.construct_array_type()
        return cls._from_sequence(data, dtype=dtype, copy=copy)

    if dtype is None:
        inferred_dtype = lib.infer_dtype(data, skipna=False)
        if inferred_dtype == 'period':
            try:
                return period_array(data, copy=copy)
            except tslibs.IncompatibleFrequency:
                # We may have a mixture of frequencies.
                # We choose to return an ndarray, rather than raising.
                pass
        elif inferred_dtype == 'interval':
            try:
                return IntervalArray(data, copy=copy)
            except ValueError:
                # We may have a mixture of `closed` here.
                # We choose to return an ndarray, rather than raising.
                pass

        elif inferred_dtype.startswith('datetime'):
            # datetime, datetime64
            try:
                return DatetimeArray._from_sequence(data, copy=copy)
            except ValueError:
                # Mixture of timezones, fall back to PandasArray
                pass

        elif inferred_dtype.startswith('timedelta'):
            # timedelta, timedelta64
            return TimedeltaArray._from_sequence(data, copy=copy)

        # TODO(BooleanArray): handle this type

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
