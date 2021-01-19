"""
Module responsible for execution of NDFrame.describe() method.

Method NDFrame.describe() delegates actual execution to function describe_ndframe().
"""
from __future__ import annotations

from typing import TYPE_CHECKING, List, Optional, Sequence, Union, cast
import warnings

import numpy as np

from pandas._libs.tslibs import Timestamp
from pandas._typing import FrameOrSeries, Hashable
from pandas.util._validators import validate_percentile

from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
    is_timedelta64_dtype,
)

from pandas.core.reshape.concat import concat

from pandas.io.formats.format import format_percentiles

if TYPE_CHECKING:
    from pandas import DataFrame, Series


def describe_ndframe(
    *,
    obj: FrameOrSeries,
    include: Optional[Union[str, Sequence[str]]],
    exclude: Optional[Union[str, Sequence[str]]],
    datetime_is_numeric: bool,
    percentiles: Optional[Sequence[float]],
) -> FrameOrSeries:
    """Describe series or dataframe.

    Called from pandas.core.generic.NDFrame.describe()

    Parameters
    ----------
    obj: DataFrame or Series
        Either dataframe or series to be described.
    include : 'all', list-like of dtypes or None (default), optional
        A white list of data types to include in the result. Ignored for ``Series``.
    exclude : list-like of dtypes or None (default), optional,
        A black list of data types to omit from the result. Ignored for ``Series``.
    datetime_is_numeric : bool, default False
        Whether to treat datetime dtypes as numeric.
    percentiles : list-like of numbers, optional
        The percentiles to include in the output. All should fall between 0 and 1.
        The default is ``[.25, .5, .75]``, which returns the 25th, 50th, and
        75th percentiles.

    Returns
    -------
    Dataframe or series description.
    """
    percentiles = refine_percentiles(percentiles)

    if obj.ndim == 1:
        result_series = describe_series(
            cast("Series", obj),
            percentiles,
            datetime_is_numeric,
        )
        return cast(FrameOrSeries, result_series)

    frame = cast("DataFrame", obj)

    if frame.ndim == 2 and frame.columns.size == 0:
        raise ValueError("Cannot describe a DataFrame without columns")

    result_frame = describe_frame(
        frame=frame,
        include=include,
        exclude=exclude,
        percentiles=percentiles,
        datetime_is_numeric=datetime_is_numeric,
    )
    return cast(FrameOrSeries, result_frame)


def describe_series(
    series: "Series",
    percentiles: Sequence[float],
    datetime_is_numeric: bool,
) -> Series:
    """Describe series.

    The reason for the delegation to ``describe_1d`` only:
    to allow for a proper stacklevel of the FutureWarning.

    Parameters
    ----------
    series : Series
        Series to be described.
    percentiles : list-like of numbers
        The percentiles to include in the output.
    datetime_is_numeric : bool, default False
        Whether to treat datetime dtypes as numeric.

    Returns
    -------
    Series
    """
    return describe_1d(
        series,
        percentiles,
        datetime_is_numeric,
        is_series=True,
    )


def describe_frame(
    frame: "DataFrame",
    include: Optional[Union[str, Sequence[str]]],
    exclude: Optional[Union[str, Sequence[str]]],
    percentiles: Sequence[float],
    datetime_is_numeric: bool,
) -> DataFrame:
    """Describe DataFrame.

    Parameters
    ----------
    frame : DataFrame
        DataFrame to be described.
    include : 'all', list-like of dtypes or None (default), optional
        A white list of data types to include in the result.
    exclude : list-like of dtypes or None (default), optional,
        A black list of data types to omit from the result.
    percentiles : list-like of numbers
        The percentiles to include in the output.
    datetime_is_numeric : bool, default False
        Whether to treat datetime dtypes as numeric.

    Returns
    -------
    DataFrame
    """
    data = select_columns(
        frame=frame,
        include=include,
        exclude=exclude,
        datetime_is_numeric=datetime_is_numeric,
    )

    ldesc = [
        describe_1d(s, percentiles, datetime_is_numeric, is_series=False)
        for _, s in data.items()
    ]

    col_names = reorder_columns(ldesc)
    d = concat(
        [x.reindex(col_names, copy=False) for x in ldesc],
        axis=1,
        sort=False,
    )
    d.columns = data.columns.copy()
    return d


def reorder_columns(ldesc: Sequence["Series"]) -> List[Hashable]:
    """Set a convenient order for rows for display."""
    names: List[Hashable] = []
    ldesc_indexes = sorted((x.index for x in ldesc), key=len)
    for idxnames in ldesc_indexes:
        for name in idxnames:
            if name not in names:
                names.append(name)
    return names


def select_columns(
    frame: "DataFrame",
    include: Optional[Union[str, Sequence[str]]],
    exclude: Optional[Union[str, Sequence[str]]],
    datetime_is_numeric: bool,
) -> DataFrame:
    """Select columns to be described."""
    if (include is None) and (exclude is None):
        # when some numerics are found, keep only numerics
        default_include = [np.number]
        if datetime_is_numeric:
            default_include.append("datetime")
        data = frame.select_dtypes(include=default_include)
        if len(data.columns) == 0:
            data = frame
    elif include == "all":
        if exclude is not None:
            msg = "exclude must be None when include is 'all'"
            raise ValueError(msg)
        data = frame
    else:
        data = frame.select_dtypes(include=include, exclude=exclude)

    return data


def describe_numeric_1d(series: "Series", percentiles: Sequence[float]) -> Series:
    """Describe series containing numerical data.

    Parameters
    ----------
    series : Series
        Series to be described.
    percentiles : list-like of numbers
        The percentiles to include in the output.
    """
    from pandas import Series

    formatted_percentiles = format_percentiles(percentiles)

    stat_index = ["count", "mean", "std", "min"] + formatted_percentiles + ["max"]
    d = (
        [series.count(), series.mean(), series.std(), series.min()]
        + series.quantile(percentiles).tolist()
        + [series.max()]
    )
    return Series(d, index=stat_index, name=series.name)


def describe_categorical_1d(data: "Series", is_series: bool) -> Series:
    """Describe series containing categorical data.

    Parameters
    ----------
    data : Series
        Series to be described.
    is_series : bool
        True if the original object is a Series.
        False if the one column of the DataFrame is described.
    """
    names = ["count", "unique"]
    objcounts = data.value_counts()
    count_unique = len(objcounts[objcounts != 0])
    result = [data.count(), count_unique]
    dtype = None
    if result[1] > 0:
        top, freq = objcounts.index[0], objcounts.iloc[0]
        if is_datetime64_any_dtype(data.dtype):
            if is_series:
                stacklevel = 6
            else:
                stacklevel = 7
            warnings.warn(
                "Treating datetime data as categorical rather than numeric in "
                "`.describe` is deprecated and will be removed in a future "
                "version of pandas. Specify `datetime_is_numeric=True` to "
                "silence this warning and adopt the future behavior now.",
                FutureWarning,
                stacklevel=stacklevel,
            )
            tz = data.dt.tz
            asint = data.dropna().values.view("i8")
            top = Timestamp(top)
            if top.tzinfo is not None and tz is not None:
                # Don't tz_localize(None) if key is already tz-aware
                top = top.tz_convert(tz)
            else:
                top = top.tz_localize(tz)
            names += ["top", "freq", "first", "last"]
            result += [
                top,
                freq,
                Timestamp(asint.min(), tz=tz),
                Timestamp(asint.max(), tz=tz),
            ]
        else:
            names += ["top", "freq"]
            result += [top, freq]

    # If the DataFrame is empty, set 'top' and 'freq' to None
    # to maintain output shape consistency
    else:
        names += ["top", "freq"]
        result += [np.nan, np.nan]
        dtype = "object"

    from pandas import Series

    return Series(result, index=names, name=data.name, dtype=dtype)


def describe_timestamp_1d(data: "Series", percentiles: Sequence[float]) -> Series:
    """Describe series containing datetime64 dtype.

    Parameters
    ----------
    data : Series
        Series to be described.
    percentiles : list-like of numbers
        The percentiles to include in the output.
    """
    # GH-30164
    from pandas import Series

    formatted_percentiles = format_percentiles(percentiles)

    stat_index = ["count", "mean", "min"] + formatted_percentiles + ["max"]
    d = (
        [data.count(), data.mean(), data.min()]
        + data.quantile(percentiles).tolist()
        + [data.max()]
    )
    return Series(d, index=stat_index, name=data.name)


def describe_1d(
    data: "Series",
    percentiles: Sequence[float],
    datetime_is_numeric: bool,
    *,
    is_series: bool,
) -> Series:
    """Describe series.

    Parameters
    ----------
    data : Series
        Series to be described.
    percentiles : list-like of numbers
        The percentiles to include in the output.
    datetime_is_numeric : bool, default False
        Whether to treat datetime dtypes as numeric.
    is_series : bool
        True if the original object is a Series.
        False if the one column of the DataFrame is described.

    Returns
    -------
    Series
    """
    if is_bool_dtype(data.dtype):
        return describe_categorical_1d(data, is_series)
    elif is_numeric_dtype(data):
        return describe_numeric_1d(data, percentiles)
    elif is_datetime64_any_dtype(data.dtype) and datetime_is_numeric:
        return describe_timestamp_1d(data, percentiles)
    elif is_timedelta64_dtype(data.dtype):
        return describe_numeric_1d(data, percentiles)
    else:
        return describe_categorical_1d(data, is_series)


def refine_percentiles(percentiles: Optional[Sequence[float]]) -> Sequence[float]:
    """Ensure that percentiles are unique and sorted.

    Parameters
    ----------
    percentiles : list-like of numbers, optional
        The percentiles to include in the output.
    """
    if percentiles is None:
        return np.array([0.25, 0.5, 0.75])

    # explicit conversion of `percentiles` to list
    percentiles = list(percentiles)

    # get them all to be in [0, 1]
    validate_percentile(percentiles)

    # median should always be included
    if 0.5 not in percentiles:
        percentiles.append(0.5)

    percentiles = np.asarray(percentiles)

    # sort and check for duplicates
    unique_pcts = np.unique(percentiles)
    assert percentiles is not None
    if len(unique_pcts) < len(percentiles):
        raise ValueError("percentiles cannot contain duplicates")

    return unique_pcts
