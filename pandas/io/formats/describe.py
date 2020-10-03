"""Module responsible for execution of NDFrame.describe() method.

Method NDFrame.describe() delegates actual execution to function describe_ndframe().

Strategy pattern is utilized.
 - The appropriate strategy is selected based on the series datatype.
 - The strategy is responsible for running proper description.
"""

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, List, Optional, Sequence, Type, Union, cast
import warnings

import numpy as np

from pandas._libs.tslibs import Timestamp
from pandas._typing import Dtype, FrameOrSeries, FrameOrSeriesUnion, Label
from pandas.util._validators import validate_percentile

from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
    is_timedelta64_dtype,
)

import pandas as pd

from pandas.io.formats.format import format_percentiles

if TYPE_CHECKING:
    from pandas import DataFrame, Series


def describe_ndframe(
    *,
    data: FrameOrSeries,
    include: Optional[Union[str, Sequence[str]]],
    exclude: Optional[Union[str, Sequence[str]]],
    datetime_is_numeric: bool,
    percentiles: Optional[Sequence[float]],
) -> FrameOrSeries:
    """Describe series or dataframe.

    Called from pandas.core.generic.NDFrame.describe()

    Parameters
    ----------
    data : FrameOrSeries
        Either dataframe or series.
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
    FrameOrSeries
        Dataframe or series description.
    """
    describer: "NDFrameDescriber"
    if data.ndim == 1:
        series = cast("Series", data)
        describer = SeriesDescriber(
            data=series,
            datetime_is_numeric=datetime_is_numeric,
        )
    else:
        dataframe = cast("DataFrame", data)
        describer = DataFrameDescriber(
            data=dataframe,
            include=include,
            exclude=exclude,
            datetime_is_numeric=datetime_is_numeric,
        )
    result = describer.describe(percentiles)
    return cast(FrameOrSeries, result)


class StrategyCreatorMixin:
    """Mixin for creating instance of appropriate strategy for describing series."""

    datetime_is_numeric: bool

    def create_strategy(
        self,
        series: "Series",
        percentiles: Optional[Sequence[float]],
    ) -> "StrategyAbstract":
        """Create strategy instance for description."""
        klass = self._select_strategy(series)
        return klass(series, percentiles)

    def _select_strategy(self, series: "Series") -> Type["StrategyAbstract"]:
        """Select strategy for description."""
        strategy: Type[StrategyAbstract] = CategoricalStrategy
        if is_bool_dtype(series.dtype):
            strategy = CategoricalStrategy
        elif is_numeric_dtype(series):
            strategy = NumericStrategy
        elif is_datetime64_any_dtype(series.dtype) and self.datetime_is_numeric:
            strategy = TimestampStrategy
        elif is_timedelta64_dtype(series.dtype):
            strategy = NumericStrategy

        if strategy == CategoricalStrategy and is_datetime64_any_dtype(series.dtype):
            strategy = TimestampAsCategoricalStrategy
            warnings.warn(
                "Treating datetime data as categorical rather than numeric in "
                "`.describe` is deprecated and will be removed in a future "
                "version of pandas. Specify `datetime_is_numeric=True` to "
                "silence this warning and adopt the future behavior now.",
                FutureWarning,
                stacklevel=6,
            )
        return strategy


class NDFrameDescriber(ABC):
    """Abstract class for describing dataframe or series."""

    @abstractmethod
    def describe(self, percentiles: Optional[Sequence[float]]) -> FrameOrSeriesUnion:
        """Do describe either series or dataframe.

        Parameters
        ----------
        percentiles : list-like of numbers, optional
            The percentiles to include in the output. All should fall between 0 and 1.
            The default is ``[.25, .5, .75]``, which returns the 25th, 50th, and
            75th percentiles.
        """


class SeriesDescriber(NDFrameDescriber, StrategyCreatorMixin):
    """Class responsible for creating series description.

    Parameters
    ----------
    data : Series
        Series to be described.
    datetime_is_numeric : bool, default False
        Whether to treat datetime dtypes as numeric.
    """

    def __init__(
        self,
        *,
        data: "Series",
        datetime_is_numeric: bool,
    ):
        self.data = data
        self.datetime_is_numeric = datetime_is_numeric

    def describe(self, percentiles: Optional[Sequence[float]]) -> "Series":
        """Do describe series."""
        strategy = self.create_strategy(self.data, percentiles)
        result = strategy.describe()
        return result


class DataFrameDescriber(NDFrameDescriber, StrategyCreatorMixin):
    """Class responsible for creating dataframe description.

    Parameters
    ----------
    data : DataFrame
        Dataframe to be described.
    include : 'all', list-like of dtypes or None (default), optional
        A white list of data types to include in the result.
    exclude : list-like of dtypes or None (default), optional,
        A black list of data types to omit from the result.
    datetime_is_numeric : bool, default False
        Whether to treat datetime dtypes as numeric.
    """

    def __init__(
        self,
        *,
        data: "DataFrame",
        include: Optional[Union[str, Sequence[str]]],
        exclude: Optional[Union[str, Sequence[str]]],
        datetime_is_numeric: bool,
    ):
        self.include = include
        self.exclude = exclude
        self.datetime_is_numeric = datetime_is_numeric
        self.data: "DataFrame" = self._initialize_data(data)

    def _initialize_data(self, data) -> "DataFrame":
        _validate_dframe_size(data)

        if self.include is None and self.exclude is None:
            # when some numerics are found, keep only numerics
            include = [np.number]
            if self.datetime_is_numeric:
                include.append("datetime")
            numeric_only = data.select_dtypes(include=include)
            if len(numeric_only.columns) == 0:
                return data
            else:
                return numeric_only

        if self.include == "all":
            if self.exclude is not None:
                msg = "exclude must be None when include is 'all'"
                raise ValueError(msg)
            return data

        return data.select_dtypes(include=self.include, exclude=self.exclude)

    def describe(self, percentiles: Optional[Sequence[float]]) -> "DataFrame":
        """Do describe dataframe."""
        ldesc: List["Series"] = []
        for _, series in self.data.items():
            strategy = self.create_strategy(series, percentiles)
            ldesc.append(strategy.describe())

        df = pd.concat(
            self._reindex_columns(ldesc),
            axis=1,
            sort=False,
        )
        df.columns = self.data.columns.copy()
        return cast("DataFrame", df)

    def _reindex_columns(self, column_data) -> List["Series"]:
        """Set a convenient order for rows."""
        names: List[Label] = []
        ldesc_indexes = sorted((x.index for x in column_data), key=len)
        for idxnames in ldesc_indexes:
            for name in idxnames:
                if name not in names:
                    names.append(name)
        return [x.reindex(names, copy=False) for x in column_data]


class StrategyAbstract(ABC):
    """Abstract strategy for describing series."""

    def __init__(
        self,
        data: "Series",
        percentiles: Optional[Sequence[float]],
    ):
        self.data = data
        self.percentiles = self._initialize_percentiles(percentiles)

    def describe(self) -> "Series":
        """Describe series."""
        return pd.Series(
            self.array,
            index=self.names,
            name=self.data.name,
            dtype=self.dtype,
        )

    @property
    @abstractmethod
    def array(self) -> List[object]:
        """Series data."""

    @property
    @abstractmethod
    def names(self) -> List[str]:
        """Series index names."""

    @property
    @abstractmethod
    def dtype(self) -> Optional[Dtype]:
        """Series dtype."""

    @property
    def formatted_percentiles(self) -> List[str]:
        """Percentiles formatted as strings, rounded."""
        return format_percentiles(self.percentiles)

    @staticmethod
    def _initialize_percentiles(
        percentiles: Optional[Sequence[float]],
    ) -> Sequence[float]:
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


class CategoricalStrategy(StrategyAbstract):
    """Strategy for series with categorical values."""

    def __init__(self, data, percentiles):
        self.data = data
        super().__init__(data, percentiles)
        self.objcounts = self.data.value_counts()

    @property
    def array(self) -> List[object]:
        top, freq = self._get_top_and_freq()
        return [
            self.count,
            self.count_unique,
            top,
            freq,
        ]

    @property
    def names(self) -> List[str]:
        return ["count", "unique", "top", "freq"]

    @property
    def dtype(self) -> Optional[Dtype]:
        if self.count_unique == 0:
            return "object"
        return None

    @property
    def count(self) -> "Series":
        return self.data.count()

    @property
    def count_unique(self) -> int:
        return len(self.objcounts[self.objcounts != 0])

    def _get_top_and_freq(self):
        if self.count_unique > 0:
            return self.objcounts.index[0], self.objcounts.iloc[0]
        return np.nan, np.nan


class TimestampAsCategoricalStrategy(CategoricalStrategy):
    """Strategy for series with timestamp values treated as categorical values."""

    @property
    def array(self) -> List[object]:
        result = [self.count, self.count_unique]
        if self.count_unique > 0:
            top, freq = self.objcounts.index[0], self.objcounts.iloc[0]
            tz = self.data.dt.tz
            asint = self.data.dropna().values.view("i8")
            top = Timestamp(top)
            if top.tzinfo is not None and tz is not None:
                # Don't tz_localize(None) if key is already tz-aware
                top = top.tz_convert(tz)
            else:
                top = top.tz_localize(tz)

            result += [
                top,
                freq,
                Timestamp(asint.min(), tz=tz),
                Timestamp(asint.max(), tz=tz),
            ]

        # If the DataFrame is empty, set 'top' and 'freq' to None
        # to maintain output shape consistency
        else:
            result += [np.nan, np.nan]
        return result

    @property
    def names(self) -> List[str]:
        names = ["count", "unique"]
        if self.count_unique > 0:
            names += ["top", "freq", "first", "last"]
        return names


class NumericStrategy(StrategyAbstract):
    """Strategy for series with numeric values."""

    @property
    def array(self) -> List[object]:
        return [
            self.data.count(),
            self.data.mean(),
            self.data.std(),
            self.data.min(),
            *self.data.quantile(self.percentiles).tolist(),
            self.data.max(),
        ]

    @property
    def names(self) -> List[str]:
        return [
            "count",
            "mean",
            "std",
            "min",
            *self.formatted_percentiles,
            "max",
        ]

    @property
    def dtype(self) -> Optional[Dtype]:
        return None


class TimestampStrategy(StrategyAbstract):
    """Strategy for series with timestamp values."""

    @property
    def array(self) -> List[object]:
        return [
            self.data.count(),
            self.data.mean(),
            self.data.min(),
            *self.data.quantile(self.percentiles).tolist(),
            self.data.max(),
        ]

    @property
    def names(self) -> List[str]:
        return [
            "count",
            "mean",
            "min",
            *self.formatted_percentiles,
            "max",
        ]

    @property
    def dtype(self) -> Optional[Dtype]:
        return None


def _validate_dframe_size(df: FrameOrSeriesUnion) -> None:
    """Validate correct size of dataframe."""
    if df.ndim == 2 and df.columns.size == 0:
        raise ValueError("Cannot describe a DataFrame without columns")
