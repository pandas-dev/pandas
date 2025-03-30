"""
This module provides Grouper objects that encapsulate the
"factorization" process - conversion of value we are grouping by
to integer codes (one per group).
"""

from __future__ import annotations

import datetime
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal, cast

import numpy as np
import pandas as pd

from xarray.coding.cftime_offsets import BaseCFTimeOffset, _new_to_legacy_freq
from xarray.core import duck_array_ops
from xarray.core.coordinates import Coordinates
from xarray.core.dataarray import DataArray
from xarray.core.groupby import T_Group, _DummyGroup
from xarray.core.indexes import safe_cast_to_index
from xarray.core.resample_cftime import CFTimeGrouper
from xarray.core.types import (
    Bins,
    DatetimeLike,
    GroupIndices,
    ResampleCompatible,
    SideOptions,
)
from xarray.core.variable import Variable

__all__ = [
    "EncodedGroups",
    "Grouper",
    "Resampler",
    "UniqueGrouper",
    "BinGrouper",
    "TimeResampler",
]

RESAMPLE_DIM = "__resample_dim__"


def _coordinates_from_variable(variable: Variable) -> Coordinates:
    from xarray.core.indexes import create_default_index_implicit

    (name,) = variable.dims
    new_index, index_vars = create_default_index_implicit(variable)
    indexes = {k: new_index for k in index_vars}
    new_vars = new_index.create_variables()
    new_vars[name].attrs = variable.attrs
    return Coordinates(new_vars, indexes)


@dataclass(init=False)
class EncodedGroups:
    """
    Dataclass for storing intermediate values for GroupBy operation.
    Returned by the ``factorize`` method on Grouper objects.

    Attributes
    ----------
    codes : DataArray
        Same shape as the DataArray to group by. Values consist of a unique integer code for each group.
    full_index : pd.Index
        Pandas Index for the group coordinate containing unique group labels.
        This can differ from ``unique_coord`` in the case of resampling and binning,
        where certain groups in the output need not be present in the input.
    group_indices : tuple of int or slice or list of int, optional
        List of indices of array elements belonging to each group. Inferred if not provided.
    unique_coord : Variable, optional
        Unique group values present in dataset. Inferred if not provided
    """

    codes: DataArray
    full_index: pd.Index
    group_indices: GroupIndices
    unique_coord: Variable | _DummyGroup
    coords: Coordinates

    def __init__(
        self,
        codes: DataArray,
        full_index: pd.Index,
        group_indices: GroupIndices | None = None,
        unique_coord: Variable | _DummyGroup | None = None,
        coords: Coordinates | None = None,
    ):
        from xarray.core.groupby import _codes_to_group_indices

        assert isinstance(codes, DataArray)
        if codes.name is None:
            raise ValueError("Please set a name on the array you are grouping by.")
        self.codes = codes
        assert isinstance(full_index, pd.Index)
        self.full_index = full_index

        if group_indices is None:
            self.group_indices = tuple(
                g
                for g in _codes_to_group_indices(codes.data.ravel(), len(full_index))
                if g
            )
        else:
            self.group_indices = group_indices

        if unique_coord is None:
            unique_values = full_index[np.unique(codes)]
            self.unique_coord = Variable(
                dims=codes.name, data=unique_values, attrs=codes.attrs
            )
        else:
            self.unique_coord = unique_coord

        if coords is None:
            assert not isinstance(self.unique_coord, _DummyGroup)
            self.coords = _coordinates_from_variable(self.unique_coord)
        else:
            self.coords = coords


class Grouper(ABC):
    """Abstract base class for Grouper objects that allow specializing GroupBy instructions."""

    @abstractmethod
    def factorize(self, group: T_Group) -> EncodedGroups:
        """
        Creates intermediates necessary for GroupBy.

        Parameters
        ----------
        group : DataArray
            DataArray we are grouping by.

        Returns
        -------
        EncodedGroups
        """
        pass


class Resampler(Grouper):
    """
    Abstract base class for Grouper objects that allow specializing resampling-type GroupBy instructions.

    Currently only used for TimeResampler, but could be used for SpaceResampler in the future.
    """

    pass


@dataclass
class UniqueGrouper(Grouper):
    """Grouper object for grouping by a categorical variable."""

    _group_as_index: pd.Index | None = field(default=None, repr=False)

    @property
    def group_as_index(self) -> pd.Index:
        """Caches the group DataArray as a pandas Index."""
        if self._group_as_index is None:
            if self.group.ndim == 1:
                self._group_as_index = self.group.to_index()
            else:
                self._group_as_index = pd.Index(np.array(self.group).ravel())
        return self._group_as_index

    def factorize(self, group: T_Group) -> EncodedGroups:
        self.group = group

        index = self.group_as_index
        is_unique_and_monotonic = isinstance(self.group, _DummyGroup) or (
            index.is_unique
            and (index.is_monotonic_increasing or index.is_monotonic_decreasing)
        )
        is_dimension = self.group.dims == (self.group.name,)
        can_squeeze = is_dimension and is_unique_and_monotonic

        if can_squeeze:
            return self._factorize_dummy()
        else:
            return self._factorize_unique()

    def _factorize_unique(self) -> EncodedGroups:
        # look through group to find the unique values
        sort = not isinstance(self.group_as_index, pd.MultiIndex)
        unique_values, codes_ = unique_value_groups(self.group_as_index, sort=sort)
        if (codes_ == -1).all():
            raise ValueError(
                "Failed to group data. Are you grouping by a variable that is all NaN?"
            )
        codes = self.group.copy(data=codes_.reshape(self.group.shape), deep=False)
        unique_coord = Variable(
            dims=codes.name, data=unique_values, attrs=self.group.attrs
        )
        full_index = pd.Index(unique_values)

        return EncodedGroups(
            codes=codes,
            full_index=full_index,
            unique_coord=unique_coord,
            coords=_coordinates_from_variable(unique_coord),
        )

    def _factorize_dummy(self) -> EncodedGroups:
        size = self.group.size
        # no need to factorize
        # use slices to do views instead of fancy indexing
        # equivalent to: group_indices = group_indices.reshape(-1, 1)
        group_indices: GroupIndices = tuple(slice(i, i + 1) for i in range(size))
        size_range = np.arange(size)
        full_index: pd.Index
        unique_coord: _DummyGroup | Variable
        if isinstance(self.group, _DummyGroup):
            codes = self.group.to_dataarray().copy(data=size_range)
            unique_coord = self.group
            full_index = pd.RangeIndex(self.group.size)
            coords = Coordinates()
        else:
            codes = self.group.copy(data=size_range, deep=False)
            unique_coord = self.group.variable.to_base_variable()
            full_index = self.group_as_index
            if isinstance(full_index, pd.MultiIndex):
                coords = Coordinates.from_pandas_multiindex(
                    full_index, dim=self.group.name
                )
            else:
                if TYPE_CHECKING:
                    assert isinstance(unique_coord, Variable)
                coords = _coordinates_from_variable(unique_coord)

        return EncodedGroups(
            codes=codes,
            group_indices=group_indices,
            full_index=full_index,
            unique_coord=unique_coord,
            coords=coords,
        )


@dataclass
class BinGrouper(Grouper):
    """
    Grouper object for binning numeric data.

    Attributes
    ----------
    bins : int, sequence of scalars, or IntervalIndex
        The criteria to bin by.

        * int : Defines the number of equal-width bins in the range of `x`. The
          range of `x` is extended by .1% on each side to include the minimum
          and maximum values of `x`.
        * sequence of scalars : Defines the bin edges allowing for non-uniform
          width. No extension of the range of `x` is done.
        * IntervalIndex : Defines the exact bins to be used. Note that
          IntervalIndex for `bins` must be non-overlapping.

    right : bool, default True
        Indicates whether `bins` includes the rightmost edge or not. If
        ``right == True`` (the default), then the `bins` ``[1, 2, 3, 4]``
        indicate (1,2], (2,3], (3,4]. This argument is ignored when
        `bins` is an IntervalIndex.
    labels : array or False, default None
        Specifies the labels for the returned bins. Must be the same length as
        the resulting bins. If False, returns only integer indicators of the
        bins. This affects the type of the output container (see below).
        This argument is ignored when `bins` is an IntervalIndex. If True,
        raises an error. When `ordered=False`, labels must be provided.
    retbins : bool, default False
        Whether to return the bins or not. Useful when bins is provided
        as a scalar.
    precision : int, default 3
        The precision at which to store and display the bins labels.
    include_lowest : bool, default False
        Whether the first interval should be left-inclusive or not.
    duplicates : {"raise", "drop"}, default: "raise"
        If bin edges are not unique, raise ValueError or drop non-uniques.
    """

    bins: Bins
    # The rest are copied from pandas
    right: bool = True
    labels: Any = None
    precision: int = 3
    include_lowest: bool = False
    duplicates: Literal["raise", "drop"] = "raise"

    def __post_init__(self) -> None:
        if duck_array_ops.isnull(self.bins).all():
            raise ValueError("All bin edges are NaN.")

    def factorize(self, group: T_Group) -> EncodedGroups:
        from xarray.core.dataarray import DataArray

        data = np.asarray(group.data)  # Cast _DummyGroup data to array

        binned, self.bins = pd.cut(  # type: ignore [call-overload]
            data.ravel(),
            bins=self.bins,
            right=self.right,
            labels=self.labels,
            precision=self.precision,
            include_lowest=self.include_lowest,
            duplicates=self.duplicates,
            retbins=True,
        )

        binned_codes = binned.codes
        if (binned_codes == -1).all():
            raise ValueError(
                f"None of the data falls within bins with edges {self.bins!r}"
            )

        new_dim_name = f"{group.name}_bins"

        full_index = binned.categories
        uniques = np.sort(pd.unique(binned_codes))
        unique_values = full_index[uniques[uniques != -1]]

        codes = DataArray(
            binned_codes.reshape(group.shape),
            getattr(group, "coords", None),
            name=new_dim_name,
        )
        unique_coord = Variable(
            dims=new_dim_name, data=unique_values, attrs=group.attrs
        )
        return EncodedGroups(
            codes=codes,
            full_index=full_index,
            unique_coord=unique_coord,
            coords=_coordinates_from_variable(unique_coord),
        )


@dataclass(repr=False)
class TimeResampler(Resampler):
    """
    Grouper object specialized to resampling the time coordinate.

    Attributes
    ----------
    freq : str, datetime.timedelta, pandas.Timestamp, or pandas.DateOffset
        Frequency to resample to. See `Pandas frequency
        aliases <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`_
        for a list of possible values.
    closed : {"left", "right"}, optional
        Side of each interval to treat as closed.
    label : {"left", "right"}, optional
        Side of each interval to use for labeling.
    origin : {'epoch', 'start', 'start_day', 'end', 'end_day'}, pandas.Timestamp, datetime.datetime, numpy.datetime64, or cftime.datetime, default 'start_day'
        The datetime on which to adjust the grouping. The timezone of origin
        must match the timezone of the index.

        If a datetime is not used, these values are also supported:
        - 'epoch': `origin` is 1970-01-01
        - 'start': `origin` is the first value of the timeseries
        - 'start_day': `origin` is the first day at midnight of the timeseries
        - 'end': `origin` is the last value of the timeseries
        - 'end_day': `origin` is the ceiling midnight of the last day
    offset : pd.Timedelta, datetime.timedelta, or str, default is None
        An offset timedelta added to the origin.
    """

    freq: ResampleCompatible
    closed: SideOptions | None = field(default=None)
    label: SideOptions | None = field(default=None)
    origin: str | DatetimeLike = field(default="start_day")
    offset: pd.Timedelta | datetime.timedelta | str | None = field(default=None)

    index_grouper: CFTimeGrouper | pd.Grouper = field(init=False, repr=False)
    group_as_index: pd.Index = field(init=False, repr=False)

    def _init_properties(self, group: T_Group) -> None:
        from xarray import CFTimeIndex

        group_as_index = safe_cast_to_index(group)
        offset = self.offset

        if not group_as_index.is_monotonic_increasing:
            # TODO: sort instead of raising an error
            raise ValueError("Index must be monotonic for resampling")

        if isinstance(group_as_index, CFTimeIndex):
            from xarray.core.resample_cftime import CFTimeGrouper

            self.index_grouper = CFTimeGrouper(
                freq=self.freq,
                closed=self.closed,
                label=self.label,
                origin=self.origin,
                offset=offset,
            )
        else:
            if isinstance(self.freq, BaseCFTimeOffset):
                raise ValueError(
                    "'BaseCFTimeOffset' resample frequencies are only supported "
                    "when resampling a 'CFTimeIndex'"
                )

            self.index_grouper = pd.Grouper(
                # TODO remove once requiring pandas >= 2.2
                freq=_new_to_legacy_freq(self.freq),
                closed=self.closed,
                label=self.label,
                origin=self.origin,
                offset=offset,
            )
        self.group_as_index = group_as_index

    def _get_index_and_items(self) -> tuple[pd.Index, pd.Series, np.ndarray]:
        first_items, codes = self.first_items()
        full_index = first_items.index
        if first_items.isnull().any():
            first_items = first_items.dropna()

        full_index = full_index.rename("__resample_dim__")
        return full_index, first_items, codes

    def first_items(self) -> tuple[pd.Series, np.ndarray]:
        from xarray.coding.cftimeindex import CFTimeIndex
        from xarray.core.resample_cftime import CFTimeGrouper

        if isinstance(self.index_grouper, CFTimeGrouper):
            return self.index_grouper.first_items(
                cast(CFTimeIndex, self.group_as_index)
            )
        else:
            s = pd.Series(np.arange(self.group_as_index.size), self.group_as_index)
            grouped = s.groupby(self.index_grouper)
            first_items = grouped.first()
            counts = grouped.count()
            # This way we generate codes for the final output index: full_index.
            # So for _flox_reduce we avoid one reindex and copy by avoiding
            # _maybe_restore_empty_groups
            codes = np.repeat(np.arange(len(first_items)), counts)
            return first_items, codes

    def factorize(self, group: T_Group) -> EncodedGroups:
        self._init_properties(group)
        full_index, first_items, codes_ = self._get_index_and_items()
        sbins = first_items.values.astype(np.int64)
        group_indices: GroupIndices = tuple(
            [slice(i, j) for i, j in zip(sbins[:-1], sbins[1:], strict=True)]
            + [slice(sbins[-1], None)]
        )

        unique_coord = Variable(
            dims=group.name, data=first_items.index, attrs=group.attrs
        )
        codes = group.copy(data=codes_.reshape(group.shape), deep=False)

        return EncodedGroups(
            codes=codes,
            group_indices=group_indices,
            full_index=full_index,
            unique_coord=unique_coord,
            coords=_coordinates_from_variable(unique_coord),
        )


def unique_value_groups(
    ar, sort: bool = True
) -> tuple[np.ndarray | pd.Index, np.ndarray]:
    """Group an array by its unique values.

    Parameters
    ----------
    ar : array-like
        Input array. This will be flattened if it is not already 1-D.
    sort : bool, default: True
        Whether or not to sort unique values.

    Returns
    -------
    values : np.ndarray
        Sorted, unique values as returned by `np.unique`.
    indices : list of lists of int
        Each element provides the integer indices in `ar` with values given by
        the corresponding value in `unique_values`.
    """
    inverse, values = pd.factorize(ar, sort=sort)
    if isinstance(values, pd.MultiIndex):
        values.names = ar.names
    return values, inverse
