"""
This module provides Grouper objects that encapsulate the
"factorization" process - conversion of value we are grouping by
to integer codes (one per group).
"""

from __future__ import annotations

import datetime
import functools
import itertools
import operator
from abc import ABC, abstractmethod
from collections import defaultdict
from collections.abc import Callable, Hashable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import partial
from itertools import chain, pairwise
from typing import TYPE_CHECKING, Any, Literal, cast

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

from xarray.coding.cftime_offsets import BaseCFTimeOffset, _new_to_legacy_freq
from xarray.coding.cftimeindex import CFTimeIndex
from xarray.compat.toolzcompat import sliding_window
from xarray.computation.apply_ufunc import apply_ufunc
from xarray.core.common import (
    _contains_cftime_datetimes,
    _contains_datetime_like_objects,
)
from xarray.core.coordinates import Coordinates, coordinates_from_variable
from xarray.core.dataarray import DataArray
from xarray.core.duck_array_ops import array_all, isnull
from xarray.core.formatting import first_n_items
from xarray.core.groupby import T_Group, _DummyGroup
from xarray.core.indexes import safe_cast_to_index
from xarray.core.resample_cftime import CFTimeGrouper
from xarray.core.types import (
    Bins,
    CFTimeDatetime,
    DatetimeLike,
    GroupIndices,
    PDDatetimeUnitOptions,
    ResampleCompatible,
    Self,
    SideOptions,
)
from xarray.core.variable import Variable
from xarray.namedarray.pycompat import is_chunked_array

__all__ = [
    "BinGrouper",
    "EncodedGroups",
    "Grouper",
    "Resampler",
    "SeasonGrouper",
    "SeasonResampler",
    "TimeResampler",
    "UniqueGrouper",
]

RESAMPLE_DIM = "__resample_dim__"


def _datetime64_via_timestamp(unit: PDDatetimeUnitOptions, **kwargs) -> np.datetime64:
    """Construct a numpy.datetime64 object through the pandas.Timestamp
    constructor with a specific resolution."""
    # TODO: when pandas 3 is our minimum requirement we will no longer need to
    # convert to np.datetime64 values prior to passing to the DatetimeIndex
    # constructor. With pandas < 3 the DatetimeIndex constructor does not
    # infer the resolution from the resolution of the Timestamp values.
    return pd.Timestamp(**kwargs).as_unit(unit).to_numpy()


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
    group_indices: GroupIndices = field(init=False, repr=False)
    unique_coord: Variable | _DummyGroup = field(init=False, repr=False)
    coords: Coordinates = field(init=False, repr=False)

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
            if not is_chunked_array(codes.data):
                self.group_indices = tuple(
                    g
                    for g in _codes_to_group_indices(
                        codes.data.ravel(), len(full_index)
                    )
                    if g
                )
            else:
                # We will not use this when grouping by a chunked array
                self.group_indices = tuple()
        else:
            self.group_indices = group_indices

        if unique_coord is None:
            unique_codes = np.sort(pd.unique(codes.data))
            # Skip the -1 sentinel
            unique_codes = unique_codes[unique_codes >= 0]
            unique_values = full_index[unique_codes]
            self.unique_coord = Variable(
                dims=codes.name, data=unique_values, attrs=codes.attrs
            )
        else:
            self.unique_coord = unique_coord

        if coords is None:
            assert not isinstance(self.unique_coord, _DummyGroup)
            self.coords = coordinates_from_variable(self.unique_coord)
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

    @abstractmethod
    def reset(self) -> Self:
        """
        Creates a new version of this Grouper clearing any caches.
        """
        pass


class Resampler(Grouper):
    """
    Abstract base class for Grouper objects that allow specializing resampling-type GroupBy instructions.

    Currently only used for TimeResampler, but could be used for SpaceResampler in the future.
    """

    def compute_chunks(self, variable: Variable, *, dim: Hashable) -> tuple[int, ...]:
        """
        Compute chunk sizes for this resampler.

        This method should be implemented by subclasses to provide appropriate
        chunking behavior for their specific resampling strategy.

        Parameters
        ----------
        variable : Variable
            The variable being chunked.
        dim : Hashable
            The name of the dimension being chunked.

        Returns
        -------
        tuple[int, ...]
            A tuple of chunk sizes for the dimension.
        """
        raise NotImplementedError("Subclasses must implement compute_chunks method")


@dataclass
class UniqueGrouper(Grouper):
    """
    Grouper object for grouping by a categorical variable.

    Parameters
    ----------
    labels: array-like, optional
        Group labels to aggregate on. This is required when grouping by a chunked array type
        (e.g. dask or cubed) since it is used to construct the coordinate on the output.
        Grouped operations will only be run on the specified group labels. Any group that is not
        present in ``labels`` will be ignored.
    """

    _group_as_index: pd.Index | None = field(default=None, repr=False, init=False)
    labels: ArrayLike | None = field(default=None)

    @property
    def group_as_index(self) -> pd.Index:
        """Caches the group DataArray as a pandas Index."""
        if self._group_as_index is None:
            if self.group.ndim == 1:
                self._group_as_index = self.group.to_index()
            else:
                self._group_as_index = pd.Index(np.array(self.group).ravel())
        return self._group_as_index

    def reset(self) -> Self:
        return type(self)()

    def factorize(self, group: T_Group) -> EncodedGroups:
        self.group = group

        if is_chunked_array(group.data) and self.labels is None:
            raise ValueError(
                "When grouping by a dask array, `labels` must be passed using "
                "a UniqueGrouper object."
            )
        if self.labels is not None:
            return self._factorize_given_labels(group)

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

    def _factorize_given_labels(self, group: T_Group) -> EncodedGroups:
        codes = apply_ufunc(
            _factorize_given_labels,
            group,
            kwargs={"labels": self.labels},
            dask="parallelized",
            output_dtypes=[np.int64],
            keep_attrs=True,
        )
        return EncodedGroups(
            codes=codes,
            full_index=pd.Index(self.labels),  # type: ignore[arg-type]
            unique_coord=Variable(
                dims=codes.name,
                data=self.labels,
                attrs=self.group.attrs,
            ),
        )

    def _factorize_unique(self) -> EncodedGroups:
        # look through group to find the unique values
        sort = not isinstance(self.group_as_index, pd.MultiIndex)
        unique_values, codes_ = unique_value_groups(self.group_as_index, sort=sort)
        if array_all(codes_ == -1):
            raise ValueError(
                "Failed to group data. Are you grouping by a variable that is all NaN?"
            )
        codes = self.group.copy(data=codes_.reshape(self.group.shape), deep=False)
        unique_coord = Variable(
            dims=codes.name, data=unique_values, attrs=self.group.attrs
        )
        full_index = (
            unique_values
            if isinstance(unique_values, pd.MultiIndex)
            else pd.Index(unique_values)
        )

        return EncodedGroups(
            codes=codes,
            full_index=full_index,
            unique_coord=unique_coord,
            coords=coordinates_from_variable(unique_coord),
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
                coords = coordinates_from_variable(unique_coord)

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
        raises an error.
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

    def reset(self) -> Self:
        return type(self)(
            bins=self.bins,
            right=self.right,
            labels=self.labels,
            precision=self.precision,
            include_lowest=self.include_lowest,
            duplicates=self.duplicates,
        )

    def __post_init__(self) -> None:
        if array_all(isnull(self.bins)):
            raise ValueError("All bin edges are NaN.")

    def _cut(self, data):
        return pd.cut(
            np.asarray(data).ravel(),
            bins=self.bins,
            right=self.right,
            labels=self.labels,
            precision=self.precision,
            include_lowest=self.include_lowest,
            duplicates=self.duplicates,
            retbins=True,
        )

    def _pandas_cut_wrapper(self, data, **kwargs):
        binned, bins = self._cut(data)
        if isinstance(self.bins, int):
            # we are running eagerly, update self.bins with actual edges instead
            self.bins = bins
        return binned.codes.reshape(data.shape)

    def factorize(self, group: T_Group) -> EncodedGroups:
        if isinstance(group, _DummyGroup):
            group = DataArray(group.data, dims=group.dims, name=group.name)
        by_is_chunked = is_chunked_array(group.data)
        if isinstance(self.bins, int) and by_is_chunked:
            raise ValueError(
                f"Bin edges must be provided when grouping by chunked arrays. Received {self.bins=!r} instead"
            )
        codes = apply_ufunc(
            self._pandas_cut_wrapper,
            group,
            dask="parallelized",
            keep_attrs=True,
            output_dtypes=[np.int64],
        )
        if not by_is_chunked and array_all(codes == -1):
            raise ValueError(
                f"None of the data falls within bins with edges {self.bins!r}"
            )

        new_dim_name = f"{group.name}_bins"
        codes.name = new_dim_name

        # This seems silly, but it lets us have Pandas handle the complexity
        # of `labels`, `precision`, and `include_lowest`, even when group is a chunked array
        # Pandas ignores labels when IntervalIndex is passed
        if self.labels is None or not isinstance(self.bins, pd.IntervalIndex):
            dummy, _ = self._cut(np.array([0]).astype(group.dtype))
            full_index = dummy.categories
        else:
            full_index = pd.Index(self.labels)

        if not by_is_chunked:
            uniques = np.sort(pd.unique(codes.data.ravel()))
            unique_values = full_index[uniques[uniques != -1]]
        else:
            unique_values = full_index

        unique_coord = Variable(
            dims=new_dim_name, data=unique_values, attrs=group.attrs
        )
        return EncodedGroups(
            codes=codes,
            full_index=full_index,
            unique_coord=unique_coord,
            coords=coordinates_from_variable(unique_coord),
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

    def reset(self) -> Self:
        return type(self)(
            freq=self.freq,
            closed=self.closed,
            label=self.label,
            origin=self.origin,
            offset=self.offset,
        )

    def _init_properties(self, group: T_Group) -> None:
        group_as_index = safe_cast_to_index(group)
        offset = self.offset

        if not group_as_index.is_monotonic_increasing:
            # TODO: sort instead of raising an error
            raise ValueError("Index must be monotonic for resampling")

        if isinstance(group_as_index, CFTimeIndex):
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
            # _maybe_reindex
            codes = np.repeat(np.arange(len(first_items)), counts)
            return first_items, codes

    def factorize(self, group: T_Group) -> EncodedGroups:
        self._init_properties(group)
        full_index, first_items, codes_ = self._get_index_and_items()
        sbins = first_items.values.astype(np.int64)
        group_indices: GroupIndices = tuple(
            list(itertools.starmap(slice, pairwise(sbins))) + [slice(sbins[-1], None)]
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
            coords=coordinates_from_variable(unique_coord),
        )

    def compute_chunks(self, variable: Variable, *, dim: Hashable) -> tuple[int, ...]:
        """
        Compute chunk sizes for this time resampler.

        This method is used during chunking operations to determine appropriate
        chunk sizes for the given variable when using this resampler.

        Parameters
        ----------
        name : Hashable
            The name of the dimension being chunked.
        variable : Variable
            The variable being chunked.

        Returns
        -------
        tuple[int, ...]
            A tuple of chunk sizes for the dimension.
        """
        if not _contains_datetime_like_objects(variable):
            raise ValueError(
                f"Computing chunks with {type(self)!r} only supported for datetime variables. "
                f"Received variable with dtype {variable.dtype!r} instead."
            )

        chunks = (
            DataArray(
                np.ones(variable.shape, dtype=int),
                dims=(dim,),
                coords={dim: variable},
            )
            .resample({dim: self})
            .sum()
        )
        # When bins (binning) or time periods are missing (resampling)
        # we can end up with NaNs. Drop them.
        if chunks.dtype.kind == "f":
            chunks = chunks.dropna(dim).astype(int)
        chunks_tuple: tuple[int, ...] = tuple(chunks.data.tolist())
        return chunks_tuple


def _factorize_given_labels(data: np.ndarray, labels: np.ndarray) -> np.ndarray:
    # Copied from flox
    sorter = np.argsort(labels)
    is_sorted = array_all(sorter == np.arange(sorter.size))
    codes = np.searchsorted(labels, data, sorter=sorter)
    mask = ~np.isin(data, labels) | isnull(data) | (codes == len(labels))
    # codes is the index in to the sorted array.
    # if we didn't want sorting, unsort it back
    if not is_sorted:
        codes[codes == len(labels)] = -1
        codes = sorter[(codes,)]
    codes[mask] = -1
    return codes


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


def season_to_month_tuple(seasons: Sequence[str]) -> tuple[tuple[int, ...], ...]:
    """
    >>> season_to_month_tuple(["DJF", "MAM", "JJA", "SON"])
    ((12, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11))
    >>> season_to_month_tuple(["DJFM", "MAMJ", "JJAS", "SOND"])
    ((12, 1, 2, 3), (3, 4, 5, 6), (6, 7, 8, 9), (9, 10, 11, 12))
    >>> season_to_month_tuple(["DJFM", "SOND"])
    ((12, 1, 2, 3), (9, 10, 11, 12))
    """
    initials = "JFMAMJJASOND"
    starts = {
        "".join(s): i + 1
        for s, i in zip(sliding_window(2, initials + "J"), range(12), strict=True)
    }
    result: list[tuple[int, ...]] = []
    for i, season in enumerate(seasons):
        if len(season) == 1:
            if i < len(seasons) - 1:
                suffix = seasons[i + 1][0]
            else:
                suffix = seasons[0][0]
        else:
            suffix = season[1]

        start = starts[season[0] + suffix]

        month_append = []
        for i in range(len(season[1:])):
            elem = start + i + 1
            month_append.append(elem - 12 * (elem > 12))
        result.append((start,) + tuple(month_append))
    return tuple(result)


def inds_to_season_string(asints: tuple[tuple[int, ...], ...]) -> tuple[str, ...]:
    inits = "JFMAMJJASOND"
    return tuple("".join([inits[i_ - 1] for i_ in t]) for t in asints)


def is_sorted_periodic(lst):
    """Used to verify that seasons provided to SeasonResampler are in order."""
    n = len(lst)

    # Find the wraparound point where the list decreases
    wrap_point = -1
    for i in range(1, n):
        if lst[i] < lst[i - 1]:
            wrap_point = i
            break

    # If no wraparound point is found, the list is already sorted
    if wrap_point == -1:
        return True

    # Check if both parts around the wrap point are sorted
    for i in range(1, wrap_point):
        if lst[i] < lst[i - 1]:
            return False
    for i in range(wrap_point + 1, n):
        if lst[i] < lst[i - 1]:
            return False

    # Check wraparound condition
    return lst[-1] <= lst[0]


@dataclass(kw_only=True, frozen=True)
class SeasonsGroup:
    seasons: tuple[str, ...]
    # tuple[integer months] corresponding to each season
    inds: tuple[tuple[int, ...], ...]
    # integer code for each season, this is not simply range(len(seasons))
    # when the seasons have overlaps
    codes: Sequence[int]


def find_independent_seasons(seasons: Sequence[str]) -> Sequence[SeasonsGroup]:
    """
    Iterates though a list of seasons e.g. ["DJF", "FMA", ...],
    and splits that into multiple sequences of non-overlapping seasons.

    >>> find_independent_seasons(
    ...     ["DJF", "FMA", "AMJ", "JJA", "ASO", "OND"]
    ... )  # doctest: +NORMALIZE_WHITESPACE
    [SeasonsGroup(seasons=('DJF', 'AMJ', 'ASO'), inds=((12, 1, 2), (4, 5, 6), (8, 9, 10)), codes=[0, 2, 4]), SeasonsGroup(seasons=('FMA', 'JJA', 'OND'), inds=((2, 3, 4), (6, 7, 8), (10, 11, 12)), codes=[1, 3, 5])]

    >>> find_independent_seasons(["DJF", "MAM", "JJA", "SON"])
    [SeasonsGroup(seasons=('DJF', 'MAM', 'JJA', 'SON'), inds=((12, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11)), codes=[0, 1, 2, 3])]
    """
    season_inds = season_to_month_tuple(seasons)
    grouped = defaultdict(list)
    codes = defaultdict(list)
    seen: set[tuple[int, ...]] = set()
    # This is quadratic, but the number of seasons is at most 12
    for i, current in enumerate(season_inds):
        # Start with a group
        if current not in seen:
            grouped[i].append(current)
            codes[i].append(i)
            seen.add(current)

        # Loop through remaining groups, and look for overlaps
        for j, second in enumerate(season_inds[i:]):
            if not (set(chain(*grouped[i])) & set(second)) and second not in seen:
                grouped[i].append(second)
                codes[i].append(j + i)
                seen.add(second)
        if len(seen) == len(seasons):
            break
        # found all non-overlapping groups for this row start over

    grouped_ints = tuple(tuple(idx) for idx in grouped.values() if idx)
    return [
        SeasonsGroup(seasons=inds_to_season_string(inds), inds=inds, codes=codes)
        for inds, codes in zip(grouped_ints, codes.values(), strict=False)
    ]


@dataclass
class SeasonGrouper(Grouper):
    """Allows grouping using a custom definition of seasons.

    Parameters
    ----------
    seasons: sequence of str
        List of strings representing seasons. E.g. ``"JF"`` or ``"JJA"`` etc.
        Overlapping seasons are allowed (e.g. ``["DJFM", "MAMJ", "JJAS", "SOND"]``)

    Examples
    --------
    >>> SeasonGrouper(["JF", "MAM", "JJAS", "OND"])
    SeasonGrouper(seasons=['JF', 'MAM', 'JJAS', 'OND'])

    The ordering is preserved

    >>> SeasonGrouper(["MAM", "JJAS", "OND", "JF"])
    SeasonGrouper(seasons=['MAM', 'JJAS', 'OND', 'JF'])

    Overlapping seasons are allowed

    >>> SeasonGrouper(["DJFM", "MAMJ", "JJAS", "SOND"])
    SeasonGrouper(seasons=['DJFM', 'MAMJ', 'JJAS', 'SOND'])
    """

    seasons: Sequence[str]
    # drop_incomplete: bool = field(default=True) # TODO

    def factorize(self, group: T_Group) -> EncodedGroups:
        if TYPE_CHECKING:
            assert not isinstance(group, _DummyGroup)
        if not _contains_datetime_like_objects(group.variable):
            raise ValueError(
                "SeasonGrouper can only be used to group by datetime-like arrays."
            )
        months = group.dt.month.data
        seasons_groups = find_independent_seasons(self.seasons)
        codes_ = np.full((len(seasons_groups),) + group.shape, -1, dtype=np.int8)
        group_indices: list[list[int]] = [[]] * len(self.seasons)
        for axis_index, seasgroup in enumerate(seasons_groups):
            for season_tuple, code in zip(
                seasgroup.inds, seasgroup.codes, strict=False
            ):
                mask = np.isin(months, season_tuple)
                codes_[axis_index, mask] = code
                (indices,) = mask.nonzero()
                group_indices[code] = indices.tolist()

        if np.all(codes_ == -1):
            raise ValueError(
                "Failed to group data. Are you grouping by a variable that is all NaN?"
            )
        needs_dummy_dim = len(seasons_groups) > 1
        codes = DataArray(
            dims=(("__season_dim__",) if needs_dummy_dim else tuple()) + group.dims,
            data=codes_ if needs_dummy_dim else codes_.squeeze(),
            attrs=group.attrs,
            name="season",
        )
        unique_coord = Variable("season", self.seasons, attrs=group.attrs)
        full_index = pd.Index(self.seasons)
        return EncodedGroups(
            codes=codes,
            group_indices=tuple(group_indices),
            unique_coord=unique_coord,
            full_index=full_index,
        )

    def reset(self) -> Self:
        return type(self)(self.seasons)


@dataclass
class SeasonResampler(Resampler):
    """Allows grouping using a custom definition of seasons.

    Parameters
    ----------
    seasons: Sequence[str]
        An ordered list of seasons.
    drop_incomplete: bool
        Whether to drop seasons that are not completely included in the data.
        For example, if a time series starts in Jan-2001, and seasons includes `"DJF"`
        then observations from Jan-2001, and Feb-2001 are ignored in the grouping
        since Dec-2000 isn't present.

    Examples
    --------
    >>> SeasonResampler(["JF", "MAM", "JJAS", "OND"])
    SeasonResampler(seasons=['JF', 'MAM', 'JJAS', 'OND'], drop_incomplete=True)

    >>> SeasonResampler(["DJFM", "AM", "JJA", "SON"])
    SeasonResampler(seasons=['DJFM', 'AM', 'JJA', 'SON'], drop_incomplete=True)
    """

    seasons: Sequence[str]
    drop_incomplete: bool = field(default=True, kw_only=True)
    season_inds: Sequence[Sequence[int]] = field(init=False, repr=False)
    season_tuples: Mapping[str, Sequence[int]] = field(init=False, repr=False)

    def __post_init__(self):
        self.season_inds = season_to_month_tuple(self.seasons)
        all_inds = functools.reduce(operator.add, self.season_inds)
        if len(all_inds) > len(set(all_inds)):
            raise ValueError(
                f"Overlapping seasons are not allowed. Received {self.seasons!r}"
            )
        self.season_tuples = dict(zip(self.seasons, self.season_inds, strict=True))

        if not is_sorted_periodic(list(itertools.chain(*self.season_inds))):
            raise ValueError(
                "Resampling is only supported with sorted seasons. "
                f"Provided seasons {self.seasons!r} are not sorted."
            )

    def factorize(self, group: T_Group) -> EncodedGroups:
        if group.ndim != 1:
            raise ValueError(
                "SeasonResampler can only be used to resample by 1D arrays."
            )
        if not isinstance(group, DataArray) or not _contains_datetime_like_objects(
            group.variable
        ):
            raise ValueError(
                "SeasonResampler can only be used to group by datetime-like DataArrays."
            )

        seasons = self.seasons
        season_inds = self.season_inds
        season_tuples = self.season_tuples

        nstr = max(len(s) for s in seasons)
        year = group.dt.year.astype(int)
        month = group.dt.month.astype(int)
        season_label = np.full(group.shape, "", dtype=f"U{nstr}")

        # offset years for seasons with December and January
        for season_str, season_ind in zip(seasons, season_inds, strict=True):
            season_label[month.isin(season_ind)] = season_str
            if "DJ" in season_str:
                after_dec = season_ind[season_str.index("D") + 1 :]
                # important: this is assuming non-overlapping seasons
                year[month.isin(after_dec)] -= 1

        # Allow users to skip one or more months?
        # present_seasons is a mask that is True for months that are requested in the output
        present_seasons = season_label != ""
        if present_seasons.all():
            # avoid copies if we can.
            present_seasons = slice(None)
        frame = pd.DataFrame(
            data={
                "index": np.arange(group[present_seasons].size),
                "month": month[present_seasons],
            },
            index=pd.MultiIndex.from_arrays(
                [year.data[present_seasons], season_label[present_seasons]],
                names=["year", "season"],
            ),
        )

        agged = (
            frame["index"]
            .groupby(["year", "season"], sort=False)
            .agg(["first", "count"])
        )
        first_items = agged["first"]
        counts = agged["count"]

        index_class: type[CFTimeIndex | pd.DatetimeIndex]
        datetime_class: CFTimeDatetime | Callable[..., np.datetime64]
        if _contains_cftime_datetimes(group.data):
            index_class = CFTimeIndex
            datetime_class = type(first_n_items(group.data, 1).item())
        else:
            index_class = pd.DatetimeIndex
            unit, _ = np.datetime_data(group.dtype)
            unit = cast(PDDatetimeUnitOptions, unit)
            datetime_class = partial(_datetime64_via_timestamp, unit)

        # these are the seasons that are present

        # TODO: when pandas 3 is our minimum requirement we will no longer need
        # to cast the list to a NumPy array prior to passing to the index
        # constructor.
        unique_coord = index_class(
            np.array(
                [
                    datetime_class(year=year, month=season_tuples[season][0], day=1)
                    for year, season in first_items.index
                ]
            )
        )

        # This sorted call is a hack. It's hard to figure out how
        # to start the iteration for arbitrary season ordering
        # for example "DJF" as first entry or last entry
        # So we construct the largest possible index and slice it to the
        # range present in the data.

        # TODO: when pandas 3 is our minimum requirement we will no longer need
        # to cast the list to a NumPy array prior to passing to the index
        # constructor.
        complete_index = index_class(
            np.array(
                sorted(
                    [
                        datetime_class(year=y, month=m, day=1)
                        for y, m in itertools.product(
                            range(year[0].item(), year[-1].item() + 1),
                            [s[0] for s in season_inds],
                        )
                    ]
                )
            )
        )

        # all years and seasons
        def get_label(year, season):
            month, *_ = season_tuples[season]
            return f"{year}-{month:02d}-01"

        unique_codes = np.arange(len(unique_coord))
        valid_season_mask = season_label != ""
        first_valid_season, last_valid_season = season_label[valid_season_mask][[0, -1]]
        first_year, last_year = year.data[[0, -1]]
        if self.drop_incomplete:
            if month.data[valid_season_mask][0] != season_tuples[first_valid_season][0]:
                if "DJ" in first_valid_season:
                    first_year += 1
                first_valid_season = seasons[
                    (seasons.index(first_valid_season) + 1) % len(seasons)
                ]
                unique_codes -= 1

            if (
                month.data[valid_season_mask][-1]
                != season_tuples[last_valid_season][-1]
            ):
                last_valid_season = seasons[seasons.index(last_valid_season) - 1]
                if "DJ" in last_valid_season:
                    last_year -= 1
                unique_codes[-1] = -1

        first_label = get_label(first_year, first_valid_season)
        last_label = get_label(last_year, last_valid_season)

        slicer = complete_index.slice_indexer(first_label, last_label)
        full_index = complete_index[slicer]

        final_codes = np.full(group.data.size, -1)
        final_codes[present_seasons] = np.repeat(unique_codes, counts)
        codes = group.copy(data=final_codes, deep=False)

        return EncodedGroups(codes=codes, full_index=full_index)

    def compute_chunks(self, variable: Variable, *, dim: Hashable) -> tuple[int, ...]:
        """
        Compute chunk sizes for this season resampler.

        This method is used during chunking operations to determine appropriate
        chunk sizes for the given variable when using this resampler.

        Parameters
        ----------
        name : Hashable
            The name of the dimension being chunked.
        variable : Variable
            The variable being chunked.

        Returns
        -------
        tuple[int, ...]
            A tuple of chunk sizes for the dimension.
        """
        if not _contains_datetime_like_objects(variable):
            raise ValueError(
                f"Computing chunks with {type(self)!r} only supported for datetime variables. "
                f"Received variable with dtype {variable.dtype!r} instead."
            )

        if len("".join(self.seasons)) != 12:
            raise ValueError(
                "Cannot rechunk with a SeasonResampler that does not cover all 12 months. "
                f"Received `seasons={self.seasons!r}`."
            )

        # Create a temporary resampler that ignores drop_incomplete for chunking
        # This prevents data from being silently dropped during chunking
        resampler_for_chunking = type(self)(seasons=self.seasons, drop_incomplete=False)

        chunks = (
            DataArray(
                np.ones(variable.shape, dtype=int),
                dims=(dim,),
                coords={dim: variable},
            )
            .resample({dim: resampler_for_chunking})
            .sum()
        )
        # When bins (binning) or time periods are missing (resampling)
        # we can end up with NaNs. Drop them.
        if chunks.dtype.kind == "f":
            chunks = chunks.dropna(dim).astype(int)
        chunks_tuple: tuple[int, ...] = tuple(chunks.data.tolist())
        return chunks_tuple

    def reset(self) -> Self:
        return type(self)(seasons=self.seasons, drop_incomplete=self.drop_incomplete)
