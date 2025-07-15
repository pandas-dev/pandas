from __future__ import annotations

import warnings
from collections.abc import Callable, Hashable, Iterable, Sequence
from typing import TYPE_CHECKING, Any, Literal

from xarray.core._aggregations import (
    DataArrayResampleAggregations,
    DatasetResampleAggregations,
)
from xarray.core.groupby import DataArrayGroupByBase, DatasetGroupByBase, GroupBy
from xarray.core.types import Dims, InterpOptions, T_Xarray

if TYPE_CHECKING:
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset
    from xarray.core.types import T_Chunks

from xarray.groupers import RESAMPLE_DIM


class Resample(GroupBy[T_Xarray]):
    """An object that extends the `GroupBy` object with additional logic
    for handling specialized re-sampling operations.

    You should create a `Resample` object by using the `DataArray.resample` or
    `Dataset.resample` methods. The dimension along re-sampling

    See Also
    --------
    DataArray.resample
    Dataset.resample

    """

    def __init__(
        self,
        *args,
        dim: Hashable | None = None,
        resample_dim: Hashable | None = None,
        **kwargs,
    ) -> None:
        if dim == resample_dim:
            raise ValueError(
                f"Proxy resampling dimension ('{resample_dim}') "
                f"cannot have the same name as actual dimension ('{dim}')!"
            )
        self._dim = dim

        super().__init__(*args, **kwargs)

    def _flox_reduce(
        self,
        dim: Dims,
        keep_attrs: bool | None = None,
        **kwargs,
    ) -> T_Xarray:
        result: T_Xarray = (
            super()
            ._flox_reduce(dim=dim, keep_attrs=keep_attrs, **kwargs)
            .rename({RESAMPLE_DIM: self._group_dim})  # type: ignore[assignment]
        )
        return result

    def shuffle_to_chunks(self, chunks: T_Chunks = None):
        """
        Sort or "shuffle" the underlying object.

        "Shuffle" means the object is sorted so that all group members occur sequentially,
        in the same chunk. Multiple groups may occur in the same chunk.
        This method is particularly useful for chunked arrays (e.g. dask, cubed).
        particularly when you need to map a function that requires all members of a group
        to be present in a single chunk. For chunked array types, the order of appearance
        is not guaranteed, but will depend on the input chunking.

        Parameters
        ----------
        chunks : int, tuple of int, "auto" or mapping of hashable to int or tuple of int, optional
            How to adjust chunks along dimensions not present in the array being grouped by.

        Returns
        -------
        DataArrayGroupBy or DatasetGroupBy

        Examples
        --------
        >>> import dask.array
        >>> da = xr.DataArray(
        ...     dims="time",
        ...     data=dask.array.arange(10, chunks=1),
        ...     coords={"time": xr.date_range("2001-01-01", freq="12h", periods=10)},
        ...     name="a",
        ... )
        >>> shuffled = da.resample(time="2D").shuffle_to_chunks()
        >>> shuffled
        <xarray.DataArray 'a' (time: 10)> Size: 80B
        dask.array<shuffle, shape=(10,), dtype=int64, chunksize=(4,), chunktype=numpy.ndarray>
        Coordinates:
          * time     (time) datetime64[ns] 80B 2001-01-01 ... 2001-01-05T12:00:00

        See Also
        --------
        dask.dataframe.DataFrame.shuffle
        dask.array.shuffle
        """
        (grouper,) = self.groupers
        return self._shuffle_obj(chunks).drop_vars(RESAMPLE_DIM)

    def _first_or_last(
        self, op: Literal["first", "last"], skipna: bool | None, keep_attrs: bool | None
    ) -> T_Xarray:
        from xarray.core.dataset import Dataset

        result = super()._first_or_last(op=op, skipna=skipna, keep_attrs=keep_attrs)
        if isinstance(result, Dataset):
            # Can't do this in the base class because group_dim is RESAMPLE_DIM
            # which is not present in the original object
            for var in result.data_vars:
                result._variables[var] = result._variables[var].transpose(
                    *self._obj._variables[var].dims
                )
        return result

    def _drop_coords(self) -> T_Xarray:
        """Drop non-dimension coordinates along the resampled dimension."""
        obj = self._obj
        for k, v in obj.coords.items():
            if k != self._dim and self._dim in v.dims:
                obj = obj.drop_vars([k])
        return obj

    def pad(self, tolerance: float | Iterable[float] | str | None = None) -> T_Xarray:
        """Forward fill new values at up-sampled frequency.

        Parameters
        ----------
        tolerance : float | Iterable[float] | str | None, default: None
            Maximum distance between original and new labels to limit
            the up-sampling method.
            Up-sampled data with indices that satisfy the equation
            ``abs(index[indexer] - target) <= tolerance`` are filled by
            new values. Data with indices that are outside the given
            tolerance are filled with ``NaN`` s.

        Returns
        -------
        padded : DataArray or Dataset
        """
        obj = self._drop_coords()
        (grouper,) = self.groupers
        return obj.reindex(
            {self._dim: grouper.full_index}, method="pad", tolerance=tolerance
        )

    ffill = pad

    def backfill(
        self, tolerance: float | Iterable[float] | str | None = None
    ) -> T_Xarray:
        """Backward fill new values at up-sampled frequency.

        Parameters
        ----------
        tolerance : float | Iterable[float] | str | None, default: None
            Maximum distance between original and new labels to limit
            the up-sampling method.
            Up-sampled data with indices that satisfy the equation
            ``abs(index[indexer] - target) <= tolerance`` are filled by
            new values. Data with indices that are outside the given
            tolerance are filled with ``NaN`` s.

        Returns
        -------
        backfilled : DataArray or Dataset
        """
        obj = self._drop_coords()
        (grouper,) = self.groupers
        return obj.reindex(
            {self._dim: grouper.full_index}, method="backfill", tolerance=tolerance
        )

    bfill = backfill

    def nearest(
        self, tolerance: float | Iterable[float] | str | None = None
    ) -> T_Xarray:
        """Take new values from nearest original coordinate to up-sampled
        frequency coordinates.

        Parameters
        ----------
        tolerance : float | Iterable[float] | str | None, default: None
            Maximum distance between original and new labels to limit
            the up-sampling method.
            Up-sampled data with indices that satisfy the equation
            ``abs(index[indexer] - target) <= tolerance`` are filled by
            new values. Data with indices that are outside the given
            tolerance are filled with ``NaN`` s.

        Returns
        -------
        upsampled : DataArray or Dataset
        """
        obj = self._drop_coords()
        (grouper,) = self.groupers
        return obj.reindex(
            {self._dim: grouper.full_index}, method="nearest", tolerance=tolerance
        )

    def interpolate(self, kind: InterpOptions = "linear", **kwargs) -> T_Xarray:
        """Interpolate up-sampled data using the original data as knots.

        Parameters
        ----------
        kind : {"linear", "nearest", "zero", "slinear", \
                "quadratic", "cubic", "polynomial"}, default: "linear"
            The method used to interpolate. The method should be supported by
            the scipy interpolator:

            - ``interp1d``: {"linear", "nearest", "zero", "slinear",
              "quadratic", "cubic", "polynomial"}
            - ``interpn``: {"linear", "nearest"}

            If ``"polynomial"`` is passed, the ``order`` keyword argument must
            also be provided.

        Returns
        -------
        interpolated : DataArray or Dataset

        See Also
        --------
        DataArray.interp
        Dataset.interp
        scipy.interpolate.interp1d

        """
        return self._interpolate(kind=kind, **kwargs)

    def _interpolate(self, kind="linear", **kwargs) -> T_Xarray:
        """Apply scipy.interpolate.interp1d along resampling dimension."""
        obj = self._drop_coords()
        (grouper,) = self.groupers
        kwargs.setdefault("bounds_error", False)
        return obj.interp(
            coords={self._dim: grouper.full_index},
            assume_sorted=True,
            method=kind,
            kwargs=kwargs,
        )


class DataArrayResample(
    Resample["DataArray"], DataArrayGroupByBase, DataArrayResampleAggregations
):
    """DataArrayGroupBy object specialized to time resampling operations over a
    specified dimension
    """

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        shortcut: bool = True,
        **kwargs: Any,
    ) -> DataArray:
        """Reduce the items in this group by applying `func` along the
        pre-defined resampling dimension.

        Parameters
        ----------
        func : callable
            Function which can be called in the form
            `func(x, axis=axis, **kwargs)` to return the result of collapsing
            an np.ndarray over an integer valued axis.
        dim : "...", str, Iterable of Hashable or None, optional
            Dimension(s) over which to apply `func`.
        keep_attrs : bool, optional
            If True, the datasets's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        **kwargs : dict
            Additional keyword arguments passed on to `func`.

        Returns
        -------
        reduced : DataArray
            Array with summarized data and the indicated dimension(s)
            removed.
        """
        return super().reduce(
            func=func,
            dim=dim,
            axis=axis,
            keep_attrs=keep_attrs,
            keepdims=keepdims,
            shortcut=shortcut,
            **kwargs,
        )

    def map(
        self,
        func: Callable[..., Any],
        args: tuple[Any, ...] = (),
        shortcut: bool | None = False,
        **kwargs: Any,
    ) -> DataArray:
        """Apply a function to each array in the group and concatenate them
        together into a new array.

        `func` is called like `func(ar, *args, **kwargs)` for each array `ar`
        in this group.

        Apply uses heuristics (like `pandas.GroupBy.apply`) to figure out how
        to stack together the array. The rule is:

        1. If the dimension along which the group coordinate is defined is
           still in the first grouped array after applying `func`, then stack
           over this dimension.
        2. Otherwise, stack over the new dimension given by name of this
           grouping (the argument to the `groupby` function).

        Parameters
        ----------
        func : callable
            Callable to apply to each array.
        shortcut : bool, optional
            Whether or not to shortcut evaluation under the assumptions that:

            (1) The action of `func` does not depend on any of the array
                metadata (attributes or coordinates) but only on the data and
                dimensions.
            (2) The action of `func` creates arrays with homogeneous metadata,
                that is, with the same dimensions and attributes.

            If these conditions are satisfied `shortcut` provides significant
            speedup. This should be the case for many common groupby operations
            (e.g., applying numpy ufuncs).
        args : tuple, optional
            Positional arguments passed on to `func`.
        **kwargs
            Used to call `func(ar, **kwargs)` for each array `ar`.

        Returns
        -------
        applied : DataArray
            The result of splitting, applying and combining this array.
        """
        # TODO: the argument order for Resample doesn't match that for its parent,
        # GroupBy
        combined = super().map(func, shortcut=shortcut, args=args, **kwargs)

        # If the aggregation function didn't drop the original resampling
        # dimension, then we need to do so before we can rename the proxy
        # dimension we used.
        if self._dim in combined.coords:
            combined = combined.drop_vars([self._dim])

        if RESAMPLE_DIM in combined.dims:
            combined = combined.rename({RESAMPLE_DIM: self._dim})

        return combined

    def apply(self, func, args=(), shortcut=None, **kwargs):
        """
        Backward compatible implementation of ``map``

        See Also
        --------
        DataArrayResample.map
        """
        warnings.warn(
            "Resample.apply may be deprecated in the future. Using Resample.map is encouraged",
            PendingDeprecationWarning,
            stacklevel=2,
        )
        return self.map(func=func, shortcut=shortcut, args=args, **kwargs)

    def asfreq(self) -> DataArray:
        """Return values of original object at the new up-sampling frequency;
        essentially a re-index with new times set to NaN.

        Returns
        -------
        resampled : DataArray
        """
        self._obj = self._drop_coords()
        return self.mean(None if self._dim is None else [self._dim])


class DatasetResample(
    Resample["Dataset"], DatasetGroupByBase, DatasetResampleAggregations
):
    """DatasetGroupBy object specialized to resampling a specified dimension"""

    def map(
        self,
        func: Callable[..., Any],
        args: tuple[Any, ...] = (),
        shortcut: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """Apply a function over each Dataset in the groups generated for
        resampling and concatenate them together into a new Dataset.

        `func` is called like `func(ds, *args, **kwargs)` for each dataset `ds`
        in this group.

        Apply uses heuristics (like `pandas.GroupBy.apply`) to figure out how
        to stack together the datasets. The rule is:

        1. If the dimension along which the group coordinate is defined is
           still in the first grouped item after applying `func`, then stack
           over this dimension.
        2. Otherwise, stack over the new dimension given by name of this
           grouping (the argument to the `groupby` function).

        Parameters
        ----------
        func : callable
            Callable to apply to each sub-dataset.
        args : tuple, optional
            Positional arguments passed on to `func`.
        **kwargs
            Used to call `func(ds, **kwargs)` for each sub-dataset `ar`.

        Returns
        -------
        applied : Dataset
            The result of splitting, applying and combining this dataset.
        """
        # ignore shortcut if set (for now)
        applied = (func(ds, *args, **kwargs) for ds in self._iter_grouped())
        combined = self._combine(applied)

        # If the aggregation function didn't drop the original resampling
        # dimension, then we need to do so before we can rename the proxy
        # dimension we used.
        if self._dim in combined.coords:
            combined = combined.drop_vars(self._dim)

        if RESAMPLE_DIM in combined.dims:
            combined = combined.rename({RESAMPLE_DIM: self._dim})

        return combined

    def apply(self, func, args=(), shortcut=None, **kwargs):
        """
        Backward compatible implementation of ``map``

        See Also
        --------
        DataSetResample.map
        """

        warnings.warn(
            "Resample.apply may be deprecated in the future. Using Resample.map is encouraged",
            PendingDeprecationWarning,
            stacklevel=2,
        )
        return self.map(func=func, shortcut=shortcut, args=args, **kwargs)

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        shortcut: bool = True,
        **kwargs: Any,
    ) -> Dataset:
        """Reduce the items in this group by applying `func` along the
        pre-defined resampling dimension.

        Parameters
        ----------
        func : callable
            Function which can be called in the form
            `func(x, axis=axis, **kwargs)` to return the result of collapsing
            an np.ndarray over an integer valued axis.
        dim : "...", str, Iterable of Hashable or None, optional
            Dimension(s) over which to apply `func`.
        keep_attrs : bool, optional
            If True, the datasets's attributes (`attrs`) will be copied from
            the original object to the new one.  If False (default), the new
            object will be returned without attributes.
        **kwargs : dict
            Additional keyword arguments passed on to `func`.

        Returns
        -------
        reduced : Dataset
            Array with summarized data and the indicated dimension(s)
            removed.
        """
        return super().reduce(
            func=func,
            dim=dim,
            axis=axis,
            keep_attrs=keep_attrs,
            keepdims=keepdims,
            shortcut=shortcut,
            **kwargs,
        )

    def asfreq(self) -> Dataset:
        """Return values of original object at the new up-sampling frequency;
        essentially a re-index with new times set to NaN.

        Returns
        -------
        resampled : Dataset
        """
        self._obj = self._drop_coords()
        return self.mean(None if self._dim is None else [self._dim])
