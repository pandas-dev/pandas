"""Mixin classes with reduction operations."""

# This file was generated using xarray.util.generate_aggregations. Do not edit manually.

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING, Any

from xarray.core import duck_array_ops
from xarray.core.options import OPTIONS
from xarray.core.types import Dims, Self
from xarray.core.utils import contains_only_chunked_or_numpy, module_available

if TYPE_CHECKING:
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

flox_available = module_available("flox")


class DataTreeAggregations:
    __slots__ = ()

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> Self:
        raise NotImplementedError()

    def count(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``count`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``count``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``count`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``count`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        pandas.DataFrame.count
        dask.dataframe.DataFrame.count
        Dataset.count
        DataArray.count
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.count()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      int64 8B 5
        """
        return self.reduce(
            duck_array_ops.count,
            dim=dim,
            numeric_only=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def all(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``all`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``all``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``all`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``all`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.all
        dask.array.all
        Dataset.all
        DataArray.all
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(
        ...             foo=(
        ...                 "time",
        ...                 np.array([True, True, True, True, True, False], dtype=bool),
        ...             )
        ...         ),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) bool 6B True True True True True False

        >>> dt.all()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      bool 1B False
        """
        return self.reduce(
            duck_array_ops.array_all,
            dim=dim,
            numeric_only=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def any(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``any`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``any``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``any`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``any`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.any
        dask.array.any
        Dataset.any
        DataArray.any
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(
        ...             foo=(
        ...                 "time",
        ...                 np.array([True, True, True, True, True, False], dtype=bool),
        ...             )
        ...         ),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) bool 6B True True True True True False

        >>> dt.any()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      bool 1B True
        """
        return self.reduce(
            duck_array_ops.array_any,
            dim=dim,
            numeric_only=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def max(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``max`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``max``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``max`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``max`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.max
        dask.array.max
        Dataset.max
        DataArray.max
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.max()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 3.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> dt.max(skipna=False)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B nan
        """
        return self.reduce(
            duck_array_ops.max,
            dim=dim,
            skipna=skipna,
            numeric_only=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def min(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``min`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``min``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``min`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``min`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.min
        dask.array.min
        Dataset.min
        DataArray.min
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.min()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 0.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> dt.min(skipna=False)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B nan
        """
        return self.reduce(
            duck_array_ops.min,
            dim=dim,
            skipna=skipna,
            numeric_only=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def mean(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``mean`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``mean``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``mean`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``mean`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.mean
        dask.array.mean
        Dataset.mean
        DataArray.mean
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.mean()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 1.6

        Use ``skipna`` to control whether NaNs are ignored.

        >>> dt.mean(skipna=False)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B nan
        """
        return self.reduce(
            duck_array_ops.mean,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def prod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``prod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``prod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``prod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``prod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.prod
        dask.array.prod
        Dataset.prod
        DataArray.prod
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.prod()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 0.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> dt.prod(skipna=False)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B nan

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> dt.prod(skipna=True, min_count=2)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 0.0
        """
        return self.reduce(
            duck_array_ops.prod,
            dim=dim,
            skipna=skipna,
            min_count=min_count,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def sum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``sum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``sum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``sum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``sum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.sum
        dask.array.sum
        Dataset.sum
        DataArray.sum
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.sum()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 8.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> dt.sum(skipna=False)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B nan

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> dt.sum(skipna=True, min_count=2)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 8.0
        """
        return self.reduce(
            duck_array_ops.sum,
            dim=dim,
            skipna=skipna,
            min_count=min_count,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def std(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``std`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``std``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``std`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``std`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.std
        dask.array.std
        Dataset.std
        DataArray.std
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.std()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 1.02

        Use ``skipna`` to control whether NaNs are ignored.

        >>> dt.std(skipna=False)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B nan

        Specify ``ddof=1`` for an unbiased estimate.

        >>> dt.std(skipna=True, ddof=1)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 1.14
        """
        return self.reduce(
            duck_array_ops.std,
            dim=dim,
            skipna=skipna,
            ddof=ddof,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def var(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``var`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``var``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``var`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``var`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.var
        dask.array.var
        Dataset.var
        DataArray.var
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.var()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 1.04

        Use ``skipna`` to control whether NaNs are ignored.

        >>> dt.var(skipna=False)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B nan

        Specify ``ddof=1`` for an unbiased estimate.

        >>> dt.var(skipna=True, ddof=1)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 1.3
        """
        return self.reduce(
            duck_array_ops.var,
            dim=dim,
            skipna=skipna,
            ddof=ddof,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def median(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``median`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``median``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``median`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``median`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.median
        dask.array.median
        Dataset.median
        DataArray.median
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.median()
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B 2.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> dt.median(skipna=False)
        <xarray.DataTree>
        Group: /
            Dimensions:  ()
            Data variables:
                foo      float64 8B nan
        """
        return self.reduce(
            duck_array_ops.median,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumsum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``cumsum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumsum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumsum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``cumsum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumsum
        dask.array.cumsum
        Dataset.cumsum
        DataArray.cumsum
        DataTree.cumulative
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.cumsum()
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Dimensions without coordinates: time
            Data variables:
                foo      (time) float64 48B 1.0 3.0 6.0 6.0 8.0 8.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> dt.cumsum(skipna=False)
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Dimensions without coordinates: time
            Data variables:
                foo      (time) float64 48B 1.0 3.0 6.0 6.0 8.0 nan
        """
        return self.reduce(
            duck_array_ops.cumsum,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumprod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataTree's data by applying ``cumprod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumprod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumprod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataTree
            New DataTree with ``cumprod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumprod
        dask.array.cumprod
        Dataset.cumprod
        DataArray.cumprod
        DataTree.cumulative
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> dt = xr.DataTree(
        ...     xr.Dataset(
        ...         data_vars=dict(foo=("time", np.array([1, 2, 3, 0, 2, np.nan]))),
        ...         coords=dict(
        ...             time=(
        ...                 "time",
        ...                 pd.date_range("2001-01-01", freq="ME", periods=6),
        ...             ),
        ...             labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...         ),
        ...     ),
        ... )
        >>> dt
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Coordinates:
              * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
                labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
            Data variables:
                foo      (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> dt.cumprod()
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Dimensions without coordinates: time
            Data variables:
                foo      (time) float64 48B 1.0 2.0 6.0 0.0 0.0 0.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> dt.cumprod(skipna=False)
        <xarray.DataTree>
        Group: /
            Dimensions:  (time: 6)
            Dimensions without coordinates: time
            Data variables:
                foo      (time) float64 48B 1.0 2.0 6.0 0.0 0.0 nan
        """
        return self.reduce(
            duck_array_ops.cumprod,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )


class DatasetAggregations:
    __slots__ = ()

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> Self:
        raise NotImplementedError()

    def count(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``count`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``count``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``count`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``count`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        pandas.DataFrame.count
        dask.dataframe.DataFrame.count
        DataArray.count
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.count()
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       int64 8B 5
        """
        return self.reduce(
            duck_array_ops.count,
            dim=dim,
            numeric_only=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def all(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``all`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``all``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``all`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``all`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.all
        dask.array.all
        DataArray.all
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 78B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) bool 6B True True True True True False

        >>> ds.all()
        <xarray.Dataset> Size: 1B
        Dimensions:  ()
        Data variables:
            da       bool 1B False
        """
        return self.reduce(
            duck_array_ops.array_all,
            dim=dim,
            numeric_only=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def any(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``any`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``any``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``any`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``any`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.any
        dask.array.any
        DataArray.any
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 78B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) bool 6B True True True True True False

        >>> ds.any()
        <xarray.Dataset> Size: 1B
        Dimensions:  ()
        Data variables:
            da       bool 1B True
        """
        return self.reduce(
            duck_array_ops.array_any,
            dim=dim,
            numeric_only=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def max(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``max`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``max``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``max`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``max`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.max
        dask.array.max
        DataArray.max
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.max()
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 3.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.max(skipna=False)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B nan
        """
        return self.reduce(
            duck_array_ops.max,
            dim=dim,
            skipna=skipna,
            numeric_only=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def min(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``min`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``min``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``min`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``min`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.min
        dask.array.min
        DataArray.min
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.min()
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 0.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.min(skipna=False)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B nan
        """
        return self.reduce(
            duck_array_ops.min,
            dim=dim,
            skipna=skipna,
            numeric_only=False,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def mean(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``mean`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``mean``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``mean`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``mean`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.mean
        dask.array.mean
        DataArray.mean
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.mean()
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 1.6

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.mean(skipna=False)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B nan
        """
        return self.reduce(
            duck_array_ops.mean,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def prod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``prod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``prod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``prod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``prod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.prod
        dask.array.prod
        DataArray.prod
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.prod()
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 0.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.prod(skipna=False)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B nan

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> ds.prod(skipna=True, min_count=2)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 0.0
        """
        return self.reduce(
            duck_array_ops.prod,
            dim=dim,
            skipna=skipna,
            min_count=min_count,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def sum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``sum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``sum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``sum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``sum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.sum
        dask.array.sum
        DataArray.sum
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.sum()
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 8.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.sum(skipna=False)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B nan

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> ds.sum(skipna=True, min_count=2)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 8.0
        """
        return self.reduce(
            duck_array_ops.sum,
            dim=dim,
            skipna=skipna,
            min_count=min_count,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def std(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``std`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``std``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``std`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``std`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.std
        dask.array.std
        DataArray.std
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.std()
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 1.02

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.std(skipna=False)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B nan

        Specify ``ddof=1`` for an unbiased estimate.

        >>> ds.std(skipna=True, ddof=1)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 1.14
        """
        return self.reduce(
            duck_array_ops.std,
            dim=dim,
            skipna=skipna,
            ddof=ddof,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def var(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``var`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``var``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``var`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``var`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.var
        dask.array.var
        DataArray.var
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.var()
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 1.04

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.var(skipna=False)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B nan

        Specify ``ddof=1`` for an unbiased estimate.

        >>> ds.var(skipna=True, ddof=1)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 1.3
        """
        return self.reduce(
            duck_array_ops.var,
            dim=dim,
            skipna=skipna,
            ddof=ddof,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def median(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``median`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``median``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``median`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``median`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.median
        dask.array.median
        DataArray.median
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.median()
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B 2.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.median(skipna=False)
        <xarray.Dataset> Size: 8B
        Dimensions:  ()
        Data variables:
            da       float64 8B nan
        """
        return self.reduce(
            duck_array_ops.median,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumsum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``cumsum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumsum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumsum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``cumsum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumsum
        dask.array.cumsum
        DataArray.cumsum
        Dataset.cumulative
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.cumsum()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 3.0 6.0 6.0 8.0 8.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.cumsum(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 3.0 6.0 6.0 8.0 nan
        """
        return self.reduce(
            duck_array_ops.cumsum,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumprod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this Dataset's data by applying ``cumprod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumprod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumprod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``cumprod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumprod
        dask.array.cumprod
        DataArray.cumprod
        Dataset.cumulative
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.cumprod()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 2.0 6.0 0.0 0.0 0.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.cumprod(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 2.0 6.0 0.0 0.0 nan
        """
        return self.reduce(
            duck_array_ops.cumprod,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )


class DataArrayAggregations:
    __slots__ = ()

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> Self:
        raise NotImplementedError()

    def count(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``count`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``count``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``count`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``count`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        pandas.DataFrame.count
        dask.dataframe.DataFrame.count
        Dataset.count
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.count()
        <xarray.DataArray ()> Size: 8B
        array(5)
        """
        return self.reduce(
            duck_array_ops.count,
            dim=dim,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def all(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``all`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``all``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``all`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``all`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.all
        dask.array.all
        Dataset.all
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 6B
        array([ True,  True,  True,  True,  True, False])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.all()
        <xarray.DataArray ()> Size: 1B
        array(False)
        """
        return self.reduce(
            duck_array_ops.array_all,
            dim=dim,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def any(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``any`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``any``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``any`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``any`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.any
        dask.array.any
        Dataset.any
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 6B
        array([ True,  True,  True,  True,  True, False])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.any()
        <xarray.DataArray ()> Size: 1B
        array(True)
        """
        return self.reduce(
            duck_array_ops.array_any,
            dim=dim,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def max(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``max`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``max``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``max`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``max`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.max
        dask.array.max
        Dataset.max
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.max()
        <xarray.DataArray ()> Size: 8B
        array(3.)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.max(skipna=False)
        <xarray.DataArray ()> Size: 8B
        array(nan)
        """
        return self.reduce(
            duck_array_ops.max,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def min(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``min`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``min``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``min`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``min`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.min
        dask.array.min
        Dataset.min
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.min()
        <xarray.DataArray ()> Size: 8B
        array(0.)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.min(skipna=False)
        <xarray.DataArray ()> Size: 8B
        array(nan)
        """
        return self.reduce(
            duck_array_ops.min,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def mean(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``mean`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``mean``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``mean`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``mean`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.mean
        dask.array.mean
        Dataset.mean
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.mean()
        <xarray.DataArray ()> Size: 8B
        array(1.6)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.mean(skipna=False)
        <xarray.DataArray ()> Size: 8B
        array(nan)
        """
        return self.reduce(
            duck_array_ops.mean,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def prod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``prod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``prod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``prod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``prod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.prod
        dask.array.prod
        Dataset.prod
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.prod()
        <xarray.DataArray ()> Size: 8B
        array(0.)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.prod(skipna=False)
        <xarray.DataArray ()> Size: 8B
        array(nan)

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> da.prod(skipna=True, min_count=2)
        <xarray.DataArray ()> Size: 8B
        array(0.)
        """
        return self.reduce(
            duck_array_ops.prod,
            dim=dim,
            skipna=skipna,
            min_count=min_count,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def sum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``sum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``sum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``sum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``sum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.sum
        dask.array.sum
        Dataset.sum
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.sum()
        <xarray.DataArray ()> Size: 8B
        array(8.)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.sum(skipna=False)
        <xarray.DataArray ()> Size: 8B
        array(nan)

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> da.sum(skipna=True, min_count=2)
        <xarray.DataArray ()> Size: 8B
        array(8.)
        """
        return self.reduce(
            duck_array_ops.sum,
            dim=dim,
            skipna=skipna,
            min_count=min_count,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def std(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``std`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``std``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``std`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``std`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.std
        dask.array.std
        Dataset.std
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.std()
        <xarray.DataArray ()> Size: 8B
        array(1.0198039)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.std(skipna=False)
        <xarray.DataArray ()> Size: 8B
        array(nan)

        Specify ``ddof=1`` for an unbiased estimate.

        >>> da.std(skipna=True, ddof=1)
        <xarray.DataArray ()> Size: 8B
        array(1.14017543)
        """
        return self.reduce(
            duck_array_ops.std,
            dim=dim,
            skipna=skipna,
            ddof=ddof,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def var(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``var`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``var``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``var`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``var`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.var
        dask.array.var
        Dataset.var
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.var()
        <xarray.DataArray ()> Size: 8B
        array(1.04)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.var(skipna=False)
        <xarray.DataArray ()> Size: 8B
        array(nan)

        Specify ``ddof=1`` for an unbiased estimate.

        >>> da.var(skipna=True, ddof=1)
        <xarray.DataArray ()> Size: 8B
        array(1.3)
        """
        return self.reduce(
            duck_array_ops.var,
            dim=dim,
            skipna=skipna,
            ddof=ddof,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def median(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``median`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``median``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``median`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``median`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.median
        dask.array.median
        Dataset.median
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.median()
        <xarray.DataArray ()> Size: 8B
        array(2.)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.median(skipna=False)
        <xarray.DataArray ()> Size: 8B
        array(nan)
        """
        return self.reduce(
            duck_array_ops.median,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumsum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``cumsum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumsum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumsum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``cumsum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumsum
        dask.array.cumsum
        Dataset.cumsum
        DataArray.cumulative
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.cumsum()
        <xarray.DataArray (time: 6)> Size: 48B
        array([1., 3., 6., 6., 8., 8.])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.cumsum(skipna=False)
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  3.,  6.,  6.,  8., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        """
        return self.reduce(
            duck_array_ops.cumsum,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumprod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this DataArray's data by applying ``cumprod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumprod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumprod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``cumprod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumprod
        dask.array.cumprod
        Dataset.cumprod
        DataArray.cumulative
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Notes
        -----
        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.cumprod()
        <xarray.DataArray (time: 6)> Size: 48B
        array([1., 2., 6., 0., 0., 0.])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.cumprod(skipna=False)
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  6.,  0.,  0., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        """
        return self.reduce(
            duck_array_ops.cumprod,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )


class DatasetGroupByAggregations:
    _obj: Dataset

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> Dataset:
        raise NotImplementedError()

    def _flox_reduce(
        self,
        dim: Dims,
        **kwargs: Any,
    ) -> Dataset:
        raise NotImplementedError()

    def count(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``count`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``count``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``count`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``count`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        pandas.DataFrame.count
        dask.dataframe.DataFrame.count
        Dataset.count
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").count()
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) int64 24B 1 2 2
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="count",
                dim=dim,
                numeric_only=False,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.count,
                dim=dim,
                numeric_only=False,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def all(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``all`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``all``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``all`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``all`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.all
        dask.array.all
        Dataset.all
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 78B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) bool 6B True True True True True False

        >>> ds.groupby("labels").all()
        <xarray.Dataset> Size: 27B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) bool 3B False True True
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="all",
                dim=dim,
                numeric_only=False,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.array_all,
                dim=dim,
                numeric_only=False,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def any(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``any`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``any``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``any`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``any`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.any
        dask.array.any
        Dataset.any
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 78B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) bool 6B True True True True True False

        >>> ds.groupby("labels").any()
        <xarray.Dataset> Size: 27B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) bool 3B True True True
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="any",
                dim=dim,
                numeric_only=False,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.array_any,
                dim=dim,
                numeric_only=False,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def max(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``max`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``max``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``max`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``max`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.max
        dask.array.max
        Dataset.max
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").max()
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B 1.0 2.0 3.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.groupby("labels").max(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 2.0 3.0
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="max",
                dim=dim,
                skipna=skipna,
                numeric_only=False,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.max,
                dim=dim,
                skipna=skipna,
                numeric_only=False,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def min(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``min`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``min``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``min`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``min`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.min
        dask.array.min
        Dataset.min
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").min()
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B 1.0 2.0 0.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.groupby("labels").min(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 2.0 0.0
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="min",
                dim=dim,
                skipna=skipna,
                numeric_only=False,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.min,
                dim=dim,
                skipna=skipna,
                numeric_only=False,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def mean(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``mean`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``mean``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``mean`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``mean`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.mean
        dask.array.mean
        Dataset.mean
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").mean()
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B 1.0 2.0 1.5

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.groupby("labels").mean(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 2.0 1.5
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="mean",
                dim=dim,
                skipna=skipna,
                numeric_only=True,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.mean,
                dim=dim,
                skipna=skipna,
                numeric_only=True,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def prod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``prod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``prod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``prod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``prod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.prod
        dask.array.prod
        Dataset.prod
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").prod()
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B 1.0 4.0 0.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.groupby("labels").prod(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 4.0 0.0

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> ds.groupby("labels").prod(skipna=True, min_count=2)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 4.0 0.0
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="prod",
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                numeric_only=True,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.prod,
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                numeric_only=True,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def sum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``sum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``sum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``sum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``sum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.sum
        dask.array.sum
        Dataset.sum
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").sum()
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B 1.0 4.0 3.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.groupby("labels").sum(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 4.0 3.0

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> ds.groupby("labels").sum(skipna=True, min_count=2)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 4.0 3.0
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="sum",
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                numeric_only=True,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.sum,
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                numeric_only=True,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def std(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``std`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``std``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``std`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``std`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.std
        dask.array.std
        Dataset.std
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").std()
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B 0.0 0.0 1.5

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.groupby("labels").std(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 0.0 1.5

        Specify ``ddof=1`` for an unbiased estimate.

        >>> ds.groupby("labels").std(skipna=True, ddof=1)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 0.0 2.121
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="std",
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                numeric_only=True,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.std,
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                numeric_only=True,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def var(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``var`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``var``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``var`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``var`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.var
        dask.array.var
        Dataset.var
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").var()
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B 0.0 0.0 2.25

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.groupby("labels").var(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 0.0 2.25

        Specify ``ddof=1`` for an unbiased estimate.

        >>> ds.groupby("labels").var(skipna=True, ddof=1)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 0.0 4.5
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="var",
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                numeric_only=True,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.var,
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                numeric_only=True,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def median(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``median`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``median``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``median`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``median`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.median
        dask.array.median
        Dataset.median
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").median()
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B 1.0 2.0 1.5

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.groupby("labels").median(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (labels: 3)
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        Data variables:
            da       (labels) float64 24B nan 2.0 1.5
        """
        return self.reduce(
            duck_array_ops.median,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumsum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``cumsum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumsum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumsum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``cumsum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumsum
        dask.array.cumsum
        Dataset.cumsum
        Dataset.cumulative
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").cumsum()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 3.0 4.0 1.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.groupby("labels").cumsum(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 3.0 4.0 nan
        """
        return self.reduce(
            duck_array_ops.cumsum,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumprod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``cumprod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumprod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumprod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``cumprod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumprod
        dask.array.cumprod
        Dataset.cumprod
        Dataset.cumulative
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.groupby("labels").cumprod()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 4.0 1.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.groupby("labels").cumprod(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 4.0 nan
        """
        return self.reduce(
            duck_array_ops.cumprod,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )


class DatasetResampleAggregations:
    _obj: Dataset

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> Dataset:
        raise NotImplementedError()

    def _flox_reduce(
        self,
        dim: Dims,
        **kwargs: Any,
    ) -> Dataset:
        raise NotImplementedError()

    def count(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``count`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``count``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``count`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``count`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        pandas.DataFrame.count
        dask.dataframe.DataFrame.count
        Dataset.count
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").count()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) int64 24B 1 3 1
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="count",
                dim=dim,
                numeric_only=False,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.count,
                dim=dim,
                numeric_only=False,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def all(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``all`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``all``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``all`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``all`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.all
        dask.array.all
        Dataset.all
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 78B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) bool 6B True True True True True False

        >>> ds.resample(time="3ME").all()
        <xarray.Dataset> Size: 27B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) bool 3B True True False
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="all",
                dim=dim,
                numeric_only=False,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.array_all,
                dim=dim,
                numeric_only=False,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def any(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``any`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``any``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``any`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``any`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.any
        dask.array.any
        Dataset.any
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 78B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) bool 6B True True True True True False

        >>> ds.resample(time="3ME").any()
        <xarray.Dataset> Size: 27B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) bool 3B True True True
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="any",
                dim=dim,
                numeric_only=False,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.array_any,
                dim=dim,
                numeric_only=False,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def max(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``max`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``max``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``max`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``max`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.max
        dask.array.max
        Dataset.max
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").max()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 3.0 2.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.resample(time="3ME").max(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 3.0 nan
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="max",
                dim=dim,
                skipna=skipna,
                numeric_only=False,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.max,
                dim=dim,
                skipna=skipna,
                numeric_only=False,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def min(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``min`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``min``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``min`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``min`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.min
        dask.array.min
        Dataset.min
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").min()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 0.0 2.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.resample(time="3ME").min(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 0.0 nan
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="min",
                dim=dim,
                skipna=skipna,
                numeric_only=False,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.min,
                dim=dim,
                skipna=skipna,
                numeric_only=False,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def mean(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``mean`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``mean``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``mean`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``mean`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.mean
        dask.array.mean
        Dataset.mean
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").mean()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 1.667 2.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.resample(time="3ME").mean(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 1.667 nan
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="mean",
                dim=dim,
                skipna=skipna,
                numeric_only=True,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.mean,
                dim=dim,
                skipna=skipna,
                numeric_only=True,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def prod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``prod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``prod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``prod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``prod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.prod
        dask.array.prod
        Dataset.prod
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").prod()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 0.0 2.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.resample(time="3ME").prod(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 0.0 nan

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> ds.resample(time="3ME").prod(skipna=True, min_count=2)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B nan 0.0 nan
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="prod",
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                numeric_only=True,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.prod,
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                numeric_only=True,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def sum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``sum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``sum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``sum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``sum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.sum
        dask.array.sum
        Dataset.sum
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").sum()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 5.0 2.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.resample(time="3ME").sum(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 5.0 nan

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> ds.resample(time="3ME").sum(skipna=True, min_count=2)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B nan 5.0 nan
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="sum",
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                numeric_only=True,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.sum,
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                numeric_only=True,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def std(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``std`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``std``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``std`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``std`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.std
        dask.array.std
        Dataset.std
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").std()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 0.0 1.247 0.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.resample(time="3ME").std(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 0.0 1.247 nan

        Specify ``ddof=1`` for an unbiased estimate.

        >>> ds.resample(time="3ME").std(skipna=True, ddof=1)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B nan 1.528 nan
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="std",
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                numeric_only=True,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.std,
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                numeric_only=True,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def var(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``var`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``var``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``var`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``var`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.var
        dask.array.var
        Dataset.var
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").var()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 0.0 1.556 0.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.resample(time="3ME").var(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 0.0 1.556 nan

        Specify ``ddof=1`` for an unbiased estimate.

        >>> ds.resample(time="3ME").var(skipna=True, ddof=1)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B nan 2.333 nan
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="var",
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                numeric_only=True,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.var,
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                numeric_only=True,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def median(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``median`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``median``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``median`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``median`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.median
        dask.array.median
        Dataset.median
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").median()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 2.0 2.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.resample(time="3ME").median(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 3)
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        Data variables:
            da       (time) float64 24B 1.0 2.0 nan
        """
        return self.reduce(
            duck_array_ops.median,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumsum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``cumsum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumsum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumsum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``cumsum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumsum
        dask.array.cumsum
        Dataset.cumsum
        Dataset.cumulative
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").cumsum()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 2.0 5.0 5.0 2.0 2.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.resample(time="3ME").cumsum(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 2.0 5.0 5.0 2.0 nan
        """
        return self.reduce(
            duck_array_ops.cumsum,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumprod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> Dataset:
        """
        Reduce this Dataset's data by applying ``cumprod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumprod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumprod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : Dataset
            New Dataset with ``cumprod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumprod
        dask.array.cumprod
        Dataset.cumprod
        Dataset.cumulative
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> ds = xr.Dataset(dict(da=da))
        >>> ds
        <xarray.Dataset> Size: 120B
        Dimensions:  (time: 6)
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Data variables:
            da       (time) float64 48B 1.0 2.0 3.0 0.0 2.0 nan

        >>> ds.resample(time="3ME").cumprod()
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 2.0 6.0 0.0 2.0 2.0

        Use ``skipna`` to control whether NaNs are ignored.

        >>> ds.resample(time="3ME").cumprod(skipna=False)
        <xarray.Dataset> Size: 48B
        Dimensions:  (time: 6)
        Dimensions without coordinates: time
        Data variables:
            da       (time) float64 48B 1.0 2.0 6.0 0.0 2.0 nan
        """
        return self.reduce(
            duck_array_ops.cumprod,
            dim=dim,
            skipna=skipna,
            numeric_only=True,
            keep_attrs=keep_attrs,
            **kwargs,
        )


class DataArrayGroupByAggregations:
    _obj: DataArray

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> DataArray:
        raise NotImplementedError()

    def _flox_reduce(
        self,
        dim: Dims,
        **kwargs: Any,
    ) -> DataArray:
        raise NotImplementedError()

    def count(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``count`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``count``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``count`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``count`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        pandas.DataFrame.count
        dask.dataframe.DataFrame.count
        DataArray.count
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").count()
        <xarray.DataArray (labels: 3)> Size: 24B
        array([1, 2, 2])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="count",
                dim=dim,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.count,
                dim=dim,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def all(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``all`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``all``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``all`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``all`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.all
        dask.array.all
        DataArray.all
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 6B
        array([ True,  True,  True,  True,  True, False])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").all()
        <xarray.DataArray (labels: 3)> Size: 3B
        array([False,  True,  True])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="all",
                dim=dim,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.array_all,
                dim=dim,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def any(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``any`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``any``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``any`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``any`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.any
        dask.array.any
        DataArray.any
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 6B
        array([ True,  True,  True,  True,  True, False])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").any()
        <xarray.DataArray (labels: 3)> Size: 3B
        array([ True,  True,  True])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="any",
                dim=dim,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.array_any,
                dim=dim,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def max(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``max`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``max``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``max`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``max`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.max
        dask.array.max
        DataArray.max
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").max()
        <xarray.DataArray (labels: 3)> Size: 24B
        array([1., 2., 3.])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.groupby("labels").max(skipna=False)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([nan,  2.,  3.])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="max",
                dim=dim,
                skipna=skipna,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.max,
                dim=dim,
                skipna=skipna,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def min(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``min`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``min``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``min`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``min`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.min
        dask.array.min
        DataArray.min
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").min()
        <xarray.DataArray (labels: 3)> Size: 24B
        array([1., 2., 0.])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.groupby("labels").min(skipna=False)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([nan,  2.,  0.])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="min",
                dim=dim,
                skipna=skipna,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.min,
                dim=dim,
                skipna=skipna,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def mean(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``mean`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``mean``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``mean`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``mean`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.mean
        dask.array.mean
        DataArray.mean
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").mean()
        <xarray.DataArray (labels: 3)> Size: 24B
        array([1. , 2. , 1.5])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.groupby("labels").mean(skipna=False)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([nan, 2. , 1.5])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="mean",
                dim=dim,
                skipna=skipna,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.mean,
                dim=dim,
                skipna=skipna,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def prod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``prod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``prod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``prod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``prod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.prod
        dask.array.prod
        DataArray.prod
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").prod()
        <xarray.DataArray (labels: 3)> Size: 24B
        array([1., 4., 0.])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.groupby("labels").prod(skipna=False)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([nan,  4.,  0.])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> da.groupby("labels").prod(skipna=True, min_count=2)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([nan,  4.,  0.])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="prod",
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.prod,
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def sum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``sum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``sum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``sum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``sum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.sum
        dask.array.sum
        DataArray.sum
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").sum()
        <xarray.DataArray (labels: 3)> Size: 24B
        array([1., 4., 3.])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.groupby("labels").sum(skipna=False)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([nan,  4.,  3.])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> da.groupby("labels").sum(skipna=True, min_count=2)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([nan,  4.,  3.])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="sum",
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.sum,
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def std(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``std`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``std``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``std`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``std`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.std
        dask.array.std
        DataArray.std
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").std()
        <xarray.DataArray (labels: 3)> Size: 24B
        array([0. , 0. , 1.5])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.groupby("labels").std(skipna=False)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([nan, 0. , 1.5])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Specify ``ddof=1`` for an unbiased estimate.

        >>> da.groupby("labels").std(skipna=True, ddof=1)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([       nan, 0.        , 2.12132034])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="std",
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.std,
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def var(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``var`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``var``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``var`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``var`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.var
        dask.array.var
        DataArray.var
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").var()
        <xarray.DataArray (labels: 3)> Size: 24B
        array([0.  , 0.  , 2.25])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.groupby("labels").var(skipna=False)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([ nan, 0.  , 2.25])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Specify ``ddof=1`` for an unbiased estimate.

        >>> da.groupby("labels").var(skipna=True, ddof=1)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([nan, 0. , 4.5])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="var",
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.var,
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def median(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``median`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``median``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``median`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``median`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.median
        dask.array.median
        DataArray.median
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").median()
        <xarray.DataArray (labels: 3)> Size: 24B
        array([1. , 2. , 1.5])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.groupby("labels").median(skipna=False)
        <xarray.DataArray (labels: 3)> Size: 24B
        array([nan, 2. , 1.5])
        Coordinates:
          * labels   (labels) object 24B 'a' 'b' 'c'
        """
        return self.reduce(
            duck_array_ops.median,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumsum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``cumsum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumsum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumsum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``cumsum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumsum
        dask.array.cumsum
        DataArray.cumsum
        DataArray.cumulative
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").cumsum()
        <xarray.DataArray (time: 6)> Size: 48B
        array([1., 2., 3., 3., 4., 1.])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.groupby("labels").cumsum(skipna=False)
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  3.,  4., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        """
        return self.reduce(
            duck_array_ops.cumsum,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumprod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``cumprod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumprod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the GroupBy dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumprod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``cumprod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumprod
        dask.array.cumprod
        DataArray.cumprod
        DataArray.cumulative
        :ref:`groupby`
            User guide on groupby operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up groupby computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.groupby("labels").cumprod()
        <xarray.DataArray (time: 6)> Size: 48B
        array([1., 2., 3., 0., 4., 1.])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.groupby("labels").cumprod(skipna=False)
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  4., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        """
        return self.reduce(
            duck_array_ops.cumprod,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )


class DataArrayResampleAggregations:
    _obj: DataArray

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keep_attrs: bool | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> DataArray:
        raise NotImplementedError()

    def _flox_reduce(
        self,
        dim: Dims,
        **kwargs: Any,
    ) -> DataArray:
        raise NotImplementedError()

    def count(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``count`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``count``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``count`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``count`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        pandas.DataFrame.count
        dask.dataframe.DataFrame.count
        DataArray.count
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").count()
        <xarray.DataArray (time: 3)> Size: 24B
        array([1, 3, 1])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="count",
                dim=dim,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.count,
                dim=dim,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def all(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``all`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``all``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``all`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``all`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.all
        dask.array.all
        DataArray.all
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 6B
        array([ True,  True,  True,  True,  True, False])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").all()
        <xarray.DataArray (time: 3)> Size: 3B
        array([ True,  True, False])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="all",
                dim=dim,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.array_all,
                dim=dim,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def any(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``any`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``any``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``any`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``any`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.any
        dask.array.any
        DataArray.any
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([True, True, True, True, True, False], dtype=bool),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 6B
        array([ True,  True,  True,  True,  True, False])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").any()
        <xarray.DataArray (time: 3)> Size: 3B
        array([ True,  True,  True])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="any",
                dim=dim,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.array_any,
                dim=dim,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def max(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``max`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``max``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``max`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``max`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.max
        dask.array.max
        DataArray.max
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").max()
        <xarray.DataArray (time: 3)> Size: 24B
        array([1., 3., 2.])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.resample(time="3ME").max(skipna=False)
        <xarray.DataArray (time: 3)> Size: 24B
        array([ 1.,  3., nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="max",
                dim=dim,
                skipna=skipna,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.max,
                dim=dim,
                skipna=skipna,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def min(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``min`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``min``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``min`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``min`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.min
        dask.array.min
        DataArray.min
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").min()
        <xarray.DataArray (time: 3)> Size: 24B
        array([1., 0., 2.])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.resample(time="3ME").min(skipna=False)
        <xarray.DataArray (time: 3)> Size: 24B
        array([ 1.,  0., nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="min",
                dim=dim,
                skipna=skipna,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.min,
                dim=dim,
                skipna=skipna,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def mean(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``mean`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``mean``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``mean`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``mean`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.mean
        dask.array.mean
        DataArray.mean
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").mean()
        <xarray.DataArray (time: 3)> Size: 24B
        array([1.        , 1.66666667, 2.        ])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.resample(time="3ME").mean(skipna=False)
        <xarray.DataArray (time: 3)> Size: 24B
        array([1.        , 1.66666667,        nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="mean",
                dim=dim,
                skipna=skipna,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.mean,
                dim=dim,
                skipna=skipna,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def prod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``prod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``prod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``prod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``prod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.prod
        dask.array.prod
        DataArray.prod
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").prod()
        <xarray.DataArray (time: 3)> Size: 24B
        array([1., 0., 2.])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.resample(time="3ME").prod(skipna=False)
        <xarray.DataArray (time: 3)> Size: 24B
        array([ 1.,  0., nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> da.resample(time="3ME").prod(skipna=True, min_count=2)
        <xarray.DataArray (time: 3)> Size: 24B
        array([nan,  0., nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="prod",
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.prod,
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def sum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``sum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``sum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        min_count : int or None, optional
            The required number of valid values to perform the operation. If
            fewer than min_count non-NA values are present the result will be
            NA. Only used if skipna is set to True or defaults to True for the
            array's dtype. Changed in version 0.17.0: if specified on an integer
            array and skipna=True, the result will be a float array.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``sum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``sum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.sum
        dask.array.sum
        DataArray.sum
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").sum()
        <xarray.DataArray (time: 3)> Size: 24B
        array([1., 5., 2.])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.resample(time="3ME").sum(skipna=False)
        <xarray.DataArray (time: 3)> Size: 24B
        array([ 1.,  5., nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> da.resample(time="3ME").sum(skipna=True, min_count=2)
        <xarray.DataArray (time: 3)> Size: 24B
        array([nan,  5., nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="sum",
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.sum,
                dim=dim,
                skipna=skipna,
                min_count=min_count,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def std(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``std`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``std``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``std`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``std`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.std
        dask.array.std
        DataArray.std
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").std()
        <xarray.DataArray (time: 3)> Size: 24B
        array([0.        , 1.24721913, 0.        ])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.resample(time="3ME").std(skipna=False)
        <xarray.DataArray (time: 3)> Size: 24B
        array([0.        , 1.24721913,        nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Specify ``ddof=1`` for an unbiased estimate.

        >>> da.resample(time="3ME").std(skipna=True, ddof=1)
        <xarray.DataArray (time: 3)> Size: 24B
        array([       nan, 1.52752523,        nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="std",
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.std,
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def var(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``var`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``var``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        ddof : int, default: 0
            “Delta Degrees of Freedom”: the divisor used in the calculation is ``N - ddof``,
            where ``N`` represents the number of elements.
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``var`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``var`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.var
        dask.array.var
        DataArray.var
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").var()
        <xarray.DataArray (time: 3)> Size: 24B
        array([0.        , 1.55555556, 0.        ])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.resample(time="3ME").var(skipna=False)
        <xarray.DataArray (time: 3)> Size: 24B
        array([0.        , 1.55555556,        nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Specify ``ddof=1`` for an unbiased estimate.

        >>> da.resample(time="3ME").var(skipna=True, ddof=1)
        <xarray.DataArray (time: 3)> Size: 24B
        array([       nan, 2.33333333,        nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        if (
            flox_available
            and OPTIONS["use_flox"]
            and contains_only_chunked_or_numpy(self._obj)
        ):
            return self._flox_reduce(
                func="var",
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                # fill_value=fill_value,
                keep_attrs=keep_attrs,
                **kwargs,
            )
        else:
            return self.reduce(
                duck_array_ops.var,
                dim=dim,
                skipna=skipna,
                ddof=ddof,
                keep_attrs=keep_attrs,
                **kwargs,
            )

    def median(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``median`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``median``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``median`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``median`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.median
        dask.array.median
        DataArray.median
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").median()
        <xarray.DataArray (time: 3)> Size: 24B
        array([1., 2., 2.])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.resample(time="3ME").median(skipna=False)
        <xarray.DataArray (time: 3)> Size: 24B
        array([ 1.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 24B 2001-01-31 2001-04-30 2001-07-31
        """
        return self.reduce(
            duck_array_ops.median,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumsum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``cumsum`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumsum``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumsum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``cumsum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumsum
        dask.array.cumsum
        DataArray.cumsum
        DataArray.cumulative
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").cumsum()
        <xarray.DataArray (time: 6)> Size: 48B
        array([1., 2., 5., 5., 2., 2.])
        Coordinates:
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Dimensions without coordinates: time

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.resample(time="3ME").cumsum(skipna=False)
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  5.,  5.,  2., nan])
        Coordinates:
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Dimensions without coordinates: time
        """
        return self.reduce(
            duck_array_ops.cumsum,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )

    def cumprod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """
        Reduce this DataArray's data by applying ``cumprod`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``cumprod``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If None, will reduce over the Resample dimensions.
            If "...", will reduce over all dimensions.
        skipna : bool or None, optional
            If True, skip missing values (as marked by NaN). By default, only
            skips missing values for float dtypes; other dtypes either do not
            have a sentinel missing value (int) or ``skipna=True`` has not been
            implemented (object, datetime64 or timedelta64).
        keep_attrs : bool or None, optional
            If True, ``attrs`` will be copied from the original
            object to the new one.  If False, the new object will be
            returned without attributes.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumprod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : DataArray
            New DataArray with ``cumprod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumprod
        dask.array.cumprod
        DataArray.cumprod
        DataArray.cumulative
        :ref:`resampling`
            User guide on resampling operations.

        Notes
        -----
        Use the ``flox`` package to significantly speed up resampling computations,
        especially with dask arrays. Xarray will use flox by default if installed.
        Pass flox-specific keyword arguments in ``**kwargs``.
        See the `flox documentation <https://flox.readthedocs.io>`_ for more.

        Non-numeric variables will be removed prior to reducing.

        Note that the methods on the ``cumulative`` method are more performant (with numbagg installed)
        and better supported. ``cumsum`` and ``cumprod`` may be deprecated
        in the future.

        Examples
        --------
        >>> da = xr.DataArray(
        ...     np.array([1, 2, 3, 0, 2, np.nan]),
        ...     dims="time",
        ...     coords=dict(
        ...         time=("time", pd.date_range("2001-01-01", freq="ME", periods=6)),
        ...         labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
        ...     ),
        ... )
        >>> da
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])
        Coordinates:
          * time     (time) datetime64[ns] 48B 2001-01-31 2001-02-28 ... 2001-06-30
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'

        >>> da.resample(time="3ME").cumprod()
        <xarray.DataArray (time: 6)> Size: 48B
        array([1., 2., 6., 0., 2., 2.])
        Coordinates:
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Dimensions without coordinates: time

        Use ``skipna`` to control whether NaNs are ignored.

        >>> da.resample(time="3ME").cumprod(skipna=False)
        <xarray.DataArray (time: 6)> Size: 48B
        array([ 1.,  2.,  6.,  0.,  2., nan])
        Coordinates:
            labels   (time) <U1 24B 'a' 'b' 'c' 'c' 'b' 'a'
        Dimensions without coordinates: time
        """
        return self.reduce(
            duck_array_ops.cumprod,
            dim=dim,
            skipna=skipna,
            keep_attrs=keep_attrs,
            **kwargs,
        )
