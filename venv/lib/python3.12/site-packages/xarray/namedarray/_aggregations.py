"""Mixin classes with reduction operations."""

# This file was generated using xarray.util.generate_aggregations. Do not edit manually.

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

from xarray.core import duck_array_ops
from xarray.core.types import Dims, Self


class NamedArrayAggregations:
    __slots__ = ()

    def reduce(
        self,
        func: Callable[..., Any],
        dim: Dims = None,
        *,
        axis: int | Sequence[int] | None = None,
        keepdims: bool = False,
        **kwargs: Any,
    ) -> Self:
        raise NotImplementedError()

    def count(
        self,
        dim: Dims = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``count`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``count``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``count`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``count`` applied to its data and the
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.count()
        <xarray.NamedArray ()> Size: 8B
        array(5)
        """
        return self.reduce(
            duck_array_ops.count,
            dim=dim,
            **kwargs,
        )

    def all(
        self,
        dim: Dims = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``all`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``all``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``all`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``all`` applied to its data and the
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray(
        ...     "x", np.array([True, True, True, True, True, False], dtype=bool)
        ... )
        >>> na
        <xarray.NamedArray (x: 6)> Size: 6B
        array([ True,  True,  True,  True,  True, False])

        >>> na.all()
        <xarray.NamedArray ()> Size: 1B
        array(False)
        """
        return self.reduce(
            duck_array_ops.array_all,
            dim=dim,
            **kwargs,
        )

    def any(
        self,
        dim: Dims = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``any`` along some dimension(s).

        Parameters
        ----------
        dim : str, Iterable of Hashable, "..." or None, default: None
            Name of dimension[s] along which to apply ``any``. For e.g. ``dim="x"``
            or ``dim=["x", "y"]``. If "..." or None, will reduce over all dimensions.
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``any`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``any`` applied to its data and the
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray(
        ...     "x", np.array([True, True, True, True, True, False], dtype=bool)
        ... )
        >>> na
        <xarray.NamedArray (x: 6)> Size: 6B
        array([ True,  True,  True,  True,  True, False])

        >>> na.any()
        <xarray.NamedArray ()> Size: 1B
        array(True)
        """
        return self.reduce(
            duck_array_ops.array_any,
            dim=dim,
            **kwargs,
        )

    def max(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``max`` along some dimension(s).

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
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``max`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``max`` applied to its data and the
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.max()
        <xarray.NamedArray ()> Size: 8B
        array(3.)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> na.max(skipna=False)
        <xarray.NamedArray ()> Size: 8B
        array(nan)
        """
        return self.reduce(
            duck_array_ops.max,
            dim=dim,
            skipna=skipna,
            **kwargs,
        )

    def min(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``min`` along some dimension(s).

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
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``min`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``min`` applied to its data and the
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.min()
        <xarray.NamedArray ()> Size: 8B
        array(0.)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> na.min(skipna=False)
        <xarray.NamedArray ()> Size: 8B
        array(nan)
        """
        return self.reduce(
            duck_array_ops.min,
            dim=dim,
            skipna=skipna,
            **kwargs,
        )

    def mean(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``mean`` along some dimension(s).

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
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``mean`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``mean`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.mean
        dask.array.mean
        Dataset.mean
        DataArray.mean
        :ref:`agg`
            User guide on reduction or aggregation operations.

        Examples
        --------
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.mean()
        <xarray.NamedArray ()> Size: 8B
        array(1.6)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> na.mean(skipna=False)
        <xarray.NamedArray ()> Size: 8B
        array(nan)
        """
        return self.reduce(
            duck_array_ops.mean,
            dim=dim,
            skipna=skipna,
            **kwargs,
        )

    def prod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``prod`` along some dimension(s).

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
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``prod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``prod`` applied to its data and the
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.prod()
        <xarray.NamedArray ()> Size: 8B
        array(0.)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> na.prod(skipna=False)
        <xarray.NamedArray ()> Size: 8B
        array(nan)

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> na.prod(skipna=True, min_count=2)
        <xarray.NamedArray ()> Size: 8B
        array(0.)
        """
        return self.reduce(
            duck_array_ops.prod,
            dim=dim,
            skipna=skipna,
            min_count=min_count,
            **kwargs,
        )

    def sum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        min_count: int | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``sum`` along some dimension(s).

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
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``sum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``sum`` applied to its data and the
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.sum()
        <xarray.NamedArray ()> Size: 8B
        array(8.)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> na.sum(skipna=False)
        <xarray.NamedArray ()> Size: 8B
        array(nan)

        Specify ``min_count`` for finer control over when NaNs are ignored.

        >>> na.sum(skipna=True, min_count=2)
        <xarray.NamedArray ()> Size: 8B
        array(8.)
        """
        return self.reduce(
            duck_array_ops.sum,
            dim=dim,
            skipna=skipna,
            min_count=min_count,
            **kwargs,
        )

    def std(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``std`` along some dimension(s).

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
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``std`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``std`` applied to its data and the
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.std()
        <xarray.NamedArray ()> Size: 8B
        array(1.0198039)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> na.std(skipna=False)
        <xarray.NamedArray ()> Size: 8B
        array(nan)

        Specify ``ddof=1`` for an unbiased estimate.

        >>> na.std(skipna=True, ddof=1)
        <xarray.NamedArray ()> Size: 8B
        array(1.14017543)
        """
        return self.reduce(
            duck_array_ops.std,
            dim=dim,
            skipna=skipna,
            ddof=ddof,
            **kwargs,
        )

    def var(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        ddof: int = 0,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``var`` along some dimension(s).

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
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``var`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``var`` applied to its data and the
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.var()
        <xarray.NamedArray ()> Size: 8B
        array(1.04)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> na.var(skipna=False)
        <xarray.NamedArray ()> Size: 8B
        array(nan)

        Specify ``ddof=1`` for an unbiased estimate.

        >>> na.var(skipna=True, ddof=1)
        <xarray.NamedArray ()> Size: 8B
        array(1.3)
        """
        return self.reduce(
            duck_array_ops.var,
            dim=dim,
            skipna=skipna,
            ddof=ddof,
            **kwargs,
        )

    def median(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``median`` along some dimension(s).

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
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``median`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``median`` applied to its data and the
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.median()
        <xarray.NamedArray ()> Size: 8B
        array(2.)

        Use ``skipna`` to control whether NaNs are ignored.

        >>> na.median(skipna=False)
        <xarray.NamedArray ()> Size: 8B
        array(nan)
        """
        return self.reduce(
            duck_array_ops.median,
            dim=dim,
            skipna=skipna,
            **kwargs,
        )

    def cumsum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``cumsum`` along some dimension(s).

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
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumsum`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``cumsum`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumsum
        dask.array.cumsum
        Dataset.cumsum
        DataArray.cumsum
        NamedArray.cumulative
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.cumsum()
        <xarray.NamedArray (x: 6)> Size: 48B
        array([1., 3., 6., 6., 8., 8.])

        Use ``skipna`` to control whether NaNs are ignored.

        >>> na.cumsum(skipna=False)
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  3.,  6.,  6.,  8., nan])
        """
        return self.reduce(
            duck_array_ops.cumsum,
            dim=dim,
            skipna=skipna,
            **kwargs,
        )

    def cumprod(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        **kwargs: Any,
    ) -> Self:
        """
        Reduce this NamedArray's data by applying ``cumprod`` along some dimension(s).

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
        **kwargs : Any
            Additional keyword arguments passed on to the appropriate array
            function for calculating ``cumprod`` on this object's data.
            These could include dask-specific kwargs like ``split_every``.

        Returns
        -------
        reduced : NamedArray
            New NamedArray with ``cumprod`` applied to its data and the
            indicated dimension(s) removed

        See Also
        --------
        numpy.cumprod
        dask.array.cumprod
        Dataset.cumprod
        DataArray.cumprod
        NamedArray.cumulative
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
        >>> from xarray.namedarray.core import NamedArray
        >>> na = NamedArray("x", np.array([1, 2, 3, 0, 2, np.nan]))
        >>> na
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  3.,  0.,  2., nan])

        >>> na.cumprod()
        <xarray.NamedArray (x: 6)> Size: 48B
        array([1., 2., 6., 0., 0., 0.])

        Use ``skipna`` to control whether NaNs are ignored.

        >>> na.cumprod(skipna=False)
        <xarray.NamedArray (x: 6)> Size: 48B
        array([ 1.,  2.,  6.,  0.,  0., nan])
        """
        return self.reduce(
            duck_array_ops.cumprod,
            dim=dim,
            skipna=skipna,
            **kwargs,
        )
