from __future__ import annotations

from collections.abc import Hashable, Iterable, Sequence
from typing import TYPE_CHECKING, Generic, Literal, cast

import numpy as np
from numpy.typing import ArrayLike

from xarray.computation.apply_ufunc import apply_ufunc
from xarray.computation.computation import dot
from xarray.core import duck_array_ops, utils
from xarray.core.types import Dims, T_DataArray, T_Xarray
from xarray.namedarray.utils import is_duck_dask_array
from xarray.structure.alignment import align, broadcast

# Weighted quantile methods are a subset of the numpy supported quantile methods.
QUANTILE_METHODS = Literal[
    "linear",
    "interpolated_inverted_cdf",
    "hazen",
    "weibull",
    "median_unbiased",
    "normal_unbiased",
]

_WEIGHTED_REDUCE_DOCSTRING_TEMPLATE = """
    Reduce this {cls}'s data by a weighted ``{fcn}`` along some dimension(s).

    Parameters
    ----------
    dim : Hashable or Iterable of Hashable, optional
        Dimension(s) over which to apply the weighted ``{fcn}``.
    skipna : bool or None, optional
        If True, skip missing values (as marked by NaN). By default, only
        skips missing values for float dtypes; other dtypes either do not
        have a sentinel missing value (int) or skipna=True has not been
        implemented (object, datetime64 or timedelta64).
    keep_attrs : bool or None, optional
        If True, the attributes (``attrs``) will be copied from the original
        object to the new one.  If False (default), the new object will be
        returned without attributes.

    Returns
    -------
    reduced : {cls}
        New {cls} object with weighted ``{fcn}`` applied to its data and
        the indicated dimension(s) removed.

    Notes
    -----
        Returns {on_zero} if the ``weights`` sum to 0.0 along the reduced
        dimension(s).
    """

_SUM_OF_WEIGHTS_DOCSTRING = """
    Calculate the sum of weights, accounting for missing values in the data.

    Parameters
    ----------
    dim : str or sequence of str, optional
        Dimension(s) over which to sum the weights.
    keep_attrs : bool, optional
        If True, the attributes (``attrs``) will be copied from the original
        object to the new one.  If False (default), the new object will be
        returned without attributes.

    Returns
    -------
    reduced : {cls}
        New {cls} object with the sum of the weights over the given dimension.
    """

_WEIGHTED_QUANTILE_DOCSTRING_TEMPLATE = """
    Apply a weighted ``quantile`` to this {cls}'s data along some dimension(s).

    Weights are interpreted as *sampling weights* (or probability weights) and
    describe how a sample is scaled to the whole population [1]_. There are
    other possible interpretations for weights, *precision weights* describing the
    precision of observations, or *frequency weights* counting the number of identical
    observations, however, they are not implemented here.

    For compatibility with NumPy's non-weighted ``quantile`` (which is used by
    ``DataArray.quantile`` and ``Dataset.quantile``), the only interpolation
    method supported by this weighted version corresponds to the default "linear"
    option of ``numpy.quantile``. This is "Type 7" option, described in Hyndman
    and Fan (1996) [2]_. The implementation is largely inspired by a blog post
    from A. Akinshin's (2023) [3]_.

    Parameters
    ----------
    q : float or sequence of float
        Quantile to compute, which must be between 0 and 1 inclusive.
    dim : str or sequence of str, optional
        Dimension(s) over which to apply the weighted ``quantile``.
    skipna : bool, optional
        If True, skip missing values (as marked by NaN). By default, only
        skips missing values for float dtypes; other dtypes either do not
        have a sentinel missing value (int) or skipna=True has not been
        implemented (object, datetime64 or timedelta64).
    keep_attrs : bool, optional
        If True, the attributes (``attrs``) will be copied from the original
        object to the new one.  If False (default), the new object will be
        returned without attributes.

    Returns
    -------
    quantiles : {cls}
        New {cls} object with weighted ``quantile`` applied to its data and
        the indicated dimension(s) removed.

    See Also
    --------
    numpy.nanquantile, pandas.Series.quantile, Dataset.quantile, DataArray.quantile

    Notes
    -----
    Returns NaN if the ``weights`` sum to 0.0 along the reduced
    dimension(s).

    References
    ----------
    .. [1] https://notstatschat.rbind.io/2020/08/04/weights-in-statistics/
    .. [2] Hyndman, R. J. & Fan, Y. (1996). Sample Quantiles in Statistical Packages.
           The American Statistician, 50(4), 361–365. https://doi.org/10.2307/2684934
    .. [3] Akinshin, A. (2023) "Weighted quantile estimators" arXiv:2304.07265 [stat.ME]
           https://arxiv.org/abs/2304.07265
    """


if TYPE_CHECKING:
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset


class Weighted(Generic[T_Xarray]):
    """An object that implements weighted operations.

    You should create a Weighted object by using the ``DataArray.weighted`` or
    ``Dataset.weighted`` methods.

    See Also
    --------
    Dataset.weighted
    DataArray.weighted
    """

    __slots__ = ("obj", "weights")

    def __init__(self, obj: T_Xarray, weights: T_DataArray) -> None:
        """
        Create a Weighted object

        Parameters
        ----------
        obj : DataArray or Dataset
            Object over which the weighted reduction operation is applied.
        weights : DataArray
            An array of weights associated with the values in the obj.
            Each value in the obj contributes to the reduction operation
            according to its associated weight.

        Notes
        -----
        ``weights`` must be a ``DataArray`` and cannot contain missing values.
        Missing values can be replaced by ``weights.fillna(0)``.
        """

        from xarray.core.dataarray import DataArray

        if not isinstance(weights, DataArray):
            raise ValueError("`weights` must be a DataArray")

        def _weight_check(w):
            # Ref https://github.com/pydata/xarray/pull/4559/files#r515968670
            if duck_array_ops.array_any(duck_array_ops.isnull(w)):
                raise ValueError(
                    "`weights` cannot contain missing values. "
                    "Missing values can be replaced by `weights.fillna(0)`."
                )
            return w

        if is_duck_dask_array(weights.data):
            # assign to copy - else the check is not triggered
            weights = weights.copy(
                data=weights.data.map_blocks(_weight_check, dtype=weights.dtype),  # type: ignore[call-arg, arg-type]
                deep=False,
            )

        else:
            _weight_check(weights.data)

        self.obj: T_Xarray = obj
        self.weights: T_DataArray = weights

    def _check_dim(self, dim: Dims):
        """raise an error if any dimension is missing"""

        dims: list[Hashable]
        if isinstance(dim, str) or not isinstance(dim, Iterable):
            dims = [dim] if dim else []
        else:
            dims = list(dim)
        all_dims = set(self.obj.dims).union(set(self.weights.dims))
        missing_dims = set(dims) - all_dims
        if missing_dims:
            raise ValueError(
                f"Dimensions {tuple(missing_dims)} not found in {self.__class__.__name__} dimensions {tuple(all_dims)}"
            )

    @staticmethod
    def _reduce(
        da: T_DataArray,
        weights: T_DataArray,
        dim: Dims = None,
        skipna: bool | None = None,
    ) -> T_DataArray:
        """reduce using dot; equivalent to (da * weights).sum(dim, skipna)

        for internal use only
        """

        # need to infer dims as we use `dot`
        if dim is None:
            dim = ...

        # need to mask invalid values in da, as `dot` does not implement skipna
        if skipna or (skipna is None and da.dtype.kind in "cfO"):
            da = da.fillna(0.0)

        # `dot` does not broadcast arrays, so this avoids creating a large
        # DataArray (if `weights` has additional dimensions)
        return dot(da, weights, dim=dim)

    def _sum_of_weights(self, da: T_DataArray, dim: Dims = None) -> T_DataArray:
        """Calculate the sum of weights, accounting for missing values"""

        # we need to mask data values that are nan; else the weights are wrong
        mask = da.notnull()

        # bool -> int, because ``xr.dot([True, True], [True, True])`` -> True
        # (and not 2); GH4074
        if self.weights.dtype == bool:
            sum_of_weights = self._reduce(
                mask,
                duck_array_ops.astype(self.weights, dtype=int),
                dim=dim,
                skipna=False,
            )
        else:
            sum_of_weights = self._reduce(mask, self.weights, dim=dim, skipna=False)

        # 0-weights are not valid
        valid_weights = sum_of_weights != 0.0

        return sum_of_weights.where(valid_weights)

    def _sum_of_squares(
        self,
        da: T_DataArray,
        dim: Dims = None,
        skipna: bool | None = None,
    ) -> T_DataArray:
        """Reduce a DataArray by a weighted ``sum_of_squares`` along some dimension(s)."""

        demeaned = da - da.weighted(self.weights).mean(dim=dim)

        # TODO: unsure why mypy complains about these being DataArray return types
        # rather than T_DataArray?
        return self._reduce((demeaned**2), self.weights, dim=dim, skipna=skipna)  # type: ignore[return-value]

    def _weighted_sum(
        self,
        da: T_DataArray,
        dim: Dims = None,
        skipna: bool | None = None,
    ) -> T_DataArray:
        """Reduce a DataArray by a weighted ``sum`` along some dimension(s)."""

        return self._reduce(da, self.weights, dim=dim, skipna=skipna)  # type: ignore[return-value]

    def _weighted_mean(
        self,
        da: T_DataArray,
        dim: Dims = None,
        skipna: bool | None = None,
    ) -> T_DataArray:
        """Reduce a DataArray by a weighted ``mean`` along some dimension(s)."""

        weighted_sum = self._weighted_sum(da, dim=dim, skipna=skipna)

        sum_of_weights = self._sum_of_weights(da, dim=dim)

        return weighted_sum / sum_of_weights

    def _weighted_var(
        self,
        da: T_DataArray,
        dim: Dims = None,
        skipna: bool | None = None,
    ) -> T_DataArray:
        """Reduce a DataArray by a weighted ``var`` along some dimension(s)."""

        sum_of_squares = self._sum_of_squares(da, dim=dim, skipna=skipna)

        sum_of_weights = self._sum_of_weights(da, dim=dim)

        return sum_of_squares / sum_of_weights

    def _weighted_std(
        self,
        da: T_DataArray,
        dim: Dims = None,
        skipna: bool | None = None,
    ) -> T_DataArray:
        """Reduce a DataArray by a weighted ``std`` along some dimension(s)."""

        return cast("T_DataArray", np.sqrt(self._weighted_var(da, dim, skipna)))

    def _weighted_quantile(
        self,
        da: T_DataArray,
        q: ArrayLike,
        dim: Dims = None,
        skipna: bool | None = None,
    ) -> T_DataArray:
        """Apply a weighted ``quantile`` to a DataArray along some dimension(s)."""

        def _get_h(n: float, q: np.ndarray, method: QUANTILE_METHODS) -> np.ndarray:
            """Return the interpolation parameter."""
            # Note that options are not yet exposed in the public API.
            h: np.ndarray
            if method == "linear":
                h = (n - 1) * q + 1
            elif method == "interpolated_inverted_cdf":
                h = n * q
            elif method == "hazen":
                h = n * q + 0.5
            elif method == "weibull":
                h = (n + 1) * q
            elif method == "median_unbiased":
                h = (n + 1 / 3) * q + 1 / 3
            elif method == "normal_unbiased":
                h = (n + 1 / 4) * q + 3 / 8
            else:
                raise ValueError(f"Invalid method: {method}.")
            return h.clip(1, n)

        def _weighted_quantile_1d(
            data: np.ndarray,
            weights: np.ndarray,
            q: np.ndarray,
            skipna: bool,
            method: QUANTILE_METHODS = "linear",
        ) -> np.ndarray:
            # This algorithm has been adapted from:
            #   https://aakinshin.net/posts/weighted-quantiles/#reference-implementation
            is_nan = np.isnan(data)
            if skipna:
                # Remove nans from data and weights
                not_nan = ~is_nan
                data = data[not_nan]
                weights = weights[not_nan]
            elif is_nan.any():
                # Return nan if data contains any nan
                return np.full(q.size, np.nan)

            # Filter out data (and weights) associated with zero weights, which also flattens them
            nonzero_weights = weights != 0
            data = data[nonzero_weights]
            weights = weights[nonzero_weights]
            n = data.size

            if n == 0:
                # Possibly empty after nan or zero weight filtering above
                return np.full(q.size, np.nan)

            # Kish's effective sample size
            nw = weights.sum() ** 2 / (weights**2).sum()

            # Sort data and weights
            sorter = np.argsort(data)
            data = data[sorter]
            weights = weights[sorter]

            # Normalize and sum the weights
            weights = weights / weights.sum()
            weights_cum = np.append(0, weights.cumsum())

            # Vectorize the computation by transposing q with respect to weights
            q = np.atleast_2d(q).T

            # Get the interpolation parameter for each q
            h = _get_h(nw, q, method)

            # Find the samples contributing to the quantile computation (at *positions* between (h-1)/nw and h/nw)
            u = np.maximum((h - 1) / nw, np.minimum(h / nw, weights_cum))

            # Compute their relative weight
            v = u * nw - h + 1
            w = np.diff(v)

            # Apply the weights
            return (data * w).sum(axis=1)

        if skipna is None and da.dtype.kind in "cfO":
            skipna = True

        q = np.atleast_1d(np.asarray(q, dtype=np.float64))

        if q.ndim > 1:
            raise ValueError("q must be a scalar or 1d")

        if np.any((q < 0) | (q > 1)):
            raise ValueError("q values must be between 0 and 1")

        if dim is None:
            dim = da.dims

        if utils.is_scalar(dim):
            dim = [dim]

        # To satisfy mypy
        dim = cast(Sequence, dim)

        # need to align *and* broadcast
        # - `_weighted_quantile_1d` requires arrays with the same shape
        # - broadcast does an outer join, which can introduce NaN to weights
        # - therefore we first need to do align(..., join="inner")

        # TODO: use broadcast(..., join="inner") once available
        # see https://github.com/pydata/xarray/issues/6304

        da, weights = align(da, self.weights, join="inner")
        da, weights = broadcast(da, weights)

        result = apply_ufunc(
            _weighted_quantile_1d,
            da,
            weights,
            input_core_dims=[dim, dim],
            output_core_dims=[["quantile"]],
            output_dtypes=[np.float64],
            dask_gufunc_kwargs=dict(output_sizes={"quantile": len(q)}),
            dask="parallelized",
            vectorize=True,
            kwargs={"q": q, "skipna": skipna},
        )

        result = result.transpose("quantile", ...)
        result = result.assign_coords(quantile=q).squeeze()
        return result

    def _implementation(self, func, dim, **kwargs):
        raise NotImplementedError("Use `Dataset.weighted` or `DataArray.weighted`")

    def sum_of_weights(
        self,
        dim: Dims = None,
        *,
        keep_attrs: bool | None = None,
    ) -> T_Xarray:
        return self._implementation(
            self._sum_of_weights, dim=dim, keep_attrs=keep_attrs
        )

    def sum_of_squares(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
    ) -> T_Xarray:
        return self._implementation(
            self._sum_of_squares, dim=dim, skipna=skipna, keep_attrs=keep_attrs
        )

    def sum(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
    ) -> T_Xarray:
        return self._implementation(
            self._weighted_sum, dim=dim, skipna=skipna, keep_attrs=keep_attrs
        )

    def mean(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
    ) -> T_Xarray:
        return self._implementation(
            self._weighted_mean, dim=dim, skipna=skipna, keep_attrs=keep_attrs
        )

    def var(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
    ) -> T_Xarray:
        return self._implementation(
            self._weighted_var, dim=dim, skipna=skipna, keep_attrs=keep_attrs
        )

    def std(
        self,
        dim: Dims = None,
        *,
        skipna: bool | None = None,
        keep_attrs: bool | None = None,
    ) -> T_Xarray:
        return self._implementation(
            self._weighted_std, dim=dim, skipna=skipna, keep_attrs=keep_attrs
        )

    def quantile(
        self,
        q: ArrayLike,
        *,
        dim: Dims = None,
        keep_attrs: bool | None = None,
        skipna: bool = True,
    ) -> T_Xarray:
        return self._implementation(
            self._weighted_quantile, q=q, dim=dim, skipna=skipna, keep_attrs=keep_attrs
        )

    def __repr__(self) -> str:
        """provide a nice str repr of our Weighted object"""

        klass = self.__class__.__name__
        weight_dims = ", ".join(map(str, self.weights.dims))
        return f"{klass} with weights along dimensions: {weight_dims}"


class DataArrayWeighted(Weighted["DataArray"]):
    def _implementation(self, func, dim, **kwargs) -> DataArray:
        self._check_dim(dim)

        dataset = self.obj._to_temp_dataset()
        dataset = dataset.map(func, dim=dim, **kwargs)
        return self.obj._from_temp_dataset(dataset)


class DatasetWeighted(Weighted["Dataset"]):
    def _implementation(self, func, dim, **kwargs) -> Dataset:
        self._check_dim(dim)
        return self.obj.map(func, dim=dim, **kwargs)


def _inject_docstring(cls, cls_name):
    cls.sum_of_weights.__doc__ = _SUM_OF_WEIGHTS_DOCSTRING.format(cls=cls_name)

    cls.sum.__doc__ = _WEIGHTED_REDUCE_DOCSTRING_TEMPLATE.format(
        cls=cls_name, fcn="sum", on_zero="0"
    )

    cls.mean.__doc__ = _WEIGHTED_REDUCE_DOCSTRING_TEMPLATE.format(
        cls=cls_name, fcn="mean", on_zero="NaN"
    )

    cls.sum_of_squares.__doc__ = _WEIGHTED_REDUCE_DOCSTRING_TEMPLATE.format(
        cls=cls_name, fcn="sum_of_squares", on_zero="0"
    )

    cls.var.__doc__ = _WEIGHTED_REDUCE_DOCSTRING_TEMPLATE.format(
        cls=cls_name, fcn="var", on_zero="NaN"
    )

    cls.std.__doc__ = _WEIGHTED_REDUCE_DOCSTRING_TEMPLATE.format(
        cls=cls_name, fcn="std", on_zero="NaN"
    )

    cls.quantile.__doc__ = _WEIGHTED_QUANTILE_DOCSTRING_TEMPLATE.format(cls=cls_name)


_inject_docstring(DataArrayWeighted, "DataArray")
_inject_docstring(DatasetWeighted, "Dataset")
