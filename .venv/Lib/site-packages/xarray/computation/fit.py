"""Fitting operations for DataArrays and Datasets."""

from __future__ import annotations

import inspect
import warnings
from collections.abc import Callable, Hashable, Iterable, Mapping, Sequence
from inspect import Parameter
from types import MappingProxyType
from typing import (
    Any,
    Literal,
    Union,
)

import numpy as np

# remove once numpy 2.0 is the oldest supported version
try:
    from numpy.exceptions import RankWarning
except ImportError:
    from numpy import RankWarning  # type: ignore[no-redef,attr-defined,unused-ignore]

from xarray.computation.apply_ufunc import apply_ufunc
from xarray.computation.computation import _ensure_numeric, where
from xarray.core.dataarray import DataArray
from xarray.core.duck_array_ops import is_duck_dask_array, least_squares
from xarray.core.types import Dims, ErrorOptions
from xarray.core.variable import Variable
from xarray.structure.alignment import broadcast


def _get_func_args(func, param_names):
    """Use `inspect.signature` to try accessing `func` args. Otherwise, ensure
    they are provided by user.
    """
    func_args: Union[dict[str, Parameter], MappingProxyType[str, Parameter]]
    try:
        func_args = inspect.signature(func).parameters
    except ValueError as err:
        func_args = {}  # type: ignore[assignment,unused-ignore]
        if not param_names:
            raise ValueError(
                "Unable to inspect `func` signature, and `param_names` was not provided."
            ) from err
    if param_names:
        params = param_names
    else:
        params = list(func_args)[1:]
        if any(
            (p.kind in [p.VAR_POSITIONAL, p.VAR_KEYWORD]) for p in func_args.values()
        ):
            raise ValueError(
                "`param_names` must be provided because `func` takes variable length arguments."
            )
    return params, func_args


def _initialize_curvefit_params(params, p0, bounds, func_args):
    """Set initial guess and bounds for curvefit.
    Priority: 1) passed args 2) func signature 3) scipy defaults
    """

    def _initialize_feasible(lb, ub):
        # Mimics functionality of scipy.optimize.minpack._initialize_feasible
        lb_finite = np.isfinite(lb)
        ub_finite = np.isfinite(ub)
        p0 = where(
            lb_finite,
            where(
                ub_finite,
                0.5 * (lb + ub),  # both bounds finite
                lb + 1,  # lower bound finite, upper infinite
            ),
            where(
                ub_finite,
                ub - 1,  # lower bound infinite, upper finite
                0,  # both bounds infinite
            ),
        )
        return p0

    param_defaults = dict.fromkeys(params, 1)
    bounds_defaults = dict.fromkeys(params, (-np.inf, np.inf))
    for p in params:
        if p in func_args and func_args[p].default is not func_args[p].empty:
            param_defaults[p] = func_args[p].default
        if p in bounds:
            lb, ub = bounds[p]
            bounds_defaults[p] = (lb, ub)
            param_defaults[p] = where(
                (param_defaults[p] < lb) | (param_defaults[p] > ub),
                _initialize_feasible(lb, ub),
                param_defaults[p],
            )
        if p in p0:
            param_defaults[p] = p0[p]
    return param_defaults, bounds_defaults


def polyfit(
    obj,
    dim: Hashable,
    deg: int,
    skipna: bool | None = None,
    rcond: np.floating[Any] | float | None = None,
    w: Hashable | Any = None,
    full: bool = False,
    cov: bool | Literal["unscaled"] = False,
):
    """
    Least squares polynomial fit.

    This replicates the behaviour of `numpy.polyfit` but differs by skipping
    invalid values when `skipna = True`.

    Parameters
    ----------
    obj : Dataset or DataArray
        Object to perform the polyfit on
    dim : hashable
        Coordinate along which to fit the polynomials.
    deg : int
        Degree of the fitting polynomial.
    skipna : bool or None, optional
        If True, removes all invalid values before fitting each 1D slices of the array.
        Default is True if data is stored in a dask.array or if there is any
        invalid values, False otherwise.
    rcond : float or None, optional
        Relative condition number to the fit.
    w : hashable or Any, optional
        Weights to apply to the y-coordinate of the sample points.
        Can be an array-like object or the name of a coordinate in the dataset.
    full : bool, default: False
        Whether to return the residuals, matrix rank and singular values in addition
        to the coefficients.
    cov : bool or "unscaled", default: False
        Whether to return to the covariance matrix in addition to the coefficients.
        The matrix is not scaled if `cov='unscaled'`.

    Returns
    -------
    Dataset
        A single dataset which contains (for each "var" in the input dataset):

        [var]_polyfit_coefficients
            The coefficients of the best fit for each variable in this dataset.
        [var]_polyfit_residuals
            The residuals of the least-square computation for each variable (only included if `full=True`)
            When the matrix rank is deficient, np.nan is returned.
        [dim]_matrix_rank
            The effective rank of the scaled Vandermonde coefficient matrix (only included if `full=True`)
            The rank is computed ignoring the NaN values that might be skipped.
        [dim]_singular_values
            The singular values of the scaled Vandermonde coefficient matrix (only included if `full=True`)
        [var]_polyfit_covariance
            The covariance matrix of the polynomial coefficient estimates (only included if `full=False` and `cov=True`)

    Warns
    -----
    RankWarning
        The rank of the coefficient matrix in the least-squares fit is deficient.
        The warning is not raised with in-memory (not dask) data and `full=True`.

    See Also
    --------
    numpy.polyfit
    numpy.polyval
    xarray.polyval
    """
    variables: dict[Hashable, Variable] = {}
    skipna_da = skipna

    x = np.asarray(_ensure_numeric(obj.coords[dim]).astype(np.float64))

    xname = f"{obj[dim].name}_"
    order = int(deg) + 1
    degree_coord_values = np.arange(order)[::-1]
    lhs = np.vander(x, order)

    if rcond is None:
        rcond = x.shape[0] * np.finfo(x.dtype).eps

    # Weights:
    if w is not None:
        if isinstance(w, Hashable):
            w = obj.coords[w]
        w = np.asarray(w)
        if w.ndim != 1:
            raise TypeError("Expected a 1-d array for weights.")
        if w.shape[0] != lhs.shape[0]:
            raise TypeError(f"Expected w and {dim} to have the same length")
        lhs *= w[:, np.newaxis]

    # Scaling
    scale = np.sqrt((lhs * lhs).sum(axis=0))
    lhs /= scale

    from xarray.core import utils

    degree_dim = utils.get_temp_dimname(obj.dims, "degree")

    rank = np.linalg.matrix_rank(lhs)

    if full:
        rank = Variable(dims=(), data=rank)
        variables[xname + "matrix_rank"] = rank
        _sing = np.linalg.svd(lhs, compute_uv=False)
        variables[xname + "singular_values"] = Variable(
            dims=(degree_dim,),
            data=np.concatenate([np.full((order - rank.data,), np.nan), _sing]),
        )

    # If we have a coordinate get its underlying dimension.
    (true_dim,) = obj.coords[dim].dims

    other_coords = {
        dim: obj._variables[dim]
        for dim in set(obj.dims) - {true_dim}
        if dim in obj._variables
    }
    present_dims: set[Hashable] = set()
    for name, var in obj._variables.items():
        if name in obj._coord_names or name in obj.dims:
            continue
        if true_dim not in var.dims:
            continue

        if is_duck_dask_array(var._data) and (rank != order or full or skipna is None):
            # Current algorithm with dask and skipna=False neither supports
            # deficient ranks nor does it output the "full" info (issue dask/dask#6516)
            skipna_da = True
        elif skipna is None:
            skipna_da = bool(np.any(var.isnull()))

        if var.ndim > 1:
            rhs = var.transpose(true_dim, ...)
            other_dims = rhs.dims[1:]
            scale_da = scale.reshape(-1, *((1,) * len(other_dims)))
        else:
            rhs = var
            scale_da = scale
            other_dims = ()

        present_dims.update(other_dims)
        if w is not None:
            rhs = rhs * w.reshape(-1, *((1,) * len(other_dims)))

        with warnings.catch_warnings():
            if full:  # Copy np.polyfit behavior
                warnings.simplefilter("ignore", RankWarning)
            else:  # Raise only once per variable
                warnings.simplefilter("once", RankWarning)

            coeffs, residuals = least_squares(
                lhs, rhs.data, rcond=rcond, skipna=skipna_da
            )

        from xarray.core.dataarray import _THIS_ARRAY

        if name is _THIS_ARRAY:
            # When polyfit is called on a DataArray, ensure the resulting
            # dataset is backwards compatible with previous behavior
            name = ""
        elif isinstance(name, str):
            name = f"{name}_"
        else:
            # For other non-string names
            name = ""

        variables[name + "polyfit_coefficients"] = Variable(
            data=coeffs / scale_da, dims=(degree_dim,) + other_dims
        )

        if full or (cov is True):
            variables[name + "polyfit_residuals"] = Variable(
                data=residuals if var.ndim > 1 else residuals.squeeze(),
                dims=other_dims,
            )

        fac: Variable | int
        if cov:
            Vbase = np.linalg.inv(np.dot(lhs.T, lhs))
            Vbase /= np.outer(scale, scale)
            if cov == "unscaled":
                fac = 1
            else:
                if x.shape[0] <= order:
                    raise ValueError(
                        "The number of data points must exceed order to scale the covariance matrix."
                    )
                fac = variables[name + "polyfit_residuals"] / (x.shape[0] - order)
            variables[name + "polyfit_covariance"] = (
                Variable(data=Vbase, dims=("cov_i", "cov_j")) * fac
            )

    return type(obj)(
        data_vars=variables,
        coords={
            degree_dim: degree_coord_values,
            **{
                name: coord
                for name, coord in other_coords.items()
                if name in present_dims
            },
        },
        attrs=obj.attrs.copy(),
    )


def curvefit(
    obj,
    coords: str | DataArray | Iterable[str | DataArray],
    func: Callable[..., Any],
    reduce_dims: Dims = None,
    skipna: bool = True,
    p0: Mapping[str, float | DataArray] | None = None,
    bounds: Mapping[str, tuple[float | DataArray, float | DataArray]] | None = None,
    param_names: Sequence[str] | None = None,
    errors: ErrorOptions = "raise",
    kwargs: dict[str, Any] | None = None,
):
    """
    Curve fitting optimization for arbitrary functions.

    Wraps `scipy.optimize.curve_fit` with `apply_ufunc`.

    Parameters
    ----------
    obj : Dataset or DataArray
        Object to perform the curvefit on
    coords : hashable, DataArray, or sequence of hashable or DataArray
        Independent coordinate(s) over which to perform the curve fitting. Must share
        at least one dimension with the calling object. When fitting multi-dimensional
        functions, supply `coords` as a sequence in the same order as arguments in
        `func`. To fit along existing dimensions of the calling object, `coords` can
        also be specified as a str or sequence of strs.
    func : callable
        User specified function in the form `f(x, *params)` which returns a numpy
        array of length `len(x)`. `params` are the fittable parameters which are optimized
        by scipy curve_fit. `x` can also be specified as a sequence containing multiple
        coordinates, e.g. `f((x0, x1), *params)`.
    reduce_dims : str, Iterable of Hashable or None, optional
        Additional dimension(s) over which to aggregate while fitting. For example,
        calling `ds.curvefit(coords='time', reduce_dims=['lat', 'lon'], ...)` will
        aggregate all lat and lon points and fit the specified function along the
        time dimension.
    skipna : bool, default: True
        Whether to skip missing values when fitting. Default is True.
    p0 : dict-like, optional
        Optional dictionary of parameter names to initial guesses passed to the
        `curve_fit` `p0` arg. If the values are DataArrays, they will be appropriately
        broadcast to the coordinates of the array. If none or only some parameters are
        passed, the rest will be assigned initial values following the default scipy
        behavior.
    bounds : dict-like, optional
        Optional dictionary of parameter names to tuples of bounding values passed to the
        `curve_fit` `bounds` arg. If any of the bounds are DataArrays, they will be
        appropriately broadcast to the coordinates of the array. If none or only some
        parameters are passed, the rest will be unbounded following the default scipy
        behavior.
    param_names : sequence of hashable, optional
        Sequence of names for the fittable parameters of `func`. If not supplied,
        this will be automatically determined by arguments of `func`. `param_names`
        should be manually supplied when fitting a function that takes a variable
        number of parameters.
    errors : {"raise", "ignore"}, default: "raise"
        If 'raise', any errors from the `scipy.optimize_curve_fit` optimization will
        raise an exception. If 'ignore', the coefficients and covariances for the
        coordinates where the fitting failed will be NaN.
    kwargs : optional
        Additional keyword arguments to passed to scipy curve_fit.

    Returns
    -------
    Dataset
        A single dataset which contains:

        [var]_curvefit_coefficients
            The coefficients of the best fit.
        [var]_curvefit_covariance
            The covariance matrix of the coefficient estimates.

    See Also
    --------
    Dataset.polyfit
    scipy.optimize.curve_fit
    """
    from scipy.optimize import curve_fit

    if p0 is None:
        p0 = {}
    if bounds is None:
        bounds = {}
    if kwargs is None:
        kwargs = {}

    reduce_dims_: list[Hashable]
    if not reduce_dims:
        reduce_dims_ = []
    elif isinstance(reduce_dims, str) or not isinstance(reduce_dims, Iterable):
        reduce_dims_ = [reduce_dims]
    else:
        reduce_dims_ = list(reduce_dims)

    if isinstance(coords, str | DataArray) or not isinstance(coords, Iterable):
        coords = [coords]
    coords_: Sequence[DataArray] = [
        obj[coord] if isinstance(coord, str) else coord for coord in coords
    ]

    # Determine whether any coords are dims on self
    for coord in coords_:
        reduce_dims_ += [c for c in obj.dims if coord.equals(obj[c])]
    reduce_dims_ = list(set(reduce_dims_))
    preserved_dims = list(set(obj.dims) - set(reduce_dims_))
    if not reduce_dims_:
        raise ValueError(
            "No arguments to `coords` were identified as a dimension on the calling "
            "object, and no dims were supplied to `reduce_dims`. This would result "
            "in fitting on scalar data."
        )

    # Check that initial guess and bounds only contain coordinates that are in preserved_dims
    for param, guess in p0.items():
        if isinstance(guess, DataArray):
            unexpected = set(guess.dims) - set(preserved_dims)
            if unexpected:
                raise ValueError(
                    f"Initial guess for '{param}' has unexpected dimensions "
                    f"{tuple(unexpected)}. It should only have dimensions that are in data "
                    f"dimensions {preserved_dims}."
                )
    for param, (lb, ub) in bounds.items():
        for label, bound in zip(("Lower", "Upper"), (lb, ub), strict=True):
            if isinstance(bound, DataArray):
                unexpected = set(bound.dims) - set(preserved_dims)
                if unexpected:
                    raise ValueError(
                        f"{label} bound for '{param}' has unexpected dimensions "
                        f"{tuple(unexpected)}. It should only have dimensions that are in data "
                        f"dimensions {preserved_dims}."
                    )

    if errors not in ["raise", "ignore"]:
        raise ValueError('errors must be either "raise" or "ignore"')

    # Broadcast all coords with each other
    coords_ = broadcast(*coords_)
    coords_ = [coord.broadcast_like(obj, exclude=preserved_dims) for coord in coords_]
    n_coords = len(coords_)

    params, func_args = _get_func_args(func, param_names)
    param_defaults, bounds_defaults = _initialize_curvefit_params(
        params, p0, bounds, func_args
    )
    n_params = len(params)

    def _wrapper(Y, *args, **kwargs):
        # Wrap curve_fit with raveled coordinates and pointwise NaN handling
        # *args contains:
        #   - the coordinates
        #   - initial guess
        #   - lower bounds
        #   - upper bounds
        coords__ = args[:n_coords]
        p0_ = args[n_coords + 0 * n_params : n_coords + 1 * n_params]
        lb = args[n_coords + 1 * n_params : n_coords + 2 * n_params]
        ub = args[n_coords + 2 * n_params :]

        x = np.vstack([c.ravel() for c in coords__])
        y = Y.ravel()
        if skipna:
            mask = np.all([np.any(~np.isnan(x), axis=0), ~np.isnan(y)], axis=0)
            x = x[:, mask]
            y = y[mask]
            if y.size == 0:
                popt = np.full([n_params], np.nan)
                pcov = np.full([n_params, n_params], np.nan)
                return popt, pcov
        x = np.squeeze(x)

        try:
            popt, pcov = curve_fit(func, x, y, p0=p0_, bounds=(lb, ub), **kwargs)
        except RuntimeError:
            if errors == "raise":
                raise
            popt = np.full([n_params], np.nan)
            pcov = np.full([n_params, n_params], np.nan)

        return popt, pcov

    from xarray.core.dataarray import _THIS_ARRAY

    result = type(obj)()
    for name, da in obj.data_vars.items():
        if name is _THIS_ARRAY:
            # When curvefit is called on a DataArray, ensure the resulting
            # dataset is backwards compatible with previous behavior
            var_name = ""
        else:
            var_name = f"{name}_"

        input_core_dims = [reduce_dims_ for _ in range(n_coords + 1)]
        input_core_dims.extend(
            [[] for _ in range(3 * n_params)]
        )  # core_dims for p0 and bounds

        popt, pcov = apply_ufunc(
            _wrapper,
            da,
            *coords_,
            *param_defaults.values(),
            *[b[0] for b in bounds_defaults.values()],
            *[b[1] for b in bounds_defaults.values()],
            vectorize=True,
            dask="parallelized",
            input_core_dims=input_core_dims,
            output_core_dims=[["param"], ["cov_i", "cov_j"]],
            dask_gufunc_kwargs={
                "output_sizes": {
                    "param": n_params,
                    "cov_i": n_params,
                    "cov_j": n_params,
                },
            },
            output_dtypes=(np.float64, np.float64),
            exclude_dims=set(reduce_dims_),
            kwargs=kwargs,
        )
        result[var_name + "curvefit_coefficients"] = popt
        result[var_name + "curvefit_covariance"] = pcov

    result = result.assign_coords({"param": params, "cov_i": params, "cov_j": params})
    result.attrs = obj.attrs.copy()

    return result
