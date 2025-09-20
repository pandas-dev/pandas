from __future__ import annotations

import functools
import itertools
import math
import warnings
from collections.abc import Callable, Hashable, Iterator, Mapping
from typing import TYPE_CHECKING, Any, Generic, TypeVar

import numpy as np

from xarray.compat import dask_array_ops
from xarray.computation.arithmetic import CoarsenArithmetic
from xarray.core import dtypes, duck_array_ops, utils
from xarray.core.options import OPTIONS, _get_keep_attrs
from xarray.core.types import CoarsenBoundaryOptions, SideOptions, T_Xarray
from xarray.core.utils import (
    either_dict_or_kwargs,
    is_duck_dask_array,
    module_available,
)
from xarray.util.deprecation_helpers import _deprecate_positional_args

try:
    import bottleneck
except ImportError:
    # use numpy methods instead
    bottleneck = None

if TYPE_CHECKING:
    from xarray.core.dataarray import DataArray
    from xarray.core.dataset import Dataset

    RollingKey = Any
    _T = TypeVar("_T")

_ROLLING_REDUCE_DOCSTRING_TEMPLATE = """\
Reduce this object's data windows by applying `{name}` along its dimension.

Parameters
----------
keep_attrs : bool, default: None
    If True, the attributes (``attrs``) will be copied from the original
    object to the new one. If False, the new object will be returned
    without attributes. If None uses the global default.
**kwargs : dict
    Additional keyword arguments passed on to `{name}`.

Returns
-------
reduced : same type as caller
    New object with `{name}` applied along its rolling dimension.
"""


class Rolling(Generic[T_Xarray]):
    """A object that implements the moving window pattern.

    See Also
    --------
    xarray.Dataset.groupby
    xarray.DataArray.groupby
    xarray.Dataset.rolling
    xarray.DataArray.rolling
    """

    __slots__ = ("center", "dim", "min_periods", "obj", "window")
    _attributes = ("window", "min_periods", "center", "dim")
    dim: list[Hashable]
    window: list[int]
    center: list[bool]
    obj: T_Xarray
    min_periods: int

    def __init__(
        self,
        obj: T_Xarray,
        windows: Mapping[Any, int],
        min_periods: int | None = None,
        center: bool | Mapping[Any, bool] = False,
    ) -> None:
        """
        Moving window object.

        Parameters
        ----------
        obj : Dataset or DataArray
            Object to window.
        windows : mapping of hashable to int
            A mapping from the name of the dimension to create the rolling
            window along (e.g. `time`) to the size of the moving window.
        min_periods : int or None, default: None
            Minimum number of observations in window required to have a value
            (otherwise result is NA). The default, None, is equivalent to
            setting min_periods equal to the size of the window.
        center : bool or dict-like Hashable to bool, default: False
            Set the labels at the center of the window. If dict-like, set this
            property per rolling dimension.

        Returns
        -------
        rolling : type of input argument
        """
        self.dim = []
        self.window = []
        for d, w in windows.items():
            self.dim.append(d)
            if w <= 0:
                raise ValueError("window must be > 0")
            self.window.append(w)

        self.center = self._mapping_to_list(center, default=False)
        self.obj = obj

        missing_dims = tuple(dim for dim in self.dim if dim not in self.obj.dims)
        if missing_dims:
            # NOTE: we raise KeyError here but ValueError in Coarsen.
            raise KeyError(
                f"Window dimensions {missing_dims} not found in {self.obj.__class__.__name__} "
                f"dimensions {tuple(self.obj.dims)}"
            )

        # attributes
        if min_periods is not None and min_periods <= 0:
            raise ValueError("min_periods must be greater than zero or None")

        self.min_periods = (
            math.prod(self.window) if min_periods is None else min_periods
        )

    def __repr__(self) -> str:
        """provide a nice str repr of our rolling object"""

        attrs = ",".join(
            f"{k}->{w}{'(center)' if c else ''}"
            for k, w, c in zip(self.dim, self.window, self.center, strict=True)
        )
        return f"{self.__class__.__name__} [{attrs}]"

    def __len__(self) -> int:
        return math.prod(self.obj.sizes[d] for d in self.dim)

    @property
    def ndim(self) -> int:
        return len(self.dim)

    def _reduce_method(  # type: ignore[misc]
        name: str,
        fillna: Any,
        rolling_agg_func: Callable | None = None,
        automatic_rechunk: bool = False,
    ) -> Callable[..., T_Xarray]:
        """Constructs reduction methods built on a numpy reduction function (e.g. sum),
        a numbagg reduction function (e.g. move_sum), a bottleneck reduction function
        (e.g. move_sum), or a Rolling reduction (_mean).

        The logic here for which function to run is quite diffuse, across this method &
        _array_reduce. Arguably we could refactor this. But one constraint is that we
        need context of xarray options, of the functions each library offers, of
        the array (e.g. dtype).

        Set automatic_rechunk=True when the reduction method makes a memory copy.
        """
        if rolling_agg_func:
            array_agg_func = None
        else:
            array_agg_func = getattr(duck_array_ops, name)

        bottleneck_move_func = getattr(bottleneck, "move_" + name, None)
        if module_available("numbagg"):
            import numbagg

            numbagg_move_func = getattr(numbagg, "move_" + name, None)
        else:
            numbagg_move_func = None

        def method(self, keep_attrs=None, **kwargs):
            keep_attrs = self._get_keep_attrs(keep_attrs)

            return self._array_reduce(
                array_agg_func=array_agg_func,
                bottleneck_move_func=bottleneck_move_func,
                numbagg_move_func=numbagg_move_func,
                rolling_agg_func=rolling_agg_func,
                keep_attrs=keep_attrs,
                fillna=fillna,
                sliding_window_view_kwargs=dict(automatic_rechunk=automatic_rechunk),
                **kwargs,
            )

        method.__name__ = name
        method.__doc__ = _ROLLING_REDUCE_DOCSTRING_TEMPLATE.format(name=name)
        return method

    def _mean(self, keep_attrs, **kwargs):
        result = self.sum(keep_attrs=False, **kwargs)
        # use dtype of result for casting of count
        # this allows for GH #7062 and GH #8864, fixes GH #10340
        result /= duck_array_ops.astype(
            self.count(keep_attrs=False), dtype=result.dtype, copy=False
        )
        if keep_attrs:
            result.attrs = self.obj.attrs

        return result

    _mean.__doc__ = _ROLLING_REDUCE_DOCSTRING_TEMPLATE.format(name="mean")

    # automatic_rechunk is set to True for reductions that make a copy.
    # std, var could be optimized after which we can set it to False
    # See #4325
    argmax = _reduce_method("argmax", dtypes.NINF, automatic_rechunk=True)
    argmin = _reduce_method("argmin", dtypes.INF, automatic_rechunk=True)
    max = _reduce_method("max", dtypes.NINF)
    min = _reduce_method("min", dtypes.INF)
    prod = _reduce_method("prod", 1)
    sum = _reduce_method("sum", 0)
    mean = _reduce_method("mean", None, _mean)
    std = _reduce_method("std", None, automatic_rechunk=True)
    var = _reduce_method("var", None, automatic_rechunk=True)
    median = _reduce_method("median", None, automatic_rechunk=True)

    def _counts(self, keep_attrs: bool | None) -> T_Xarray:
        raise NotImplementedError()

    def count(self, keep_attrs: bool | None = None) -> T_Xarray:
        keep_attrs = self._get_keep_attrs(keep_attrs)
        rolling_count = self._counts(keep_attrs=keep_attrs)
        enough_periods = rolling_count >= self.min_periods
        return rolling_count.where(enough_periods)

    count.__doc__ = _ROLLING_REDUCE_DOCSTRING_TEMPLATE.format(name="count")

    def _mapping_to_list(
        self,
        arg: _T | Mapping[Any, _T],
        default: _T | None = None,
        allow_default: bool = True,
        allow_allsame: bool = True,
    ) -> list[_T]:
        if utils.is_dict_like(arg):
            if allow_default:
                return [arg.get(d, default) for d in self.dim]
            for d in self.dim:
                if d not in arg:
                    raise KeyError(f"Argument has no dimension key {d}.")
            return [arg[d] for d in self.dim]
        if allow_allsame:  # for single argument
            return [arg] * self.ndim  # type: ignore[list-item]  # no check for negatives
        if self.ndim == 1:
            return [arg]  # type: ignore[list-item]  # no check for negatives
        raise ValueError(f"Mapping argument is necessary for {self.ndim}d-rolling.")

    def _get_keep_attrs(self, keep_attrs):
        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)

        return keep_attrs


class DataArrayRolling(Rolling["DataArray"]):
    __slots__ = ("window_labels",)

    def __init__(
        self,
        obj: DataArray,
        windows: Mapping[Any, int],
        min_periods: int | None = None,
        center: bool | Mapping[Any, bool] = False,
    ) -> None:
        """
        Moving window object for DataArray.
        You should use DataArray.rolling() method to construct this object
        instead of the class constructor.

        Parameters
        ----------
        obj : DataArray
            Object to window.
        windows : mapping of hashable to int
            A mapping from the name of the dimension to create the rolling
            exponential window along (e.g. `time`) to the size of the moving window.
        min_periods : int, default: None
            Minimum number of observations in window required to have a value
            (otherwise result is NA). The default, None, is equivalent to
            setting min_periods equal to the size of the window.
        center : bool, default: False
            Set the labels at the center of the window. The default, False,
            sets the labels at the right edge of the window.

        Returns
        -------
        rolling : type of input argument

        See Also
        --------
        xarray.DataArray.rolling
        xarray.DataArray.groupby
        xarray.Dataset.rolling
        xarray.Dataset.groupby
        """
        super().__init__(obj, windows, min_periods=min_periods, center=center)

        # TODO legacy attribute
        self.window_labels = self.obj[self.dim[0]]

    def __iter__(self) -> Iterator[tuple[DataArray, DataArray]]:
        if self.ndim > 1:
            raise ValueError("__iter__ is only supported for 1d-rolling")

        dim0 = self.dim[0]
        window0 = int(self.window[0])
        offset = (window0 + 1) // 2 if self.center[0] else 1
        stops = np.arange(offset, self.obj.sizes[dim0] + offset)
        starts = stops - window0
        starts[: window0 - offset] = 0

        for label, start, stop in zip(self.window_labels, starts, stops, strict=True):
            window = self.obj.isel({dim0: slice(start, stop)})

            counts = window.count(dim=[dim0])
            window = window.where(counts >= self.min_periods)

            yield (label, window)

    @_deprecate_positional_args("v2024.11.0")
    def construct(
        self,
        window_dim: Hashable | Mapping[Any, Hashable] | None = None,
        *,
        stride: int | Mapping[Any, int] = 1,
        fill_value: Any = dtypes.NA,
        keep_attrs: bool | None = None,
        sliding_window_view_kwargs: Mapping[Any, Any] | None = None,
        **window_dim_kwargs: Hashable,
    ) -> DataArray:
        """
        Convert this rolling object to xr.DataArray,
        where the window dimension is stacked as a new dimension

        Parameters
        ----------
        window_dim : Hashable or dict-like to Hashable, optional
            A mapping from dimension name to the new window dimension names.
        stride : int or mapping of int, default: 1
            Size of stride for the rolling window.
        fill_value : default: dtypes.NA
            Filling value to match the dimension size.
        keep_attrs : bool, default: None
            If True, the attributes (``attrs``) will be copied from the original
            object to the new one. If False, the new object will be returned
            without attributes. If None uses the global default.
        sliding_window_view_kwargs : Mapping
            Keyword arguments that should be passed to the underlying array type's
            ``sliding_window_view`` function.
        **window_dim_kwargs : Hashable, optional
            The keyword arguments form of ``window_dim`` {dim: new_name, ...}.

        Returns
        -------
        DataArray
            a view of the original array. By default, the returned array is not writeable.
            For numpy arrays, one can pass ``writeable=True`` in ``sliding_window_view_kwargs``.

        See Also
        --------
        numpy.lib.stride_tricks.sliding_window_view
        dask.array.lib.stride_tricks.sliding_window_view

        Notes
        -----
        With dask arrays, it's possible to pass the ``automatic_rechunk`` kwarg as
        ``sliding_window_view_kwargs={"automatic_rechunk": True}``. This controls
        whether dask should automatically rechunk the output to avoid
        exploding chunk sizes. Automatically rechunking is the default behaviour.
        Importantly, each chunk will be a view of the data so large chunk sizes are
        only safe if *no* copies are made later.

        Examples
        --------
        >>> da = xr.DataArray(np.arange(8).reshape(2, 4), dims=("a", "b"))

        >>> rolling = da.rolling(b=3)
        >>> rolling.construct("window_dim")
        <xarray.DataArray (a: 2, b: 4, window_dim: 3)> Size: 192B
        array([[[nan, nan,  0.],
                [nan,  0.,  1.],
                [ 0.,  1.,  2.],
                [ 1.,  2.,  3.]],
        <BLANKLINE>
               [[nan, nan,  4.],
                [nan,  4.,  5.],
                [ 4.,  5.,  6.],
                [ 5.,  6.,  7.]]])
        Dimensions without coordinates: a, b, window_dim

        >>> rolling = da.rolling(b=3, center=True)
        >>> rolling.construct("window_dim")
        <xarray.DataArray (a: 2, b: 4, window_dim: 3)> Size: 192B
        array([[[nan,  0.,  1.],
                [ 0.,  1.,  2.],
                [ 1.,  2.,  3.],
                [ 2.,  3., nan]],
        <BLANKLINE>
               [[nan,  4.,  5.],
                [ 4.,  5.,  6.],
                [ 5.,  6.,  7.],
                [ 6.,  7., nan]]])
        Dimensions without coordinates: a, b, window_dim

        """

        if sliding_window_view_kwargs is None:
            sliding_window_view_kwargs = {}
        return self._construct(
            self.obj,
            window_dim=window_dim,
            stride=stride,
            fill_value=fill_value,
            keep_attrs=keep_attrs,
            sliding_window_view_kwargs=sliding_window_view_kwargs,
            **window_dim_kwargs,
        )

    def _construct(
        self,
        obj: DataArray,
        *,
        window_dim: Hashable | Mapping[Any, Hashable] | None = None,
        stride: int | Mapping[Any, int] = 1,
        fill_value: Any = dtypes.NA,
        keep_attrs: bool | None = None,
        sliding_window_view_kwargs: Mapping[Any, Any] | None = None,
        **window_dim_kwargs: Hashable,
    ) -> DataArray:
        from xarray.core.dataarray import DataArray

        if sliding_window_view_kwargs is None:
            sliding_window_view_kwargs = {}

        keep_attrs = self._get_keep_attrs(keep_attrs)

        if window_dim is None:
            if len(window_dim_kwargs) == 0:
                raise ValueError(
                    "Either window_dim or window_dim_kwargs need to be specified."
                )
            window_dim = {d: window_dim_kwargs[str(d)] for d in self.dim}

        window_dims = self._mapping_to_list(
            window_dim, allow_default=False, allow_allsame=False
        )
        strides = self._mapping_to_list(stride, default=1)

        window = obj.variable.rolling_window(
            self.dim,
            self.window,
            window_dims,
            center=self.center,
            fill_value=fill_value,
            **sliding_window_view_kwargs,
        )

        attrs = obj.attrs if keep_attrs else {}

        result = DataArray(
            window,
            dims=obj.dims + tuple(window_dims),
            coords=obj.coords,
            attrs=attrs,
            name=obj.name,
        )
        return result.isel(
            {d: slice(None, None, s) for d, s in zip(self.dim, strides, strict=True)}
        )

    def reduce(
        self,
        func: Callable,
        keep_attrs: bool | None = None,
        *,
        sliding_window_view_kwargs: Mapping[Any, Any] | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """Reduce each window by applying `func`.

        Equivalent to ``.construct(...).reduce(func, ...)``.

        Parameters
        ----------
        func : callable
            Function which can be called in the form
            `func(x, **kwargs)` to return the result of collapsing an
            np.ndarray over an the rolling dimension.
        keep_attrs : bool, default: None
            If True, the attributes (``attrs``) will be copied from the original
            object to the new one. If False, the new object will be returned
            without attributes. If None uses the global default.
        sliding_window_view_kwargs
            Keyword arguments that should be passed to the underlying array type's
            ``sliding_window_view`` function.
        **kwargs : dict
            Additional keyword arguments passed on to `func`.

        Returns
        -------
        reduced : DataArray
            Array with summarized data.

        See Also
        --------
        numpy.lib.stride_tricks.sliding_window_view
        dask.array.lib.stride_tricks.sliding_window_view

        Notes
        -----
        With dask arrays, it's possible to pass the ``automatic_rechunk`` kwarg as
        ``sliding_window_view_kwargs={"automatic_rechunk": True}``. This controls
        whether dask should automatically rechunk the output to avoid
        exploding chunk sizes. Automatically rechunking is the default behaviour.
        Importantly, each chunk will be a view of the data so large chunk sizes are
        only safe if *no* copies are made later.

        Examples
        --------
        >>> da = xr.DataArray(np.arange(8).reshape(2, 4), dims=("a", "b"))
        >>> rolling = da.rolling(b=3)
        >>> rolling.construct("window_dim")
        <xarray.DataArray (a: 2, b: 4, window_dim: 3)> Size: 192B
        array([[[nan, nan,  0.],
                [nan,  0.,  1.],
                [ 0.,  1.,  2.],
                [ 1.,  2.,  3.]],
        <BLANKLINE>
               [[nan, nan,  4.],
                [nan,  4.,  5.],
                [ 4.,  5.,  6.],
                [ 5.,  6.,  7.]]])
        Dimensions without coordinates: a, b, window_dim

        >>> rolling.reduce(np.sum)
        <xarray.DataArray (a: 2, b: 4)> Size: 64B
        array([[nan, nan,  3.,  6.],
               [nan, nan, 15., 18.]])
        Dimensions without coordinates: a, b

        >>> rolling = da.rolling(b=3, min_periods=1)
        >>> rolling.reduce(np.nansum)
        <xarray.DataArray (a: 2, b: 4)> Size: 64B
        array([[ 0.,  1.,  3.,  6.],
               [ 4.,  9., 15., 18.]])
        Dimensions without coordinates: a, b
        """

        keep_attrs = self._get_keep_attrs(keep_attrs)

        rolling_dim = {
            d: utils.get_temp_dimname(self.obj.dims, f"_rolling_dim_{d}")
            for d in self.dim
        }

        # save memory with reductions GH4325
        fillna = kwargs.pop("fillna", dtypes.NA)
        if fillna is not dtypes.NA:
            obj = self.obj.fillna(fillna)
        else:
            obj = self.obj
        windows = self._construct(
            obj,
            window_dim=rolling_dim,
            keep_attrs=keep_attrs,
            fill_value=fillna,
            sliding_window_view_kwargs=sliding_window_view_kwargs,
        )

        dim = list(rolling_dim.values())
        result = windows.reduce(func, dim=dim, keep_attrs=keep_attrs, **kwargs)

        # Find valid windows based on count.
        counts = self._counts(keep_attrs=False)
        return result.where(counts >= self.min_periods)

    def _counts(self, keep_attrs: bool | None) -> DataArray:
        """Number of non-nan entries in each rolling window."""

        rolling_dim = {
            d: utils.get_temp_dimname(self.obj.dims, f"_rolling_dim_{d}")
            for d in self.dim
        }
        # We use False as the fill_value instead of np.nan, since boolean
        # array is faster to be reduced than object array.
        # The use of skipna==False is also faster since it does not need to
        # copy the strided array.
        dim = list(rolling_dim.values())
        counts = (
            self.obj.notnull(keep_attrs=keep_attrs)
            .rolling(
                dict(zip(self.dim, self.window, strict=True)),
                center={d: self.center[i] for i, d in enumerate(self.dim)},
            )
            .construct(rolling_dim, fill_value=False, keep_attrs=keep_attrs)
            .sum(dim=dim, skipna=False, keep_attrs=keep_attrs)
        )
        return counts

    def _numbagg_reduce(self, func, keep_attrs, **kwargs):
        # Some of this is copied from `_bottleneck_reduce`, we could reduce this as part
        # of a wider refactor.

        axis = self.obj.get_axis_num(self.dim[0])

        padded = self.obj.variable
        if self.center[0]:
            if is_duck_dask_array(padded.data):
                # workaround to make the padded chunk size larger than
                # self.window - 1
                shift = -(self.window[0] + 1) // 2
                offset = (self.window[0] - 1) // 2
                valid = (slice(None),) * axis + (
                    slice(offset, offset + self.obj.shape[axis]),
                )
            else:
                shift = (-self.window[0] // 2) + 1
                valid = (slice(None),) * axis + (slice(-shift, None),)
            padded = padded.pad({self.dim[0]: (0, -shift)}, mode="constant")

        if is_duck_dask_array(padded.data) and False:
            raise AssertionError("should not be reachable")
        else:
            values = func(
                padded.data,
                window=self.window[0],
                min_count=self.min_periods,
                axis=axis,
            )

        if self.center[0]:
            values = values[valid]

        attrs = self.obj.attrs if keep_attrs else {}

        return self.obj.__class__(
            values, self.obj.coords, attrs=attrs, name=self.obj.name
        )

    def _bottleneck_reduce(self, func, keep_attrs, **kwargs):
        # bottleneck doesn't allow min_count to be 0, although it should
        # work the same as if min_count = 1
        # Note bottleneck only works with 1d-rolling.
        if self.min_periods == 0:
            min_count = 1
        else:
            min_count = self.min_periods

        axis = self.obj.get_axis_num(self.dim[0])

        padded = self.obj.variable
        if self.center[0]:
            if is_duck_dask_array(padded.data):
                # workaround to make the padded chunk size larger than
                # self.window - 1
                shift = -(self.window[0] + 1) // 2
                offset = (self.window[0] - 1) // 2
                valid = (slice(None),) * axis + (
                    slice(offset, offset + self.obj.shape[axis]),
                )
            else:
                shift = (-self.window[0] // 2) + 1
                valid = (slice(None),) * axis + (slice(-shift, None),)
            padded = padded.pad({self.dim[0]: (0, -shift)}, mode="constant")

        if is_duck_dask_array(padded.data):
            values = dask_array_ops.dask_rolling_wrapper(
                func, padded, axis=axis, window=self.window[0], min_count=min_count
            )
        else:
            values = func(
                padded.data, window=self.window[0], min_count=min_count, axis=axis
            )
        # index 0 is at the rightmost edge of the window
        # need to reverse index here
        # see GH #8541
        if func in [bottleneck.move_argmin, bottleneck.move_argmax]:
            values = self.window[0] - 1 - values

        if self.center[0]:
            values = values[valid]

        attrs = self.obj.attrs if keep_attrs else {}

        return self.obj.__class__(
            values, self.obj.coords, attrs=attrs, name=self.obj.name
        )

    def _array_reduce(
        self,
        array_agg_func,
        bottleneck_move_func,
        numbagg_move_func,
        rolling_agg_func,
        keep_attrs,
        fillna,
        **kwargs,
    ):
        if "dim" in kwargs:
            warnings.warn(
                f"Reductions are applied along the rolling dimension(s) "
                f"'{self.dim}'. Passing the 'dim' kwarg to reduction "
                f"operations has no effect.",
                DeprecationWarning,
                stacklevel=3,
            )
            del kwargs["dim"]

        xp = duck_array_ops.get_array_namespace(self.obj.data)
        if (
            OPTIONS["use_numbagg"]
            and module_available("numbagg")
            and numbagg_move_func is not None
            # TODO: we could at least allow this for the equivalent of `apply_ufunc`'s
            # "parallelized". `rolling_exp` does this, as an example (but rolling_exp is
            # much simpler)
            and not is_duck_dask_array(self.obj.data)
            # Numbagg doesn't handle object arrays and generally has dtype consistency,
            # so doesn't deal well with bool arrays which are expected to change type.
            and self.obj.data.dtype.kind not in "ObMm"
            # TODO: we could also allow this, probably as part of a refactoring of this
            # module, so we can use the machinery in `self.reduce`.
            and self.ndim == 1
            and xp is np
        ):
            import numbagg

            # Numbagg has a default ddof of 1. I (@max-sixty) think we should make
            # this the default in xarray too, but until we do, don't use numbagg for
            # std and var unless ddof is set to 1.
            if (
                numbagg_move_func not in [numbagg.move_std, numbagg.move_var]
                or kwargs.get("ddof") == 1
            ):
                return self._numbagg_reduce(
                    numbagg_move_func, keep_attrs=keep_attrs, **kwargs
                )

        if (
            OPTIONS["use_bottleneck"]
            and bottleneck_move_func is not None
            and (
                not is_duck_dask_array(self.obj.data)
                or module_available("dask", "2024.11.0")
            )
            and self.ndim == 1
            and xp is np
        ):
            return self._bottleneck_reduce(
                bottleneck_move_func, keep_attrs=keep_attrs, **kwargs
            )

        if rolling_agg_func:
            return rolling_agg_func(self, keep_attrs=self._get_keep_attrs(keep_attrs))

        if fillna is not None:
            if fillna is dtypes.INF:
                fillna = dtypes.get_pos_infinity(self.obj.dtype, max_for_int=True)
            elif fillna is dtypes.NINF:
                fillna = dtypes.get_neg_infinity(self.obj.dtype, min_for_int=True)
            kwargs.setdefault("skipna", False)
            kwargs.setdefault("fillna", fillna)

        return self.reduce(array_agg_func, keep_attrs=keep_attrs, **kwargs)


class DatasetRolling(Rolling["Dataset"]):
    __slots__ = ("rollings",)

    def __init__(
        self,
        obj: Dataset,
        windows: Mapping[Any, int],
        min_periods: int | None = None,
        center: bool | Mapping[Any, bool] = False,
    ) -> None:
        """
        Moving window object for Dataset.
        You should use Dataset.rolling() method to construct this object
        instead of the class constructor.

        Parameters
        ----------
        obj : Dataset
            Object to window.
        windows : mapping of hashable to int
            A mapping from the name of the dimension to create the rolling
            exponential window along (e.g. `time`) to the size of the moving window.
        min_periods : int, default: None
            Minimum number of observations in window required to have a value
            (otherwise result is NA). The default, None, is equivalent to
            setting min_periods equal to the size of the window.
        center : bool or mapping of hashable to bool, default: False
            Set the labels at the center of the window. The default, False,
            sets the labels at the right edge of the window.

        Returns
        -------
        rolling : type of input argument

        See Also
        --------
        xarray.Dataset.rolling
        xarray.DataArray.rolling
        xarray.Dataset.groupby
        xarray.DataArray.groupby
        """
        super().__init__(obj, windows, min_periods, center)

        # Keep each Rolling object as a dictionary
        self.rollings = {}
        for key, da in self.obj.data_vars.items():
            # keeps rollings only for the dataset depending on self.dim
            dims, center = [], {}
            for i, d in enumerate(self.dim):
                if d in da.dims:
                    dims.append(d)
                    center[d] = self.center[i]

            if dims:
                w = {d: windows[d] for d in dims}
                self.rollings[key] = DataArrayRolling(da, w, min_periods, center)

    def _dataset_implementation(self, func, keep_attrs, **kwargs):
        from xarray.core.dataset import Dataset

        keep_attrs = self._get_keep_attrs(keep_attrs)

        reduced = {}
        for key, da in self.obj.data_vars.items():
            if any(d in da.dims for d in self.dim):
                reduced[key] = func(self.rollings[key], keep_attrs=keep_attrs, **kwargs)
            else:
                reduced[key] = self.obj[key].copy()
                # we need to delete the attrs of the copied DataArray
                if not keep_attrs:
                    reduced[key].attrs = {}

        attrs = self.obj.attrs if keep_attrs else {}
        return Dataset(reduced, coords=self.obj.coords, attrs=attrs)

    def reduce(
        self,
        func: Callable,
        keep_attrs: bool | None = None,
        sliding_window_view_kwargs: Mapping[Any, Any] | None = None,
        **kwargs: Any,
    ) -> DataArray:
        """Reduce the items in this group by applying `func` along some
        dimension(s).

        Parameters
        ----------
        func : callable
            Function which can be called in the form
            `func(x, **kwargs)` to return the result of collapsing an
            np.ndarray over an the rolling dimension.
        keep_attrs : bool, default: None
            If True, the attributes (``attrs``) will be copied from the original
            object to the new one. If False, the new object will be returned
            without attributes. If None uses the global default.
        sliding_window_view_kwargs : Mapping
            Keyword arguments that should be passed to the underlying array type's
            ``sliding_window_view`` function.
        **kwargs : dict
            Additional keyword arguments passed on to `func`.

        Returns
        -------
        reduced : DataArray
            Array with summarized data.

        See Also
        --------
        numpy.lib.stride_tricks.sliding_window_view
        dask.array.lib.stride_tricks.sliding_window_view

        Notes
        -----
        With dask arrays, it's possible to pass the ``automatic_rechunk`` kwarg as
        ``sliding_window_view_kwargs={"automatic_rechunk": True}``. This controls
        whether dask should automatically rechunk the output to avoid
        exploding chunk sizes. Automatically rechunking is the default behaviour.
        Importantly, each chunk will be a view of the data so large chunk sizes are
        only safe if *no* copies are made later.
        """
        return self._dataset_implementation(
            functools.partial(DataArrayRolling.reduce, func=func),
            keep_attrs=keep_attrs,
            sliding_window_view_kwargs=sliding_window_view_kwargs,
            **kwargs,
        )

    def _counts(self, keep_attrs: bool | None) -> Dataset:
        return self._dataset_implementation(
            DataArrayRolling._counts, keep_attrs=keep_attrs
        )

    def _array_reduce(
        self,
        array_agg_func,
        bottleneck_move_func,
        rolling_agg_func,
        keep_attrs,
        **kwargs,
    ):
        return self._dataset_implementation(
            functools.partial(
                DataArrayRolling._array_reduce,
                array_agg_func=array_agg_func,
                bottleneck_move_func=bottleneck_move_func,
                rolling_agg_func=rolling_agg_func,
            ),
            keep_attrs=keep_attrs,
            **kwargs,
        )

    @_deprecate_positional_args("v2024.11.0")
    def construct(
        self,
        window_dim: Hashable | Mapping[Any, Hashable] | None = None,
        *,
        stride: int | Mapping[Any, int] = 1,
        fill_value: Any = dtypes.NA,
        keep_attrs: bool | None = None,
        sliding_window_view_kwargs: Mapping[Any, Any] | None = None,
        **window_dim_kwargs: Hashable,
    ) -> Dataset:
        """
        Convert this rolling object to xr.Dataset,
        where the window dimension is stacked as a new dimension

        Parameters
        ----------
        window_dim : str or mapping, optional
            A mapping from dimension name to the new window dimension names.
            Just a string can be used for 1d-rolling.
        stride : int, optional
            size of stride for the rolling window.
        fill_value : Any, default: dtypes.NA
            Filling value to match the dimension size.
        sliding_window_view_kwargs
            Keyword arguments that should be passed to the underlying array type's
            ``sliding_window_view`` function.
        **window_dim_kwargs : {dim: new_name, ...}, optional
            The keyword arguments form of ``window_dim``.

        Returns
        -------
        Dataset
            Dataset with views of the original arrays. By default, the returned arrays are not writeable.
            For numpy arrays, one can pass ``writeable=True`` in ``sliding_window_view_kwargs``.

        See Also
        --------
        numpy.lib.stride_tricks.sliding_window_view
        dask.array.lib.stride_tricks.sliding_window_view

        Notes
        -----
        With dask arrays, it's possible to pass the ``automatic_rechunk`` kwarg as
        ``sliding_window_view_kwargs={"automatic_rechunk": True}``. This controls
        whether dask should automatically rechunk the output to avoid
        exploding chunk sizes. Automatically rechunking is the default behaviour.
        Importantly, each chunk will be a view of the data so large chunk sizes are
        only safe if *no* copies are made later.
        """

        from xarray.core.dataset import Dataset

        keep_attrs = self._get_keep_attrs(keep_attrs)

        if window_dim is None:
            if len(window_dim_kwargs) == 0:
                raise ValueError(
                    "Either window_dim or window_dim_kwargs need to be specified."
                )
            window_dim = {d: window_dim_kwargs[str(d)] for d in self.dim}

        window_dims = self._mapping_to_list(
            window_dim, allow_default=False, allow_allsame=False
        )
        strides = self._mapping_to_list(stride, default=1)

        dataset = {}
        for key, da in self.obj.data_vars.items():
            # keeps rollings only for the dataset depending on self.dim
            dims = [d for d in self.dim if d in da.dims]
            if dims:
                wi = {d: window_dims[i] for i, d in enumerate(self.dim) if d in da.dims}
                st = {d: strides[i] for i, d in enumerate(self.dim) if d in da.dims}

                dataset[key] = self.rollings[key].construct(
                    window_dim=wi,
                    fill_value=fill_value,
                    stride=st,
                    keep_attrs=keep_attrs,
                    sliding_window_view_kwargs=sliding_window_view_kwargs,
                )
            else:
                dataset[key] = da.copy()

            # as the DataArrays can be copied we need to delete the attrs
            if not keep_attrs:
                dataset[key].attrs = {}

        # Need to stride coords as well. TODO: is there a better way?
        coords = self.obj.isel(
            {d: slice(None, None, s) for d, s in zip(self.dim, strides, strict=True)}
        ).coords

        attrs = self.obj.attrs if keep_attrs else {}

        return Dataset(dataset, coords=coords, attrs=attrs)


class Coarsen(CoarsenArithmetic, Generic[T_Xarray]):
    """A object that implements the coarsen.

    See Also
    --------
    Dataset.coarsen
    DataArray.coarsen
    """

    __slots__ = (
        "boundary",
        "coord_func",
        "obj",
        "side",
        "trim_excess",
        "windows",
    )
    _attributes = ("windows", "side", "trim_excess")
    obj: T_Xarray
    windows: Mapping[Hashable, int]
    side: SideOptions | Mapping[Hashable, SideOptions]
    boundary: CoarsenBoundaryOptions
    coord_func: Mapping[Hashable, str | Callable]

    def __init__(
        self,
        obj: T_Xarray,
        windows: Mapping[Any, int],
        boundary: CoarsenBoundaryOptions,
        side: SideOptions | Mapping[Any, SideOptions],
        coord_func: str | Callable | Mapping[Any, str | Callable],
    ) -> None:
        """
        Moving window object.

        Parameters
        ----------
        obj : Dataset or DataArray
            Object to window.
        windows : mapping of hashable to int
            A mapping from the name of the dimension to create the rolling
            exponential window along (e.g. `time`) to the size of the moving window.
        boundary : {"exact", "trim", "pad"}
            If 'exact', a ValueError will be raised if dimension size is not a
            multiple of window size. If 'trim', the excess indexes are trimmed.
            If 'pad', NA will be padded.
        side : 'left' or 'right' or mapping from dimension to 'left' or 'right'
        coord_func : function (name) or mapping from coordinate name to function (name).

        Returns
        -------
        coarsen

        """
        self.obj = obj
        self.windows = windows
        self.side = side
        self.boundary = boundary

        missing_dims = tuple(dim for dim in windows.keys() if dim not in self.obj.dims)
        if missing_dims:
            raise ValueError(
                f"Window dimensions {missing_dims} not found in {self.obj.__class__.__name__} "
                f"dimensions {tuple(self.obj.dims)}"
            )

        if utils.is_dict_like(coord_func):
            coord_func_map = coord_func
        else:
            coord_func_map = dict.fromkeys(self.obj.dims, coord_func)
        for c in self.obj.coords:
            if c not in coord_func_map:
                coord_func_map[c] = duck_array_ops.mean  # type: ignore[index]
        self.coord_func = coord_func_map

    def _get_keep_attrs(self, keep_attrs):
        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)

        return keep_attrs

    def __repr__(self) -> str:
        """provide a nice str repr of our coarsen object"""

        attrs = ",".join(
            f"{k}->{getattr(self, k)}"
            for k in self._attributes
            if getattr(self, k, None) is not None
        )
        return f"{self.__class__.__name__} [{attrs}]"

    def construct(
        self,
        window_dim=None,
        keep_attrs=None,
        **window_dim_kwargs,
    ) -> T_Xarray:
        """
        Convert this Coarsen object to a DataArray or Dataset,
        where the coarsening dimension is split or reshaped to two
        new dimensions.

        Parameters
        ----------
        window_dim: mapping
            A mapping from existing dimension name to new dimension names.
            The size of the second dimension will be the length of the
            coarsening window.
        keep_attrs: bool, optional
            Preserve attributes if True
        **window_dim_kwargs : {dim: new_name, ...}
            The keyword arguments form of ``window_dim``.

        Returns
        -------
        Dataset or DataArray with reshaped dimensions

        Examples
        --------
        >>> da = xr.DataArray(np.arange(24), dims="time")
        >>> da.coarsen(time=12).construct(time=("year", "month"))
        <xarray.DataArray (year: 2, month: 12)> Size: 192B
        array([[ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11],
               [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]])
        Dimensions without coordinates: year, month

        See Also
        --------
        DataArrayRolling.construct
        DatasetRolling.construct
        """

        from xarray.core.dataarray import DataArray
        from xarray.core.dataset import Dataset

        window_dim = either_dict_or_kwargs(
            window_dim, window_dim_kwargs, "Coarsen.construct"
        )
        if not window_dim:
            raise ValueError(
                "Either window_dim or window_dim_kwargs need to be specified."
            )

        bad_new_dims = tuple(
            win
            for win, dims in window_dim.items()
            if len(dims) != 2 or isinstance(dims, str)
        )
        if bad_new_dims:
            raise ValueError(
                f"Please provide exactly two dimension names for the following coarsening dimensions: {bad_new_dims}"
            )

        if keep_attrs is None:
            keep_attrs = _get_keep_attrs(default=True)

        missing_dims = set(window_dim) - set(self.windows)
        if missing_dims:
            raise ValueError(
                f"'window_dim' must contain entries for all dimensions to coarsen. Missing {missing_dims}"
            )
        extra_windows = set(self.windows) - set(window_dim)
        if extra_windows:
            raise ValueError(
                f"'window_dim' includes dimensions that will not be coarsened: {extra_windows}"
            )

        reshaped = Dataset()
        if isinstance(self.obj, DataArray):
            obj = self.obj._to_temp_dataset()
        else:
            obj = self.obj

        reshaped.attrs = obj.attrs if keep_attrs else {}

        for key, var in obj.variables.items():
            reshaped_dims = tuple(
                itertools.chain(*[window_dim.get(dim, [dim]) for dim in list(var.dims)])
            )
            if reshaped_dims != var.dims:
                windows = {w: self.windows[w] for w in window_dim if w in var.dims}
                reshaped_var, _ = var.coarsen_reshape(windows, self.boundary, self.side)
                attrs = var.attrs if keep_attrs else {}
                reshaped[key] = (reshaped_dims, reshaped_var, attrs)
            else:
                reshaped[key] = var

        # should handle window_dim being unindexed
        should_be_coords = (set(window_dim) & set(self.obj.coords)) | set(
            self.obj.coords
        )
        result = reshaped.set_coords(should_be_coords)
        if isinstance(self.obj, DataArray):
            return self.obj._from_temp_dataset(result)
        else:
            return result


class DataArrayCoarsen(Coarsen["DataArray"]):
    __slots__ = ()

    _reduce_extra_args_docstring = """"""

    @classmethod
    def _reduce_method(
        cls, func: Callable, include_skipna: bool = False, numeric_only: bool = False
    ) -> Callable[..., DataArray]:
        """
        Return a wrapped function for injecting reduction methods.
        see ops.inject_reduce_methods
        """
        kwargs: dict[str, Any] = {}
        if include_skipna:
            kwargs["skipna"] = None

        def wrapped_func(
            self: DataArrayCoarsen, keep_attrs: bool | None = None, **kwargs
        ) -> DataArray:
            from xarray.core.dataarray import DataArray

            keep_attrs = self._get_keep_attrs(keep_attrs)

            reduced = self.obj.variable.coarsen(
                self.windows, func, self.boundary, self.side, keep_attrs, **kwargs
            )
            coords = {}
            for c, v in self.obj.coords.items():
                if c == self.obj.name:
                    coords[c] = reduced
                elif any(d in self.windows for d in v.dims):
                    coords[c] = v.variable.coarsen(
                        self.windows,
                        self.coord_func[c],
                        self.boundary,
                        self.side,
                        keep_attrs,
                        **kwargs,
                    )
                else:
                    coords[c] = v
            return DataArray(
                reduced, dims=self.obj.dims, coords=coords, name=self.obj.name
            )

        return wrapped_func

    def reduce(
        self, func: Callable, keep_attrs: bool | None = None, **kwargs
    ) -> DataArray:
        """Reduce the items in this group by applying `func` along some
        dimension(s).

        Parameters
        ----------
        func : callable
            Function which can be called in the form `func(x, axis, **kwargs)`
            to return the result of collapsing an np.ndarray over the coarsening
            dimensions.  It must be possible to provide the `axis` argument
            with a tuple of integers.
        keep_attrs : bool, default: None
            If True, the attributes (``attrs``) will be copied from the original
            object to the new one. If False, the new object will be returned
            without attributes. If None uses the global default.
        **kwargs : dict
            Additional keyword arguments passed on to `func`.

        Returns
        -------
        reduced : DataArray
            Array with summarized data.

        Examples
        --------
        >>> da = xr.DataArray(np.arange(8).reshape(2, 4), dims=("a", "b"))
        >>> coarsen = da.coarsen(b=2)
        >>> coarsen.reduce(np.sum)
        <xarray.DataArray (a: 2, b: 2)> Size: 32B
        array([[ 1,  5],
               [ 9, 13]])
        Dimensions without coordinates: a, b
        """
        wrapped_func = self._reduce_method(func)
        return wrapped_func(self, keep_attrs=keep_attrs, **kwargs)


class DatasetCoarsen(Coarsen["Dataset"]):
    __slots__ = ()

    _reduce_extra_args_docstring = """"""

    @classmethod
    def _reduce_method(
        cls, func: Callable, include_skipna: bool = False, numeric_only: bool = False
    ) -> Callable[..., Dataset]:
        """
        Return a wrapped function for injecting reduction methods.
        see ops.inject_reduce_methods
        """
        kwargs: dict[str, Any] = {}
        if include_skipna:
            kwargs["skipna"] = None

        def wrapped_func(
            self: DatasetCoarsen, keep_attrs: bool | None = None, **kwargs
        ) -> Dataset:
            from xarray.core.dataset import Dataset

            keep_attrs = self._get_keep_attrs(keep_attrs)

            if keep_attrs:
                attrs = self.obj.attrs
            else:
                attrs = {}

            reduced = {}
            for key, da in self.obj.data_vars.items():
                reduced[key] = da.variable.coarsen(
                    self.windows,
                    func,
                    self.boundary,
                    self.side,
                    keep_attrs=keep_attrs,
                    **kwargs,
                )

            coords = {}
            for c, v in self.obj.coords.items():
                # variable.coarsen returns variables not containing the window dims
                # unchanged (maybe removes attrs)
                coords[c] = v.variable.coarsen(
                    self.windows,
                    self.coord_func[c],
                    self.boundary,
                    self.side,
                    keep_attrs=keep_attrs,
                    **kwargs,
                )

            return Dataset(reduced, coords=coords, attrs=attrs)

        return wrapped_func

    def reduce(self, func: Callable, keep_attrs=None, **kwargs) -> Dataset:
        """Reduce the items in this group by applying `func` along some
        dimension(s).

        Parameters
        ----------
        func : callable
            Function which can be called in the form `func(x, axis, **kwargs)`
            to return the result of collapsing an np.ndarray over the coarsening
            dimensions.  It must be possible to provide the `axis` argument with
            a tuple of integers.
        keep_attrs : bool, default: None
            If True, the attributes (``attrs``) will be copied from the original
            object to the new one. If False, the new object will be returned
            without attributes. If None uses the global default.
        **kwargs : dict
            Additional keyword arguments passed on to `func`.

        Returns
        -------
        reduced : Dataset
            Arrays with summarized data.
        """
        wrapped_func = self._reduce_method(func)
        return wrapped_func(self, keep_attrs=keep_attrs, **kwargs)
