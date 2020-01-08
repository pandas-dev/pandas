"""
Provide a generic structure to support window functions,
similar to how we have a Groupby object.
"""
from datetime import timedelta
from functools import partial
import inspect
from textwrap import dedent
from typing import Callable, Dict, List, Optional, Set, Tuple, Union

import numpy as np

import pandas._libs.window.aggregations as window_aggregations
from pandas._typing import Axis, FrameOrSeries, Scalar
from pandas.compat._optional import import_optional_dependency
from pandas.compat.numpy import function as nv
from pandas.util._decorators import Appender, Substitution, cache_readonly

from pandas.core.dtypes.common import (
    ensure_float64,
    is_bool,
    is_float_dtype,
    is_integer,
    is_integer_dtype,
    is_list_like,
    is_scalar,
    needs_i8_conversion,
)
from pandas.core.dtypes.generic import (
    ABCDataFrame,
    ABCDateOffset,
    ABCDatetimeIndex,
    ABCPeriodIndex,
    ABCSeries,
    ABCTimedeltaIndex,
)

from pandas.core.base import DataError, PandasObject, SelectionMixin, ShallowMixin
import pandas.core.common as com
from pandas.core.indexes.api import Index, ensure_index
from pandas.core.window.common import (
    WindowGroupByMixin,
    _doc_template,
    _flex_binary_moment,
    _shared_docs,
    calculate_center_offset,
    calculate_min_periods,
    get_weighted_roll_func,
    zsqrt,
)
from pandas.core.window.indexers import (
    BaseIndexer,
    FixedWindowIndexer,
    VariableWindowIndexer,
)
from pandas.core.window.numba_ import generate_numba_apply_func


class _Window(PandasObject, ShallowMixin, SelectionMixin):
    _attributes: List[str] = [
        "window",
        "min_periods",
        "center",
        "win_type",
        "axis",
        "on",
        "closed",
    ]
    exclusions: Set[str] = set()

    def __init__(
        self,
        obj,
        window=None,
        min_periods: Optional[int] = None,
        center: Optional[bool] = False,
        win_type: Optional[str] = None,
        axis: Axis = 0,
        on: Optional[Union[str, Index]] = None,
        closed: Optional[str] = None,
        **kwargs,
    ):

        self.__dict__.update(kwargs)
        self.obj = obj
        self.on = on
        self.closed = closed
        self.window = window
        self.min_periods = min_periods
        self.center = center
        self.win_type = win_type
        self.win_freq = None
        self.axis = obj._get_axis_number(axis) if axis is not None else None
        self.validate()
        self._numba_func_cache: Dict[Optional[str], Callable] = dict()

    @property
    def _constructor(self):
        return Window

    @property
    def is_datetimelike(self) -> Optional[bool]:
        return None

    @property
    def _on(self):
        return None

    @property
    def is_freq_type(self) -> bool:
        return self.win_type == "freq"

    def validate(self) -> None:
        if self.center is not None and not is_bool(self.center):
            raise ValueError("center must be a boolean")
        if self.min_periods is not None and not is_integer(self.min_periods):
            raise ValueError("min_periods must be an integer")
        if self.closed is not None and self.closed not in [
            "right",
            "both",
            "left",
            "neither",
        ]:
            raise ValueError("closed must be 'right', 'left', 'both' or 'neither'")
        if not isinstance(self.obj, (ABCSeries, ABCDataFrame)):
            raise TypeError(f"invalid type: {type(self)}")
        if isinstance(self.window, BaseIndexer):
            self._validate_get_window_bounds_signature(self.window)

    @staticmethod
    def _validate_get_window_bounds_signature(window: BaseIndexer) -> None:
        """
        Validate that the passed BaseIndexer subclass has
        a get_window_bounds with the correct signature.
        """
        get_window_bounds_signature = inspect.signature(
            window.get_window_bounds
        ).parameters.keys()
        expected_signature = inspect.signature(
            BaseIndexer().get_window_bounds
        ).parameters.keys()
        if get_window_bounds_signature != expected_signature:
            raise ValueError(
                f"{type(window).__name__} does not implement the correct signature for "
                f"get_window_bounds"
            )

    def _create_blocks(self):
        """
        Split data into blocks & return conformed data.
        """

        obj = self._selected_obj

        # filter out the on from the object
        if self.on is not None and not isinstance(self.on, Index):
            if obj.ndim == 2:
                obj = obj.reindex(columns=obj.columns.difference([self.on]), copy=False)
        blocks = obj._to_dict_of_blocks(copy=False).values()

        return blocks, obj

    def _gotitem(self, key, ndim, subset=None):
        """
        Sub-classes to define. Return a sliced object.

        Parameters
        ----------
        key : str / list of selections
        ndim : 1,2
            requested ndim of result
        subset : object, default None
            subset to act on
        """

        # create a new object to prevent aliasing
        if subset is None:
            subset = self.obj
        self = self._shallow_copy(subset)
        self._reset_cache()
        if subset.ndim == 2:
            if is_scalar(key) and key in subset or is_list_like(key):
                self._selection = key
        return self

    def __getattr__(self, attr: str):
        if attr in self._internal_names_set:
            return object.__getattribute__(self, attr)
        if attr in self.obj:
            return self[attr]

        raise AttributeError(
            f"'{type(self).__name__}' object has no attribute '{attr}'"
        )

    def _dir_additions(self):
        return self.obj._dir_additions()

    def _get_win_type(self, kwargs: Dict):
        """
        Exists for compatibility, overriden by subclass Window.

        Parameters
        ----------
        kwargs : dict
            ignored, exists for compatibility

        Returns
        -------
        None
        """
        return None

    def _get_window(self, other=None, win_type: Optional[str] = None) -> int:
        """
        Return window length.

        Parameters
        ----------
        other :
            ignored, exists for compatibility
        win_type :
            ignored, exists for compatibility

        Returns
        -------
        window : int
        """
        if isinstance(self.window, BaseIndexer):
            return self.min_periods or 0
        return self.window

    @property
    def _window_type(self) -> str:
        return type(self).__name__

    def __repr__(self) -> str:
        """
        Provide a nice str repr of our rolling object.
        """

        attrs_list = (
            f"{attr_name}={getattr(self, attr_name)}"
            for attr_name in self._attributes
            if getattr(self, attr_name, None) is not None
        )
        attrs = ",".join(attrs_list)
        return f"{self._window_type} [{attrs}]"

    def __iter__(self):
        url = "https://github.com/pandas-dev/pandas/issues/11704"
        raise NotImplementedError(f"See issue #11704 {url}")

    def _prep_values(self, values: Optional[np.ndarray] = None) -> np.ndarray:
        """Convert input to numpy arrays for Cython routines"""
        if values is None:
            values = getattr(self._selected_obj, "values", self._selected_obj)

        # GH #12373 : rolling functions error on float32 data
        # make sure the data is coerced to float64
        if is_float_dtype(values.dtype):
            values = ensure_float64(values)
        elif is_integer_dtype(values.dtype):
            values = ensure_float64(values)
        elif needs_i8_conversion(values.dtype):
            raise NotImplementedError(
                f"ops for {self._window_type} for this "
                f"dtype {values.dtype} are not implemented"
            )
        else:
            try:
                values = ensure_float64(values)
            except (ValueError, TypeError):
                raise TypeError(f"cannot handle this type -> {values.dtype}")

        # Convert inf to nan for C funcs
        inf = np.isinf(values)
        if inf.any():
            values = np.where(inf, np.nan, values)

        return values

    def _wrap_result(self, result, block=None, obj=None):
        """
        Wrap a single result.
        """

        if obj is None:
            obj = self._selected_obj
        index = obj.index

        if isinstance(result, np.ndarray):

            if result.ndim == 1:
                from pandas import Series

                return Series(result, index, name=obj.name)

            return type(obj)(result, index=index, columns=block.columns)
        return result

    def _wrap_results(self, results, blocks, obj, exclude=None) -> FrameOrSeries:
        """
        Wrap the results.

        Parameters
        ----------
        results : list of ndarrays
        blocks : list of blocks
        obj : conformed data (may be resampled)
        exclude: list of columns to exclude, default to None
        """

        from pandas import Series, concat

        final = []
        for result, block in zip(results, blocks):

            result = self._wrap_result(result, block=block, obj=obj)
            if result.ndim == 1:
                return result
            final.append(result)

        # if we have an 'on' column
        # we want to put it back into the results
        # in the same location
        columns = self._selected_obj.columns
        if self.on is not None and not self._on.equals(obj.index):

            name = self._on.name
            final.append(Series(self._on, index=obj.index, name=name))

            if self._selection is not None:

                selection = ensure_index(self._selection)

                # need to reorder to include original location of
                # the on column (if its not already there)
                if name not in selection:
                    columns = self.obj.columns
                    indexer = columns.get_indexer(selection.tolist() + [name])
                    columns = columns.take(sorted(indexer))

        # exclude nuisance columns so that they are not reindexed
        if exclude is not None and exclude:
            columns = [c for c in columns if c not in exclude]

            if not columns:
                raise DataError("No numeric types to aggregate")

        if not len(final):
            return obj.astype("float64")
        return concat(final, axis=1).reindex(columns=columns, copy=False)

    def _center_window(self, result, window) -> np.ndarray:
        """
        Center the result in the window.
        """
        if self.axis > result.ndim - 1:
            raise ValueError("Requested axis is larger then no. of argument dimensions")

        offset = calculate_center_offset(window)
        if offset > 0:
            lead_indexer = [slice(None)] * result.ndim
            lead_indexer[self.axis] = slice(offset, None)
            result = np.copy(result[tuple(lead_indexer)])
        return result

    def _get_roll_func(self, func_name: str) -> Callable:
        """
        Wrap rolling function to check values passed.

        Parameters
        ----------
        func_name : str
            Cython function used to calculate rolling statistics

        Returns
        -------
        func : callable
        """
        window_func = getattr(window_aggregations, func_name, None)
        if window_func is None:
            raise ValueError(
                f"we do not support this function in window_aggregations.{func_name}"
            )
        return window_func

    def _get_cython_func_type(self, func: str) -> Callable:
        """
        Return a variable or fixed cython function type.

        Variable algorithms do not use window while fixed do.
        """
        if self.is_freq_type or isinstance(self.window, BaseIndexer):
            return self._get_roll_func(f"{func}_variable")
        return partial(self._get_roll_func(f"{func}_fixed"), win=self._get_window())

    def _get_window_indexer(self, window: int) -> BaseIndexer:
        """
        Return an indexer class that will compute the window start and end bounds
        """
        if isinstance(self.window, BaseIndexer):
            return self.window
        if self.is_freq_type:
            return VariableWindowIndexer(index_array=self._on.asi8, window_size=window)
        return FixedWindowIndexer(window_size=window)

    def _apply(
        self,
        func: Callable,
        center: bool,
        require_min_periods: int = 0,
        floor: int = 1,
        is_weighted: bool = False,
        name: Optional[str] = None,
        use_numba_cache: bool = False,
        **kwargs,
    ):
        """
        Rolling statistical measure using supplied function.

        Designed to be used with passed-in Cython array-based functions.

        Parameters
        ----------
        func : callable function to apply
        center : bool
        require_min_periods : int
        floor : int
        is_weighted : bool
        name : str,
            compatibility with groupby.rolling
        use_numba_cache : bool
            whether to cache a numba compiled function. Only available for numba
            enabled methods (so far only apply)
        **kwargs
            additional arguments for rolling function and window function

        Returns
        -------
        y : type of input
        """
        win_type = self._get_win_type(kwargs)
        window = self._get_window(win_type=win_type)

        blocks, obj = self._create_blocks()
        block_list = list(blocks)
        window_indexer = self._get_window_indexer(window)

        results = []
        exclude: List[Scalar] = []
        for i, b in enumerate(blocks):
            try:
                values = self._prep_values(b.values)

            except (TypeError, NotImplementedError):
                if isinstance(obj, ABCDataFrame):
                    exclude.extend(b.columns)
                    del block_list[i]
                    continue
                else:
                    raise DataError("No numeric types to aggregate")

            if values.size == 0:
                results.append(values.copy())
                continue

            # calculation function
            offset = calculate_center_offset(window) if center else 0
            additional_nans = np.array([np.nan] * offset)

            if not is_weighted:

                def calc(x):
                    x = np.concatenate((x, additional_nans))
                    if not isinstance(window, BaseIndexer):
                        min_periods = calculate_min_periods(
                            window, self.min_periods, len(x), require_min_periods, floor
                        )
                    else:
                        min_periods = calculate_min_periods(
                            self.min_periods or 1,
                            self.min_periods,
                            len(x),
                            require_min_periods,
                            floor,
                        )
                    start, end = window_indexer.get_window_bounds(
                        num_values=len(x),
                        min_periods=self.min_periods,
                        center=self.center,
                        closed=self.closed,
                    )
                    return func(x, start, end, min_periods)

            else:

                def calc(x):
                    x = np.concatenate((x, additional_nans))
                    return func(x, window, self.min_periods)

            with np.errstate(all="ignore"):
                if values.ndim > 1:
                    result = np.apply_along_axis(calc, self.axis, values)
                else:
                    result = calc(values)
                    result = np.asarray(result)

            if use_numba_cache:
                self._numba_func_cache[name] = func

            if center:
                result = self._center_window(result, window)

            results.append(result)

        return self._wrap_results(results, block_list, obj, exclude)

    def aggregate(self, func, *args, **kwargs):
        result, how = self._aggregate(func, *args, **kwargs)
        if result is None:
            return self.apply(func, raw=False, args=args, kwargs=kwargs)
        return result

    agg = aggregate

    _shared_docs["sum"] = dedent(
        """
    Calculate %(name)s sum of given DataFrame or Series.

    Parameters
    ----------
    *args, **kwargs
        For compatibility with other %(name)s methods. Has no effect
        on the computed value.

    Returns
    -------
    Series or DataFrame
        Same type as the input, with the same index, containing the
        %(name)s sum.

    See Also
    --------
    Series.sum : Reducing sum for Series.
    DataFrame.sum : Reducing sum for DataFrame.

    Examples
    --------
    >>> s = pd.Series([1, 2, 3, 4, 5])
    >>> s
    0    1
    1    2
    2    3
    3    4
    4    5
    dtype: int64

    >>> s.rolling(3).sum()
    0     NaN
    1     NaN
    2     6.0
    3     9.0
    4    12.0
    dtype: float64

    >>> s.expanding(3).sum()
    0     NaN
    1     NaN
    2     6.0
    3    10.0
    4    15.0
    dtype: float64

    >>> s.rolling(3, center=True).sum()
    0     NaN
    1     6.0
    2     9.0
    3    12.0
    4     NaN
    dtype: float64

    For DataFrame, each %(name)s sum is computed column-wise.

    >>> df = pd.DataFrame({"A": s, "B": s ** 2})
    >>> df
       A   B
    0  1   1
    1  2   4
    2  3   9
    3  4  16
    4  5  25

    >>> df.rolling(3).sum()
          A     B
    0   NaN   NaN
    1   NaN   NaN
    2   6.0  14.0
    3   9.0  29.0
    4  12.0  50.0
    """
    )

    _shared_docs["mean"] = dedent(
        """
    Calculate the %(name)s mean of the values.

    Parameters
    ----------
    *args
        Under Review.
    **kwargs
        Under Review.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the %(name)s
        calculation.

    See Also
    --------
    Series.%(name)s : Calling object with Series data.
    DataFrame.%(name)s : Calling object with DataFrames.
    Series.mean : Equivalent method for Series.
    DataFrame.mean : Equivalent method for DataFrame.

    Examples
    --------
    The below examples will show rolling mean calculations with window sizes of
    two and three, respectively.

    >>> s = pd.Series([1, 2, 3, 4])
    >>> s.rolling(2).mean()
    0    NaN
    1    1.5
    2    2.5
    3    3.5
    dtype: float64

    >>> s.rolling(3).mean()
    0    NaN
    1    NaN
    2    2.0
    3    3.0
    dtype: float64
    """
    )

    _shared_docs["var"] = dedent(
        """
    Calculate unbiased %(name)s variance.
    %(versionadded)s
    Normalized by N-1 by default. This can be changed using the `ddof`
    argument.

    Parameters
    ----------
    ddof : int, default 1
        Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.
    *args, **kwargs
        For NumPy compatibility. No additional arguments are used.

    Returns
    -------
    Series or DataFrame
        Returns the same object type as the caller of the %(name)s calculation.

    See Also
    --------
    Series.%(name)s : Calling object with Series data.
    DataFrame.%(name)s : Calling object with DataFrames.
    Series.var : Equivalent method for Series.
    DataFrame.var : Equivalent method for DataFrame.
    numpy.var : Equivalent method for Numpy array.

    Notes
    -----
    The default `ddof` of 1 used in :meth:`Series.var` is different than the
    default `ddof` of 0 in :func:`numpy.var`.

    A minimum of 1 period is required for the rolling calculation.

    Examples
    --------
    >>> s = pd.Series([5, 5, 6, 7, 5, 5, 5])
    >>> s.rolling(3).var()
    0         NaN
    1         NaN
    2    0.333333
    3    1.000000
    4    1.000000
    5    1.333333
    6    0.000000
    dtype: float64

    >>> s.expanding(3).var()
    0         NaN
    1         NaN
    2    0.333333
    3    0.916667
    4    0.800000
    5    0.700000
    6    0.619048
    dtype: float64
    """
    )

    _shared_docs["std"] = dedent(
        """
    Calculate %(name)s standard deviation.
    %(versionadded)s
    Normalized by N-1 by default. This can be changed using the `ddof`
    argument.

    Parameters
    ----------
    ddof : int, default 1
        Delta Degrees of Freedom.  The divisor used in calculations
        is ``N - ddof``, where ``N`` represents the number of elements.
    *args, **kwargs
        For NumPy compatibility. No additional arguments are used.

    Returns
    -------
    Series or DataFrame
        Returns the same object type as the caller of the %(name)s calculation.

    See Also
    --------
    Series.%(name)s : Calling object with Series data.
    DataFrame.%(name)s : Calling object with DataFrames.
    Series.std : Equivalent method for Series.
    DataFrame.std : Equivalent method for DataFrame.
    numpy.std : Equivalent method for Numpy array.

    Notes
    -----
    The default `ddof` of 1 used in Series.std is different than the default
    `ddof` of 0 in numpy.std.

    A minimum of one period is required for the rolling calculation.

    Examples
    --------
    >>> s = pd.Series([5, 5, 6, 7, 5, 5, 5])
    >>> s.rolling(3).std()
    0         NaN
    1         NaN
    2    0.577350
    3    1.000000
    4    1.000000
    5    1.154701
    6    0.000000
    dtype: float64

    >>> s.expanding(3).std()
    0         NaN
    1         NaN
    2    0.577350
    3    0.957427
    4    0.894427
    5    0.836660
    6    0.786796
    dtype: float64
    """
    )


class Window(_Window):
    """
    Provide rolling window calculations.

    Parameters
    ----------
    window : int, offset, or BaseIndexer subclass
        Size of the moving window. This is the number of observations used for
        calculating the statistic. Each window will be a fixed size.

        If its an offset then this will be the time period of each window. Each
        window will be a variable sized based on the observations included in
        the time-period. This is only valid for datetimelike indexes.

        If a BaseIndexer subclass is passed, calculates the window boundaries
        based on the defined ``get_window_bounds`` method. Additional rolling
        keyword arguments, namely `min_periods`, `center`, and
        `closed` will be passed to `get_window_bounds`.
    min_periods : int, default None
        Minimum number of observations in window required to have a value
        (otherwise result is NA). For a window that is specified by an offset,
        `min_periods` will default to 1. Otherwise, `min_periods` will default
        to the size of the window.
    center : bool, default False
        Set the labels at the center of the window.
    win_type : str, default None
        Provide a window type. If ``None``, all points are evenly weighted.
        See the notes below for further information.
    on : str, optional
        For a DataFrame, a datetime-like column or MultiIndex level on which
        to calculate the rolling window, rather than the DataFrame's index.
        Provided integer column is ignored and excluded from result since
        an integer index is not used to calculate the rolling window.
    axis : int or str, default 0
    closed : str, default None
        Make the interval closed on the 'right', 'left', 'both' or
        'neither' endpoints.
        For offset-based windows, it defaults to 'right'.
        For fixed windows, defaults to 'both'. Remaining cases not implemented
        for fixed windows.

    Returns
    -------
    a Window or Rolling sub-classed for the particular operation

    See Also
    --------
    expanding : Provides expanding transformations.
    ewm : Provides exponential weighted functions.

    Notes
    -----
    By default, the result is set to the right edge of the window. This can be
    changed to the center of the window by setting ``center=True``.

    To learn more about the offsets & frequency strings, please see `this link
    <https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html#offset-aliases>`__.

    The recognized win_types are:

    * ``boxcar``
    * ``triang``
    * ``blackman``
    * ``hamming``
    * ``bartlett``
    * ``parzen``
    * ``bohman``
    * ``blackmanharris``
    * ``nuttall``
    * ``barthann``
    * ``kaiser`` (needs beta)
    * ``gaussian`` (needs std)
    * ``general_gaussian`` (needs power, width)
    * ``slepian`` (needs width)
    * ``exponential`` (needs tau), center is set to None.

    If ``win_type=None`` all points are evenly weighted. To learn more about
    different window types see `scipy.signal window functions
    <https://docs.scipy.org/doc/scipy/reference/signal.html#window-functions>`__.

    Examples
    --------

    >>> df = pd.DataFrame({'B': [0, 1, 2, np.nan, 4]})
    >>> df
         B
    0  0.0
    1  1.0
    2  2.0
    3  NaN
    4  4.0

    Rolling sum with a window length of 2, using the 'triang'
    window type.

    >>> df.rolling(2, win_type='triang').sum()
         B
    0  NaN
    1  0.5
    2  1.5
    3  NaN
    4  NaN

    Rolling sum with a window length of 2, using the 'gaussian'
    window type (note how we need to specify std).

    >>> df.rolling(2, win_type='gaussian').sum(std=3)
              B
    0       NaN
    1  0.986207
    2  2.958621
    3       NaN
    4       NaN

    Rolling sum with a window length of 2, min_periods defaults
    to the window length.

    >>> df.rolling(2).sum()
         B
    0  NaN
    1  1.0
    2  3.0
    3  NaN
    4  NaN

    Same as above, but explicitly set the min_periods

    >>> df.rolling(2, min_periods=1).sum()
         B
    0  0.0
    1  1.0
    2  3.0
    3  2.0
    4  4.0

    A ragged (meaning not-a-regular frequency), time-indexed DataFrame

    >>> df = pd.DataFrame({'B': [0, 1, 2, np.nan, 4]},
    ...                   index = [pd.Timestamp('20130101 09:00:00'),
    ...                            pd.Timestamp('20130101 09:00:02'),
    ...                            pd.Timestamp('20130101 09:00:03'),
    ...                            pd.Timestamp('20130101 09:00:05'),
    ...                            pd.Timestamp('20130101 09:00:06')])

    >>> df
                           B
    2013-01-01 09:00:00  0.0
    2013-01-01 09:00:02  1.0
    2013-01-01 09:00:03  2.0
    2013-01-01 09:00:05  NaN
    2013-01-01 09:00:06  4.0

    Contrasting to an integer rolling window, this will roll a variable
    length window corresponding to the time period.
    The default for min_periods is 1.

    >>> df.rolling('2s').sum()
                           B
    2013-01-01 09:00:00  0.0
    2013-01-01 09:00:02  1.0
    2013-01-01 09:00:03  3.0
    2013-01-01 09:00:05  NaN
    2013-01-01 09:00:06  4.0
    """

    def validate(self):
        super().validate()

        window = self.window
        if isinstance(window, BaseIndexer):
            raise NotImplementedError(
                "BaseIndexer subclasses not implemented with win_types."
            )
        elif isinstance(window, (list, tuple, np.ndarray)):
            pass
        elif is_integer(window):
            if window <= 0:
                raise ValueError("window must be > 0 ")
            import_optional_dependency(
                "scipy", extra="Scipy is required to generate window weight."
            )
            import scipy.signal as sig

            if not isinstance(self.win_type, str):
                raise ValueError(f"Invalid win_type {self.win_type}")
            if getattr(sig, self.win_type, None) is None:
                raise ValueError(f"Invalid win_type {self.win_type}")
        else:
            raise ValueError(f"Invalid window {window}")

    def _get_win_type(self, kwargs: Dict) -> Union[str, Tuple]:
        """
        Extract arguments for the window type, provide validation for it
        and return the validated window type.

        Parameters
        ----------
        kwargs : dict

        Returns
        -------
        win_type : str, or tuple
        """
        # the below may pop from kwargs
        def _validate_win_type(win_type, kwargs):
            arg_map = {
                "kaiser": ["beta"],
                "gaussian": ["std"],
                "general_gaussian": ["power", "width"],
                "slepian": ["width"],
                "exponential": ["tau"],
            }

            if win_type in arg_map:
                win_args = _pop_args(win_type, arg_map[win_type], kwargs)
                if win_type == "exponential":
                    # exponential window requires the first arg (center)
                    # to be set to None (necessary for symmetric window)
                    win_args.insert(0, None)

                return tuple([win_type] + win_args)

            return win_type

        def _pop_args(win_type, arg_names, kwargs):
            all_args = []
            for n in arg_names:
                if n not in kwargs:
                    raise ValueError(f"{win_type} window requires {n}")
                all_args.append(kwargs.pop(n))
            return all_args

        return _validate_win_type(self.win_type, kwargs)

    def _get_window(
        self, other=None, win_type: Optional[Union[str, Tuple]] = None
    ) -> np.ndarray:
        """
        Get the window, weights.

        Parameters
        ----------
        other :
            ignored, exists for compatibility
        win_type : str, or tuple
            type of window to create

        Returns
        -------
        window : ndarray
            the window, weights
        """

        window = self.window
        if isinstance(window, (list, tuple, np.ndarray)):
            return com.asarray_tuplesafe(window).astype(float)
        elif is_integer(window):
            import scipy.signal as sig

            # GH #15662. `False` makes symmetric window, rather than periodic.
            return sig.get_window(win_type, window, False).astype(float)

    _agg_see_also_doc = dedent(
        """
    See Also
    --------
    pandas.DataFrame.rolling.aggregate
    pandas.DataFrame.aggregate
    """
    )

    _agg_examples_doc = dedent(
        """
    Examples
    --------

    >>> df = pd.DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'])
    >>> df
              A         B         C
    0 -2.385977 -0.102758  0.438822
    1 -1.004295  0.905829 -0.954544
    2  0.735167 -0.165272 -1.619346
    3 -0.702657 -1.340923 -0.706334
    4 -0.246845  0.211596 -0.901819
    5  2.463718  3.157577 -1.380906
    6 -1.142255  2.340594 -0.039875
    7  1.396598 -1.647453  1.677227
    8 -0.543425  1.761277 -0.220481
    9 -0.640505  0.289374 -1.550670

    >>> df.rolling(3, win_type='boxcar').agg('mean')
              A         B         C
    0       NaN       NaN       NaN
    1       NaN       NaN       NaN
    2 -0.885035  0.212600 -0.711689
    3 -0.323928 -0.200122 -1.093408
    4 -0.071445 -0.431533 -1.075833
    5  0.504739  0.676083 -0.996353
    6  0.358206  1.903256 -0.774200
    7  0.906020  1.283573  0.085482
    8 -0.096361  0.818139  0.472290
    9  0.070889  0.134399 -0.031308
    """
    )

    @Substitution(
        see_also=_agg_see_also_doc,
        examples=_agg_examples_doc,
        versionadded="",
        klass="Series/DataFrame",
        axis="",
    )
    @Appender(_shared_docs["aggregate"])
    def aggregate(self, func, *args, **kwargs):
        result, how = self._aggregate(func, *args, **kwargs)
        if result is None:

            # these must apply directly
            result = func(self)

        return result

    agg = aggregate

    @Substitution(name="window")
    @Appender(_shared_docs["sum"])
    def sum(self, *args, **kwargs):
        nv.validate_window_func("sum", args, kwargs)
        window_func = self._get_roll_func("roll_weighted_sum")
        window_func = get_weighted_roll_func(window_func)
        return self._apply(
            window_func, center=self.center, is_weighted=True, name="sum", **kwargs
        )

    @Substitution(name="window")
    @Appender(_shared_docs["mean"])
    def mean(self, *args, **kwargs):
        nv.validate_window_func("mean", args, kwargs)
        window_func = self._get_roll_func("roll_weighted_mean")
        window_func = get_weighted_roll_func(window_func)
        return self._apply(
            window_func, center=self.center, is_weighted=True, name="mean", **kwargs
        )

    @Substitution(name="window", versionadded="\n.. versionadded:: 1.0.0\n")
    @Appender(_shared_docs["var"])
    def var(self, ddof=1, *args, **kwargs):
        nv.validate_window_func("var", args, kwargs)
        window_func = partial(self._get_roll_func("roll_weighted_var"), ddof=ddof)
        window_func = get_weighted_roll_func(window_func)
        kwargs.pop("name", None)
        return self._apply(
            window_func, center=self.center, is_weighted=True, name="var", **kwargs
        )

    @Substitution(name="window", versionadded="\n.. versionadded:: 1.0.0\n")
    @Appender(_shared_docs["std"])
    def std(self, ddof=1, *args, **kwargs):
        nv.validate_window_func("std", args, kwargs)
        return zsqrt(self.var(ddof=ddof, name="std", **kwargs))


class _Rolling(_Window):
    @property
    def _constructor(self):
        return Rolling


class _Rolling_and_Expanding(_Rolling):

    _shared_docs["count"] = dedent(
        r"""
    The %(name)s count of any non-NaN observations inside the window.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the %(name)s
        calculation.

    See Also
    --------
    Series.%(name)s : Calling object with Series data.
    DataFrame.%(name)s : Calling object with DataFrames.
    DataFrame.count : Count of the full DataFrame.

    Examples
    --------
    >>> s = pd.Series([2, 3, np.nan, 10])
    >>> s.rolling(2).count()
    0    1.0
    1    2.0
    2    1.0
    3    1.0
    dtype: float64
    >>> s.rolling(3).count()
    0    1.0
    1    2.0
    2    2.0
    3    2.0
    dtype: float64
    >>> s.rolling(4).count()
    0    1.0
    1    2.0
    2    2.0
    3    3.0
    dtype: float64
    """
    )

    def count(self):

        blocks, obj = self._create_blocks()

        window = self._get_window()
        window = min(window, len(obj)) if not self.center else window

        results = []
        for b in blocks:
            result = b.notna().astype(int)
            result = self._constructor(
                result,
                window=window,
                min_periods=0,
                center=self.center,
                axis=self.axis,
                closed=self.closed,
            ).sum()
            results.append(result)

        return self._wrap_results(results, blocks, obj)

    _shared_docs["apply"] = dedent(
        r"""
    The %(name)s function's apply function.

    Parameters
    ----------
    func : function
        Must produce a single value from an ndarray input if ``raw=True``
        or a single value from a Series if ``raw=False``. Can also accept a
        Numba JIT function with ``engine='numba'`` specified.

        .. versionchanged:: 1.0.0

    raw : bool, default None
        * ``False`` : passes each row or column as a Series to the
          function.
        * ``True`` : the passed function will receive ndarray
          objects instead.
          If you are just applying a NumPy reduction function this will
          achieve much better performance.
    engine : str, default 'cython'
        * ``'cython'`` : Runs rolling apply through C-extensions from cython.
        * ``'numba'`` : Runs rolling apply through JIT compiled code from numba.
          Only available when ``raw`` is set to ``True``.

          .. versionadded:: 1.0.0

    engine_kwargs : dict, default None
        * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
        * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
          and ``parallel`` dictionary keys. The values must either be ``True`` or
          ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
          ``{'nopython': True, 'nogil': False, 'parallel': False}`` and will be
          applied to both the ``func`` and the ``apply`` rolling aggregation.

          .. versionadded:: 1.0.0

    args : tuple, default None
        Positional arguments to be passed into func.
    kwargs : dict, default None
        Keyword arguments to be passed into func.

    Returns
    -------
    Series or DataFrame
        Return type is determined by the caller.

    See Also
    --------
    Series.%(name)s : Series %(name)s.
    DataFrame.%(name)s : DataFrame %(name)s.

    Notes
    -----
    See :ref:`stats.rolling_apply` for extended documentation and performance
    considerations for the Numba engine.
    """
    )

    def apply(
        self,
        func,
        raw: bool = False,
        engine: str = "cython",
        engine_kwargs: Optional[Dict] = None,
        args: Optional[Tuple] = None,
        kwargs: Optional[Dict] = None,
    ):
        if args is None:
            args = ()
        if kwargs is None:
            kwargs = {}
        kwargs.pop("_level", None)
        kwargs.pop("floor", None)
        window = self._get_window()
        offset = calculate_center_offset(window) if self.center else 0
        if not is_bool(raw):
            raise ValueError("raw parameter must be `True` or `False`")

        if engine == "cython":
            if engine_kwargs is not None:
                raise ValueError("cython engine does not accept engine_kwargs")
            apply_func = self._generate_cython_apply_func(
                args, kwargs, raw, offset, func
            )
        elif engine == "numba":
            if raw is False:
                raise ValueError("raw must be `True` when using the numba engine")
            if func in self._numba_func_cache:
                # Return an already compiled version of roll_apply if available
                apply_func = self._numba_func_cache[func]
            else:
                apply_func = generate_numba_apply_func(
                    args, kwargs, func, engine_kwargs
                )
        else:
            raise ValueError("engine must be either 'numba' or 'cython'")

        # TODO: Why do we always pass center=False?
        # name=func for WindowGroupByMixin._apply
        return self._apply(
            apply_func,
            center=False,
            floor=0,
            name=func,
            use_numba_cache=engine == "numba",
        )

    def _generate_cython_apply_func(self, args, kwargs, raw, offset, func):
        from pandas import Series

        window_func = partial(
            self._get_cython_func_type("roll_generic"),
            args=args,
            kwargs=kwargs,
            raw=raw,
            offset=offset,
            func=func,
        )

        def apply_func(values, begin, end, min_periods, raw=raw):
            if not raw:
                values = Series(values, index=self.obj.index)
            return window_func(values, begin, end, min_periods)

        return apply_func

    def sum(self, *args, **kwargs):
        nv.validate_window_func("sum", args, kwargs)
        window_func = self._get_cython_func_type("roll_sum")
        kwargs.pop("floor", None)
        return self._apply(
            window_func, center=self.center, floor=0, name="sum", **kwargs
        )

    _shared_docs["max"] = dedent(
        """
    Calculate the %(name)s maximum.

    Parameters
    ----------
    *args, **kwargs
        Arguments and keyword arguments to be passed into func.
    """
    )

    def max(self, *args, **kwargs):
        nv.validate_window_func("max", args, kwargs)
        window_func = self._get_cython_func_type("roll_max")
        return self._apply(window_func, center=self.center, name="max", **kwargs)

    _shared_docs["min"] = dedent(
        """
    Calculate the %(name)s minimum.

    Parameters
    ----------
    **kwargs
        Under Review.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the %(name)s
        calculation.

    See Also
    --------
    Series.%(name)s : Calling object with a Series.
    DataFrame.%(name)s : Calling object with a DataFrame.
    Series.min : Similar method for Series.
    DataFrame.min : Similar method for DataFrame.

    Examples
    --------
    Performing a rolling minimum with a window size of 3.

    >>> s = pd.Series([4, 3, 5, 2, 6])
    >>> s.rolling(3).min()
    0    NaN
    1    NaN
    2    3.0
    3    2.0
    4    2.0
    dtype: float64
    """
    )

    def min(self, *args, **kwargs):
        nv.validate_window_func("min", args, kwargs)
        window_func = self._get_cython_func_type("roll_min")
        return self._apply(window_func, center=self.center, name="min", **kwargs)

    def mean(self, *args, **kwargs):
        nv.validate_window_func("mean", args, kwargs)
        window_func = self._get_cython_func_type("roll_mean")
        return self._apply(window_func, center=self.center, name="mean", **kwargs)

    _shared_docs["median"] = dedent(
        """
    Calculate the %(name)s median.

    Parameters
    ----------
    **kwargs
        For compatibility with other %(name)s methods. Has no effect
        on the computed median.

    Returns
    -------
    Series or DataFrame
        Returned type is the same as the original object.

    See Also
    --------
    Series.%(name)s : Calling object with Series data.
    DataFrame.%(name)s : Calling object with DataFrames.
    Series.median : Equivalent method for Series.
    DataFrame.median : Equivalent method for DataFrame.

    Examples
    --------
    Compute the rolling median of a series with a window size of 3.

    >>> s = pd.Series([0, 1, 2, 3, 4])
    >>> s.rolling(3).median()
    0    NaN
    1    NaN
    2    1.0
    3    2.0
    4    3.0
    dtype: float64
    """
    )

    def median(self, **kwargs):
        window_func = self._get_roll_func("roll_median_c")
        window_func = partial(window_func, win=self._get_window())
        return self._apply(window_func, center=self.center, name="median", **kwargs)

    def std(self, ddof=1, *args, **kwargs):
        nv.validate_window_func("std", args, kwargs)
        kwargs.pop("require_min_periods", None)
        window_func = self._get_cython_func_type("roll_var")

        def zsqrt_func(values, begin, end, min_periods):
            return zsqrt(window_func(values, begin, end, min_periods, ddof=ddof))

        # ddof passed again for compat with groupby.rolling
        return self._apply(
            zsqrt_func,
            center=self.center,
            require_min_periods=1,
            name="std",
            ddof=ddof,
            **kwargs,
        )

    def var(self, ddof=1, *args, **kwargs):
        nv.validate_window_func("var", args, kwargs)
        kwargs.pop("require_min_periods", None)
        window_func = partial(self._get_cython_func_type("roll_var"), ddof=ddof)
        # ddof passed again for compat with groupby.rolling
        return self._apply(
            window_func,
            center=self.center,
            require_min_periods=1,
            name="var",
            ddof=ddof,
            **kwargs,
        )

    _shared_docs[
        "skew"
    ] = """
    Unbiased %(name)s skewness.

    Parameters
    ----------
    **kwargs
        Keyword arguments to be passed into func.
    """

    def skew(self, **kwargs):
        window_func = self._get_cython_func_type("roll_skew")
        kwargs.pop("require_min_periods", None)
        return self._apply(
            window_func,
            center=self.center,
            require_min_periods=3,
            name="skew",
            **kwargs,
        )

    _shared_docs["kurt"] = dedent(
        """
    Calculate unbiased %(name)s kurtosis.

    This function uses Fisher's definition of kurtosis without bias.

    Parameters
    ----------
    **kwargs
        Under Review.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the %(name)s
        calculation.

    See Also
    --------
    Series.%(name)s : Calling object with Series data.
    DataFrame.%(name)s : Calling object with DataFrames.
    Series.kurt : Equivalent method for Series.
    DataFrame.kurt : Equivalent method for DataFrame.
    scipy.stats.skew : Third moment of a probability density.
    scipy.stats.kurtosis : Reference SciPy method.

    Notes
    -----
    A minimum of 4 periods is required for the %(name)s calculation.
    """
    )

    def kurt(self, **kwargs):
        window_func = self._get_cython_func_type("roll_kurt")
        kwargs.pop("require_min_periods", None)
        return self._apply(
            window_func,
            center=self.center,
            require_min_periods=4,
            name="kurt",
            **kwargs,
        )

    _shared_docs["quantile"] = dedent(
        """
    Calculate the %(name)s quantile.

    Parameters
    ----------
    quantile : float
        Quantile to compute. 0 <= quantile <= 1.
    interpolation : {'linear', 'lower', 'higher', 'midpoint', 'nearest'}
        .. versionadded:: 0.23.0

        This optional parameter specifies the interpolation method to use,
        when the desired quantile lies between two data points `i` and `j`:

            * linear: `i + (j - i) * fraction`, where `fraction` is the
              fractional part of the index surrounded by `i` and `j`.
            * lower: `i`.
            * higher: `j`.
            * nearest: `i` or `j` whichever is nearest.
            * midpoint: (`i` + `j`) / 2.
    **kwargs
        For compatibility with other %(name)s methods. Has no effect on
        the result.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the %(name)s
        calculation.

    See Also
    --------
    Series.quantile : Computes value at the given quantile over all data
        in Series.
    DataFrame.quantile : Computes values at the given quantile over
        requested axis in DataFrame.

    Examples
    --------
    >>> s = pd.Series([1, 2, 3, 4])
    >>> s.rolling(2).quantile(.4, interpolation='lower')
    0    NaN
    1    1.0
    2    2.0
    3    3.0
    dtype: float64

    >>> s.rolling(2).quantile(.4, interpolation='midpoint')
    0    NaN
    1    1.5
    2    2.5
    3    3.5
    dtype: float64
    """
    )

    def quantile(self, quantile, interpolation="linear", **kwargs):
        if quantile == 1.0:
            window_func = self._get_cython_func_type("roll_max")
        elif quantile == 0.0:
            window_func = self._get_cython_func_type("roll_min")
        else:
            window_func = partial(
                self._get_roll_func("roll_quantile"),
                win=self._get_window(),
                quantile=quantile,
                interpolation=interpolation,
            )

        # Pass through for groupby.rolling
        kwargs["quantile"] = quantile
        kwargs["interpolation"] = interpolation
        return self._apply(window_func, center=self.center, name="quantile", **kwargs)

    _shared_docs[
        "cov"
    ] = """
        Calculate the %(name)s sample covariance.

        Parameters
        ----------
        other : Series, DataFrame, or ndarray, optional
            If not supplied then will default to self and produce pairwise
            output.
        pairwise : bool, default None
            If False then only matching columns between self and other will be
            used and the output will be a DataFrame.
            If True then all pairwise combinations will be calculated and the
            output will be a MultiIndexed DataFrame in the case of DataFrame
            inputs. In the case of missing elements, only complete pairwise
            observations will be used.
        ddof : int, default 1
            Delta Degrees of Freedom.  The divisor used in calculations
            is ``N - ddof``, where ``N`` represents the number of elements.
        **kwargs
            Keyword arguments to be passed into func.
    """

    def cov(self, other=None, pairwise=None, ddof=1, **kwargs):
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)

        # GH 16058: offset window
        if self.is_freq_type:
            window = self.win_freq
        else:
            window = self._get_window(other)

        def _get_cov(X, Y):
            # GH #12373 : rolling functions error on float32 data
            # to avoid potential overflow, cast the data to float64
            X = X.astype("float64")
            Y = Y.astype("float64")
            mean = lambda x: x.rolling(
                window, self.min_periods, center=self.center
            ).mean(**kwargs)
            count = (X + Y).rolling(window=window, center=self.center).count(**kwargs)
            bias_adj = count / (count - ddof)
            return (mean(X * Y) - mean(X) * mean(Y)) * bias_adj

        return _flex_binary_moment(
            self._selected_obj, other._selected_obj, _get_cov, pairwise=bool(pairwise)
        )

    _shared_docs["corr"] = dedent(
        """
    Calculate %(name)s correlation.

    Parameters
    ----------
    other : Series, DataFrame, or ndarray, optional
        If not supplied then will default to self.
    pairwise : bool, default None
        Calculate pairwise combinations of columns within a
        DataFrame. If `other` is not specified, defaults to `True`,
        otherwise defaults to `False`.
        Not relevant for :class:`~pandas.Series`.
    **kwargs
        Unused.

    Returns
    -------
    Series or DataFrame
        Returned object type is determined by the caller of the
        %(name)s calculation.

    See Also
    --------
    Series.%(name)s : Calling object with Series data.
    DataFrame.%(name)s : Calling object with DataFrames.
    Series.corr : Equivalent method for Series.
    DataFrame.corr : Equivalent method for DataFrame.
    %(name)s.cov : Similar method to calculate covariance.
    numpy.corrcoef : NumPy Pearson's correlation calculation.

    Notes
    -----
    This function uses Pearson's definition of correlation
    (https://en.wikipedia.org/wiki/Pearson_correlation_coefficient).

    When `other` is not specified, the output will be self correlation (e.g.
    all 1's), except for :class:`~pandas.DataFrame` inputs with `pairwise`
    set to `True`.

    Function will return ``NaN`` for correlations of equal valued sequences;
    this is the result of a 0/0 division error.

    When `pairwise` is set to `False`, only matching columns between `self` and
    `other` will be used.

    When `pairwise` is set to `True`, the output will be a MultiIndex DataFrame
    with the original index on the first level, and the `other` DataFrame
    columns on the second level.

    In the case of missing elements, only complete pairwise observations
    will be used.

    Examples
    --------
    The below example shows a rolling calculation with a window size of
    four matching the equivalent function call using :meth:`numpy.corrcoef`.

    >>> v1 = [3, 3, 3, 5, 8]
    >>> v2 = [3, 4, 4, 4, 8]
    >>> # numpy returns a 2X2 array, the correlation coefficient
    >>> # is the number at entry [0][1]
    >>> print(f"{np.corrcoef(v1[:-1], v2[:-1])[0][1]:.6f}")
    0.333333
    >>> print(f"{np.corrcoef(v1[1:], v2[1:])[0][1]:.6f}")
    0.916949
    >>> s1 = pd.Series(v1)
    >>> s2 = pd.Series(v2)
    >>> s1.rolling(4).corr(s2)
    0         NaN
    1         NaN
    2         NaN
    3    0.333333
    4    0.916949
    dtype: float64

    The below example shows a similar rolling calculation on a
    DataFrame using the pairwise option.

    >>> matrix = np.array([[51., 35.], [49., 30.], [47., 32.],\
    [46., 31.], [50., 36.]])
    >>> print(np.corrcoef(matrix[:-1,0], matrix[:-1,1]).round(7))
    [[1.         0.6263001]
     [0.6263001  1.       ]]
    >>> print(np.corrcoef(matrix[1:,0], matrix[1:,1]).round(7))
    [[1.         0.5553681]
     [0.5553681  1.        ]]
    >>> df = pd.DataFrame(matrix, columns=['X','Y'])
    >>> df
          X     Y
    0  51.0  35.0
    1  49.0  30.0
    2  47.0  32.0
    3  46.0  31.0
    4  50.0  36.0
    >>> df.rolling(4).corr(pairwise=True)
                X         Y
    0 X       NaN       NaN
      Y       NaN       NaN
    1 X       NaN       NaN
      Y       NaN       NaN
    2 X       NaN       NaN
      Y       NaN       NaN
    3 X  1.000000  0.626300
      Y  0.626300  1.000000
    4 X  1.000000  0.555368
      Y  0.555368  1.000000
    """
    )

    def corr(self, other=None, pairwise=None, **kwargs):
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)
        window = self._get_window(other)

        def _get_corr(a, b):
            a = a.rolling(
                window=window, min_periods=self.min_periods, center=self.center
            )
            b = b.rolling(
                window=window, min_periods=self.min_periods, center=self.center
            )

            return a.cov(b, **kwargs) / (a.std(**kwargs) * b.std(**kwargs))

        return _flex_binary_moment(
            self._selected_obj, other._selected_obj, _get_corr, pairwise=bool(pairwise)
        )


class Rolling(_Rolling_and_Expanding):
    @cache_readonly
    def is_datetimelike(self) -> bool:
        return isinstance(
            self._on, (ABCDatetimeIndex, ABCTimedeltaIndex, ABCPeriodIndex)
        )

    @cache_readonly
    def _on(self) -> Index:
        if self.on is None:
            if self.axis == 0:
                return self.obj.index
            else:
                # i.e. self.axis == 1
                return self.obj.columns
        elif isinstance(self.on, Index):
            return self.on
        elif isinstance(self.obj, ABCDataFrame) and self.on in self.obj.columns:
            return Index(self.obj[self.on])
        else:
            raise ValueError(
                f"invalid on specified as {self.on}, "
                "must be a column (of DataFrame), an Index "
                "or None"
            )

    def validate(self):
        super().validate()

        # we allow rolling on a datetimelike index
        if (self.obj.empty or self.is_datetimelike) and isinstance(
            self.window, (str, ABCDateOffset, timedelta)
        ):

            self._validate_monotonic()
            freq = self._validate_freq()

            # we don't allow center
            if self.center:
                raise NotImplementedError(
                    "center is not implemented "
                    "for datetimelike and offset "
                    "based windows"
                )

            # this will raise ValueError on non-fixed freqs
            self.win_freq = self.window
            self.window = freq.nanos
            self.win_type = "freq"

            # min_periods must be an integer
            if self.min_periods is None:
                self.min_periods = 1

        elif isinstance(self.window, BaseIndexer):
            # Passed BaseIndexer subclass should handle all other rolling kwargs
            return
        elif not is_integer(self.window):
            raise ValueError("window must be an integer")
        elif self.window < 0:
            raise ValueError("window must be non-negative")

        if not self.is_datetimelike and self.closed is not None:
            raise ValueError(
                "closed only implemented for datetimelike and offset based windows"
            )

    def _validate_monotonic(self):
        """
        Validate monotonic (increasing or decreasing).
        """
        if not (self._on.is_monotonic_increasing or self._on.is_monotonic_decreasing):
            formatted = self.on
            if self.on is None:
                formatted = "index"
            raise ValueError(f"{formatted} must be monotonic")

    def _validate_freq(self):
        """
        Validate & return window frequency.
        """
        from pandas.tseries.frequencies import to_offset

        try:
            return to_offset(self.window)
        except (TypeError, ValueError):
            raise ValueError(
                f"passed window {self.window} is not "
                "compatible with a datetimelike "
                "index"
            )

    _agg_see_also_doc = dedent(
        """
    See Also
    --------
    Series.rolling
    DataFrame.rolling
    """
    )

    _agg_examples_doc = dedent(
        """
    Examples
    --------

    >>> df = pd.DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'])
    >>> df
              A         B         C
    0 -2.385977 -0.102758  0.438822
    1 -1.004295  0.905829 -0.954544
    2  0.735167 -0.165272 -1.619346
    3 -0.702657 -1.340923 -0.706334
    4 -0.246845  0.211596 -0.901819
    5  2.463718  3.157577 -1.380906
    6 -1.142255  2.340594 -0.039875
    7  1.396598 -1.647453  1.677227
    8 -0.543425  1.761277 -0.220481
    9 -0.640505  0.289374 -1.550670

    >>> df.rolling(3).sum()
              A         B         C
    0       NaN       NaN       NaN
    1       NaN       NaN       NaN
    2 -2.655105  0.637799 -2.135068
    3 -0.971785 -0.600366 -3.280224
    4 -0.214334 -1.294599 -3.227500
    5  1.514216  2.028250 -2.989060
    6  1.074618  5.709767 -2.322600
    7  2.718061  3.850718  0.256446
    8 -0.289082  2.454418  1.416871
    9  0.212668  0.403198 -0.093924

    >>> df.rolling(3).agg({'A':'sum', 'B':'min'})
              A         B
    0       NaN       NaN
    1       NaN       NaN
    2 -2.655105 -0.165272
    3 -0.971785 -1.340923
    4 -0.214334 -1.340923
    5  1.514216 -1.340923
    6  1.074618  0.211596
    7  2.718061 -1.647453
    8 -0.289082 -1.647453
    9  0.212668 -1.647453
    """
    )

    @Substitution(
        see_also=_agg_see_also_doc,
        examples=_agg_examples_doc,
        versionadded="",
        klass="Series/Dataframe",
        axis="",
    )
    @Appender(_shared_docs["aggregate"])
    def aggregate(self, func, *args, **kwargs):
        return super().aggregate(func, *args, **kwargs)

    agg = aggregate

    @Substitution(name="rolling")
    @Appender(_shared_docs["count"])
    def count(self):

        # different impl for freq counting
        if self.is_freq_type:
            window_func = self._get_roll_func("roll_count")
            return self._apply(window_func, center=self.center, name="count")

        return super().count()

    @Substitution(name="rolling")
    @Appender(_shared_docs["apply"])
    def apply(
        self,
        func,
        raw=False,
        engine="cython",
        engine_kwargs=None,
        args=None,
        kwargs=None,
    ):
        return super().apply(
            func,
            raw=raw,
            engine=engine,
            engine_kwargs=engine_kwargs,
            args=args,
            kwargs=kwargs,
        )

    @Substitution(name="rolling")
    @Appender(_shared_docs["sum"])
    def sum(self, *args, **kwargs):
        nv.validate_rolling_func("sum", args, kwargs)
        return super().sum(*args, **kwargs)

    @Substitution(name="rolling")
    @Appender(_doc_template)
    @Appender(_shared_docs["max"])
    def max(self, *args, **kwargs):
        nv.validate_rolling_func("max", args, kwargs)
        return super().max(*args, **kwargs)

    @Substitution(name="rolling")
    @Appender(_shared_docs["min"])
    def min(self, *args, **kwargs):
        nv.validate_rolling_func("min", args, kwargs)
        return super().min(*args, **kwargs)

    @Substitution(name="rolling")
    @Appender(_shared_docs["mean"])
    def mean(self, *args, **kwargs):
        nv.validate_rolling_func("mean", args, kwargs)
        return super().mean(*args, **kwargs)

    @Substitution(name="rolling")
    @Appender(_shared_docs["median"])
    def median(self, **kwargs):
        return super().median(**kwargs)

    @Substitution(name="rolling", versionadded="")
    @Appender(_shared_docs["std"])
    def std(self, ddof=1, *args, **kwargs):
        nv.validate_rolling_func("std", args, kwargs)
        return super().std(ddof=ddof, **kwargs)

    @Substitution(name="rolling", versionadded="")
    @Appender(_shared_docs["var"])
    def var(self, ddof=1, *args, **kwargs):
        nv.validate_rolling_func("var", args, kwargs)
        return super().var(ddof=ddof, **kwargs)

    @Substitution(name="rolling")
    @Appender(_doc_template)
    @Appender(_shared_docs["skew"])
    def skew(self, **kwargs):
        return super().skew(**kwargs)

    _agg_doc = dedent(
        """
    Examples
    --------

    The example below will show a rolling calculation with a window size of
    four matching the equivalent function call using `scipy.stats`.

    >>> arr = [1, 2, 3, 4, 999]
    >>> import scipy.stats
    >>> print(f"{scipy.stats.kurtosis(arr[:-1], bias=False):.6f}")
    -1.200000
    >>> print(f"{scipy.stats.kurtosis(arr[1:], bias=False):.6f}")
    3.999946
    >>> s = pd.Series(arr)
    >>> s.rolling(4).kurt()
    0         NaN
    1         NaN
    2         NaN
    3   -1.200000
    4    3.999946
    dtype: float64
    """
    )

    @Appender(_agg_doc)
    @Substitution(name="rolling")
    @Appender(_shared_docs["kurt"])
    def kurt(self, **kwargs):
        return super().kurt(**kwargs)

    @Substitution(name="rolling")
    @Appender(_shared_docs["quantile"])
    def quantile(self, quantile, interpolation="linear", **kwargs):
        return super().quantile(
            quantile=quantile, interpolation=interpolation, **kwargs
        )

    @Substitution(name="rolling")
    @Appender(_doc_template)
    @Appender(_shared_docs["cov"])
    def cov(self, other=None, pairwise=None, ddof=1, **kwargs):
        return super().cov(other=other, pairwise=pairwise, ddof=ddof, **kwargs)

    @Substitution(name="rolling")
    @Appender(_shared_docs["corr"])
    def corr(self, other=None, pairwise=None, **kwargs):
        return super().corr(other=other, pairwise=pairwise, **kwargs)


Rolling.__doc__ = Window.__doc__


class RollingGroupby(WindowGroupByMixin, Rolling):
    """
    Provide a rolling groupby implementation.
    """

    @property
    def _constructor(self):
        return Rolling

    def _gotitem(self, key, ndim, subset=None):

        # we are setting the index on the actual object
        # here so our index is carried thru to the selected obj
        # when we do the splitting for the groupby
        if self.on is not None:
            self._groupby.obj = self._groupby.obj.set_index(self._on)
            self.on = None
        return super()._gotitem(key, ndim, subset=subset)

    def _validate_monotonic(self):
        """
        Validate that on is monotonic;
        we don't care for groupby.rolling
        because we have already validated at a higher
        level.
        """
        pass
