import datetime
from functools import partial
from textwrap import dedent
from typing import Optional, Union

import numpy as np

from pandas._libs.tslibs import Timedelta
import pandas._libs.window.aggregations as window_aggregations
from pandas._typing import FrameOrSeries, TimedeltaConvertibleTypes
from pandas.compat.numpy import function as nv
from pandas.util._decorators import Appender, Substitution, doc

from pandas.core.dtypes.common import is_datetime64_ns_dtype
from pandas.core.dtypes.generic import ABCDataFrame

from pandas.core.base import DataError
import pandas.core.common as common
from pandas.core.window.common import _doc_template, _shared_docs, zsqrt
from pandas.core.window.rolling import _flex_binary_moment, _Rolling

_bias_template = """
        Parameters
        ----------
        bias : bool, default False
            Use a standard estimation bias correction.
        *args, **kwargs
            Arguments and keyword arguments to be passed into func.
"""


def get_center_of_mass(
    comass: Optional[float],
    span: Optional[float],
    halflife: Optional[float],
    alpha: Optional[float],
) -> float:
    valid_count = common.count_not_none(comass, span, halflife, alpha)
    if valid_count > 1:
        raise ValueError("comass, span, halflife, and alpha are mutually exclusive")

    # Convert to center of mass; domain checks ensure 0 < alpha <= 1
    if comass is not None:
        if comass < 0:
            raise ValueError("comass must satisfy: comass >= 0")
    elif span is not None:
        if span < 1:
            raise ValueError("span must satisfy: span >= 1")
        comass = (span - 1) / 2.0
    elif halflife is not None:
        if halflife <= 0:
            raise ValueError("halflife must satisfy: halflife > 0")
        decay = 1 - np.exp(np.log(0.5) / halflife)
        comass = 1 / decay - 1
    elif alpha is not None:
        if alpha <= 0 or alpha > 1:
            raise ValueError("alpha must satisfy: 0 < alpha <= 1")
        comass = (1.0 - alpha) / alpha
    else:
        raise ValueError("Must pass one of comass, span, halflife, or alpha")

    return float(comass)


class ExponentialMovingWindow(_Rolling):
    r"""
    Provide exponential weighted (EW) functions.

    Available EW functions: ``mean()``, ``var()``, ``std()``, ``corr()``, ``cov()``.

    Exactly one parameter: ``com``, ``span``, ``halflife``, or ``alpha`` must be
    provided.

    Parameters
    ----------
    com : float, optional
        Specify decay in terms of center of mass,
        :math:`\alpha = 1 / (1 + com)`, for :math:`com \geq 0`.
    span : float, optional
        Specify decay in terms of span,
        :math:`\alpha = 2 / (span + 1)`, for :math:`span \geq 1`.
    halflife : float, str, timedelta, optional
        Specify decay in terms of half-life,
        :math:`\alpha = 1 - \exp\left(-\ln(2) / halflife\right)`, for
        :math:`halflife > 0`.

        If ``times`` is specified, the time unit (str or timedelta) over which an
        observation decays to half its value. Only applicable to ``mean()``
        and halflife value will not apply to the other functions.

        .. versionadded:: 1.1.0

    alpha : float, optional
        Specify smoothing factor :math:`\alpha` directly,
        :math:`0 < \alpha \leq 1`.
    min_periods : int, default 0
        Minimum number of observations in window required to have a value
        (otherwise result is NA).
    adjust : bool, default True
        Divide by decaying adjustment factor in beginning periods to account
        for imbalance in relative weightings (viewing EWMA as a moving average).

        - When ``adjust=True`` (default), the EW function is calculated using weights
          :math:`w_i = (1 - \alpha)^i`. For example, the EW moving average of the series
          [:math:`x_0, x_1, ..., x_t`] would be:

        .. math::
            y_t = \frac{x_t + (1 - \alpha)x_{t-1} + (1 - \alpha)^2 x_{t-2} + ... + (1 -
            \alpha)^t x_0}{1 + (1 - \alpha) + (1 - \alpha)^2 + ... + (1 - \alpha)^t}

        - When ``adjust=False``, the exponentially weighted function is calculated
          recursively:

        .. math::
            \begin{split}
                y_0 &= x_0\\
                y_t &= (1 - \alpha) y_{t-1} + \alpha x_t,
            \end{split}
    ignore_na : bool, default False
        Ignore missing values when calculating weights; specify ``True`` to reproduce
        pre-0.15.0 behavior.

        - When ``ignore_na=False`` (default), weights are based on absolute positions.
          For example, the weights of :math:`x_0` and :math:`x_2` used in calculating
          the final weighted average of [:math:`x_0`, None, :math:`x_2`] are
          :math:`(1-\alpha)^2` and :math:`1` if ``adjust=True``, and
          :math:`(1-\alpha)^2` and :math:`\alpha` if ``adjust=False``.

        - When ``ignore_na=True`` (reproducing pre-0.15.0 behavior), weights are based
          on relative positions. For example, the weights of :math:`x_0` and :math:`x_2`
          used in calculating the final weighted average of
          [:math:`x_0`, None, :math:`x_2`] are :math:`1-\alpha` and :math:`1` if
          ``adjust=True``, and :math:`1-\alpha` and :math:`\alpha` if ``adjust=False``.
    axis : {0, 1}, default 0
        The axis to use. The value 0 identifies the rows, and 1
        identifies the columns.
    times : str, np.ndarray, Series, default None

        .. versionadded:: 1.1.0

        Times corresponding to the observations. Must be monotonically increasing and
        ``datetime64[ns]`` dtype.

        If str, the name of the column in the DataFrame representing the times.

        If 1-D array like, a sequence with the same shape as the observations.

        Only applicable to ``mean()``.

    Returns
    -------
    DataFrame
        A Window sub-classed for the particular operation.

    See Also
    --------
    rolling : Provides rolling window calculations.
    expanding : Provides expanding transformations.

    Notes
    -----

    More details can be found at:
    :ref:`Exponentially weighted windows <stats.moments.exponentially_weighted>`.

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

    >>> df.ewm(com=0.5).mean()
              B
    0  0.000000
    1  0.750000
    2  1.615385
    3  1.615385
    4  3.670213

    Specifying ``times`` with a timedelta ``halflife`` when computing mean.

    >>> times = ['2020-01-01', '2020-01-03', '2020-01-10', '2020-01-15', '2020-01-17']
    >>> df.ewm(halflife='4 days', times=pd.DatetimeIndex(times)).mean()
              B
    0  0.000000
    1  0.585786
    2  1.523889
    3  1.523889
    4  3.233686
    """

    _attributes = ["com", "min_periods", "adjust", "ignore_na", "axis"]

    def __init__(
        self,
        obj,
        com: Optional[float] = None,
        span: Optional[float] = None,
        halflife: Optional[Union[float, TimedeltaConvertibleTypes]] = None,
        alpha: Optional[float] = None,
        min_periods: int = 0,
        adjust: bool = True,
        ignore_na: bool = False,
        axis: int = 0,
        times: Optional[Union[str, np.ndarray, FrameOrSeries]] = None,
    ):
        self.com: Optional[float]
        self.obj = obj
        self.min_periods = max(int(min_periods), 1)
        self.adjust = adjust
        self.ignore_na = ignore_na
        self.axis = axis
        self.on = None
        if times is not None:
            if isinstance(times, str):
                times = self._selected_obj[times]
            if not is_datetime64_ns_dtype(times):
                raise ValueError("times must be datetime64[ns] dtype.")
            if len(times) != len(obj):
                raise ValueError("times must be the same length as the object.")
            if not isinstance(halflife, (str, datetime.timedelta)):
                raise ValueError(
                    "halflife must be a string or datetime.timedelta object"
                )
            self.times = np.asarray(times.astype(np.int64))
            self.halflife = Timedelta(halflife).value
            # Halflife is no longer applicable when calculating COM
            # But allow COM to still be calculated if the user passes other decay args
            if common.count_not_none(com, span, alpha) > 0:
                self.com = get_center_of_mass(com, span, None, alpha)
            else:
                self.com = None
        else:
            if halflife is not None and isinstance(halflife, (str, datetime.timedelta)):
                raise ValueError(
                    "halflife can only be a timedelta convertible argument if "
                    "times is not None."
                )
            self.times = None
            self.halflife = None
            self.com = get_center_of_mass(com, span, halflife, alpha)

    @property
    def _constructor(self):
        return ExponentialMovingWindow

    _agg_see_also_doc = dedent(
        """
    See Also
    --------
    pandas.DataFrame.rolling.aggregate
    """
    )

    _agg_examples_doc = dedent(
        """
    Examples
    --------
    >>> df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6], "C": [7, 8, 9]})
    >>> df
       A  B  C
    0  1  4  7
    1  2  5  8
    2  3  6  9

    >>> df.ewm(alpha=0.5).mean()
              A         B         C
    0  1.000000  4.000000  7.000000
    1  1.666667  4.666667  7.666667
    2  2.428571  5.428571  8.428571
    """
    )

    @doc(
        _shared_docs["aggregate"],
        see_also=_agg_see_also_doc,
        examples=_agg_examples_doc,
        versionadded="",
        klass="Series/Dataframe",
        axis="",
    )
    def aggregate(self, func, *args, **kwargs):
        return super().aggregate(func, *args, **kwargs)

    agg = aggregate

    def _apply(self, func):
        """
        Rolling statistical measure using supplied function. Designed to be
        used with passed-in Cython array-based functions.

        Parameters
        ----------
        func : str/callable to apply

        Returns
        -------
        y : same type as input argument
        """
        blocks, obj = self._create_blocks(self._selected_obj)
        block_list = list(blocks)

        results = []
        exclude = []
        for i, b in enumerate(blocks):
            try:
                values = self._prep_values(b.values)

            except (TypeError, NotImplementedError) as err:
                if isinstance(obj, ABCDataFrame):
                    exclude.extend(b.columns)
                    del block_list[i]
                    continue
                else:
                    raise DataError("No numeric types to aggregate") from err

            if values.size == 0:
                results.append(values.copy())
                continue

            results.append(np.apply_along_axis(func, self.axis, values))

        return self._wrap_results(results, block_list, obj, exclude)

    @Substitution(name="ewm", func_name="mean")
    @Appender(_doc_template)
    def mean(self, *args, **kwargs):
        """
        Exponential weighted moving average.

        Parameters
        ----------
        *args, **kwargs
            Arguments and keyword arguments to be passed into func.
        """
        nv.validate_window_func("mean", args, kwargs)
        if self.times is not None:
            window_func = self._get_roll_func("ewma_time")
            window_func = partial(
                window_func,
                minp=self.min_periods,
                times=self.times,
                halflife=self.halflife,
            )
        else:
            window_func = self._get_roll_func("ewma")
            window_func = partial(
                window_func,
                com=self.com,
                adjust=self.adjust,
                ignore_na=self.ignore_na,
                minp=self.min_periods,
            )
        return self._apply(window_func)

    @Substitution(name="ewm", func_name="std")
    @Appender(_doc_template)
    @Appender(_bias_template)
    def std(self, bias: bool = False, *args, **kwargs):
        """
        Exponential weighted moving stddev.
        """
        nv.validate_window_func("std", args, kwargs)
        return zsqrt(self.var(bias=bias, **kwargs))

    vol = std

    @Substitution(name="ewm", func_name="var")
    @Appender(_doc_template)
    @Appender(_bias_template)
    def var(self, bias: bool = False, *args, **kwargs):
        """
        Exponential weighted moving variance.
        """
        nv.validate_window_func("var", args, kwargs)

        def f(arg):
            return window_aggregations.ewmcov(
                arg, arg, self.com, self.adjust, self.ignore_na, self.min_periods, bias,
            )

        return self._apply(f)

    @Substitution(name="ewm", func_name="cov")
    @Appender(_doc_template)
    def cov(
        self,
        other: Optional[Union[np.ndarray, FrameOrSeries]] = None,
        pairwise: Optional[bool] = None,
        bias: bool = False,
        **kwargs,
    ):
        """
        Exponential weighted sample covariance.

        Parameters
        ----------
        other : Series, DataFrame, or ndarray, optional
            If not supplied then will default to self and produce pairwise
            output.
        pairwise : bool, default None
            If False then only matching columns between self and other will be
            used and the output will be a DataFrame.
            If True then all pairwise combinations will be calculated and the
            output will be a MultiIndex DataFrame in the case of DataFrame
            inputs. In the case of missing elements, only complete pairwise
            observations will be used.
        bias : bool, default False
            Use a standard estimation bias correction.
        **kwargs
           Keyword arguments to be passed into func.
        """
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)

        def _get_cov(X, Y):
            X = self._shallow_copy(X)
            Y = self._shallow_copy(Y)
            cov = window_aggregations.ewmcov(
                X._prep_values(),
                Y._prep_values(),
                self.com,
                self.adjust,
                self.ignore_na,
                self.min_periods,
                bias,
            )
            return X._wrap_result(cov)

        return _flex_binary_moment(
            self._selected_obj, other._selected_obj, _get_cov, pairwise=bool(pairwise)
        )

    @Substitution(name="ewm", func_name="corr")
    @Appender(_doc_template)
    def corr(
        self,
        other: Optional[Union[np.ndarray, FrameOrSeries]] = None,
        pairwise: Optional[bool] = None,
        **kwargs,
    ):
        """
        Exponential weighted sample correlation.

        Parameters
        ----------
        other : Series, DataFrame, or ndarray, optional
            If not supplied then will default to self and produce pairwise
            output.
        pairwise : bool, default None
            If False then only matching columns between self and other will be
            used and the output will be a DataFrame.
            If True then all pairwise combinations will be calculated and the
            output will be a MultiIndex DataFrame in the case of DataFrame
            inputs. In the case of missing elements, only complete pairwise
            observations will be used.
        **kwargs
           Keyword arguments to be passed into func.
        """
        if other is None:
            other = self._selected_obj
            # only default unset
            pairwise = True if pairwise is None else pairwise
        other = self._shallow_copy(other)

        def _get_corr(X, Y):
            X = self._shallow_copy(X)
            Y = self._shallow_copy(Y)

            def _cov(x, y):
                return window_aggregations.ewmcov(
                    x, y, self.com, self.adjust, self.ignore_na, self.min_periods, 1,
                )

            x_values = X._prep_values()
            y_values = Y._prep_values()
            with np.errstate(all="ignore"):
                cov = _cov(x_values, y_values)
                x_var = _cov(x_values, x_values)
                y_var = _cov(y_values, y_values)
                corr = cov / zsqrt(x_var * y_var)
            return X._wrap_result(corr)

        return _flex_binary_moment(
            self._selected_obj, other._selected_obj, _get_corr, pairwise=bool(pairwise)
        )
