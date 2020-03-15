from textwrap import dedent

import numpy as np

import pandas._libs.window.aggregations as window_aggregations
from pandas.compat.numpy import function as nv
from pandas.util._decorators import Appender, Substitution

from pandas.core.dtypes.generic import ABCDataFrame

from pandas.core.base import DataError
from pandas.core.window.common import (
    _doc_template,
    _get_center_of_mass,
    _shared_docs,
    zsqrt,
)
from pandas.core.window.rolling import _flex_binary_moment, _Rolling

_bias_template = """
        Parameters
        ----------
        bias : bool, default False
            Use a standard estimation bias correction.
        *args, **kwargs
            Arguments and keyword arguments to be passed into func.
"""


class EWM(_Rolling):
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
    halflife : float, optional
        Specify decay in terms of half-life,
        :math:`\alpha = 1 - \exp\left(-\ln(2) / halflife\right)`, for
        :math:`halflife > 0`.
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
    axis : {0 or 'index', 1 or 'columns'}, default 0
        The axis to use. The value 0 identifies the rows, and 1
        identifies the columns.

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
    """

    _attributes = ["com", "min_periods", "adjust", "ignore_na", "axis"]

    def __init__(
        self,
        obj,
        com=None,
        span=None,
        halflife=None,
        alpha=None,
        min_periods=0,
        adjust=True,
        ignore_na=False,
        axis=0,
    ):
        self.obj = obj
        self.com = _get_center_of_mass(com, span, halflife, alpha)
        self.min_periods = min_periods
        self.adjust = adjust
        self.ignore_na = ignore_na
        self.axis = axis
        self.on = None

    @property
    def _constructor(self):
        return EWM

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

    >>> df.ewm(alpha=0.5).mean()
              A         B         C
    0 -2.385977 -0.102758  0.438822
    1 -1.464856  0.569633 -0.490089
    2 -0.207700  0.149687 -1.135379
    3 -0.471677 -0.645305 -0.906555
    4 -0.355635 -0.203033 -0.904111
    5  1.076417  1.503943 -1.146293
    6 -0.041654  1.925562 -0.588728
    7  0.680292  0.132049  0.548693
    8  0.067236  0.948257  0.163353
    9 -0.286980  0.618493 -0.694496
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

    def _apply(self, func, **kwargs):
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
        blocks, obj = self._create_blocks()
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

            # if we have a string function name, wrap it
            if isinstance(func, str):
                cfunc = getattr(window_aggregations, func, None)
                if cfunc is None:
                    raise ValueError(
                        f"we do not support this function in window_aggregations.{func}"
                    )

                def func(arg):
                    return cfunc(
                        arg,
                        self.com,
                        int(self.adjust),
                        int(self.ignore_na),
                        int(self.min_periods),
                    )

            results.append(np.apply_along_axis(func, self.axis, values))

        return self._wrap_results(results, block_list, obj, exclude)

    @Substitution(name="ewm")
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
        return self._apply("ewma", **kwargs)

    @Substitution(name="ewm")
    @Appender(_doc_template)
    @Appender(_bias_template)
    def std(self, bias=False, *args, **kwargs):
        """
        Exponential weighted moving stddev.
        """
        nv.validate_window_func("std", args, kwargs)
        return zsqrt(self.var(bias=bias, **kwargs))

    vol = std

    @Substitution(name="ewm")
    @Appender(_doc_template)
    @Appender(_bias_template)
    def var(self, bias=False, *args, **kwargs):
        """
        Exponential weighted moving variance.
        """
        nv.validate_window_func("var", args, kwargs)

        def f(arg):
            return window_aggregations.ewmcov(
                arg,
                arg,
                self.com,
                int(self.adjust),
                int(self.ignore_na),
                int(self.min_periods),
                int(bias),
            )

        return self._apply(f, **kwargs)

    @Substitution(name="ewm")
    @Appender(_doc_template)
    def cov(self, other=None, pairwise=None, bias=False, **kwargs):
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
                int(self.adjust),
                int(self.ignore_na),
                int(self.min_periods),
                int(bias),
            )
            return X._wrap_result(cov)

        return _flex_binary_moment(
            self._selected_obj, other._selected_obj, _get_cov, pairwise=bool(pairwise)
        )

    @Substitution(name="ewm")
    @Appender(_doc_template)
    def corr(self, other=None, pairwise=None, **kwargs):
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
                    x,
                    y,
                    self.com,
                    int(self.adjust),
                    int(self.ignore_na),
                    int(self.min_periods),
                    1,
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
