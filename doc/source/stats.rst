.. currentmodule:: pandas.stats.api

.. _stats:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   import pandas.util.testing as tm
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   import matplotlib.pyplot as plt
   plt.close('all')

**********************************
Moving window stats and regression
**********************************

.. currentmodule:: pandas

.. currentmodule:: pandas.stats.api

.. _stats.moments:

Moving (rolling) statistics / moments
-------------------------------------

For working with time series data, a number of functions are provided for
computing common *moving* or *rolling* statistics. Among these are count, sum,
mean, median, correlation, variance, covariance, standard deviation, skewness,
and kurtosis. All of these methods are in the :mod:`pandas` namespace, but
otherwise they can be found in :mod:`pandas.stats.moments`.

.. csv-table::
    :header: "Function", "Description"
    :widths: 20, 80

    ``rolling_count``, Number of non-null observations
    ``rolling_sum``, Sum of values
    ``rolling_mean``, Mean of values
    ``rolling_median``, Arithmetic median of values
    ``rolling_min``, Minimum
    ``rolling_max``, Maximum
    ``rolling_std``, Unbiased standard deviation
    ``rolling_var``, Unbiased variance
    ``rolling_skew``, Unbiased skewness (3rd moment)
    ``rolling_kurt``, Unbiased kurtosis (4th moment)
    ``rolling_quantile``, Sample quantile (value at %)
    ``rolling_apply``, Generic apply
    ``rolling_cov``, Unbiased covariance (binary)
    ``rolling_corr``, Correlation (binary)

Generally these methods all have the same interface. The binary operators
(e.g. ``rolling_corr``) take two Series or DataFrames. Otherwise, they all
accept the following arguments:

  - ``window``: size of moving window
  - ``min_periods``: threshold of non-null data points to require (otherwise
    result is NA)
  - ``time_rule``: optionally specify a :ref:`time rule <timeseries.timerule>`
    to pre-conform the data to

These functions can be applied to ndarrays or Series objects:

.. ipython:: python

   ts = Series(randn(1000), index=DateRange('1/1/2000', periods=1000))
   ts = ts.cumsum()

   ts.plot(style='k--')

   @savefig rolling_mean_ex.png width=4.5in
   rolling_mean(ts, 60).plot(style='k')

They can also be applied to DataFrame objects. This is really just syntactic
sugar for applying the moving window operator to all of the DataFrame's columns:

.. ipython:: python
   :suppress:

   plt.close('all')

.. ipython:: python

   df = DataFrame(randn(1000, 4), index=ts.index,
                  columns=['A', 'B', 'C', 'D'])
   df = df.cumsum()

   @savefig rolling_mean_frame.png width=4.5in
   rolling_sum(df, 60).plot(subplots=True)

Exponentially weighted moment functions
---------------------------------------

A related set of functions are exponentially weighted versions of many of the
above statistics. A number of EW (exponentially weighted) functions are
provided using the blending method. For example, where :math:`y_t` is the
result and :math:`x_t` the input, we compute an exponentially weighted moving
average as

.. math::

    y_t = (1-\alpha) y_{t-1} + \alpha x_t

One must have :math:`0 < \alpha \leq 1`, but rather than pass :math:`\alpha`
directly, it's easier to think about either the **span** or **center of mass
(com)** of an EW moment:

.. math::

   \alpha =
    \begin{cases}
	\frac{2}{s + 1}, s = \text{span}\\
	\frac{1}{c + 1}, c = \text{center of mass}
    \end{cases}

You can pass one or the other to these functions but not both. **Span**
corresponds to what is commonly called a "20-day EW moving average" for
example. **Center of mass** has a more physical interpretation. For example,
**span** = 20 corresponds to **com** = 9.5. Here is the list of functions
available:

.. csv-table::
    :header: "Function", "Description"
    :widths: 20, 80

    ``ewma``, EW moving average
    ``ewvar``, EW moving variance
    ``ewstd``, EW moving standard deviation
    ``ewmcorr``, EW moving correlation
    ``ewmcov``, EW moving covariance

Here are an example for a univariate time series:

.. ipython:: python

   plt.close('all')
   ts.plot(style='k--')

   @savefig ewma_ex.png width=4.5in
   ewma(ts, span=20).plot(style='k')

.. _stats.ols:

Linear and panel regression
---------------------------

.. note::

   We plan to move this functionality to `statsmodels
   <http://statsmodels.sourceforge.net>`__ for the next release. Some of the
   result attributes may change names in order to foster naming consistency
   with the rest of statsmodels. We will provide every effort to provide
   compatibility with older versions of pandas, however.

We have implemented a very fast set of *moving-window linear regression*
classes in pandas. Two different types of regressions are supported:

  - Standard ordinary least squares (OLS) multiple regression
  - Multiple regression (OLS-based) on `panel data
    <http://en.wikipedia.org/wiki/Panel_data>`__ including with fixed-effects
    (also known as entity or individual effects) or time-effects.

Both kinds of linear models are accessed through the ``ols`` function in the
pandas namespace. They all take the following arguments to specify either a
static (full sample) or dynamic (moving window) regression:

  - ``window_type``: ``'full sample'`` (default), ``'expanding'``, or
    ``rolling``
  - ``window``: size of the moving window in the ``window_type='rolling'``
    case. If ``window`` is specified, ``window_type`` will be automatically set
    to ``'rolling'``
  - ``min_periods``: minimum number of time periods to require to compute the
    regression coefficients

Generally speaking, the ``ols`` works by being given a ``y`` (response) object
and an ``x`` (predictors) object. These can take many forms:

  - ``y``: a Series, ndarray, or DataFrame (panel model)
  - ``x``: Series, DataFrame, dict of Series, dict of DataFrame, Panel,
    LongPanel

Based on the types of ``y`` and ``x``, the model will be inferred to either a
panel model or a regular linear model. If the ``y`` variable is a DataFrame,
the result will be a panel model. In this case, the ``x`` variable must either
be a Panel, LongPanel, or a dict of DataFrame (which will be coerced into a
Panel).

Standard OLS regression
~~~~~~~~~~~~~~~~~~~~~~~

Let's pull in some sample data:

.. ipython:: python

   from pandas.io.data import DataReader
   symbols = ['MSFT', 'GOOG', 'AAPL']
   data = dict((sym, DataReader(sym, "yahoo"))
               for sym in symbols)
   panel = Panel(data).swapaxes('items', 'minor')
   close_px = panel['close']

   # convert closing prices to returns
   rets = close_px / close_px.shift(1) - 1
   rets.info()

Let's do a static regression of ``AAPL`` returns on ``GOOG`` returns:

.. ipython:: python

   model = ols(y=rets['AAPL'], x=rets.ix[:, ['GOOG']])
   model
   model.beta

If we had passed a Series instead of a DataFrame with the single ``GOOG``
column, the model would have assigned the generic name ``x`` to the sole
right-hand side variable.

We can do a moving window regression to see how the relationship changes over
time:

.. ipython:: python
   :suppress:

   plt.close('all')

.. ipython:: python

   model = ols(y=rets['AAPL'], x=rets.ix[:, ['GOOG']],
               window=250)

   # just plot the coefficient for GOOG
   @savefig moving_lm_ex.png width=5in
   model.beta['GOOG'].plot()

It looks like there are some outliers rolling in and out of the window in the
above regression, influencing the results. We could perform a simple
`winsorization <http://en.wikipedia.org/wiki/Winsorising>`__ at the 3 STD level
to trim the impact of outliers:

.. ipython:: python
   :suppress:

   plt.close('all')

.. ipython:: python

   winz = rets.copy()
   std_1year = rolling_std(rets, 250, min_periods=20)

   # cap at 3 * 1 year standard deviation
   cap_level = 3 * np.sign(winz) * std_1year
   winz[np.abs(winz) > 3 * std_1year] = cap_level

   winz_model = ols(y=winz['AAPL'], x=winz.ix[:, ['GOOG']],
               window=250)

   model.beta['GOOG'].plot(label="With outliers")

   @savefig moving_lm_winz.png width=5in
   winz_model.beta['GOOG'].plot(label="Winsorized"); plt.legend(loc='best')

So in this simple example we see the impact of winsorization is actually quite
significant. Note the correlation after winsorization remains high:

.. ipython:: python

   winz.corrwith(rets)

Multiple regressions can be run by passing a DataFrame with multiple columns
for the predictors ``x``:

.. ipython:: python

   ols(y=winz['AAPL'], x=winz.drop(['AAPL'], axis=1))

Panel regression
~~~~~~~~~~~~~~~~

We've implemented moving window panel regression on potentially unbalanced
panel data (see `this article <http://en.wikipedia.org/wiki/Panel_data>`__ if
this means nothing to you). Suppose we wanted to model the relationship between
the magnitude of the daily return and trading volume among a group of stocks,
and we want to pool all the data together to run one big regression. This is
actually quite easy:

.. ipython:: python

   # make the units somewhat comparable
   volume = panel['volume'] / 1e8
   model = ols(y=volume, x={'return' : np.abs(rets)})
   model

In a panel model, we can insert dummy (0-1) variables for the "entities"
involved (here, each of the stocks) to account the a entity-specific effect
(intercept):

.. ipython:: python

   fe_model = ols(y=volume, x={'return' : np.abs(rets)},
                  entity_effects=True)
   fe_model

Because we ran the regression with an intercept, one of the dummy variables
must be dropped or the design matrix will not be full rank. If we do not use an
intercept, all of the dummy variables will be included:

.. ipython:: python

   fe_model = ols(y=volume, x={'return' : np.abs(rets)},
                  entity_effects=True, intercept=False)
   fe_model

We can also include *time effects*, which demeans the data cross-sectionally at
each point in time (equivalent to including dummy variables for each
date). More mathematical care must be taken to properly compute the standard
errors in this case:

.. ipython:: python

   te_model = ols(y=volume, x={'return' : np.abs(rets)},
                  time_effects=True, entity_effects=True)
   te_model

Here the intercept (the mean term) is dropped by default because it will be 0
according to the model assumptions, having subtracted off the group means.

Result fields and tests
~~~~~~~~~~~~~~~~~~~~~~~

We'll leave it to the user to explore the docstrings and source, especially as
we'll be moving this code into statsmodels in the near future.
