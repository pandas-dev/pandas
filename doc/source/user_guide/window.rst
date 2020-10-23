.. _window:

{{ header }}

********************
Windowing Operations
********************

pandas contains a compact set of APIs for performing windowing operations - an operation that performs
an aggregation over a sliding partition of values. The API functions similarly to the ``groupby`` API
in that :class:`Series` and :class:`DataFrame` call the windowing method with
necessary parameters and then subsequently call the aggregation function.

.. ipython:: python

   s = pd.Series(range(5))
   s.rolling(window=2).sum()

The windows are comprised by looking back the length of the window from the current observation.
The result above can be derived by taking the sum of the following, windowed partitions of data:

.. ipython:: python

   for window in s.rolling(window=2):
       print(window)

Overview
--------

pandas supports 4 types of windowing operations:

# . Rolling window: Generic fixed or variable sliding window over the values.
# . Weighted window: Weighted, non-rectangular window supplied by the ``scipy.signal`` library.
# . Expanding window: Accumulating window over the values.
# . Exponentially Weighted window: Accumulating and exponentially weighted window over the values.

=============================   =================  ===========================   ===========================  ========================
Concept                         Method             Returned Object               Supports time-based windows  Supports chained groupby
=============================   =================  ===========================   ===========================  ========================
Rolling window                  ``rolling``        ``Rolling``                   Yes                          Yes
Weighted window                 ``rolling``        ``Window``                    No                           No
Expanding window                ``expanding``      ``Expanding``                 No                           Yes
Exponentially Weighted window   ``ewm``            ``ExponentialMovingWindow``   No                           No
=============================   =================  ===========================   ===========================  ========================

As noted above, some methods support specifying a window based on a timespan:

.. ipython:: python

   s = pd.Series(range(5), index=pd.date_range('2020-01-01', periods=5, freq='1D'))
   s.rolling(window='2D').sum()

Additionally, some methods support chaining a ``groupby`` operation with a windowing operation
which will first group the data by the specified keys and then perform a windowing operation per group.

.. ipython:: python

   df = pd.DataFrame({'A': ['a', 'b', 'a', 'b', 'a'], 'B': range(5))
   df.groupby('A').expanding().sum()

.. note::

   Windowing operations currently only support operation on numeric data (integer and float)
   and will always return ``float64`` values.

.. warning::

    Some windowing aggregation, ``mean``, ``sum``, ``var`` and ``std`` methods may suffer from numerical
    imprecision due to the underlying windowing algorithms accumulating sums. When values differ
    with magnitude :math:`1/np.finfo(np.double).eps` this results in truncation. It must be
    noted, that large values may have an impact on windows, which do not include these values. `Kahan summation
    <https://en.wikipedia.org/wiki/Kahan_summation_algorithm>`__ is used
    to compute the rolling sums to preserve accuracy as much as possible.


All windowing operations support a ``min_periods`` argument that dictates the minimum amount of
non-``np.nan`` values a window must have; otherwise, the resulting value is ``np.nan``.
``min_peridos`` defaults to 1 for time-based windows and ``window`` for fixed windows

.. ipython:: python

   s = pd.Series([np.nan, 1, 2, np.nan, np.nan, 3])
   # Equivalent to min_periods=3
   s.rolling(window=3, min_periods=None).sum()
   s.rolling(window=3, min_periods=2).sum()


Rolling window
--------------

Generic rolling windows support specifying windows as a fixed number of observations or variable
number of observations based on an offset. If a time based offset is provided, the corresponding
time based index must be monotonic.

.. ipython:: python

   times = ['2020-01-01', '2020-01-03', '2020-01-04', '2020-01-05, '2020-01-29']
   s = pd.Series(range(5), index=pd.DatetimeIndex(times))
   s
   # Window with 2 observations
   s.rolling(window=2).sum()
   # Window with 2 days worth of observations
   s.rolling(window='2D').sum()

For all supported aggregation functions, see :ref:`api.functions_rolling`.


Centering windows
~~~~~~~~~~~~~~~~~

By default the labels are set to the right edge of the window, but a
``center`` keyword is available so the labels can be set at the center.

.. ipython:: python

   s = pd.Series(range(10))
   s.rolling(window=5).mean()
   s.rolling(window=5, center=True).mean()


.. _window.endpoints:

Rolling window endpoints
~~~~~~~~~~~~~~~~~~~~~~~~

The inclusion of the interval endpoints in rolling window calculations can be specified with the ``closed``
parameter:

.. csv-table::
    :header: "``closed``", "Description", "Default for"
    :widths: 20, 30, 30

    ``right``, close right endpoint,
    ``left``, close left endpoint,
    ``both``, close both endpoints,
    ``neither``, open endpoints,

For example, having the right endpoint open is useful in many problems that require that there is no contamination
from present information back to past information. This allows the rolling window to compute statistics
"up to that point in time", but not including that point in time.

.. ipython:: python

   df = pd.DataFrame(
       {"x": 1},
       index=[
           pd.Timestamp("20130101 09:00:01"),
           pd.Timestamp("20130101 09:00:02"),
           pd.Timestamp("20130101 09:00:03"),
           pd.Timestamp("20130101 09:00:04"),
           pd.Timestamp("20130101 09:00:06"),
       ],
   )

   df["right"] = df.rolling("2s", closed="right").x.sum()  # default
   df["both"] = df.rolling("2s", closed="both").x.sum()
   df["left"] = df.rolling("2s", closed="left").x.sum()
   df["neither"] = df.rolling("2s", closed="neither").x.sum()

   df


.. _window.custom_rolling_window:

Custom window rolling
~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 1.0

In addition to accepting an integer or offset as a ``window`` argument, ``rolling`` also accepts
a ``BaseIndexer`` subclass that allows a user to define a custom method for calculating window bounds.
The ``BaseIndexer`` subclass will need to define a ``get_window_bounds`` method that returns
a tuple of two arrays, the first being the starting indices of the windows and second being the
ending indices of the windows. Additionally, ``num_values``, ``min_periods``, ``center``, ``closed``
and will automatically be passed to ``get_window_bounds`` and the defined method must
always accept these arguments.

For example, if we have the following ``DataFrame``:

.. ipython:: python

   use_expanding = [True, False, True, False, True]
   use_expanding
   df = pd.DataFrame({"values": range(5)})
   df

and we want to use an expanding window where ``use_expanding`` is ``True`` otherwise a window of size
1, we can create the following ``BaseIndexer`` subclass:

.. code-block:: ipython

   In [2]: from pandas.api.indexers import BaseIndexer
   ...:
   ...: class CustomIndexer(BaseIndexer):
   ...:
   ...:    def get_window_bounds(self, num_values, min_periods, center, closed):
   ...:        start = np.empty(num_values, dtype=np.int64)
   ...:        end = np.empty(num_values, dtype=np.int64)
   ...:        for i in range(num_values):
   ...:            if self.use_expanding[i]:
   ...:                start[i] = 0
   ...:                end[i] = i + 1
   ...:            else:
   ...:                start[i] = i
   ...:                end[i] = i + self.window_size
   ...:        return start, end
   ...:

   In [3]: indexer = CustomIndexer(window_size=1, use_expanding=use_expanding)

   In [4]: df.rolling(indexer).sum()
   Out[4]:
       values
   0     0.0
   1     1.0
   2     3.0
   3     3.0
   4    10.0

You can view other examples of ``BaseIndexer`` subclasses `here <https://github.com/pandas-dev/pandas/blob/master/pandas/core/window/indexers.py>`__

.. versionadded:: 1.1

One subclass of note within those examples is the ``VariableOffsetWindowIndexer`` that allows
rolling operations over a non-fixed offset like a ``BusinessDay``.

.. ipython:: python

   from pandas.api.indexers import VariableOffsetWindowIndexer

   df = pd.DataFrame(range(10), index=pd.date_range("2020", periods=10))
   offset = pd.offsets.BDay(1)
   indexer = VariableOffsetWindowIndexer(index=df.index, offset=offset)
   df
   df.rolling(indexer).sum()

For some problems knowledge of the future is available for analysis. For example, this occurs when
each data point is a full time series read from an experiment, and the task is to extract underlying
conditions. In these cases it can be useful to perform forward-looking rolling window computations.
:func:`FixedForwardWindowIndexer <pandas.api.indexers.FixedForwardWindowIndexer>` class is available for this purpose.
This :func:`BaseIndexer <pandas.api.indexers.BaseIndexer>` subclass implements a closed fixed-width
forward-looking rolling window, and we can use it as follows:

.. ipython:: ipython

   from pandas.api.indexers import FixedForwardWindowIndexer
   indexer = FixedForwardWindowIndexer(window_size=2)
   df.rolling(indexer, min_periods=1).sum()


.. _window.rolling_apply:

Rolling apply
~~~~~~~~~~~~~

The :meth:`~Rolling.apply` function takes an extra ``func`` argument and performs
generic rolling computations. The ``func`` argument should be a single function
that produces a single value from an ndarray input. ``raw`` specifies whether
the windows are cast as :class:`Series` objects ``raw=False`` or ndarray objects ``raw=True``.

.. ipython:: python

   def mad(x):
       return np.fabs(x - x.mean()).mean()

   s = pd.Series(range(10))
   s.rolling(window=4).apply(mad, raw=True)


Using the Numba engine
~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 1.0

Additionally, :meth:`~Rolling.apply` can leverage `Numba <https://numba.pydata.org/>`__
if installed as an optional dependency. The apply aggregation can be executed using Numba by specifying
``engine='numba'`` and ``engine_kwargs`` arguments (``raw`` must also be set to ``True``).
Numba will be applied in potentially two routines:

1. If ``func`` is a standard Python function, the engine will `JIT <https://numba.pydata.org/numba-doc/latest/user/overview.html>`__
the passed function. ``func`` can also be a JITed function in which case the engine will not JIT the function again.

2. The engine will JIT the for loop where the apply function is applied to each window.

The ``engine_kwargs`` argument is a dictionary of keyword arguments that will be passed into the
`numba.jit decorator <https://numba.pydata.org/numba-doc/latest/reference/jit-compilation.html#numba.jit>`__.
These keyword arguments will be applied to *both* the passed function (if a standard Python function)
and the apply for loop over each window. Currently only ``nogil``, ``nopython``, and ``parallel`` are supported,
and their default values are set to ``False``, ``True`` and ``False`` respectively.

.. note::

   In terms of performance, **the first time a function is run using the Numba engine will be slow**
   as Numba will have some function compilation overhead. However, the compiled functions are cached,
   and subsequent calls will be fast. In general, the Numba engine is performant with
   a larger amount of data points (e.g. 1+ million).

.. code-block:: ipython

   In [1]: data = pd.Series(range(1_000_000))

   In [2]: roll = data.rolling(10)

   In [3]: def f(x):
      ...:     return np.sum(x) + 5
   # Run the first time, compilation time will affect performance
   In [4]: %timeit -r 1 -n 1 roll.apply(f, engine='numba', raw=True)  # noqa: E225
   1.23 s ± 0 ns per loop (mean ± std. dev. of 1 run, 1 loop each)
   # Function is cached and performance will improve
   In [5]: %timeit roll.apply(f, engine='numba', raw=True)
   188 ms ± 1.93 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)

   In [6]: %timeit roll.apply(f, engine='cython', raw=True)
   3.92 s ± 59 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

Weighted window
---------------

The ``win_type`` argument in ``.rolling`` generates a weighted windows that are commonly used in filtering
and spectral estimation. ``win_type`` must be string that corresponds to a `scipy.signal window function
<https://docs.scipy.org/doc/scipy/reference/signal.windows.html#module-scipy.signal.windows>`__.
Scipy must be installed in order to use these windows, and supplementary arguments
that the Scipy window methods take must be specified in the aggregation function.


.. ipython:: python

   s = pd.Series(range(10))
   s.rolling(window=5).mean()
   s.rolling(window=5, win_type="triang").mean()
   # Supplementary Scipy arguments passed in the aggregation function
   s.rolling(window=5, win_type="gaussian").mean(std=0.1)

For all supported aggregation functions, see :ref:`api.functions_window`.

Expanding window
----------------

For all supported aggregation functions, see :ref:`api.functions_expanding`.

Exponentially Weighted window
-----------------------------

For all supported aggregation functions, see :ref:`api.functions_ewm`.


.. _stats.moments.binary:

Binary window functions
~~~~~~~~~~~~~~~~~~~~~~~

:meth:`~Rolling.cov` and :meth:`~Rolling.corr` can compute moving window statistics about
two ``Series`` or any combination of ``DataFrame/Series`` or
``DataFrame/DataFrame``. Here is the behavior in each case:

* two ``Series``: compute the statistic for the pairing.
* ``DataFrame/Series``: compute the statistics for each column of the DataFrame
  with the passed Series, thus returning a DataFrame.
* ``DataFrame/DataFrame``: by default compute the statistic for matching column
  names, returning a DataFrame. If the keyword argument ``pairwise=True`` is
  passed then computes the statistic for each pair of columns, returning a
  ``MultiIndexed DataFrame`` whose ``index`` are the dates in question (see :ref:`the next section
  <stats.moments.corr_pairwise>`).

For example:

.. ipython:: python

   df = pd.DataFrame(
       np.random.randn(10, 4),
       index=pd.date_range("2020-01-01", periods=10),
       columns=["A", "B", "C", "D"],
   )
   df = df.cumsum()

   df2 = df[:4]
   df2.rolling(window=2).corr(df2["B"])

.. _stats.moments.corr_pairwise:

Computing rolling pairwise covariances and correlations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In financial data analysis and other fields it's common to compute covariance
and correlation matrices for a collection of time series. Often one is also
interested in moving-window covariance and correlation matrices. This can be
done by passing the ``pairwise`` keyword argument, which in the case of
``DataFrame`` inputs will yield a MultiIndexed ``DataFrame`` whose ``index`` are the dates in
question. In the case of a single DataFrame argument the ``pairwise`` argument
can even be omitted:

.. note::

    Missing values are ignored and each entry is computed using the pairwise
    complete observations.  Please see the :ref:`covariance section
    <computation.covariance>` for :ref:`caveats
    <computation.covariance.caveats>` associated with this method of
    calculating covariance and correlation matrices.

.. ipython:: python

   covs = (
       df[["B", "C", "D"]]
       .rolling(window=50)
       .cov(df[["A", "B", "C"]], pairwise=True)
   )
   covs.loc["2002-09-22":]

.. ipython:: python

   correls = df.rolling(window=50).corr()
   correls.loc["2002-09-22":]

You can efficiently retrieve the time series of correlations between two
columns by reshaping and indexing:

.. ipython:: python
   :suppress:

   plt.close("all")

.. ipython:: python

   @savefig rolling_corr_pairwise_ex.png
   correls.unstack(1)[("A", "C")].plot()

.. _stats.aggregate:

Aggregation
-----------

Once the ``Rolling``, ``Expanding`` or ``ExponentialMovingWindow`` objects have been created, several methods are available to
perform multiple computations on the data. These operations are similar to the :ref:`aggregating API <basics.aggregate>`,
:ref:`groupby API <groupby.aggregate>`, and :ref:`resample API <timeseries.aggregate>`.


.. ipython:: python

   dfa = pd.DataFrame(
       np.random.randn(1000, 3),
       index=pd.date_range("1/1/2000", periods=1000),
       columns=["A", "B", "C"],
   )
   r = dfa.rolling(window=60, min_periods=1)
   r

We can aggregate by passing a function to the entire DataFrame, or select a
Series (or multiple Series) via standard ``__getitem__``.

.. ipython:: python

   r.aggregate(np.sum)

   r["A"].aggregate(np.sum)

   r[["A", "B"]].aggregate(np.sum)

As you can see, the result of the aggregation will have the selected columns, or all
columns if none are selected.

.. _stats.aggregate.multifunc:

Applying multiple functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

With windowed ``Series`` you can also pass a list of functions to do
aggregation with, outputting a DataFrame:

.. ipython:: python

   r["A"].agg([np.sum, np.mean, np.std])

On a windowed DataFrame, you can pass a list of functions to apply to each
column, which produces an aggregated result with a hierarchical index:

.. ipython:: python

   r.agg([np.sum, np.mean])

Passing a dict of functions has different behavior by default, see the next
section.

Applying different functions to DataFrame columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By passing a dict to ``aggregate`` you can apply a different aggregation to the
columns of a ``DataFrame``:

.. ipython:: python

   r.agg({"A": np.sum, "B": lambda x: np.std(x, ddof=1)})

The function names can also be strings. In order for a string to be valid it
must be implemented on the windowed object

.. ipython:: python

   r.agg({"A": "sum", "B": "std"})

Furthermore you can pass a nested dict to indicate different aggregations on different columns.

.. ipython:: python

   r.agg({"A": ["sum", "std"], "B": ["mean", "std"]})


.. _stats.moments.expanding:

Expanding windows
-----------------

A common alternative to rolling statistics is to use an *expanding* window,
which yields the value of the statistic with all the data available up to that
point in time.

These follow a similar interface to ``.rolling``, with the ``.expanding`` method
returning an :class:`~pandas.core.window.Expanding` object.

As these calculations are a special case of rolling statistics,
they are implemented in pandas such that the following two calls are equivalent:

.. ipython:: python

   df.rolling(window=len(df), min_periods=1).mean()[:5]

   df.expanding(min_periods=1).mean()[:5]

These have a similar set of methods to ``.rolling`` methods.

Method summary
~~~~~~~~~~~~~~

.. currentmodule:: pandas.core.window

.. csv-table::
    :header: "Function", "Description"
    :widths: 20, 80

    :meth:`~Expanding.count`, Number of non-null observations
    :meth:`~Expanding.sum`, Sum of values
    :meth:`~Expanding.mean`, Mean of values
    :meth:`~Expanding.median`, Arithmetic median of values
    :meth:`~Expanding.min`, Minimum
    :meth:`~Expanding.max`, Maximum
    :meth:`~Expanding.std`, Sample standard deviation
    :meth:`~Expanding.var`, Sample variance
    :meth:`~Expanding.skew`, Sample skewness (3rd moment)
    :meth:`~Expanding.kurt`, Sample kurtosis (4th moment)
    :meth:`~Expanding.quantile`, Sample quantile (value at %)
    :meth:`~Expanding.apply`, Generic apply
    :meth:`~Expanding.cov`, Sample covariance (binary)
    :meth:`~Expanding.corr`, Sample correlation (binary)
    :meth:`~Expanding.sem`, Standard error of mean

.. note::

   Using sample variance formulas for :meth:`~Expanding.std` and
   :meth:`~Expanding.var` comes with the same caveats as using them with rolling
   windows. See :ref:`this section <computation.window_variance.caveats>` for more
   information.

   The same caveats apply to using any supported statistical sample methods.

.. currentmodule:: pandas

Aside from not having a ``window`` parameter, these functions have the same
interfaces as their ``.rolling`` counterparts. Like above, the parameters they
all accept are:

* ``min_periods``: threshold of non-null data points to require. Defaults to
  minimum needed to compute statistic. No ``NaNs`` will be output once
  ``min_periods`` non-null data points have been seen.
* ``center``: boolean, whether to set the labels at the center (default is False).

.. _stats.moments.expanding.note:
.. note::

   The output of the ``.rolling`` and ``.expanding`` methods do not return a
   ``NaN`` if there are at least ``min_periods`` non-null values in the current
   window. For example:

   .. ipython:: python

        sn = pd.Series([1, 2, np.nan, 3, np.nan, 4])
        sn
        sn.rolling(2).max()
        sn.rolling(2, min_periods=1).max()

   In case of expanding functions, this differs from :meth:`~DataFrame.cumsum`,
   :meth:`~DataFrame.cumprod`, :meth:`~DataFrame.cummax`,
   and :meth:`~DataFrame.cummin`, which return ``NaN`` in the output wherever
   a ``NaN`` is encountered in the input. In order to match the output of ``cumsum``
   with ``expanding``, use :meth:`~DataFrame.fillna`:

   .. ipython:: python

        sn.expanding().sum()
        sn.cumsum()
        sn.cumsum().fillna(method="ffill")


An expanding window statistic will be more stable (and less responsive) than
its rolling window counterpart as the increasing window size decreases the
relative impact of an individual data point. As an example, here is the
:meth:`~core.window.Expanding.mean` output for the previous time series dataset:

.. ipython:: python
   :suppress:

   plt.close("all")

.. ipython:: python

   s.plot(style="k--")

   @savefig expanding_mean_frame.png
   s.expanding().mean().plot(style="k")


.. _stats.moments.exponentially_weighted:

Exponentially weighted windows
------------------------------

.. currentmodule:: pandas.core.window

A related set of functions are exponentially weighted versions of several of
the above statistics. A similar interface to ``.rolling`` and ``.expanding`` is accessed
through the ``.ewm`` method to receive an :class:`~ExponentialMovingWindow` object.
A number of expanding EW (exponentially weighted)
methods are provided:


.. csv-table::
    :header: "Function", "Description"
    :widths: 20, 80

    :meth:`~ExponentialMovingWindow.mean`, EW moving average
    :meth:`~ExponentialMovingWindow.var`, EW moving variance
    :meth:`~ExponentialMovingWindow.std`, EW moving standard deviation
    :meth:`~ExponentialMovingWindow.corr`, EW moving correlation
    :meth:`~ExponentialMovingWindow.cov`, EW moving covariance

In general, a weighted moving average is calculated as

.. math::

    y_t = \frac{\sum_{i=0}^t w_i x_{t-i}}{\sum_{i=0}^t w_i},

where :math:`x_t` is the input, :math:`y_t` is the result and the :math:`w_i`
are the weights.

The EW functions support two variants of exponential weights.
The default, ``adjust=True``, uses the weights :math:`w_i = (1 - \alpha)^i`
which gives

.. math::

    y_t = \frac{x_t + (1 - \alpha)x_{t-1} + (1 - \alpha)^2 x_{t-2} + ...
    + (1 - \alpha)^t x_{0}}{1 + (1 - \alpha) + (1 - \alpha)^2 + ...
    + (1 - \alpha)^t}

When ``adjust=False`` is specified, moving averages are calculated as

.. math::

    y_0 &= x_0 \\
    y_t &= (1 - \alpha) y_{t-1} + \alpha x_t,

which is equivalent to using weights

.. math::

    w_i = \begin{cases}
        \alpha (1 - \alpha)^i & \text{if } i < t \\
        (1 - \alpha)^i        & \text{if } i = t.
    \end{cases}

.. note::

   These equations are sometimes written in terms of :math:`\alpha' = 1 - \alpha`, e.g.

   .. math::

      y_t = \alpha' y_{t-1} + (1 - \alpha') x_t.

The difference between the above two variants arises because we are
dealing with series which have finite history. Consider a series of infinite
history, with ``adjust=True``:

.. math::

    y_t = \frac{x_t + (1 - \alpha)x_{t-1} + (1 - \alpha)^2 x_{t-2} + ...}
    {1 + (1 - \alpha) + (1 - \alpha)^2 + ...}

Noting that the denominator is a geometric series with initial term equal to 1
and a ratio of :math:`1 - \alpha` we have

.. math::

    y_t &= \frac{x_t + (1 - \alpha)x_{t-1} + (1 - \alpha)^2 x_{t-2} + ...}
    {\frac{1}{1 - (1 - \alpha)}}\\
    &= [x_t + (1 - \alpha)x_{t-1} + (1 - \alpha)^2 x_{t-2} + ...] \alpha \\
    &= \alpha x_t + [(1-\alpha)x_{t-1} + (1 - \alpha)^2 x_{t-2} + ...]\alpha \\
    &= \alpha x_t + (1 - \alpha)[x_{t-1} + (1 - \alpha) x_{t-2} + ...]\alpha\\
    &= \alpha x_t + (1 - \alpha) y_{t-1}

which is the same expression as ``adjust=False`` above and therefore
shows the equivalence of the two variants for infinite series.
When ``adjust=False``, we have :math:`y_0 = x_0` and
:math:`y_t = \alpha x_t + (1 - \alpha) y_{t-1}`.
Therefore, there is an assumption that :math:`x_0` is not an ordinary value
but rather an exponentially weighted moment of the infinite series up to that
point.

One must have :math:`0 < \alpha \leq 1`, and while it is possible to pass
:math:`\alpha` directly, it's often easier to think about either the
**span**, **center  of mass (com)** or **half-life** of an EW moment:

.. math::

   \alpha =
    \begin{cases}
        \frac{2}{s + 1},               & \text{for span}\ s \geq 1\\
        \frac{1}{1 + c},               & \text{for center of mass}\ c \geq 0\\
        1 - \exp^{\frac{\log 0.5}{h}}, & \text{for half-life}\ h > 0
    \end{cases}

One must specify precisely one of **span**, **center of mass**, **half-life**
and **alpha** to the EW functions:

* **Span** corresponds to what is commonly called an "N-day EW moving average".
* **Center of mass** has a more physical interpretation and can be thought of
  in terms of span: :math:`c = (s - 1) / 2`.
* **Half-life** is the period of time for the exponential weight to reduce to
  one half.
* **Alpha** specifies the smoothing factor directly.

.. versionadded:: 1.1.0

You can also specify ``halflife`` in terms of a timedelta convertible unit to specify the amount of
time it takes for an observation to decay to half its value when also specifying a sequence
of ``times``.

.. ipython:: python

    df = pd.DataFrame({"B": [0, 1, 2, np.nan, 4]})
    df
    times = ["2020-01-01", "2020-01-03", "2020-01-10", "2020-01-15", "2020-01-17"]
    df.ewm(halflife="4 days", times=pd.DatetimeIndex(times)).mean()

The following formula is used to compute exponentially weighted mean with an input vector of times:

.. math::

    y_t = \frac{\sum_{i=0}^t 0.5^\frac{t_{t} - t_{i}}{\lambda} x_{t-i}}{0.5^\frac{t_{t} - t_{i}}{\lambda}},

Here is an example for a univariate time series:

.. ipython:: python

   s.plot(style="k--")

   @savefig ewma_ex.png
   s.ewm(span=20).mean().plot(style="k")

ExponentialMovingWindow has a ``min_periods`` argument, which has the same
meaning it does for all the ``.expanding`` and ``.rolling`` methods:
no output values will be set until at least ``min_periods`` non-null values
are encountered in the (expanding) window.

ExponentialMovingWindow also has an ``ignore_na`` argument, which determines how
intermediate null values affect the calculation of the weights.
When ``ignore_na=False`` (the default), weights are calculated based on absolute
positions, so that intermediate null values affect the result.
When ``ignore_na=True``,
weights are calculated by ignoring intermediate null values.
For example, assuming ``adjust=True``, if ``ignore_na=False``, the weighted
average of ``3, NaN, 5`` would be calculated as

.. math::

	\frac{(1-\alpha)^2 \cdot 3 + 1 \cdot 5}{(1-\alpha)^2 + 1}.

Whereas if ``ignore_na=True``, the weighted average would be calculated as

.. math::

	\frac{(1-\alpha) \cdot 3 + 1 \cdot 5}{(1-\alpha) + 1}.

The :meth:`~Ewm.var`, :meth:`~Ewm.std`, and :meth:`~Ewm.cov` functions have a ``bias`` argument,
specifying whether the result should contain biased or unbiased statistics.
For example, if ``bias=True``, ``ewmvar(x)`` is calculated as
``ewmvar(x) = ewma(x**2) - ewma(x)**2``;
whereas if ``bias=False`` (the default), the biased variance statistics
are scaled by debiasing factors

.. math::

    \frac{\left(\sum_{i=0}^t w_i\right)^2}{\left(\sum_{i=0}^t w_i\right)^2 - \sum_{i=0}^t w_i^2}.

(For :math:`w_i = 1`, this reduces to the usual :math:`N / (N - 1)` factor,
with :math:`N = t + 1`.)
See `Weighted Sample Variance <https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance>`__
on Wikipedia for further details.
