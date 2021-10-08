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
The result above can be derived by taking the sum of the following windowed partitions of data:

.. ipython:: python

   for window in s.rolling(window=2):
       print(window)


.. _window.overview:

Overview
--------

pandas supports 4 types of windowing operations:

#. Rolling window: Generic fixed or variable sliding window over the values.
#. Weighted window: Weighted, non-rectangular window supplied by the ``scipy.signal`` library.
#. Expanding window: Accumulating window over the values.
#. Exponentially Weighted window: Accumulating and exponentially weighted window over the values.

=============================   =================  ===========================   ===========================  ========================  ===================================  ===========================
Concept                         Method             Returned Object               Supports time-based windows  Supports chained groupby  Supports table method                Supports online operations
=============================   =================  ===========================   ===========================  ========================  ===================================  ===========================
Rolling window                  ``rolling``        ``Rolling``                   Yes                          Yes                       Yes (as of version 1.3)              No
Weighted window                 ``rolling``        ``Window``                    No                           No                        No                                   No
Expanding window                ``expanding``      ``Expanding``                 No                           Yes                       Yes (as of version 1.3)              No
Exponentially Weighted window   ``ewm``            ``ExponentialMovingWindow``   No                           Yes (as of version 1.2)   No                                   Yes (as of version 1.3)
=============================   =================  ===========================   ===========================  ========================  ===================================  ===========================

As noted above, some operations support specifying a window based on a time offset:

.. ipython:: python

   s = pd.Series(range(5), index=pd.date_range('2020-01-01', periods=5, freq='1D'))
   s.rolling(window='2D').sum()

Additionally, some methods support chaining a ``groupby`` operation with a windowing operation
which will first group the data by the specified keys and then perform a windowing operation per group.

.. ipython:: python

   df = pd.DataFrame({'A': ['a', 'b', 'a', 'b', 'a'], 'B': range(5)})
   df.groupby('A').expanding().sum()

.. note::

   Windowing operations currently only support numeric data (integer and float)
   and will always return ``float64`` values.

.. warning::

    Some windowing aggregation, ``mean``, ``sum``, ``var`` and ``std`` methods may suffer from numerical
    imprecision due to the underlying windowing algorithms accumulating sums. When values differ
    with magnitude :math:`1/np.finfo(np.double).eps` this results in truncation. It must be
    noted, that large values may have an impact on windows, which do not include these values. `Kahan summation
    <https://en.wikipedia.org/wiki/Kahan_summation_algorithm>`__ is used
    to compute the rolling sums to preserve accuracy as much as possible.


.. versionadded:: 1.3.0

Some windowing operations also support the ``method='table'`` option in the constructor which
performs the windowing operation over an entire :class:`DataFrame` instead of a single column or row at a time.
This can provide a useful performance benefit for a :class:`DataFrame` with many columns or rows
(with the corresponding ``axis`` argument) or the ability to utilize other columns during the windowing
operation. The ``method='table'`` option can only be used if ``engine='numba'`` is specified
in the corresponding method call.

For example, a `weighted mean <https://en.wikipedia.org/wiki/Weighted_arithmetic_mean>`__ calculation can
be calculated with :meth:`~Rolling.apply` by specifying a separate column of weights.

.. ipython:: python

   def weighted_mean(x):
       arr = np.ones((1, x.shape[1]))
       arr[:, :2] = (x[:, :2] * x[:, 2]).sum(axis=0) / x[:, 2].sum()
       return arr

   df = pd.DataFrame([[1, 2, 0.6], [2, 3, 0.4], [3, 4, 0.2], [4, 5, 0.7]])
   df.rolling(2, method="table", min_periods=0).apply(weighted_mean, raw=True, engine="numba")  # noqa:E501

.. versionadded:: 1.3

Some windowing operations also support an ``online`` method after constructing a windowing object
which returns a new object that supports passing in new :class:`DataFrame` or :class:`Series` objects
to continue the windowing calculation with the new values (i.e. online calculations).

The methods on this new windowing objects must call the aggregation method first to "prime" the initial
state of the online calculation. Then, new :class:`DataFrame` or :class:`Series` objects can be passed in
the ``update`` argument to continue the windowing calculation.

.. ipython:: python

   df = pd.DataFrame([[1, 2, 0.6], [2, 3, 0.4], [3, 4, 0.2], [4, 5, 0.7]])
   df.ewm(0.5).mean()

.. ipython:: python

   online_ewm = df.head(2).ewm(0.5).online()
   online_ewm.mean()
   online_ewm.mean(update=df.tail(1))

All windowing operations support a ``min_periods`` argument that dictates the minimum amount of
non-``np.nan`` values a window must have; otherwise, the resulting value is ``np.nan``.
``min_periods`` defaults to 1 for time-based windows and ``window`` for fixed windows

.. ipython:: python

   s = pd.Series([np.nan, 1, 2, np.nan, np.nan, 3])
   s.rolling(window=3, min_periods=1).sum()
   s.rolling(window=3, min_periods=2).sum()
   # Equivalent to min_periods=3
   s.rolling(window=3, min_periods=None).sum()


Additionally, all windowing operations supports the ``aggregate`` method for returning a result
of multiple aggregations applied to a window.

.. ipython:: python

   df = pd.DataFrame({"A": range(5), "B": range(10, 15)})
   df.expanding().agg([np.sum, np.mean, np.std])


.. _window.generic:

Rolling window
--------------

Generic rolling windows support specifying windows as a fixed number of observations or variable
number of observations based on an offset. If a time based offset is provided, the corresponding
time based index must be monotonic.

.. ipython:: python

   times = ['2020-01-01', '2020-01-03', '2020-01-04', '2020-01-05', '2020-01-29']
   s = pd.Series(range(5), index=pd.DatetimeIndex(times))
   s
   # Window with 2 observations
   s.rolling(window=2).sum()
   # Window with 2 days worth of observations
   s.rolling(window='2D').sum()

For all supported aggregation functions, see :ref:`api.functions_rolling`.

.. _window.center:

Centering windows
~~~~~~~~~~~~~~~~~

By default the labels are set to the right edge of the window, but a
``center`` keyword is available so the labels can be set at the center.

.. ipython:: python

   s = pd.Series(range(10))
   s.rolling(window=5).mean()
   s.rolling(window=5, center=True).mean()


This can also be applied to datetime-like indices.

.. versionadded:: 1.3.0

.. ipython:: python

    df = pd.DataFrame(
        {"A": [0, 1, 2, 3, 4]}, index=pd.date_range("2020", periods=5, freq="1D")
    )
    df
    df.rolling("2D", center=False).mean()
    df.rolling("2D", center=True).mean()


.. _window.endpoints:

Rolling window endpoints
~~~~~~~~~~~~~~~~~~~~~~~~

The inclusion of the interval endpoints in rolling window calculations can be specified with the ``closed``
parameter:

=============  ====================
Value          Behavior
=============  ====================
``'right'``     close right endpoint
``'left'``     close left endpoint
``'both'``     close both endpoints
``'neither'``  open endpoints
=============  ====================

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

For example, if we have the following :class:`DataFrame`

.. ipython:: python

   use_expanding = [True, False, True, False, True]
   use_expanding
   df = pd.DataFrame({"values": range(5)})
   df

and we want to use an expanding window where ``use_expanding`` is ``True`` otherwise a window of size
1, we can create the following ``BaseIndexer`` subclass:

.. code-block:: ipython

   In [2]: from pandas.api.indexers import BaseIndexer

   In [3]: class CustomIndexer(BaseIndexer):
      ...:     def get_window_bounds(self, num_values, min_periods, center, closed):
      ...:         start = np.empty(num_values, dtype=np.int64)
      ...:         end = np.empty(num_values, dtype=np.int64)
      ...:         for i in range(num_values):
      ...:             if self.use_expanding[i]:
      ...:                 start[i] = 0
      ...:                 end[i] = i + 1
      ...:             else:
      ...:                 start[i] = i
      ...:                 end[i] = i + self.window_size
      ...:         return start, end

   In [4]: indexer = CustomIndexer(window_size=1, use_expanding=use_expanding)

   In [5]: df.rolling(indexer).sum()
   Out[5]:
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

.. ipython:: python

   from pandas.api.indexers import FixedForwardWindowIndexer
   indexer = FixedForwardWindowIndexer(window_size=2)
   df.rolling(indexer, min_periods=1).sum()

We can also achieve this by using slicing, applying rolling aggregation, and then flipping the result as shown in example below:

.. ipython:: python

   df = pd.DataFrame(
       data=[
           [pd.Timestamp("2018-01-01 00:00:00"), 100],
           [pd.Timestamp("2018-01-01 00:00:01"), 101],
           [pd.Timestamp("2018-01-01 00:00:03"), 103],
           [pd.Timestamp("2018-01-01 00:00:04"), 111],
       ],
       columns=["time", "value"],
   ).set_index("time")
   df

   reversed_df = df[::-1].rolling("2s").sum()[::-1]
   reversed_df

.. _window.rolling_apply:

Rolling apply
~~~~~~~~~~~~~

The :meth:`~Rolling.apply` function takes an extra ``func`` argument and performs
generic rolling computations. The ``func`` argument should be a single function
that produces a single value from an ndarray input. ``raw`` specifies whether
the windows are cast as :class:`Series` objects (``raw=False``) or ndarray objects (``raw=True``).

.. ipython:: python

   def mad(x):
       return np.fabs(x - x.mean()).mean()

   s = pd.Series(range(10))
   s.rolling(window=4).apply(mad, raw=True)

.. _window.numba_engine:

Numba engine
~~~~~~~~~~~~

.. versionadded:: 1.0

Additionally, :meth:`~Rolling.apply` can leverage `Numba <https://numba.pydata.org/>`__
if installed as an optional dependency. The apply aggregation can be executed using Numba by specifying
``engine='numba'`` and ``engine_kwargs`` arguments (``raw`` must also be set to ``True``).
See :ref:`enhancing performance with Numba <enhancingperf.numba>` for general usage of the arguments and performance considerations.

Numba will be applied in potentially two routines:

#. If ``func`` is a standard Python function, the engine will `JIT <https://numba.pydata.org/numba-doc/latest/user/overview.html>`__ the passed function. ``func`` can also be a JITed function in which case the engine will not JIT the function again.
#. The engine will JIT the for loop where the apply function is applied to each window.

The ``engine_kwargs`` argument is a dictionary of keyword arguments that will be passed into the
`numba.jit decorator <https://numba.pydata.org/numba-doc/latest/reference/jit-compilation.html#numba.jit>`__.
These keyword arguments will be applied to *both* the passed function (if a standard Python function)
and the apply for loop over each window.

.. versionadded:: 1.3.0

``mean``, ``median``, ``max``, ``min``, and ``sum`` also support the ``engine`` and ``engine_kwargs`` arguments.

.. _window.cov_corr:

Binary window functions
~~~~~~~~~~~~~~~~~~~~~~~

:meth:`~Rolling.cov` and :meth:`~Rolling.corr` can compute moving window statistics about
two :class:`Series` or any combination of :class:`DataFrame`/:class:`Series` or
:class:`DataFrame`/:class:`DataFrame`. Here is the behavior in each case:

* two :class:`Series`: compute the statistic for the pairing.
* :class:`DataFrame`/:class:`Series`: compute the statistics for each column of the DataFrame
  with the passed Series, thus returning a DataFrame.
* :class:`DataFrame`/:class:`DataFrame`: by default compute the statistic for matching column
  names, returning a DataFrame. If the keyword argument ``pairwise=True`` is
  passed then computes the statistic for each pair of columns, returning a :class:`DataFrame` with a
  :class:`MultiIndex` whose values are the dates in question (see :ref:`the next section
  <window.corr_pairwise>`).

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

.. _window.corr_pairwise:

Computing rolling pairwise covariances and correlations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In financial data analysis and other fields it's common to compute covariance
and correlation matrices for a collection of time series. Often one is also
interested in moving-window covariance and correlation matrices. This can be
done by passing the ``pairwise`` keyword argument, which in the case of
:class:`DataFrame` inputs will yield a MultiIndexed :class:`DataFrame` whose ``index`` are the dates in
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
       .rolling(window=4)
       .cov(df[["A", "B", "C"]], pairwise=True)
   )
   covs


.. _window.weighted:

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

.. _window.expanding:

Expanding window
----------------

An expanding window yields the value of an aggregation statistic with all the data available up to that
point in time. Since these calculations are a special case of rolling statistics,
they are implemented in pandas such that the following two calls are equivalent:

.. ipython:: python

   df = pd.DataFrame(range(5))
   df.rolling(window=len(df), min_periods=1).mean()
   df.expanding(min_periods=1).mean()

For all supported aggregation functions, see :ref:`api.functions_expanding`.


.. _window.exponentially_weighted:

Exponentially Weighted window
-----------------------------

An exponentially weighted window is similar to an expanding window but with each prior point
being exponentially weighted down relative to the current point.

In general, a weighted moving average is calculated as

.. math::

    y_t = \frac{\sum_{i=0}^t w_i x_{t-i}}{\sum_{i=0}^t w_i},

where :math:`x_t` is the input, :math:`y_t` is the result and the :math:`w_i`
are the weights.

For all supported aggregation functions, see :ref:`api.functions_ewm`.

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

    y_t = \frac{\sum_{i=0}^t 0.5^\frac{t_{t} - t_{i}}{\lambda} x_{t-i}}{\sum_{i=0}^t 0.5^\frac{t_{t} - t_{i}}{\lambda}},


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
