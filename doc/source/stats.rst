.. currentmodule:: pandas.stats.api

.. _stats:

**********************************
Built-in statistical functionality
**********************************

.. currentmodule:: pandas

.. currentmodule:: pandas.stats.api

.. _stats.moments:

Moving (rolling) statistics / moments
-------------------------------------

For working with time series data, a number of functions are provided
for computing common *moving* or *rolling* statistics. Among these are
count, sum, mean, median, correlation, variance, covariance, standard
deviation, skewness, and kurtosis. All of these methods are in the
:mod:`pandas` namespace, but otherwise they can be found in
:mod:`pandas.stats.moments`.

Each of these methods observes the same interface (with relevant
methods accepting two Series arguments instead of one):

.. function:: Unary moment functions

   :Parameters:
       **arg** : ndarray, Series, DataFrame, etc.
	       If a DataFrame is passed, function will be applied to **columns**

       **window** : int
           Number of periods to include in window

       **min_periods** : int or None
	       Number of periods to require to compute a value (defaults to
	       **window**)

       **time_rule** : string or DateOffset
            Frequency to pre-convert data to

.. function:: Binary moment functions

   :Parameters:
       **arg1** : ndarray, Series
	       If a DataFrame is passed, function will be applied to **columns**

       **arg2** : ndarray, Series
	       Must be same type as **arg1**

       **window** : int
           Number of periods to include in window

       **min_periods** : int or None
	       Number of periods to require to compute a value (defaults to
	       **window**)

       **time_rule** : string or DateOffset
            Frequency to pre-convert data to

::

    >>> ts
    2000-01-31 00:00:00    -0.550139282247
    2000-02-01 00:00:00    0.0950636484432
    2000-02-02 00:00:00    0.0621763420914
    2000-02-03 00:00:00    0.125698607137
    2000-02-04 00:00:00    0.222288320816
    2000-02-07 00:00:00    0.903314747152
    2000-02-08 00:00:00    -0.391449402196
    2000-02-09 00:00:00    -0.726137553115
    2000-02-10 00:00:00    -0.89302167539
    2000-02-11 00:00:00    0.228509179513

    >>> rolling_sum(ts, 5, min_periods=3)
    2000-01-31 00:00:00    NaN
    2000-02-01 00:00:00    NaN
    2000-02-02 00:00:00    -0.0913037710365
    2000-02-03 00:00:00    0.798752592168
    2000-02-04 00:00:00    1.39432346651
    2000-02-07 00:00:00    2.44074916551
    2000-02-08 00:00:00    2.77458564938
    2000-02-09 00:00:00    1.87181399193
    2000-02-10 00:00:00    2.48549563273
    2000-02-11 00:00:00    1.81285272663

If passing a DataFrame argument, the statistics will be applied independently to
the columns:

::

    >>> df
			   A              B              C              D
    2000-01-31 00:00:00    NaN            NaN            0.03752        -0.3952
    2000-02-01 00:00:00    NaN            NaN            -1.511         -0.1126
    2000-02-02 00:00:00    1.136          NaN            0.777          -0.3502
    2000-02-03 00:00:00    0.8901         NaN            1.196          0.7456
    2000-02-04 00:00:00    0.5956         0.7684         0.9042         0.4984
    2000-02-07 00:00:00    -0.3502        1.015          0.5366         0.6628
    2000-02-08 00:00:00    0.5036         1.825          0.8682         -1.69
    2000-02-09 00:00:00    0.2327         -0.3899        0.4493         -0.1267
    2000-02-10 00:00:00    1.504          0.3904         -0.06148       1.717
    2000-02-11 00:00:00    -0.07707       0.2286         -1.039         0.1438

    >>> rolling_mean(df, 5, min_periods=3)
			   A              B              C              D
    2000-01-31 00:00:00    NaN            NaN            NaN            NaN
    2000-02-01 00:00:00    NaN            NaN            NaN            NaN
    2000-02-02 00:00:00    NaN            NaN            -0.2321        -0.286
    2000-02-03 00:00:00    NaN            NaN            0.125          -0.02811
    2000-02-04 00:00:00    0.8737         NaN            0.2809         0.07718
    2000-02-07 00:00:00    0.5677         NaN            0.3807         0.2888
    2000-02-08 00:00:00    0.5549         1.203          0.8565         -0.0267
    2000-02-09 00:00:00    0.3744         0.8047         0.7909         0.018
    2000-02-10 00:00:00    0.4971         0.7219         0.5394         0.2123
    2000-02-11 00:00:00    0.3626         0.6139         0.1507         0.1414

Each of these methods can optionally accept a **time_rule** argument
(see :ref:`time rules <datetools.timerules>`) which is provided as a
convenience when the user wishes to guarantee that the window of the
statistic.

Here are some plots of the unary moment functions:

.. plot:: plots/stats/moments_rolling.py

And the binary moment functions:

.. plot:: plots/stats/moments_rolling_binary.py

Exponentially weighted moment functions
---------------------------------------

It's also quite common to want to do non-equally weighted moving statistics,
such as exponentially weighted (EW) moving average or EW moving standard
deviation. A number of EW functions are provided using the blending method. For
example, where :math:`y_t` is the result and :math:`x_t` the input, we compute
an exponentially weighted moving average as

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
**span** = 20 corresponds to **com** = 9.5.

Here are some examples for a univariate time series:

.. plot:: plots/stats/moments_expw.py

The binary `emwcov` and `ewmcorr` are similar to their equal-weighted
counterparts above.

.. autofunction:: pandas.stats.moments.ewma


.. _stats.ols:

Linear and panel regression
---------------------------

.. automodule:: pandas.stats.interface
   :members:
