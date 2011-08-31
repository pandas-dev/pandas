.. currentmodule:: pandas
.. _series:

Iterating
---------

Series iterates by default over its values as though it were a regular
ndarray.

Otherwise, methods providing dict-like iteration are available:

::

    >>> for x in ts:
            print x
    0.0
    1.0
    2.0
    3.0
    4.0

    >>> for index, value in ts.iteritems():
            print index, value
    2009-01-01 00:00:00 0.0
    2009-01-02 00:00:00 1.0
    2009-01-05 00:00:00 2.0
    2009-01-06 00:00:00 3.0
    2009-01-07 00:00:00 4.0


.. autosummary::
   :toctree: generated/

   Series.values
   Series.iteritems

.. _series.statistics:

Basic statistical functions
---------------------------

There are many built-in ndarray methods providing basic descriptive
statistics. Since these do not handle missing observations (which are
represented in our case as NaN), we've overridden these methods to do
the appropriate handling.

For example:

::

    >>> s
    a    0.0
    b    1.0
    c    nan
    d    3.0
    e    4.0

    >>> s.count()
    4

    >>> s.std()
    1.8257418583505536

    >>> s.cumsum()
    a    0.0
    b    1.0
    c    nan
    d    4.0
    e    8.0

Due to the way the numpy.{sum, mean, var, std} are implemented, they
can be used safely:

::

    >>> np.mean(s)
    2.0


TimeSeries-oriented methods
---------------------------

.. seealso::
    :ref:`Reindexing methods <series.reindexing>`;
    :ref:`DateRange and date offsets / time rules <datetools>`

.. note::

    While pandas does not force you to sort your dates, many of these
    methods may have unexpected or incorrect behavior in that case. In
    other words, *be careful*.

When working with time series data, a number of different
time-oriented operations may be useful. The first is **frequency
conversion**, which has similar options to :func:`Series.reindex`:

::

    >>> dr = DateRange('1/1/2010', periods=10,
                       offset=datetools.BMonthEnd())
    >>> ts = Series(np.arange(10.), index=dr)
    >>> ts
    2010-01-29 00:00:00    0.0
    2010-02-26 00:00:00    1.0
    2010-03-31 00:00:00    2.0
    2010-04-30 00:00:00    3.0
    2010-05-31 00:00:00    4.0
    2010-06-30 00:00:00    5.0
    2010-07-30 00:00:00    6.0
    2010-08-31 00:00:00    7.0
    2010-09-30 00:00:00    8.0
    2010-10-29 00:00:00    9.0

    >>> ts.asfreq('WEEKDAY', method='pad')
    2010-01-29 00:00:00    0.0
    2010-02-01 00:00:00    0.0
    2010-02-02 00:00:00    0.0
    2010-02-03 00:00:00    0.0
    2010-02-04 00:00:00    0.0
    <snip>
    2010-10-22 00:00:00    8.0
    2010-10-25 00:00:00    8.0
    2010-10-26 00:00:00    8.0
    2010-10-27 00:00:00    8.0
    2010-10-28 00:00:00    8.0
    2010-10-29 00:00:00    9.0

We often will also want to **shift** or *lag* a TimeSeries:

::

    >>> ts.shift(1)
    2010-01-29 00:00:00    NaN
    2010-02-26 00:00:00    0.0
    2010-03-31 00:00:00    1.0
    2010-04-30 00:00:00    2.0
    2010-05-31 00:00:00    3.0
    2010-06-30 00:00:00    4.0
    2010-07-30 00:00:00    5.0
    2010-08-31 00:00:00    6.0
    2010-09-30 00:00:00    7.0
    2010-10-29 00:00:00    8.0

    >>> ts.shift(5, offset=datetools.bday)
    2010-02-05 00:00:00    0.0
    2010-03-05 00:00:00    1.0
    2010-04-07 00:00:00    2.0
    2010-05-07 00:00:00    3.0
    2010-06-07 00:00:00    4.0
    2010-07-07 00:00:00    5.0
    2010-08-06 00:00:00    6.0
    2010-09-07 00:00:00    7.0
    2010-10-07 00:00:00    8.0
    2010-11-05 00:00:00    9.0

In the presence of missing data with sorted dates

A convenience method for selecting weekdays, similar to
:mod:`scikits.timeseries` is also provided:

::

    >>> dr = DateRange('1/1/2010', periods=10, offset=datetools.bday)
    >>> ts = Series(np.arange(10.), index=dr)
    >>> ts[ts.weekday == 2]
    2010-01-06 00:00:00    3.0
    2010-01-13 00:00:00    8.0
