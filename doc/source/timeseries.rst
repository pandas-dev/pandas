.. currentmodule:: pandas
.. _timeseries:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   from dateutil import relativedelta
   from pandas.core.datetools import *

********************************
Time Series / Date functionality
********************************

pandas has proven very successful as a tool for working with time series data,
especially in the financial data analysis space. Over the coming year we will
be looking to consolidate the various Python libraries for time series data,
e.g. ``scikits.timeseries``, using the new NumPy ``datetime64`` dtype, to
create a very nice integrated solution. Everything in pandas at the moment is
based on using Python ``datetime`` objects.

In working with time series data, we will frequently seek to:

  - generate sequences of fixed-frequency dates
  - conform or convert time series to a particular frequency
  - compute "relative" dates based on various non-standard time increments
    (e.g. 5 business days before the last business day of the year), or "roll"
    dates forward or backward

pandas provides a relatively compact and self-contained set of tools for
performing the above tasks.

.. note::

   This area of pandas has gotten less development attention recently, though
   this should change in the near future.

.. _timeseries.offsets:

DateOffset objects
------------------

A ``DateOffset`` instance represents a frequency increment. Different offset
logic via subclasses:

.. csv-table::
    :header: "Class name", "Description"
    :widths: 15, 65

    DateOffset, "Generic offset class, defaults to 1 calendar day"
    BDay, "business day (weekday)"
    Week, "one week, optionally anchored on a day of the week"
    MonthEnd, "calendar month end"
    BMonthEnd, "business month end"
    BQuarterEnd, "business quarter end"
    YearEnd, "calendar year end"
    YearBegin, "calendar year begin"
    BYearEnd, "business year end"
    Hour, "one hour"
    Minute, "one minute"
    Second, "one second"

The basic ``DateOffset`` takes the same arguments as
``dateutil.relativedelta``, which works like:

.. ipython:: python

   d = datetime(2008, 8, 18)
   d + relativedelta(months=4, days=5)

We could have done the same thing with ``DateOffset``:

.. ipython:: python

   d + DateOffset(months=4, days=5)

The key features of a ``DateOffset`` object are:

  - it can be added / subtracted to/from a datetime object to obtain a
    shifted date
  - it can be multiplied by an integer (positive or negative) so that the
    increment will be applied multiple times
  - it has ``rollforward`` and ``rollback`` methods for moving a date forward
    or backward to the next or previous "offset date"

Subclasses of ``DateOffset`` define the ``apply`` function which dictates
custom date increment logic, such as adding business days:

.. code-block:: python

    class BDay(DateOffset):
	"""DateOffset increments between business days"""
        def apply(self, other):
            ...

.. ipython:: python

   d - 5 * BDay()
   d + BMonthEnd()

The ``rollforward`` and ``rollback`` methods do exactly what you would expect:

.. ipython:: python

   d
   offset = BMonthEnd()
   offset.rollforward(d)
   offset.rollback(d)

It's definitely worth exploring the ``pandas.core.datetools`` module and the
various docstrings for the classes.

Parametric offsets
~~~~~~~~~~~~~~~~~~

Some of the offsets can be "parameterized" when created to result in different
behavior. For example, the ``Week`` offset for generating weekly data accepts a
``weekday`` parameter which results in the generated dates always lying on a
particular day of the week:

.. ipython:: python

   d + Week()
   d + Week(weekday=4)
   (d + Week(weekday=4)).weekday()

Time rules
~~~~~~~~~~

A number of string aliases are given to useful common time series
frequencies. We will refer to these aliases as *time rules*.

.. csv-table::
    :header: "Rule name", "Description"
    :widths: 15, 65

    "WEEKDAY", "business day frequency"
    "EOM", "business month end frequency"
    "W\@MON", "weekly frequency (mondays)"
    "W\@TUE", "weekly frequency (tuesdays)"
    "W\@WED", "weekly frequency (wednesdays)"
    "W\@THU", "weekly frequency (thursdays)"
    "W\@FRI", "weekly frequency (fridays)"
    "Q\@JAN", "quarterly frequency, starting January"
    "Q\@FEB", "quarterly frequency, starting February"
    "Q\@MAR", "quarterly frequency, starting March"
    "A\@DEC", "annual frequency, year end (December)"
    "A\@JAN", "annual frequency, anchored end of January"
    "A\@FEB", "annual frequency, anchored end of February"
    "A\@MAR", "annual frequency, anchored end of March"
    "A\@APR", "annual frequency, anchored end of April"
    "A\@MAY", "annual frequency, anchored end of May"
    "A\@JUN", "annual frequency, anchored end of June"
    "A\@JUL", "annual frequency, anchored end of July"
    "A\@AUG", "annual frequency, anchored end of August"
    "A\@SEP", "annual frequency, anchored end of September"
    "A\@OCT", "annual frequency, anchored end of October"
    "A\@NOV", "annual frequency, anchored end of November"

These can be used as arguments to ``DateRange`` and various other time
series-related functions in pandas.

.. _timeseries.daterange:

Generating date ranges (DateRange)
----------------------------------

The ``DateRange`` class utilizes these offsets (and any ones that we might add)
to generate fixed-frequency date ranges:

.. ipython:: python

   start = datetime(2009, 1, 1)
   end = datetime(2010, 1, 1)

   rng = DateRange(start, end, offset=BDay())
   rng
   DateRange(start, end, offset=BMonthEnd())

**Business day frequency** is the default for ``DateRange``. You can also
strictly generate a ``DateRange`` of a certain length by providing either a
start or end date and a ``periods`` argument:

.. ipython:: python

   DateRange(start, periods=20)
   DateRange(end=end, periods=20)

The start and end dates are strictly inclusive. So it will not generate any
dates outside of those dates if specified.

DateRange is a valid Index
~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the main uses for ``DateRange`` is as an index for pandas objects. When
working with a lot of time series data, there are several reasons to use
``DateRange`` objects when possible:

  - A large range of dates for various offsets are pre-computed and cached
    under the hood in order to make generating subsequent date ranges very fast
    (just have to grab a slice)
  - Fast shifting using the ``shift`` method on pandas objects
  - Unioning of overlapping DateRange objects with the same frequency is very
    fast (important for fast data alignment)

The ``DateRange`` is a valid index and can even be intelligent when doing
slicing, etc.

.. ipython:: python

   rng = DateRange(start, end, offset=BMonthEnd())
   ts = Series(randn(len(rng)), index=rng)
   ts.index
   ts[:5].index
   ts[::2].index

More complicated fancy indexing will result in an ``Index`` that is no longer a
``DateRange``, however:

.. ipython:: python

   ts[[0, 2, 6]].index

Time series-related instance methods
------------------------------------

.. seealso::
    :ref:`Reindexing methods <basics.reindexing>`

.. note::

    While pandas does not force you to have a sorted date index, some of these
    methods may have unexpected or incorrect behavior if the dates are
    unsorted. So please be careful.

Shifting / lagging
~~~~~~~~~~~~~~~~~~

One may want to *shift* or *lag* the values in a TimeSeries back and forward in
time. The method for this is ``shift``, which is available on all of the pandas
objects. In DataFrame, ``shift`` will currently only shift along the ``index``
and in Panel along the ``major_axis``.

.. ipython:: python

   ts = ts[:5]
   ts.shift(1)
   ts.shift(5, offset=datetools.bday)

Frequency conversion
~~~~~~~~~~~~~~~~~~~~

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

Filling forward / backward
~~~~~~~~~~~~~~~~~~~~~~~~~~

TimeSeries-oriented methods
---------------------------
