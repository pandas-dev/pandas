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
    QuarterEnd, "calendar quarter end"
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

   from pandas.core.datetools import *
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

.. _timeseries.freq:

Offset Aliases
~~~~~~~~~~~~~~

A number of string aliases are given to useful common time series
frequencies. We will refer to these aliases as *offset aliases*
(referred to as *time rules* prior to v0.8.0).

.. csv-table::
    :header: "Alias", "Description"
    :widths: 15, 65

    "B", "business day frequency"
    "D", "calendar day frequency"
    "W", "weekly frequency"
    "M", "month end frequency"
    "BM", "business month end frequency"
    "MS", "month start frequency"
    "BMS", "business month start frequency"
    "Q", "quarter end frequency"
    "BQ", "business quarter endfrequency"
    "QS", "quarter start frequency"
    "BQS", "business quarter start frequency"
    "A", "year end frequency"
    "BA", "business year end frequency"
    "AS", "year start frequency"
    "BAS", "business year start frequency"
    "H", "hourly frequency"
    "T", "minutely frequency"
    "S", "secondly frequency"
    "L", "milliseonds"
    "U", "microseconds"

For some frequencies you can specify an anchoring suffix:

.. csv-table::
    :header: "Alias", "Description"
    :widths: 15, 65

    "W\-SUN", "weekly frequency (sundays). Same as 'W'"
    "W\-MON", "weekly frequency (mondays)"
    "W\-TUE", "weekly frequency (tuesdays)"
    "W\-WED", "weekly frequency (wednesdays)"
    "W\-THU", "weekly frequency (thursdays)"
    "W\-FRI", "weekly frequency (fridays)"
    "W\-SAT", "weekly frequency (saturdays)"
    "(B)Q(S)\-DEC", "quarterly frequency, year ends in December. Same as 'Q'"
    "(B)Q(S)\-JAN", "quarterly frequency, year ends in January"
    "(B)Q(S)\-FEB", "quarterly frequency, year ends in February"
    "(B)Q(S)\-MAR", "quarterly frequency, year ends in March"
    "(B)Q(S)\-APR", "quarterly frequency, year ends in April"
    "(B)Q(S)\-MAY", "quarterly frequency, year ends in May"
    "(B)Q(S)\-JUN", "quarterly frequency, year ends in June"
    "(B)Q(S)\-JUL", "quarterly frequency, year ends in July"
    "(B)Q(S)\-AUG", "quarterly frequency, year ends in August"
    "(B)Q(S)\-SEP", "quarterly frequency, year ends in September"
    "(B)Q(S)\-OCT", "quarterly frequency, year ends in October"
    "(B)Q(S)\-NOV", "quarterly frequency, year ends in November"
    "(B)A(S)\-DEC", "annual frequency, anchored end of December. Same as 'A'"
    "(B)A(S)\-JAN", "annual frequency, anchored end of January"
    "(B)A(S)\-FEB", "annual frequency, anchored end of February"
    "(B)A(S)\-MAR", "annual frequency, anchored end of March"
    "(B)A(S)\-APR", "annual frequency, anchored end of April"
    "(B)A(S)\-MAY", "annual frequency, anchored end of May"
    "(B)A(S)\-JUN", "annual frequency, anchored end of June"
    "(B)A(S)\-JUL", "annual frequency, anchored end of July"
    "(B)A(S)\-AUG", "annual frequency, anchored end of August"
    "(B)A(S)\-SEP", "annual frequency, anchored end of September"
    "(B)A(S)\-OCT", "annual frequency, anchored end of October"
    "(B)A(S)\-NOV", "annual frequency, anchored end of November"

These can be used as arguments to ``date_range``, ``period_range``, constructors
for ``PeriodIndex`` and ``DatetimeIndex``, as well as various other time
series-related functions in pandas.

Note that prior to v0.8.0, time rules had a slightly different look. Pandas
will continue to support the legacy time rules for the time being but it is
strongly recommended that you switch to using the new offset aliases.

.. csv-table::
    :header: "Legacy Time Rule", "Offset Alias"
    :widths: 15, 15

    "WEEKDAY", "B"
    "EOM", "BM"
    "W\@MON", "W\-MON"
    "W\@TUE", "W\-TUE"
    "W\@WED", "W\-WED"
    "W\@THU", "W\-THU"
    "W\@FRI", "W\-FRI"
    "W\@SAT", "W\-SAT"
    "W\@SUN", "W\-SUN"
    "Q\@JAN", "BQ\-JAN"
    "Q\@FEB", "BQ\-FEB"
    "Q\@MAR", "BQ\-MAR"
    "A\@JAN", "BA\-JAN"
    "A\@FEB", "BA\-FEB"
    "A\@MAR", "BA\-MAR"
    "A\@APR", "BA\-APR"
    "A\@MAY", "BA\-MAY"
    "A\@JUN", "BA\-JUN"
    "A\@JUL", "BA\-JUL"
    "A\@AUG", "BA\-AUG"
    "A\@SEP", "BA\-SEP"
    "A\@OCT", "BA\-OCT"
    "A\@NOV", "BA\-NOV"
    "A\@DEC", "BA\-DEC"
    "min", "T"
    "ms", "L"
    "us": "U"

Note that the legacy quarterly and annual frequencies are business quarter and
business year ends. Also note the legacy time rule for milliseconds ``ms``
versus the new offset alias for month start ``MS``. This means that offset
alias parsing is case sensitive.

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

The shift method accepts an ``offset`` argument which can accept a
``DateOffset`` class or other ``timedelta``-like object or also a :ref:`time
rule <timeseries.timerule>`:

.. ipython:: python

   ts.shift(5, offset=datetools.bday)
   ts.shift(5, offset='EOM')

Frequency conversion
~~~~~~~~~~~~~~~~~~~~

The primary function for changing frequencies is the ``asfreq`` function. This
is basically just a thin, but convenient wrapper around ``reindex`` which
generates a ``DateRange`` and calls ``reindex``.

.. ipython:: python

   dr = DateRange('1/1/2010', periods=3, offset=3 * datetools.bday)
   ts = Series(randn(3), index=dr)
   ts
   ts.asfreq(BDay())
   ts.asfreq(BDay(), method='pad')

Filling forward / backward
~~~~~~~~~~~~~~~~~~~~~~~~~~

Related to ``asfreq`` and ``reindex`` is the ``fillna`` function documented in
the :ref:`missing data section <missing_data.fillna>`.

Up- and downsampling
--------------------

We plan to add some efficient methods for doing resampling during frequency
conversion. For example, converting secondly data into 5-minutely data. This is
extremely common in, but not limited to, financial applications.

Until then, your best bet is a clever (or kludgy, depending on your point of
view) application of GroupBy. Carry out the following steps:

1. Generate the target ``DateRange`` of interest

.. code-block:: python

   dr1hour = DateRange(start, end, offset=Hour())
   dr5day = DateRange(start, end, offset=5 * datetools.day)
   dr10day = DateRange(start, end, offset=10 * datetools.day)


2. Use the ``asof`` function ("as of") of the DateRange to do a groupby
   expression

.. code-block:: python

   grouped = data.groupby(dr5day.asof)
   means = grouped.mean()

Here is a fully-worked example:

.. ipython:: python

   # some minutely data
   minutely = DateRange('1/3/2000 00:00:00', '1/3/2000 12:00:00',
                        offset=datetools.Minute())
   ts = Series(randn(len(minutely)), index=minutely)
   ts.index

   hourly = DateRange('1/3/2000', '1/4/2000', offset=datetools.Hour())

   grouped = ts.groupby(hourly.asof)
   grouped.mean()

Some things to note about this method:

  - This is rather inefficient because we haven't exploited the orderedness of
    the data at all. Calling the ``asof`` function on every date in the
    minutely time series is not strictly necessary. We'll be writing some
    significantly more efficient methods in the near future
  - The dates in the result mark the **beginning of the period**. Be careful
    about which convention you use; you don't want to end up misaligning data
    because you used the wrong upsampling convention
