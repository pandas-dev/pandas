.. currentmodule:: pandas
.. _timeseries:

.. ipython:: python
   :suppress:

   from datetime import datetime
   import numpy as np
   np.random.seed(123456)
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   from dateutil import relativedelta
   from pandas.tseries.api import *

********************************
Time Series / Date functionality
********************************

pandas has proven very successful as a tool for working with time series data,
especially in the financial data analysis space. With the 0.8 release, we have
further improved the time series API in pandas by leaps and bounds. Using the
new NumPy ``datetime64`` dtype, we have consolidated a large number of features
from other Python libraries like ``scikits.timeseries`` as well as created
a tremendous amount of new functionality for manipulating time series data.

In working with time series data, we will frequently seek to:

  - generate sequences of fixed-frequency dates and time spans
  - conform or convert time series to a particular frequency
  - compute "relative" dates based on various non-standard time increments
    (e.g. 5 business days before the last business day of the year), or "roll"
    dates forward or backward

pandas provides a relatively compact and self-contained set of tools for
performing the above tasks.

.. _timeseries.representation:

Time Stamps vs. Time Spans
--------------------------

Time-stamped data is the most basic type of timeseries data that associates
values with points in time. For pandas objects it means using the points in
time to create the index

.. ipython:: python

   dates = [datetime(2012, 5, 1), datetime(2012, 5, 2), datetime(2012, 5, 3)]
   ts = Series(np.random.randn(3), dates)

   type(ts.index)

   ts

However, in many cases it is more natural to associate things like change
variables with a time span instead.

For example:

.. ipython:: python

   periods = PeriodIndex([Period('2012-01'), Period('2012-02'),
                          Period('2012-03')])

   ts = Series(np.random.randn(3), periods)

   type(ts.index)

   ts

Starting with 0.8, pandas allows you to capture both representations and
convert between them. Under the hood, pandas represents timestamps using
instances of ``Timestamp`` and sequences of timestamps using instances of
``DatetimeIndex``. For regular time spans, pandas uses ``Period`` objects for
scalar values and ``PeriodIndex`` for sequences of spans. Better support for
irregular intervals with arbitrary start and end points are forth-coming in
future releases.

.. _timeseries.daterange:

Generating Ranges of Timestamps
-------------------------------

To generate an index with time stamps, you can use either the DatetimeIndex or
Index constructor and pass in a list of datetime objects:

.. ipython:: python

   dates = [datetime(2012, 5, 1), datetime(2012, 5, 2), datetime(2012, 5, 3)]
   index = DatetimeIndex(dates)
   index # Note the frequency information

   index = Index(dates)
   index # Automatically converted to DatetimeIndex

Practically, this becomes very cumbersome because we often need a very long
index with a large number of timestamps. If we need timestamps on a regular
frequency, we can use the pandas functions ``date_range`` and ``bdate_range``
to create timestamp indexes.

.. ipython:: python

   index = date_range('2000-1-1', periods=1000, freq='M')
   index

   index = bdate_range('2012-1-1', periods=250)
   index

Convenience functions like ``date_range`` and ``bdate_range`` utilizes a
variety of frequency aliases. The default frequency for ``date_range`` is a
**calendar day** while the default for ``bdate_range`` is a **business day**

.. ipython:: python

   start = datetime(2011, 1, 1)
   end = datetime(2012, 1, 1)

   rng = date_range(start, end)
   rng

   rng = bdate_range(start, end)
   rng

``date_range`` and ``bdate_range`` makes it easy to generate a range of dates
using various combinations of its parameters like ``start``, ``end``,
``periods``, and ``freq``:

   date_range(start, end, freq='BM')

   date_range(start, end, freq='W')

   bdate_range(end=end, periods=20)

   bdate_range(start=start, periods=20)

The start and end dates are strictly inclusive. So it will not generate any
dates outside of those dates if specified.

.. _timeseries.datetimeindex:

DatetimeIndex
~~~~~~~~~~~~~

One of the main uses for ``DatetimeIndex`` is as an index for pandas objects.
The ``DatetimeIndex`` class contains many timeseries related optimizations:

  - A large range of dates for various offsets are pre-computed and cached
    under the hood in order to make generating subsequent date ranges very fast
    (just have to grab a slice)
  - Fast shifting using the ``shift`` and ``tshift`` method on pandas objects
  - Unioning of overlapping DatetimeIndex objects with the same frequency is
    very fast (important for fast data alignment)
  - Quick access to date fields via properties such as ``year``, ``month``, etc.
  - Regularization functions like ``snap`` and very fast ``asof`` logic

``DatetimeIndex`` can be used like a regular index and offers all of its
intelligent functionality like selection, slicing, etc.

.. ipython:: python

   rng = date_range(start, end, freq='BM')
   ts = Series(randn(len(rng)), index=rng)
   ts.index
   ts[:5].index
   ts[::2].index

You can pass in dates and strings that parses to dates as indexing parameters:

.. ipython:: python

   ts['1/7/2011']

   ts[datetime(2011, 12, 25):]

   ts['10/31/2011':'11/30/2011']

A ``truncate`` convenience function is provided that is equivalent to slicing:

.. ipython:: python

   ts.truncate(after='10/31/2011', before='11/30/2011')

To provide convenience for accessing longer time series, you can also pass in
the year or year and month as strings:

.. ipython:: python

   ts['2011']

   ts['2011-6']

Even complicated fancy indexing that breaks the DatetimeIndex's frequency
regularity will result in a ``DatetimeIndex`` (but frequency is lost):

.. ipython:: python

   ts[[0, 2, 6]].index

DatetimeIndex objects has all the basic functionality of regular Index objects
and a smorgasbord of advanced timeseries-specific methods for easy frequency
processing.

.. seealso::
    :ref:`Reindexing methods <basics.reindexing>`

.. note::

    While pandas does not force you to have a sorted date index, some of these
    methods may have unexpected or incorrect behavior if the dates are
    unsorted. So please be careful.


.. _timeseries.offsets:

DateOffset objects
------------------

In the preceding examples, we created DatetimeIndex objects at various
frequencies by passing in frequency strings like 'M', 'W', and 'BM to the
``freq`` keyword. Under the hood, these frequency strings are being translated
into an instance of pandas ``DateOffset``, which represents a regular
frequency increment. Specific offset logic like "month", "business day", or
"one hour" is represented in its various subclasses.

.. csv-table::
    :header: "Class name", "Description"
    :widths: 15, 65

    DateOffset, "Generic offset class, defaults to 1 calendar day"
    BDay, "business day (weekday)"
    Week, "one week, optionally anchored on a day of the week"
    WeekOfMonth, "the x-th day of the y-th week of each month"
    MonthEnd, "calendar month end"
    MonthBegin, "calendar month begin"
    BMonthEnd, "business month end"
    BMonthBegin, "business month begin"
    QuarterEnd, "calendar quarter end"
    QuarterBegin, "calendar quarter begin"
    BQuarterEnd, "business quarter end"
    BQuarterBegin, "business quarter begin"
    YearEnd, "calendar year end"
    YearBegin, "calendar year begin"
    BYearEnd, "business year end"
    BYearBegin, "business year begin"
    Hour, "one hour"
    Minute, "one minute"
    Second, "one second"
    Milli, "one millisecond"
    Micro, "one microsecond"

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

It's definitely worth exploring the ``pandas.tseries.offsets`` module and the
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

Another example is parameterizing ``YearEnd`` with the specific ending month:

.. ipython:: python

   d + YearEnd()
   d + YearEnd(month=6)

.. _timeseries.alias:

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

These can be used as arguments to ``date_range``, ``bdate_range``, constructors
for ``DatetimeIndex``, as well as various other timeseries-related functions
in pandas.

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

As you can see, legacy quarterly and annual frequencies are business quarter
and business year ends. Please also note the legacy time rule for milliseconds
``ms`` versus the new offset alias for month start ``MS``. This means that
offset alias parsing is case sensitive.

.. _timeseries.advanced_datetime:

Time series-related instance methods
------------------------------------

Shifting / lagging
~~~~~~~~~~~~~~~~~~

One may want to *shift* or *lag* the values in a TimeSeries back and forward in
time. The method for this is ``shift``, which is available on all of the pandas
objects. In DataFrame, ``shift`` will currently only shift along the ``index``
and in Panel along the ``major_axis``.

.. ipython:: python

   ts = ts[:5]
   ts.shift(1)

The shift method accepts an ``freq`` argument which can accept a
``DateOffset`` class or other ``timedelta``-like object or also a :ref:`offset alias <timeseries.alias>`:

.. ipython:: python

   ts.shift(5, freq=datetools.bday)
   ts.shift(5, freq='BM')

Rather than changing the alignment of the data and the index, ``DataFrame`` and
``TimeSeries`` objects also have a ``tshift`` convenience method that changes
all the dates in the index by a specified number of offsets:

.. ipython:: python

   ts.tshift(5, freq='D')

Note that with ``tshift``, the leading entry is no longer NaN because the data
is not being realigned.

Frequency conversion
~~~~~~~~~~~~~~~~~~~~

The primary function for changing frequencies is the ``asfreq`` function.
For a ``DatetimeIndex``, this is basically just a thin, but convenient wrapper
around ``reindex`` which generates a ``date_range`` and calls ``reindex``.

.. ipython:: python

   dr = date_range('1/1/2010', periods=3, freq=3 * datetools.bday)
   ts = Series(randn(3), index=dr)
   ts
   ts.asfreq(BDay())

``asfreq`` provides a further convenience so you can specify an interpolation
method for any gaps that may appear after the frequency conversion

.. ipython:: python

   ts.asfreq(BDay(), method='pad')

Filling forward / backward
~~~~~~~~~~~~~~~~~~~~~~~~~~

Related to ``asfreq`` and ``reindex`` is the ``fillna`` function documented in
the :ref:`missing data section <missing_data.fillna>`.

.. _timeseries.resampling:

Up- and downsampling
--------------------

With 0.8, pandas introduces some efficient and easy-to-use methods for
performing resampling operations during frequency conversion (e.g., converting
secondly data into 5-minutely data). This is extremely common in, but not
limited to, financial applications.

.. ipython:: python



Until then, your best bet is a clever (or kludgy, depending on your point of
view) application of GroupBy. Carry out the following steps:

1. Generate the target ``date_range`` of interest

.. code-block:: python

   dr1hour = date_range(start, end, freq=Hour())
   dr5day = date_range(start, end, freq=5 * datetools.day)
   dr10day = date_range(start, end, freq=10 * datetools.day)


2. Use the ``asof`` function ("as of") of the date_range to do a groupby
   expression

.. code-block:: python

   grouped = data.groupby(dr5day.asof)
   means = grouped.mean()

Here is a fully-worked example:

.. ipython:: python

   # some minutely data
   minutely = date_range('1/3/2000 00:00:00', '1/3/2000 12:00:00',
                        freq=datetools.Minute())
   ts = Series(randn(len(minutely)), index=minutely)
   ts.index

   hourly = date_range('1/3/2000', '1/4/2000', freq=datetools.Hour())

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
