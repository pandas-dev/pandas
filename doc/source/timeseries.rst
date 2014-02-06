.. currentmodule:: pandas
.. _timeseries:

.. ipython:: python
   :suppress:

   from datetime import datetime
   import numpy as np
   np.random.seed(123456)
   from pandas import *
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)
   options.display.max_rows=15
   from dateutil.relativedelta import relativedelta
   from pandas.tseries.api import *
   from pandas.tseries.offsets import *

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

Create a range of dates:

.. ipython:: python

   # 72 hours starting with midnight Jan 1st, 2011
   rng = date_range('1/1/2011', periods=72, freq='H')
   rng[:5]

Index pandas objects with dates:

.. ipython:: python

   ts = Series(randn(len(rng)), index=rng)
   ts.head()

Change frequency and fill gaps:

.. ipython:: python

   # to 45 minute frequency and forward fill
   converted = ts.asfreq('45Min', method='pad')
   converted.head()

Resample:

.. ipython:: python

   # Daily means
   ts.resample('D', how='mean')


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


.. _timeseries.converting:

Converting to Timestamps
------------------------

To convert a Series or list-like object of date-like objects e.g. strings,
epochs, or a mixture, you can use the ``to_datetime`` function. When passed
a Series, this returns a Series (with the same index), while a list-like
is converted to a DatetimeIndex:

.. ipython:: python

    to_datetime(Series(['Jul 31, 2009', '2010-01-10', None]))

    to_datetime(['2005/11/23', '2010.12.31'])

If you use dates which start with the day first (i.e. European style),
you can pass the ``dayfirst`` flag:

.. ipython:: python

    to_datetime(['04-01-2012 10:00'], dayfirst=True)

    to_datetime(['14-01-2012', '01-14-2012'], dayfirst=True)

.. warning::

   You see in the above example that ``dayfirst`` isn't strict, so if a date
   can't be parsed with the day being first it will be parsed as if
   ``dayfirst`` were False.

.. note::
   Specifying a ``format`` argument will potentially speed up the conversion
   considerably and on versions later then 0.13.0 explicitly specifying
   a format string of '%Y%m%d' takes a faster path still.


Invalid Data
~~~~~~~~~~~~

Pass ``coerce=True`` to convert invalid data to ``NaT`` (not a time):

.. ipython:: python

   to_datetime(['2009-07-31', 'asd'])

   to_datetime(['2009-07-31', 'asd'], coerce=True)


Take care, ``to_datetime`` may not act as you expect on mixed data:

.. ipython:: python

   to_datetime([1, '1'])

Epoch Timestamps
~~~~~~~~~~~~~~~~

It's also possible to convert integer or float epoch times. The default unit
for these is nanoseconds (since these are how Timestamps are stored). However,
often epochs are stored in another ``unit`` which can be specified:

Typical epoch stored units

.. ipython:: python

   to_datetime([1349720105, 1349806505, 1349892905,
                1349979305, 1350065705], unit='s')

   to_datetime([1349720105100, 1349720105200, 1349720105300,
                1349720105400, 1349720105500 ], unit='ms')

These *work*, but the results may be unexpected.

.. ipython:: python

   to_datetime([1])

   to_datetime([1, 3.14], unit='s')

.. note::

   Epoch times will be rounded to the nearest nanosecond.

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

Convenience functions like ``date_range`` and ``bdate_range`` utilize a
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
using various combinations of parameters like ``start``, ``end``,
``periods``, and ``freq``:

.. ipython:: python

   date_range(start, end, freq='BM')

   date_range(start, end, freq='W')

   bdate_range(end=end, periods=20)

   bdate_range(start=start, periods=20)

The start and end dates are strictly inclusive. So it will not generate any
dates outside of those dates if specified.

.. _timeseries.datetimeindex:

DatetimeIndex
-------------

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

DatetimeIndex objects has all the basic functionality of regular Index objects
and a smorgasbord of advanced timeseries-specific methods for easy frequency
processing.

.. seealso::
    :ref:`Reindexing methods <basics.reindexing>`

.. note::

    While pandas does not force you to have a sorted date index, some of these
    methods may have unexpected or incorrect behavior if the dates are
    unsorted. So please be careful.

``DatetimeIndex`` can be used like a regular index and offers all of its
intelligent functionality like selection, slicing, etc.

.. ipython:: python

   rng = date_range(start, end, freq='BM')
   ts = Series(randn(len(rng)), index=rng)
   ts.index
   ts[:5].index
   ts[::2].index

.. _timeseries.partialindexing:

DatetimeIndex Partial String Indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can pass in dates and strings that parse to dates as indexing parameters:

.. ipython:: python

   ts['1/31/2011']

   ts[datetime(2011, 12, 25):]

   ts['10/31/2011':'12/31/2011']

To provide convenience for accessing longer time series, you can also pass in
the year or year and month as strings:

.. ipython:: python

   ts['2011']

   ts['2011-6']

This type of slicing will work on a DataFrame with a ``DateTimeIndex`` as well. Since the
partial string selection is a form of label slicing, the endpoints **will be** included. This
would include matching times on an included date. Here's an example:

.. ipython:: python

   dft = DataFrame(randn(100000,1),columns=['A'],index=date_range('20130101',periods=100000,freq='T'))
   dft
   dft['2013']

This starts on the very first time in the month, and includes the last date & time for the month

.. ipython:: python

   dft['2013-1':'2013-2']

This specifies a stop time **that includes all of the times on the last day**

.. ipython:: python

   dft['2013-1':'2013-2-28']

This specifies an **exact** stop time (and is not the same as the above)

.. ipython:: python

   dft['2013-1':'2013-2-28 00:00:00']

We are stopping on the included end-point as its part of the index

.. ipython:: python

   dft['2013-1-15':'2013-1-15 12:30:00']

.. warning::

   The following selection will raise a ``KeyError``; otherwise this selection methodology
   would be inconsistent with other selection methods in pandas (as this is not a *slice*, nor does it
   resolve to one)

   .. code-block:: python

      dft['2013-1-15 12:30:00']

   To select a single row, use ``.loc``

   .. ipython:: python

      dft.loc['2013-1-15 12:30:00']


Datetime Indexing
~~~~~~~~~~~~~~~~~

Indexing a ``DateTimeIndex`` with a partial string depends on the "accuracy" of the period, in other words how specific the interval is in relation to the frequency of the index. In contrast, indexing with datetime objects is exact, because the objects have exact meaning. These also follow the sematics of *including both endpoints*.

These ``datetime`` objects  are specific ``hours, minutes,`` and ``seconds`` even though they were not explicity specified (they are ``0``).

.. ipython:: python

   dft[datetime(2013, 1, 1):datetime(2013,2,28)]

With no defaults.

.. ipython:: python

   dft[datetime(2013, 1, 1, 10, 12, 0):datetime(2013, 2, 28, 10, 12, 0)]


Truncating & Fancy Indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A ``truncate`` convenience function is provided that is equivalent to slicing:

.. ipython:: python

   ts.truncate(before='10/31/2011', after='12/31/2011')

Even complicated fancy indexing that breaks the DatetimeIndex's frequency
regularity will result in a ``DatetimeIndex`` (but frequency is lost):

.. ipython:: python

   ts[[0, 2, 6]].index

.. _timeseries.offsets:

Time/Date Components
~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are several time/date properties that one can access from ``Timestamp`` or a collection of timestamps like a ``DateTimeIndex``.

.. csv-table::
    :header: "Property", "Description"
    :widths: 15, 65

    year, "The year of the datetime"
    month,"The month of the datetime"
    day,"The days of the datetime"
    hour,"The hour of the datetime"
    minute,"The minutes of the datetime"
    second,"The seconds of the datetime"
    microsecond,"The microseconds of the datetime"
    nanosecond,"The nanoseconds of the datetime"
    date,"Returns datetime.date"
    time,"Returns datetime.time"
    dayofyear,"The ordinal day of year"
    weekofyear,"The week ordinal of the year"
    week,"The week ordinal of the year"
    dayofweek,"The day of the week with Monday=0, Sunday=6"
    weekday,"The day of the week with Monday=0, Sunday=6"
    quarter,"Quarter of the date: Jan=Mar = 1, Apr-Jun = 2, etc."
    is_month_start,"Logical indicating if first day of month (defined by frequency)"
    is_month_end,"Logical indicating if last day of month (defined by frequency)"
    is_quarter_start,"Logical indicating if first day of quarter (defined by frequency)"
    is_quarter_end,"Logical indicating if last day of quarter (defined by frequency)"
    is_year_start,"Logical indicating if first day of year (defined by frequency)"
    is_year_end,"Logical indicating if last day of year (defined by frequency)"


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
    CDay, "custom business day (experimental)"
    Week, "one week, optionally anchored on a day of the week"
    WeekOfMonth, "the x-th day of the y-th week of each month"
    LastWeekOfMonth, "the x-th day of the last week of each month"
    MonthEnd, "calendar month end"
    MonthBegin, "calendar month begin"
    BMonthEnd, "business month end"
    BMonthBegin, "business month begin"
    CBMonthEnd, "custom business month end"
    CBMonthBegin, "custom business month begin"
    QuarterEnd, "calendar quarter end"
    QuarterBegin, "calendar quarter begin"
    BQuarterEnd, "business quarter end"
    BQuarterBegin, "business quarter begin"
    FY5253Quarter, "retail (aka 52-53 week) quarter"
    YearEnd, "calendar year end"
    YearBegin, "calendar year begin"
    BYearEnd, "business year end"
    BYearBegin, "business year begin"
    FY5253, "retail (aka 52-53 week) year"
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

   from pandas.tseries.offsets import *
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

Custom Business Days (Experimental)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``CDay`` or ``CustomBusinessDay`` class provides a parametric
``BusinessDay`` class which can be used to create customized business day
calendars which account for local holidays and local weekend conventions.

.. ipython:: python

    from pandas.tseries.offsets import CustomBusinessDay
    # As an interesting example, let's look at Egypt where
    # a Friday-Saturday weekend is observed.
    weekmask_egypt = 'Sun Mon Tue Wed Thu'
    # They also observe International Workers' Day so let's
    # add that for a couple of years
    holidays = ['2012-05-01', datetime(2013, 5, 1), np.datetime64('2014-05-01')]
    bday_egypt = CustomBusinessDay(holidays=holidays, weekmask=weekmask_egypt)
    dt = datetime(2013, 4, 30)
    dt + 2 * bday_egypt
    dts = date_range(dt, periods=5, freq=bday_egypt)
    Series(dts.weekday, dts).map(Series('Mon Tue Wed Thu Fri Sat Sun'.split()))

As of v0.14 holiday calendars can be used to provide the list of holidays.  See the
:ref:`holiday calendar<timeseries.holiday>` section for more information.

.. ipython:: python

    from pandas.tseries.holiday import USFederalHolidayCalendar
    bday_us = CustomBusinessDay(calendar=USFederalHolidayCalendar())
    # Friday before MLK Day
    dt = datetime(2014, 1, 17)
    # Tuesday after MLK Day (Monday is skipped because it's a holiday)
    dt + bday_us

Monthly offsets that respect a certain holiday calendar can be defined
in the usual way.

.. ipython:: python

    from pandas.tseries.offsets import CustomBusinessMonthBegin
    bmth_us = CustomBusinessMonthBegin(calendar=USFederalHolidayCalendar())
    # Skip new years
    dt = datetime(2013, 12, 17)
    dt + bmth_us

    # Define date index with custom offset
    from pandas import DatetimeIndex
    DatetimeIndex(start='20100101',end='20120101',freq=bmth_us)

.. note::

    The frequency string 'C' is used to indicate that a CustomBusinessDay
    DateOffset is used, it is important to note that since CustomBusinessDay is
    a parameterised type, instances of CustomBusinessDay may differ and this is
    not detectable from the 'C' frequency string. The user therefore needs to
    ensure that the 'C' frequency string is used consistently within the user's
    application.


.. note::

    This uses the ``numpy.busdaycalendar`` API introduced in Numpy 1.7 and
    therefore requires Numpy 1.7.0 or newer.

.. warning::

    There are known problems with the timezone handling in Numpy 1.7 and users
    should therefore use this **experimental(!)** feature with caution and at
    their own risk.

    To the extent that the ``datetime64`` and ``busdaycalendar`` APIs in Numpy
    have to change to fix the timezone issues, the behaviour of the
    ``CustomBusinessDay`` class may have to change in future versions.

Offset Aliases
~~~~~~~~~~~~~~

A number of string aliases are given to useful common time series
frequencies. We will refer to these aliases as *offset aliases*
(referred to as *time rules* prior to v0.8.0).

.. csv-table::
    :header: "Alias", "Description"
    :widths: 15, 100

    "B", "business day frequency"
    "C", "custom business day frequency (experimental)"
    "D", "calendar day frequency"
    "W", "weekly frequency"
    "M", "month end frequency"
    "BM", "business month end frequency"
    "CBM", "custom business month end frequency"
    "MS", "month start frequency"
    "BMS", "business month start frequency"
    "CBMS", "custom business month start frequency"
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

Combining Aliases
~~~~~~~~~~~~~~~~~

As we have seen previously, the alias and the offset instance are fungible in
most functions:

.. ipython:: python

   date_range(start, periods=5, freq='B')

   date_range(start, periods=5, freq=BDay())

You can combine together day and intraday offsets:

.. ipython:: python

   date_range(start, periods=10, freq='2h20min')

   date_range(start, periods=10, freq='1D10U')

Anchored Offsets
~~~~~~~~~~~~~~~~

For some frequencies you can specify an anchoring suffix:

.. csv-table::
    :header: "Alias", "Description"
    :widths: 15, 100

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

Legacy Aliases
~~~~~~~~~~~~~~
Note that prior to v0.8.0, time rules had a slightly different look. pandas
will continue to support the legacy time rules for the time being but it is
strongly recommended that you switch to using the new offset aliases.

.. csv-table::
    :header: "Legacy Time Rule", "Offset Alias"
    :widths: 15, 65

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
    "us", "U"

As you can see, legacy quarterly and annual frequencies are business quarter
and business year ends. Please also note the legacy time rule for milliseconds
``ms`` versus the new offset alias for month start ``MS``. This means that
offset alias parsing is case sensitive.

.. _timeseries.holiday:

Holidays / Holiday Calendars
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Holidays and calendars provide a simple way to define holiday rules to be used
with ``CustomBusinessDay`` or in other analysis that requires a predefined
set of holidays.  The ``AbstractHolidayCalendar`` class provides all the necessary
methods to return a list of holidays and only ``rules`` need to be defined
in a specific holiday calendar class.  Further, ``start_date`` and ``end_date``
class attributes determine over what date range holidays are generated.  These
should be overwritten on the ``AbstractHolidayCalendar`` class to have the range
apply to all calendar subclasses.  ``USFederalHolidayCalendar`` is the
only calendar that exists and primarily serves as an example for developing
other calendars.

For holidays that occur on fixed dates (e.g., US Memorial Day or July 4th) an
observance rule determines when that holiday is observed if it falls on a weekend
or some other non-observed day.  Defined observance rules are:

.. csv-table::
    :header: "Rule", "Description"
    :widths: 15, 70

    "nearest_workday", "move Saturday to Friday and Sunday to Monday"
    "sunday_to_monday", "move Sunday to following Monday"
    "next_monday_or_tuesday", "move Saturday to Monday and Sunday/Monday to Tuesday"
    "previous_friday", move Saturday and Sunday to previous Friday"
    "next_monday", "move Saturday and Sunday to following Monday"

An example of how holidays and holiday calendars are defined:

.. ipython:: python

    from pandas.tseries.holiday import Holiday, USMemorialDay,\
        AbstractHolidayCalendar, nearest_workday, MO
    class ExampleCalendar(AbstractHolidayCalendar):
        rules = [
            USMemorialDay,
            Holiday('July 4th', month=7, day=4, observance=nearest_workday),
            Holiday('Columbus Day', month=10, day=1,
                offset=DateOffset(weekday=MO(2))), #same as 2*Week(weekday=2)
            ]
    cal = ExampleCalendar()
    cal.holidays(datetime(2012, 1, 1), datetime(2012, 12, 31))

Using this calendar, creating an index or doing offset arithmetic skips weekends
and holidays (i.e., Memorial Day/July 4th).

.. ipython:: python

    DatetimeIndex(start='7/1/2012', end='7/10/2012',
        freq=CDay(calendar=cal)).to_pydatetime()
    offset = CustomBusinessDay(calendar=cal)
    datetime(2012, 5, 25) + offset
    datetime(2012, 7, 3) + offset
    datetime(2012, 7, 3) + 2 * offset
    datetime(2012, 7, 6) + offset

Ranges are defined by the ``start_date`` and ``end_date`` class attributes
of ``AbstractHolidayCalendar``.  The defaults are below.

.. ipython:: python

    AbstractHolidayCalendar.start_date
    AbstractHolidayCalendar.end_date

These dates can be overwritten by setting the attributes as
datetime/Timestamp/string.

.. ipython:: python

    AbstractHolidayCalendar.start_date = datetime(2012, 1, 1)
    AbstractHolidayCalendar.end_date = datetime(2012, 12, 31)
    cal.holidays()

Every calendar class is accessible by name using the ``get_calendar`` function
which returns a holiday class instance.  Any imported calendar class will
automatically be available by this function.  Also, ``HolidayCalendarFactory``
provides an easy interface to create calendars that are combinations of calendars
or calendars with additional rules.

.. ipython:: python

    from pandas.tseries.holiday import get_calendar, HolidayCalendarFactory,\
        USLaborDay
    cal = get_calendar('ExampleCalendar')
    cal.rules
    new_cal = HolidayCalendarFactory('NewExampleCalendar', cal, USLaborDay)
    new_cal.rules

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

Converting to Python datetimes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``DatetimeIndex`` can be converted to an array of Python native datetime.datetime objects using the
``to_pydatetime`` method.

.. _timeseries.resampling:

Up- and downsampling
--------------------

With 0.8, pandas introduces simple, powerful, and efficient functionality for
performing resampling operations during frequency conversion (e.g., converting
secondly data into 5-minutely data). This is extremely common in, but not
limited to, financial applications.

See some :ref:`cookbook examples <cookbook.resample>` for some advanced strategies

.. ipython:: python

   rng = date_range('1/1/2012', periods=100, freq='S')

   ts = Series(randint(0, 500, len(rng)), index=rng)

   ts.resample('5Min', how='sum')

The ``resample`` function is very flexible and allows you to specify many
different parameters to control the frequency conversion and resampling
operation.

The ``how`` parameter can be a function name or numpy array function that takes
an array and produces aggregated values:

.. ipython:: python

   ts.resample('5Min') # default is mean

   ts.resample('5Min', how='ohlc')

   ts.resample('5Min', how=np.max)

Any function available via :ref:`dispatching <groupby.dispatch>` can be given to
the ``how`` parameter by name, including ``sum``, ``mean``, ``std``, ``max``,
``min``, ``median``, ``first``, ``last``, ``ohlc``.

For downsampling, ``closed`` can be set to 'left' or 'right' to specify which
end of the interval is closed:

.. ipython:: python

   ts.resample('5Min', closed='right')

   ts.resample('5Min', closed='left')

For upsampling, the ``fill_method`` and ``limit`` parameters can be specified
to interpolate over the gaps that are created:

.. ipython:: python

   # from secondly to every 250 milliseconds

   ts[:2].resample('250L')

   ts[:2].resample('250L', fill_method='pad')

   ts[:2].resample('250L', fill_method='pad', limit=2)

Parameters like ``label`` and ``loffset`` are used to manipulate the resulting
labels. ``label`` specifies whether the result is labeled with the beginning or
the end of the interval. ``loffset`` performs a time adjustment on the output
labels.

.. ipython:: python

   ts.resample('5Min') # by default label='right'

   ts.resample('5Min', label='left')

   ts.resample('5Min', label='left', loffset='1s')

The ``axis`` parameter can be set to 0 or 1 and allows you to resample the
specified axis for a DataFrame.

``kind`` can be set to 'timestamp' or 'period' to convert the resulting index
to/from time-stamp and time-span representations. By default ``resample``
retains the input representation.

``convention`` can be set to 'start' or 'end' when resampling period data
(detail below). It specifies how low frequency periods are converted to higher
frequency periods.

Note that 0.8 marks a watershed in the timeseries functionality in pandas. In
previous versions, resampling had to be done using a combination of
``date_range``, ``groupby`` with ``asof``, and then calling an aggregation
function on the grouped object. This was not nearly convenient or performant as
the new pandas timeseries API.

.. _timeseries.periods:

Time Span Representation
------------------------

Regular intervals of time are represented by ``Period`` objects in pandas while
sequences of ``Period`` objects are collected in a ``PeriodIndex``, which can
be created with the convenience function ``period_range``.

Period
~~~~~~
A ``Period`` represents a span of time (e.g., a day, a month, a quarter, etc).
It can be created using a frequency alias:

.. ipython:: python

   Period('2012', freq='A-DEC')

   Period('2012-1-1', freq='D')

   Period('2012-1-1 19:00', freq='H')

Unlike time stamped data, pandas does not support frequencies at multiples of
DateOffsets (e.g., '3Min') for periods.

Adding and subtracting integers from periods shifts the period by its own
frequency.

.. ipython:: python

   p = Period('2012', freq='A-DEC')

   p + 1

   p - 3

Taking the difference of ``Period`` instances with the same frequency will
return the number of frequency units between them:

.. ipython:: python

   Period('2012', freq='A-DEC') - Period('2002', freq='A-DEC')

PeriodIndex and period_range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Regular sequences of ``Period`` objects can be collected in a ``PeriodIndex``,
which can be constructed using the ``period_range`` convenience function:

.. ipython:: python

   prng = period_range('1/1/2011', '1/1/2012', freq='M')
   prng

The ``PeriodIndex`` constructor can also be used directly:

.. ipython:: python

   PeriodIndex(['2011-1', '2011-2', '2011-3'], freq='M')

Just like ``DatetimeIndex``, a ``PeriodIndex`` can also be used to index pandas
objects:

.. ipython:: python

   ps = Series(randn(len(prng)), prng)
   ps

PeriodIndex Partial String Indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can pass in dates and strings to `Series` and `DataFrame` with `PeriodIndex`, as the same manner as `DatetimeIndex`. For details, refer to :ref:`DatetimeIndex Partial String Indexing <timeseries.partialindexing>`.

.. ipython:: python

   ps['2011-01']

   ps[datetime(2011, 12, 25):]

   ps['10/31/2011':'12/31/2011']

Passing string represents lower frequency than `PeriodIndex` returns partial sliced data.

.. ipython:: python

   ps['2011']

   dfp = DataFrame(randn(600,1), columns=['A'],
                   index=period_range('2013-01-01 9:00', periods=600, freq='T'))
   dfp
   dfp['2013-01-01 10H']

As the same as `DatetimeIndex`, the endpoints will be included in the result. Below example slices data starting from 10:00 to 11:59.

.. ipython:: python

   dfp['2013-01-01 10H':'2013-01-01 11H']

Frequency Conversion and Resampling with PeriodIndex
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The frequency of Periods and PeriodIndex can be converted via the ``asfreq``
method. Let's start with the fiscal year 2011, ending in December:

.. ipython:: python

   p = Period('2011', freq='A-DEC')
   p

We can convert it to a monthly frequency. Using the ``how`` parameter, we can
specify whether to return the starting or ending month:

.. ipython:: python

   p.asfreq('M', how='start')

   p.asfreq('M', how='end')

The shorthands 's' and 'e' are provided for convenience:

.. ipython:: python

   p.asfreq('M', 's')
   p.asfreq('M', 'e')

Converting to a "super-period" (e.g., annual frequency is a super-period of
quarterly frequency) automatically returns the super-period that includes the
input period:

.. ipython:: python

   p = Period('2011-12', freq='M')

   p.asfreq('A-NOV')

Note that since we converted to an annual frequency that ends the year in
November, the monthly period of December 2011 is actually in the 2012 A-NOV
period.

.. _timeseries.quarterly:

Period conversions with anchored frequencies are particularly useful for
working with various quarterly data common to economics, business, and other
fields. Many organizations define quarters relative to the month in which their
fiscal year start and ends. Thus, first quarter of 2011 could start in 2010 or
a few months into 2011. Via anchored frequencies, pandas works all quarterly
frequencies ``Q-JAN`` through ``Q-DEC``.

``Q-DEC`` define regular calendar quarters:

.. ipython:: python

   p = Period('2012Q1', freq='Q-DEC')

   p.asfreq('D', 's')

   p.asfreq('D', 'e')

``Q-MAR`` defines fiscal year end in March:

.. ipython:: python

   p = Period('2011Q4', freq='Q-MAR')

   p.asfreq('D', 's')

   p.asfreq('D', 'e')

.. _timeseries.interchange:

Converting between Representations
----------------------------------

Timestamped data can be converted to PeriodIndex-ed data using ``to_period``
and vice-versa using ``to_timestamp``:

.. ipython:: python

   rng = date_range('1/1/2012', periods=5, freq='M')

   ts = Series(randn(len(rng)), index=rng)

   ts

   ps = ts.to_period()

   ps

   ps.to_timestamp()

Remember that 's' and 'e' can be used to return the timestamps at the start or
end of the period:

.. ipython:: python

   ps.to_timestamp('D', how='s')

Converting between period and timestamp enables some convenient arithmetic
functions to be used. In the following example, we convert a quarterly
frequency with year ending in November to 9am of the end of the month following
the quarter end:

.. ipython:: python

   prng = period_range('1990Q1', '2000Q4', freq='Q-NOV')

   ts = Series(randn(len(prng)), prng)

   ts.index = (prng.asfreq('M', 'e') + 1).asfreq('H', 's') + 9

   ts.head()

.. _timeseries.timezone:

Time Zone Handling
------------------

Pandas provides rich support for working with timestamps in different time zones using ``pytz`` and ``dateutil`` libraries. 
``dateutil`` support is new [in 0.14.1] and currently only supported for fixed offset and tzfile zones. The default library is ``pytz``.
Support for ``dateutil`` is provided for compatibility with other applications e.g. if you use ``dateutil`` in other python packages.

By default, pandas objects are time zone unaware:

.. ipython:: python

   rng = date_range('3/6/2012 00:00', periods=15, freq='D')
   rng.tz is None

To supply the time zone, you can use the ``tz`` keyword to ``date_range`` and
other functions. Dateutil time zone strings are distinguished from ``pytz``
time zones by starting with ``dateutil/``. In ``pytz`` you can find a list of
common (and less common) time zones using ``from pytz import common_timezones, all_timezones``. 
``dateutil`` uses the OS timezones so there isn't a fixed list available. For
common zones, the names are the same as ``pytz``.

.. ipython:: python
   
   # pytz
   rng_utc = date_range('3/6/2012 00:00', periods=10, freq='D', tz='UTC')
   rng_utc.tz

   # dateutil
   rng_utc_dateutil = date_range('3/6/2012 00:00', periods=10, freq='D',
                                 tz='dateutil/UTC')
   rng_utc_dateutil.tz

You can also construct the timezone explicitly first, which gives you more control over which
time zone is used:

.. ipython:: python

   # pytz
   import pytz
   tz_pytz = pytz.timezone('UTC')
   rng_utc = date_range('3/6/2012 00:00', periods=10, freq='D', tz=tz_pytz)
   rng_utc.tz

   # dateutil
   import dateutil
   tz_dateutil = dateutil.tz.gettz('UTC')
   rng_utc_dateutil = date_range('3/6/2012 00:00', periods=10, freq='D',
                                 tz=tz_dateutil)
   rng_utc_dateutil.tz


Timestamps, like Python's ``datetime.datetime`` object can be either time zone
naive or time zone aware. Naive time series and DatetimeIndex objects can be
*localized* using ``tz_localize``:

.. ipython:: python

   ts = Series(randn(len(rng)), rng)

   ts_utc = ts.tz_localize('UTC')
   ts_utc

Again, you can explicitly construct the timezone object first.
You can use the ``tz_convert`` method to convert pandas objects to convert
tz-aware data to another time zone:

.. ipython:: python

   ts_utc.tz_convert('US/Eastern')

.. warning::
    Be very wary of conversions between libraries as ``pytz`` and ``dateutil``
    may have different definitions of the time zones. This is more of a problem for
    unusual timezones than for 'standard' zones like ``US/Eastern``. 

Under the hood, all timestamps are stored in UTC. Scalar values from a
``DatetimeIndex`` with a time zone will have their fields (day, hour, minute)
localized to the time zone. However, timestamps with the same UTC value are
still considered to be equal even if they are in different time zones:

.. ipython:: python

   rng_eastern = rng_utc.tz_convert('US/Eastern')
   rng_berlin = rng_utc.tz_convert('Europe/Berlin')

   rng_eastern[5]
   rng_berlin[5]
   rng_eastern[5] == rng_berlin[5]

Like Series, DataFrame, and DatetimeIndex, Timestamps can be converted to other
time zones using ``tz_convert``:

.. ipython:: python

   rng_eastern[5]
   rng_berlin[5]
   rng_eastern[5].tz_convert('Europe/Berlin')

Localization of Timestamps functions just like DatetimeIndex and TimeSeries:

.. ipython:: python

   rng[5]
   rng[5].tz_localize('Asia/Shanghai')


Operations between TimeSeries in different time zones will yield UTC
TimeSeries, aligning the data on the UTC timestamps:

.. ipython:: python

   eastern = ts_utc.tz_convert('US/Eastern')
   berlin = ts_utc.tz_convert('Europe/Berlin')
   result = eastern + berlin
   result
   result.index

.. _timeseries.timedeltas:

In some cases, localize cannot determine the DST and non-DST hours when there are
duplicates.  This often happens when reading files that simply duplicate the hours.
The infer_dst argument in tz_localize will attempt
to determine the right offset.

.. ipython:: python
   :okexcept:

   rng_hourly = DatetimeIndex(['11/06/2011 00:00', '11/06/2011 01:00',
                               '11/06/2011 01:00', '11/06/2011 02:00',
                               '11/06/2011 03:00'])
   rng_hourly.tz_localize('US/Eastern')
   rng_hourly_eastern = rng_hourly.tz_localize('US/Eastern', infer_dst=True)
   rng_hourly_eastern.values

Time Deltas
-----------

Timedeltas are differences in times, expressed in difference units, e.g. days,hours,minutes,seconds.
They can be both positive and negative. :ref:`DateOffsets<timeseries.offsets>` that are absolute in nature
(``Day, Hour, Minute, Second, Milli, Micro, Nano``) can be used as ``timedeltas``.

.. ipython:: python

   from datetime import datetime, timedelta
   s = Series(date_range('2012-1-1', periods=3, freq='D'))
   td = Series([ timedelta(days=i) for i in range(3) ])
   df = DataFrame(dict(A = s, B = td))
   df
   df['C'] = df['A'] + df['B']
   df
   df.dtypes

   s - s.max()
   s - datetime(2011,1,1,3,5)
   s + timedelta(minutes=5)
   s + Minute(5)
   s + Minute(5) + Milli(5)

Getting scalar results from a ``timedelta64[ns]`` series

.. ipython:: python

   y = s - s[0]
   y

Series of timedeltas with ``NaT`` values are supported

.. ipython:: python

   y = s - s.shift()
   y

Elements can be set to ``NaT`` using ``np.nan`` analagously to datetimes

.. ipython:: python

   y[1] = np.nan
   y

Operands can also appear in a reversed order (a singluar object operated with a Series)

.. ipython:: python

   s.max() - s
   datetime(2011,1,1,3,5) - s
   timedelta(minutes=5) + s

Some timedelta numeric like operations are supported.

.. ipython:: python

   td - timedelta(minutes=5, seconds=5, microseconds=5)

``min, max`` and the corresponding ``idxmin, idxmax`` operations are supported on frames

.. ipython:: python

   A = s - Timestamp('20120101') - timedelta(minutes=5, seconds=5)
   B = s - Series(date_range('2012-1-2', periods=3, freq='D'))

   df = DataFrame(dict(A=A, B=B))
   df

   df.min()
   df.min(axis=1)

   df.idxmin()
   df.idxmax()

``min, max`` operations are supported on series; these return a single element
``timedelta64[ns]`` Series (this avoids having to deal with numpy timedelta64
issues). ``idxmin, idxmax`` are supported as well.

.. ipython:: python

   df.min().max()
   df.min(axis=1).min()

   df.min().idxmax()
   df.min(axis=1).idxmin()

You can fillna on timedeltas. Integers will be interpreted as seconds. You can
pass a timedelta to get a particular value.

.. ipython:: python

   y.fillna(0)
   y.fillna(10)
   y.fillna(timedelta(days=-1,seconds=5))

.. _timeseries.timedeltas_reductions:

Time Deltas & Reductions
------------------------

.. warning::

   A numeric reduction operation for ``timedelta64[ns]`` can return a single-element ``Series`` of
   dtype ``timedelta64[ns]``.

You can do numeric reduction operations on timedeltas.

.. ipython:: python

   y2 = y.fillna(timedelta(days=-1,seconds=5))
   y2
   y2.mean()
   y2.quantile(.1)

.. _timeseries.timedeltas_convert:

Time Deltas & Conversions
-------------------------

.. versionadded:: 0.13

**string/integer conversion**

Using the top-level ``to_timedelta``, you can convert a scalar or array from the standard
timedelta format (produced by ``to_csv``) into a timedelta type (``np.timedelta64`` in ``nanoseconds``).
It can also construct Series.

.. warning::

   This requires ``numpy >= 1.7``

.. ipython:: python

   to_timedelta('1 days 06:05:01.00003')
   to_timedelta('15.5us')
   to_timedelta(['1 days 06:05:01.00003','15.5us','nan'])
   to_timedelta(np.arange(5),unit='s')
   to_timedelta(np.arange(5),unit='d')

**frequency conversion**

Timedeltas can be converted to other 'frequencies' by dividing by another timedelta,
or by astyping to a specific timedelta type. These operations yield ``float64`` dtyped Series.

.. ipython:: python

   td = Series(date_range('20130101',periods=4))-Series(date_range('20121201',periods=4))
   td[2] += np.timedelta64(timedelta(minutes=5,seconds=3))
   td[3] = np.nan
   td

   # to days
   td / np.timedelta64(1,'D')
   td.astype('timedelta64[D]')

   # to seconds
   td / np.timedelta64(1,'s')
   td.astype('timedelta64[s]')

Dividing or multiplying a ``timedelta64[ns]`` Series by an integer or integer Series
yields another ``timedelta64[ns]`` dtypes Series.

.. ipython:: python

   td * -1
   td * Series([1,2,3,4])

Numpy < 1.7 Compatibility
~~~~~~~~~~~~~~~~~~~~~~~~~

Numpy < 1.7 has a broken ``timedelta64`` type that does not correctly work
for arithmetic. pandas bypasses this, but for frequency conversion as above,
you need to create the divisor yourself. The ``np.timetimedelta64`` type only
has 1 argument, the number of **micro** seconds.

The following are equivalent statements in the two versions of numpy.

.. code-block:: python

   from distutils.version import LooseVersion
   if LooseVersion(np.__version__) <= '1.6.2':
       y / np.timedelta(86400*int(1e6))
       y / np.timedelta(int(1e6))
   else:
       y / np.timedelta64(1,'D')
       y / np.timedelta64(1,'s')
