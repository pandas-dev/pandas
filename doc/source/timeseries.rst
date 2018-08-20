.. currentmodule:: pandas
.. _timeseries:

.. ipython:: python
   :suppress:

   from datetime import datetime, timedelta, time
   import numpy as np
   import pandas as pd
   from pandas import offsets
   np.random.seed(123456)
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)
   pd.options.display.max_rows=15
   import dateutil
   import pytz
   from dateutil.relativedelta import relativedelta

********************************
Time Series / Date functionality
********************************

pandas has proven very successful as a tool for working with time series data,
especially in the financial data analysis space. Using the NumPy ``datetime64`` and ``timedelta64`` dtypes,
we have consolidated a large number of features from other Python libraries like ``scikits.timeseries`` as well as created
a tremendous amount of new functionality for manipulating time series data.

In working with time series data, we will frequently seek to:

* generate sequences of fixed-frequency dates and time spans
* conform or convert time series to a particular frequency
* compute "relative" dates based on various non-standard time increments
  (e.g. 5 business days before the last business day of the year), or "roll"
  dates forward or backward

pandas provides a relatively compact and self-contained set of tools for
performing the above tasks.

Create a range of dates:

.. ipython:: python

   # 72 hours starting with midnight Jan 1st, 2011
   rng = pd.date_range('1/1/2011', periods=72, freq='H')
   rng[:5]

Index pandas objects with dates:

.. ipython:: python

   ts = pd.Series(np.random.randn(len(rng)), index=rng)
   ts.head()

Change frequency and fill gaps:

.. ipython:: python

   # to 45 minute frequency and forward fill
   converted = ts.asfreq('45Min', method='pad')
   converted.head()

Resample the series to a daily frequency:

.. ipython:: python

   # Daily means
   ts.resample('D').mean()


.. _timeseries.overview:

Overview
--------

The following table shows the type of time-related classes pandas can handle and
how to create them.

=================  =============================== ===================================================================
Class              Remarks                         How to create
=================  =============================== ===================================================================
``Timestamp``      Represents a single timestamp   ``to_datetime``, ``Timestamp``
``DatetimeIndex``  Index of ``Timestamp``          ``to_datetime``, ``date_range``, ``bdate_range``, ``DatetimeIndex``
``Period``         Represents a single time span   ``Period``
``PeriodIndex``    Index of ``Period``             ``period_range``, ``PeriodIndex``
=================  =============================== ===================================================================

.. _timeseries.representation:

Timestamps vs. Time Spans
-------------------------

Timestamped data is the most basic type of time series data that associates
values with points in time. For pandas objects it means using the points in
time.

.. ipython:: python

   pd.Timestamp(datetime(2012, 5, 1))
   pd.Timestamp('2012-05-01')
   pd.Timestamp(2012, 5, 1)

However, in many cases it is more natural to associate things like change
variables with a time span instead. The span represented by ``Period`` can be
specified explicitly, or inferred from datetime string format.

For example:

.. ipython:: python

   pd.Period('2011-01')

   pd.Period('2012-05', freq='D')

:class:`Timestamp` and :class:`Period` can serve as an index. Lists of 
``Timestamp`` and ``Period`` are automatically coerced to :class:`DatetimeIndex`
and :class:`PeriodIndex` respectively.

.. ipython:: python

   dates = [pd.Timestamp('2012-05-01'), pd.Timestamp('2012-05-02'), pd.Timestamp('2012-05-03')]
   ts = pd.Series(np.random.randn(3), dates)

   type(ts.index)
   ts.index

   ts

   periods = [pd.Period('2012-01'), pd.Period('2012-02'), pd.Period('2012-03')]

   ts = pd.Series(np.random.randn(3), periods)

   type(ts.index)
   ts.index

   ts

pandas allows you to capture both representations and
convert between them. Under the hood, pandas represents timestamps using
instances of ``Timestamp`` and sequences of timestamps using instances of
``DatetimeIndex``. For regular time spans, pandas uses ``Period`` objects for
scalar values and ``PeriodIndex`` for sequences of spans. Better support for
irregular intervals with arbitrary start and end points are forth-coming in
future releases.


.. _timeseries.converting:

Converting to Timestamps
------------------------

To convert a :class:`Series` or list-like object of date-like objects e.g. strings,
epochs, or a mixture, you can use the ``to_datetime`` function. When passed
a ``Series``, this returns a ``Series`` (with the same index), while a list-like
is converted to a ``DatetimeIndex``:

.. ipython:: python

    pd.to_datetime(pd.Series(['Jul 31, 2009', '2010-01-10', None]))

    pd.to_datetime(['2005/11/23', '2010.12.31'])

If you use dates which start with the day first (i.e. European style),
you can pass the ``dayfirst`` flag:

.. ipython:: python

    pd.to_datetime(['04-01-2012 10:00'], dayfirst=True)

    pd.to_datetime(['14-01-2012', '01-14-2012'], dayfirst=True)

.. warning::

   You see in the above example that ``dayfirst`` isn't strict, so if a date
   can't be parsed with the day being first it will be parsed as if
   ``dayfirst`` were False.

If you pass a single string to ``to_datetime``, it returns a single ``Timestamp``. 
``Timestamp`` can also accept string input, but it doesn't accept string parsing
options like ``dayfirst`` or ``format``, so use ``to_datetime`` if these are required.

.. ipython:: python

    pd.to_datetime('2010/11/12')

    pd.Timestamp('2010/11/12')

You can also use the ``DatetimeIndex`` constructor directly:

.. ipython:: python

    pd.DatetimeIndex(['2018-01-01', '2018-01-03', '2018-01-05'])

The string 'infer' can be passed in order to set the frequency of the index as the
inferred frequency upon creation:

.. ipython:: python

    pd.DatetimeIndex(['2018-01-01', '2018-01-03', '2018-01-05'], freq='infer')

Providing a Format Argument
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the required datetime string, a ``format`` argument can be passed to ensure specific parsing.
This could also potentially speed up the conversion considerably.

.. ipython:: python

    pd.to_datetime('2010/11/12', format='%Y/%m/%d')

    pd.to_datetime('12-11-2010 00:00', format='%d-%m-%Y %H:%M')

For more information on the choices available when specifying the ``format`` 
option, see the Python `datetime documentation`_.

.. _datetime documentation: https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior

Assembling Datetime from Multiple DataFrame Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.18.1

You can also pass a ``DataFrame`` of integer or string columns to assemble into a ``Series`` of ``Timestamps``.

.. ipython:: python

   df = pd.DataFrame({'year': [2015, 2016],
                      'month': [2, 3],
                      'day': [4, 5],
                      'hour': [2, 3]})
   pd.to_datetime(df)


You can pass only the columns that you need to assemble.

.. ipython:: python

   pd.to_datetime(df[['year', 'month', 'day']])

``pd.to_datetime`` looks for standard designations of the datetime component in the column names, including:

* required: ``year``, ``month``, ``day``
* optional: ``hour``, ``minute``, ``second``, ``millisecond``, ``microsecond``, ``nanosecond``

Invalid Data
~~~~~~~~~~~~

The default behavior, ``errors='raise'``, is to raise when unparseable:

.. code-block:: ipython

    In [2]: pd.to_datetime(['2009/07/31', 'asd'], errors='raise')
    ValueError: Unknown string format

Pass ``errors='ignore'`` to return the original input when unparseable:

.. ipython:: python

   pd.to_datetime(['2009/07/31', 'asd'], errors='ignore')

Pass ``errors='coerce'`` to convert unparseable data to ``NaT`` (not a time):

.. ipython:: python

   pd.to_datetime(['2009/07/31', 'asd'], errors='coerce')


.. _timeseries.converting.epoch:

Epoch Timestamps
~~~~~~~~~~~~~~~~

pandas supports converting integer or float epoch times to ``Timestamp`` and
``DatetimeIndex``. The default unit is nanoseconds, since that is how ``Timestamp``
objects are stored internally. However, epochs are often stored in another ``unit``
which can be specified. These are computed from the starting point specified by the
``origin`` parameter.

.. ipython:: python

   pd.to_datetime([1349720105, 1349806505, 1349892905,
                   1349979305, 1350065705], unit='s')

   pd.to_datetime([1349720105100, 1349720105200, 1349720105300,
                   1349720105400, 1349720105500 ], unit='ms')

.. note::

   Epoch times will be rounded to the nearest nanosecond.

.. warning::

   Conversion of float epoch times can lead to inaccurate and unexpected results.
   :ref:`Python floats <python:tut-fp-issues>` have about 15 digits precision in
   decimal. Rounding during conversion from float to high precision ``Timestamp`` is
   unavoidable. The only way to achieve exact precision is to use a fixed-width
   types (e.g. an int64).

   .. ipython:: python

      pd.to_datetime([1490195805.433, 1490195805.433502912], unit='s')
      pd.to_datetime(1490195805433502912, unit='ns')

.. seealso::

   :ref:`timeseries.origin`

.. _timeseries.converting.epoch_inverse:

From Timestamps to Epoch
~~~~~~~~~~~~~~~~~~~~~~~~

To invert the operation from above, namely, to convert from a ``Timestamp`` to a 'unix' epoch:

.. ipython:: python

   stamps = pd.date_range('2012-10-08 18:15:05', periods=4, freq='D')
   stamps

We subtract the epoch (midnight at January 1, 1970 UTC) and then floor divide by the
"unit" (1 second).

.. ipython:: python

   (stamps - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')

.. _timeseries.origin:

Using the ``origin`` Parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.20.0

Using the ``origin`` parameter, one can specify an alternative starting point for creation
of a ``DatetimeIndex``. For example, to use 1960-01-01 as the starting date:

.. ipython:: python

   pd.to_datetime([1, 2, 3], unit='D', origin=pd.Timestamp('1960-01-01'))

The default is set at ``origin='unix'``, which defaults to ``1970-01-01 00:00:00``.
Commonly called 'unix epoch' or POSIX time.

.. ipython:: python

   pd.to_datetime([1, 2, 3], unit='D')

.. _timeseries.daterange:

Generating Ranges of Timestamps
-------------------------------

To generate an index with timestamps, you can use either the ``DatetimeIndex`` or
``Index`` constructor and pass in a list of datetime objects:

.. ipython:: python

   dates = [datetime(2012, 5, 1), datetime(2012, 5, 2), datetime(2012, 5, 3)]

   # Note the frequency information
   index = pd.DatetimeIndex(dates)
   index

   # Automatically converted to DatetimeIndex
   index = pd.Index(dates)
   index

In practice this becomes very cumbersome because we often need a very long
index with a large number of timestamps. If we need timestamps on a regular
frequency, we can use the :func:`date_range` and :func:`bdate_range` functions
to create a ``DatetimeIndex``. The default frequency for ``date_range`` is a
**day** while the default for ``bdate_range`` is a **business day**:

.. ipython:: python

   start = datetime(2011, 1, 1)
   end = datetime(2012, 1, 1)

   index = pd.date_range(start, end)
   index

   index = pd.bdate_range(start, end)
   index

Convenience functions like ``date_range`` and ``bdate_range`` can utilize a
variety of :ref:`frequency aliases <timeseries.offset_aliases>`:

.. ipython:: python

   pd.date_range(start, periods=1000, freq='M')

   pd.bdate_range(start, periods=250, freq='BQS')

``date_range`` and ``bdate_range`` make it easy to generate a range of dates
using various combinations of parameters like ``start``, ``end``, ``periods``,
and ``freq``. The start and end dates are strictly inclusive, so dates outside
of those specified will not be generated:

.. ipython:: python

   pd.date_range(start, end, freq='BM')

   pd.date_range(start, end, freq='W')

   pd.bdate_range(end=end, periods=20)

   pd.bdate_range(start=start, periods=20)

.. versionadded:: 0.23.0

Specifying ``start``, ``end``, and ``periods`` will generate a range of evenly spaced
dates from ``start`` to ``end`` inclusively, with ``periods`` number of elements in the
resulting ``DatetimeIndex``:

.. ipython:: python

   pd.date_range('2018-01-01', '2018-01-05', periods=5)

   pd.date_range('2018-01-01', '2018-01-05', periods=10)

.. _timeseries.custom-freq-ranges:

Custom Frequency Ranges
~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   This functionality was originally exclusive to ``cdate_range``, which is
   deprecated as of version 0.21.0 in favor of ``bdate_range``.  Note that
   ``cdate_range`` only utilizes the ``weekmask`` and ``holidays`` parameters
   when custom business day, 'C', is passed as the frequency string. Support has 
   been expanded with ``bdate_range`` to work with any custom frequency string.

.. versionadded:: 0.21.0

``bdate_range`` can also generate a range of custom frequency dates by using
the ``weekmask`` and ``holidays`` parameters.  These parameters will only be
used if a custom frequency string is passed.

.. ipython:: python

   weekmask = 'Mon Wed Fri'

   holidays = [datetime(2011, 1, 5), datetime(2011, 3, 14)]

   pd.bdate_range(start, end, freq='C', weekmask=weekmask, holidays=holidays)

   pd.bdate_range(start, end, freq='CBMS', weekmask=weekmask)

.. seealso::

   :ref:`timeseries.custombusinessdays`

.. _timeseries.timestamp-limits:

Timestamp Limitations
---------------------

Since pandas represents timestamps in nanosecond resolution, the time span that
can be represented using a 64-bit integer is limited to approximately 584 years:

.. ipython:: python

   pd.Timestamp.min
   pd.Timestamp.max

.. seealso::

   :ref:`timeseries.oob`

.. _timeseries.datetimeindex:

Indexing
--------

One of the main uses for ``DatetimeIndex`` is as an index for pandas objects.
The ``DatetimeIndex`` class contains many time series related optimizations:

* A large range of dates for various offsets are pre-computed and cached
  under the hood in order to make generating subsequent date ranges very fast
  (just have to grab a slice).
* Fast shifting using the ``shift`` and ``tshift`` method on pandas objects.
* Unioning of overlapping ``DatetimeIndex`` objects with the same frequency is
  very fast (important for fast data alignment).
* Quick access to date fields via properties such as ``year``, ``month``, etc.
* Regularization functions like ``snap`` and very fast ``asof`` logic.

``DatetimeIndex`` objects have all the basic functionality of regular ``Index``
objects, and a smorgasbord of advanced time series specific methods for easy
frequency processing.

.. seealso::
    :ref:`Reindexing methods <basics.reindexing>`

.. note::

    While pandas does not force you to have a sorted date index, some of these
    methods may have unexpected or incorrect behavior if the dates are unsorted.

``DatetimeIndex`` can be used like a regular index and offers all of its
intelligent functionality like selection, slicing, etc.

.. ipython:: python

   rng = pd.date_range(start, end, freq='BM')
   ts = pd.Series(np.random.randn(len(rng)), index=rng)
   ts.index
   ts[:5].index
   ts[::2].index

.. _timeseries.partialindexing:

Partial String Indexing
~~~~~~~~~~~~~~~~~~~~~~~

Dates and strings that parse to timestamps can be passed as indexing parameters:

.. ipython:: python

   ts['1/31/2011']

   ts[datetime(2011, 12, 25):]

   ts['10/31/2011':'12/31/2011']

To provide convenience for accessing longer time series, you can also pass in
the year or year and month as strings:

.. ipython:: python

   ts['2011']

   ts['2011-6']

This type of slicing will work on a ``DataFrame`` with a ``DatetimeIndex`` as well. Since the
partial string selection is a form of label slicing, the endpoints **will be** included. This
would include matching times on an included date:

.. ipython:: python

   dft = pd.DataFrame(randn(100000,1),
                      columns=['A'],
                      index=pd.date_range('20130101',periods=100000,freq='T'))
   dft
   dft['2013']

This starts on the very first time in the month, and includes the last date and 
time for the month:

.. ipython:: python

   dft['2013-1':'2013-2']

This specifies a stop time **that includes all of the times on the last day**:

.. ipython:: python

   dft['2013-1':'2013-2-28']

This specifies an **exact** stop time (and is not the same as the above):

.. ipython:: python

   dft['2013-1':'2013-2-28 00:00:00']

We are stopping on the included end-point as it is part of the index:

.. ipython:: python

   dft['2013-1-15':'2013-1-15 12:30:00']

.. versionadded:: 0.18.0

``DatetimeIndex`` partial string indexing also works on a ``DataFrame`` with a ``MultiIndex``:

.. ipython:: python

   dft2 = pd.DataFrame(np.random.randn(20, 1),
                       columns=['A'],
                       index=pd.MultiIndex.from_product([pd.date_range('20130101',
                                                                       periods=10,
                                                                       freq='12H'),
                                                        ['a', 'b']]))
   dft2
   dft2.loc['2013-01-05']
   idx = pd.IndexSlice
   dft2 = dft2.swaplevel(0, 1).sort_index()
   dft2.loc[idx[:, '2013-01-05'], :]

.. _timeseries.slice_vs_exact_match:

Slice vs. Exact Match
~~~~~~~~~~~~~~~~~~~~~

.. versionchanged:: 0.20.0

The same string used as an indexing parameter can be treated either as a slice or as an exact match depending on the resolution of the index. If the string is less accurate than the index, it will be treated as a slice, otherwise as an exact match.

Consider a ``Series`` object with a minute resolution index:

.. ipython:: python

    series_minute = pd.Series([1, 2, 3],
                              pd.DatetimeIndex(['2011-12-31 23:59:00',
                                                '2012-01-01 00:00:00',
                                                '2012-01-01 00:02:00']))
    series_minute.index.resolution

A timestamp string less accurate than a minute gives a ``Series`` object.

.. ipython:: python

    series_minute['2011-12-31 23']

A timestamp string with minute resolution (or more accurate), gives a scalar instead, i.e. it is not casted to a slice.

.. ipython:: python

    series_minute['2011-12-31 23:59']
    series_minute['2011-12-31 23:59:00']

If index resolution is second, then the minute-accurate timestamp gives a 
``Series``.

.. ipython:: python

    series_second = pd.Series([1, 2, 3],
                              pd.DatetimeIndex(['2011-12-31 23:59:59',
                                                '2012-01-01 00:00:00',
                                                '2012-01-01 00:00:01']))
    series_second.index.resolution
    series_second['2011-12-31 23:59']

If the timestamp string is treated as a slice, it can be used to index ``DataFrame`` with ``[]`` as well.

.. ipython:: python

    dft_minute = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]},
                               index=series_minute.index)
    dft_minute['2011-12-31 23']


.. warning::

   However, if the string is treated as an exact match, the selection in ``DataFrame``'s ``[]`` will be column-wise and not row-wise, see :ref:`Indexing Basics <indexing.basics>`. For example ``dft_minute['2011-12-31 23:59']`` will raise ``KeyError`` as ``'2012-12-31 23:59'`` has the same resolution as the index and there is no column with such name:

   To *always* have unambiguous selection, whether the row is treated as a slice or a single selection, use ``.loc``.

   .. ipython:: python

     dft_minute.loc['2011-12-31 23:59']

Note also that ``DatetimeIndex`` resolution cannot be less precise than day.

.. ipython:: python

    series_monthly = pd.Series([1, 2, 3],
                              pd.DatetimeIndex(['2011-12',
                                                '2012-01',
                                                '2012-02']))
    series_monthly.index.resolution
    series_monthly['2011-12'] # returns Series


Exact Indexing
~~~~~~~~~~~~~~

As discussed in previous section, indexing a ``DatetimeIndex`` with a partial string depends on the "accuracy" of the period, in other words how specific the interval is in relation to the resolution of the index. In contrast, indexing with ``Timestamp`` or ``datetime`` objects is exact, because the objects have exact meaning. These also follow the semantics of *including both endpoints*.

These ``Timestamp`` and ``datetime`` objects have exact ``hours, minutes,`` and ``seconds``, even though they were not explicitly specified (they are ``0``).

.. ipython:: python

   dft[datetime(2013, 1, 1):datetime(2013,2,28)]

With no defaults.

.. ipython:: python

   dft[datetime(2013, 1, 1, 10, 12, 0):datetime(2013, 2, 28, 10, 12, 0)]


Truncating & Fancy Indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :meth:`~DataFrame.truncate` convenience function is provided that is similar 
to slicing. Note that ``truncate`` assumes a 0 value for any unspecified date 
component in a ``DatetimeIndex`` in contrast to slicing which returns any 
partially matching dates:

.. ipython:: python

   rng2 = pd.date_range('2011-01-01', '2012-01-01', freq='W')
   ts2 = pd.Series(np.random.randn(len(rng2)), index=rng2)

   ts2.truncate(before='2011-11', after='2011-12')
   ts2['2011-11':'2011-12']

Even complicated fancy indexing that breaks the ``DatetimeIndex`` frequency
regularity will result in a ``DatetimeIndex``, although frequency is lost:

.. ipython:: python

   ts2[[0, 2, 6]].index

.. _timeseries.components:

Time/Date Components
--------------------

There are several time/date properties that one can access from ``Timestamp`` or a collection of timestamps like a ``DatetimeIndex``.

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
    date,"Returns datetime.date (does not contain timezone information)"
    time,"Returns datetime.time (does not contain timezone information)"
    timetz,"Returns datetime.time as local time with timezone information"
    dayofyear,"The ordinal day of year"
    weekofyear,"The week ordinal of the year"
    week,"The week ordinal of the year"
    dayofweek,"The number of the day of the week with Monday=0, Sunday=6"
    weekday,"The number of the day of the week with Monday=0, Sunday=6"
    weekday_name,"The name of the day in a week (ex: Friday)"
    quarter,"Quarter of the date: Jan-Mar = 1, Apr-Jun = 2, etc."
    days_in_month,"The number of days in the month of the datetime"
    is_month_start,"Logical indicating if first day of month (defined by frequency)"
    is_month_end,"Logical indicating if last day of month (defined by frequency)"
    is_quarter_start,"Logical indicating if first day of quarter (defined by frequency)"
    is_quarter_end,"Logical indicating if last day of quarter (defined by frequency)"
    is_year_start,"Logical indicating if first day of year (defined by frequency)"
    is_year_end,"Logical indicating if last day of year (defined by frequency)"
    is_leap_year,"Logical indicating if the date belongs to a leap year"

Furthermore, if you have a ``Series`` with datetimelike values, then you can 
access these properties via the ``.dt`` accessor, as detailed in the section
on :ref:`.dt accessors<basics.dt_accessors>`.

.. _timeseries.offsets:

DateOffset Objects
------------------

In the preceding examples, we created ``DatetimeIndex`` objects at various
frequencies by passing in :ref:`frequency strings <timeseries.offset_aliases>`
like 'M', 'W', and 'BM' to the ``freq`` keyword. Under the hood, these frequency
strings are being translated into an instance of :class:`DateOffset`,
which represents a regular frequency increment. Specific offset logic like
"month", "business day", or "one hour" is represented in its various subclasses.

.. csv-table::
    :header: "Class name", "Description"
    :widths: 15, 65

    DateOffset, "Generic offset class, defaults to 1 calendar day"
    BDay, "business day (weekday)"
    CDay, "custom business day"
    Week, "one week, optionally anchored on a day of the week"
    WeekOfMonth, "the x-th day of the y-th week of each month"
    LastWeekOfMonth, "the x-th day of the last week of each month"
    MonthEnd, "calendar month end"
    MonthBegin, "calendar month begin"
    BMonthEnd, "business month end"
    BMonthBegin, "business month begin"
    CBMonthEnd, "custom business month end"
    CBMonthBegin, "custom business month begin"
    SemiMonthEnd, "15th (or other day_of_month) and calendar month end"
    SemiMonthBegin, "15th (or other day_of_month) and calendar month begin"
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
    BusinessHour, "business hour"
    CustomBusinessHour, "custom business hour"
    Hour, "one hour"
    Minute, "one minute"
    Second, "one second"
    Milli, "one millisecond"
    Micro, "one microsecond"
    Nano, "one nanosecond"

The basic ``DateOffset`` takes the same arguments as
``dateutil.relativedelta``, which works as follows:

.. ipython:: python

   d = datetime(2008, 8, 18, 9, 0)
   d + relativedelta(months=4, days=5)

We could have done the same thing with ``DateOffset``:

.. ipython:: python

   from pandas.tseries.offsets import *
   d + DateOffset(months=4, days=5)

The key features of a ``DateOffset`` object are:

* It can be added / subtracted to/from a datetime object to obtain a
  shifted date.
* It can be multiplied by an integer (positive or negative) so that the
  increment will be applied multiple times.
* It has :meth:`~pandas.DateOffset.rollforward` and
  :meth:`~pandas.DateOffset.rollback` methods for moving a date forward or 
  backward to the next or previous "offset date".

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

These operations (``apply``, ``rollforward`` and ``rollback``) preserve time 
(hour, minute, etc) information by default. To reset time, use ``normalize=True`` 
when creating the offset instance. If ``normalize=True``, the result is 
normalized after the function is applied.


.. ipython:: python

   day = Day()
   day.apply(pd.Timestamp('2014-01-01 09:00'))

   day = Day(normalize=True)
   day.apply(pd.Timestamp('2014-01-01 09:00'))

   hour = Hour()
   hour.apply(pd.Timestamp('2014-01-01 22:00'))

   hour = Hour(normalize=True)
   hour.apply(pd.Timestamp('2014-01-01 22:00'))
   hour.apply(pd.Timestamp('2014-01-01 23:00'))


.. _timeseries.dayvscalendarday:

Day vs. CalendarDay
~~~~~~~~~~~~~~~~~~~

:class:`Day` (``'D'``) is a timedelta-like offset that respects absolute time
arithmetic and essentially is an alias for 24 :class:`Hour`. This offset is the default
argument to many pandas time related function like :func:`date_range` and :func:`timedelta_range`.

:class:`CalendarDay` (``'CD'``) is a relativedelta-like offset that respects
calendar time arithmetic. :class:`CalendarDay` is useful preserving calendar day
semantics with date times with have day light savings transitions.

.. ipython:: python

   ts = pd.Timestamp('2016-10-30 00:00:00', tz='Europe/Helsinki')
   ts + pd.offsets.Day(1)
   ts + pd.offsets.CalendarDay(1)


Parametric Offsets
~~~~~~~~~~~~~~~~~~

Some of the offsets can be "parameterized" when created to result in different
behaviors. For example, the ``Week`` offset for generating weekly data accepts a
``weekday`` parameter which results in the generated dates always lying on a
particular day of the week:

.. ipython:: python

   d
   d + Week()
   d + Week(weekday=4)
   (d + Week(weekday=4)).weekday()

   d - Week()

The ``normalize`` option will be effective for addition and subtraction.

.. ipython:: python

   d + Week(normalize=True)
   d - Week(normalize=True)


Another example is parameterizing ``YearEnd`` with the specific ending month:

.. ipython:: python

   d + YearEnd()
   d + YearEnd(month=6)


.. _timeseries.offsetseries:

Using Offsets with ``Series`` / ``DatetimeIndex``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Offsets can be used with either a ``Series`` or ``DatetimeIndex`` to
apply the offset to each element.

.. ipython:: python

   rng = pd.date_range('2012-01-01', '2012-01-03')
   s = pd.Series(rng)
   rng
   rng + DateOffset(months=2)
   s + DateOffset(months=2)
   s - DateOffset(months=2)

If the offset class maps directly to a ``Timedelta`` (``Day``, ``Hour``,
``Minute``, ``Second``, ``Micro``, ``Milli``, ``Nano``) it can be
used exactly like a ``Timedelta`` - see the
:ref:`Timedelta section<timedeltas.operations>` for more examples.

.. ipython:: python

   s - Day(2)
   td = s - pd.Series(pd.date_range('2011-12-29', '2011-12-31'))
   td
   td + Minute(15)

Note that some offsets (such as ``BQuarterEnd``) do not have a
vectorized implementation.  They can still be used but may
calculate significantly slower and will show a ``PerformanceWarning``

.. ipython:: python
   :okwarning:

   rng + BQuarterEnd()


.. _timeseries.custombusinessdays:

Custom Business Days
~~~~~~~~~~~~~~~~~~~~

The ``CDay`` or ``CustomBusinessDay`` class provides a parametric
``BusinessDay`` class which can be used to create customized business day
calendars which account for local holidays and local weekend conventions.

As an interesting example, let's look at Egypt where a Friday-Saturday weekend is observed.

.. ipython:: python

    from pandas.tseries.offsets import CustomBusinessDay
    weekmask_egypt = 'Sun Mon Tue Wed Thu'

    # They also observe International Workers' Day so let's
    # add that for a couple of years

    holidays = ['2012-05-01', datetime(2013, 5, 1), np.datetime64('2014-05-01')]
    bday_egypt = CustomBusinessDay(holidays=holidays, weekmask=weekmask_egypt)
    dt = datetime(2013, 4, 30)
    dt + 2 * bday_egypt

Let's map to the weekday names:

.. ipython:: python

    dts = pd.date_range(dt, periods=5, freq=bday_egypt)

    pd.Series(dts.weekday, dts).map(pd.Series('Mon Tue Wed Thu Fri Sat Sun'.split()))

Holiday calendars can be used to provide the list of holidays.  See the
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
    pd.DatetimeIndex(start='20100101',end='20120101',freq=bmth_us)

.. note::

    The frequency string 'C' is used to indicate that a CustomBusinessDay
    DateOffset is used, it is important to note that since CustomBusinessDay is
    a parameterised type, instances of CustomBusinessDay may differ and this is
    not detectable from the 'C' frequency string. The user therefore needs to
    ensure that the 'C' frequency string is used consistently within the user's
    application.

.. _timeseries.businesshour:

Business Hour
~~~~~~~~~~~~~

The ``BusinessHour`` class provides a business hour representation on ``BusinessDay``,
allowing to use specific start and end times.

By default, ``BusinessHour`` uses 9:00 - 17:00 as business hours.
Adding ``BusinessHour`` will increment ``Timestamp`` by hourly frequency.
If target ``Timestamp`` is out of business hours, move to the next business hour 
then increment it. If the result exceeds the business hours end, the remaining 
hours are added to the next business day.

.. ipython:: python

    bh = BusinessHour()
    bh

    # 2014-08-01 is Friday
    pd.Timestamp('2014-08-01 10:00').weekday()
    pd.Timestamp('2014-08-01 10:00') + bh

    # Below example is the same as: pd.Timestamp('2014-08-01 09:00') + bh
    pd.Timestamp('2014-08-01 08:00') + bh

    # If the results is on the end time, move to the next business day
    pd.Timestamp('2014-08-01 16:00') + bh

    # Remainings are added to the next day
    pd.Timestamp('2014-08-01 16:30') + bh

    # Adding 2 business hours
    pd.Timestamp('2014-08-01 10:00') + BusinessHour(2)

    # Subtracting 3 business hours
    pd.Timestamp('2014-08-01 10:00') + BusinessHour(-3)

You can also specify ``start`` and ``end`` time by keywords. The argument must 
be a ``str`` with an ``hour:minute`` representation or a ``datetime.time`` 
instance. Specifying seconds, microseconds and nanoseconds as business hour 
results in ``ValueError``.

.. ipython:: python

    bh = BusinessHour(start='11:00', end=time(20, 0))
    bh

    pd.Timestamp('2014-08-01 13:00') + bh
    pd.Timestamp('2014-08-01 09:00') + bh
    pd.Timestamp('2014-08-01 18:00') + bh

Passing ``start`` time later than ``end`` represents midnight business hour.
In this case, business hour exceeds midnight and overlap to the next day.
Valid business hours are distinguished by whether it started from valid ``BusinessDay``.

.. ipython:: python

    bh = BusinessHour(start='17:00', end='09:00')
    bh

    pd.Timestamp('2014-08-01 17:00') + bh
    pd.Timestamp('2014-08-01 23:00') + bh

    # Although 2014-08-02 is Satuaday,
    # it is valid because it starts from 08-01 (Friday).
    pd.Timestamp('2014-08-02 04:00') + bh

    # Although 2014-08-04 is Monday,
    # it is out of business hours because it starts from 08-03 (Sunday).
    pd.Timestamp('2014-08-04 04:00') + bh

Applying ``BusinessHour.rollforward`` and ``rollback`` to out of business hours results in
the next business hour start or previous day's end. Different from other offsets, ``BusinessHour.rollforward``
may output different results from ``apply`` by definition.

This is because one day's business hour end is equal to next day's business hour start. For example,
under the default business hours (9:00 - 17:00), there is no gap (0 minutes) between ``2014-08-01 17:00`` and
``2014-08-04 09:00``.

.. ipython:: python

    # This adjusts a Timestamp to business hour edge
    BusinessHour().rollback(pd.Timestamp('2014-08-02 15:00'))
    BusinessHour().rollforward(pd.Timestamp('2014-08-02 15:00'))

    # It is the same as BusinessHour().apply(pd.Timestamp('2014-08-01 17:00')).
    # And it is the same as BusinessHour().apply(pd.Timestamp('2014-08-04 09:00'))
    BusinessHour().apply(pd.Timestamp('2014-08-02 15:00'))

    # BusinessDay results (for reference)
    BusinessHour().rollforward(pd.Timestamp('2014-08-02'))

    # It is the same as BusinessDay().apply(pd.Timestamp('2014-08-01'))
    # The result is the same as rollworward because BusinessDay never overlap.
    BusinessHour().apply(pd.Timestamp('2014-08-02'))

``BusinessHour`` regards Saturday and Sunday as holidays. To use arbitrary 
holidays, you can use ``CustomBusinessHour`` offset, as explained in the 
following subsection.

.. _timeseries.custombusinesshour:

Custom Business Hour
~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.18.1

The ``CustomBusinessHour`` is a mixture of ``BusinessHour`` and ``CustomBusinessDay`` which
allows you to specify arbitrary holidays. ``CustomBusinessHour`` works as the same
as ``BusinessHour`` except that it skips specified custom holidays.

.. ipython:: python

    from pandas.tseries.holiday import USFederalHolidayCalendar
    bhour_us = CustomBusinessHour(calendar=USFederalHolidayCalendar())
    # Friday before MLK Day
    dt = datetime(2014, 1, 17, 15)

    dt + bhour_us

    # Tuesday after MLK Day (Monday is skipped because it's a holiday)
    dt + bhour_us * 2

You can use keyword arguments supported by either ``BusinessHour`` and ``CustomBusinessDay``.

.. ipython:: python

    bhour_mon = CustomBusinessHour(start='10:00', weekmask='Tue Wed Thu Fri')

    # Monday is skipped because it's a holiday, business hour starts from 10:00
    dt + bhour_mon * 2

.. _timeseries.offset_aliases:

Offset Aliases
~~~~~~~~~~~~~~

A number of string aliases are given to useful common time series
frequencies. We will refer to these aliases as *offset aliases*.

.. csv-table::
    :header: "Alias", "Description"
    :widths: 15, 100

    "B", "business day frequency"
    "C", "custom business day frequency"
    "D", "day frequency"
    "CD", "calendar day frequency"
    "W", "weekly frequency"
    "M", "month end frequency"
    "SM", "semi-month end frequency (15th and end of month)"
    "BM", "business month end frequency"
    "CBM", "custom business month end frequency"
    "MS", "month start frequency"
    "SMS", "semi-month start frequency (1st and 15th)"
    "BMS", "business month start frequency"
    "CBMS", "custom business month start frequency"
    "Q", "quarter end frequency"
    "BQ", "business quarter end frequency"
    "QS", "quarter start frequency"
    "BQS", "business quarter start frequency"
    "A, Y", "year end frequency"
    "BA, BY", "business year end frequency"
    "AS, YS", "year start frequency"
    "BAS, BYS", "business year start frequency"
    "BH", "business hour frequency"
    "H", "hourly frequency"
    "T, min", "minutely frequency"
    "S", "secondly frequency"
    "L, ms", "milliseconds"
    "U, us", "microseconds"
    "N", "nanoseconds"

Combining Aliases
~~~~~~~~~~~~~~~~~

As we have seen previously, the alias and the offset instance are fungible in
most functions:

.. ipython:: python

   pd.date_range(start, periods=5, freq='B')

   pd.date_range(start, periods=5, freq=BDay())

You can combine together day and intraday offsets:

.. ipython:: python

   pd.date_range(start, periods=10, freq='2h20min')

   pd.date_range(start, periods=10, freq='1D10U')

Anchored Offsets
~~~~~~~~~~~~~~~~

For some frequencies you can specify an anchoring suffix:

.. csv-table::
    :header: "Alias", "Description"
    :widths: 15, 100

    "W\-SUN", "weekly frequency (Sundays). Same as 'W'"
    "W\-MON", "weekly frequency (Mondays)"
    "W\-TUE", "weekly frequency (Tuesdays)"
    "W\-WED", "weekly frequency (Wednesdays)"
    "W\-THU", "weekly frequency (Thursdays)"
    "W\-FRI", "weekly frequency (Fridays)"
    "W\-SAT", "weekly frequency (Saturdays)"
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

Anchored Offset Semantics
~~~~~~~~~~~~~~~~~~~~~~~~~

For those offsets that are anchored to the start or end of specific
frequency (``MonthEnd``, ``MonthBegin``, ``WeekEnd``, etc), the following
rules apply to rolling forward and backwards.

When ``n`` is not 0, if the given date is not on an anchor point, it snapped to the next(previous)
anchor point, and moved ``|n|-1`` additional steps forwards or backwards.

.. ipython:: python

   pd.Timestamp('2014-01-02') + MonthBegin(n=1)
   pd.Timestamp('2014-01-02') + MonthEnd(n=1)

   pd.Timestamp('2014-01-02') - MonthBegin(n=1)
   pd.Timestamp('2014-01-02') - MonthEnd(n=1)

   pd.Timestamp('2014-01-02') + MonthBegin(n=4)
   pd.Timestamp('2014-01-02') - MonthBegin(n=4)

If the given date *is* on an anchor point, it is moved ``|n|`` points forwards
or backwards.

.. ipython:: python

   pd.Timestamp('2014-01-01') + MonthBegin(n=1)
   pd.Timestamp('2014-01-31') + MonthEnd(n=1)

   pd.Timestamp('2014-01-01') - MonthBegin(n=1)
   pd.Timestamp('2014-01-31') - MonthEnd(n=1)

   pd.Timestamp('2014-01-01') + MonthBegin(n=4)
   pd.Timestamp('2014-01-31') - MonthBegin(n=4)

For the case when ``n=0``, the date is not moved if on an anchor point, otherwise
it is rolled forward to the next anchor point.

.. ipython:: python

   pd.Timestamp('2014-01-02') + MonthBegin(n=0)
   pd.Timestamp('2014-01-02') + MonthEnd(n=0)

   pd.Timestamp('2014-01-01') + MonthBegin(n=0)
   pd.Timestamp('2014-01-31') + MonthEnd(n=0)

.. _timeseries.holiday:

Holidays / Holiday Calendars
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Holidays and calendars provide a simple way to define holiday rules to be used
with ``CustomBusinessDay`` or in other analysis that requires a predefined
set of holidays.  The ``AbstractHolidayCalendar`` class provides all the necessary
methods to return a list of holidays and only ``rules`` need to be defined
in a specific holiday calendar class. Furthermore, the ``start_date`` and ``end_date``
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
and holidays (i.e., Memorial Day/July 4th).  For example, the below defines
a custom business day offset using the ``ExampleCalendar``.  Like any other offset,
it can be used to create a ``DatetimeIndex`` or added to ``datetime``
or ``Timestamp`` objects.

.. ipython:: python

    from pandas.tseries.offsets import CDay
    pd.DatetimeIndex(start='7/1/2012', end='7/10/2012',
        freq=CDay(calendar=cal)).to_pydatetime()
    offset = CustomBusinessDay(calendar=cal)
    datetime(2012, 5, 25) + offset
    datetime(2012, 7, 3) + offset
    datetime(2012, 7, 3) + 2 * offset
    datetime(2012, 7, 6) + offset

Ranges are defined by the ``start_date`` and ``end_date`` class attributes
of ``AbstractHolidayCalendar``.  The defaults are shown below.

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

Time Series-Related Instance Methods
------------------------------------

Shifting / Lagging
~~~~~~~~~~~~~~~~~~

One may want to *shift* or *lag* the values in a time series back and forward in
time. The method for this is :meth:`~Series.shift`, which is available on all of 
the pandas objects.

.. ipython:: python

   ts = ts[:5]
   ts.shift(1)

The ``shift`` method accepts an ``freq`` argument which can accept a
``DateOffset`` class or other ``timedelta``-like object or also an 
:ref:`offset alias <timeseries.offset_aliases>`:

.. ipython:: python

   ts.shift(5, freq=offsets.BDay())
   ts.shift(5, freq='BM')

Rather than changing the alignment of the data and the index, ``DataFrame`` and
``Series`` objects also have a :meth:`~Series.tshift` convenience method that 
changes all the dates in the index by a specified number of offsets:

.. ipython:: python

   ts.tshift(5, freq='D')

Note that with ``tshift``, the leading entry is no longer NaN because the data
is not being realigned.

Frequency Conversion
~~~~~~~~~~~~~~~~~~~~

The primary function for changing frequencies is the :meth:`~Series.asfreq` 
method. For a ``DatetimeIndex``, this is basically just a thin, but convenient 
wrapper around :meth:`~Series.reindex`  which generates a ``date_range`` and 
calls ``reindex``.

.. ipython:: python

   dr = pd.date_range('1/1/2010', periods=3, freq=3 * offsets.BDay())
   ts = pd.Series(randn(3), index=dr)
   ts
   ts.asfreq(BDay())

``asfreq`` provides a further convenience so you can specify an interpolation
method for any gaps that may appear after the frequency conversion.

.. ipython:: python

   ts.asfreq(BDay(), method='pad')

Filling Forward / Backward
~~~~~~~~~~~~~~~~~~~~~~~~~~

Related to ``asfreq`` and ``reindex`` is :meth:`~Series.fillna`, which is 
documented in the :ref:`missing data section <missing_data.fillna>`.

Converting to Python Datetimes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``DatetimeIndex`` can be converted to an array of Python native 
:py:class:`datetime.datetime` objects using the ``to_pydatetime`` method.

.. _timeseries.resampling:

Resampling
----------

.. warning::

   The interface to ``.resample`` has changed in 0.18.0 to be more groupby-like and hence more flexible.
   See the :ref:`whatsnew docs <whatsnew_0180.breaking.resample>` for a comparison with prior versions.

Pandas has a simple, powerful, and efficient functionality for performing 
resampling operations during frequency conversion (e.g., converting secondly 
data into 5-minutely data). This is extremely common in, but not limited to, 
financial applications.

:meth:`~Series.resample` is a time-based groupby, followed by a reduction method 
on each of its groups. See some :ref:`cookbook examples <cookbook.resample>` for 
some advanced strategies.

Starting in version 0.18.1, the ``resample()`` function can be used directly from
``DataFrameGroupBy`` objects, see the :ref:`groupby docs <groupby.transform.window_resample>`.

.. note::

   ``.resample()`` is similar to using a :meth:`~Series.rolling` operation with 
   a time-based offset, see a discussion :ref:`here <stats.moments.ts-versus-resampling>`.

Basics
~~~~~~

.. ipython:: python

   rng = pd.date_range('1/1/2012', periods=100, freq='S')

   ts = pd.Series(np.random.randint(0, 500, len(rng)), index=rng)

   ts.resample('5Min').sum()

The ``resample`` function is very flexible and allows you to specify many
different parameters to control the frequency conversion and resampling
operation.

Any function available via :ref:`dispatching <groupby.dispatch>` is available as
a method of the returned object, including ``sum``, ``mean``, ``std``, ``sem``,
``max``, ``min``, ``median``, ``first``, ``last``, ``ohlc``:

.. ipython:: python

   ts.resample('5Min').mean()

   ts.resample('5Min').ohlc()

   ts.resample('5Min').max()


For downsampling, ``closed`` can be set to 'left' or 'right' to specify which
end of the interval is closed:

.. ipython:: python

   ts.resample('5Min', closed='right').mean()

   ts.resample('5Min', closed='left').mean()

Parameters like ``label`` and ``loffset`` are used to manipulate the resulting
labels. ``label`` specifies whether the result is labeled with the beginning or
the end of the interval. ``loffset`` performs a time adjustment on the output
labels.

.. ipython:: python

   ts.resample('5Min').mean() # by default label='left'

   ts.resample('5Min', label='left').mean()

   ts.resample('5Min', label='left', loffset='1s').mean()

.. note::

    The default values for ``label`` and ``closed`` is 'left' for all 
    frequency offsets except for 'M', 'A', 'Q', 'BM', 'BA', 'BQ', and 'W' 
    which all have a default of 'right'.

    .. ipython:: python

       rng2 = pd.date_range('1/1/2012', end='3/31/2012', freq='D')
       ts2 = pd.Series(range(len(rng2)), index=rng2)

       # default: label='right', closed='right'
       ts2.resample('M').max()

       # default: label='left', closed='left'
       ts2.resample('SM').max()

       ts2.resample('SM', label='right', closed='right').max()

The ``axis`` parameter can be set to 0 or 1 and allows you to resample the
specified axis for a ``DataFrame``.

``kind`` can be set to 'timestamp' or 'period' to convert the resulting index
to/from timestamp and time span representations. By default ``resample``
retains the input representation.

``convention`` can be set to 'start' or 'end' when resampling period data
(detail below). It specifies how low frequency periods are converted to higher
frequency periods.


Upsampling
~~~~~~~~~~

For upsampling, you can specify a way to upsample and the ``limit`` parameter to interpolate over the gaps that are created:

.. ipython:: python

   # from secondly to every 250 milliseconds

   ts[:2].resample('250L').asfreq()

   ts[:2].resample('250L').ffill()

   ts[:2].resample('250L').ffill(limit=2)

Sparse Resampling
~~~~~~~~~~~~~~~~~

Sparse timeseries are the ones where you have a lot fewer points relative
to the amount of time you are looking to resample. Naively upsampling a sparse 
series can potentially generate lots of intermediate values. When you don't want 
to use a method to fill these values, e.g. ``fill_method`` is ``None``, then 
intermediate values will be filled with ``NaN``.

Since ``resample`` is a time-based groupby, the following is a method to efficiently
resample only the groups that are not all ``NaN``.

.. ipython:: python

    rng = pd.date_range('2014-1-1', periods=100, freq='D') + pd.Timedelta('1s')
    ts = pd.Series(range(100), index=rng)

If we want to resample to the full range of the series:

.. ipython:: python

    ts.resample('3T').sum()

We can instead only resample those groups where we have points as follows:

.. ipython:: python

    from functools import partial
    from pandas.tseries.frequencies import to_offset

    def round(t, freq):
        # round a Timestamp to a specified freq
        freq = to_offset(freq)
        return pd.Timestamp((t.value // freq.delta.value) * freq.delta.value)

    ts.groupby(partial(round, freq='3T')).sum()

.. _timeseries.aggregate:

Aggregation
~~~~~~~~~~~

Similar to the :ref:`aggregating API <basics.aggregate>`, :ref:`groupby API <groupby.aggregate>`, and the :ref:`window functions API <stats.aggregate>`,
a ``Resampler`` can be selectively resampled.

Resampling a ``DataFrame``, the default will be to act on all columns with the same function.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 3),
                     index=pd.date_range('1/1/2012', freq='S', periods=1000),
                     columns=['A', 'B', 'C'])
   r = df.resample('3T')
   r.mean()

We can select a specific column or columns using standard getitem.

.. ipython:: python

   r['A'].mean()

   r[['A','B']].mean()

You can pass a list or dict of functions to do aggregation with, outputting a ``DataFrame``:

.. ipython:: python

   r['A'].agg([np.sum, np.mean, np.std])

On a resampled ``DataFrame``, you can pass a list of functions to apply to each
column, which produces an aggregated result with a hierarchical index:

.. ipython:: python

   r.agg([np.sum, np.mean])

By passing a dict to ``aggregate`` you can apply a different aggregation to the
columns of a ``DataFrame``:

.. ipython:: python
   :okexcept:

   r.agg({'A' : np.sum,
          'B' : lambda x: np.std(x, ddof=1)})

The function names can also be strings. In order for a string to be valid it
must be implemented on the resampled object:

.. ipython:: python

   r.agg({'A' : 'sum', 'B' : 'std'})

Furthermore, you can also specify multiple aggregation functions for each column separately.

.. ipython:: python

   r.agg({'A' : ['sum','std'], 'B' : ['mean','std'] })


If a ``DataFrame`` does not have a datetimelike index, but instead you want
to resample based on datetimelike column in the frame, it can passed to the
``on`` keyword.

.. ipython:: python

   df = pd.DataFrame({'date': pd.date_range('2015-01-01', freq='W', periods=5),
                      'a': np.arange(5)},
                     index=pd.MultiIndex.from_arrays([
                              [1,2,3,4,5],
                              pd.date_range('2015-01-01', freq='W', periods=5)],
                          names=['v','d']))
   df
   df.resample('M', on='date').sum()

Similarly, if you instead want to resample by a datetimelike
level of ``MultiIndex``, its name or location can be passed to the
``level`` keyword.

.. ipython:: python

   df.resample('M', level='d').sum()


.. _timeseries.periods:

Time Span Representation
------------------------

Regular intervals of time are represented by ``Period`` objects in pandas while
sequences of ``Period`` objects are collected in a ``PeriodIndex``, which can
be created with the convenience function ``period_range``.

Period
~~~~~~

A ``Period`` represents a span of time (e.g., a day, a month, a quarter, etc).
You can specify the span via ``freq`` keyword using a frequency alias like below.
Because ``freq`` represents a span of ``Period``, it cannot be negative like "-3D".

.. ipython:: python

   pd.Period('2012', freq='A-DEC')

   pd.Period('2012-1-1', freq='D')

   pd.Period('2012-1-1 19:00', freq='H')

   pd.Period('2012-1-1 19:00', freq='5H')

Adding and subtracting integers from periods shifts the period by its own
frequency. Arithmetic is not allowed between ``Period`` with different ``freq`` (span).

.. ipython:: python

   p = pd.Period('2012', freq='A-DEC')
   p + 1
   p - 3
   p = pd.Period('2012-01', freq='2M')
   p + 2
   p - 1
   @okexcept
   p == pd.Period('2012-01', freq='3M')


If ``Period`` freq is daily or higher (``D``, ``H``, ``T``, ``S``, ``L``, ``U``, ``N``), ``offsets`` and ``timedelta``-like can be added if the result can have the same freq. Otherwise, ``ValueError`` will be raised.

.. ipython:: python

   p = pd.Period('2014-07-01 09:00', freq='H')
   p + Hour(2)
   p + timedelta(minutes=120)
   p + np.timedelta64(7200, 's')

.. code-block:: ipython

   In [1]: p + Minute(5)
   Traceback
      ...
   ValueError: Input has different freq from Period(freq=H)

If ``Period`` has other frequencies, only the same ``offsets`` can be added. Otherwise, ``ValueError`` will be raised.

.. ipython:: python

   p = pd.Period('2014-07', freq='M')
   p + MonthEnd(3)

.. code-block:: ipython

   In [1]: p + MonthBegin(3)
   Traceback
      ...
   ValueError: Input has different freq from Period(freq=M)

Taking the difference of ``Period`` instances with the same frequency will
return the number of frequency units between them:

.. ipython:: python

   pd.Period('2012', freq='A-DEC') - pd.Period('2002', freq='A-DEC')

PeriodIndex and period_range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Regular sequences of ``Period`` objects can be collected in a ``PeriodIndex``,
which can be constructed using the ``period_range`` convenience function:

.. ipython:: python

   prng = pd.period_range('1/1/2011', '1/1/2012', freq='M')
   prng

The ``PeriodIndex`` constructor can also be used directly:

.. ipython:: python

   pd.PeriodIndex(['2011-1', '2011-2', '2011-3'], freq='M')

Passing multiplied frequency outputs a sequence of ``Period`` which
has multiplied span.

.. ipython:: python

   pd.PeriodIndex(start='2014-01', freq='3M', periods=4)

If ``start`` or ``end`` are ``Period`` objects, they will be used as anchor
endpoints for a ``PeriodIndex`` with frequency matching that of the
``PeriodIndex`` constructor.

.. ipython:: python

   pd.PeriodIndex(start=pd.Period('2017Q1', freq='Q'),
                  end=pd.Period('2017Q2', freq='Q'), freq='M')

Just like ``DatetimeIndex``, a ``PeriodIndex`` can also be used to index pandas
objects:

.. ipython:: python

   ps = pd.Series(np.random.randn(len(prng)), prng)
   ps

``PeriodIndex`` supports addition and subtraction with the same rule as ``Period``.

.. ipython:: python

   idx = pd.period_range('2014-07-01 09:00', periods=5, freq='H')
   idx
   idx + Hour(2)

   idx = pd.period_range('2014-07', periods=5, freq='M')
   idx
   idx + MonthEnd(3)

``PeriodIndex`` has its own dtype named ``period``, refer to :ref:`Period Dtypes <timeseries.period_dtype>`.

.. _timeseries.period_dtype:

Period Dtypes
~~~~~~~~~~~~~

.. versionadded:: 0.19.0

``PeriodIndex`` has a custom ``period`` dtype. This is a pandas extension
dtype similar to the :ref:`timezone aware dtype <timeseries.timezone_series>` (``datetime64[ns, tz]``).

The ``period`` dtype holds the ``freq`` attribute and is represented with
``period[freq]`` like ``period[D]`` or ``period[M]``, using :ref:`frequency strings <timeseries.offset_aliases>`.

.. ipython:: python

   pi = pd.period_range('2016-01-01', periods=3, freq='M')
   pi
   pi.dtype

The ``period`` dtype can be used in ``.astype(...)``. It allows one to change the
``freq`` of a ``PeriodIndex`` like ``.asfreq()`` and convert a
``DatetimeIndex`` to ``PeriodIndex`` like ``to_period()``:

.. ipython:: python

   # change monthly freq to daily freq
   pi.astype('period[D]')

   # convert to DatetimeIndex
   pi.astype('datetime64[ns]')

   # convert to PeriodIndex
   dti = pd.date_range('2011-01-01', freq='M', periods=3)
   dti
   dti.astype('period[M]')


PeriodIndex Partial String Indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can pass in dates and strings to ``Series`` and ``DataFrame`` with ``PeriodIndex``, in the same manner as ``DatetimeIndex``. For details, refer to :ref:`DatetimeIndex Partial String Indexing <timeseries.partialindexing>`.

.. ipython:: python

   ps['2011-01']

   ps[datetime(2011, 12, 25):]

   ps['10/31/2011':'12/31/2011']

Passing a string representing a lower frequency than ``PeriodIndex`` returns partial sliced data.

.. ipython:: python

   ps['2011']

   dfp = pd.DataFrame(np.random.randn(600,1),
                      columns=['A'],
                      index=pd.period_range('2013-01-01 9:00', periods=600, freq='T'))
   dfp
   dfp['2013-01-01 10H']

As with ``DatetimeIndex``, the endpoints will be included in the result. The example below slices data starting from 10:00 to 11:59.

.. ipython:: python

   dfp['2013-01-01 10H':'2013-01-01 11H']

Frequency Conversion and Resampling with PeriodIndex
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The frequency of ``Period`` and ``PeriodIndex`` can be converted via the ``asfreq``
method. Let's start with the fiscal year 2011, ending in December:

.. ipython:: python

   p = pd.Period('2011', freq='A-DEC')
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

   p = pd.Period('2011-12', freq='M')

   p.asfreq('A-NOV')

Note that since we converted to an annual frequency that ends the year in
November, the monthly period of December 2011 is actually in the 2012 A-NOV
period.

.. _timeseries.quarterly:

Period conversions with anchored frequencies are particularly useful for
working with various quarterly data common to economics, business, and other
fields. Many organizations define quarters relative to the month in which their
fiscal year starts and ends. Thus, first quarter of 2011 could start in 2010 or
a few months into 2011. Via anchored frequencies, pandas works for all quarterly
frequencies ``Q-JAN`` through ``Q-DEC``.

``Q-DEC`` define regular calendar quarters:

.. ipython:: python

   p = pd.Period('2012Q1', freq='Q-DEC')

   p.asfreq('D', 's')

   p.asfreq('D', 'e')

``Q-MAR`` defines fiscal year end in March:

.. ipython:: python

   p = pd.Period('2011Q4', freq='Q-MAR')

   p.asfreq('D', 's')

   p.asfreq('D', 'e')

.. _timeseries.interchange:

Converting Between Representations
----------------------------------

Timestamped data can be converted to PeriodIndex-ed data using ``to_period``
and vice-versa using ``to_timestamp``:

.. ipython:: python

   rng = pd.date_range('1/1/2012', periods=5, freq='M')

   ts = pd.Series(np.random.randn(len(rng)), index=rng)

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

   prng = pd.period_range('1990Q1', '2000Q4', freq='Q-NOV')

   ts = pd.Series(np.random.randn(len(prng)), prng)

   ts.index = (prng.asfreq('M', 'e') + 1).asfreq('H', 's') + 9

   ts.head()

.. _timeseries.oob:

Representing Out-of-Bounds Spans
--------------------------------

If you have data that is outside of the ``Timestamp`` bounds, see :ref:`Timestamp limitations <timeseries.timestamp-limits>`,
then you can use a ``PeriodIndex`` and/or ``Series`` of ``Periods`` to do computations.

.. ipython:: python

   span = pd.period_range('1215-01-01', '1381-01-01', freq='D')
   span

To convert from an ``int64`` based YYYYMMDD representation.

.. ipython:: python

   s = pd.Series([20121231, 20141130, 99991231])
   s

   def conv(x):
       return pd.Period(year = x // 10000, month = x//100 % 100, day = x%100, freq='D')

   s.apply(conv)
   s.apply(conv)[2]

These can easily be converted to a ``PeriodIndex``:

.. ipython:: python

   span = pd.PeriodIndex(s.apply(conv))
   span

.. _timeseries.timezone:

Time Zone Handling
------------------

Pandas provides rich support for working with timestamps in different time
zones using ``pytz`` and ``dateutil`` libraries. ``dateutil`` currently is only
supported for fixed offset and tzfile zones. The default library is ``pytz``.
Support for ``dateutil`` is provided for compatibility with other
applications e.g. if you use ``dateutil`` in other Python packages.

Working with Time Zones
~~~~~~~~~~~~~~~~~~~~~~~

By default, pandas objects are time zone unaware:

.. ipython:: python

   rng = pd.date_range('3/6/2012 00:00', periods=15, freq='D')
   rng.tz is None

To supply the time zone, you can use the ``tz`` keyword to ``date_range`` and
other functions. Dateutil time zone strings are distinguished from ``pytz``
time zones by starting with ``dateutil/``.

* In ``pytz`` you can find a list of common (and less common) time zones using
  ``from pytz import common_timezones, all_timezones``.
* ``dateutil`` uses the OS timezones so there isn't a fixed list available. For
  common zones, the names are the same as ``pytz``.

.. ipython:: python

   # pytz
   rng_pytz = pd.date_range('3/6/2012 00:00', periods=10, freq='D',
                            tz='Europe/London')
   rng_pytz.tz

   # dateutil
   rng_dateutil = pd.date_range('3/6/2012 00:00', periods=10, freq='D',
                                tz='dateutil/Europe/London')
   rng_dateutil.tz

   # dateutil - utc special case
   rng_utc = pd.date_range('3/6/2012 00:00', periods=10, freq='D',
                           tz=dateutil.tz.tzutc())
   rng_utc.tz

Note that the ``UTC`` timezone is a special case in ``dateutil`` and should be constructed explicitly
as an instance of ``dateutil.tz.tzutc``. You can also construct other timezones explicitly first,
which gives you more control over which time zone is used:

.. ipython:: python

   # pytz
   tz_pytz = pytz.timezone('Europe/London')
   rng_pytz = pd.date_range('3/6/2012 00:00', periods=10, freq='D',
                            tz=tz_pytz)
   rng_pytz.tz == tz_pytz

   # dateutil
   tz_dateutil = dateutil.tz.gettz('Europe/London')
   rng_dateutil = pd.date_range('3/6/2012 00:00', periods=10, freq='D',
                                tz=tz_dateutil)
   rng_dateutil.tz == tz_dateutil

Timestamps, like Python's ``datetime.datetime`` object can be either time zone
naive or time zone aware. Naive time series and ``DatetimeIndex`` objects can be
*localized* using ``tz_localize``:

.. ipython:: python

   ts = pd.Series(np.random.randn(len(rng)), rng)

   ts_utc = ts.tz_localize('UTC')
   ts_utc

Again, you can explicitly construct the timezone object first.
You can use the ``tz_convert`` method to convert pandas objects to convert
tz-aware data to another time zone:

.. ipython:: python

   ts_utc.tz_convert('US/Eastern')

.. warning::

	Be wary of conversions between libraries. For some zones ``pytz`` and ``dateutil`` have different
	definitions of the zone. This is more of a problem for unusual timezones than for
	'standard' zones like ``US/Eastern``.

.. warning::

       Be aware that a timezone definition across versions of timezone libraries may not
       be considered equal.  This may cause problems when working with stored data that
       is localized using one version and operated on with a different version.
       See :ref:`here<io.hdf5-notes>` for how to handle such a situation.

.. warning::

       It is incorrect to pass a timezone directly into the ``datetime.datetime`` constructor (e.g.,
       ``datetime.datetime(2011, 1, 1, tz=timezone('US/Eastern'))``.  Instead, the datetime
       needs to be localized using the localize method on the timezone.

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

Like ``Series``, ``DataFrame``, and ``DatetimeIndex``; ``Timestamp`` objects
can be converted to other time zones using ``tz_convert``:

.. ipython:: python

   rng_eastern[5]
   rng_berlin[5]
   rng_eastern[5].tz_convert('Europe/Berlin')

Localization of ``Timestamp`` functions just like ``DatetimeIndex`` and ``Series``:

.. ipython:: python

   rng[5]
   rng[5].tz_localize('Asia/Shanghai')


Operations between ``Series`` in different time zones will yield UTC
``Series``, aligning the data on the UTC timestamps:

.. ipython:: python

   eastern = ts_utc.tz_convert('US/Eastern')
   berlin = ts_utc.tz_convert('Europe/Berlin')
   result = eastern + berlin
   result
   result.index

To remove timezone from tz-aware ``DatetimeIndex``, use ``tz_localize(None)`` or ``tz_convert(None)``.
``tz_localize(None)`` will remove timezone holding local time representations.
``tz_convert(None)`` will remove timezone after converting to UTC time.

.. ipython:: python

   didx = pd.DatetimeIndex(start='2014-08-01 09:00', freq='H', periods=10, tz='US/Eastern')
   didx
   didx.tz_localize(None)
   didx.tz_convert(None)

   # tz_convert(None) is identical with tz_convert('UTC').tz_localize(None)
   didx.tz_convert('UCT').tz_localize(None)

.. _timeseries.timezone_ambiguous:

Ambiguous Times when Localizing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, localize cannot determine the DST and non-DST hours when there are
duplicates.  This often happens when reading files or database records that simply
duplicate the hours.  Passing ``ambiguous='infer'`` into ``tz_localize`` will
attempt to determine the right offset.  Below the top example will fail as it
contains ambiguous times and the bottom will infer the right offset.

.. ipython:: python

   rng_hourly = pd.DatetimeIndex(['11/06/2011 00:00', '11/06/2011 01:00',
                                  '11/06/2011 01:00', '11/06/2011 02:00',
                                  '11/06/2011 03:00'])

This will fail as there are ambiguous times

.. code-block:: ipython

   In [2]: rng_hourly.tz_localize('US/Eastern')
   AmbiguousTimeError: Cannot infer dst time from Timestamp('2011-11-06 01:00:00'), try using the 'ambiguous' argument

Infer the ambiguous times

.. ipython:: python

   rng_hourly_eastern = rng_hourly.tz_localize('US/Eastern', ambiguous='infer')
   rng_hourly_eastern.tolist()

In addition to 'infer', there are several other arguments supported.  Passing
an array-like of bools or 0s/1s where True represents a DST hour and False a
non-DST hour, allows for distinguishing more than one DST
transition (e.g., if you have multiple records in a database each with their
own DST transition).  Or passing 'NaT' will fill in transition times
with not-a-time values.  These methods are available in the ``DatetimeIndex``
constructor as well as ``tz_localize``.

.. ipython:: python

   rng_hourly_dst = np.array([1, 1, 0, 0, 0])
   rng_hourly.tz_localize('US/Eastern', ambiguous=rng_hourly_dst).tolist()
   rng_hourly.tz_localize('US/Eastern', ambiguous='NaT').tolist()

   didx = pd.DatetimeIndex(start='2014-08-01 09:00', freq='H', periods=10, tz='US/Eastern')
   didx
   didx.tz_localize(None)
   didx.tz_convert(None)

   # tz_convert(None) is identical with tz_convert('UTC').tz_localize(None)
   didx.tz_convert('UCT').tz_localize(None)

.. _timeseries.timezone_series:

TZ Aware Dtypes
~~~~~~~~~~~~~~~

``Series/DatetimeIndex`` with a timezone **naive** value are represented with a dtype of ``datetime64[ns]``.

.. ipython:: python

   s_naive = pd.Series(pd.date_range('20130101',periods=3))
   s_naive

``Series/DatetimeIndex`` with a timezone **aware** value are represented with a dtype of ``datetime64[ns, tz]``.

.. ipython:: python

   s_aware = pd.Series(pd.date_range('20130101',periods=3,tz='US/Eastern'))
   s_aware

Both of these ``Series`` can be manipulated via the ``.dt`` accessor, see :ref:`here <basics.dt_accessors>`.

For example, to localize and convert a naive stamp to timezone aware.

.. ipython:: python

   s_naive.dt.tz_localize('UTC').dt.tz_convert('US/Eastern')


Further more you can ``.astype(...)`` timezone aware (and naive). This operation is effectively a localize AND convert on a naive stamp, and
a convert on an aware stamp.

.. ipython:: python

   # localize and convert a naive timezone
   s_naive.astype('datetime64[ns, US/Eastern]')

   # make an aware tz naive
   s_aware.astype('datetime64[ns]')

   # convert to a new timezone
   s_aware.astype('datetime64[ns, CET]')

.. note::

   Using the ``.values`` accessor on a ``Series``, returns an NumPy array of the data.
   These values are converted to UTC, as NumPy does not currently support timezones (even though it is *printing* in the local timezone!).

   .. ipython:: python

      s_naive.values
      s_aware.values

   Further note that once converted to a NumPy array these would lose the tz tenor.

   .. ipython:: python

      pd.Series(s_aware.values)

   However, these can be easily converted:

   .. ipython:: python

      pd.Series(s_aware.values).dt.tz_localize('UTC').dt.tz_convert('US/Eastern')
