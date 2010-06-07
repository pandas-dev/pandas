.. currentmodule:: pandas

.. _datetools:

************************
Date and frequency tools
************************

In dealing with the economic realities of time series data, we will
frequently seek to:

* generate sequences of fixed-frequency dates

* conform or convert time series to a particular frequency

* compute "relative" dates based on various non-standard time
  increments (e.g. 5 business days before the last business day of the
  year), or "roll" dates forward or backward

pandas provides a relatively compact and self-contained set of tools
for performing the above tasks. The functionality revolves around the
**DateOffset** class and its subclasses, which are objects which can
be added and subtracted from regular Python **datetime** objects.

.. _datetools.offsets:

DateOffset and subclasses
-------------------------

A **DateOffset** instance represents a frequency increment. Different
offset logic via subclasses. The vanilla version of the DateOffset
works identically to an object in the dateutil package,
dateutil.relativedelta, which works like:

::

    >>> d = datetime(2008, 8, 18)
    >>> d + relativedelta(months=4, days=5)
    datetime.datetime(2008, 12, 23, 0, 0)

The basic DateOffset class accepts the same set of arguments that
relativedelta does:

::

    >>> d + DateOffset(months=4, days=5)
    datetime.datetime(2008, 12, 23, 0, 0)

The key characteristics of a DateOffset are:

* it can be added / subtracted to/from a datetime object to obtain a
  shifted date

* it can be multiplied by an integer to multiple the effect as
  desired

Subclasses of DateOffset define specialized logic for various
increments of interest, such a weekday / business day.

::

    class BDay(DateOffset):
	"""DateOffset increments between business days"""
	def apply(self, other):
	    ...

::

    >>> d - 5 * BDay()
    datetime.datetime(2008, 8, 11, 0, 0)

A subclass need only redefine the method 'apply' which describes the
interaction with a datetime object.

::

    >>> d + BMonthEnd(0)
    datetime.datetime(2008, 8, 29, 0, 0)

An offset is also *callable* as functions to avoid those times where
you might otherwise create a lambda, for example:

::

    map(BDay(), dateList) versus
    map(lambda x: x + BDay(), dateList)

.. csv-table::
    :header: "Class name", "Description"
    :widths: 15, 65

    DateOffset, "Generic offset, arguments as **datetime.relativedelta**"
    BDay, "business day (weekday)"
    Week, "one week, optionally anchored on a day of the week"
    MonthEnd, "calendar month end"
    BMonthEnd, "business month end"
    YearEnd, "calendar year end"
    YearBegin, "calendar year begin"
    BYearEnd, "business year end"
    Hour, "one hour"
    Minute, "one minute"
    Second, "one second"

.. _daterange:

Creating date ranges (DateRange)
--------------------------------

The DateRange class utilizes these offsets (and any ones that we might
add) to generate lists of dates for general purposes:

::

    In [327]: DateRange(fromDate=d, nPeriods=10, offset=daterange.bmonthEnd)
    Out[327]:
    [datetime.datetime(2008, 8, 29, 0, 0),
     datetime.datetime(2008, 9, 30, 0, 0),
     datetime.datetime(2008, 10, 31, 0, 0),
     datetime.datetime(2008, 11, 28, 0, 0),
     datetime.datetime(2008, 12, 31, 0, 0),
     datetime.datetime(2009, 1, 30, 0, 0),
     datetime.datetime(2009, 2, 27, 0, 0),
     datetime.datetime(2009, 3, 31, 0, 0),
     datetime.datetime(2009, 4, 30, 0, 0),
     datetime.datetime(2009, 5, 29, 0, 0)]


One can also specify a 'toDate', the endpoints are included if and
only if they lie on the offset (so here, d was today, and is not
included)

There's a lot of general usefulness in these date ranges and offsets,
among others:

1) you can reindex a TimeSeries or DataFrame with them

2) you can use them to do groupby statistics (use a 0-parameter offset
to 'roll' all your dates to the next offset, e.g. EOM date)

3) you can use them to shift the index of your TimeSeries or DataFrame
by some amount:

::

    In [334]: ts
    Out[334]:
    2007-09-28 00:00:00     0.01025
    2007-10-31 00:00:00     0.0059684
    2007-11-30 00:00:00     -0.063186
    2007-12-31 00:00:00     -0.037049
    2008-01-31 00:00:00     -0.10556
    2008-02-29 00:00:00     -0.020221
    2008-03-31 00:00:00     -0.046154
    2008-04-30 00:00:00     0.095594
    2008-05-30 00:00:00     0.031957
    2008-06-30 00:00:00     -0.058455

    In [335]: ts.shift(5, offset=daterange.bday)
    Out[335]:
    2007-10-05 00:00:00     0.01025
    2007-11-07 00:00:00     0.0059684
    2007-12-07 00:00:00     -0.063186
    2008-01-07 00:00:00     -0.037049
    2008-02-07 00:00:00     -0.10556
    2008-03-07 00:00:00     -0.020221
    2008-04-07 00:00:00     -0.046154
    2008-05-07 00:00:00     0.095594
    2008-06-06 00:00:00     0.031957
    2008-07-07 00:00:00     -0.058455

Tips for custom frequencies
---------------------------

.. _datetools.timerules:

Time rules
----------

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

