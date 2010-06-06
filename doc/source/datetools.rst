.. currentmodule:: pandas

.. _datetools:

************************
Date and frequency tools
************************



.. _datetools.offsets:

DateOffset and subclasses
-------------------------

It is often necessary to do non-standard date logic, adding business
days, rolling to the next month or quarter end, rolling to the next
business month end, etc. One might wish to aggregate more granular
data based on these fixed-frequency periods, or maybe you want to
shift a TimeSeries in some particular way without resorting to a hack.

In searching for the most "pythonic" approach to this problem (and I
don't know if I have found it or not), I chose to create objects
representing the frequency increments in the form of the DateOffset
object specifying the mechanics of the "increment", and one can then
define different offset logic via subclasses. The vanilla version of
the DateOffset works identically to an object in the dateutil package,
dateutil.relativedelta, which works like:

::

    In [319]: d = datetime(2008, 8, 18)

    In [321]: d + relativedelta(months=4, days=5)
    Out[321]: datetime.datetime(2008, 12, 23, 0, 0)


You could write this using DateOffset in nearly the same way (here I
have not imported the DateOffset name in this example, so you could do
"from daterange import DateOffset" to avoid having to write the module
name):

::

	In [323]: d + daterange.DateOffset(months=4, days=5)
	Out[323]: datetime.datetime(2008, 12, 23, 0, 0)


The key characteristics of the DateOffset are:

  - it can be added / subtracted to/from a datetime object to obtain a
    shifted date
  - it can be multiplied by an integer to multiple the effect as
    desired

::
    In [325]: d - 5 * daterange.BDay()
    Out[325]: datetime.datetime(2008, 8, 11, 0, 0)

  - in the subclass one only need redefine the method 'apply' which
    describes the interaction with a datetime object.

::

    class BMonthEnd(DateOffset):
	"""DateOffset increments between business EOM dates"""
	def apply(self, other):
	    ...

Turns out each of these 'apply' methods for various standard offsets
are pretty easy to write (the longest is 10 lines), and I've written a
few so far:

business day (BDay)
business month end aka EOM (BMonthEnd)
month end (MonthEnd), quarter end (QuarterEnd), year end/begin (YearEnd / YearBegin)

Each of these when initialized with the parameter 0 (e.g. BMonthEnd(0)
) have the effect of rolling forward to the nearest date lying on a
particular offset:

::

    In [324]: d + daterange.BMonthEnd(0)
    Out[324]: datetime.datetime(2008, 8, 29, 0, 0)


For convenience, the daterange module has shorthand instances of these
classes available for a number of reasons:

bday, monthEnd, bmonthEnd, yearEnd, etc.

::

    In [326]: d - 5 * daterange.bday
    Out[326]: datetime.datetime(2008, 8, 11, 0, 0)


These offsets are also *callable* as functions to avoid those times
where you might otherwise create a lambda, for example:

    map(daterange.bday, dateList) versus
    map(lambda x: x + daterange.bday, dateList)



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

.. csv-table::
    :header: "Class name", "Description"
    :widths: 15, 65

    DateOffset, ""
    BDay, ""
    Week, ""
    MonthEnd, ""
    BMonthEnd, ""
    YearEnd, ""
    YearBegin, ""
    BYearEnd, ""
    Hour, ""
    Minute, ""
    Second, ""

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

Defining custom frequencies
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

