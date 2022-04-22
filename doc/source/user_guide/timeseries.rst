.. _timeseries:

{{ header }}

********************************
Time series / date functionality
********************************

pandas contains extensive capabilities and features for working with time series data for all domains.
Using the NumPy ``datetime64`` and ``timedelta64`` dtypes, pandas has consolidated a large number of
features from other Python libraries like ``scikits.timeseries`` as well as created
a tremendous amount of new functionality for manipulating time series data.

For example, pandas supports:

Parsing time series information from various sources and formats

.. ipython:: python

   import datetime

   dti = pd.to_datetime(
       ["1/1/2018", np.datetime64("2018-01-01"), datetime.datetime(2018, 1, 1)]
   )
   dti

Generate sequences of fixed-frequency dates and time spans

.. ipython:: python

   dti = pd.date_range("2018-01-01", periods=3, freq="H")
   dti

Manipulating and converting date times with timezone information

.. ipython:: python

   dti = dti.tz_localize("UTC")
   dti
   dti.tz_convert("US/Pacific")

Resampling or converting a time series to a particular frequency

.. ipython:: python

   idx = pd.date_range("2018-01-01", periods=5, freq="H")
   ts = pd.Series(range(len(idx)), index=idx)
   ts
   ts.resample("2H").mean()

Performing date and time arithmetic with absolute or relative time increments

.. ipython:: python

    friday = pd.Timestamp("2018-01-05")
    friday.day_name()
    # Add 1 day
    saturday = friday + pd.Timedelta("1 day")
    saturday.day_name()
    # Add 1 business day (Friday --> Monday)
    monday = friday + pd.offsets.BDay()
    monday.day_name()

pandas provides a relatively compact and self-contained set of tools for
performing the above tasks and more.


.. _timeseries.overview:

Overview
--------

pandas captures 4 general time related concepts:

#. Date times: A specific date and time with timezone support. Similar to ``datetime.datetime`` from the standard library.
#. Time deltas: An absolute time duration. Similar to ``datetime.timedelta`` from the standard library.
#. Time spans: A span of time defined by a point in time and its associated frequency.
#. Date offsets: A relative time duration that respects calendar arithmetic. Similar to ``dateutil.relativedelta.relativedelta`` from the ``dateutil`` package.

=====================   =================  ===================   ============================================  ========================================
Concept                 Scalar Class       Array Class           pandas Data Type                              Primary Creation Method
=====================   =================  ===================   ============================================  ========================================
Date times              ``Timestamp``      ``DatetimeIndex``     ``datetime64[ns]`` or ``datetime64[ns, tz]``  ``to_datetime`` or ``date_range``
Time deltas             ``Timedelta``      ``TimedeltaIndex``    ``timedelta64[ns]``                           ``to_timedelta`` or ``timedelta_range``
Time spans              ``Period``         ``PeriodIndex``       ``period[freq]``                              ``Period`` or ``period_range``
Date offsets            ``DateOffset``     ``None``              ``None``                                      ``DateOffset``
=====================   =================  ===================   ============================================  ========================================

For time series data, it's conventional to represent the time component in the index of a :class:`Series` or :class:`DataFrame`
so manipulations can be performed with respect to the time element.

.. ipython:: python

   pd.Series(range(3), index=pd.date_range("2000", freq="D", periods=3))

However, :class:`Series` and :class:`DataFrame` can directly also support the time component as data itself.

.. ipython:: python

   pd.Series(pd.date_range("2000", freq="D", periods=3))

:class:`Series` and :class:`DataFrame` have extended data type support and functionality for ``datetime``, ``timedelta``
and ``Period`` data when passed into those constructors. ``DateOffset``
data however will be stored as ``object`` data.

.. ipython:: python

   pd.Series(pd.period_range("1/1/2011", freq="M", periods=3))
   pd.Series([pd.DateOffset(1), pd.DateOffset(2)])
   pd.Series(pd.date_range("1/1/2011", freq="M", periods=3))

Lastly, pandas represents null date times, time deltas, and time spans as ``NaT`` which
is useful for representing missing or null date like values and behaves similar
as ``np.nan`` does for float data.

.. ipython:: python

   pd.Timestamp(pd.NaT)
   pd.Timedelta(pd.NaT)
   pd.Period(pd.NaT)
   # Equality acts as np.nan would
   pd.NaT == pd.NaT

.. _timeseries.representation:

Timestamps vs. time spans
-------------------------

Timestamped data is the most basic type of time series data that associates
values with points in time. For pandas objects it means using the points in
time.

.. ipython:: python

   pd.Timestamp(datetime.datetime(2012, 5, 1))
   pd.Timestamp("2012-05-01")
   pd.Timestamp(2012, 5, 1)

However, in many cases it is more natural to associate things like change
variables with a time span instead. The span represented by ``Period`` can be
specified explicitly, or inferred from datetime string format.

For example:

.. ipython:: python

   pd.Period("2011-01")

   pd.Period("2012-05", freq="D")

:class:`Timestamp` and :class:`Period` can serve as an index. Lists of
``Timestamp`` and ``Period`` are automatically coerced to :class:`DatetimeIndex`
and :class:`PeriodIndex` respectively.

.. ipython:: python

   dates = [
       pd.Timestamp("2012-05-01"),
       pd.Timestamp("2012-05-02"),
       pd.Timestamp("2012-05-03"),
   ]
   ts = pd.Series(np.random.randn(3), dates)

   type(ts.index)
   ts.index

   ts

   periods = [pd.Period("2012-01"), pd.Period("2012-02"), pd.Period("2012-03")]

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

Converting to timestamps
------------------------

To convert a :class:`Series` or list-like object of date-like objects e.g. strings,
epochs, or a mixture, you can use the ``to_datetime`` function. When passed
a ``Series``, this returns a ``Series`` (with the same index), while a list-like
is converted to a ``DatetimeIndex``:

.. ipython:: python

    pd.to_datetime(pd.Series(["Jul 31, 2009", "2010-01-10", None]))

    pd.to_datetime(["2005/11/23", "2010.12.31"])

If you use dates which start with the day first (i.e. European style),
you can pass the ``dayfirst`` flag:

.. ipython:: python
   :okwarning:

    pd.to_datetime(["04-01-2012 10:00"], dayfirst=True)

    pd.to_datetime(["14-01-2012", "01-14-2012"], dayfirst=True)

.. warning::

   You see in the above example that ``dayfirst`` isn't strict. If a date
   can't be parsed with the day being first it will be parsed as if
   ``dayfirst`` were False, and in the case of parsing delimited date strings
   (e.g. ``31-12-2012``) then a warning will also be raised.

If you pass a single string to ``to_datetime``, it returns a single ``Timestamp``.
``Timestamp`` can also accept string input, but it doesn't accept string parsing
options like ``dayfirst`` or ``format``, so use ``to_datetime`` if these are required.

.. ipython:: python

    pd.to_datetime("2010/11/12")

    pd.Timestamp("2010/11/12")

You can also use the ``DatetimeIndex`` constructor directly:

.. ipython:: python

    pd.DatetimeIndex(["2018-01-01", "2018-01-03", "2018-01-05"])

The string 'infer' can be passed in order to set the frequency of the index as the
inferred frequency upon creation:

.. ipython:: python

    pd.DatetimeIndex(["2018-01-01", "2018-01-03", "2018-01-05"], freq="infer")

.. _timeseries.converting.format:

Providing a format argument
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the required datetime string, a ``format`` argument can be passed to ensure specific parsing.
This could also potentially speed up the conversion considerably.

.. ipython:: python

    pd.to_datetime("2010/11/12", format="%Y/%m/%d")

    pd.to_datetime("12-11-2010 00:00", format="%d-%m-%Y %H:%M")

For more information on the choices available when specifying the ``format``
option, see the Python `datetime documentation`_.

.. _datetime documentation: https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior

Assembling datetime from multiple DataFrame columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also pass a ``DataFrame`` of integer or string columns to assemble into a ``Series`` of ``Timestamps``.

.. ipython:: python

   df = pd.DataFrame(
       {"year": [2015, 2016], "month": [2, 3], "day": [4, 5], "hour": [2, 3]}
   )
   pd.to_datetime(df)


You can pass only the columns that you need to assemble.

.. ipython:: python

   pd.to_datetime(df[["year", "month", "day"]])

``pd.to_datetime`` looks for standard designations of the datetime component in the column names, including:

* required: ``year``, ``month``, ``day``
* optional: ``hour``, ``minute``, ``second``, ``millisecond``, ``microsecond``, ``nanosecond``

Invalid data
~~~~~~~~~~~~

The default behavior, ``errors='raise'``, is to raise when unparsable:

.. code-block:: ipython

    In [2]: pd.to_datetime(['2009/07/31', 'asd'], errors='raise')
    ValueError: Unknown string format

Pass ``errors='ignore'`` to return the original input when unparsable:

.. ipython:: python

   pd.to_datetime(["2009/07/31", "asd"], errors="ignore")

Pass ``errors='coerce'`` to convert unparsable data to ``NaT`` (not a time):

.. ipython:: python

   pd.to_datetime(["2009/07/31", "asd"], errors="coerce")


.. _timeseries.converting.epoch:

Epoch timestamps
~~~~~~~~~~~~~~~~

pandas supports converting integer or float epoch times to ``Timestamp`` and
``DatetimeIndex``. The default unit is nanoseconds, since that is how ``Timestamp``
objects are stored internally. However, epochs are often stored in another ``unit``
which can be specified. These are computed from the starting point specified by the
``origin`` parameter.

.. ipython:: python

   pd.to_datetime(
       [1349720105, 1349806505, 1349892905, 1349979305, 1350065705], unit="s"
   )

   pd.to_datetime(
       [1349720105100, 1349720105200, 1349720105300, 1349720105400, 1349720105500],
       unit="ms",
   )

.. note::

   The ``unit`` parameter does not use the same strings as the ``format`` parameter
   that was discussed :ref:`above<timeseries.converting.format>`). The
   available units are listed on the documentation for :func:`pandas.to_datetime`.

.. versionchanged:: 1.0.0

Constructing a :class:`Timestamp` or :class:`DatetimeIndex` with an epoch timestamp
with the ``tz`` argument specified will raise a ValueError. If you have
epochs in wall time in another timezone, you can read the epochs
as timezone-naive timestamps and then localize to the appropriate timezone:

.. ipython:: python

   pd.Timestamp(1262347200000000000).tz_localize("US/Pacific")
   pd.DatetimeIndex([1262347200000000000]).tz_localize("US/Pacific")

.. note::

   Epoch times will be rounded to the nearest nanosecond.

.. warning::

   Conversion of float epoch times can lead to inaccurate and unexpected results.
   :ref:`Python floats <python:tut-fp-issues>` have about 15 digits precision in
   decimal. Rounding during conversion from float to high precision ``Timestamp`` is
   unavoidable. The only way to achieve exact precision is to use a fixed-width
   types (e.g. an int64).

   .. ipython:: python

      pd.to_datetime([1490195805.433, 1490195805.433502912], unit="s")
      pd.to_datetime(1490195805433502912, unit="ns")

.. seealso::

   :ref:`timeseries.origin`

.. _timeseries.converting.epoch_inverse:

From timestamps to epoch
~~~~~~~~~~~~~~~~~~~~~~~~

To invert the operation from above, namely, to convert from a ``Timestamp`` to a 'unix' epoch:

.. ipython:: python

   stamps = pd.date_range("2012-10-08 18:15:05", periods=4, freq="D")
   stamps

We subtract the epoch (midnight at January 1, 1970 UTC) and then floor divide by the
"unit" (1 second).

.. ipython:: python

   (stamps - pd.Timestamp("1970-01-01")) // pd.Timedelta("1s")

.. _timeseries.origin:

Using the ``origin`` Parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Using the ``origin`` parameter, one can specify an alternative starting point for creation
of a ``DatetimeIndex``. For example, to use 1960-01-01 as the starting date:

.. ipython:: python

   pd.to_datetime([1, 2, 3], unit="D", origin=pd.Timestamp("1960-01-01"))

The default is set at ``origin='unix'``, which defaults to ``1970-01-01 00:00:00``.
Commonly called 'unix epoch' or POSIX time.

.. ipython:: python

   pd.to_datetime([1, 2, 3], unit="D")

.. _timeseries.daterange:

Generating ranges of timestamps
-------------------------------

To generate an index with timestamps, you can use either the ``DatetimeIndex`` or
``Index`` constructor and pass in a list of datetime objects:

.. ipython:: python

   dates = [
       datetime.datetime(2012, 5, 1),
       datetime.datetime(2012, 5, 2),
       datetime.datetime(2012, 5, 3),
   ]

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
**calendar day** while the default for ``bdate_range`` is a **business day**:

.. ipython:: python

   start = datetime.datetime(2011, 1, 1)
   end = datetime.datetime(2012, 1, 1)

   index = pd.date_range(start, end)
   index

   index = pd.bdate_range(start, end)
   index

Convenience functions like ``date_range`` and ``bdate_range`` can utilize a
variety of :ref:`frequency aliases <timeseries.offset_aliases>`:

.. ipython:: python

   pd.date_range(start, periods=1000, freq="M")

   pd.bdate_range(start, periods=250, freq="BQS")

``date_range`` and ``bdate_range`` make it easy to generate a range of dates
using various combinations of parameters like ``start``, ``end``, ``periods``,
and ``freq``. The start and end dates are strictly inclusive, so dates outside
of those specified will not be generated:

.. ipython:: python

   pd.date_range(start, end, freq="BM")

   pd.date_range(start, end, freq="W")

   pd.bdate_range(end=end, periods=20)

   pd.bdate_range(start=start, periods=20)

Specifying ``start``, ``end``, and ``periods`` will generate a range of evenly spaced
dates from ``start`` to ``end`` inclusively, with ``periods`` number of elements in the
resulting ``DatetimeIndex``:

.. ipython:: python

   pd.date_range("2018-01-01", "2018-01-05", periods=5)

   pd.date_range("2018-01-01", "2018-01-05", periods=10)

.. _timeseries.custom-freq-ranges:

Custom frequency ranges
~~~~~~~~~~~~~~~~~~~~~~~

``bdate_range`` can also generate a range of custom frequency dates by using
the ``weekmask`` and ``holidays`` parameters.  These parameters will only be
used if a custom frequency string is passed.

.. ipython:: python

   weekmask = "Mon Wed Fri"

   holidays = [datetime.datetime(2011, 1, 5), datetime.datetime(2011, 3, 14)]

   pd.bdate_range(start, end, freq="C", weekmask=weekmask, holidays=holidays)

   pd.bdate_range(start, end, freq="CBMS", weekmask=weekmask)

.. seealso::

   :ref:`timeseries.custombusinessdays`

.. _timeseries.timestamp-limits:

Timestamp limitations
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
* Fast shifting using the ``shift`` method on pandas objects.
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

   rng = pd.date_range(start, end, freq="BM")
   ts = pd.Series(np.random.randn(len(rng)), index=rng)
   ts.index
   ts[:5].index
   ts[::2].index

.. _timeseries.partialindexing:

Partial string indexing
~~~~~~~~~~~~~~~~~~~~~~~

Dates and strings that parse to timestamps can be passed as indexing parameters:

.. ipython:: python

   ts["1/31/2011"]

   ts[datetime.datetime(2011, 12, 25):]

   ts["10/31/2011":"12/31/2011"]

To provide convenience for accessing longer time series, you can also pass in
the year or year and month as strings:

.. ipython:: python

   ts["2011"]

   ts["2011-6"]

This type of slicing will work on a ``DataFrame`` with a ``DatetimeIndex`` as well. Since the
partial string selection is a form of label slicing, the endpoints **will be** included. This
would include matching times on an included date:

.. warning::

   Indexing ``DataFrame`` rows with a *single* string with getitem (e.g. ``frame[dtstring]``)
   is deprecated starting with pandas 1.2.0 (given the ambiguity whether it is indexing
   the rows or selecting a column) and will be removed in a future version. The equivalent
   with ``.loc`` (e.g. ``frame.loc[dtstring]``) is still supported.

.. ipython:: python

   dft = pd.DataFrame(
       np.random.randn(100000, 1),
       columns=["A"],
       index=pd.date_range("20130101", periods=100000, freq="T"),
   )
   dft
   dft.loc["2013"]

This starts on the very first time in the month, and includes the last date and
time for the month:

.. ipython:: python

   dft["2013-1":"2013-2"]

This specifies a stop time **that includes all of the times on the last day**:

.. ipython:: python

   dft["2013-1":"2013-2-28"]

This specifies an **exact** stop time (and is not the same as the above):

.. ipython:: python

   dft["2013-1":"2013-2-28 00:00:00"]

We are stopping on the included end-point as it is part of the index:

.. ipython:: python

   dft["2013-1-15":"2013-1-15 12:30:00"]

``DatetimeIndex`` partial string indexing also works on a ``DataFrame`` with a ``MultiIndex``:

.. ipython:: python

   dft2 = pd.DataFrame(
       np.random.randn(20, 1),
       columns=["A"],
       index=pd.MultiIndex.from_product(
           [pd.date_range("20130101", periods=10, freq="12H"), ["a", "b"]]
       ),
   )
   dft2
   dft2.loc["2013-01-05"]
   idx = pd.IndexSlice
   dft2 = dft2.swaplevel(0, 1).sort_index()
   dft2.loc[idx[:, "2013-01-05"], :]

.. versionadded:: 0.25.0

Slicing with string indexing also honors UTC offset.

.. ipython:: python

    df = pd.DataFrame([0], index=pd.DatetimeIndex(["2019-01-01"], tz="US/Pacific"))
    df
    df["2019-01-01 12:00:00+04:00":"2019-01-01 13:00:00+04:00"]

.. _timeseries.slice_vs_exact_match:

Slice vs. exact match
~~~~~~~~~~~~~~~~~~~~~

The same string used as an indexing parameter can be treated either as a slice or as an exact match depending on the resolution of the index. If the string is less accurate than the index, it will be treated as a slice, otherwise as an exact match.

Consider a ``Series`` object with a minute resolution index:

.. ipython:: python

    series_minute = pd.Series(
        [1, 2, 3],
        pd.DatetimeIndex(
            ["2011-12-31 23:59:00", "2012-01-01 00:00:00", "2012-01-01 00:02:00"]
        ),
    )
    series_minute.index.resolution

A timestamp string less accurate than a minute gives a ``Series`` object.

.. ipython:: python

    series_minute["2011-12-31 23"]

A timestamp string with minute resolution (or more accurate), gives a scalar instead, i.e. it is not casted to a slice.

.. ipython:: python

    series_minute["2011-12-31 23:59"]
    series_minute["2011-12-31 23:59:00"]

If index resolution is second, then the minute-accurate timestamp gives a
``Series``.

.. ipython:: python

    series_second = pd.Series(
        [1, 2, 3],
        pd.DatetimeIndex(
            ["2011-12-31 23:59:59", "2012-01-01 00:00:00", "2012-01-01 00:00:01"]
        ),
    )
    series_second.index.resolution
    series_second["2011-12-31 23:59"]

If the timestamp string is treated as a slice, it can be used to index ``DataFrame`` with ``.loc[]`` as well.

.. ipython:: python

    dft_minute = pd.DataFrame(
        {"a": [1, 2, 3], "b": [4, 5, 6]}, index=series_minute.index
    )
    dft_minute.loc["2011-12-31 23"]


.. warning::

   However, if the string is treated as an exact match, the selection in ``DataFrame``'s ``[]`` will be column-wise and not row-wise, see :ref:`Indexing Basics <indexing.basics>`. For example ``dft_minute['2011-12-31 23:59']`` will raise ``KeyError`` as ``'2012-12-31 23:59'`` has the same resolution as the index and there is no column with such name:

   To *always* have unambiguous selection, whether the row is treated as a slice or a single selection, use ``.loc``.

   .. ipython:: python

      dft_minute.loc["2011-12-31 23:59"]

Note also that ``DatetimeIndex`` resolution cannot be less precise than day.

.. ipython:: python

    series_monthly = pd.Series(
        [1, 2, 3], pd.DatetimeIndex(["2011-12", "2012-01", "2012-02"])
    )
    series_monthly.index.resolution
    series_monthly["2011-12"]  # returns Series


Exact indexing
~~~~~~~~~~~~~~

As discussed in previous section, indexing a ``DatetimeIndex`` with a partial string depends on the "accuracy" of the period, in other words how specific the interval is in relation to the resolution of the index. In contrast, indexing with ``Timestamp`` or ``datetime`` objects is exact, because the objects have exact meaning. These also follow the semantics of *including both endpoints*.

These ``Timestamp`` and ``datetime`` objects have exact ``hours, minutes,`` and ``seconds``, even though they were not explicitly specified (they are ``0``).

.. ipython:: python

   dft[datetime.datetime(2013, 1, 1): datetime.datetime(2013, 2, 28)]

With no defaults.

.. ipython:: python

   dft[
       datetime.datetime(2013, 1, 1, 10, 12, 0): datetime.datetime(
           2013, 2, 28, 10, 12, 0
       )
   ]

Truncating & fancy indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :meth:`~DataFrame.truncate` convenience function is provided that is similar
to slicing. Note that ``truncate`` assumes a 0 value for any unspecified date
component in a ``DatetimeIndex`` in contrast to slicing which returns any
partially matching dates:

.. ipython:: python

   rng2 = pd.date_range("2011-01-01", "2012-01-01", freq="W")
   ts2 = pd.Series(np.random.randn(len(rng2)), index=rng2)

   ts2.truncate(before="2011-11", after="2011-12")
   ts2["2011-11":"2011-12"]

Even complicated fancy indexing that breaks the ``DatetimeIndex`` frequency
regularity will result in a ``DatetimeIndex``, although frequency is lost:

.. ipython:: python

   ts2[[0, 2, 6]].index

.. _timeseries.components:

Time/date components
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
    day_of_year,"The ordinal day of year"
    weekofyear,"The week ordinal of the year"
    week,"The week ordinal of the year"
    dayofweek,"The number of the day of the week with Monday=0, Sunday=6"
    day_of_week,"The number of the day of the week with Monday=0, Sunday=6"
    weekday,"The number of the day of the week with Monday=0, Sunday=6"
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

.. versionadded:: 1.1.0

You may obtain the year, week and day components of the ISO year from the ISO 8601 standard:

.. ipython:: python

   idx = pd.date_range(start="2019-12-29", freq="D", periods=4)
   idx.isocalendar()
   idx.to_series().dt.isocalendar()

.. _timeseries.offsets:

DateOffset objects
------------------

In the preceding examples, frequency strings (e.g. ``'D'``) were used to specify
a frequency that defined:

* how the date times in :class:`DatetimeIndex` were spaced when using :meth:`date_range`
* the frequency of a :class:`Period` or :class:`PeriodIndex`

These frequency strings map to a :class:`DateOffset` object and its subclasses. A :class:`DateOffset`
is similar to a :class:`Timedelta` that represents a duration of time but follows specific calendar duration rules.
For example, a :class:`Timedelta` day will always increment ``datetimes`` by 24 hours, while a :class:`DateOffset` day
will increment ``datetimes`` to the same time the next day whether a day represents 23, 24 or 25 hours due to daylight
savings time. However, all :class:`DateOffset` subclasses that are an hour or smaller
(``Hour``, ``Minute``, ``Second``, ``Milli``, ``Micro``, ``Nano``) behave like
:class:`Timedelta` and respect absolute time.

The basic :class:`DateOffset` acts similar to ``dateutil.relativedelta`` (`relativedelta documentation`_)
that shifts a date time by the corresponding calendar duration specified. The
arithmetic operator (``+``) can be used to perform the shift.

.. ipython:: python

   # This particular day contains a day light savings time transition
   ts = pd.Timestamp("2016-10-30 00:00:00", tz="Europe/Helsinki")
   # Respects absolute time
   ts + pd.Timedelta(days=1)
   # Respects calendar time
   ts + pd.DateOffset(days=1)
   friday = pd.Timestamp("2018-01-05")
   friday.day_name()
   # Add 2 business days (Friday --> Tuesday)
   two_business_days = 2 * pd.offsets.BDay()
   friday + two_business_days
   (friday + two_business_days).day_name()


Most ``DateOffsets`` have associated frequencies strings, or offset aliases, that can be passed
into ``freq`` keyword arguments. The available date offsets and associated frequency strings can be found below:

.. csv-table::
    :header: "Date Offset", "Frequency String", "Description"
    :widths: 15, 15, 65

    :class:`~pandas.tseries.offsets.DateOffset`, None, "Generic offset class, defaults to absolute 24 hours"
    :class:`~pandas.tseries.offsets.BDay` or :class:`~pandas.tseries.offsets.BusinessDay`, ``'B'``,"business day (weekday)"
    :class:`~pandas.tseries.offsets.CDay` or :class:`~pandas.tseries.offsets.CustomBusinessDay`, ``'C'``, "custom business day"
    :class:`~pandas.tseries.offsets.Week`, ``'W'``, "one week, optionally anchored on a day of the week"
    :class:`~pandas.tseries.offsets.WeekOfMonth`, ``'WOM'``, "the x-th day of the y-th week of each month"
    :class:`~pandas.tseries.offsets.LastWeekOfMonth`, ``'LWOM'``, "the x-th day of the last week of each month"
    :class:`~pandas.tseries.offsets.MonthEnd`, ``'M'``, "calendar month end"
    :class:`~pandas.tseries.offsets.MonthBegin`, ``'MS'``, "calendar month begin"
    :class:`~pandas.tseries.offsets.BMonthEnd` or :class:`~pandas.tseries.offsets.BusinessMonthEnd`, ``'BM'``, "business month end"
    :class:`~pandas.tseries.offsets.BMonthBegin` or :class:`~pandas.tseries.offsets.BusinessMonthBegin`, ``'BMS'``, "business month begin"
    :class:`~pandas.tseries.offsets.CBMonthEnd` or :class:`~pandas.tseries.offsets.CustomBusinessMonthEnd`, ``'CBM'``, "custom business month end"
    :class:`~pandas.tseries.offsets.CBMonthBegin` or :class:`~pandas.tseries.offsets.CustomBusinessMonthBegin`, ``'CBMS'``, "custom business month begin"
    :class:`~pandas.tseries.offsets.SemiMonthEnd`, ``'SM'``, "15th (or other day_of_month) and calendar month end"
    :class:`~pandas.tseries.offsets.SemiMonthBegin`, ``'SMS'``, "15th (or other day_of_month) and calendar month begin"
    :class:`~pandas.tseries.offsets.QuarterEnd`, ``'Q'``, "calendar quarter end"
    :class:`~pandas.tseries.offsets.QuarterBegin`, ``'QS'``, "calendar quarter begin"
    :class:`~pandas.tseries.offsets.BQuarterEnd`, ``'BQ``, "business quarter end"
    :class:`~pandas.tseries.offsets.BQuarterBegin`, ``'BQS'``, "business quarter begin"
    :class:`~pandas.tseries.offsets.FY5253Quarter`, ``'REQ'``, "retail (aka 52-53 week) quarter"
    :class:`~pandas.tseries.offsets.YearEnd`, ``'A'``, "calendar year end"
    :class:`~pandas.tseries.offsets.YearBegin`, ``'AS'`` or ``'BYS'``,"calendar year begin"
    :class:`~pandas.tseries.offsets.BYearEnd`, ``'BA'``, "business year end"
    :class:`~pandas.tseries.offsets.BYearBegin`, ``'BAS'``, "business year begin"
    :class:`~pandas.tseries.offsets.FY5253`, ``'RE'``, "retail (aka 52-53 week) year"
    :class:`~pandas.tseries.offsets.Easter`, None, "Easter holiday"
    :class:`~pandas.tseries.offsets.BusinessHour`, ``'BH'``, "business hour"
    :class:`~pandas.tseries.offsets.CustomBusinessHour`, ``'CBH'``, "custom business hour"
    :class:`~pandas.tseries.offsets.Day`, ``'D'``, "one absolute day"
    :class:`~pandas.tseries.offsets.Hour`, ``'H'``, "one hour"
    :class:`~pandas.tseries.offsets.Minute`, ``'T'`` or ``'min'``,"one minute"
    :class:`~pandas.tseries.offsets.Second`, ``'S'``, "one second"
    :class:`~pandas.tseries.offsets.Milli`, ``'L'`` or ``'ms'``, "one millisecond"
    :class:`~pandas.tseries.offsets.Micro`, ``'U'`` or ``'us'``, "one microsecond"
    :class:`~pandas.tseries.offsets.Nano`, ``'N'``, "one nanosecond"

``DateOffsets`` additionally have :meth:`rollforward` and :meth:`rollback`
methods for moving a date forward or backward respectively to a valid offset
date relative to the offset. For example, business offsets will roll dates
that land on the weekends (Saturday and Sunday) forward to Monday since
business offsets operate on the weekdays.

.. ipython:: python

   ts = pd.Timestamp("2018-01-06 00:00:00")
   ts.day_name()
   # BusinessHour's valid offset dates are Monday through Friday
   offset = pd.offsets.BusinessHour(start="09:00")
   # Bring the date to the closest offset date (Monday)
   offset.rollforward(ts)
   # Date is brought to the closest offset date first and then the hour is added
   ts + offset

These operations preserve time (hour, minute, etc) information by default.
To reset time to midnight, use :meth:`normalize` before or after applying
the operation (depending on whether you want the time information included
in the operation).

.. ipython:: python

   ts = pd.Timestamp("2014-01-01 09:00")
   day = pd.offsets.Day()
   day + ts
   (day + ts).normalize()

   ts = pd.Timestamp("2014-01-01 22:00")
   hour = pd.offsets.Hour()
   hour + ts
   (hour + ts).normalize()
   (hour + pd.Timestamp("2014-01-01 23:30")).normalize()

.. _relativedelta documentation: https://dateutil.readthedocs.io/en/stable/relativedelta.html


Parametric offsets
~~~~~~~~~~~~~~~~~~

Some of the offsets can be "parameterized" when created to result in different
behaviors. For example, the ``Week`` offset for generating weekly data accepts a
``weekday`` parameter which results in the generated dates always lying on a
particular day of the week:

.. ipython:: python

   d = datetime.datetime(2008, 8, 18, 9, 0)
   d
   d + pd.offsets.Week()
   d + pd.offsets.Week(weekday=4)
   (d + pd.offsets.Week(weekday=4)).weekday()

   d - pd.offsets.Week()

The ``normalize`` option will be effective for addition and subtraction.

.. ipython:: python

   d + pd.offsets.Week(normalize=True)
   d - pd.offsets.Week(normalize=True)


Another example is parameterizing ``YearEnd`` with the specific ending month:

.. ipython:: python

   d + pd.offsets.YearEnd()
   d + pd.offsets.YearEnd(month=6)


.. _timeseries.offsetseries:

Using offsets with ``Series`` / ``DatetimeIndex``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Offsets can be used with either a ``Series`` or ``DatetimeIndex`` to
apply the offset to each element.

.. ipython:: python

   rng = pd.date_range("2012-01-01", "2012-01-03")
   s = pd.Series(rng)
   rng
   rng + pd.DateOffset(months=2)
   s + pd.DateOffset(months=2)
   s - pd.DateOffset(months=2)

If the offset class maps directly to a ``Timedelta`` (``Day``, ``Hour``,
``Minute``, ``Second``, ``Micro``, ``Milli``, ``Nano``) it can be
used exactly like a ``Timedelta`` - see the
:ref:`Timedelta section<timedeltas.operations>` for more examples.

.. ipython:: python

   s - pd.offsets.Day(2)
   td = s - pd.Series(pd.date_range("2011-12-29", "2011-12-31"))
   td
   td + pd.offsets.Minute(15)

Note that some offsets (such as ``BQuarterEnd``) do not have a
vectorized implementation.  They can still be used but may
calculate significantly slower and will show a ``PerformanceWarning``

.. ipython:: python
   :okwarning:

   rng + pd.offsets.BQuarterEnd()


.. _timeseries.custombusinessdays:

Custom business days
~~~~~~~~~~~~~~~~~~~~

The ``CDay`` or ``CustomBusinessDay`` class provides a parametric
``BusinessDay`` class which can be used to create customized business day
calendars which account for local holidays and local weekend conventions.

As an interesting example, let's look at Egypt where a Friday-Saturday weekend is observed.

.. ipython:: python

    weekmask_egypt = "Sun Mon Tue Wed Thu"

    # They also observe International Workers' Day so let's
    # add that for a couple of years

    holidays = [
        "2012-05-01",
        datetime.datetime(2013, 5, 1),
        np.datetime64("2014-05-01"),
    ]
    bday_egypt = pd.offsets.CustomBusinessDay(
        holidays=holidays,
        weekmask=weekmask_egypt,
    )
    dt = datetime.datetime(2013, 4, 30)
    dt + 2 * bday_egypt

Let's map to the weekday names:

.. ipython:: python

    dts = pd.date_range(dt, periods=5, freq=bday_egypt)

    pd.Series(dts.weekday, dts).map(pd.Series("Mon Tue Wed Thu Fri Sat Sun".split()))

Holiday calendars can be used to provide the list of holidays.  See the
:ref:`holiday calendar<timeseries.holiday>` section for more information.

.. ipython:: python

    from pandas.tseries.holiday import USFederalHolidayCalendar

    bday_us = pd.offsets.CustomBusinessDay(calendar=USFederalHolidayCalendar())

    # Friday before MLK Day
    dt = datetime.datetime(2014, 1, 17)

    # Tuesday after MLK Day (Monday is skipped because it's a holiday)
    dt + bday_us

Monthly offsets that respect a certain holiday calendar can be defined
in the usual way.

.. ipython:: python

    bmth_us = pd.offsets.CustomBusinessMonthBegin(calendar=USFederalHolidayCalendar())

    # Skip new years
    dt = datetime.datetime(2013, 12, 17)
    dt + bmth_us

    # Define date index with custom offset
    pd.date_range(start="20100101", end="20120101", freq=bmth_us)

.. note::

    The frequency string 'C' is used to indicate that a CustomBusinessDay
    DateOffset is used, it is important to note that since CustomBusinessDay is
    a parameterised type, instances of CustomBusinessDay may differ and this is
    not detectable from the 'C' frequency string. The user therefore needs to
    ensure that the 'C' frequency string is used consistently within the user's
    application.

.. _timeseries.businesshour:

Business hour
~~~~~~~~~~~~~

The ``BusinessHour`` class provides a business hour representation on ``BusinessDay``,
allowing to use specific start and end times.

By default, ``BusinessHour`` uses 9:00 - 17:00 as business hours.
Adding ``BusinessHour`` will increment ``Timestamp`` by hourly frequency.
If target ``Timestamp`` is out of business hours, move to the next business hour
then increment it. If the result exceeds the business hours end, the remaining
hours are added to the next business day.

.. ipython:: python

    bh = pd.offsets.BusinessHour()
    bh

    # 2014-08-01 is Friday
    pd.Timestamp("2014-08-01 10:00").weekday()
    pd.Timestamp("2014-08-01 10:00") + bh

    # Below example is the same as: pd.Timestamp('2014-08-01 09:00') + bh
    pd.Timestamp("2014-08-01 08:00") + bh

    # If the results is on the end time, move to the next business day
    pd.Timestamp("2014-08-01 16:00") + bh

    # Remainings are added to the next day
    pd.Timestamp("2014-08-01 16:30") + bh

    # Adding 2 business hours
    pd.Timestamp("2014-08-01 10:00") + pd.offsets.BusinessHour(2)

    # Subtracting 3 business hours
    pd.Timestamp("2014-08-01 10:00") + pd.offsets.BusinessHour(-3)

You can also specify ``start`` and ``end`` time by keywords. The argument must
be a ``str`` with an ``hour:minute`` representation or a ``datetime.time``
instance. Specifying seconds, microseconds and nanoseconds as business hour
results in ``ValueError``.

.. ipython:: python

    bh = pd.offsets.BusinessHour(start="11:00", end=datetime.time(20, 0))
    bh

    pd.Timestamp("2014-08-01 13:00") + bh
    pd.Timestamp("2014-08-01 09:00") + bh
    pd.Timestamp("2014-08-01 18:00") + bh

Passing ``start`` time later than ``end`` represents midnight business hour.
In this case, business hour exceeds midnight and overlap to the next day.
Valid business hours are distinguished by whether it started from valid ``BusinessDay``.

.. ipython:: python

    bh = pd.offsets.BusinessHour(start="17:00", end="09:00")
    bh

    pd.Timestamp("2014-08-01 17:00") + bh
    pd.Timestamp("2014-08-01 23:00") + bh

    # Although 2014-08-02 is Saturday,
    # it is valid because it starts from 08-01 (Friday).
    pd.Timestamp("2014-08-02 04:00") + bh

    # Although 2014-08-04 is Monday,
    # it is out of business hours because it starts from 08-03 (Sunday).
    pd.Timestamp("2014-08-04 04:00") + bh

Applying ``BusinessHour.rollforward`` and ``rollback`` to out of business hours results in
the next business hour start or previous day's end. Different from other offsets, ``BusinessHour.rollforward``
may output different results from ``apply`` by definition.

This is because one day's business hour end is equal to next day's business hour start. For example,
under the default business hours (9:00 - 17:00), there is no gap (0 minutes) between ``2014-08-01 17:00`` and
``2014-08-04 09:00``.

.. ipython:: python

    # This adjusts a Timestamp to business hour edge
    pd.offsets.BusinessHour().rollback(pd.Timestamp("2014-08-02 15:00"))
    pd.offsets.BusinessHour().rollforward(pd.Timestamp("2014-08-02 15:00"))

    # It is the same as BusinessHour() + pd.Timestamp('2014-08-01 17:00').
    # And it is the same as BusinessHour() + pd.Timestamp('2014-08-04 09:00')
    pd.offsets.BusinessHour() + pd.Timestamp("2014-08-02 15:00")

    # BusinessDay results (for reference)
    pd.offsets.BusinessHour().rollforward(pd.Timestamp("2014-08-02"))

    # It is the same as BusinessDay() + pd.Timestamp('2014-08-01')
    # The result is the same as rollworward because BusinessDay never overlap.
    pd.offsets.BusinessHour() + pd.Timestamp("2014-08-02")

``BusinessHour`` regards Saturday and Sunday as holidays. To use arbitrary
holidays, you can use ``CustomBusinessHour`` offset, as explained in the
following subsection.

.. _timeseries.custombusinesshour:

Custom business hour
~~~~~~~~~~~~~~~~~~~~

The ``CustomBusinessHour`` is a mixture of ``BusinessHour`` and ``CustomBusinessDay`` which
allows you to specify arbitrary holidays. ``CustomBusinessHour`` works as the same
as ``BusinessHour`` except that it skips specified custom holidays.

.. ipython:: python

    from pandas.tseries.holiday import USFederalHolidayCalendar

    bhour_us = pd.offsets.CustomBusinessHour(calendar=USFederalHolidayCalendar())
    # Friday before MLK Day
    dt = datetime.datetime(2014, 1, 17, 15)

    dt + bhour_us

    # Tuesday after MLK Day (Monday is skipped because it's a holiday)
    dt + bhour_us * 2

You can use keyword arguments supported by either ``BusinessHour`` and ``CustomBusinessDay``.

.. ipython:: python

    bhour_mon = pd.offsets.CustomBusinessHour(start="10:00", weekmask="Tue Wed Thu Fri")

    # Monday is skipped because it's a holiday, business hour starts from 10:00
    dt + bhour_mon * 2

.. _timeseries.offset_aliases:

Offset aliases
~~~~~~~~~~~~~~

A number of string aliases are given to useful common time series
frequencies. We will refer to these aliases as *offset aliases*.

.. csv-table::
    :header: "Alias", "Description"
    :widths: 15, 100

    "B", "business day frequency"
    "C", "custom business day frequency"
    "D", "calendar day frequency"
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

.. note::

    When using the offset aliases above, it should be noted that functions
    such as :func:`date_range`, :func:`bdate_range`, will only return
    timestamps that are in the interval defined by ``start_date`` and
    ``end_date``. If the ``start_date`` does not correspond to the frequency,
    the returned timestamps will start at the next valid timestamp, same for
    ``end_date``, the returned timestamps will stop at the previous valid
    timestamp.

   For example, for the offset ``MS``, if the ``start_date`` is not the first
   of the month, the returned timestamps will start with the first day of the
   next month. If ``end_date`` is not the first day of a month, the last
   returned timestamp will be the first day of the corresponding month.

   .. ipython:: python

       dates_lst_1 = pd.date_range("2020-01-06", "2020-04-03", freq="MS")
       dates_lst_1

       dates_lst_2 = pd.date_range("2020-01-01", "2020-04-01", freq="MS")
       dates_lst_2

   We can see in the above example :func:`date_range` and
   :func:`bdate_range` will only return the valid timestamps between the
   ``start_date`` and ``end_date``. If these are not valid timestamps for the
   given frequency it will roll to the next value for ``start_date``
   (respectively previous for the ``end_date``)


Combining aliases
~~~~~~~~~~~~~~~~~

As we have seen previously, the alias and the offset instance are fungible in
most functions:

.. ipython:: python

   pd.date_range(start, periods=5, freq="B")

   pd.date_range(start, periods=5, freq=pd.offsets.BDay())

You can combine together day and intraday offsets:

.. ipython:: python

   pd.date_range(start, periods=10, freq="2h20min")

   pd.date_range(start, periods=10, freq="1D10U")

Anchored offsets
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

Anchored offset semantics
~~~~~~~~~~~~~~~~~~~~~~~~~

For those offsets that are anchored to the start or end of specific
frequency (``MonthEnd``, ``MonthBegin``, ``WeekEnd``, etc), the following
rules apply to rolling forward and backwards.

When ``n`` is not 0, if the given date is not on an anchor point, it snapped to the next(previous)
anchor point, and moved ``|n|-1`` additional steps forwards or backwards.

.. ipython:: python

   pd.Timestamp("2014-01-02") + pd.offsets.MonthBegin(n=1)
   pd.Timestamp("2014-01-02") + pd.offsets.MonthEnd(n=1)

   pd.Timestamp("2014-01-02") - pd.offsets.MonthBegin(n=1)
   pd.Timestamp("2014-01-02") - pd.offsets.MonthEnd(n=1)

   pd.Timestamp("2014-01-02") + pd.offsets.MonthBegin(n=4)
   pd.Timestamp("2014-01-02") - pd.offsets.MonthBegin(n=4)

If the given date *is* on an anchor point, it is moved ``|n|`` points forwards
or backwards.

.. ipython:: python

   pd.Timestamp("2014-01-01") + pd.offsets.MonthBegin(n=1)
   pd.Timestamp("2014-01-31") + pd.offsets.MonthEnd(n=1)

   pd.Timestamp("2014-01-01") - pd.offsets.MonthBegin(n=1)
   pd.Timestamp("2014-01-31") - pd.offsets.MonthEnd(n=1)

   pd.Timestamp("2014-01-01") + pd.offsets.MonthBegin(n=4)
   pd.Timestamp("2014-01-31") - pd.offsets.MonthBegin(n=4)

For the case when ``n=0``, the date is not moved if on an anchor point, otherwise
it is rolled forward to the next anchor point.

.. ipython:: python

   pd.Timestamp("2014-01-02") + pd.offsets.MonthBegin(n=0)
   pd.Timestamp("2014-01-02") + pd.offsets.MonthEnd(n=0)

   pd.Timestamp("2014-01-01") + pd.offsets.MonthBegin(n=0)
   pd.Timestamp("2014-01-31") + pd.offsets.MonthEnd(n=0)

.. _timeseries.holiday:

Holidays / holiday calendars
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

    from pandas.tseries.holiday import (
        Holiday,
        USMemorialDay,
        AbstractHolidayCalendar,
        nearest_workday,
        MO,
    )

    class ExampleCalendar(AbstractHolidayCalendar):
        rules = [
            USMemorialDay,
            Holiday("July 4th", month=7, day=4, observance=nearest_workday),
            Holiday(
                "Columbus Day",
                month=10,
                day=1,
                offset=pd.DateOffset(weekday=MO(2)),
            ),
        ]

    cal = ExampleCalendar()
    cal.holidays(datetime.datetime(2012, 1, 1), datetime.datetime(2012, 12, 31))

:hint:
   **weekday=MO(2)** is same as **2 * Week(weekday=2)**

Using this calendar, creating an index or doing offset arithmetic skips weekends
and holidays (i.e., Memorial Day/July 4th).  For example, the below defines
a custom business day offset using the ``ExampleCalendar``.  Like any other offset,
it can be used to create a ``DatetimeIndex`` or added to ``datetime``
or ``Timestamp`` objects.

.. ipython:: python

    pd.date_range(
        start="7/1/2012", end="7/10/2012", freq=pd.offsets.CDay(calendar=cal)
    ).to_pydatetime()
    offset = pd.offsets.CustomBusinessDay(calendar=cal)
    datetime.datetime(2012, 5, 25) + offset
    datetime.datetime(2012, 7, 3) + offset
    datetime.datetime(2012, 7, 3) + 2 * offset
    datetime.datetime(2012, 7, 6) + offset

Ranges are defined by the ``start_date`` and ``end_date`` class attributes
of ``AbstractHolidayCalendar``.  The defaults are shown below.

.. ipython:: python

    AbstractHolidayCalendar.start_date
    AbstractHolidayCalendar.end_date

These dates can be overwritten by setting the attributes as
datetime/Timestamp/string.

.. ipython:: python

    AbstractHolidayCalendar.start_date = datetime.datetime(2012, 1, 1)
    AbstractHolidayCalendar.end_date = datetime.datetime(2012, 12, 31)
    cal.holidays()

Every calendar class is accessible by name using the ``get_calendar`` function
which returns a holiday class instance.  Any imported calendar class will
automatically be available by this function.  Also, ``HolidayCalendarFactory``
provides an easy interface to create calendars that are combinations of calendars
or calendars with additional rules.

.. ipython:: python

    from pandas.tseries.holiday import get_calendar, HolidayCalendarFactory, USLaborDay

    cal = get_calendar("ExampleCalendar")
    cal.rules
    new_cal = HolidayCalendarFactory("NewExampleCalendar", cal, USLaborDay)
    new_cal.rules

.. _timeseries.advanced_datetime:

Time series-related instance methods
------------------------------------

Shifting / lagging
~~~~~~~~~~~~~~~~~~

One may want to *shift* or *lag* the values in a time series back and forward in
time. The method for this is :meth:`~Series.shift`, which is available on all of
the pandas objects.

.. ipython:: python

   ts = pd.Series(range(len(rng)), index=rng)
   ts = ts[:5]
   ts.shift(1)

The ``shift`` method accepts an ``freq`` argument which can accept a
``DateOffset`` class or other ``timedelta``-like object or also an
:ref:`offset alias <timeseries.offset_aliases>`.

When ``freq`` is specified, ``shift`` method changes all the dates in the index
rather than changing the alignment of the data and the index:

.. ipython:: python

   ts.shift(5, freq="D")
   ts.shift(5, freq=pd.offsets.BDay())
   ts.shift(5, freq="BM")

Note that with when ``freq`` is specified, the leading entry is no longer NaN
because the data is not being realigned.

Frequency conversion
~~~~~~~~~~~~~~~~~~~~

The primary function for changing frequencies is the :meth:`~Series.asfreq`
method. For a ``DatetimeIndex``, this is basically just a thin, but convenient
wrapper around :meth:`~Series.reindex`  which generates a ``date_range`` and
calls ``reindex``.

.. ipython:: python

   dr = pd.date_range("1/1/2010", periods=3, freq=3 * pd.offsets.BDay())
   ts = pd.Series(np.random.randn(3), index=dr)
   ts
   ts.asfreq(pd.offsets.BDay())

``asfreq`` provides a further convenience so you can specify an interpolation
method for any gaps that may appear after the frequency conversion.

.. ipython:: python

   ts.asfreq(pd.offsets.BDay(), method="pad")

Filling forward / backward
~~~~~~~~~~~~~~~~~~~~~~~~~~

Related to ``asfreq`` and ``reindex`` is :meth:`~Series.fillna`, which is
documented in the :ref:`missing data section <missing_data.fillna>`.

Converting to Python datetimes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``DatetimeIndex`` can be converted to an array of Python native
:py:class:`datetime.datetime` objects using the ``to_pydatetime`` method.

.. _timeseries.resampling:

Resampling
----------

pandas has a simple, powerful, and efficient functionality for performing
resampling operations during frequency conversion (e.g., converting secondly
data into 5-minutely data). This is extremely common in, but not limited to,
financial applications.

:meth:`~Series.resample` is a time-based groupby, followed by a reduction method
on each of its groups. See some :ref:`cookbook examples <cookbook.resample>` for
some advanced strategies.

The ``resample()`` method can be used directly from ``DataFrameGroupBy`` objects,
see the :ref:`groupby docs <groupby.transform.window_resample>`.

Basics
~~~~~~

.. ipython:: python

   rng = pd.date_range("1/1/2012", periods=100, freq="S")

   ts = pd.Series(np.random.randint(0, 500, len(rng)), index=rng)

   ts.resample("5Min").sum()

The ``resample`` function is very flexible and allows you to specify many
different parameters to control the frequency conversion and resampling
operation.

Any function available via :ref:`dispatching <groupby.dispatch>` is available as
a method of the returned object, including ``sum``, ``mean``, ``std``, ``sem``,
``max``, ``min``, ``median``, ``first``, ``last``, ``ohlc``:

.. ipython:: python

   ts.resample("5Min").mean()

   ts.resample("5Min").ohlc()

   ts.resample("5Min").max()


For downsampling, ``closed`` can be set to 'left' or 'right' to specify which
end of the interval is closed:

.. ipython:: python

   ts.resample("5Min", closed="right").mean()

   ts.resample("5Min", closed="left").mean()

Parameters like ``label`` are used to manipulate the resulting labels.
``label`` specifies whether the result is labeled with the beginning or
the end of the interval.

.. ipython:: python

   ts.resample("5Min").mean()  # by default label='left'

   ts.resample("5Min", label="left").mean()

.. warning::

    The default values for ``label`` and ``closed`` is '**left**' for all
    frequency offsets except for 'M', 'A', 'Q', 'BM', 'BA', 'BQ', and 'W'
    which all have a default of 'right'.

    This might unintendedly lead to looking ahead, where the value for a later
    time is pulled back to a previous time as in the following example with
    the :class:`~pandas.tseries.offsets.BusinessDay` frequency:

    .. ipython:: python

        s = pd.date_range("2000-01-01", "2000-01-05").to_series()
        s.iloc[2] = pd.NaT
        s.dt.day_name()

        # default: label='left', closed='left'
        s.resample("B").last().dt.day_name()

    Notice how the value for Sunday got pulled back to the previous Friday.
    To get the behavior where the value for Sunday is pushed to Monday, use
    instead

    .. ipython:: python

        s.resample("B", label="right", closed="right").last().dt.day_name()

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

   ts[:2].resample("250L").asfreq()

   ts[:2].resample("250L").ffill()

   ts[:2].resample("250L").ffill(limit=2)

Sparse resampling
~~~~~~~~~~~~~~~~~

Sparse timeseries are the ones where you have a lot fewer points relative
to the amount of time you are looking to resample. Naively upsampling a sparse
series can potentially generate lots of intermediate values. When you don't want
to use a method to fill these values, e.g. ``fill_method`` is ``None``, then
intermediate values will be filled with ``NaN``.

Since ``resample`` is a time-based groupby, the following is a method to efficiently
resample only the groups that are not all ``NaN``.

.. ipython:: python

    rng = pd.date_range("2014-1-1", periods=100, freq="D") + pd.Timedelta("1s")
    ts = pd.Series(range(100), index=rng)

If we want to resample to the full range of the series:

.. ipython:: python

    ts.resample("3T").sum()

We can instead only resample those groups where we have points as follows:

.. ipython:: python

    from functools import partial
    from pandas.tseries.frequencies import to_offset

    def round(t, freq):
        # round a Timestamp to a specified freq
        freq = to_offset(freq)
        return pd.Timestamp((t.value // freq.delta.value) * freq.delta.value)

    ts.groupby(partial(round, freq="3T")).sum()

.. _timeseries.aggregate:

Aggregation
~~~~~~~~~~~

Similar to the :ref:`aggregating API <basics.aggregate>`, :ref:`groupby API <groupby.aggregate>`, and the :ref:`window API <window.overview>`,
a ``Resampler`` can be selectively resampled.

Resampling a ``DataFrame``, the default will be to act on all columns with the same function.

.. ipython:: python

   df = pd.DataFrame(
       np.random.randn(1000, 3),
       index=pd.date_range("1/1/2012", freq="S", periods=1000),
       columns=["A", "B", "C"],
   )
   r = df.resample("3T")
   r.mean()

We can select a specific column or columns using standard getitem.

.. ipython:: python

   r["A"].mean()

   r[["A", "B"]].mean()

You can pass a list or dict of functions to do aggregation with, outputting a ``DataFrame``:

.. ipython:: python

   r["A"].agg([np.sum, np.mean, np.std])

On a resampled ``DataFrame``, you can pass a list of functions to apply to each
column, which produces an aggregated result with a hierarchical index:

.. ipython:: python

   r.agg([np.sum, np.mean])

By passing a dict to ``aggregate`` you can apply a different aggregation to the
columns of a ``DataFrame``:

.. ipython:: python
   :okexcept:

   r.agg({"A": np.sum, "B": lambda x: np.std(x, ddof=1)})

The function names can also be strings. In order for a string to be valid it
must be implemented on the resampled object:

.. ipython:: python

   r.agg({"A": "sum", "B": "std"})

Furthermore, you can also specify multiple aggregation functions for each column separately.

.. ipython:: python

   r.agg({"A": ["sum", "std"], "B": ["mean", "std"]})


If a ``DataFrame`` does not have a datetimelike index, but instead you want
to resample based on datetimelike column in the frame, it can passed to the
``on`` keyword.

.. ipython:: python

   df = pd.DataFrame(
       {"date": pd.date_range("2015-01-01", freq="W", periods=5), "a": np.arange(5)},
       index=pd.MultiIndex.from_arrays(
           [[1, 2, 3, 4, 5], pd.date_range("2015-01-01", freq="W", periods=5)],
           names=["v", "d"],
       ),
   )
   df
   df.resample("M", on="date").sum()

Similarly, if you instead want to resample by a datetimelike
level of ``MultiIndex``, its name or location can be passed to the
``level`` keyword.

.. ipython:: python

   df.resample("M", level="d").sum()

.. _timeseries.iterating-label:

Iterating through groups
~~~~~~~~~~~~~~~~~~~~~~~~

With the ``Resampler`` object in hand, iterating through the grouped data is very
natural and functions similarly to :py:func:`itertools.groupby`:

.. ipython:: python

   small = pd.Series(
       range(6),
       index=pd.to_datetime(
           [
               "2017-01-01T00:00:00",
               "2017-01-01T00:30:00",
               "2017-01-01T00:31:00",
               "2017-01-01T01:00:00",
               "2017-01-01T03:00:00",
               "2017-01-01T03:05:00",
           ]
       ),
   )
   resampled = small.resample("H")

   for name, group in resampled:
       print("Group: ", name)
       print("-" * 27)
       print(group, end="\n\n")

See :ref:`groupby.iterating-label` or :class:`Resampler.__iter__` for more.

.. _timeseries.adjust-the-start-of-the-bins:

Use ``origin`` or ``offset`` to adjust the start of the bins
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 1.1.0

The bins of the grouping are adjusted based on the beginning of the day of the time series starting point. This works well with frequencies that are multiples of a day (like ``30D``) or that divide a day evenly (like ``90s`` or ``1min``). This can create inconsistencies with some frequencies that do not meet this criteria. To change this behavior you can specify a fixed Timestamp with the argument ``origin``.

For example:

.. ipython:: python

    start, end = "2000-10-01 23:30:00", "2000-10-02 00:30:00"
    middle = "2000-10-02 00:00:00"
    rng = pd.date_range(start, end, freq="7min")
    ts = pd.Series(np.arange(len(rng)) * 3, index=rng)
    ts

Here we can see that, when using ``origin`` with its default value (``'start_day'``), the result after ``'2000-10-02 00:00:00'`` are not identical depending on the start of time series:

.. ipython:: python

    ts.resample("17min", origin="start_day").sum()
    ts[middle:end].resample("17min", origin="start_day").sum()


Here we can see that, when setting ``origin`` to ``'epoch'``, the result after ``'2000-10-02 00:00:00'`` are identical depending on the start of time series:

.. ipython:: python

   ts.resample("17min", origin="epoch").sum()
   ts[middle:end].resample("17min", origin="epoch").sum()


If needed you can use a custom timestamp for ``origin``:

.. ipython:: python

   ts.resample("17min", origin="2001-01-01").sum()
   ts[middle:end].resample("17min", origin=pd.Timestamp("2001-01-01")).sum()

If needed you can just adjust the bins with an ``offset`` Timedelta that would be added to the default ``origin``.
Those two examples are equivalent for this time series:

.. ipython:: python

    ts.resample("17min", origin="start").sum()
    ts.resample("17min", offset="23h30min").sum()


Note the use of ``'start'`` for ``origin`` on the last example. In that case, ``origin`` will be set to the first value of the timeseries.

Backward resample
~~~~~~~~~~~~~~~~~

.. versionadded:: 1.3.0

Instead of adjusting the beginning of bins, sometimes we need to fix the end of the bins to make a backward resample with a given ``freq``. The backward resample sets ``closed`` to ``'right'`` by default since the last value should be considered as the edge point for the last bin.

We can set ``origin`` to ``'end'``. The value for a specific ``Timestamp`` index stands for the resample result from the current ``Timestamp`` minus ``freq`` to the current ``Timestamp`` with a right close.

.. ipython:: python

   ts.resample('17min', origin='end').sum()

Besides, in contrast with the ``'start_day'`` option, ``end_day`` is supported. This will set the origin as the ceiling midnight of the largest ``Timestamp``.

.. ipython:: python

   ts.resample('17min', origin='end_day').sum()

The above result uses ``2000-10-02 00:29:00`` as the last bin's right edge since the following computation.

.. ipython:: python

   ceil_mid = rng.max().ceil('D')
   freq = pd.offsets.Minute(17)
   bin_res = ceil_mid - freq * ((ceil_mid - rng.max()) // freq)
   bin_res

.. _timeseries.periods:

Time span representation
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

   pd.Period("2012", freq="A-DEC")

   pd.Period("2012-1-1", freq="D")

   pd.Period("2012-1-1 19:00", freq="H")

   pd.Period("2012-1-1 19:00", freq="5H")

Adding and subtracting integers from periods shifts the period by its own
frequency. Arithmetic is not allowed between ``Period`` with different ``freq`` (span).

.. ipython:: python

   p = pd.Period("2012", freq="A-DEC")
   p + 1
   p - 3
   p = pd.Period("2012-01", freq="2M")
   p + 2
   p - 1
   @okexcept
   p == pd.Period("2012-01", freq="3M")


If ``Period`` freq is daily or higher (``D``, ``H``, ``T``, ``S``, ``L``, ``U``, ``N``), ``offsets`` and ``timedelta``-like can be added if the result can have the same freq. Otherwise, ``ValueError`` will be raised.

.. ipython:: python

   p = pd.Period("2014-07-01 09:00", freq="H")
   p + pd.offsets.Hour(2)
   p + datetime.timedelta(minutes=120)
   p + np.timedelta64(7200, "s")

.. code-block:: ipython

   In [1]: p + pd.offsets.Minute(5)
   Traceback
      ...
   ValueError: Input has different freq from Period(freq=H)

If ``Period`` has other frequencies, only the same ``offsets`` can be added. Otherwise, ``ValueError`` will be raised.

.. ipython:: python

   p = pd.Period("2014-07", freq="M")
   p + pd.offsets.MonthEnd(3)

.. code-block:: ipython

   In [1]: p + pd.offsets.MonthBegin(3)
   Traceback
      ...
   ValueError: Input has different freq from Period(freq=M)

Taking the difference of ``Period`` instances with the same frequency will
return the number of frequency units between them:

.. ipython:: python

   pd.Period("2012", freq="A-DEC") - pd.Period("2002", freq="A-DEC")

PeriodIndex and period_range
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Regular sequences of ``Period`` objects can be collected in a ``PeriodIndex``,
which can be constructed using the ``period_range`` convenience function:

.. ipython:: python

   prng = pd.period_range("1/1/2011", "1/1/2012", freq="M")
   prng

The ``PeriodIndex`` constructor can also be used directly:

.. ipython:: python

   pd.PeriodIndex(["2011-1", "2011-2", "2011-3"], freq="M")

Passing multiplied frequency outputs a sequence of ``Period`` which
has multiplied span.

.. ipython:: python

   pd.period_range(start="2014-01", freq="3M", periods=4)

If ``start`` or ``end`` are ``Period`` objects, they will be used as anchor
endpoints for a ``PeriodIndex`` with frequency matching that of the
``PeriodIndex`` constructor.

.. ipython:: python

   pd.period_range(
       start=pd.Period("2017Q1", freq="Q"), end=pd.Period("2017Q2", freq="Q"), freq="M"
   )

Just like ``DatetimeIndex``, a ``PeriodIndex`` can also be used to index pandas
objects:

.. ipython:: python

   ps = pd.Series(np.random.randn(len(prng)), prng)
   ps

``PeriodIndex`` supports addition and subtraction with the same rule as ``Period``.

.. ipython:: python

   idx = pd.period_range("2014-07-01 09:00", periods=5, freq="H")
   idx
   idx + pd.offsets.Hour(2)

   idx = pd.period_range("2014-07", periods=5, freq="M")
   idx
   idx + pd.offsets.MonthEnd(3)

``PeriodIndex`` has its own dtype named ``period``, refer to :ref:`Period Dtypes <timeseries.period_dtype>`.

.. _timeseries.period_dtype:

Period dtypes
~~~~~~~~~~~~~

``PeriodIndex`` has a custom ``period`` dtype. This is a pandas extension
dtype similar to the :ref:`timezone aware dtype <timeseries.timezone_series>` (``datetime64[ns, tz]``).

The ``period`` dtype holds the ``freq`` attribute and is represented with
``period[freq]`` like ``period[D]`` or ``period[M]``, using :ref:`frequency strings <timeseries.offset_aliases>`.

.. ipython:: python

   pi = pd.period_range("2016-01-01", periods=3, freq="M")
   pi
   pi.dtype

The ``period`` dtype can be used in ``.astype(...)``. It allows one to change the
``freq`` of a ``PeriodIndex`` like ``.asfreq()`` and convert a
``DatetimeIndex`` to ``PeriodIndex`` like ``to_period()``:

.. ipython:: python

   # change monthly freq to daily freq
   pi.astype("period[D]")

   # convert to DatetimeIndex
   pi.astype("datetime64[ns]")

   # convert to PeriodIndex
   dti = pd.date_range("2011-01-01", freq="M", periods=3)
   dti
   dti.astype("period[M]")

PeriodIndex partial string indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PeriodIndex now supports partial string slicing with non-monotonic indexes.

.. versionadded:: 1.1.0

You can pass in dates and strings to ``Series`` and ``DataFrame`` with ``PeriodIndex``, in the same manner as ``DatetimeIndex``. For details, refer to :ref:`DatetimeIndex Partial String Indexing <timeseries.partialindexing>`.

.. ipython:: python

   ps["2011-01"]

   ps[datetime.datetime(2011, 12, 25):]

   ps["10/31/2011":"12/31/2011"]

Passing a string representing a lower frequency than ``PeriodIndex`` returns partial sliced data.

.. ipython:: python

   ps["2011"]

   dfp = pd.DataFrame(
       np.random.randn(600, 1),
       columns=["A"],
       index=pd.period_range("2013-01-01 9:00", periods=600, freq="T"),
   )
   dfp
   dfp.loc["2013-01-01 10H"]

As with ``DatetimeIndex``, the endpoints will be included in the result. The example below slices data starting from 10:00 to 11:59.

.. ipython:: python

   dfp["2013-01-01 10H":"2013-01-01 11H"]


Frequency conversion and resampling with PeriodIndex
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The frequency of ``Period`` and ``PeriodIndex`` can be converted via the ``asfreq``
method. Let's start with the fiscal year 2011, ending in December:

.. ipython:: python

   p = pd.Period("2011", freq="A-DEC")
   p

We can convert it to a monthly frequency. Using the ``how`` parameter, we can
specify whether to return the starting or ending month:

.. ipython:: python

   p.asfreq("M", how="start")

   p.asfreq("M", how="end")

The shorthands 's' and 'e' are provided for convenience:

.. ipython:: python

   p.asfreq("M", "s")
   p.asfreq("M", "e")

Converting to a "super-period" (e.g., annual frequency is a super-period of
quarterly frequency) automatically returns the super-period that includes the
input period:

.. ipython:: python

   p = pd.Period("2011-12", freq="M")

   p.asfreq("A-NOV")

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

   p = pd.Period("2012Q1", freq="Q-DEC")

   p.asfreq("D", "s")

   p.asfreq("D", "e")

``Q-MAR`` defines fiscal year end in March:

.. ipython:: python

   p = pd.Period("2011Q4", freq="Q-MAR")

   p.asfreq("D", "s")

   p.asfreq("D", "e")

.. _timeseries.interchange:

Converting between representations
----------------------------------

Timestamped data can be converted to PeriodIndex-ed data using ``to_period``
and vice-versa using ``to_timestamp``:

.. ipython:: python

   rng = pd.date_range("1/1/2012", periods=5, freq="M")

   ts = pd.Series(np.random.randn(len(rng)), index=rng)

   ts

   ps = ts.to_period()

   ps

   ps.to_timestamp()

Remember that 's' and 'e' can be used to return the timestamps at the start or
end of the period:

.. ipython:: python

   ps.to_timestamp("D", how="s")

Converting between period and timestamp enables some convenient arithmetic
functions to be used. In the following example, we convert a quarterly
frequency with year ending in November to 9am of the end of the month following
the quarter end:

.. ipython:: python

   prng = pd.period_range("1990Q1", "2000Q4", freq="Q-NOV")

   ts = pd.Series(np.random.randn(len(prng)), prng)

   ts.index = (prng.asfreq("M", "e") + 1).asfreq("H", "s") + 9

   ts.head()

.. _timeseries.oob:

Representing out-of-bounds spans
--------------------------------

If you have data that is outside of the ``Timestamp`` bounds, see :ref:`Timestamp limitations <timeseries.timestamp-limits>`,
then you can use a ``PeriodIndex`` and/or ``Series`` of ``Periods`` to do computations.

.. ipython:: python

   span = pd.period_range("1215-01-01", "1381-01-01", freq="D")
   span

To convert from an ``int64`` based YYYYMMDD representation.

.. ipython:: python

   s = pd.Series([20121231, 20141130, 99991231])
   s

   def conv(x):
       return pd.Period(year=x // 10000, month=x // 100 % 100, day=x % 100, freq="D")

   s.apply(conv)
   s.apply(conv)[2]

These can easily be converted to a ``PeriodIndex``:

.. ipython:: python

   span = pd.PeriodIndex(s.apply(conv))
   span

.. _timeseries.timezone:

Time zone handling
------------------

pandas provides rich support for working with timestamps in different time
zones using the ``pytz`` and ``dateutil`` libraries or :class:`datetime.timezone`
objects from the standard library.


Working with time zones
~~~~~~~~~~~~~~~~~~~~~~~

By default, pandas objects are time zone unaware:

.. ipython:: python

   rng = pd.date_range("3/6/2012 00:00", periods=15, freq="D")
   rng.tz is None

To localize these dates to a time zone (assign a particular time zone to a naive date),
you can use the ``tz_localize`` method or the ``tz`` keyword argument in
:func:`date_range`, :class:`Timestamp`, or :class:`DatetimeIndex`.
You can either pass ``pytz`` or ``dateutil`` time zone objects or Olson time zone database strings.
Olson time zone strings will return ``pytz`` time zone objects by default.
To return ``dateutil`` time zone objects, append ``dateutil/`` before the string.

* In ``pytz`` you can find a list of common (and less common) time zones using
  ``from pytz import common_timezones, all_timezones``.
* ``dateutil`` uses the OS time zones so there isn't a fixed list available. For
  common zones, the names are the same as ``pytz``.

.. ipython:: python

   import dateutil

   # pytz
   rng_pytz = pd.date_range("3/6/2012 00:00", periods=3, freq="D", tz="Europe/London")
   rng_pytz.tz

   # dateutil
   rng_dateutil = pd.date_range("3/6/2012 00:00", periods=3, freq="D")
   rng_dateutil = rng_dateutil.tz_localize("dateutil/Europe/London")
   rng_dateutil.tz

   # dateutil - utc special case
   rng_utc = pd.date_range(
       "3/6/2012 00:00",
       periods=3,
       freq="D",
       tz=dateutil.tz.tzutc(),
   )
   rng_utc.tz

.. versionadded:: 0.25.0

.. ipython:: python

   # datetime.timezone
   rng_utc = pd.date_range(
       "3/6/2012 00:00",
       periods=3,
       freq="D",
       tz=datetime.timezone.utc,
   )
   rng_utc.tz

Note that the ``UTC`` time zone is a special case in ``dateutil`` and should be constructed explicitly
as an instance of ``dateutil.tz.tzutc``. You can also construct other time
zones objects explicitly first.

.. ipython:: python

   import pytz

   # pytz
   tz_pytz = pytz.timezone("Europe/London")
   rng_pytz = pd.date_range("3/6/2012 00:00", periods=3, freq="D")
   rng_pytz = rng_pytz.tz_localize(tz_pytz)
   rng_pytz.tz == tz_pytz

   # dateutil
   tz_dateutil = dateutil.tz.gettz("Europe/London")
   rng_dateutil = pd.date_range("3/6/2012 00:00", periods=3, freq="D", tz=tz_dateutil)
   rng_dateutil.tz == tz_dateutil

To convert a time zone aware pandas object from one time zone to another,
you can use the ``tz_convert`` method.

.. ipython:: python

   rng_pytz.tz_convert("US/Eastern")

.. note::

    When using ``pytz`` time zones, :class:`DatetimeIndex` will construct a different
    time zone object than a :class:`Timestamp` for the same time zone input. A :class:`DatetimeIndex`
    can hold a collection of :class:`Timestamp` objects that may have different UTC offsets and cannot be
    succinctly represented by one ``pytz`` time zone instance while one :class:`Timestamp`
    represents one point in time with a specific UTC offset.

    .. ipython:: python

       dti = pd.date_range("2019-01-01", periods=3, freq="D", tz="US/Pacific")
       dti.tz
       ts = pd.Timestamp("2019-01-01", tz="US/Pacific")
       ts.tz

.. warning::

        Be wary of conversions between libraries. For some time zones, ``pytz`` and ``dateutil`` have different
        definitions of the zone. This is more of a problem for unusual time zones than for
        'standard' zones like ``US/Eastern``.

.. warning::

    Be aware that a time zone definition across versions of time zone libraries may not
    be considered equal.  This may cause problems when working with stored data that
    is localized using one version and operated on with a different version.
    See :ref:`here<io.hdf5-notes>` for how to handle such a situation.

.. warning::

    For ``pytz`` time zones, it is incorrect to pass a time zone object directly into
    the ``datetime.datetime`` constructor
    (e.g., ``datetime.datetime(2011, 1, 1, tzinfo=pytz.timezone('US/Eastern'))``.
    Instead, the datetime needs to be localized using the ``localize`` method
    on the ``pytz`` time zone object.

.. warning::

    Be aware that for times in the future, correct conversion between time zones
    (and UTC) cannot be guaranteed by any time zone library because a timezone's
    offset from UTC may be changed by the respective government.

.. warning::

    If you are using dates beyond 2038-01-18, due to current deficiencies
    in the underlying libraries caused by the year 2038 problem, daylight saving time (DST) adjustments
    to timezone aware dates will not be applied. If and when the underlying libraries are fixed,
    the DST transitions will be applied.

    For example, for two dates that are in British Summer Time (and so would normally be GMT+1), both the following asserts evaluate as true:

    .. ipython:: python

       d_2037 = "2037-03-31T010101"
       d_2038 = "2038-03-31T010101"
       DST = "Europe/London"
       assert pd.Timestamp(d_2037, tz=DST) != pd.Timestamp(d_2037, tz="GMT")
       assert pd.Timestamp(d_2038, tz=DST) == pd.Timestamp(d_2038, tz="GMT")

Under the hood, all timestamps are stored in UTC. Values from a time zone aware
:class:`DatetimeIndex` or :class:`Timestamp` will have their fields (day, hour, minute, etc.)
localized to the time zone. However, timestamps with the same UTC value are
still considered to be equal even if they are in different time zones:

.. ipython:: python

   rng_eastern = rng_utc.tz_convert("US/Eastern")
   rng_berlin = rng_utc.tz_convert("Europe/Berlin")

   rng_eastern[2]
   rng_berlin[2]
   rng_eastern[2] == rng_berlin[2]

Operations between :class:`Series` in different time zones will yield UTC
:class:`Series`, aligning the data on the UTC timestamps:

.. ipython:: python

   ts_utc = pd.Series(range(3), pd.date_range("20130101", periods=3, tz="UTC"))
   eastern = ts_utc.tz_convert("US/Eastern")
   berlin = ts_utc.tz_convert("Europe/Berlin")
   result = eastern + berlin
   result
   result.index

To remove time zone information, use ``tz_localize(None)`` or ``tz_convert(None)``.
``tz_localize(None)`` will remove the time zone yielding the local time representation.
``tz_convert(None)`` will remove the time zone after converting to UTC time.

.. ipython:: python

   didx = pd.date_range(start="2014-08-01 09:00", freq="H", periods=3, tz="US/Eastern")
   didx
   didx.tz_localize(None)
   didx.tz_convert(None)

   # tz_convert(None) is identical to tz_convert('UTC').tz_localize(None)
   didx.tz_convert("UTC").tz_localize(None)

.. _timeseries.fold:

Fold
~~~~

.. versionadded:: 1.1.0

For ambiguous times, pandas supports explicitly specifying the keyword-only fold argument.
Due to daylight saving time, one wall clock time can occur twice when shifting
from summer to winter time; fold describes whether the datetime-like corresponds
to the first (0) or the second time (1) the wall clock hits the ambiguous time.
Fold is supported only for constructing from naive ``datetime.datetime``
(see `datetime documentation <https://docs.python.org/3/library/datetime.html>`__ for details) or from :class:`Timestamp`
or for constructing from components (see below). Only ``dateutil`` timezones are supported
(see `dateutil documentation <https://dateutil.readthedocs.io/en/stable/tz.html#dateutil.tz.enfold>`__
for ``dateutil`` methods that deal with ambiguous datetimes) as ``pytz``
timezones do not support fold (see `pytz documentation <http://pytz.sourceforge.net/index.html>`__
for details on how ``pytz`` deals with ambiguous datetimes). To localize an ambiguous datetime
with ``pytz``, please use :meth:`Timestamp.tz_localize`. In general, we recommend to rely
on :meth:`Timestamp.tz_localize` when localizing ambiguous datetimes if you need direct
control over how they are handled.

.. ipython:: python

   pd.Timestamp(
       datetime.datetime(2019, 10, 27, 1, 30, 0, 0),
       tz="dateutil/Europe/London",
       fold=0,
   )
   pd.Timestamp(
       year=2019,
       month=10,
       day=27,
       hour=1,
       minute=30,
       tz="dateutil/Europe/London",
       fold=1,
   )

.. _timeseries.timezone_ambiguous:

Ambiguous times when localizing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``tz_localize`` may not be able to determine the UTC offset of a timestamp
because daylight savings time (DST) in a local time zone causes some times to occur
twice within one day ("clocks fall back"). The following options are available:

* ``'raise'``: Raises a ``pytz.AmbiguousTimeError`` (the default behavior)
* ``'infer'``: Attempt to determine the correct offset base on the monotonicity of the timestamps
* ``'NaT'``: Replaces ambiguous times with ``NaT``
* ``bool``: ``True`` represents a DST time, ``False`` represents non-DST time. An array-like of ``bool`` values is supported for a sequence of times.

.. ipython:: python

   rng_hourly = pd.DatetimeIndex(
       ["11/06/2011 00:00", "11/06/2011 01:00", "11/06/2011 01:00", "11/06/2011 02:00"]
   )

This will fail as there are ambiguous times (``'11/06/2011 01:00'``)

.. code-block:: ipython

   In [2]: rng_hourly.tz_localize('US/Eastern')
   AmbiguousTimeError: Cannot infer dst time from Timestamp('2011-11-06 01:00:00'), try using the 'ambiguous' argument

Handle these ambiguous times by specifying the following.

.. ipython:: python

   rng_hourly.tz_localize("US/Eastern", ambiguous="infer")
   rng_hourly.tz_localize("US/Eastern", ambiguous="NaT")
   rng_hourly.tz_localize("US/Eastern", ambiguous=[True, True, False, False])

.. _timeseries.timezone_nonexistent:

Nonexistent times when localizing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A DST transition may also shift the local time ahead by 1 hour creating nonexistent
local times ("clocks spring forward"). The behavior of localizing a timeseries with nonexistent times
can be controlled by the ``nonexistent`` argument. The following options are available:

* ``'raise'``: Raises a ``pytz.NonExistentTimeError`` (the default behavior)
* ``'NaT'``: Replaces nonexistent times with ``NaT``
* ``'shift_forward'``: Shifts nonexistent times forward to the closest real time
* ``'shift_backward'``: Shifts nonexistent times backward to the closest real time
* timedelta object: Shifts nonexistent times by the timedelta duration

.. ipython:: python

    dti = pd.date_range(start="2015-03-29 02:30:00", periods=3, freq="H")
    # 2:30 is a nonexistent time

Localization of nonexistent times will raise an error by default.

.. code-block:: ipython

   In [2]: dti.tz_localize('Europe/Warsaw')
   NonExistentTimeError: 2015-03-29 02:30:00

Transform nonexistent times to ``NaT`` or shift the times.

.. ipython:: python

    dti
    dti.tz_localize("Europe/Warsaw", nonexistent="shift_forward")
    dti.tz_localize("Europe/Warsaw", nonexistent="shift_backward")
    dti.tz_localize("Europe/Warsaw", nonexistent=pd.Timedelta(1, unit="H"))
    dti.tz_localize("Europe/Warsaw", nonexistent="NaT")


.. _timeseries.timezone_series:

Time zone series operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :class:`Series` with time zone **naive** values is
represented with a dtype of ``datetime64[ns]``.

.. ipython:: python

   s_naive = pd.Series(pd.date_range("20130101", periods=3))
   s_naive

A :class:`Series` with a time zone **aware** values is
represented with a dtype of ``datetime64[ns, tz]`` where ``tz`` is the time zone

.. ipython:: python

   s_aware = pd.Series(pd.date_range("20130101", periods=3, tz="US/Eastern"))
   s_aware

Both of these :class:`Series` time zone information
can be manipulated via the ``.dt`` accessor, see :ref:`the dt accessor section <basics.dt_accessors>`.

For example, to localize and convert a naive stamp to time zone aware.

.. ipython:: python

   s_naive.dt.tz_localize("UTC").dt.tz_convert("US/Eastern")

Time zone information can also be manipulated using the ``astype`` method.
This method can convert between different timezone-aware dtypes.

.. ipython:: python

   # convert to a new time zone
   s_aware.astype("datetime64[ns, CET]")

.. note::

   Using :meth:`Series.to_numpy` on a ``Series``, returns a NumPy array of the data.
   NumPy does not currently support time zones (even though it is *printing* in the local time zone!),
   therefore an object array of Timestamps is returned for time zone aware data:

   .. ipython:: python

      s_naive.to_numpy()
      s_aware.to_numpy()

   By converting to an object array of Timestamps, it preserves the time zone
   information. For example, when converting back to a Series:

   .. ipython:: python

      pd.Series(s_aware.to_numpy())

   However, if you want an actual NumPy ``datetime64[ns]`` array (with the values
   converted to UTC) instead of an array of objects, you can specify the
   ``dtype`` argument:

   .. ipython:: python

      s_aware.to_numpy(dtype="datetime64[ns]")
