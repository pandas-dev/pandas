{{ header }}

.. _api.scalars:

=======
Scalars
=======
.. currentmodule:: pandas

Period
------
.. autosummary::
    :toctree: generated/

    Period

Properties
~~~~~~~~~~
.. autosummary::
    :toctree: generated/

    Period.day
    Period.dayofweek
    Period.dayofyear
    Period.days_in_month
    Period.daysinmonth
    Period.end_time
    Period.freq
    Period.freqstr
    Period.hour
    Period.is_leap_year
    Period.minute
    Period.month
    Period.ordinal
    Period.quarter
    Period.qyear
    Period.second
    Period.start_time
    Period.week
    Period.weekday
    Period.weekofyear
    Period.year

Methods
~~~~~~~
.. autosummary::
    :toctree: generated/

    Period.asfreq
    Period.now
    Period.strftime
    Period.to_timestamp

Timestamp
---------
.. autosummary::
    :toctree: generated/

    Timestamp

Properties
~~~~~~~~~~
.. autosummary::
    :toctree: generated/

    Timestamp.asm8
    Timestamp.day
    Timestamp.dayofweek
    Timestamp.dayofyear
    Timestamp.days_in_month
    Timestamp.daysinmonth
    Timestamp.fold
    Timestamp.hour
    Timestamp.is_leap_year
    Timestamp.is_month_end
    Timestamp.is_month_start
    Timestamp.is_quarter_end
    Timestamp.is_quarter_start
    Timestamp.is_year_end
    Timestamp.is_year_start
    Timestamp.max
    Timestamp.microsecond
    Timestamp.min
    Timestamp.minute
    Timestamp.month
    Timestamp.nanosecond
    Timestamp.quarter
    Timestamp.resolution
    Timestamp.second
    Timestamp.tz
    Timestamp.tzinfo
    Timestamp.value
    Timestamp.week
    Timestamp.weekofyear
    Timestamp.year

Methods
~~~~~~~
.. autosummary::
    :toctree: generated/

    Timestamp.astimezone
    Timestamp.ceil
    Timestamp.combine
    Timestamp.ctime
    Timestamp.date
    Timestamp.day_name
    Timestamp.dst
    Timestamp.floor
    Timestamp.freq
    Timestamp.freqstr
    Timestamp.fromordinal
    Timestamp.fromtimestamp
    Timestamp.isocalendar
    Timestamp.isoformat
    Timestamp.isoweekday
    Timestamp.month_name
    Timestamp.normalize
    Timestamp.now
    Timestamp.replace
    Timestamp.round
    Timestamp.strftime
    Timestamp.strptime
    Timestamp.time
    Timestamp.timestamp
    Timestamp.timetuple
    Timestamp.timetz
    Timestamp.to_datetime64
    Timestamp.to_julian_date
    Timestamp.to_period
    Timestamp.to_pydatetime
    Timestamp.today
    Timestamp.toordinal
    Timestamp.tz_convert
    Timestamp.tz_localize
    Timestamp.tzname
    Timestamp.utcfromtimestamp
    Timestamp.utcnow
    Timestamp.utcoffset
    Timestamp.utctimetuple
    Timestamp.weekday

Interval
--------
.. autosummary::
    :toctree: generated/

    Interval

Properties
~~~~~~~~~~
.. autosummary::
    :toctree: generated/

    Interval.closed
    Interval.closed_left
    Interval.closed_right
    Interval.left
    Interval.length
    Interval.mid
    Interval.open_left
    Interval.open_right
    Interval.overlaps
    Interval.right

Timedelta
---------
.. autosummary::
    :toctree: generated/

    Timedelta

Properties
~~~~~~~~~~
.. autosummary::
    :toctree: generated/

    Timedelta.asm8
    Timedelta.components
    Timedelta.days
    Timedelta.delta
    Timedelta.freq
    Timedelta.is_populated
    Timedelta.max
    Timedelta.microseconds
    Timedelta.min
    Timedelta.nanoseconds
    Timedelta.resolution
    Timedelta.seconds
    Timedelta.value
    Timedelta.view

Methods
~~~~~~~
.. autosummary::
    :toctree: generated/

    Timedelta.ceil
    Timedelta.floor
    Timedelta.isoformat
    Timedelta.round
    Timedelta.to_pytimedelta
    Timedelta.to_timedelta64
    Timedelta.total_seconds
