{{ header }}

.. _api.dateoffsets:

============
Date offsets
============
.. currentmodule:: pandas.tseries.offsets

DateOffset
----------
.. autosummary::
   :toctree: api/

    DateOffset

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    DateOffset.freqstr
    DateOffset.kwds
    DateOffset.name
    DateOffset.nanos
    DateOffset.normalize
    DateOffset.rule_code
    DateOffset.n
    DateOffset.is_month_start
    DateOffset.is_month_end

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    DateOffset.copy
    DateOffset.is_anchored
    DateOffset.is_on_offset
    DateOffset.is_month_start
    DateOffset.is_month_end
    DateOffset.is_quarter_start
    DateOffset.is_quarter_end
    DateOffset.is_year_start
    DateOffset.is_year_end

BusinessDay
-----------

.. autosummary::
   :toctree: api/

    BusinessDay

Alias:

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   BDay

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessDay.freqstr
    BusinessDay.kwds
    BusinessDay.name
    BusinessDay.nanos
    BusinessDay.normalize
    BusinessDay.rule_code
    BusinessDay.n
    BusinessDay.weekmask
    BusinessDay.holidays
    BusinessDay.calendar

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessDay.copy
    BusinessDay.is_anchored
    BusinessDay.is_on_offset
    BusinessDay.is_month_start
    BusinessDay.is_month_end
    BusinessDay.is_quarter_start
    BusinessDay.is_quarter_end
    BusinessDay.is_year_start
    BusinessDay.is_year_end

BusinessHour
------------
.. autosummary::
   :toctree: api/

    BusinessHour

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessHour.freqstr
    BusinessHour.kwds
    BusinessHour.name
    BusinessHour.nanos
    BusinessHour.normalize
    BusinessHour.rule_code
    BusinessHour.n
    BusinessHour.start
    BusinessHour.end
    BusinessHour.weekmask
    BusinessHour.holidays
    BusinessHour.calendar

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessHour.copy
    BusinessHour.is_anchored
    BusinessHour.is_on_offset
    BusinessHour.is_month_start
    BusinessHour.is_month_end
    BusinessHour.is_quarter_start
    BusinessHour.is_quarter_end
    BusinessHour.is_year_start
    BusinessHour.is_year_end

CustomBusinessDay
-----------------

.. autosummary::
   :toctree: api/

    CustomBusinessDay

Alias:

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   CDay

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessDay.freqstr
    CustomBusinessDay.kwds
    CustomBusinessDay.name
    CustomBusinessDay.nanos
    CustomBusinessDay.normalize
    CustomBusinessDay.rule_code
    CustomBusinessDay.n
    CustomBusinessDay.weekmask
    CustomBusinessDay.calendar
    CustomBusinessDay.holidays

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessDay.copy
    CustomBusinessDay.is_anchored
    CustomBusinessDay.is_on_offset
    CustomBusinessDay.is_month_start
    CustomBusinessDay.is_month_end
    CustomBusinessDay.is_quarter_start
    CustomBusinessDay.is_quarter_end
    CustomBusinessDay.is_year_start
    CustomBusinessDay.is_year_end

CustomBusinessHour
------------------
.. autosummary::
   :toctree: api/

    CustomBusinessHour

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessHour.freqstr
    CustomBusinessHour.kwds
    CustomBusinessHour.name
    CustomBusinessHour.nanos
    CustomBusinessHour.normalize
    CustomBusinessHour.rule_code
    CustomBusinessHour.n
    CustomBusinessHour.weekmask
    CustomBusinessHour.calendar
    CustomBusinessHour.holidays
    CustomBusinessHour.start
    CustomBusinessHour.end

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessHour.copy
    CustomBusinessHour.is_anchored
    CustomBusinessHour.is_on_offset
    CustomBusinessHour.is_month_start
    CustomBusinessHour.is_month_end
    CustomBusinessHour.is_quarter_start
    CustomBusinessHour.is_quarter_end
    CustomBusinessHour.is_year_start
    CustomBusinessHour.is_year_end

MonthEnd
--------
.. autosummary::
   :toctree: api/

    MonthEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    MonthEnd.freqstr
    MonthEnd.kwds
    MonthEnd.name
    MonthEnd.nanos
    MonthEnd.normalize
    MonthEnd.rule_code
    MonthEnd.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    MonthEnd.copy
    MonthEnd.is_anchored
    MonthEnd.is_on_offset
    MonthEnd.is_month_start
    MonthEnd.is_month_end
    MonthEnd.is_quarter_start
    MonthEnd.is_quarter_end
    MonthEnd.is_year_start
    MonthEnd.is_year_end

MonthBegin
----------
.. autosummary::
   :toctree: api/

    MonthBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    MonthBegin.freqstr
    MonthBegin.kwds
    MonthBegin.name
    MonthBegin.nanos
    MonthBegin.normalize
    MonthBegin.rule_code
    MonthBegin.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    MonthBegin.copy
    MonthBegin.is_anchored
    MonthBegin.is_on_offset
    MonthBegin.is_month_start
    MonthBegin.is_month_end
    MonthBegin.is_quarter_start
    MonthBegin.is_quarter_end
    MonthBegin.is_year_start
    MonthBegin.is_year_end

BusinessMonthEnd
----------------

.. autosummary::
   :toctree: api/

    BusinessMonthEnd

Alias:

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   BMonthEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessMonthEnd.freqstr
    BusinessMonthEnd.kwds
    BusinessMonthEnd.name
    BusinessMonthEnd.nanos
    BusinessMonthEnd.normalize
    BusinessMonthEnd.rule_code
    BusinessMonthEnd.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessMonthEnd.copy
    BusinessMonthEnd.is_anchored
    BusinessMonthEnd.is_on_offset
    BusinessMonthEnd.is_month_start
    BusinessMonthEnd.is_month_end
    BusinessMonthEnd.is_quarter_start
    BusinessMonthEnd.is_quarter_end
    BusinessMonthEnd.is_year_start
    BusinessMonthEnd.is_year_end

BusinessMonthBegin
------------------

.. autosummary::
   :toctree: api/

    BusinessMonthBegin

Alias:

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   BMonthBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessMonthBegin.freqstr
    BusinessMonthBegin.kwds
    BusinessMonthBegin.name
    BusinessMonthBegin.nanos
    BusinessMonthBegin.normalize
    BusinessMonthBegin.rule_code
    BusinessMonthBegin.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessMonthBegin.copy
    BusinessMonthBegin.is_anchored
    BusinessMonthBegin.is_on_offset
    BusinessMonthBegin.is_month_start
    BusinessMonthBegin.is_month_end
    BusinessMonthBegin.is_quarter_start
    BusinessMonthBegin.is_quarter_end
    BusinessMonthBegin.is_year_start
    BusinessMonthBegin.is_year_end

CustomBusinessMonthEnd
----------------------

.. autosummary::
   :toctree: api/

    CustomBusinessMonthEnd

Alias:

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   CBMonthEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessMonthEnd.freqstr
    CustomBusinessMonthEnd.kwds
    CustomBusinessMonthEnd.m_offset
    CustomBusinessMonthEnd.name
    CustomBusinessMonthEnd.nanos
    CustomBusinessMonthEnd.normalize
    CustomBusinessMonthEnd.rule_code
    CustomBusinessMonthEnd.n
    CustomBusinessMonthEnd.weekmask
    CustomBusinessMonthEnd.calendar
    CustomBusinessMonthEnd.holidays

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessMonthEnd.copy
    CustomBusinessMonthEnd.is_anchored
    CustomBusinessMonthEnd.is_on_offset
    CustomBusinessMonthEnd.is_month_start
    CustomBusinessMonthEnd.is_month_end
    CustomBusinessMonthEnd.is_quarter_start
    CustomBusinessMonthEnd.is_quarter_end
    CustomBusinessMonthEnd.is_year_start
    CustomBusinessMonthEnd.is_year_end

CustomBusinessMonthBegin
------------------------

.. autosummary::
   :toctree: api/

    CustomBusinessMonthBegin

Alias:

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   CBMonthBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessMonthBegin.freqstr
    CustomBusinessMonthBegin.kwds
    CustomBusinessMonthBegin.m_offset
    CustomBusinessMonthBegin.name
    CustomBusinessMonthBegin.nanos
    CustomBusinessMonthBegin.normalize
    CustomBusinessMonthBegin.rule_code
    CustomBusinessMonthBegin.n
    CustomBusinessMonthBegin.weekmask
    CustomBusinessMonthBegin.calendar
    CustomBusinessMonthBegin.holidays

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessMonthBegin.copy
    CustomBusinessMonthBegin.is_anchored
    CustomBusinessMonthBegin.is_on_offset
    CustomBusinessMonthBegin.is_month_start
    CustomBusinessMonthBegin.is_month_end
    CustomBusinessMonthBegin.is_quarter_start
    CustomBusinessMonthBegin.is_quarter_end
    CustomBusinessMonthBegin.is_year_start
    CustomBusinessMonthBegin.is_year_end

SemiMonthEnd
------------
.. autosummary::
   :toctree: api/

    SemiMonthEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    SemiMonthEnd.freqstr
    SemiMonthEnd.kwds
    SemiMonthEnd.name
    SemiMonthEnd.nanos
    SemiMonthEnd.normalize
    SemiMonthEnd.rule_code
    SemiMonthEnd.n
    SemiMonthEnd.day_of_month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    SemiMonthEnd.copy
    SemiMonthEnd.is_anchored
    SemiMonthEnd.is_on_offset
    SemiMonthEnd.is_month_start
    SemiMonthEnd.is_month_end
    SemiMonthEnd.is_quarter_start
    SemiMonthEnd.is_quarter_end
    SemiMonthEnd.is_year_start
    SemiMonthEnd.is_year_end

SemiMonthBegin
--------------
.. autosummary::
   :toctree: api/

    SemiMonthBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    SemiMonthBegin.freqstr
    SemiMonthBegin.kwds
    SemiMonthBegin.name
    SemiMonthBegin.nanos
    SemiMonthBegin.normalize
    SemiMonthBegin.rule_code
    SemiMonthBegin.n
    SemiMonthBegin.day_of_month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    SemiMonthBegin.copy
    SemiMonthBegin.is_anchored
    SemiMonthBegin.is_on_offset
    SemiMonthBegin.is_month_start
    SemiMonthBegin.is_month_end
    SemiMonthBegin.is_quarter_start
    SemiMonthBegin.is_quarter_end
    SemiMonthBegin.is_year_start
    SemiMonthBegin.is_year_end

Week
----
.. autosummary::
   :toctree: api/

    Week

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Week.freqstr
    Week.kwds
    Week.name
    Week.nanos
    Week.normalize
    Week.rule_code
    Week.n
    Week.weekday

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Week.copy
    Week.is_anchored
    Week.is_on_offset
    Week.is_month_start
    Week.is_month_end
    Week.is_quarter_start
    Week.is_quarter_end
    Week.is_year_start
    Week.is_year_end

WeekOfMonth
-----------
.. autosummary::
   :toctree: api/

    WeekOfMonth

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    WeekOfMonth.freqstr
    WeekOfMonth.kwds
    WeekOfMonth.name
    WeekOfMonth.nanos
    WeekOfMonth.normalize
    WeekOfMonth.rule_code
    WeekOfMonth.n
    WeekOfMonth.week

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    WeekOfMonth.copy
    WeekOfMonth.is_anchored
    WeekOfMonth.is_on_offset
    WeekOfMonth.weekday
    WeekOfMonth.is_month_start
    WeekOfMonth.is_month_end
    WeekOfMonth.is_quarter_start
    WeekOfMonth.is_quarter_end
    WeekOfMonth.is_year_start
    WeekOfMonth.is_year_end

LastWeekOfMonth
---------------
.. autosummary::
   :toctree: api/

    LastWeekOfMonth

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    LastWeekOfMonth.freqstr
    LastWeekOfMonth.kwds
    LastWeekOfMonth.name
    LastWeekOfMonth.nanos
    LastWeekOfMonth.normalize
    LastWeekOfMonth.rule_code
    LastWeekOfMonth.n
    LastWeekOfMonth.weekday
    LastWeekOfMonth.week

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    LastWeekOfMonth.copy
    LastWeekOfMonth.is_anchored
    LastWeekOfMonth.is_on_offset
    LastWeekOfMonth.is_month_start
    LastWeekOfMonth.is_month_end
    LastWeekOfMonth.is_quarter_start
    LastWeekOfMonth.is_quarter_end
    LastWeekOfMonth.is_year_start
    LastWeekOfMonth.is_year_end

BQuarterEnd
-----------
.. autosummary::
   :toctree: api/

    BQuarterEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BQuarterEnd.freqstr
    BQuarterEnd.kwds
    BQuarterEnd.name
    BQuarterEnd.nanos
    BQuarterEnd.normalize
    BQuarterEnd.rule_code
    BQuarterEnd.n
    BQuarterEnd.startingMonth

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BQuarterEnd.copy
    BQuarterEnd.is_anchored
    BQuarterEnd.is_on_offset
    BQuarterEnd.is_month_start
    BQuarterEnd.is_month_end
    BQuarterEnd.is_quarter_start
    BQuarterEnd.is_quarter_end
    BQuarterEnd.is_year_start
    BQuarterEnd.is_year_end

BQuarterBegin
-------------
.. autosummary::
   :toctree: api/

    BQuarterBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BQuarterBegin.freqstr
    BQuarterBegin.kwds
    BQuarterBegin.name
    BQuarterBegin.nanos
    BQuarterBegin.normalize
    BQuarterBegin.rule_code
    BQuarterBegin.n
    BQuarterBegin.startingMonth

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BQuarterBegin.copy
    BQuarterBegin.is_anchored
    BQuarterBegin.is_on_offset
    BQuarterBegin.is_month_start
    BQuarterBegin.is_month_end
    BQuarterBegin.is_quarter_start
    BQuarterBegin.is_quarter_end
    BQuarterBegin.is_year_start
    BQuarterBegin.is_year_end

QuarterEnd
----------
.. autosummary::
   :toctree: api/

    QuarterEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    QuarterEnd.freqstr
    QuarterEnd.kwds
    QuarterEnd.name
    QuarterEnd.nanos
    QuarterEnd.normalize
    QuarterEnd.rule_code
    QuarterEnd.n
    QuarterEnd.startingMonth

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    QuarterEnd.copy
    QuarterEnd.is_anchored
    QuarterEnd.is_on_offset
    QuarterEnd.is_month_start
    QuarterEnd.is_month_end
    QuarterEnd.is_quarter_start
    QuarterEnd.is_quarter_end
    QuarterEnd.is_year_start
    QuarterEnd.is_year_end

QuarterBegin
------------
.. autosummary::
   :toctree: api/

    QuarterBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    QuarterBegin.freqstr
    QuarterBegin.kwds
    QuarterBegin.name
    QuarterBegin.nanos
    QuarterBegin.normalize
    QuarterBegin.rule_code
    QuarterBegin.n
    QuarterBegin.startingMonth

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    QuarterBegin.copy
    QuarterBegin.is_anchored
    QuarterBegin.is_on_offset
    QuarterBegin.is_month_start
    QuarterBegin.is_month_end
    QuarterBegin.is_quarter_start
    QuarterBegin.is_quarter_end
    QuarterBegin.is_year_start
    QuarterBegin.is_year_end

BYearEnd
--------
.. autosummary::
   :toctree: api/

    BYearEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BYearEnd.freqstr
    BYearEnd.kwds
    BYearEnd.name
    BYearEnd.nanos
    BYearEnd.normalize
    BYearEnd.rule_code
    BYearEnd.n
    BYearEnd.month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BYearEnd.copy
    BYearEnd.is_anchored
    BYearEnd.is_on_offset
    BYearEnd.is_month_start
    BYearEnd.is_month_end
    BYearEnd.is_quarter_start
    BYearEnd.is_quarter_end
    BYearEnd.is_year_start
    BYearEnd.is_year_end

BYearBegin
----------
.. autosummary::
   :toctree: api/

    BYearBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BYearBegin.freqstr
    BYearBegin.kwds
    BYearBegin.name
    BYearBegin.nanos
    BYearBegin.normalize
    BYearBegin.rule_code
    BYearBegin.n
    BYearBegin.month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BYearBegin.copy
    BYearBegin.is_anchored
    BYearBegin.is_on_offset
    BYearBegin.is_month_start
    BYearBegin.is_month_end
    BYearBegin.is_quarter_start
    BYearBegin.is_quarter_end
    BYearBegin.is_year_start
    BYearBegin.is_year_end

YearEnd
-------
.. autosummary::
   :toctree: api/

    YearEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    YearEnd.freqstr
    YearEnd.kwds
    YearEnd.name
    YearEnd.nanos
    YearEnd.normalize
    YearEnd.rule_code
    YearEnd.n
    YearEnd.month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    YearEnd.copy
    YearEnd.is_anchored
    YearEnd.is_on_offset
    YearEnd.is_month_start
    YearEnd.is_month_end
    YearEnd.is_quarter_start
    YearEnd.is_quarter_end
    YearEnd.is_year_start
    YearEnd.is_year_end

YearBegin
---------
.. autosummary::
   :toctree: api/

    YearBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    YearBegin.freqstr
    YearBegin.kwds
    YearBegin.name
    YearBegin.nanos
    YearBegin.normalize
    YearBegin.rule_code
    YearBegin.n
    YearBegin.month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    YearBegin.copy
    YearBegin.is_anchored
    YearBegin.is_on_offset
    YearBegin.is_month_start
    YearBegin.is_month_end
    YearBegin.is_quarter_start
    YearBegin.is_quarter_end
    YearBegin.is_year_start
    YearBegin.is_year_end

FY5253
------
.. autosummary::
   :toctree: api/

    FY5253

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    FY5253.freqstr
    FY5253.kwds
    FY5253.name
    FY5253.nanos
    FY5253.normalize
    FY5253.rule_code
    FY5253.n
    FY5253.startingMonth
    FY5253.variation
    FY5253.weekday

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    FY5253.copy
    FY5253.get_rule_code_suffix
    FY5253.get_year_end
    FY5253.is_anchored
    FY5253.is_on_offset
    FY5253.is_month_start
    FY5253.is_month_end
    FY5253.is_quarter_start
    FY5253.is_quarter_end
    FY5253.is_year_start
    FY5253.is_year_end

FY5253Quarter
-------------
.. autosummary::
   :toctree: api/

    FY5253Quarter

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    FY5253Quarter.freqstr
    FY5253Quarter.kwds
    FY5253Quarter.name
    FY5253Quarter.nanos
    FY5253Quarter.normalize
    FY5253Quarter.rule_code
    FY5253Quarter.n
    FY5253Quarter.qtr_with_extra_week
    FY5253Quarter.startingMonth
    FY5253Quarter.variation
    FY5253Quarter.weekday

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    FY5253Quarter.copy
    FY5253Quarter.get_rule_code_suffix
    FY5253Quarter.get_weeks
    FY5253Quarter.is_anchored
    FY5253Quarter.is_on_offset
    FY5253Quarter.year_has_extra_week
    FY5253Quarter.is_month_start
    FY5253Quarter.is_month_end
    FY5253Quarter.is_quarter_start
    FY5253Quarter.is_quarter_end
    FY5253Quarter.is_year_start
    FY5253Quarter.is_year_end

Easter
------
.. autosummary::
   :toctree: api/

    Easter

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Easter.freqstr
    Easter.kwds
    Easter.name
    Easter.nanos
    Easter.normalize
    Easter.rule_code
    Easter.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Easter.copy
    Easter.is_anchored
    Easter.is_on_offset
    Easter.is_month_start
    Easter.is_month_end
    Easter.is_quarter_start
    Easter.is_quarter_end
    Easter.is_year_start
    Easter.is_year_end

Tick
----
.. autosummary::
   :toctree: api/

    Tick

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Tick.delta
    Tick.freqstr
    Tick.kwds
    Tick.name
    Tick.nanos
    Tick.normalize
    Tick.rule_code
    Tick.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Tick.copy
    Tick.is_anchored
    Tick.is_on_offset
    Tick.is_month_start
    Tick.is_month_end
    Tick.is_quarter_start
    Tick.is_quarter_end
    Tick.is_year_start
    Tick.is_year_end

Day
---
.. autosummary::
   :toctree: api/

    Day

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Day.delta
    Day.freqstr
    Day.kwds
    Day.name
    Day.nanos
    Day.normalize
    Day.rule_code
    Day.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Day.copy
    Day.is_anchored
    Day.is_on_offset
    Day.is_month_start
    Day.is_month_end
    Day.is_quarter_start
    Day.is_quarter_end
    Day.is_year_start
    Day.is_year_end

Hour
----
.. autosummary::
   :toctree: api/

    Hour

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Hour.delta
    Hour.freqstr
    Hour.kwds
    Hour.name
    Hour.nanos
    Hour.normalize
    Hour.rule_code
    Hour.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Hour.copy
    Hour.is_anchored
    Hour.is_on_offset
    Hour.is_month_start
    Hour.is_month_end
    Hour.is_quarter_start
    Hour.is_quarter_end
    Hour.is_year_start
    Hour.is_year_end

Minute
------
.. autosummary::
   :toctree: api/

    Minute

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Minute.delta
    Minute.freqstr
    Minute.kwds
    Minute.name
    Minute.nanos
    Minute.normalize
    Minute.rule_code
    Minute.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Minute.copy
    Minute.is_anchored
    Minute.is_on_offset
    Minute.is_month_start
    Minute.is_month_end
    Minute.is_quarter_start
    Minute.is_quarter_end
    Minute.is_year_start
    Minute.is_year_end

Second
------
.. autosummary::
   :toctree: api/

    Second

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Second.delta
    Second.freqstr
    Second.kwds
    Second.name
    Second.nanos
    Second.normalize
    Second.rule_code
    Second.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Second.copy
    Second.is_anchored
    Second.is_on_offset
    Second.is_month_start
    Second.is_month_end
    Second.is_quarter_start
    Second.is_quarter_end
    Second.is_year_start
    Second.is_year_end

Milli
-----
.. autosummary::
   :toctree: api/

    Milli

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Milli.delta
    Milli.freqstr
    Milli.kwds
    Milli.name
    Milli.nanos
    Milli.normalize
    Milli.rule_code
    Milli.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Milli.copy
    Milli.is_anchored
    Milli.is_on_offset
    Milli.is_month_start
    Milli.is_month_end
    Milli.is_quarter_start
    Milli.is_quarter_end
    Milli.is_year_start
    Milli.is_year_end

Micro
-----
.. autosummary::
   :toctree: api/

    Micro

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Micro.delta
    Micro.freqstr
    Micro.kwds
    Micro.name
    Micro.nanos
    Micro.normalize
    Micro.rule_code
    Micro.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Micro.copy
    Micro.is_anchored
    Micro.is_on_offset
    Micro.is_month_start
    Micro.is_month_end
    Micro.is_quarter_start
    Micro.is_quarter_end
    Micro.is_year_start
    Micro.is_year_end

Nano
----
.. autosummary::
   :toctree: api/

    Nano

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Nano.delta
    Nano.freqstr
    Nano.kwds
    Nano.name
    Nano.nanos
    Nano.normalize
    Nano.rule_code
    Nano.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Nano.copy
    Nano.is_anchored
    Nano.is_on_offset
    Nano.is_month_start
    Nano.is_month_end
    Nano.is_quarter_start
    Nano.is_quarter_end
    Nano.is_year_start
    Nano.is_year_end

.. _api.frequencies:

===========
Frequencies
===========
.. currentmodule:: pandas.tseries.frequencies

.. _api.offsets:

.. autosummary::
   :toctree: api/

   to_offset
