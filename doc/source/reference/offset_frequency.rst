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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    DateOffset.apply
    DateOffset.copy
    DateOffset.isAnchored
    DateOffset.onOffset
    DateOffset.is_anchored
    DateOffset.is_on_offset

BusinessDay
-----------
.. autosummary::
   :toctree: api/

    BusinessDay

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessDay.apply
    BusinessDay.apply_index
    BusinessDay.copy
    BusinessDay.isAnchored
    BusinessDay.onOffset
    BusinessDay.is_anchored
    BusinessDay.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessHour.apply
    BusinessHour.copy
    BusinessHour.isAnchored
    BusinessHour.onOffset
    BusinessHour.is_anchored
    BusinessHour.is_on_offset

CustomBusinessDay
-----------------
.. autosummary::
   :toctree: api/

    CustomBusinessDay

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessDay.apply
    CustomBusinessDay.copy
    CustomBusinessDay.isAnchored
    CustomBusinessDay.onOffset
    CustomBusinessDay.is_anchored
    CustomBusinessDay.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessHour.apply
    CustomBusinessHour.copy
    CustomBusinessHour.isAnchored
    CustomBusinessHour.onOffset
    CustomBusinessHour.is_anchored
    CustomBusinessHour.is_on_offset

MonthOffset
-----------
.. autosummary::
   :toctree: api/

    MonthOffset

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    MonthOffset.freqstr
    MonthOffset.kwds
    MonthOffset.name
    MonthOffset.nanos
    MonthOffset.normalize
    MonthOffset.rule_code

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    MonthOffset.apply
    MonthOffset.apply_index
    MonthOffset.copy
    MonthOffset.isAnchored
    MonthOffset.onOffset
    MonthOffset.is_anchored
    MonthOffset.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    MonthEnd.apply
    MonthEnd.apply_index
    MonthEnd.copy
    MonthEnd.isAnchored
    MonthEnd.onOffset
    MonthEnd.is_anchored
    MonthEnd.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    MonthBegin.apply
    MonthBegin.apply_index
    MonthBegin.copy
    MonthBegin.isAnchored
    MonthBegin.onOffset
    MonthBegin.is_anchored
    MonthBegin.is_on_offset

BusinessMonthEnd
----------------
.. autosummary::
   :toctree: api/

    BusinessMonthEnd

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessMonthEnd.apply
    BusinessMonthEnd.apply_index
    BusinessMonthEnd.copy
    BusinessMonthEnd.isAnchored
    BusinessMonthEnd.onOffset
    BusinessMonthEnd.is_anchored
    BusinessMonthEnd.is_on_offset

BusinessMonthBegin
------------------
.. autosummary::
   :toctree: api/

    BusinessMonthBegin

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessMonthBegin.apply
    BusinessMonthBegin.apply_index
    BusinessMonthBegin.copy
    BusinessMonthBegin.isAnchored
    BusinessMonthBegin.onOffset
    BusinessMonthBegin.is_anchored
    BusinessMonthBegin.is_on_offset

CustomBusinessMonthEnd
----------------------
.. autosummary::
   :toctree: api/

    CustomBusinessMonthEnd

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessMonthEnd.apply
    CustomBusinessMonthEnd.copy
    CustomBusinessMonthEnd.isAnchored
    CustomBusinessMonthEnd.onOffset
    CustomBusinessMonthEnd.is_anchored
    CustomBusinessMonthEnd.is_on_offset

CustomBusinessMonthBegin
------------------------
.. autosummary::
   :toctree: api/

    CustomBusinessMonthBegin

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessMonthBegin.apply
    CustomBusinessMonthBegin.copy
    CustomBusinessMonthBegin.isAnchored
    CustomBusinessMonthBegin.onOffset
    CustomBusinessMonthBegin.is_anchored
    CustomBusinessMonthBegin.is_on_offset

SemiMonthOffset
---------------
.. autosummary::
   :toctree: api/

    SemiMonthOffset

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    SemiMonthOffset.freqstr
    SemiMonthOffset.kwds
    SemiMonthOffset.name
    SemiMonthOffset.nanos
    SemiMonthOffset.normalize
    SemiMonthOffset.rule_code

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    SemiMonthOffset.apply
    SemiMonthOffset.apply_index
    SemiMonthOffset.copy
    SemiMonthOffset.isAnchored
    SemiMonthOffset.onOffset
    SemiMonthOffset.is_anchored
    SemiMonthOffset.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    SemiMonthEnd.apply
    SemiMonthEnd.apply_index
    SemiMonthEnd.copy
    SemiMonthEnd.isAnchored
    SemiMonthEnd.onOffset
    SemiMonthEnd.is_anchored
    SemiMonthEnd.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    SemiMonthBegin.apply
    SemiMonthBegin.apply_index
    SemiMonthBegin.copy
    SemiMonthBegin.isAnchored
    SemiMonthBegin.onOffset
    SemiMonthBegin.is_anchored
    SemiMonthBegin.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Week.apply
    Week.apply_index
    Week.copy
    Week.isAnchored
    Week.onOffset
    Week.is_anchored
    Week.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    WeekOfMonth.apply
    WeekOfMonth.copy
    WeekOfMonth.isAnchored
    WeekOfMonth.onOffset
    WeekOfMonth.is_anchored
    WeekOfMonth.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    LastWeekOfMonth.apply
    LastWeekOfMonth.copy
    LastWeekOfMonth.isAnchored
    LastWeekOfMonth.onOffset
    LastWeekOfMonth.is_anchored
    LastWeekOfMonth.is_on_offset

QuarterOffset
-------------
.. autosummary::
   :toctree: api/

    QuarterOffset

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    QuarterOffset.freqstr
    QuarterOffset.kwds
    QuarterOffset.name
    QuarterOffset.nanos
    QuarterOffset.normalize
    QuarterOffset.rule_code

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    QuarterOffset.apply
    QuarterOffset.apply_index
    QuarterOffset.copy
    QuarterOffset.isAnchored
    QuarterOffset.onOffset
    QuarterOffset.is_anchored
    QuarterOffset.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BQuarterEnd.apply
    BQuarterEnd.apply_index
    BQuarterEnd.copy
    BQuarterEnd.isAnchored
    BQuarterEnd.onOffset
    BQuarterEnd.is_anchored
    BQuarterEnd.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BQuarterBegin.apply
    BQuarterBegin.apply_index
    BQuarterBegin.copy
    BQuarterBegin.isAnchored
    BQuarterBegin.onOffset
    BQuarterBegin.is_anchored
    BQuarterBegin.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    QuarterEnd.apply
    QuarterEnd.apply_index
    QuarterEnd.copy
    QuarterEnd.isAnchored
    QuarterEnd.onOffset
    QuarterEnd.is_anchored
    QuarterEnd.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    QuarterBegin.apply
    QuarterBegin.apply_index
    QuarterBegin.copy
    QuarterBegin.isAnchored
    QuarterBegin.onOffset
    QuarterBegin.is_anchored
    QuarterBegin.is_on_offset

YearOffset
----------
.. autosummary::
   :toctree: api/

    YearOffset

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    YearOffset.freqstr
    YearOffset.kwds
    YearOffset.name
    YearOffset.nanos
    YearOffset.normalize
    YearOffset.rule_code

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    YearOffset.apply
    YearOffset.apply_index
    YearOffset.copy
    YearOffset.isAnchored
    YearOffset.onOffset
    YearOffset.is_anchored
    YearOffset.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BYearEnd.apply
    BYearEnd.apply_index
    BYearEnd.copy
    BYearEnd.isAnchored
    BYearEnd.onOffset
    BYearEnd.is_anchored
    BYearEnd.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BYearBegin.apply
    BYearBegin.apply_index
    BYearBegin.copy
    BYearBegin.isAnchored
    BYearBegin.onOffset
    BYearBegin.is_anchored
    BYearBegin.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    YearEnd.apply
    YearEnd.apply_index
    YearEnd.copy
    YearEnd.isAnchored
    YearEnd.onOffset
    YearEnd.is_anchored
    YearEnd.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    YearBegin.apply
    YearBegin.apply_index
    YearBegin.copy
    YearBegin.isAnchored
    YearBegin.onOffset
    YearBegin.is_anchored
    YearBegin.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    FY5253.apply
    FY5253.copy
    FY5253.get_rule_code_suffix
    FY5253.get_year_end
    FY5253.isAnchored
    FY5253.onOffset
    FY5253.is_anchored
    FY5253.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    FY5253Quarter.apply
    FY5253Quarter.copy
    FY5253Quarter.get_weeks
    FY5253Quarter.isAnchored
    FY5253Quarter.onOffset
    FY5253Quarter.is_anchored
    FY5253Quarter.is_on_offset
    FY5253Quarter.year_has_extra_week

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Easter.apply
    Easter.copy
    Easter.isAnchored
    Easter.onOffset
    Easter.is_anchored
    Easter.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Tick.copy
    Tick.isAnchored
    Tick.onOffset
    Tick.is_anchored
    Tick.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Day.copy
    Day.isAnchored
    Day.onOffset
    Day.is_anchored
    Day.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Hour.copy
    Hour.isAnchored
    Hour.onOffset
    Hour.is_anchored
    Hour.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Minute.copy
    Minute.isAnchored
    Minute.onOffset
    Minute.is_anchored
    Minute.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Second.copy
    Second.isAnchored
    Second.onOffset
    Second.is_anchored
    Second.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Milli.copy
    Milli.isAnchored
    Milli.onOffset
    Milli.is_anchored
    Milli.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Micro.copy
    Micro.isAnchored
    Micro.onOffset
    Micro.is_anchored
    Micro.is_on_offset

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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    Nano.copy
    Nano.isAnchored
    Nano.onOffset
    Nano.is_anchored
    Nano.is_on_offset

BDay
----
.. autosummary::
   :toctree: api/

    BDay

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BDay.base
    BDay.freqstr
    BDay.kwds
    BDay.name
    BDay.nanos
    BDay.normalize
    BDay.offset
    BDay.rule_code

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BDay.apply
    BDay.apply_index
    BDay.copy
    BDay.isAnchored
    BDay.onOffset
    BDay.is_anchored
    BDay.is_on_offset
    BDay.rollback
    BDay.rollforward

BMonthEnd
---------
.. autosummary::
   :toctree: api/

    BMonthEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BMonthEnd.base
    BMonthEnd.freqstr
    BMonthEnd.kwds
    BMonthEnd.name
    BMonthEnd.nanos
    BMonthEnd.normalize
    BMonthEnd.rule_code

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BMonthEnd.apply
    BMonthEnd.apply_index
    BMonthEnd.copy
    BMonthEnd.isAnchored
    BMonthEnd.onOffset
    BMonthEnd.is_anchored
    BMonthEnd.is_on_offset
    BMonthEnd.rollback
    BMonthEnd.rollforward

BMonthBegin
-----------
.. autosummary::
   :toctree: api/

    BMonthBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BMonthBegin.base
    BMonthBegin.freqstr
    BMonthBegin.kwds
    BMonthBegin.name
    BMonthBegin.nanos
    BMonthBegin.normalize
    BMonthBegin.rule_code

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BMonthBegin.apply
    BMonthBegin.apply_index
    BMonthBegin.copy
    BMonthBegin.isAnchored
    BMonthBegin.onOffset
    BMonthBegin.is_anchored
    BMonthBegin.is_on_offset
    BMonthBegin.rollback
    BMonthBegin.rollforward

CBMonthEnd
----------
.. autosummary::
   :toctree: api/

    CBMonthEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    CBMonthEnd.base
    CBMonthEnd.cbday_roll
    CBMonthEnd.freqstr
    CBMonthEnd.kwds
    CBMonthEnd.m_offset
    CBMonthEnd.month_roll
    CBMonthEnd.name
    CBMonthEnd.nanos
    CBMonthEnd.normalize
    CBMonthEnd.offset
    CBMonthEnd.rule_code

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CBMonthEnd.apply
    CBMonthEnd.apply_index
    CBMonthEnd.copy
    CBMonthEnd.isAnchored
    CBMonthEnd.onOffset
    CBMonthEnd.is_anchored
    CBMonthEnd.is_on_offset
    CBMonthEnd.rollback
    CBMonthEnd.rollforward

CBMonthBegin
------------
.. autosummary::
   :toctree: api/

    CBMonthBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    CBMonthBegin.base
    CBMonthBegin.cbday_roll
    CBMonthBegin.freqstr
    CBMonthBegin.kwds
    CBMonthBegin.m_offset
    CBMonthBegin.month_roll
    CBMonthBegin.name
    CBMonthBegin.nanos
    CBMonthBegin.normalize
    CBMonthBegin.offset
    CBMonthBegin.rule_code

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CBMonthBegin.apply
    CBMonthBegin.apply_index
    CBMonthBegin.copy
    CBMonthBegin.isAnchored
    CBMonthBegin.onOffset
    CBMonthBegin.is_anchored
    CBMonthBegin.is_on_offset
    CBMonthBegin.rollback
    CBMonthBegin.rollforward

CDay
----
.. autosummary::
   :toctree: api/

    CDay

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    CDay.base
    CDay.freqstr
    CDay.kwds
    CDay.name
    CDay.nanos
    CDay.normalize
    CDay.offset
    CDay.rule_code

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    CDay.apply
    CDay.apply_index
    CDay.copy
    CDay.isAnchored
    CDay.onOffset
    CDay.is_anchored
    CDay.is_on_offset
    CDay.rollback
    CDay.rollforward

.. _api.frequencies:

===========
Frequencies
===========
.. currentmodule:: pandas.tseries.frequencies

.. _api.offsets:

.. autosummary::
   :toctree: api/

   to_offset
