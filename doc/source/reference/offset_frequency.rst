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
    CDay.is_anchored
    CDay.is_on_offset
    CDay.rollback
    CDay.rollforward

.. This is to prevent warnings in the doc build. We don't want to encourage
.. these methods.

..
    .. toctree::

        ..api/pandas.tseries.offsets.DateOffset.isAnchored
        ..api/pandas.tseries.offsets.BusinessDay.isAnchored
        ..api/pandas.tseries.offsets.CDay.isAnchored
        ..api/pandas.tseries.offsets.CBMonthBegin.isAnchored
        ..api/pandas.tseries.offsets.CBMonthEnd.isAnchored
        ..api/pandas.tseries.offsets.BMonthBegin.isAnchored
        ..api/pandas.tseries.offsets.BMonthEnd.isAnchored
        ..api/pandas.tseries.offsets.BDay.isAnchored
        ..api/pandas.tseries.offsets.Nano.isAnchored
        ..api/pandas.tseries.offsets.Micro.isAnchored
        ..api/pandas.tseries.offsets.Milli.isAnchored
        ..api/pandas.tseries.offsets.Second.isAnchored
        ..api/pandas.tseries.offsets.Minute.isAnchored
        ..api/pandas.tseries.offsets.Hour.isAnchored
        ..api/pandas.tseries.offsets.Day.isAnchored
        ..api/pandas.tseries.offsets.Tick.isAnchored
        ..api/pandas.tseries.offsets.Easter.isAnchored
        ..api/pandas.tseries.offsets.FY5253Quarter.isAnchored
        ..api/pandas.tseries.offsets.FY5253.isAnchored
        ..api/pandas.tseries.offsets.YearBegin.isAnchored
        ..api/pandas.tseries.offsets.YearEnd.isAnchored
        ..api/pandas.tseries.offsets.BYearBegin.isAnchored
        ..api/pandas.tseries.offsets.BYearEnd.isAnchored
        ..api/pandas.tseries.offsets.YearOffset.isAnchored
        ..api/pandas.tseries.offsets.QuarterBegin.isAnchored
        ..api/pandas.tseries.offsets.QuarterEnd.isAnchored
        ..api/pandas.tseries.offsets.BQuarterBegin.isAnchored
        ..api/pandas.tseries.offsets.BQuarterEnd.isAnchored
        ..api/pandas.tseries.offsets.QuarterOffset.isAnchored
        ..api/pandas.tseries.offsets.LastWeekOfMonth.isAnchored
        ..api/pandas.tseries.offsets.WeekOfMonth.isAnchored
        ..api/pandas.tseries.offsets.Week.isAnchored
        ..api/pandas.tseries.offsets.SemiMonthBegin.isAnchored
        ..api/pandas.tseries.offsets.SemiMonthEnd.isAnchored
        ..api/pandas.tseries.offsets.SemiMonthOffset.isAnchored
        ..api/pandas.tseries.offsets.CustomBusinessMonthBegin.isAnchored
        ..api/pandas.tseries.offsets.CustomBusinessMonthEnd.isAnchored
        ..api/pandas.tseries.offsets.BusinessMonthBegin.isAnchored
        ..api/pandas.tseries.offsets.BusinessMonthEnd.isAnchored
        ..api/pandas.tseries.offsets.MonthBegin.isAnchored
        ..api/pandas.tseries.offsets.MonthEnd.isAnchored
        ..api/pandas.tseries.offsets.MonthOffset.isAnchored
        ..api/pandas.tseries.offsets.CustomBusinessHour.isAnchored
        ..api/pandas.tseries.offsets.CustomBusinessDay.isAnchored
        ..api/pandas.tseries.offsets.BusinessHour.isAnchored
        ..api/pandas.tseries.offsets.DateOffset.onOffset
        ..api/pandas.tseries.offsets.BusinessDay.onOffset
        ..api/pandas.tseries.offsets.CDay.onOffset
        ..api/pandas.tseries.offsets.CBMonthBegin.onOffset
        ..api/pandas.tseries.offsets.CBMonthEnd.onOffset
        ..api/pandas.tseries.offsets.BMonthBegin.onOffset
        ..api/pandas.tseries.offsets.BMonthEnd.onOffset
        ..api/pandas.tseries.offsets.BDay.onOffset
        ..api/pandas.tseries.offsets.Nano.onOffset
        ..api/pandas.tseries.offsets.Micro.onOffset
        ..api/pandas.tseries.offsets.Milli.onOffset
        ..api/pandas.tseries.offsets.Second.onOffset
        ..api/pandas.tseries.offsets.Minute.onOffset
        ..api/pandas.tseries.offsets.Hour.onOffset
        ..api/pandas.tseries.offsets.Day.onOffset
        ..api/pandas.tseries.offsets.Tick.onOffset
        ..api/pandas.tseries.offsets.Easter.onOffset
        ..api/pandas.tseries.offsets.FY5253Quarter.onOffset
        ..api/pandas.tseries.offsets.FY5253.onOffset
        ..api/pandas.tseries.offsets.YearBegin.onOffset
        ..api/pandas.tseries.offsets.YearEnd.onOffset
        ..api/pandas.tseries.offsets.BYearBegin.onOffset
        ..api/pandas.tseries.offsets.BYearEnd.onOffset
        ..api/pandas.tseries.offsets.YearOffset.onOffset
        ..api/pandas.tseries.offsets.QuarterBegin.onOffset
        ..api/pandas.tseries.offsets.QuarterEnd.onOffset
        ..api/pandas.tseries.offsets.BQuarterBegin.onOffset
        ..api/pandas.tseries.offsets.BQuarterEnd.onOffset
        ..api/pandas.tseries.offsets.QuarterOffset.onOffset
        ..api/pandas.tseries.offsets.LastWeekOfMonth.onOffset
        ..api/pandas.tseries.offsets.WeekOfMonth.onOffset
        ..api/pandas.tseries.offsets.Week.onOffset
        ..api/pandas.tseries.offsets.SemiMonthBegin.onOffset
        ..api/pandas.tseries.offsets.SemiMonthEnd.onOffset
        ..api/pandas.tseries.offsets.SemiMonthOffset.onOffset
        ..api/pandas.tseries.offsets.CustomBusinessMonthBegin.onOffset
        ..api/pandas.tseries.offsets.CustomBusinessMonthEnd.onOffset
        ..api/pandas.tseries.offsets.BusinessMonthBegin.onOffset
        ..api/pandas.tseries.offsets.BusinessMonthEnd.onOffset
        ..api/pandas.tseries.offsets.MonthBegin.onOffset
        ..api/pandas.tseries.offsets.MonthEnd.onOffset
        ..api/pandas.tseries.offsets.MonthOffset.onOffset
        ..api/pandas.tseries.offsets.CustomBusinessHour.onOffset
        ..api/pandas.tseries.offsets.CustomBusinessDay.onOffset
        ..api/pandas.tseries.offsets.BusinessHour.onOffset
        api/pandas.tseries.offsets.DateOffset.isAnchored
        api/pandas.tseries.offsets.BusinessDay.isAnchored
        api/pandas.tseries.offsets.CDay.isAnchored
        api/pandas.tseries.offsets.CBMonthBegin.isAnchored
        api/pandas.tseries.offsets.CBMonthEnd.isAnchored
        api/pandas.tseries.offsets.BMonthBegin.isAnchored
        api/pandas.tseries.offsets.BMonthEnd.isAnchored
        api/pandas.tseries.offsets.BDay.isAnchored
        api/pandas.tseries.offsets.Nano.isAnchored
        api/pandas.tseries.offsets.Micro.isAnchored
        api/pandas.tseries.offsets.Milli.isAnchored
        api/pandas.tseries.offsets.Second.isAnchored
        api/pandas.tseries.offsets.Minute.isAnchored
        api/pandas.tseries.offsets.Hour.isAnchored
        api/pandas.tseries.offsets.Day.isAnchored
        api/pandas.tseries.offsets.Tick.isAnchored
        api/pandas.tseries.offsets.Easter.isAnchored
        api/pandas.tseries.offsets.FY5253Quarter.isAnchored
        api/pandas.tseries.offsets.FY5253.isAnchored
        api/pandas.tseries.offsets.YearBegin.isAnchored
        api/pandas.tseries.offsets.YearEnd.isAnchored
        api/pandas.tseries.offsets.BYearBegin.isAnchored
        api/pandas.tseries.offsets.BYearEnd.isAnchored
        api/pandas.tseries.offsets.YearOffset.isAnchored
        api/pandas.tseries.offsets.QuarterBegin.isAnchored
        api/pandas.tseries.offsets.QuarterEnd.isAnchored
        api/pandas.tseries.offsets.BQuarterBegin.isAnchored
        api/pandas.tseries.offsets.BQuarterEnd.isAnchored
        api/pandas.tseries.offsets.QuarterOffset.isAnchored
        api/pandas.tseries.offsets.LastWeekOfMonth.isAnchored
        api/pandas.tseries.offsets.WeekOfMonth.isAnchored
        api/pandas.tseries.offsets.Week.isAnchored
        api/pandas.tseries.offsets.SemiMonthBegin.isAnchored
        api/pandas.tseries.offsets.SemiMonthEnd.isAnchored
        api/pandas.tseries.offsets.SemiMonthOffset.isAnchored
        api/pandas.tseries.offsets.CustomBusinessMonthBegin.isAnchored
        api/pandas.tseries.offsets.CustomBusinessMonthEnd.isAnchored
        api/pandas.tseries.offsets.BusinessMonthBegin.isAnchored
        api/pandas.tseries.offsets.BusinessMonthEnd.isAnchored
        api/pandas.tseries.offsets.MonthBegin.isAnchored
        api/pandas.tseries.offsets.MonthEnd.isAnchored
        api/pandas.tseries.offsets.MonthOffset.isAnchored
        api/pandas.tseries.offsets.CustomBusinessHour.isAnchored
        api/pandas.tseries.offsets.CustomBusinessDay.isAnchored
        api/pandas.tseries.offsets.BusinessHour.isAnchored
        api/pandas.tseries.offsets.DateOffset.onOffset
        api/pandas.tseries.offsets.BusinessDay.onOffset
        api/pandas.tseries.offsets.CDay.onOffset
        api/pandas.tseries.offsets.CBMonthBegin.onOffset
        api/pandas.tseries.offsets.CBMonthEnd.onOffset
        api/pandas.tseries.offsets.BMonthBegin.onOffset
        api/pandas.tseries.offsets.BMonthEnd.onOffset
        api/pandas.tseries.offsets.BDay.onOffset
        api/pandas.tseries.offsets.Nano.onOffset
        api/pandas.tseries.offsets.Micro.onOffset
        api/pandas.tseries.offsets.Milli.onOffset
        api/pandas.tseries.offsets.Second.onOffset
        api/pandas.tseries.offsets.Minute.onOffset
        api/pandas.tseries.offsets.Hour.onOffset
        api/pandas.tseries.offsets.Day.onOffset
        api/pandas.tseries.offsets.Tick.onOffset
        api/pandas.tseries.offsets.Easter.onOffset
        api/pandas.tseries.offsets.FY5253Quarter.onOffset
        api/pandas.tseries.offsets.FY5253.onOffset
        api/pandas.tseries.offsets.YearBegin.onOffset
        api/pandas.tseries.offsets.YearEnd.onOffset
        api/pandas.tseries.offsets.BYearBegin.onOffset
        api/pandas.tseries.offsets.BYearEnd.onOffset
        api/pandas.tseries.offsets.YearOffset.onOffset
        api/pandas.tseries.offsets.QuarterBegin.onOffset
        api/pandas.tseries.offsets.QuarterEnd.onOffset
        api/pandas.tseries.offsets.BQuarterBegin.onOffset
        api/pandas.tseries.offsets.BQuarterEnd.onOffset
        api/pandas.tseries.offsets.QuarterOffset.onOffset
        api/pandas.tseries.offsets.LastWeekOfMonth.onOffset
        api/pandas.tseries.offsets.WeekOfMonth.onOffset
        api/pandas.tseries.offsets.Week.onOffset
        api/pandas.tseries.offsets.SemiMonthBegin.onOffset
        api/pandas.tseries.offsets.SemiMonthEnd.onOffset
        api/pandas.tseries.offsets.SemiMonthOffset.onOffset
        api/pandas.tseries.offsets.CustomBusinessMonthBegin.onOffset
        api/pandas.tseries.offsets.CustomBusinessMonthEnd.onOffset
        api/pandas.tseries.offsets.BusinessMonthBegin.onOffset
        api/pandas.tseries.offsets.BusinessMonthEnd.onOffset
        api/pandas.tseries.offsets.MonthBegin.onOffset
        api/pandas.tseries.offsets.MonthEnd.onOffset
        api/pandas.tseries.offsets.MonthOffset.onOffset
        api/pandas.tseries.offsets.CustomBusinessHour.onOffset
        api/pandas.tseries.offsets.CustomBusinessDay.onOffset
        api/pandas.tseries.offsets.BusinessHour.onOffset
        ..api/pandas.offsets.DateOffset.isAnchored
        ..api/pandas.offsets.BusinessDay.isAnchored
        ..api/pandas.offsets.CDay.isAnchored
        ..api/pandas.offsets.CBMonthBegin.isAnchored
        ..api/pandas.offsets.CBMonthEnd.isAnchored
        ..api/pandas.offsets.BMonthBegin.isAnchored
        ..api/pandas.offsets.BMonthEnd.isAnchored
        ..api/pandas.offsets.BDay.isAnchored
        ..api/pandas.offsets.Nano.isAnchored
        ..api/pandas.offsets.Micro.isAnchored
        ..api/pandas.offsets.Milli.isAnchored
        ..api/pandas.offsets.Second.isAnchored
        ..api/pandas.offsets.Minute.isAnchored
        ..api/pandas.offsets.Hour.isAnchored
        ..api/pandas.offsets.Day.isAnchored
        ..api/pandas.offsets.Tick.isAnchored
        ..api/pandas.offsets.Easter.isAnchored
        ..api/pandas.offsets.FY5253Quarter.isAnchored
        ..api/pandas.offsets.FY5253.isAnchored
        ..api/pandas.offsets.YearBegin.isAnchored
        ..api/pandas.offsets.YearEnd.isAnchored
        ..api/pandas.offsets.BYearBegin.isAnchored
        ..api/pandas.offsets.BYearEnd.isAnchored
        ..api/pandas.offsets.YearOffset.isAnchored
        ..api/pandas.offsets.QuarterBegin.isAnchored
        ..api/pandas.offsets.QuarterEnd.isAnchored
        ..api/pandas.offsets.BQuarterBegin.isAnchored
        ..api/pandas.offsets.BQuarterEnd.isAnchored
        ..api/pandas.offsets.QuarterOffset.isAnchored
        ..api/pandas.offsets.LastWeekOfMonth.isAnchored
        ..api/pandas.offsets.WeekOfMonth.isAnchored
        ..api/pandas.offsets.Week.isAnchored
        ..api/pandas.offsets.SemiMonthBegin.isAnchored
        ..api/pandas.offsets.SemiMonthEnd.isAnchored
        ..api/pandas.offsets.SemiMonthOffset.isAnchored
        ..api/pandas.offsets.CustomBusinessMonthBegin.isAnchored
        ..api/pandas.offsets.CustomBusinessMonthEnd.isAnchored
        ..api/pandas.offsets.BusinessMonthBegin.isAnchored
        ..api/pandas.offsets.BusinessMonthEnd.isAnchored
        ..api/pandas.offsets.MonthBegin.isAnchored
        ..api/pandas.offsets.MonthEnd.isAnchored
        ..api/pandas.offsets.MonthOffset.isAnchored
        ..api/pandas.offsets.CustomBusinessHour.isAnchored
        ..api/pandas.offsets.CustomBusinessDay.isAnchored
        ..api/pandas.offsets.BusinessHour.isAnchored
        ..api/pandas.offsets.DateOffset.onOffset
        ..api/pandas.offsets.BusinessDay.onOffset
        ..api/pandas.offsets.CDay.onOffset
        ..api/pandas.offsets.CBMonthBegin.onOffset
        ..api/pandas.offsets.CBMonthEnd.onOffset
        ..api/pandas.offsets.BMonthBegin.onOffset
        ..api/pandas.offsets.BMonthEnd.onOffset
        ..api/pandas.offsets.BDay.onOffset
        ..api/pandas.offsets.Nano.onOffset
        ..api/pandas.offsets.Micro.onOffset
        ..api/pandas.offsets.Milli.onOffset
        ..api/pandas.offsets.Second.onOffset
        ..api/pandas.offsets.Minute.onOffset
        ..api/pandas.offsets.Hour.onOffset
        ..api/pandas.offsets.Day.onOffset
        ..api/pandas.offsets.Tick.onOffset
        ..api/pandas.offsets.Easter.onOffset
        ..api/pandas.offsets.FY5253Quarter.onOffset
        ..api/pandas.offsets.FY5253.onOffset
        ..api/pandas.offsets.YearBegin.onOffset
        ..api/pandas.offsets.YearEnd.onOffset
        ..api/pandas.offsets.BYearBegin.onOffset
        ..api/pandas.offsets.BYearEnd.onOffset
        ..api/pandas.offsets.YearOffset.onOffset
        ..api/pandas.offsets.QuarterBegin.onOffset
        ..api/pandas.offsets.QuarterEnd.onOffset
        ..api/pandas.offsets.BQuarterBegin.onOffset
        ..api/pandas.offsets.BQuarterEnd.onOffset
        ..api/pandas.offsets.QuarterOffset.onOffset
        ..api/pandas.offsets.LastWeekOfMonth.onOffset
        ..api/pandas.offsets.WeekOfMonth.onOffset
        ..api/pandas.offsets.Week.onOffset
        ..api/pandas.offsets.SemiMonthBegin.onOffset
        ..api/pandas.offsets.SemiMonthEnd.onOffset
        ..api/pandas.offsets.SemiMonthOffset.onOffset
        ..api/pandas.offsets.CustomBusinessMonthBegin.onOffset
        ..api/pandas.offsets.CustomBusinessMonthEnd.onOffset
        ..api/pandas.offsets.BusinessMonthBegin.onOffset
        ..api/pandas.offsets.BusinessMonthEnd.onOffset
        ..api/pandas.offsets.MonthBegin.onOffset
        ..api/pandas.offsets.MonthEnd.onOffset
        ..api/pandas.offsets.MonthOffset.onOffset
        ..api/pandas.offsets.CustomBusinessHour.onOffset
        ..api/pandas.offsets.CustomBusinessDay.onOffset
        ..api/pandas.offsets.BusinessHour.onOffset
        api/pandas.offsets.DateOffset.isAnchored
        api/pandas.offsets.BusinessDay.isAnchored
        api/pandas.offsets.CDay.isAnchored
        api/pandas.offsets.CBMonthBegin.isAnchored
        api/pandas.offsets.CBMonthEnd.isAnchored
        api/pandas.offsets.BMonthBegin.isAnchored
        api/pandas.offsets.BMonthEnd.isAnchored
        api/pandas.offsets.BDay.isAnchored
        api/pandas.offsets.Nano.isAnchored
        api/pandas.offsets.Micro.isAnchored
        api/pandas.offsets.Milli.isAnchored
        api/pandas.offsets.Second.isAnchored
        api/pandas.offsets.Minute.isAnchored
        api/pandas.offsets.Hour.isAnchored
        api/pandas.offsets.Day.isAnchored
        api/pandas.offsets.Tick.isAnchored
        api/pandas.offsets.Easter.isAnchored
        api/pandas.offsets.FY5253Quarter.isAnchored
        api/pandas.offsets.FY5253.isAnchored
        api/pandas.offsets.YearBegin.isAnchored
        api/pandas.offsets.YearEnd.isAnchored
        api/pandas.offsets.BYearBegin.isAnchored
        api/pandas.offsets.BYearEnd.isAnchored
        api/pandas.offsets.YearOffset.isAnchored
        api/pandas.offsets.QuarterBegin.isAnchored
        api/pandas.offsets.QuarterEnd.isAnchored
        api/pandas.offsets.BQuarterBegin.isAnchored
        api/pandas.offsets.BQuarterEnd.isAnchored
        api/pandas.offsets.QuarterOffset.isAnchored
        api/pandas.offsets.LastWeekOfMonth.isAnchored
        api/pandas.offsets.WeekOfMonth.isAnchored
        api/pandas.offsets.Week.isAnchored
        api/pandas.offsets.SemiMonthBegin.isAnchored
        api/pandas.offsets.SemiMonthEnd.isAnchored
        api/pandas.offsets.SemiMonthOffset.isAnchored
        api/pandas.offsets.CustomBusinessMonthBegin.isAnchored
        api/pandas.offsets.CustomBusinessMonthEnd.isAnchored
        api/pandas.offsets.BusinessMonthBegin.isAnchored
        api/pandas.offsets.BusinessMonthEnd.isAnchored
        api/pandas.offsets.MonthBegin.isAnchored
        api/pandas.offsets.MonthEnd.isAnchored
        api/pandas.offsets.MonthOffset.isAnchored
        api/pandas.offsets.CustomBusinessHour.isAnchored
        api/pandas.offsets.CustomBusinessDay.isAnchored
        api/pandas.offsets.BusinessHour.isAnchored
        api/pandas.offsets.DateOffset.onOffset
        api/pandas.offsets.BusinessDay.onOffset
        api/pandas.offsets.CDay.onOffset
        api/pandas.offsets.CBMonthBegin.onOffset
        api/pandas.offsets.CBMonthEnd.onOffset
        api/pandas.offsets.BMonthBegin.onOffset
        api/pandas.offsets.BMonthEnd.onOffset
        api/pandas.offsets.BDay.onOffset
        api/pandas.offsets.Nano.onOffset
        api/pandas.offsets.Micro.onOffset
        api/pandas.offsets.Milli.onOffset
        api/pandas.offsets.Second.onOffset
        api/pandas.offsets.Minute.onOffset
        api/pandas.offsets.Hour.onOffset
        api/pandas.offsets.Day.onOffset
        api/pandas.offsets.Tick.onOffset
        api/pandas.offsets.Easter.onOffset
        api/pandas.offsets.FY5253Quarter.onOffset
        api/pandas.offsets.FY5253.onOffset
        api/pandas.offsets.YearBegin.onOffset
        api/pandas.offsets.YearEnd.onOffset
        api/pandas.offsets.BYearBegin.onOffset
        api/pandas.offsets.BYearEnd.onOffset
        api/pandas.offsets.YearOffset.onOffset
        api/pandas.offsets.QuarterBegin.onOffset
        api/pandas.offsets.QuarterEnd.onOffset
        api/pandas.offsets.BQuarterBegin.onOffset
        api/pandas.offsets.BQuarterEnd.onOffset
        api/pandas.offsets.QuarterOffset.onOffset
        api/pandas.offsets.LastWeekOfMonth.onOffset
        api/pandas.offsets.WeekOfMonth.onOffset
        api/pandas.offsets.Week.onOffset
        api/pandas.offsets.SemiMonthBegin.onOffset
        api/pandas.offsets.SemiMonthEnd.onOffset
        api/pandas.offsets.SemiMonthOffset.onOffset
        api/pandas.offsets.CustomBusinessMonthBegin.onOffset
        api/pandas.offsets.CustomBusinessMonthEnd.onOffset
        api/pandas.offsets.BusinessMonthBegin.onOffset
        api/pandas.offsets.BusinessMonthEnd.onOffset
        api/pandas.offsets.MonthBegin.onOffset
        api/pandas.offsets.MonthEnd.onOffset
        api/pandas.offsets.MonthOffset.onOffset
        api/pandas.offsets.CustomBusinessHour.onOffset
        api/pandas.offsets.CustomBusinessDay.onOffset
        api/pandas.offsets.BusinessHour.onOffset
        ..api/pandas.DateOffset.isAnchored
        ..api/pandas.BusinessDay.isAnchored
        ..api/pandas.CDay.isAnchored
        ..api/pandas.CBMonthBegin.isAnchored
        ..api/pandas.CBMonthEnd.isAnchored
        ..api/pandas.BMonthBegin.isAnchored
        ..api/pandas.BMonthEnd.isAnchored
        ..api/pandas.BDay.isAnchored
        ..api/pandas.Nano.isAnchored
        ..api/pandas.Micro.isAnchored
        ..api/pandas.Milli.isAnchored
        ..api/pandas.Second.isAnchored
        ..api/pandas.Minute.isAnchored
        ..api/pandas.Hour.isAnchored
        ..api/pandas.Day.isAnchored
        ..api/pandas.Tick.isAnchored
        ..api/pandas.Easter.isAnchored
        ..api/pandas.FY5253Quarter.isAnchored
        ..api/pandas.FY5253.isAnchored
        ..api/pandas.YearBegin.isAnchored
        ..api/pandas.YearEnd.isAnchored
        ..api/pandas.BYearBegin.isAnchored
        ..api/pandas.BYearEnd.isAnchored
        ..api/pandas.YearOffset.isAnchored
        ..api/pandas.QuarterBegin.isAnchored
        ..api/pandas.QuarterEnd.isAnchored
        ..api/pandas.BQuarterBegin.isAnchored
        ..api/pandas.BQuarterEnd.isAnchored
        ..api/pandas.QuarterOffset.isAnchored
        ..api/pandas.LastWeekOfMonth.isAnchored
        ..api/pandas.WeekOfMonth.isAnchored
        ..api/pandas.Week.isAnchored
        ..api/pandas.SemiMonthBegin.isAnchored
        ..api/pandas.SemiMonthEnd.isAnchored
        ..api/pandas.SemiMonthOffset.isAnchored
        ..api/pandas.CustomBusinessMonthBegin.isAnchored
        ..api/pandas.CustomBusinessMonthEnd.isAnchored
        ..api/pandas.BusinessMonthBegin.isAnchored
        ..api/pandas.BusinessMonthEnd.isAnchored
        ..api/pandas.MonthBegin.isAnchored
        ..api/pandas.MonthEnd.isAnchored
        ..api/pandas.MonthOffset.isAnchored
        ..api/pandas.CustomBusinessHour.isAnchored
        ..api/pandas.CustomBusinessDay.isAnchored
        ..api/pandas.BusinessHour.isAnchored
        ..api/pandas.DateOffset.onOffset
        ..api/pandas.BusinessDay.onOffset
        ..api/pandas.CDay.onOffset
        ..api/pandas.CBMonthBegin.onOffset
        ..api/pandas.CBMonthEnd.onOffset
        ..api/pandas.BMonthBegin.onOffset
        ..api/pandas.BMonthEnd.onOffset
        ..api/pandas.BDay.onOffset
        ..api/pandas.Nano.onOffset
        ..api/pandas.Micro.onOffset
        ..api/pandas.Milli.onOffset
        ..api/pandas.Second.onOffset
        ..api/pandas.Minute.onOffset
        ..api/pandas.Hour.onOffset
        ..api/pandas.Day.onOffset
        ..api/pandas.Tick.onOffset
        ..api/pandas.Easter.onOffset
        ..api/pandas.FY5253Quarter.onOffset
        ..api/pandas.FY5253.onOffset
        ..api/pandas.YearBegin.onOffset
        ..api/pandas.YearEnd.onOffset
        ..api/pandas.BYearBegin.onOffset
        ..api/pandas.BYearEnd.onOffset
        ..api/pandas.YearOffset.onOffset
        ..api/pandas.QuarterBegin.onOffset
        ..api/pandas.QuarterEnd.onOffset
        ..api/pandas.BQuarterBegin.onOffset
        ..api/pandas.BQuarterEnd.onOffset
        ..api/pandas.QuarterOffset.onOffset
        ..api/pandas.LastWeekOfMonth.onOffset
        ..api/pandas.WeekOfMonth.onOffset
        ..api/pandas.Week.onOffset
        ..api/pandas.SemiMonthBegin.onOffset
        ..api/pandas.SemiMonthEnd.onOffset
        ..api/pandas.SemiMonthOffset.onOffset
        ..api/pandas.CustomBusinessMonthBegin.onOffset
        ..api/pandas.CustomBusinessMonthEnd.onOffset
        ..api/pandas.BusinessMonthBegin.onOffset
        ..api/pandas.BusinessMonthEnd.onOffset
        ..api/pandas.MonthBegin.onOffset
        ..api/pandas.MonthEnd.onOffset
        ..api/pandas.MonthOffset.onOffset
        ..api/pandas.CustomBusinessHour.onOffset
        ..api/pandas.CustomBusinessDay.onOffset
        ..api/pandas.BusinessHour.onOffset
        api/pandas.DateOffset.isAnchored
        api/pandas.BusinessDay.isAnchored
        api/pandas.CDay.isAnchored
        api/pandas.CBMonthBegin.isAnchored
        api/pandas.CBMonthEnd.isAnchored
        api/pandas.BMonthBegin.isAnchored
        api/pandas.BMonthEnd.isAnchored
        api/pandas.BDay.isAnchored
        api/pandas.Nano.isAnchored
        api/pandas.Micro.isAnchored
        api/pandas.Milli.isAnchored
        api/pandas.Second.isAnchored
        api/pandas.Minute.isAnchored
        api/pandas.Hour.isAnchored
        api/pandas.Day.isAnchored
        api/pandas.Tick.isAnchored
        api/pandas.Easter.isAnchored
        api/pandas.FY5253Quarter.isAnchored
        api/pandas.FY5253.isAnchored
        api/pandas.YearBegin.isAnchored
        api/pandas.YearEnd.isAnchored
        api/pandas.BYearBegin.isAnchored
        api/pandas.BYearEnd.isAnchored
        api/pandas.YearOffset.isAnchored
        api/pandas.QuarterBegin.isAnchored
        api/pandas.QuarterEnd.isAnchored
        api/pandas.BQuarterBegin.isAnchored
        api/pandas.BQuarterEnd.isAnchored
        api/pandas.QuarterOffset.isAnchored
        api/pandas.LastWeekOfMonth.isAnchored
        api/pandas.WeekOfMonth.isAnchored
        api/pandas.Week.isAnchored
        api/pandas.SemiMonthBegin.isAnchored
        api/pandas.SemiMonthEnd.isAnchored
        api/pandas.SemiMonthOffset.isAnchored
        api/pandas.CustomBusinessMonthBegin.isAnchored
        api/pandas.CustomBusinessMonthEnd.isAnchored
        api/pandas.BusinessMonthBegin.isAnchored
        api/pandas.BusinessMonthEnd.isAnchored
        api/pandas.MonthBegin.isAnchored
        api/pandas.MonthEnd.isAnchored
        api/pandas.MonthOffset.isAnchored
        api/pandas.CustomBusinessHour.isAnchored
        api/pandas.CustomBusinessDay.isAnchored
        api/pandas.BusinessHour.isAnchored
        api/pandas.DateOffset.onOffset
        api/pandas.BusinessDay.onOffset
        api/pandas.CDay.onOffset
        api/pandas.CBMonthBegin.onOffset
        api/pandas.CBMonthEnd.onOffset
        api/pandas.BMonthBegin.onOffset
        api/pandas.BMonthEnd.onOffset
        api/pandas.BDay.onOffset
        api/pandas.Nano.onOffset
        api/pandas.Micro.onOffset
        api/pandas.Milli.onOffset
        api/pandas.Second.onOffset
        api/pandas.Minute.onOffset
        api/pandas.Hour.onOffset
        api/pandas.Day.onOffset
        api/pandas.Tick.onOffset
        api/pandas.Easter.onOffset
        api/pandas.FY5253Quarter.onOffset
        api/pandas.FY5253.onOffset
        api/pandas.YearBegin.onOffset
        api/pandas.YearEnd.onOffset
        api/pandas.BYearBegin.onOffset
        api/pandas.BYearEnd.onOffset
        api/pandas.YearOffset.onOffset
        api/pandas.QuarterBegin.onOffset
        api/pandas.QuarterEnd.onOffset
        api/pandas.BQuarterBegin.onOffset
        api/pandas.BQuarterEnd.onOffset
        api/pandas.QuarterOffset.onOffset
        api/pandas.LastWeekOfMonth.onOffset
        api/pandas.WeekOfMonth.onOffset
        api/pandas.Week.onOffset
        api/pandas.SemiMonthBegin.onOffset
        api/pandas.SemiMonthEnd.onOffset
        api/pandas.SemiMonthOffset.onOffset
        api/pandas.CustomBusinessMonthBegin.onOffset
        api/pandas.CustomBusinessMonthEnd.onOffset
        api/pandas.BusinessMonthBegin.onOffset
        api/pandas.BusinessMonthEnd.onOffset
        api/pandas.MonthBegin.onOffset
        api/pandas.MonthEnd.onOffset
        api/pandas.MonthOffset.onOffset
        api/pandas.CustomBusinessHour.onOffset
        api/pandas.CustomBusinessDay.onOffset
        api/pandas.BusinessHour.onOffset

.. _api.frequencies:

===========
Frequencies
===========
.. currentmodule:: pandas.tseries.frequencies

.. _api.offsets:

.. autosummary::
   :toctree: api/

   to_offset
