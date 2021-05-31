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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    DateOffset.apply
    DateOffset.apply_index
    DateOffset.copy
    DateOffset.isAnchored
    DateOffset.onOffset
    DateOffset.is_anchored
    DateOffset.is_on_offset
    DateOffset.__call__

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

    BusinessDay.weekmask
    BusinessDay.holidays
    BusinessDay.calendar

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


BusinessHour
------------
.. autosummary::
   :toctree: api/

    BusinessHour

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BusinessHour.start
    BusinessHour.end
    BusinessHour.weekmask
    BusinessHour.holidays
    BusinessHour.calendar

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

    CustomBusinessDay.weekmask
    CustomBusinessDay.calendar
    CustomBusinessDay.holidays

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


CustomBusinessHour
------------------
.. autosummary::
   :toctree: api/

    CustomBusinessHour

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    CustomBusinessHour.weekmask
    CustomBusinessHour.calendar
    CustomBusinessHour.holidays
    CustomBusinessHour.start
    CustomBusinessHour.end

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


MonthEnd
--------
.. autosummary::
   :toctree: api/

    MonthEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/


Methods
~~~~~~~
.. autosummary::
   :toctree: api/


MonthBegin
----------
.. autosummary::
   :toctree: api/

    MonthBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/


Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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


Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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


Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

    CustomBusinessMonthEnd.m_offset
    CustomBusinessMonthEnd.weekmask
    CustomBusinessMonthEnd.calendar
    CustomBusinessMonthEnd.holidays

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

    CustomBusinessMonthBegin.m_offset
    CustomBusinessMonthBegin.weekmask
    CustomBusinessMonthBegin.calendar
    CustomBusinessMonthBegin.holidays

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


SemiMonthEnd
------------
.. autosummary::
   :toctree: api/

    SemiMonthEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    SemiMonthEnd.day_of_month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


SemiMonthBegin
--------------
.. autosummary::
   :toctree: api/

    SemiMonthBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    SemiMonthBegin.day_of_month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


Week
----
.. autosummary::
   :toctree: api/

    Week

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    Week.weekday

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


WeekOfMonth
-----------
.. autosummary::
   :toctree: api/

    WeekOfMonth

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    WeekOfMonth.week

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    WeekOfMonth.weekday

LastWeekOfMonth
---------------
.. autosummary::
   :toctree: api/

    LastWeekOfMonth

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    LastWeekOfMonth.weekday
    LastWeekOfMonth.week

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


BQuarterEnd
-----------
.. autosummary::
   :toctree: api/

    BQuarterEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BQuarterEnd.startingMonth

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


BQuarterBegin
-------------
.. autosummary::
   :toctree: api/

    BQuarterBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BQuarterBegin.startingMonth

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


QuarterEnd
----------
.. autosummary::
   :toctree: api/

    QuarterEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    QuarterEnd.startingMonth

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


QuarterBegin
------------
.. autosummary::
   :toctree: api/

    QuarterBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    QuarterBegin.startingMonth

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


BYearEnd
--------
.. autosummary::
   :toctree: api/

    BYearEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BYearEnd.month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


BYearBegin
----------
.. autosummary::
   :toctree: api/

    BYearBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BYearBegin.month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


YearEnd
-------
.. autosummary::
   :toctree: api/

    YearEnd

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    YearEnd.month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


YearBegin
---------
.. autosummary::
   :toctree: api/

    YearBegin

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    YearBegin.month

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


FY5253
------
.. autosummary::
   :toctree: api/

    FY5253

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    FY5253.startingMonth
    FY5253.variation
    FY5253.weekday

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    FY5253.get_rule_code_suffix
    FY5253.get_year_end

FY5253Quarter
-------------
.. autosummary::
   :toctree: api/

    FY5253Quarter

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    FY5253Quarter.qtr_with_extra_week
    FY5253Quarter.startingMonth
    FY5253Quarter.variation
    FY5253Quarter.weekday

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    FY5253Quarter.get_rule_code_suffix
    FY5253Quarter.get_weeks
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


Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


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

Methods
~~~~~~~
.. autosummary::
   :toctree: api/


.. _api.frequencies:

===========
Frequencies
===========
.. currentmodule:: pandas.tseries.frequencies

.. _api.offsets:

.. autosummary::
   :toctree: api/

   to_offset
