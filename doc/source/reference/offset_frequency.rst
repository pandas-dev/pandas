{{ header }}

.. _api.dateoffsets:

============
Date offsets
============
.. currentmodule:: pandas.tseries.offsets

BaseOffset
----------
.. autosummary::
   :toctree: api/

    BaseOffset

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

    BaseOffset.freqstr
    BaseOffset.kwds
    DateOffset.name
    BaseOffset.nanos
    BaseOffset.normalize
    BaseOffset.rule_code
    BaseOffset.n

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

    BaseOffset.copy
    BaseOffset.is_on_offset
    BaseOffset.is_month_start
    BaseOffset.is_month_end
    BaseOffset.is_quarter_start
    BaseOffset.is_quarter_end
    BaseOffset.is_year_start
    BaseOffset.is_year_end
    BaseOffset.rollback
    BaseOffset.rollforward

DateOffset
----------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    DateOffset

BusinessDay
-----------

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    BusinessDay
    BDay

Custom properties:

.. autosummary::
   :toctree: api/

    BusinessDay.offset
    BusinessDay.weekmask
    BusinessDay.holidays
    BusinessDay.calendar

BusinessHour
------------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    BusinessHour

Custom properties:

.. autosummary::
   :toctree: api/

    BusinessHour.offset
    BusinessHour.start
    BusinessHour.end
    BusinessHour.weekmask
    BusinessHour.holidays
    BusinessHour.calendar

CustomBusinessDay
-----------------

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    CustomBusinessDay
    CDay

Custom properties:

.. autosummary::
   :toctree: api/

    CustomBusinessDay.weekmask
    CustomBusinessDay.calendar
    CustomBusinessDay.holidays
    CustomBusinessDay.offset

CustomBusinessHour
------------------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    CustomBusinessHour

Custom properties:

.. autosummary::
   :toctree: api/

    CustomBusinessHour.weekmask
    CustomBusinessHour.calendar
    CustomBusinessHour.holidays
    CustomBusinessHour.start
    CustomBusinessHour.end
    CustomBusinessHour.offset

Month offsets
-------------

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    MonthEnd
    MonthBegin

Business variants:

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    BusinessMonthEnd
    BMonthEnd
    BusinessMonthBegin
    BMonthBegin

CustomBusinessMonthEnd
----------------------

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    CustomBusinessMonthEnd
    CBMonthEnd

Custom properties:

.. autosummary::
   :toctree: api/

    CustomBusinessMonthEnd.m_offset
    CustomBusinessMonthEnd.weekmask
    CustomBusinessMonthEnd.calendar
    CustomBusinessMonthEnd.holidays
    CustomBusinessMonthEnd.offset

CustomBusinessMonthBegin
------------------------

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    CustomBusinessMonthBegin
    CBMonthBegin

Custom properties:

.. autosummary::
   :toctree: api/

    CustomBusinessMonthBegin.m_offset
    CustomBusinessMonthBegin.weekmask
    CustomBusinessMonthBegin.calendar
    CustomBusinessMonthBegin.holidays
    CustomBusinessMonthBegin.offset

SemiMonthEnd
------------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    SemiMonthEnd
    SemiMonthBegin

Custom properties:

.. autosummary::
   :toctree: api/

    SemiMonthEnd.day_of_month
    SemiMonthBegin.day_of_month

Week
----
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    Week

Custom properties:

.. autosummary::
   :toctree: api/

    Week.weekday

WeekOfMonth
-----------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    WeekOfMonth

Custom properties:

.. autosummary::
   :toctree: api/

    WeekOfMonth.week
    WeekOfMonth.weekday

LastWeekOfMonth
---------------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    LastWeekOfMonth

Custom properties:

.. autosummary::
   :toctree: api/

    LastWeekOfMonth.weekday

Quarter offsets
---------------

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    QuarterEnd
    QuarterBegin
    BQuarterEnd
    BQuarterBegin

Custom properties:

.. autosummary::
   :toctree: api/

    BQuarterEnd.startingMonth
    BQuarterBegin.startingMonth
    QuarterEnd.startingMonth
    QuarterBegin.startingMonth

HalfYear offsets
----------------

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    HalfYearEnd
    HalfYearBegin
    BHalfYearEnd
    BHalfYearBegin

Custom properties:

.. autosummary::
   :toctree: api/

    BHalfYearEnd.startingMonth
    BHalfYearBegin.startingMonth
    HalfYearEnd.startingMonth
    HalfYearBegin.startingMonth

Year offsets
------------

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    YearEnd
    YearBegin
    BYearEnd
    BYearBegin

Custom properties:

.. autosummary::
   :toctree: api/

    YearEnd.month
    YearBegin.month
    BYearEnd.month
    BYearBegin.month

FY5253
------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    FY5253

Custom properties:

.. autosummary::
   :toctree: api/

    FY5253.startingMonth
    FY5253.variation
    FY5253.weekday

Custom methods:

.. autosummary::
   :toctree: api/

    FY5253.get_rule_code_suffix
    FY5253.get_year_end

FY5253Quarter
-------------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    FY5253Quarter

Custom properties:

.. autosummary::
   :toctree: api/

    FY5253Quarter.qtr_with_extra_week
    FY5253Quarter.startingMonth
    FY5253Quarter.variation
    FY5253Quarter.weekday

Custom methods:

.. autosummary::
   :toctree: api/

    FY5253Quarter.get_rule_code_suffix
    FY5253Quarter.get_weeks
    FY5253Quarter.year_has_extra_week

Easter
------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    Easter

Custom properties:

.. autosummary::
   :toctree: api/

    Easter.method

Ticks
-----

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

    Tick
    Day
    Hour
    Minute
    Second
    Milli
    Micro
    Nano

.. _api.frequencies:

===========
Frequencies
===========
.. currentmodule:: pandas.tseries.frequencies

.. _api.offsets:

.. autosummary::
   :toctree: api/

   to_offset

.. _api.holiday:

==================
Holiday calendars
==================
.. currentmodule:: pandas.tseries.holiday

.. autosummary::
   :toctree: api/

   AbstractHolidayCalendar
   Holiday
