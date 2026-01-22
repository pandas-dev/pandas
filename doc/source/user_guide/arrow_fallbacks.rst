.. _arrow-fallbacks:

{{ header }}

***************
Arrow Fallbacks
***************

This document shows the runtime behavior of pandas methods on
Arrow-backed arrays. Results are determined by actually running
each operation and observing the outcome.

Legend
======

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Status
     - Description
   * - |arrow|
     - Operation returns Arrow-backed result
   * - |numpy|
     - Operation falls back to NumPy
   * - |elementwise|
     - Uses element-wise Python processing
   * - |object|
     - Returns object dtype
   * - |notimpl|
     - Raises NotImplementedError
   * - |typeerror|
     - Not supported for this dtype
   * - |error|
     - Other error

.. |arrow| replace:: ✓ Arrow
.. |numpy| replace:: → NumPy
.. |elementwise| replace:: ⟳ Elem
.. |object| replace:: → Object
.. |notimpl| replace:: ✗ N/I
.. |typeerror| replace:: ✗ Type
.. |error| replace:: ✗ Err


String Methods (Series.str.*)
=============================

.. list-table::
   :widths: 20 10 10
   :header-rows: 1

   * - Method
     - large_string
     - string
   * - ``capitalize``
     - |arrow|
     - |arrow|
   * - ``casefold``
     - |elementwise|
     - |numpy|
   * - ``center``
     - |arrow|
     - |arrow|
   * - ``contains``
     - |arrow|
     - |numpy|
   * - ``count``
     - |arrow|
     - |numpy|
   * - ``encode``
     - |elementwise|
     - |numpy|
   * - ``endswith``
     - |arrow|
     - |numpy|
   * - ``extract``
     - |error|
     - |numpy|
   * - ``extractall``
     - |arrow|
     - |arrow|
   * - ``find``
     - |arrow|
     - |numpy|
   * - ``findall``
     - |elementwise|
     - |numpy|
   * - ``fullmatch``
     - |arrow|
     - |numpy|
   * - ``get``
     - |arrow|
     - |arrow|
   * - ``get_dummies``
     - |numpy|
     - |numpy|
   * - ``index``
     - |error|
     - |error|
   * - ``isalnum``
     - |arrow|
     - |numpy|
   * - ``isalpha``
     - |arrow|
     - |numpy|
   * - ``isascii``
     - |arrow|
     - |numpy|
   * - ``isdecimal``
     - |arrow|
     - |numpy|
   * - ``isdigit``
     - |arrow|
     - |numpy|
   * - ``islower``
     - |arrow|
     - |numpy|
   * - ``isnumeric``
     - |arrow|
     - |numpy|
   * - ``isspace``
     - |arrow|
     - |numpy|
   * - ``istitle``
     - |arrow|
     - |numpy|
   * - ``isupper``
     - |arrow|
     - |numpy|
   * - ``join``
     - |elementwise|
     - |numpy|
   * - ``len``
     - |arrow|
     - |numpy|
   * - ``ljust``
     - |arrow|
     - |arrow|
   * - ``lower``
     - |arrow|
     - |arrow|
   * - ``lstrip``
     - |arrow|
     - |arrow|
   * - ``match``
     - |arrow|
     - |numpy|
   * - ``normalize``
     - |elementwise|
     - |numpy|
   * - ``pad``
     - |arrow|
     - |arrow|
   * - ``partition``
     - |elementwise|
     - |numpy|
   * - ``removeprefix``
     - |arrow|
     - |arrow|
   * - ``removesuffix``
     - |arrow|
     - |arrow|
   * - ``repeat``
     - |arrow|
     - |arrow|
   * - ``replace``
     - |arrow|
     - |arrow|
   * - ``rfind``
     - |elementwise|
     - |numpy|
   * - ``rindex``
     - |error|
     - |error|
   * - ``rjust``
     - |arrow|
     - |arrow|
   * - ``rpartition``
     - |elementwise|
     - |numpy|
   * - ``rsplit``
     - |arrow|
     - |numpy|
   * - ``rstrip``
     - |arrow|
     - |arrow|
   * - ``slice``
     - |arrow|
     - |arrow|
   * - ``slice_replace``
     - |arrow|
     - |arrow|
   * - ``split``
     - |arrow|
     - |numpy|
   * - ``startswith``
     - |arrow|
     - |numpy|
   * - ``strip``
     - |arrow|
     - |arrow|
   * - ``swapcase``
     - |arrow|
     - |arrow|
   * - ``title``
     - |arrow|
     - |arrow|
   * - ``translate``
     - |elementwise|
     - |numpy|
   * - ``upper``
     - |arrow|
     - |arrow|
   * - ``wrap``
     - |elementwise|
     - |numpy|
   * - ``zfill``
     - |arrow|
     - |numpy|

Datetime Methods (Series.dt.*)
==============================

.. list-table::
   :widths: 20 10 10 10
   :header-rows: 1

   * - Method
     - timestamp[ns]
     - timestamp[us_tz]
     - timestamp[us]
   * - ``as_unit``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``ceil``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``date``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``day``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``day_name``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``day_of_week``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``day_of_year``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``dayofweek``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``dayofyear``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``days_in_month``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``daysinmonth``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``floor``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``hour``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``is_leap_year``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``is_month_end``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``is_month_start``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``is_quarter_end``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``is_quarter_start``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``is_year_end``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``is_year_start``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``isocalendar``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``microsecond``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``minute``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``month``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``month_name``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``nanosecond``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``normalize``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``quarter``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``round``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``second``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``strftime``
     - |typeerror|
     - |typeerror|
     - |typeerror|
   * - ``time``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``to_pydatetime``
     - |object|
     - |object|
     - |object|
   * - ``tz``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``tz_convert``
     - |error|
     - |arrow|
     - |error|
   * - ``tz_localize``
     - |arrow|
     - |error|
     - |arrow|
   * - ``unit``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``weekday``
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``year``
     - |arrow|
     - |arrow|
     - |arrow|

Timedelta Methods (Series.dt.*)
===============================

.. list-table::
   :widths: 20 10 10
   :header-rows: 1

   * - Method
     - duration[ns]
     - duration[us]
   * - ``as_unit``
     - |arrow|
     - |arrow|
   * - ``days``
     - |numpy|
     - |numpy|
   * - ``microseconds``
     - |numpy|
     - |numpy|
   * - ``nanoseconds``
     - |numpy|
     - |numpy|
   * - ``seconds``
     - |numpy|
     - |numpy|
   * - ``to_pytimedelta``
     - |arrow|
     - |arrow|
   * - ``total_seconds``
     - |arrow|
     - |arrow|

Aggregation Methods
===================

.. list-table::
   :widths: 20 10 10 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Method
     - double
     - float
     - int16
     - int32
     - int64
     - int8
     - uint16
     - uint32
     - uint64
     - uint8
   * - ``all``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``any``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``count``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``kurt``
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
   * - ``max``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``mean``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``median``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``min``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``prod``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``sem``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``skew``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``std``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``sum``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``var``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|

Array Methods
=============

.. list-table::
   :widths: 20 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Method
     - binary
     - bool
     - date32[day]
     - date64[ms]
     - double
     - duration[ns]
     - duration[us]
     - float
     - int16
     - int32
     - int64
     - int8
     - large_binary
     - large_string
     - string
     - time64[us]
     - timestamp[ns]
     - timestamp[us_tz]
     - timestamp[us]
     - uint16
     - uint32
     - uint64
     - uint8
   * - ``abs``
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``argsort``
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
   * - ``bfill``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``clip``
     - |error|
     - |error|
     - |typeerror|
     - |typeerror|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |error|
     - |error|
     - |error|
     - |numpy|
     - |numpy|
     - |typeerror|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
   * - ``cummax``
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``cummin``
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``cumprod``
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``cumsum``
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``diff``
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``drop_duplicates``
     - |arrow|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |arrow|
     - |arrow|
     - |arrow|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
   * - ``dropna``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``duplicated``
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
   * - ``factorize``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``ffill``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``fillna``
     - |arrow|
     - |arrow|
     - |error|
     - |error|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``interpolate``
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``isna``
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
   * - ``notna``
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
   * - ``rank``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |numpy|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``round``
     - |error|
     - |error|
     - |error|
     - |error|
     - |arrow|
     - |error|
     - |error|
     - |arrow|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
   * - ``searchsorted``
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
     - |error|
   * - ``shift``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``sort_values``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``unique``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``value_counts``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|

Arithmetic Operations
=====================

.. list-table::
   :widths: 20 10 10 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Method
     - double
     - float
     - int16
     - int32
     - int64
     - int8
     - uint16
     - uint32
     - uint64
     - uint8
   * - ``add``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``floordiv``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``mod``
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
   * - ``mul``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``pow``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``sub``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``truediv``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|

Comparison Operations
=====================

.. list-table::
   :widths: 20 10 10 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Method
     - double
     - float
     - int16
     - int32
     - int64
     - int8
     - uint16
     - uint32
     - uint64
     - uint8
   * - ``eq``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``ge``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``gt``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``le``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``lt``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``ne``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
