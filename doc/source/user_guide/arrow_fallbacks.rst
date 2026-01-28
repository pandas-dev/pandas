.. _arrow-fallbacks:

{{ header }}

***************
Arrow Fallbacks
***************

This document shows the runtime behavior of pandas methods on
Arrow-backed arrays. Results are determined by actually running
each operation and observing the outcome.

.. note::

   Results may vary based on method parameters. For example, some string
   methods like ``split(expand=True)`` use Arrow, while
   ``split(expand=False)`` may fall back to NumPy. This table shows
   results for default or common parameter values only.

.. note::

   Some methods have version-gated behavior that depends on your PyArrow
   version. See :ref:`arrow-version-gated` below for details.

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
   * - |mixed|
     - Mixed results within dtype group

.. |arrow| replace:: ✓ Arrow
.. |numpy| replace:: → NumPy
.. |elementwise| replace:: ⟳ Elem
.. |object| replace:: → Object
.. |notimpl| replace:: ✗ N/I
.. |typeerror| replace:: ✗ Type
.. |error| replace:: ✗ Err
.. |mixed| replace:: ~ Mixed


String Methods (Series.str.*)
=============================

.. list-table::
   :widths: 25 10
   :header-rows: 1

   * - Method
     - string
   * - ``capitalize``
     - |arrow|
   * - ``casefold``
     - |mixed|
   * - ``center``
     - |arrow|
   * - ``contains(flags=0)``
     - |mixed|
   * - ``contains(flags=re.I)``
     - |mixed|
   * - ``count(flags=0)``
     - |mixed|
   * - ``count(flags=re.I)``
     - |mixed|
   * - ``encode``
     - |mixed|
   * - ``endswith``
     - |mixed|
   * - ``extract(expand=False)``
     - |arrow|
   * - ``extract(expand=True)``
     - |arrow|
   * - ``extractall``
     - |arrow|
   * - ``find``
     - |mixed|
   * - ``findall(flags=0)``
     - |mixed|
   * - ``findall(flags=re.I)``
     - |mixed|
   * - ``fullmatch(flags=0)``
     - |mixed|
   * - ``fullmatch(flags=re.I)``
     - |mixed|
   * - ``get``
     - |arrow|
   * - ``get_dummies``
     - |mixed|
   * - ``index``
     - |error|
   * - ``isalnum``
     - |mixed|
   * - ``isalpha``
     - |mixed|
   * - ``isascii``
     - |mixed|
   * - ``isdecimal``
     - |mixed|
   * - ``isdigit``
     - |mixed|
   * - ``islower``
     - |mixed|
   * - ``isnumeric``
     - |mixed|
   * - ``isspace``
     - |mixed|
   * - ``istitle``
     - |mixed|
   * - ``isupper``
     - |mixed|
   * - ``join``
     - |mixed|
   * - ``len``
     - |mixed|
   * - ``ljust``
     - |arrow|
   * - ``lower``
     - |arrow|
   * - ``lstrip``
     - |arrow|
   * - ``match(flags=0)``
     - |mixed|
   * - ``match(flags=re.I)``
     - |mixed|
   * - ``normalize``
     - |mixed|
   * - ``pad(side=both)``
     - |arrow|
   * - ``pad(side=left)``
     - |arrow|
   * - ``pad(side=right)``
     - |arrow|
   * - ``partition``
     - |mixed|
   * - ``removeprefix``
     - |arrow|
   * - ``removesuffix``
     - |arrow|
   * - ``repeat``
     - |arrow|
   * - ``replace(case=False)``
     - |mixed|
   * - ``replace(case=True)``
     - |arrow|
   * - ``replace(repl=callable)``
     - |mixed|
   * - ``rfind``
     - |mixed|
   * - ``rindex``
     - |error|
   * - ``rjust``
     - |arrow|
   * - ``rpartition``
     - |mixed|
   * - ``rsplit(expand=False)``
     - |mixed|
   * - ``rsplit(expand=True)``
     - |arrow|
   * - ``rstrip``
     - |arrow|
   * - ``slice``
     - |arrow|
   * - ``slice_replace``
     - |arrow|
   * - ``split(expand=False)``
     - |mixed|
   * - ``split(expand=True)``
     - |arrow|
   * - ``startswith``
     - |mixed|
   * - ``strip``
     - |arrow|
   * - ``swapcase``
     - |arrow|
   * - ``title``
     - |arrow|
   * - ``translate``
     - |mixed|
   * - ``upper``
     - |arrow|
   * - ``wrap``
     - |mixed|
   * - ``zfill``
     - |arrow|

Datetime Methods (Series.dt.*)
==============================

.. list-table::
   :widths: 25 10
   :header-rows: 1

   * - Method
     - timestamp
   * - ``as_unit``
     - |arrow|
   * - ``ceil``
     - |arrow|
   * - ``date``
     - |arrow|
   * - ``day``
     - |arrow|
   * - ``day_name``
     - |arrow|
   * - ``day_of_week``
     - |arrow|
   * - ``day_of_year``
     - |arrow|
   * - ``dayofweek``
     - |arrow|
   * - ``dayofyear``
     - |arrow|
   * - ``days_in_month``
     - |arrow|
   * - ``daysinmonth``
     - |arrow|
   * - ``floor``
     - |arrow|
   * - ``hour``
     - |arrow|
   * - ``is_leap_year``
     - |arrow|
   * - ``is_month_end``
     - |arrow|
   * - ``is_month_start``
     - |arrow|
   * - ``is_quarter_end``
     - |arrow|
   * - ``is_quarter_start``
     - |arrow|
   * - ``is_year_end``
     - |arrow|
   * - ``is_year_start``
     - |arrow|
   * - ``isocalendar``
     - |arrow|
   * - ``microsecond``
     - |arrow|
   * - ``minute``
     - |arrow|
   * - ``month``
     - |arrow|
   * - ``month_name``
     - |arrow|
   * - ``nanosecond``
     - |arrow|
   * - ``normalize``
     - |arrow|
   * - ``quarter``
     - |arrow|
   * - ``round``
     - |arrow|
   * - ``second``
     - |arrow|
   * - ``strftime``
     - |typeerror|
   * - ``time``
     - |arrow|
   * - ``to_pydatetime``
     - |object|
   * - ``tz``
     - |arrow|
   * - ``tz_convert``
     - |arrow|
   * - ``tz_localize(ambiguous=NaT)``
     - |notimpl|
   * - ``tz_localize(ambiguous=raise)``
     - |arrow|
   * - ``unit``
     - |arrow|
   * - ``weekday``
     - |arrow|
   * - ``year``
     - |arrow|

Timedelta Methods (Series.dt.*)
===============================

.. list-table::
   :widths: 25 10
   :header-rows: 1

   * - Method
     - duration
   * - ``as_unit``
     - |arrow|
   * - ``days``
     - |arrow|
   * - ``microseconds``
     - |arrow|
   * - ``nanoseconds``
     - |arrow|
   * - ``seconds``
     - |arrow|
   * - ``to_pytimedelta``
     - |arrow|
   * - ``total_seconds``
     - |arrow|

Aggregation Methods
===================

.. list-table::
   :widths: 25 10 10
   :header-rows: 1

   * - Method
     - integer
     - float
   * - ``all``
     - |arrow|
     - |arrow|
   * - ``any``
     - |arrow|
     - |arrow|
   * - ``count``
     - |arrow|
     - |arrow|
   * - ``kurt``
     - |typeerror|
     - |typeerror|
   * - ``max``
     - |arrow|
     - |arrow|
   * - ``mean``
     - |arrow|
     - |arrow|
   * - ``median``
     - |arrow|
     - |arrow|
   * - ``min``
     - |arrow|
     - |arrow|
   * - ``prod``
     - |arrow|
     - |arrow|
   * - ``sem``
     - |arrow|
     - |arrow|
   * - ``skew``
     - |arrow|
     - |arrow|
   * - ``std``
     - |arrow|
     - |arrow|
   * - ``sum``
     - |arrow|
     - |arrow|
   * - ``var``
     - |arrow|
     - |arrow|

Array Methods
=============

.. list-table::
   :widths: 25 10 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Method
     - string
     - integer
     - float
     - bool
     - timestamp
     - date
     - duration
     - time
     - binary
   * - ``abs``
     -
     - |arrow|
     - |arrow|
     -
     - |notimpl|
     - |notimpl|
     - |arrow|
     - |notimpl|
     -
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
   * - ``clip``
     -
     - |arrow|
     - |arrow|
     -
     - |mixed|
     - |typeerror|
     - |arrow|
     - |typeerror|
     -
   * - ``cummax``
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
   * - ``cummin``
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
   * - ``cumprod``
     -
     - |arrow|
     - |arrow|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     -
   * - ``cumsum``
     -
     - |arrow|
     - |arrow|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |arrow|
     - |typeerror|
     -
   * - ``diff``
     -
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     -
   * - ``drop_duplicates``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
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
   * - ``fillna(limit=1)``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |error|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``fillna(limit=None)``
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |error|
     - |arrow|
     - |arrow|
     - |arrow|
   * - ``interpolate(method=linear)``
     -
     - |arrow|
     - |arrow|
     -
     - |typeerror|
     - |typeerror|
     - |typeerror|
     - |typeerror|
     -
   * - ``interpolate(method=pad)``
     -
     - |error|
     - |error|
     -
     - |error|
     - |error|
     - |error|
     - |error|
     -
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
   * - ``round``
     -
     -
     - |arrow|
     -
     -
     -
     -
     -
     -
   * - ``searchsorted``
     - |numpy|
     - |numpy|
     - |numpy|
     - |numpy|
     - |mixed|
     - |typeerror|
     - |numpy|
     - |typeerror|
     - |numpy|
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

Algorithmic Methods
===================

These methods have behavior that may vary across versions.
Explicit kwargs are used for deterministic results.

.. list-table::
   :widths: 25 10 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Method
     - string
     - integer
     - float
     - bool
     - timestamp
     - date
     - duration
     - time
     - binary
   * - ``factorize``
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
     - |mixed|
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

Arithmetic Operations
=====================

.. list-table::
   :widths: 25 10 10
   :header-rows: 1

   * - Method
     - integer
     - float
   * - ``add``
     - |arrow|
     - |arrow|
   * - ``floordiv``
     - |arrow|
     - |arrow|
   * - ``mod``
     - |notimpl|
     - |notimpl|
   * - ``mul``
     - |arrow|
     - |arrow|
   * - ``pow``
     - |arrow|
     - |arrow|
   * - ``sub``
     - |arrow|
     - |arrow|
   * - ``truediv``
     - |arrow|
     - |arrow|

Comparison Operations
=====================

.. list-table::
   :widths: 25 10 10
   :header-rows: 1

   * - Method
     - integer
     - float
   * - ``eq``
     - |arrow|
     - |arrow|
   * - ``ge``
     - |arrow|
     - |arrow|
   * - ``gt``
     - |arrow|
     - |arrow|
   * - ``le``
     - |arrow|
     - |arrow|
   * - ``lt``
     - |arrow|
     - |arrow|
   * - ``ne``
     - |arrow|
     - |arrow|

.. _arrow-version-gated:

Version-Gated Methods
=====================

The following methods have behavior that depends on your installed PyArrow
version. The table above shows results for the current environment; if you
have an older PyArrow version, some methods may fall back to slower paths.

.. list-table::
   :widths: 25 15 60
   :header-rows: 1

   * - Method
     - Min PyArrow
     - Notes
   * - ``str.center``
     - 17.0
     - Uses pc.utf8_center for side='both' padding. Below 17.0: object; 17.0+: arrow.
   * - ``str.isdigit``
     - 21.0
     - PyArrow < 21.0 has incorrect utf8_is_digit for some digits. Below 21.0: elementwise; 21.0+: arrow.
   * - ``str.pad(side=both)``
     - 17.0
     - Uses pc.utf8_center for side='both' padding. Below 17.0: object; 17.0+: arrow.
   * - ``str.zfill``
     - 21.0
     - pc.utf8_zfill was added in PyArrow 21.0. Below 21.0: elementwise; 21.0+: arrow.
