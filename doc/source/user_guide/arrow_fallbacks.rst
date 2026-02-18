.. _arrow-fallbacks:

{{ header }}

***************
Arrow Fallbacks
***************

This document shows the runtime behavior of pandas methods on
Arrow-backed arrays. Results are determined by actually running
each operation and observing the outcome.

.. note::

   ``string[pyarrow]`` (StringDtype) and ``large_string[pyarrow]``
   (ArrowDtype) use different implementations and are shown in separate
   columns. ``string[pyarrow]`` is the default string dtype; most users
   should focus on that column.

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
   :widths: 25 10 10
   :header-rows: 1

   * - Method
     - string
     - large_string
   * - ``capitalize``
     - |arrow|
     - |arrow|
   * - ``casefold``
     - |arrow|
     - |elementwise|
   * - ``center``
     - |arrow|
     - |arrow|
   * - ``contains(flags=0)``
     - |numpy|
     - |arrow|
   * - ``contains(flags=re.I)``
     - |numpy|
     - |notimpl|
   * - ``count(flags=0)``
     - |numpy|
     - |arrow|
   * - ``count(flags=re.I)``
     - |numpy|
     - |notimpl|
   * - ``encode``
     - |object|
     - |elementwise|
   * - ``endswith``
     - |numpy|
     - |arrow|
   * - ``extract(expand=False)``
     - |arrow|
     - |arrow|
   * - ``extract(expand=True)``
     - |arrow|
     - |arrow|
   * - ``extractall``
     - |arrow|
     - |arrow|
   * - ``find``
     - |numpy|
     - |arrow|
   * - ``findall(flags=0)``
     - |object|
     - |elementwise|
   * - ``findall(flags=re.I)``
     - |object|
     - |elementwise|
   * - ``fullmatch(flags=0)``
     - |numpy|
     - |arrow|
   * - ``fullmatch(flags=re.I)``
     - |numpy|
     - |notimpl|
   * - ``get``
     - |arrow|
     - |arrow|
   * - ``get_dummies``
     - |numpy|
     - |arrow|
   * - ``index``
     - |numpy|
     - |elementwise|
   * - ``isalnum``
     - |numpy|
     - |arrow|
   * - ``isalpha``
     - |numpy|
     - |arrow|
   * - ``isascii``
     - |numpy|
     - |arrow|
   * - ``isdecimal``
     - |numpy|
     - |arrow|
   * - ``isdigit``
     - |numpy|
     - |arrow|
   * - ``islower``
     - |numpy|
     - |arrow|
   * - ``isnumeric``
     - |numpy|
     - |arrow|
   * - ``isspace``
     - |numpy|
     - |arrow|
   * - ``istitle``
     - |numpy|
     - |arrow|
   * - ``isupper``
     - |numpy|
     - |arrow|
   * - ``join``
     - |arrow|
     - |elementwise|
   * - ``len``
     - |numpy|
     - |arrow|
   * - ``ljust``
     - |arrow|
     - |arrow|
   * - ``lower``
     - |arrow|
     - |arrow|
   * - ``lstrip``
     - |arrow|
     - |arrow|
   * - ``match(flags=0)``
     - |numpy|
     - |error|
   * - ``match(flags=re.I)``
     - |numpy|
     - |error|
   * - ``normalize``
     - |arrow|
     - |elementwise|
   * - ``pad(side=both)``
     - |arrow|
     - |arrow|
   * - ``pad(side=left)``
     - |arrow|
     - |arrow|
   * - ``pad(side=right)``
     - |arrow|
     - |arrow|
   * - ``partition``
     - |arrow|
     - |elementwise|
   * - ``removeprefix``
     - |arrow|
     - |arrow|
   * - ``removesuffix``
     - |arrow|
     - |arrow|
   * - ``repeat``
     - |arrow|
     - |arrow|
   * - ``replace(case=False)``
     - |arrow|
     - |notimpl|
   * - ``replace(case=True)``
     - |arrow|
     - |arrow|
   * - ``replace(repl=callable)``
     - |arrow|
     - |notimpl|
   * - ``rfind``
     - |numpy|
     - |elementwise|
   * - ``rindex``
     - |numpy|
     - |elementwise|
   * - ``rjust``
     - |arrow|
     - |arrow|
   * - ``rpartition``
     - |arrow|
     - |elementwise|
   * - ``rsplit(expand=False)``
     - |object|
     - |arrow|
   * - ``rsplit(expand=True)``
     - |arrow|
     - |arrow|
   * - ``rstrip``
     - |arrow|
     - |arrow|
   * - ``slice``
     - |arrow|
     - |arrow|
   * - ``slice_replace``
     - |arrow|
     - |arrow|
   * - ``split(expand=False)``
     - |object|
     - |arrow|
   * - ``split(expand=True)``
     - |arrow|
     - |arrow|
   * - ``startswith``
     - |numpy|
     - |arrow|
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
     - |arrow|
     - |elementwise|
   * - ``upper``
     - |arrow|
     - |arrow|
   * - ``wrap``
     - |arrow|
     - |elementwise|
   * - ``zfill``
     - |arrow|
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
   :widths: 25 10 10 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Method
     - string
     - large_string
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
   * - ``clip``
     -
     -
     - |arrow|
     - |arrow|
     -
     - |arrow|
     - |typeerror|
     - |arrow|
     - |typeerror|
     -
   * - ``cummax``
     - |arrow|
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
     - |arrow|
     - |typeerror|
     - |arrow|
     - |arrow|
     - |arrow|
     - |arrow|
     - |typeerror|
   * - ``cumprod``
     -
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
   * - ``fillna(limit=1)``
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
   * - ``fillna(limit=None)``
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
   * - ``interpolate(method=linear)``
     -
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
     -
     - |notimpl|
     - |notimpl|
     -
     - |notimpl|
     - |notimpl|
     - |notimpl|
     - |notimpl|
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
   * - ``round``
     -
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
     - |numpy|
     - |numpy|
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

Algorithmic Methods
===================

These methods have behavior that may vary across versions.
Explicit kwargs are used for deterministic results.

.. list-table::
   :widths: 25 10 10 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Method
     - string
     - large_string
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
     - |numpy|
   * - ``rank``
     - |numpy|
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
