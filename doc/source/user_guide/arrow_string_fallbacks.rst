.. _arrow-fallbacks:

{{ header }}

******************************
Arrow Method Support Reference
******************************

This document provides a comprehensive reference of which pandas methods
use native PyArrow compute functions vs. falling back to Python/NumPy
implementations when using Arrow-backed arrays.

Minimum PyArrow version: **13.0.0**

Legend
======

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Status
     - Description
   * - |arrow|
     - Uses native PyArrow compute kernel
   * - |conditional|
     - Uses Arrow by default, falls back to Python for certain parameters
   * - |version|
     - Behavior depends on PyArrow version
   * - |pyfallback|
     - Python fallback (combines |elementwise| + |object| in summary table)
   * - |elementwise|
     - Processes elements one-by-one in Python (not vectorized)
   * - |object|
     - Converts to Python objects, then processes element-by-element
   * - |notimpl|
     - Not implemented for Arrow arrays

.. |arrow| replace:: Arrow
.. |conditional| replace:: Conditional
.. |version| replace:: Version-gated
.. |pyfallback| replace:: Py fallback
.. |elementwise| replace:: Element-wise
.. |object| replace:: Object fallback
.. |notimpl| replace:: Not implemented


Summary
=======

.. list-table::
   :widths: 30 12 12 12 12 12 10
   :header-rows: 1

   * - Category
     - |arrow|
     - |conditional|
     - |version|
     - |pyfallback|
     - |notimpl|
     - Total
   * - String methods
     - 25
     - 4
     - 2
     - 24
     - 0
     - 55
   * - Arithmetic & comparison
     - 15
     - 6
     - 0
     - 0
     - 4
     - 25
   * - Datetime methods
     - 35
     - 0
     - 0
     - 9
     - 0
     - 44
   * - Aggregation methods
     - 9
     - 2
     - 1
     - 0
     - 1
     - 13
   * - Array methods
     - 9
     - 7
     - 0
     - 3
     - 0
     - 19
   * - List accessor
     - 2
     - 1
     - 0
     - 0
     - 0
     - 3
   * - Struct accessor
     - 3
     - 0
     - 0
     - 0
     - 0
     - 3
   * - **Total**
     - **98**
     - **20**
     - **3**
     - **36**
     - **5**
     - **162**

**98** methods (60%) use native Arrow.
**20** methods (12%) use Arrow with conditional fallbacks.
Overall, **73%** have Arrow-native code paths.


String Methods (Series.str.*)
=============================

.. list-table::
   :widths: 20 12 25 43
   :header-rows: 1

   * - Method
     - Status
     - Arrow Function
     - Notes
   * - ``str.accumulate``
     - |object|
     - ``pc.all``
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.capitalize``
     - |arrow|
     - ``pc.utf8_capitalize``
     - -
   * - ``str.casefold``
     - |object|
     - -
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.contains``
     - |conditional|
     - -
     - Falls back: case-insensitive matching; compiled regex pattern (+1 more)
   * - ``str.count``
     - |conditional|
     - ``pc.count_substring_regex``
     - Falls back: regex flags used
   * - ``str.encode``
     - |object|
     - -
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.endswith``
     - |arrow|
     - ``pc.ends_with``
     - -
   * - ``str.extract``
     - |object|
     - ``pc.extract_regex``
     - Falls back: regex flags used | Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.find``
     - |conditional|
     - ``pc.find_substring``
     - Falls back: start/end parameters with edge cases | _apply_elementwise() fallback; Related: GH#56792
   * - ``str.find_``
     - |object|
     - -
     - Object fallback
   * - ``str.findall``
     - |object|
     - -
     - Falls back: regex flags used | Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.fullmatch``
     - |object|
     - -
     - Falls back: compiled regex pattern; regex flags used | via ObjectStringArrayMixin
   * - ``str.get``
     - |arrow|
     - ``pc.utf8_length``
     - -
   * - ``str.get_dummies``
     - |object|
     - ``pc.split_pattern``
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.getitem``
     - |object|
     - -
     - Object fallback
   * - ``str.index``
     - |object|
     - -
     - Falls back: start/end parameters with edge cases | Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.isalnum``
     - |arrow|
     - ``pc.utf8_is_alnum``
     - -
   * - ``str.isalpha``
     - |arrow|
     - ``pc.utf8_is_alpha``
     - -
   * - ``str.isascii``
     - |arrow|
     - ``pc.string_is_ascii``
     - -
   * - ``str.isdecimal``
     - |arrow|
     - ``pc.utf8_is_decimal``
     - -
   * - ``str.isdigit``
     - |version|
     - ``pc.utf8_is_digit``
     - Requires PyArrow >= 21.0.0 | Falls back: PyArrow < 21.0.0
   * - ``str.islower``
     - |arrow|
     - ``pc.utf8_is_lower``
     - -
   * - ``str.isnumeric``
     - |arrow|
     - ``pc.utf8_is_numeric``
     - -
   * - ``str.isspace``
     - |arrow|
     - ``pc.utf8_is_space``
     - -
   * - ``str.istitle``
     - |arrow|
     - ``pc.utf8_is_title``
     - -
   * - ``str.isupper``
     - |arrow|
     - ``pc.utf8_is_upper``
     - -
   * - ``str.join``
     - |object|
     - ``pc.binary_join``
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.len``
     - |arrow|
     - ``pc.utf8_length``
     - -
   * - ``str.lower``
     - |arrow|
     - ``pc.utf8_lower``
     - -
   * - ``str.lstrip``
     - |arrow|
     - ``pc.utf8_ltrim_whitespace``
     - -
   * - ``str.map``
     - |object|
     - -
     - Object fallback
   * - ``str.match``
     - |object|
     - -
     - Falls back: compiled regex pattern; regex flags used | via ObjectStringArrayMixin
   * - ``str.normalize``
     - |object|
     - -
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.pad/center/ljust/rjust``
     - |version|
     - -
     - Requires PyArrow >= 17.0.0 | Falls back: PyArrow < 17.0.0 | Converts to object dtype; Related: GH#59624
   * - ``str.partition``
     - |object|
     - -
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.removeprefix``
     - |arrow|
     - ``pc.starts_with``
     - -
   * - ``str.removesuffix``
     - |arrow|
     - ``pc.ends_with``
     - -
   * - ``str.repeat``
     - |object|
     - -
     - Falls back: array-like repeats | via ObjectStringArrayMixin
   * - ``str.replace``
     - |conditional|
     - -
     - Falls back: callable replacement; case-insensitive matching (+3 more) | Raises: replace is not supported with a re.Pattern, cal...
   * - ``str.rfind``
     - |object|
     - -
     - Falls back: start/end parameters with edge cases | Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.rindex``
     - |object|
     - -
     - Falls back: start/end parameters with edge cases | Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.rpartition``
     - |object|
     - -
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.rsplit``
     - |object|
     - ``pc.utf8_split_whitespace``
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.rstrip``
     - |arrow|
     - ``pc.utf8_rtrim_whitespace``
     - -
   * - ``str.slice``
     - |arrow|
     - ``pc.utf8_slice_codeunits``
     - Related: GH#59710
   * - ``str.slice_replace``
     - |arrow|
     - ``pc.utf8_replace_slice``
     - -
   * - ``str.split``
     - |object|
     - -
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.startswith``
     - |arrow|
     - ``pc.starts_with``
     - -
   * - ``str.strip``
     - |arrow|
     - ``pc.utf8_trim_whitespace``
     - -
   * - ``str.swapcase``
     - |arrow|
     - ``pc.utf8_swapcase``
     - -
   * - ``str.title``
     - |arrow|
     - ``pc.utf8_title``
     - -
   * - ``str.translate``
     - |object|
     - -
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.upper``
     - |arrow|
     - ``pc.utf8_upper``
     - -
   * - ``str.wrap``
     - |object|
     - -
     - Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()
   * - ``str.zfill``
     - |object|
     - ``pc.utf8_zfill``
     - Falls back: PyArrow < 21.0.0 | Falls back for dtype='string[pyarrow]'; ArrowDtype uses _apply_elementwise()

Arithmetic & Comparison Operations
==================================

.. list-table::
   :widths: 20 12 25 43
   :header-rows: 1

   * - Method
     - Status
     - Arrow Function
     - Notes
   * - ``!=``
     - |conditional|
     - ``pc.not_equal``
     - Falls back: ArrowNotImplementedError for unsupported type combinations
   * - ``%``
     - |notimpl|
     - -
     - Not supported for Arrow arrays, raises TypeError
   * - ``% (right)``
     - |notimpl|
     - -
     - Not supported for Arrow arrays, raises TypeError
   * - ``&``
     - |arrow|
     - ``pc.and_kleene``
     - Uses Kleene logic for boolean, bit_wise for integers
   * - ``*``
     - |arrow|
     - ``pc.multiply_checked``
     - -
   * - ``* (right)``
     - |arrow|
     - ``pc.multiply_checked``
     - -
   * - ``**``
     - |arrow|
     - ``pc.power_checked``
     - -
   * - ``** (right)``
     - |arrow|
     - ``pc.power_checked``
     - -
   * - ``+``
     - |arrow|
     - ``pc.add_checked``
     - -
   * - ``+ (right)``
     - |arrow|
     - ``pc.add_checked``
     - -
   * - ``-``
     - |arrow|
     - ``pc.subtract_checked``
     - -
   * - ``- (right)``
     - |arrow|
     - ``pc.subtract_checked``
     - -
   * - ``/``
     - |arrow|
     - ``pc.divide``
     - -
   * - ``/ (right)``
     - |arrow|
     - ``pc.divide``
     - -
   * - ``//``
     - |arrow|
     - ``floordiv_compat``
     - -
   * - ``// (right)``
     - |arrow|
     - ``floordiv_compat``
     - -
   * - ``<``
     - |conditional|
     - ``pc.less``
     - Falls back: ArrowNotImplementedError for unsupported type combinations
   * - ``<=``
     - |conditional|
     - ``pc.less_equal``
     - Falls back: ArrowNotImplementedError for unsupported type combinations
   * - ``==``
     - |conditional|
     - ``pc.equal``
     - Falls back: ArrowNotImplementedError for unsupported type combinations
   * - ``>``
     - |conditional|
     - ``pc.greater``
     - Falls back: ArrowNotImplementedError for unsupported type combinations
   * - ``>=``
     - |conditional|
     - ``pc.greater_equal``
     - Falls back: ArrowNotImplementedError for unsupported type combinations
   * - ``^``
     - |arrow|
     - ``pc.xor``
     - Uses Kleene logic for boolean, bit_wise for integers
   * - ``divmod``
     - |notimpl|
     - -
     - Not supported for Arrow arrays, raises TypeError
   * - ``divmod (right)``
     - |notimpl|
     - -
     - Not supported for Arrow arrays, raises TypeError
   * - ``|``
     - |arrow|
     - ``pc.or_kleene``
     - Uses Kleene logic for boolean, bit_wise for integers

Datetime Methods (Series.dt.*)
==============================

.. list-table::
   :widths: 20 12 25 43
   :header-rows: 1

   * - Method
     - Status
     - Arrow Function
     - Notes
   * - ``dt.as_unit``
     - |arrow|
     - ``pc.cast``
     - -
   * - ``dt.ceil``
     - |arrow|
     - -
     - -
   * - ``dt.date``
     - |arrow|
     - ``ChunkedArray.cast``
     - -
   * - ``dt.day``
     - |arrow|
     - ``pc.day``
     - -
   * - ``dt.day_name``
     - |arrow|
     - ``pc.strftime``
     - -
   * - ``dt.day_of_week/dayofweek/weekday``
     - |arrow|
     - ``pc.day_of_week``
     - -
   * - ``dt.day_of_year/dayofyear``
     - |arrow|
     - ``pc.day_of_year``
     - -
   * - ``dt.days``
     - |object|
     - -
     - via to_numpy(); via TimedeltaArray; .components.days
   * - ``dt.days_in_month/daysinmonth``
     - |arrow|
     - ``pc.days_between``
     - -
   * - ``dt.floor``
     - |arrow|
     - -
     - -
   * - ``dt.hour``
     - |arrow|
     - ``pc.hour``
     - -
   * - ``dt.hours (td)``
     - |object|
     - -
     - via to_numpy(); via TimedeltaArray; .components.hours
   * - ``dt.is_leap_year``
     - |arrow|
     - ``pc.is_leap_year``
     - -
   * - ``dt.is_month_end``
     - |arrow|
     - ``pc.equal``
     - -
   * - ``dt.is_month_start``
     - |arrow|
     - ``pc.equal``
     - -
   * - ``dt.is_quarter_end``
     - |arrow|
     - ``pc.equal``
     - -
   * - ``dt.is_quarter_start``
     - |arrow|
     - ``pc.equal``
     - -
   * - ``dt.is_year_end``
     - |arrow|
     - ``pc.and_``
     - -
   * - ``dt.is_year_start``
     - |arrow|
     - ``pc.and_``
     - -
   * - ``dt.isocalendar``
     - |arrow|
     - ``pc.iso_calendar``
     - -
   * - ``dt.microsecond``
     - |arrow|
     - ``pc.microsecond``
     - -
   * - ``dt.microseconds (td)``
     - |object|
     - -
     - via to_numpy(); via TimedeltaArray; .components.microseconds
   * - ``dt.milliseconds``
     - |object|
     - -
     - via to_numpy(); via TimedeltaArray; .components.milliseconds
   * - ``dt.minute``
     - |arrow|
     - ``pc.minute``
     - -
   * - ``dt.minutes (td)``
     - |object|
     - -
     - via to_numpy(); via TimedeltaArray; .components.minutes
   * - ``dt.month``
     - |arrow|
     - ``pc.month``
     - -
   * - ``dt.month_name``
     - |arrow|
     - ``pc.strftime``
     - -
   * - ``dt.nanosecond``
     - |arrow|
     - ``pc.nanosecond``
     - -
   * - ``dt.nanoseconds (td)``
     - |object|
     - -
     - via to_numpy(); via TimedeltaArray; .components.nanoseconds
   * - ``dt.normalize``
     - |arrow|
     - ``pc.floor_temporal``
     - -
   * - ``dt.quarter``
     - |arrow|
     - ``pc.quarter``
     - -
   * - ``dt.round``
     - |arrow|
     - -
     - -
   * - ``dt.second``
     - |arrow|
     - ``pc.second``
     - -
   * - ``dt.seconds (td)``
     - |object|
     - -
     - via to_numpy(); via TimedeltaArray; .components.seconds
   * - ``dt.strftime``
     - |arrow|
     - ``pc.strftime``
     - -
   * - ``dt.time``
     - |arrow|
     - ``ChunkedArray.cast``
     - -
   * - ``dt.to_pydatetime``
     - |object|
     - -
     - via to_numpy(); via pylist
   * - ``dt.to_pytimedelta``
     - |object|
     - -
     - via to_numpy(); via pylist
   * - ``dt.total_seconds``
     - |arrow|
     - ``pc.divide``
     - -
   * - ``dt.tz``
     - |arrow|
     - -
     - -
   * - ``dt.tz_convert``
     - |arrow|
     - ``ChunkedArray.cast``
     - -
   * - ``dt.tz_localize``
     - |arrow|
     - ``pc.local_timestamp``
     - -
   * - ``dt.unit``
     - |arrow|
     - -
     - -
   * - ``dt.year``
     - |arrow|
     - ``pc.year``
     - -

Series Aggregation Methods (Series.sum(), etc.)
===============================================

.. list-table::
   :widths: 20 12 25 43
   :header-rows: 1

   * - Method
     - Status
     - Arrow Function
     - Notes
   * - ``all``
     - |conditional|
     - ``pc.all``
     - Series-level only; converts non-bool to bool first; special handling for temporal types
   * - ``any``
     - |conditional|
     - ``pc.any``
     - Series-level only; converts non-bool to bool first; special handling for temporal types
   * - ``kurt``
     - |notimpl|
     - -
     - Not supported for Arrow arrays
   * - ``max``
     - |arrow|
     - ``pc.max``
     - Series-level only; special handling for temporal types
   * - ``mean``
     - |arrow|
     - ``pc.mean``
     - Series-level only; special handling for temporal types
   * - ``median``
     - |arrow|
     - ``pc.quantile``
     - Series-level only; special handling for temporal types
   * - ``min``
     - |arrow|
     - ``pc.min``
     - Series-level only; special handling for temporal types
   * - ``prod``
     - |arrow|
     - ``pc.product``
     - Series-level only; special handling for temporal types
   * - ``sem``
     - |arrow|
     - ``pc.stddev / pc.sqrt_checked / pc.count / pc.divide_checked``
     - Series-level only; computed from stddev, sqrt_checked, count, divide_checked; special handling for temporal types
   * - ``skew``
     - |version|
     - ``pc.skew``
     - Series-level only; requires PyArrow >= 20.0
   * - ``std``
     - |arrow|
     - ``pc.stddev``
     - Series-level only; special handling for temporal types
   * - ``sum``
     - |arrow|
     - ``pc.sum``
     - Series-level only; special handling for temporal types
   * - ``var``
     - |arrow|
     - ``pc.variance``
     - Series-level only; special handling for temporal types

.. note::

   The aggregation methods above apply to **Series-level** operations
   only (e.g., ``ser.sum()``). As of this pandas version, **GroupBy**
   **aggregations** (``df.groupby(...).sum()``) use the standard
   pandas implementation rather than Arrow compute kernels. This may
   change in future versions.


General Array Methods
=====================

.. list-table::
   :widths: 20 12 25 43
   :header-rows: 1

   * - Method
     - Status
     - Arrow Function
     - Notes
   * - ``argsort()``
     - |arrow|
     - ``pc.array_sort_indices``
     - -
   * - ``bfill() / backfill()``
     - |conditional|
     - ``pc.fill_null_backward``
     - Falls back: limit parameter specified; limit_area parameter specified | Falls back if limit parameter specified or limit_area parameter specified
   * - ``copy()``
     - |arrow|
     - -
     - -
   * - ``cummax()``
     - |conditional|
     - ``pc.cumulative_max``
     - String dtype uses NumPy fallback
   * - ``cummin()``
     - |conditional|
     - ``pc.cumulative_min``
     - String dtype uses NumPy fallback
   * - ``cumprod()``
     - |conditional|
     - ``pc.cumulative_prod_checked``
     - Not supported for string dtype
   * - ``cumsum()``
     - |conditional|
     - ``pc.cumulative_sum_checked``
     - String dtype uses NumPy fallback
   * - ``dropna()``
     - |arrow|
     - ``pc.drop_null``
     - -
   * - ``duplicated()``
     - |object|
     - -
     - Falls back: temporal types have special handling | via to_numpy()
   * - ``factorize()``
     - |arrow|
     - ``pc.fill_null``
     - -
   * - ``ffill() / pad()``
     - |conditional|
     - ``pc.fill_null_forward``
     - Falls back: limit parameter specified; limit_area parameter specified | Falls back if limit parameter specified or limit_area parameter specified
   * - ``fillna(value)``
     - |conditional|
     - ``pc.fill_null``
     - Falls back: falls back to base implementation for some cases; limit parameter specified
   * - ``interpolate()``
     - |arrow|
     - ``pc.fill_null_backward / pc.pairwise_diff_checked``
     - Falls back: limit_area parameter specified
   * - ``isna() / isnull()``
     - |object|
     - -
     - via to_numpy(); fast null path
   * - ``round(decimals)``
     - |arrow|
     - ``pc.round``
     - -
   * - ``searchsorted(value)``
     - |object|
     - -
     - Falls back: duration types have special handling | via to_numpy()
   * - ``take(indices)``
     - |arrow|
     - ``pc.fill_null``
     - -
   * - ``unique()``
     - |arrow|
     - ``pc.unique``
     - -
   * - ``value_counts()``
     - |arrow|
     - ``ChunkedArray.value_counts``
     - -

List Accessor Methods (Series.list.*)
=====================================

.. list-table::
   :widths: 20 12 25 43
   :header-rows: 1

   * - Method
     - Status
     - Arrow Function
     - Notes
   * - ``list.[index] / [slice]``
     - |conditional|
     - ``pc.add / pc.list_value_length``
     - Negative indices not supported; Different handling for int vs slice
   * - ``list.flatten()``
     - |arrow|
     - ``pa.compute.list_flatten``
     - -
   * - ``list.len()``
     - |arrow|
     - ``pc.list_value_length``
     - -

Struct Accessor Methods (Series.struct.*)
=========================================

.. list-table::
   :widths: 20 12 25 43
   :header-rows: 1

   * - Method
     - Status
     - Arrow Function
     - Notes
   * - ``struct.dtypes (property)``
     - |arrow|
     - -
     - Returns dtype info (property)
   * - ``struct.explode()``
     - |arrow|
     - -
     - Returns all fields as DataFrame
   * - ``struct.field(name_or_index)``
     - |arrow|
     - ``pc.field / pc.struct_field``
     - Supports nested field access via list
