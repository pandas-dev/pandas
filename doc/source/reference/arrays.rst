{{ header }}

.. _api.arrays:

======================================
pandas arrays, scalars, and data types
======================================

*******
Objects
*******

.. currentmodule:: pandas

For most data types, pandas uses NumPy arrays as the concrete
objects contained with a :class:`Index`, :class:`Series`, or
:class:`DataFrame`.

For some data types, pandas extends NumPy's type system. String aliases for these types
can be found at :ref:`basics.dtypes`.

=================== ========================== ============================= =============================
Kind of Data        pandas Data Type           Scalar                        Array
=================== ========================== ============================= =============================
TZ-aware datetime   :class:`DatetimeTZDtype`   :class:`Timestamp`            :ref:`api.arrays.datetime`
Timedeltas          (none)                     :class:`Timedelta`            :ref:`api.arrays.timedelta`
Period (time spans) :class:`PeriodDtype`       :class:`Period`               :ref:`api.arrays.period`
Intervals           :class:`IntervalDtype`     :class:`Interval`             :ref:`api.arrays.interval`
Nullable Integer    :class:`Int64Dtype`, ...   (none)                        :ref:`api.arrays.integer_na`
Nullable Float      :class:`Float64Dtype`, ... (none)                        :ref:`api.arrays.float_na`
Categorical         :class:`CategoricalDtype`  (none)                        :ref:`api.arrays.categorical`
Sparse              :class:`SparseDtype`       (none)                        :ref:`api.arrays.sparse`
Strings             :class:`StringDtype`       :class:`str`                  :ref:`api.arrays.string`
Nullable Boolean    :class:`BooleanDtype`      :class:`bool`                 :ref:`api.arrays.bool`
PyArrow             :class:`ArrowDtype`        Python Scalars or :class:`NA` :ref:`api.arrays.arrow`
=================== ========================== ============================= =============================

pandas and third-party libraries can extend NumPy's type system (see :ref:`extending.extension-types`).
The top-level :meth:`array` method can be used to create a new array, which may be
stored in a :class:`Series`, :class:`Index`, or as a column in a :class:`DataFrame`.

.. autosummary::
   :toctree: api/

   array

.. _api.arrays.arrow:

PyArrow
-------

.. warning::

    This feature is experimental, and the API can change in a future release without warning.

The :class:`arrays.ArrowExtensionArray` is backed by a :external+pyarrow:py:class:`pyarrow.ChunkedArray` with a
:external+pyarrow:py:class:`pyarrow.DataType` instead of a NumPy array and data type. The ``.dtype`` of a :class:`arrays.ArrowExtensionArray`
is an :class:`ArrowDtype`.

`Pyarrow <https://arrow.apache.org/docs/python/index.html>`__ provides similar array and `data type <https://arrow.apache.org/docs/python/api/datatypes.html>`__
support as NumPy including first-class nullability support for all data types, immutability and more.

The table below shows the equivalent pyarrow-backed (``pa``), pandas extension, and numpy (``np``) types that are recognized by pandas.
Pyarrow-backed types below need to be passed into :class:`ArrowDtype` to be recognized by pandas e.g. ``pd.ArrowDtype(pa.bool_())``

=============================================== ========================== ===================
PyArrow type                                    pandas extension type      NumPy type
=============================================== ========================== ===================
:external+pyarrow:py:func:`pyarrow.bool_`       :class:`BooleanDtype`      ``np.bool_``
:external+pyarrow:py:func:`pyarrow.int8`        :class:`Int8Dtype`         ``np.int8``
:external+pyarrow:py:func:`pyarrow.int16`       :class:`Int16Dtype`        ``np.int16``
:external+pyarrow:py:func:`pyarrow.int32`       :class:`Int32Dtype`        ``np.int32``
:external+pyarrow:py:func:`pyarrow.int64`       :class:`Int64Dtype`        ``np.int64``
:external+pyarrow:py:func:`pyarrow.uint8`       :class:`UInt8Dtype`        ``np.uint8``
:external+pyarrow:py:func:`pyarrow.uint16`      :class:`UInt16Dtype`       ``np.uint16``
:external+pyarrow:py:func:`pyarrow.uint32`      :class:`UInt32Dtype`       ``np.uint32``
:external+pyarrow:py:func:`pyarrow.uint64`      :class:`UInt64Dtype`       ``np.uint64``
:external+pyarrow:py:func:`pyarrow.float32`     :class:`Float32Dtype`      ``np.float32``
:external+pyarrow:py:func:`pyarrow.float64`     :class:`Float64Dtype`      ``np.float64``
:external+pyarrow:py:func:`pyarrow.time32`      (none)                     (none)
:external+pyarrow:py:func:`pyarrow.time64`      (none)                     (none)
:external+pyarrow:py:func:`pyarrow.timestamp`   :class:`DatetimeTZDtype`   ``np.datetime64``
:external+pyarrow:py:func:`pyarrow.date32`      (none)                     (none)
:external+pyarrow:py:func:`pyarrow.date64`      (none)                     (none)
:external+pyarrow:py:func:`pyarrow.duration`    (none)                     ``np.timedelta64``
:external+pyarrow:py:func:`pyarrow.binary`      (none)                     (none)
:external+pyarrow:py:func:`pyarrow.string`      :class:`StringDtype`       ``np.str_``
:external+pyarrow:py:func:`pyarrow.decimal128`  (none)                     (none)
:external+pyarrow:py:func:`pyarrow.list_`       (none)                     (none)
:external+pyarrow:py:func:`pyarrow.map_`        (none)                     (none)
:external+pyarrow:py:func:`pyarrow.dictionary`  :class:`CategoricalDtype`  (none)
=============================================== ========================== ===================

.. note::

    Pyarrow-backed string support is provided by both ``pd.StringDtype("pyarrow")`` and ``pd.ArrowDtype(pa.string())``.
    ``pd.StringDtype("pyarrow")`` is described below in the :ref:`string section <api.arrays.string>`
    and will be returned if the string alias ``"string[pyarrow]"`` is specified. ``pd.ArrowDtype(pa.string())``
    generally has better interoperability with :class:`ArrowDtype` of different types.

While individual values in an :class:`arrays.ArrowExtensionArray` are stored as a PyArrow objects, scalars are **returned**
as Python scalars corresponding to the data type, e.g. a PyArrow int64 will be returned as Python int, or :class:`NA` for missing
values.

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   arrays.ArrowExtensionArray

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   ArrowDtype

For more information, please see the :ref:`PyArrow user guide <pyarrow>`

.. _api.arrays.datetime:

Datetimes
---------

NumPy cannot natively represent timezone-aware datetimes. pandas supports this
with the :class:`arrays.DatetimeArray` extension array, which can hold timezone-naive
or timezone-aware values.

:class:`Timestamp`, a subclass of :class:`datetime.datetime`, is pandas'
scalar type for timezone-naive or timezone-aware datetime data. :class:`NaT`
is the missing value for datetime data.

.. autosummary::
   :toctree: api/

   Timestamp

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Timestamp.asm8
   Timestamp.day
   Timestamp.dayofweek
   Timestamp.day_of_week
   Timestamp.dayofyear
   Timestamp.day_of_year
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
   Timestamp.unit
   Timestamp.value
   Timestamp.week
   Timestamp.weekofyear
   Timestamp.year

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

   Timestamp.as_unit
   Timestamp.astimezone
   Timestamp.ceil
   Timestamp.combine
   Timestamp.ctime
   Timestamp.date
   Timestamp.day_name
   Timestamp.dst
   Timestamp.floor
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
   Timestamp.to_numpy
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

A collection of timestamps may be stored in a :class:`arrays.DatetimeArray`.
For timezone-aware data, the ``.dtype`` of a :class:`arrays.DatetimeArray` is a
:class:`DatetimeTZDtype`. For timezone-naive data, ``np.dtype("datetime64[ns]")``
is used.

If the data are timezone-aware, then every value in the array must have the same timezone.

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   arrays.DatetimeArray

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   DatetimeTZDtype

.. _api.arrays.timedelta:

Timedeltas
----------

NumPy can natively represent timedeltas. pandas provides :class:`Timedelta`
for symmetry with :class:`Timestamp`. :class:`NaT`
is the missing value for timedelta data.

.. autosummary::
   :toctree: api/

   Timedelta

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Timedelta.asm8
   Timedelta.components
   Timedelta.days
   Timedelta.max
   Timedelta.microseconds
   Timedelta.min
   Timedelta.nanoseconds
   Timedelta.resolution
   Timedelta.seconds
   Timedelta.unit
   Timedelta.value
   Timedelta.view

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

   Timedelta.as_unit
   Timedelta.ceil
   Timedelta.floor
   Timedelta.isoformat
   Timedelta.round
   Timedelta.to_pytimedelta
   Timedelta.to_timedelta64
   Timedelta.to_numpy
   Timedelta.total_seconds

A collection of :class:`Timedelta` may be stored in a :class:`TimedeltaArray`.

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   arrays.TimedeltaArray

.. _api.arrays.period:

Periods
-------

pandas represents spans of times as :class:`Period` objects.

Period
------
.. autosummary::
   :toctree: api/

   Period

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Period.day
   Period.dayofweek
   Period.day_of_week
   Period.dayofyear
   Period.day_of_year
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
   :toctree: api/

   Period.asfreq
   Period.now
   Period.strftime
   Period.to_timestamp

A collection of :class:`Period` may be stored in a :class:`arrays.PeriodArray`.
Every period in a :class:`arrays.PeriodArray` must have the same ``freq``.

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   arrays.PeriodArray

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   PeriodDtype

.. _api.arrays.interval:

Intervals
---------

Arbitrary intervals can be represented as :class:`Interval` objects.

.. autosummary::
   :toctree: api/

    Interval

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Interval.closed
   Interval.closed_left
   Interval.closed_right
   Interval.is_empty
   Interval.left
   Interval.length
   Interval.mid
   Interval.open_left
   Interval.open_right
   Interval.overlaps
   Interval.right

A collection of intervals may be stored in an :class:`arrays.IntervalArray`.

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   arrays.IntervalArray

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   IntervalDtype


.. Those attributes and methods are included in the API because the docstrings
.. of IntervalIndex and IntervalArray are shared. Including it here to make
.. sure a docstring page is built for them to avoid warnings

..
    .. autosummary::
      :toctree: api/

      arrays.IntervalArray.left
      arrays.IntervalArray.right
      arrays.IntervalArray.closed
      arrays.IntervalArray.mid
      arrays.IntervalArray.length
      arrays.IntervalArray.is_empty
      arrays.IntervalArray.is_non_overlapping_monotonic
      arrays.IntervalArray.from_arrays
      arrays.IntervalArray.from_tuples
      arrays.IntervalArray.from_breaks
      arrays.IntervalArray.contains
      arrays.IntervalArray.overlaps
      arrays.IntervalArray.set_closed
      arrays.IntervalArray.to_tuples


.. _api.arrays.integer_na:

Nullable integer
----------------

:class:`numpy.ndarray` cannot natively represent integer-data with missing values.
pandas provides this through :class:`arrays.IntegerArray`.

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   arrays.IntegerArray

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   Int8Dtype
   Int16Dtype
   Int32Dtype
   Int64Dtype
   UInt8Dtype
   UInt16Dtype
   UInt32Dtype
   UInt64Dtype

.. _api.arrays.float_na:

Nullable float
--------------

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   arrays.FloatingArray

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   Float32Dtype
   Float64Dtype

.. _api.arrays.categorical:

Categoricals
------------

pandas defines a custom data type for representing data that can take only a
limited, fixed set of values. The dtype of a :class:`Categorical` can be described by
a :class:`CategoricalDtype`.

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   CategoricalDtype

.. autosummary::
   :toctree: api/

   CategoricalDtype.categories
   CategoricalDtype.ordered

Categorical data can be stored in a :class:`pandas.Categorical`

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   Categorical

The alternative :meth:`Categorical.from_codes` constructor can be used when you
have the categories and integer codes already:

.. autosummary::
   :toctree: api/

   Categorical.from_codes

The dtype information is available on the :class:`Categorical`

.. autosummary::
   :toctree: api/

   Categorical.dtype
   Categorical.categories
   Categorical.ordered
   Categorical.codes

``np.asarray(categorical)`` works by implementing the array interface. Be aware, that this converts
the :class:`Categorical` back to a NumPy array, so categories and order information is not preserved!

.. autosummary::
   :toctree: api/

   Categorical.__array__

A :class:`Categorical` can be stored in a :class:`Series` or :class:`DataFrame`.
To create a Series of dtype ``category``, use ``cat = s.astype(dtype)`` or
``Series(..., dtype=dtype)`` where ``dtype`` is either

* the string ``'category'``
* an instance of :class:`CategoricalDtype`.

If the :class:`Series` is of dtype :class:`CategoricalDtype`, ``Series.cat`` can be used to change the categorical
data. See :ref:`api.series.cat` for more.

.. _api.arrays.sparse:

Sparse
------

Data where a single value is repeated many times (e.g. ``0`` or ``NaN``) may
be stored efficiently as a :class:`arrays.SparseArray`.

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   arrays.SparseArray

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   SparseDtype

The ``Series.sparse`` accessor may be used to access sparse-specific attributes
and methods if the :class:`Series` contains sparse values. See
:ref:`api.series.sparse` and :ref:`the user guide <sparse>` for more.


.. _api.arrays.string:

Strings
-------

When working with text data, where each valid element is a string or missing,
we recommend using :class:`StringDtype` (with the alias ``"string"``).

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   arrays.StringArray
   arrays.ArrowStringArray

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   StringDtype

The ``Series.str`` accessor is available for :class:`Series` backed by a :class:`arrays.StringArray`.
See :ref:`api.series.str` for more.


.. _api.arrays.bool:

Nullable Boolean
----------------

The boolean dtype (with the alias ``"boolean"``) provides support for storing
boolean data (``True``, ``False``) with missing values, which is not possible
with a bool :class:`numpy.ndarray`.

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   arrays.BooleanArray

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   BooleanDtype


.. Dtype attributes which are manually listed in their docstrings: including
.. it here to make sure a docstring page is built for them

..
    .. autosummary::
      :toctree: api/

      DatetimeTZDtype.unit
      DatetimeTZDtype.tz
      PeriodDtype.freq
      IntervalDtype.subtype

*********
Utilities
*********

Constructors
------------
.. autosummary::
   :toctree: api/

   api.types.union_categoricals
   api.types.infer_dtype
   api.types.pandas_dtype

Data type introspection
~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

    api.types.is_any_real_numeric_dtype
    api.types.is_bool_dtype
    api.types.is_categorical_dtype
    api.types.is_complex_dtype
    api.types.is_datetime64_any_dtype
    api.types.is_datetime64_dtype
    api.types.is_datetime64_ns_dtype
    api.types.is_datetime64tz_dtype
    api.types.is_extension_array_dtype
    api.types.is_float_dtype
    api.types.is_int64_dtype
    api.types.is_integer_dtype
    api.types.is_interval_dtype
    api.types.is_numeric_dtype
    api.types.is_object_dtype
    api.types.is_period_dtype
    api.types.is_signed_integer_dtype
    api.types.is_string_dtype
    api.types.is_timedelta64_dtype
    api.types.is_timedelta64_ns_dtype
    api.types.is_unsigned_integer_dtype
    api.types.is_sparse

Iterable introspection
~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

    api.types.is_dict_like
    api.types.is_file_like
    api.types.is_list_like
    api.types.is_named_tuple
    api.types.is_iterator

Scalar introspection
~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

    api.types.is_bool
    api.types.is_complex
    api.types.is_float
    api.types.is_hashable
    api.types.is_integer
    api.types.is_interval
    api.types.is_number
    api.types.is_re
    api.types.is_re_compilable
    api.types.is_scalar
