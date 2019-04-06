{{ header }}

.. _api.arrays:

=============
Pandas Arrays
=============

.. currentmodule:: pandas

For most data types, pandas uses NumPy arrays as the concrete
objects contained with a :class:`Index`, :class:`Series`, or
:class:`DataFrame`.

For some data types, pandas extends NumPy's type system.

=================== ========================= ================== =============================
Kind of Data        Pandas Data Type          Scalar             Array
=================== ========================= ================== =============================
TZ-aware datetime   :class:`DatetimeTZDtype`  :class:`Timestamp` :ref:`api.arrays.datetime`
Timedeltas          (none)                    :class:`Timedelta` :ref:`api.arrays.timedelta`
Period (time spans) :class:`PeriodDtype`      :class:`Period`    :ref:`api.arrays.period`
Intervals           :class:`IntervalDtype`    :class:`Interval`  :ref:`api.arrays.interval`
Nullable Integer    :class:`Int64Dtype`, ...  (none)             :ref:`api.arrays.integer_na`
Categorical         :class:`CategoricalDtype` (none)             :ref:`api.arrays.categorical`
Sparse              :class:`SparseDtype`      (none)             :ref:`api.arrays.sparse`
=================== ========================= ================== =============================

Pandas and third-party libraries can extend NumPy's type system (see :ref:`extending.extension-types`).
The top-level :meth:`array` method can be used to create a new array, which may be
stored in a :class:`Series`, :class:`Index`, or as a column in a :class:`DataFrame`.

.. autosummary::
   :toctree: api/

   array

.. _api.arrays.datetime:

Datetime Data
-------------

NumPy cannot natively represent timezone-aware datetimes. Pandas supports this
with the :class:`arrays.DatetimeArray` extension array, which can hold timezone-naive
or timezone-aware values.

:class:`Timestamp`, a subclass of :class:`datetime.datetime`, is pandas'
scalar type for timezone-naive or timezone-aware datetime data.

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
   Timestamp.dayofyear
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
   Timestamp.value
   Timestamp.week
   Timestamp.weekofyear
   Timestamp.year

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

   Timestamp.astimezone
   Timestamp.ceil
   Timestamp.combine
   Timestamp.ctime
   Timestamp.date
   Timestamp.day_name
   Timestamp.dst
   Timestamp.floor
   Timestamp.freq
   Timestamp.freqstr
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
For timezone-aware data, the ``.dtype`` of a ``DatetimeArray`` is a
:class:`DatetimeTZDtype`. For timezone-naive data, ``np.dtype("datetime64[ns]")``
is used.

If the data are tz-aware, then every value in the array must have the same timezone.

.. autosummary::
   :toctree: api/

   arrays.DatetimeArray
   DatetimeTZDtype

.. _api.arrays.timedelta:

Timedelta Data
--------------

NumPy can natively represent timedeltas. Pandas provides :class:`Timedelta`
for symmetry with :class:`Timestamp`.

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
   Timedelta.delta
   Timedelta.freq
   Timedelta.is_populated
   Timedelta.max
   Timedelta.microseconds
   Timedelta.min
   Timedelta.nanoseconds
   Timedelta.resolution
   Timedelta.seconds
   Timedelta.value
   Timedelta.view

Methods
~~~~~~~
.. autosummary::
   :toctree: api/

   Timedelta.ceil
   Timedelta.floor
   Timedelta.isoformat
   Timedelta.round
   Timedelta.to_pytimedelta
   Timedelta.to_timedelta64
   Timedelta.to_numpy
   Timedelta.total_seconds

A collection of timedeltas may be stored in a :class:`TimedeltaArray`.

.. autosummary::
   :toctree: api/

   arrays.TimedeltaArray

.. _api.arrays.period:

Timespan Data
-------------

Pandas represents spans of times as :class:`Period` objects.

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
   Period.dayofyear
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

A collection of timedeltas may be stored in a :class:`arrays.PeriodArray`.
Every period in a ``PeriodArray`` must have the same ``freq``.

.. autosummary::
   :toctree: api/

   arrays.PeriodArray
   PeriodDtype

.. _api.arrays.interval:

Interval Data
-------------

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

   arrays.IntervalArray
   IntervalDtype

.. _api.arrays.integer_na:

Nullable Integer
----------------

:class:`numpy.ndarray` cannot natively represent integer-data with missing values.
Pandas provides this through :class:`arrays.IntegerArray`.

.. autosummary::
   :toctree: api/

   arrays.IntegerArray
   Int8Dtype
   Int16Dtype
   Int32Dtype
   Int64Dtype
   UInt8Dtype
   UInt16Dtype
   UInt32Dtype
   UInt64Dtype

.. _api.arrays.categorical:

Categorical Data
----------------

Pandas defines a custom data type for representing data that can take only a
limited, fixed set of values. The dtype of a ``Categorical`` can be described by
a :class:`pandas.api.types.CategoricalDtype`.

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

The dtype information is available on the ``Categorical``

.. autosummary::
   :toctree: api/

   Categorical.dtype
   Categorical.categories
   Categorical.ordered
   Categorical.codes

``np.asarray(categorical)`` works by implementing the array interface. Be aware, that this converts
the Categorical back to a NumPy array, so categories and order information is not preserved!

.. autosummary::
   :toctree: api/

   Categorical.__array__

A ``Categorical`` can be stored in a ``Series`` or ``DataFrame``.
To create a Series of dtype ``category``, use ``cat = s.astype(dtype)`` or
``Series(..., dtype=dtype)`` where ``dtype`` is either

* the string ``'category'``
* an instance of :class:`~pandas.api.types.CategoricalDtype`.

If the Series is of dtype ``CategoricalDtype``, ``Series.cat`` can be used to change the categorical
data. See :ref:`api.series.cat` for more.

.. _api.arrays.sparse:

Sparse Data
-----------

Data where a single value is repeated many times (e.g. ``0`` or ``NaN``) may
be stored efficiently as a :class:`SparseArray`.

.. autosummary::
   :toctree: api/

   SparseArray
   SparseDtype

The ``Series.sparse`` accessor may be used to access sparse-specific attributes
and methods if the :class:`Series` contains sparse values. See
:ref:`api.series.sparse` for more.
