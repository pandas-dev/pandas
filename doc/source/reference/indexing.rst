{{ header }}

.. _api.indexing:

=============
Index objects
=============

Index
-----
.. currentmodule:: pandas

**Many of these methods or variants thereof are available on the objects
that contain an index (Series/DataFrame) and those should most likely be
used before calling these methods directly.**

.. autosummary::
   :toctree: api/

   Index

Properties
~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Index.values
   Index.is_monotonic
   Index.is_monotonic_increasing
   Index.is_monotonic_decreasing
   Index.is_unique
   Index.has_duplicates
   Index.hasnans
   Index.dtype
   Index.inferred_type
   Index.is_all_dates
   Index.shape
   Index.name
   Index.names
   Index.nbytes
   Index.ndim
   Index.size
   Index.empty
   Index.T
   Index.memory_usage

Modifying and computations
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Index.all
   Index.any
   Index.argmin
   Index.argmax
   Index.copy
   Index.delete
   Index.drop
   Index.drop_duplicates
   Index.duplicated
   Index.equals
   Index.factorize
   Index.identical
   Index.insert
   Index.is_
   Index.is_boolean
   Index.is_categorical
   Index.is_floating
   Index.is_integer
   Index.is_interval
   Index.is_mixed
   Index.is_numeric
   Index.is_object
   Index.min
   Index.max
   Index.reindex
   Index.rename
   Index.repeat
   Index.where
   Index.take
   Index.putmask
   Index.unique
   Index.nunique
   Index.value_counts

Compatibility with MultiIndex
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Index.set_names
   Index.droplevel

Missing values
~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Index.fillna
   Index.dropna
   Index.isna
   Index.notna

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Index.astype
   Index.item
   Index.map
   Index.ravel
   Index.to_list
   Index.to_native_types
   Index.to_series
   Index.to_frame
   Index.view

Sorting
~~~~~~~
.. autosummary::
   :toctree: api/

   Index.argsort
   Index.searchsorted
   Index.sort_values

Time-specific operations
~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Index.shift

Combining / joining / set operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   Index.append
   Index.join
   Index.intersection
   Index.union
   Index.difference
   Index.symmetric_difference

Selecting
~~~~~~~~~
.. autosummary::
   :toctree: api/

   Index.asof
   Index.asof_locs
   Index.get_indexer
   Index.get_indexer_for
   Index.get_indexer_non_unique
   Index.get_level_values
   Index.get_loc
   Index.get_slice_bound
   Index.get_value
   Index.isin
   Index.slice_indexer
   Index.slice_locs

.. _api.numericindex:

Numeric Index
-------------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   RangeIndex
   Int64Index
   UInt64Index
   Float64Index

.. We need this autosummary so that the methods are generated.
.. Separate block, since they aren't classes.

.. autosummary::
   :toctree: api/

   RangeIndex.start
   RangeIndex.stop
   RangeIndex.step
   RangeIndex.from_range

.. _api.categoricalindex:

CategoricalIndex
----------------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   CategoricalIndex

Categorical components
~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   CategoricalIndex.codes
   CategoricalIndex.categories
   CategoricalIndex.ordered
   CategoricalIndex.rename_categories
   CategoricalIndex.reorder_categories
   CategoricalIndex.add_categories
   CategoricalIndex.remove_categories
   CategoricalIndex.remove_unused_categories
   CategoricalIndex.set_categories
   CategoricalIndex.as_ordered
   CategoricalIndex.as_unordered

Modifying and computations
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   CategoricalIndex.map
   CategoricalIndex.equals

.. _api.intervalindex:

IntervalIndex
-------------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   IntervalIndex

IntervalIndex components
~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   IntervalIndex.from_arrays
   IntervalIndex.from_tuples
   IntervalIndex.from_breaks
   IntervalIndex.left
   IntervalIndex.right
   IntervalIndex.mid
   IntervalIndex.closed
   IntervalIndex.length
   IntervalIndex.values
   IntervalIndex.is_empty
   IntervalIndex.is_non_overlapping_monotonic
   IntervalIndex.is_overlapping
   IntervalIndex.get_loc
   IntervalIndex.get_indexer
   IntervalIndex.set_closed
   IntervalIndex.contains
   IntervalIndex.overlaps
   IntervalIndex.to_tuples

.. _api.multiindex:

MultiIndex
----------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   MultiIndex

.. autosummary::
   :toctree: api/

   IndexSlice

MultiIndex constructors
~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   MultiIndex.from_arrays
   MultiIndex.from_tuples
   MultiIndex.from_product
   MultiIndex.from_frame

MultiIndex properties
~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   MultiIndex.names
   MultiIndex.levels
   MultiIndex.codes
   MultiIndex.nlevels
   MultiIndex.levshape

MultiIndex components
~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   MultiIndex.set_levels
   MultiIndex.set_codes
   MultiIndex.to_flat_index
   MultiIndex.to_frame
   MultiIndex.is_lexsorted
   MultiIndex.sortlevel
   MultiIndex.droplevel
   MultiIndex.swaplevel
   MultiIndex.reorder_levels
   MultiIndex.remove_unused_levels

MultiIndex selecting
~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   MultiIndex.get_loc
   MultiIndex.get_locs
   MultiIndex.get_loc_level
   MultiIndex.get_indexer
   MultiIndex.get_level_values

.. _api.datetimeindex:

DatetimeIndex
-------------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   DatetimeIndex

Time/Date components
~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   DatetimeIndex.year
   DatetimeIndex.month
   DatetimeIndex.day
   DatetimeIndex.hour
   DatetimeIndex.minute
   DatetimeIndex.second
   DatetimeIndex.microsecond
   DatetimeIndex.nanosecond
   DatetimeIndex.date
   DatetimeIndex.time
   DatetimeIndex.timetz
   DatetimeIndex.dayofyear
   DatetimeIndex.weekofyear
   DatetimeIndex.week
   DatetimeIndex.dayofweek
   DatetimeIndex.weekday
   DatetimeIndex.quarter
   DatetimeIndex.tz
   DatetimeIndex.freq
   DatetimeIndex.freqstr
   DatetimeIndex.is_month_start
   DatetimeIndex.is_month_end
   DatetimeIndex.is_quarter_start
   DatetimeIndex.is_quarter_end
   DatetimeIndex.is_year_start
   DatetimeIndex.is_year_end
   DatetimeIndex.is_leap_year
   DatetimeIndex.inferred_freq

Selecting
~~~~~~~~~
.. autosummary::
   :toctree: api/

   DatetimeIndex.indexer_at_time
   DatetimeIndex.indexer_between_time


Time-specific operations
~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: api/

   DatetimeIndex.normalize
   DatetimeIndex.strftime
   DatetimeIndex.snap
   DatetimeIndex.tz_convert
   DatetimeIndex.tz_localize
   DatetimeIndex.round
   DatetimeIndex.floor
   DatetimeIndex.ceil
   DatetimeIndex.month_name
   DatetimeIndex.day_name

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: api/

   DatetimeIndex.to_period
   DatetimeIndex.to_perioddelta
   DatetimeIndex.to_pydatetime
   DatetimeIndex.to_series
   DatetimeIndex.to_frame

Methods
~~~~~~~
.. autosummary::
    :toctree: api/

    DatetimeIndex.mean

TimedeltaIndex
--------------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   TimedeltaIndex

Components
~~~~~~~~~~
.. autosummary::
   :toctree: api/

   TimedeltaIndex.days
   TimedeltaIndex.seconds
   TimedeltaIndex.microseconds
   TimedeltaIndex.nanoseconds
   TimedeltaIndex.components
   TimedeltaIndex.inferred_freq

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: api/

   TimedeltaIndex.to_pytimedelta
   TimedeltaIndex.to_series
   TimedeltaIndex.round
   TimedeltaIndex.floor
   TimedeltaIndex.ceil
   TimedeltaIndex.to_frame

Methods
~~~~~~~
.. autosummary::
    :toctree: api/

    TimedeltaIndex.mean

.. currentmodule:: pandas

PeriodIndex
-----------
.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   PeriodIndex

Properties
~~~~~~~~~~
.. autosummary::
    :toctree: api/

    PeriodIndex.day
    PeriodIndex.dayofweek
    PeriodIndex.dayofyear
    PeriodIndex.days_in_month
    PeriodIndex.daysinmonth
    PeriodIndex.end_time
    PeriodIndex.freq
    PeriodIndex.freqstr
    PeriodIndex.hour
    PeriodIndex.is_leap_year
    PeriodIndex.minute
    PeriodIndex.month
    PeriodIndex.quarter
    PeriodIndex.qyear
    PeriodIndex.second
    PeriodIndex.start_time
    PeriodIndex.week
    PeriodIndex.weekday
    PeriodIndex.weekofyear
    PeriodIndex.year

Methods
~~~~~~~
.. autosummary::
    :toctree: api/

    PeriodIndex.asfreq
    PeriodIndex.strftime
    PeriodIndex.to_timestamp
