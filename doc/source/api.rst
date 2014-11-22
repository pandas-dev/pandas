.. currentmodule:: pandas
.. _api:

*************
API Reference
*************

.. _api.functions:

Input/Output
------------

Pickling
~~~~~~~~

.. autosummary::
   :toctree: generated/

   read_pickle

Flat File
~~~~~~~~~

.. autosummary::
   :toctree: generated/

   read_table
   read_csv
   read_fwf

Clipboard
~~~~~~~~~

.. autosummary::
   :toctree: generated/

   read_clipboard

Excel
~~~~~

.. autosummary::
   :toctree: generated/

   read_excel
   ExcelFile.parse

JSON
~~~~

.. autosummary::
   :toctree: generated/

   read_json

HTML
~~~~

.. autosummary::
   :toctree: generated/

   read_html

HDFStore: PyTables (HDF5)
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   read_hdf
   HDFStore.put
   HDFStore.append
   HDFStore.get
   HDFStore.select

SQL
~~~

.. autosummary::
   :toctree: generated/

   read_sql_table
   read_sql_query
   read_sql

Google BigQuery
~~~~~~~~~~~~~~~
.. currentmodule:: pandas.io.gbq

.. autosummary::
   :toctree: generated/

   read_gbq
   to_gbq

.. currentmodule:: pandas


STATA
~~~~~

.. autosummary::
   :toctree: generated/

   read_stata

.. currentmodule:: pandas.io.stata

.. autosummary::
   :toctree: generated/

   StataReader.data
   StataReader.data_label
   StataReader.value_labels
   StataReader.variable_labels
   StataWriter.write_file

.. currentmodule:: pandas

General functions
-----------------

Data manipulations
~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   melt
   pivot
   pivot_table
   crosstab
   cut
   qcut
   merge
   concat
   get_dummies
   factorize

Top-level missing data
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   isnull
   notnull

Top-level dealing with datetimelike
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   to_datetime
   to_timedelta
   date_range
   bdate_range
   period_range
   timedelta_range

Top-level evaluation
~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   eval

Standard moving window functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   rolling_count
   rolling_sum
   rolling_mean
   rolling_median
   rolling_var
   rolling_std
   rolling_min
   rolling_max
   rolling_corr
   rolling_corr_pairwise
   rolling_cov
   rolling_skew
   rolling_kurt
   rolling_apply
   rolling_quantile
   rolling_window

.. _api.functions_expanding:

Standard expanding window functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   expanding_count
   expanding_sum
   expanding_mean
   expanding_median
   expanding_var
   expanding_std
   expanding_min
   expanding_max
   expanding_corr
   expanding_corr_pairwise
   expanding_cov
   expanding_skew
   expanding_kurt
   expanding_apply
   expanding_quantile

Exponentially-weighted moving window functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   ewma
   ewmstd
   ewmvar
   ewmcorr
   ewmcov

.. _api.series:

Series
------

Constructor
~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series

Attributes
~~~~~~~~~~
**Axes**
  * **index**: axis labels

.. autosummary::
   :toctree: generated/

   Series.values
   Series.dtype
   Series.ftype
   Series.shape
   Series.nbytes
   Series.ndim
   Series.size
   Series.strides
   Series.itemsize
   Series.base
   Series.T

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.astype
   Series.copy
   Series.isnull
   Series.notnull

Indexing, iteration
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.get
   Series.at
   Series.iat
   Series.ix
   Series.loc
   Series.iloc
   Series.__iter__
   Series.iteritems

For more information on ``.at``, ``.iat``, ``.ix``, ``.loc``, and
``.iloc``,  see the :ref:`indexing documentation <indexing>`.

Binary operator functions
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.add
   Series.sub
   Series.mul
   Series.div
   Series.truediv
   Series.floordiv
   Series.mod
   Series.pow
   Series.radd
   Series.rsub
   Series.rmul
   Series.rdiv
   Series.rtruediv
   Series.rfloordiv
   Series.rmod
   Series.rpow
   Series.combine
   Series.combine_first
   Series.round
   Series.lt
   Series.gt
   Series.le
   Series.ge
   Series.ne
   Series.eq

Function application, GroupBy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.apply
   Series.map
   Series.groupby

.. _api.series.stats:

Computations / Descriptive Stats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.abs
   Series.all
   Series.any
   Series.autocorr
   Series.between
   Series.clip
   Series.clip_lower
   Series.clip_upper
   Series.corr
   Series.count
   Series.cov
   Series.cummax
   Series.cummin
   Series.cumprod
   Series.cumsum
   Series.describe
   Series.diff
   Series.factorize
   Series.kurt
   Series.mad
   Series.max
   Series.mean
   Series.median
   Series.min
   Series.mode
   Series.pct_change
   Series.prod
   Series.quantile
   Series.rank
   Series.sem
   Series.skew
   Series.std
   Series.sum
   Series.var
   Series.unique
   Series.nunique
   Series.value_counts

Reindexing / Selection / Label manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.align
   Series.drop
   Series.drop_duplicates
   Series.duplicated
   Series.equals
   Series.first
   Series.head
   Series.idxmax
   Series.idxmin
   Series.isin
   Series.last
   Series.reindex
   Series.reindex_like
   Series.rename
   Series.reset_index
   Series.select
   Series.take
   Series.tail
   Series.truncate
   Series.where
   Series.mask

Missing data handling
~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.dropna
   Series.fillna
   Series.interpolate

Reshaping, sorting
~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.argsort
   Series.order
   Series.reorder_levels
   Series.sort
   Series.sort_index
   Series.sortlevel
   Series.swaplevel
   Series.unstack
   Series.searchsorted

Combining / joining / merging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.append
   Series.replace
   Series.update

Time series-related
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.asfreq
   Series.asof
   Series.shift
   Series.first_valid_index
   Series.last_valid_index
   Series.resample
   Series.tz_convert
   Series.tz_localize

Datetimelike Properties
~~~~~~~~~~~~~~~~~~~~~~~

``Series.dt`` can be used to access the values of the series as
datetimelike and return several properties.
Due to implementation details the methods show up here as methods of the
``DatetimeProperties/PeriodProperties/TimedeltaProperties`` classes. These can be accessed like ``Series.dt.<property>``.

.. currentmodule:: pandas.tseries.common

**Datetime Properties**

.. autosummary::
   :toctree: generated/

   DatetimeProperties.date
   DatetimeProperties.time
   DatetimeProperties.year
   DatetimeProperties.month
   DatetimeProperties.day
   DatetimeProperties.hour
   DatetimeProperties.minute
   DatetimeProperties.second
   DatetimeProperties.microsecond
   DatetimeProperties.nanosecond
   DatetimeProperties.second
   DatetimeProperties.weekofyear
   DatetimeProperties.dayofweek
   DatetimeProperties.weekday
   DatetimeProperties.dayofyear
   DatetimeProperties.quarter
   DatetimeProperties.is_month_start
   DatetimeProperties.is_month_end
   DatetimeProperties.is_quarter_start
   DatetimeProperties.is_quarter_end
   DatetimeProperties.is_year_start
   DatetimeProperties.is_year_end

**Datetime Methods**

.. autosummary::
   :toctree: generated/

   DatetimeProperties.to_period
   DatetimeProperties.to_pydatetime
   DatetimeProperties.tz_localize
   DatetimeProperties.tz_convert

**Timedelta Properties**

.. autosummary::
   :toctree: generated/

   TimedeltaProperties.days
   TimedeltaProperties.hours
   TimedeltaProperties.minutes
   TimedeltaProperties.seconds
   TimedeltaProperties.milliseconds
   TimedeltaProperties.microseconds
   TimedeltaProperties.nanoseconds
   TimedeltaProperties.components

**Timedelta Methods**

.. autosummary::
   :toctree: generated/

   TimedeltaProperties.to_pytimedelta

String handling
~~~~~~~~~~~~~~~
``Series.str`` can be used to access the values of the series as
strings and apply several methods to it. Due to implementation
details the methods show up here as methods of the
``StringMethods`` class. These can be acccessed like ``Series.str.<function/property>``.

.. currentmodule:: pandas.core.strings

.. autosummary::
   :toctree: generated/

   StringMethods.cat
   StringMethods.center
   StringMethods.contains
   StringMethods.count
   StringMethods.decode
   StringMethods.encode
   StringMethods.endswith
   StringMethods.extract
   StringMethods.findall
   StringMethods.get
   StringMethods.join
   StringMethods.len
   StringMethods.lower
   StringMethods.lstrip
   StringMethods.match
   StringMethods.pad
   StringMethods.repeat
   StringMethods.replace
   StringMethods.rstrip
   StringMethods.slice
   StringMethods.slice_replace
   StringMethods.split
   StringMethods.startswith
   StringMethods.strip
   StringMethods.title
   StringMethods.upper
   StringMethods.get_dummies

.. _api.categorical:

Categorical
~~~~~~~~~~~

.. currentmodule:: pandas.core.categorical

If the Series is of dtype ``category``, ``Series.cat`` can be used to change the the categorical
data. This accessor is similar to the ``Series.dt`` or ``Series.str`` and has the
following usable methods and properties (all available as ``Series.cat.<method_or_property>``).

.. autosummary::
   :toctree: generated/

   Categorical.categories
   Categorical.ordered
   Categorical.rename_categories
   Categorical.reorder_categories
   Categorical.add_categories
   Categorical.remove_categories
   Categorical.remove_unused_categories
   Categorical.set_categories
   Categorical.codes

To create a Series of dtype ``category``, use ``cat = s.astype("category")``.

The following two ``Categorical`` constructors are considered API but should only be used when
adding ordering information or special categories is need at creation time of the categorical data:

.. autosummary::
   :toctree: generated/

   Categorical
   Categorical.from_codes

``np.asarray(categorical)`` works by implementing the array interface. Be aware, that this converts
the Categorical back to a numpy array, so levels and order information is not preserved!

.. autosummary::
   :toctree: generated/

   Categorical.__array__

Plotting
~~~~~~~~
.. currentmodule:: pandas

.. autosummary::
   :toctree: generated/

   Series.hist
   Series.plot

Serialization / IO / Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.from_csv
   Series.to_pickle
   Series.to_csv
   Series.to_dict
   Series.to_frame
   Series.to_hdf
   Series.to_sql
   Series.to_msgpack
   Series.to_json
   Series.to_sparse
   Series.to_dense
   Series.to_string
   Series.to_clipboard

.. _api.dataframe:

DataFrame
---------

Constructor
~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame

Attributes and underlying data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Axes**

  * **index**: row labels
  * **columns**: column labels

.. autosummary::
   :toctree: generated/

   DataFrame.as_matrix
   DataFrame.dtypes
   DataFrame.ftypes
   DataFrame.get_dtype_counts
   DataFrame.get_ftype_counts
   DataFrame.select_dtypes
   DataFrame.values
   DataFrame.axes
   DataFrame.ndim
   DataFrame.size
   DataFrame.shape

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.astype
   DataFrame.convert_objects
   DataFrame.copy
   DataFrame.isnull
   DataFrame.notnull

Indexing, iteration
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.head
   DataFrame.at
   DataFrame.iat
   DataFrame.ix
   DataFrame.loc
   DataFrame.iloc
   DataFrame.insert
   DataFrame.__iter__
   DataFrame.iteritems
   DataFrame.iterrows
   DataFrame.itertuples
   DataFrame.lookup
   DataFrame.pop
   DataFrame.tail
   DataFrame.xs
   DataFrame.isin
   DataFrame.where
   DataFrame.mask
   DataFrame.query

For more information on ``.at``, ``.iat``, ``.ix``, ``.loc``, and
``.iloc``,  see the :ref:`indexing documentation <indexing>`.


Binary operator functions
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.add
   DataFrame.sub
   DataFrame.mul
   DataFrame.div
   DataFrame.truediv
   DataFrame.floordiv
   DataFrame.mod
   DataFrame.pow
   DataFrame.radd
   DataFrame.rsub
   DataFrame.rmul
   DataFrame.rdiv
   DataFrame.rtruediv
   DataFrame.rfloordiv
   DataFrame.rmod
   DataFrame.rpow
   DataFrame.lt
   DataFrame.gt
   DataFrame.le
   DataFrame.ge
   DataFrame.ne
   DataFrame.eq
   DataFrame.combine
   DataFrame.combineAdd
   DataFrame.combine_first
   DataFrame.combineMult

Function application, GroupBy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.apply
   DataFrame.applymap
   DataFrame.groupby

.. _api.dataframe.stats:

Computations / Descriptive Stats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.abs
   DataFrame.all
   DataFrame.any
   DataFrame.clip
   DataFrame.clip_lower
   DataFrame.clip_upper
   DataFrame.corr
   DataFrame.corrwith
   DataFrame.count
   DataFrame.cov
   DataFrame.cummax
   DataFrame.cummin
   DataFrame.cumprod
   DataFrame.cumsum
   DataFrame.describe
   DataFrame.diff
   DataFrame.eval
   DataFrame.kurt
   DataFrame.mad
   DataFrame.max
   DataFrame.mean
   DataFrame.median
   DataFrame.min
   DataFrame.mode
   DataFrame.pct_change
   DataFrame.prod
   DataFrame.quantile
   DataFrame.rank
   DataFrame.sem
   DataFrame.skew
   DataFrame.sum
   DataFrame.std
   DataFrame.var

Reindexing / Selection / Label manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.add_prefix
   DataFrame.add_suffix
   DataFrame.align
   DataFrame.drop
   DataFrame.drop_duplicates
   DataFrame.duplicated
   DataFrame.equals
   DataFrame.filter
   DataFrame.first
   DataFrame.head
   DataFrame.idxmax
   DataFrame.idxmin
   DataFrame.last
   DataFrame.reindex
   DataFrame.reindex_axis
   DataFrame.reindex_like
   DataFrame.rename
   DataFrame.reset_index
   DataFrame.select
   DataFrame.set_index
   DataFrame.tail
   DataFrame.take
   DataFrame.truncate

.. _api.dataframe.missing:

Missing data handling
~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.dropna
   DataFrame.fillna
   DataFrame.replace

Reshaping, sorting, transposing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.pivot
   DataFrame.reorder_levels
   DataFrame.sort
   DataFrame.sort_index
   DataFrame.sortlevel
   DataFrame.swaplevel
   DataFrame.stack
   DataFrame.unstack
   DataFrame.T
   DataFrame.to_panel
   DataFrame.transpose

Combining / joining / merging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.append
   DataFrame.join
   DataFrame.merge
   DataFrame.update

Time series-related
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.asfreq
   DataFrame.shift
   DataFrame.first_valid_index
   DataFrame.last_valid_index
   DataFrame.resample
   DataFrame.to_period
   DataFrame.to_timestamp
   DataFrame.tz_convert
   DataFrame.tz_localize

Plotting
~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.boxplot
   DataFrame.hist
   DataFrame.plot

Serialization / IO / Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.from_csv
   DataFrame.from_dict
   DataFrame.from_items
   DataFrame.from_records
   DataFrame.info
   DataFrame.to_pickle
   DataFrame.to_csv
   DataFrame.to_hdf
   DataFrame.to_sql
   DataFrame.to_dict
   DataFrame.to_excel
   DataFrame.to_json
   DataFrame.to_html
   DataFrame.to_latex
   DataFrame.to_stata
   DataFrame.to_msgpack
   DataFrame.to_gbq
   DataFrame.to_records
   DataFrame.to_sparse
   DataFrame.to_dense
   DataFrame.to_string
   DataFrame.to_clipboard

.. _api.panel:

Panel
------

Constructor
~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel

Attributes and underlying data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Axes**

  * **items**: axis 0; each item corresponds to a DataFrame contained inside
  * **major_axis**: axis 1; the index (rows) of each of the DataFrames
  * **minor_axis**: axis 2; the columns of each of the DataFrames

.. autosummary::
   :toctree: generated/

   Panel.values
   Panel.axes
   Panel.ndim
   Panel.size
   Panel.shape
   Panel.dtypes
   Panel.ftypes
   Panel.get_dtype_counts
   Panel.get_ftype_counts

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.astype
   Panel.copy
   Panel.isnull
   Panel.notnull

Getting and setting
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.get_value
   Panel.set_value

Indexing, iteration, slicing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.at
   Panel.iat
   Panel.ix
   Panel.loc
   Panel.iloc
   Panel.__iter__
   Panel.iteritems
   Panel.pop
   Panel.xs
   Panel.major_xs
   Panel.minor_xs

For more information on ``.at``, ``.iat``, ``.ix``, ``.loc``, and
``.iloc``,  see the :ref:`indexing documentation <indexing>`.

Binary operator functions
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.add
   Panel.sub
   Panel.mul
   Panel.div
   Panel.truediv
   Panel.floordiv
   Panel.mod
   Panel.pow
   Panel.radd
   Panel.rsub
   Panel.rmul
   Panel.rdiv
   Panel.rtruediv
   Panel.rfloordiv
   Panel.rmod
   Panel.rpow
   Panel.lt
   Panel.gt
   Panel.le
   Panel.ge
   Panel.ne
   Panel.eq

Function application, GroupBy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.apply
   Panel.groupby

.. _api.panel.stats:

Computations / Descriptive Stats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.abs
   Panel.clip
   Panel.clip_lower
   Panel.clip_upper
   Panel.count
   Panel.cummax
   Panel.cummin
   Panel.cumprod
   Panel.cumsum
   Panel.max
   Panel.mean
   Panel.median
   Panel.min
   Panel.pct_change
   Panel.prod
   Panel.sem
   Panel.skew
   Panel.sum
   Panel.std
   Panel.var

Reindexing / Selection / Label manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.add_prefix
   Panel.add_suffix
   Panel.drop
   Panel.equals
   Panel.filter
   Panel.first
   Panel.last
   Panel.reindex
   Panel.reindex_axis
   Panel.reindex_like
   Panel.rename
   Panel.select
   Panel.take
   Panel.truncate

Missing data handling
~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.dropna
   Panel.fillna

Reshaping, sorting, transposing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.sort_index
   Panel.swaplevel
   Panel.transpose
   Panel.swapaxes
   Panel.conform

Combining / joining / merging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.join
   Panel.update

Time series-related
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.asfreq
   Panel.shift
   Panel.resample
   Panel.tz_convert
   Panel.tz_localize

Serialization / IO / Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.from_dict
   Panel.to_pickle
   Panel.to_excel
   Panel.to_hdf
   Panel.to_json
   Panel.to_sparse
   Panel.to_frame
   Panel.to_clipboard

.. _api.panel4d:

Panel4D
-------

Constructor
~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel4D

Attributes and underlying data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Axes**

  * **labels**: axis 1; each label corresponds to a Panel contained inside
  * **items**: axis 2; each item corresponds to a DataFrame contained inside
  * **major_axis**: axis 3; the index (rows) of each of the DataFrames
  * **minor_axis**: axis 4; the columns of each of the DataFrames

.. autosummary::
   :toctree: generated/

   Panel4D.values
   Panel4D.axes
   Panel4D.ndim
   Panel4D.size
   Panel4D.shape
   Panel4D.dtypes
   Panel4D.ftypes
   Panel4D.get_dtype_counts
   Panel4D.get_ftype_counts

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel4D.astype
   Panel4D.copy
   Panel4D.isnull
   Panel4D.notnull

.. _api.index:

Index
-----

**Many of these methods or variants thereof are available on the objects
that contain an index (Series/Dataframe) and those should most likely be
used before calling these methods directly.**

.. autosummary::
   :toctree: generated/

   Index

Attributes
~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Index.values
   Index.is_monotonic
   Index.is_monotonic_increasing
   Index.is_monotonic_decreasing
   Index.is_unique
   Index.dtype
   Index.inferred_type
   Index.is_all_dates
   Index.shape
   Index.nbytes
   Index.ndim
   Index.size
   Index.strides
   Index.itemsize
   Index.base
   Index.T

Modifying and Computations
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.all
   Index.any
   Index.argmin
   Index.argmax
   Index.copy
   Index.delete
   Index.diff
   Index.sym_diff
   Index.drop
   Index.drop_duplicates
   Index.duplicated
   Index.equals
   Index.factorize
   Index.identical
   Index.insert
   Index.min
   Index.max
   Index.order
   Index.reindex
   Index.repeat
   Index.take
   Index.putmask
   Index.set_names
   Index.unique
   Index.nunique
   Index.value_counts

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.astype
   Index.tolist
   Index.to_datetime
   Index.to_series

Sorting
~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.argsort
   Index.order
   Index.sort

Time-specific operations
~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.shift

Combining / joining / merging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.append
   Index.intersection
   Index.join
   Index.union

Selecting
~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.get_indexer
   Index.get_indexer_non_unique
   Index.get_level_values
   Index.get_loc
   Index.get_value
   Index.isin
   Index.slice_indexer
   Index.slice_locs

.. _api.datetimeindex:

DatetimeIndex
-------------

.. autosummary::
   :toctree: generated/

   DatetimeIndex

Time/Date Components
~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

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

Selecting
~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DatetimeIndex.indexer_at_time
   DatetimeIndex.indexer_between_time


Time-specific operations
~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DatetimeIndex.normalize
   DatetimeIndex.snap
   DatetimeIndex.tz_convert
   DatetimeIndex.tz_localize


Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DatetimeIndex.to_datetime
   DatetimeIndex.to_period
   DatetimeIndex.to_pydatetime
   DatetimeIndex.to_series

TimedeltaIndex
--------------

.. autosummary::
   :toctree: generated/

   TimedeltaIndex

Components
~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   TimedeltaIndex.days
   TimedeltaIndex.hours
   TimedeltaIndex.minutes
   TimedeltaIndex.seconds
   TimedeltaIndex.milliseconds
   TimedeltaIndex.microseconds
   TimedeltaIndex.nanoseconds
   TimedeltaIndex.components

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   TimedeltaIndex.to_pytimedelta
   TimedeltaIndex.to_series

GroupBy
-------
.. currentmodule:: pandas.core.groupby

GroupBy objects are returned by groupby calls: :func:`pandas.DataFrame.groupby`, :func:`pandas.Series.groupby`, etc.

Indexing, iteration
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   GroupBy.__iter__
   GroupBy.groups
   GroupBy.indices
   GroupBy.get_group

.. currentmodule:: pandas

.. autosummary::
   :toctree: generated/

   Grouper

.. currentmodule:: pandas.core.groupby

Function application
~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   GroupBy.apply
   GroupBy.aggregate
   GroupBy.transform

Computations / Descriptive Stats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   GroupBy.count
   GroupBy.cumcount
   GroupBy.first
   GroupBy.head
   GroupBy.last
   GroupBy.max
   GroupBy.mean
   GroupBy.median
   GroupBy.min
   GroupBy.nth
   GroupBy.ohlc
   GroupBy.prod
   GroupBy.size
   GroupBy.sem
   GroupBy.std
   GroupBy.sum
   GroupBy.var
   GroupBy.tail

The following methods are available in both ``SeriesGroupBy`` and
``DataFrameGroupBy`` objects, but may differ slightly, usually in that
the ``DataFrameGroupBy`` version usually permits the specification of an
axis argument, and often an argument indicating whether to restrict
application to columns of a specific data type.

.. autosummary::
   :toctree: generated/

   DataFrameGroupBy.bfill
   DataFrameGroupBy.cummax
   DataFrameGroupBy.cummin
   DataFrameGroupBy.cumprod
   DataFrameGroupBy.cumsum
   DataFrameGroupBy.describe
   DataFrameGroupBy.all
   DataFrameGroupBy.any
   DataFrameGroupBy.corr
   DataFrameGroupBy.cov
   DataFrameGroupBy.diff
   DataFrameGroupBy.ffill
   DataFrameGroupBy.fillna
   DataFrameGroupBy.hist
   DataFrameGroupBy.idxmax
   DataFrameGroupBy.idxmin
   DataFrameGroupBy.irow
   DataFrameGroupBy.mad
   DataFrameGroupBy.pct_change
   DataFrameGroupBy.plot
   DataFrameGroupBy.quantile
   DataFrameGroupBy.rank
   DataFrameGroupBy.resample
   DataFrameGroupBy.shift
   DataFrameGroupBy.skew
   DataFrameGroupBy.take
   DataFrameGroupBy.tshift

The following methods are available only for ``SeriesGroupBy`` objects.

.. autosummary::
   :toctree: generated/

   SeriesGroupBy.nlargest
   SeriesGroupBy.nsmallest
   SeriesGroupBy.nunique
   SeriesGroupBy.unique
   SeriesGroupBy.value_counts

The following methods are available only for ``DataFrameGroupBy`` objects.

.. autosummary::
   :toctree: generated/

   DataFrameGroupBy.corrwith
   DataFrameGroupBy.boxplot

.. currentmodule:: pandas

General utility functions
-------------------------

Working with options
~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   describe_option
   reset_option
   get_option
   set_option
   option_context


..
    HACK - see github issue #4539. To ensure old links remain valid, include
    here the autosummaries with previous currentmodules as a comment and add
    them to a hidden toctree (to avoid warnings):

.. toctree::
   :hidden:

   generated/pandas.core.common.isnull
   generated/pandas.core.common.notnull
   generated/pandas.core.reshape.get_dummies
   generated/pandas.io.clipboard.read_clipboard
   generated/pandas.io.excel.ExcelFile.parse
   generated/pandas.io.excel.read_excel
   generated/pandas.io.html.read_html
   generated/pandas.io.json.read_json
   generated/pandas.io.parsers.read_csv
   generated/pandas.io.parsers.read_fwf
   generated/pandas.io.parsers.read_table
   generated/pandas.io.pickle.read_pickle
   generated/pandas.io.pytables.HDFStore.append
   generated/pandas.io.pytables.HDFStore.get
   generated/pandas.io.pytables.HDFStore.put
   generated/pandas.io.pytables.HDFStore.select
   generated/pandas.io.pytables.read_hdf
   generated/pandas.io.sql.read_sql
   generated/pandas.io.sql.read_frame
   generated/pandas.io.sql.write_frame
   generated/pandas.io.stata.read_stata
   generated/pandas.stats.moments.ewma
   generated/pandas.stats.moments.ewmcorr
   generated/pandas.stats.moments.ewmcov
   generated/pandas.stats.moments.ewmstd
   generated/pandas.stats.moments.ewmvar
   generated/pandas.stats.moments.expanding_apply
   generated/pandas.stats.moments.expanding_corr
   generated/pandas.stats.moments.expanding_count
   generated/pandas.stats.moments.expanding_cov
   generated/pandas.stats.moments.expanding_kurt
   generated/pandas.stats.moments.expanding_mean
   generated/pandas.stats.moments.expanding_median
   generated/pandas.stats.moments.expanding_quantile
   generated/pandas.stats.moments.expanding_skew
   generated/pandas.stats.moments.expanding_std
   generated/pandas.stats.moments.expanding_sum
   generated/pandas.stats.moments.expanding_var
   generated/pandas.stats.moments.rolling_apply
   generated/pandas.stats.moments.rolling_corr
   generated/pandas.stats.moments.rolling_count
   generated/pandas.stats.moments.rolling_cov
   generated/pandas.stats.moments.rolling_kurt
   generated/pandas.stats.moments.rolling_mean
   generated/pandas.stats.moments.rolling_median
   generated/pandas.stats.moments.rolling_quantile
   generated/pandas.stats.moments.rolling_skew
   generated/pandas.stats.moments.rolling_std
   generated/pandas.stats.moments.rolling_sum
   generated/pandas.stats.moments.rolling_var
   generated/pandas.tools.merge.concat
   generated/pandas.tools.merge.merge
   generated/pandas.tools.pivot.pivot_table
   generated/pandas.tseries.tools.to_datetime

..
    .. currentmodule:: pandas.io.pickle

    .. autosummary::
       :toctree: generated/

       read_pickle

    .. currentmodule:: pandas.io.parsers

    .. autosummary::
       :toctree: generated/

       read_table
       read_csv
       read_fwf

    .. currentmodule:: pandas.io.clipboard

    .. autosummary::
       :toctree: generated/

       read_clipboard

    .. currentmodule:: pandas.io.excel

    .. autosummary::
       :toctree: generated/

       read_excel
       ExcelFile.parse

    .. currentmodule:: pandas.io.json

    .. autosummary::
       :toctree: generated/

       read_json

    .. currentmodule:: pandas.io.html

    .. autosummary::
       :toctree: generated/

       read_html

    .. currentmodule:: pandas.io.pytables

    .. autosummary::
       :toctree: generated/

       read_hdf
       HDFStore.put
       HDFStore.append
       HDFStore.get
       HDFStore.select

    .. currentmodule:: pandas.io.sql

    .. autosummary::
       :toctree: generated/

       read_sql
       read_frame
       write_frame

    .. currentmodule:: pandas.io.stata

    .. autosummary::
       :toctree: generated/

       read_stata
       StataReader.data
       StataReader.data_label
       StataReader.value_labels
       StataReader.variable_labels
       StataWriter.write_file

    .. currentmodule:: pandas.tools.pivot

    .. autosummary::
       :toctree: generated/

       pivot_table

    .. currentmodule:: pandas.tools.merge

    .. autosummary::
       :toctree: generated/

       merge
       concat

    .. currentmodule:: pandas.core.reshape

    .. autosummary::
       :toctree: generated/

       get_dummies

    .. currentmodule:: pandas.core.common

    .. autosummary::
       :toctree: generated/

       isnull
       notnull

    .. currentmodule:: pandas.tseries.tools

    .. autosummary::
       :toctree: generated/

       to_datetime


    .. currentmodule:: pandas.stats.moments

    .. autosummary::
       :toctree: generated/

       rolling_count
       rolling_sum
       rolling_mean
       rolling_median
       rolling_var
       rolling_std
       rolling_corr
       rolling_cov
       rolling_skew
       rolling_kurt
       rolling_apply
       rolling_quantile


    .. currentmodule:: pandas.stats.moments

    .. autosummary::
       :toctree: generated/

       expanding_count
       expanding_sum
       expanding_mean
       expanding_median
       expanding_var
       expanding_std
       expanding_corr
       expanding_cov
       expanding_skew
       expanding_kurt
       expanding_apply
       expanding_quantile


    .. autosummary::
       :toctree: generated/

       ewma
       ewmstd
       ewmvar
       ewmcorr
       ewmcov
