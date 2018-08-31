.. currentmodule:: pandas
.. _api:

*************
API Reference
*************

This page gives an overview of all public pandas objects, functions and
methods. All classes and functions exposed in ``pandas.*`` namespace are public.

Some subpackages are public which include ``pandas.errors``,
``pandas.plotting``, and ``pandas.testing``. Public functions in
``pandas.io`` and ``pandas.tseries`` submodules are mentioned in
the documentation. ``pandas.api.types`` subpackage holds some
public functions related to data types in pandas.


.. warning::

    The ``pandas.core``, ``pandas.compat``, and ``pandas.util`` top-level modules are PRIVATE. Stable functionality in such modules is not guaranteed.


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
   read_msgpack

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

.. currentmodule:: pandas.io.json

.. autosummary::
   :toctree: generated/

   json_normalize
   build_table_schema

.. currentmodule:: pandas

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
   HDFStore.info
   HDFStore.keys
   HDFStore.walk

Feather
~~~~~~~

.. autosummary::
   :toctree: generated/

   read_feather

Parquet
~~~~~~~

.. autosummary::
   :toctree: generated/

   read_parquet

SAS
~~~

.. autosummary::
   :toctree: generated/

   read_sas

SQL
~~~

.. autosummary::
   :toctree: generated/

   read_sql_table
   read_sql_query
   read_sql

Google BigQuery
~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   read_gbq


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
   merge_ordered
   merge_asof
   concat
   get_dummies
   factorize
   unique
   wide_to_long

Top-level missing data
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   isna
   isnull
   notna
   notnull

Top-level conversions
~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   to_numeric

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
   infer_freq

Top-level dealing with intervals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   interval_range

Top-level evaluation
~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   eval

Testing
~~~~~~~

.. autosummary::
   :toctree: generated/

   test

.. _api.series:

Series
------

Constructor
~~~~~~~~~~~

.. currentmodule:: pandas

.. autosummary::
   :toctree: generated/

   Series

Attributes
~~~~~~~~~~
**Axes**

.. autosummary::
   :toctree: generated/

   Series.index

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
   Series.memory_usage
   Series.hasnans
   Series.flags
   Series.empty
   Series.dtypes
   Series.ftypes
   Series.data
   Series.is_copy
   Series.name
   Series.put

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.astype
   Series.infer_objects
   Series.convert_objects
   Series.copy
   Series.bool
   Series.to_period
   Series.to_timestamp
   Series.tolist
   Series.get_values


Indexing, iteration
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.get
   Series.at
   Series.iat
   Series.loc
   Series.iloc
   Series.__iter__
   Series.iteritems
   Series.items
   Series.keys
   Series.pop
   Series.item
   Series.xs

For more information on ``.at``, ``.iat``, ``.loc``, and
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
   Series.product
   Series.dot

Function application, GroupBy & Window
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.apply
   Series.agg
   Series.aggregate
   Series.transform
   Series.map
   Series.groupby
   Series.rolling
   Series.expanding
   Series.ewm
   Series.pipe

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
   Series.nlargest
   Series.nsmallest
   Series.pct_change
   Series.prod
   Series.quantile
   Series.rank
   Series.sem
   Series.skew
   Series.std
   Series.sum
   Series.var
   Series.kurtosis
   Series.unique
   Series.nunique
   Series.is_unique
   Series.is_monotonic
   Series.is_monotonic_increasing
   Series.is_monotonic_decreasing
   Series.value_counts
   Series.compound
   Series.nonzero


Reindexing / Selection / Label manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.align
   Series.drop
   Series.droplevel
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
   Series.rename_axis
   Series.reset_index
   Series.sample
   Series.select
   Series.set_axis
   Series.take
   Series.tail
   Series.truncate
   Series.where
   Series.mask
   Series.add_prefix
   Series.add_suffix
   Series.filter

Missing data handling
~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.isna
   Series.notna
   Series.dropna
   Series.fillna
   Series.interpolate

Reshaping, sorting
~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.argsort
   Series.argmin
   Series.argmax
   Series.reorder_levels
   Series.sort_values
   Series.sort_index
   Series.swaplevel
   Series.unstack
   Series.searchsorted
   Series.ravel
   Series.repeat
   Series.squeeze
   Series.view
   Series.sortlevel


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
   Series.at_time
   Series.between_time
   Series.tshift
   Series.slice_shift

Datetimelike Properties
~~~~~~~~~~~~~~~~~~~~~~~

``Series.dt`` can be used to access the values of the series as
datetimelike and return several properties.
These can be accessed like ``Series.dt.<property>``.

**Datetime Properties**

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_attribute.rst

   Series.dt.date
   Series.dt.time
   Series.dt.timetz
   Series.dt.year
   Series.dt.month
   Series.dt.day
   Series.dt.hour
   Series.dt.minute
   Series.dt.second
   Series.dt.microsecond
   Series.dt.nanosecond
   Series.dt.week
   Series.dt.weekofyear
   Series.dt.dayofweek
   Series.dt.weekday
   Series.dt.dayofyear
   Series.dt.quarter
   Series.dt.is_month_start
   Series.dt.is_month_end
   Series.dt.is_quarter_start
   Series.dt.is_quarter_end
   Series.dt.is_year_start
   Series.dt.is_year_end
   Series.dt.is_leap_year
   Series.dt.daysinmonth
   Series.dt.days_in_month
   Series.dt.tz
   Series.dt.freq

**Datetime Methods**

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   Series.dt.to_period
   Series.dt.to_pydatetime
   Series.dt.tz_localize
   Series.dt.tz_convert
   Series.dt.is_dst
   Series.dt.normalize
   Series.dt.strftime
   Series.dt.round
   Series.dt.floor
   Series.dt.ceil
   Series.dt.month_name
   Series.dt.day_name

**Timedelta Properties**

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_attribute.rst

   Series.dt.days
   Series.dt.seconds
   Series.dt.microseconds
   Series.dt.nanoseconds
   Series.dt.components

**Timedelta Methods**

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   Series.dt.to_pytimedelta
   Series.dt.total_seconds

String handling
~~~~~~~~~~~~~~~
``Series.str`` can be used to access the values of the series as
strings and apply several methods to it. These can be accessed like
``Series.str.<function/property>``.

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   Series.str.capitalize
   Series.str.cat
   Series.str.center
   Series.str.contains
   Series.str.count
   Series.str.decode
   Series.str.encode
   Series.str.endswith
   Series.str.extract
   Series.str.extractall
   Series.str.find
   Series.str.findall
   Series.str.get
   Series.str.index
   Series.str.join
   Series.str.len
   Series.str.ljust
   Series.str.lower
   Series.str.lstrip
   Series.str.match
   Series.str.normalize
   Series.str.pad
   Series.str.partition
   Series.str.repeat
   Series.str.replace
   Series.str.rfind
   Series.str.rindex
   Series.str.rjust
   Series.str.rpartition
   Series.str.rstrip
   Series.str.slice
   Series.str.slice_replace
   Series.str.split
   Series.str.rsplit
   Series.str.startswith
   Series.str.strip
   Series.str.swapcase
   Series.str.title
   Series.str.translate
   Series.str.upper
   Series.str.wrap
   Series.str.zfill
   Series.str.isalnum
   Series.str.isalpha
   Series.str.isdigit
   Series.str.isspace
   Series.str.islower
   Series.str.isupper
   Series.str.istitle
   Series.str.isnumeric
   Series.str.isdecimal
   Series.str.get_dummies

..
    The following is needed to ensure the generated pages are created with the
    correct template (otherwise they would be created in the Series/Index class page)

..
    .. autosummary::
       :toctree: generated/
       :template: autosummary/accessor.rst

       Series.str
       Series.cat
       Series.dt
       Index.str

.. _api.categorical:

Categorical
~~~~~~~~~~~

Pandas defines a custom data type for representing data that can take only a
limited, fixed set of values. The dtype of a ``Categorical`` can be described by
a :class:`pandas.api.types.CategoricalDtype`.

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_autosummary.rst

   api.types.CategoricalDtype

.. autosummary::
   :toctree: generated/

   api.types.CategoricalDtype.categories
   api.types.CategoricalDtype.ordered

Categorical data can be stored in a :class:`pandas.Categorical`

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_autosummary.rst

   Categorical


The alternative :meth:`Categorical.from_codes` constructor can be used when you
have the categories and integer codes already:

.. autosummary::
   :toctree: generated/

   Categorical.from_codes

The dtype information is available on the ``Categorical``

.. autosummary::
   :toctree: generated/

   Categorical.dtype
   Categorical.categories
   Categorical.ordered
   Categorical.codes

``np.asarray(categorical)`` works by implementing the array interface. Be aware, that this converts
the Categorical back to a NumPy array, so categories and order information is not preserved!

.. autosummary::
   :toctree: generated/

   Categorical.__array__

A ``Categorical`` can be stored in a ``Series`` or ``DataFrame``.
To create a Series of dtype ``category``, use ``cat = s.astype(dtype)`` or
``Series(..., dtype=dtype)`` where ``dtype`` is either

* the string ``'category'``
* an instance of :class:`~pandas.api.types.CategoricalDtype`.

If the Series is of dtype ``CategoricalDtype``, ``Series.cat`` can be used to change the categorical
data. This accessor is similar to the ``Series.dt`` or ``Series.str`` and has the
following usable methods and properties:

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_attribute.rst

   Series.cat.categories
   Series.cat.ordered
   Series.cat.codes

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   Series.cat.rename_categories
   Series.cat.reorder_categories
   Series.cat.add_categories
   Series.cat.remove_categories
   Series.cat.remove_unused_categories
   Series.cat.set_categories
   Series.cat.as_ordered
   Series.cat.as_unordered

Plotting
~~~~~~~~

``Series.plot`` is both a callable method and a namespace attribute for
specific plotting methods of the form ``Series.plot.<kind>``.

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_callable.rst

   Series.plot

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   Series.plot.area
   Series.plot.bar
   Series.plot.barh
   Series.plot.box
   Series.plot.density
   Series.plot.hist
   Series.plot.kde
   Series.plot.line
   Series.plot.pie

.. autosummary::
   :toctree: generated/

   Series.hist

Serialization / IO / Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.to_pickle
   Series.to_csv
   Series.to_dict
   Series.to_excel
   Series.to_frame
   Series.to_xarray
   Series.to_hdf
   Series.to_sql
   Series.to_msgpack
   Series.to_json
   Series.to_sparse
   Series.to_dense
   Series.to_string
   Series.to_clipboard
   Series.to_latex

Sparse
~~~~~~
.. autosummary::
   :toctree: generated/

   SparseSeries.to_coo
   SparseSeries.from_coo

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

.. autosummary::
   :toctree: generated/

   DataFrame.index
   DataFrame.columns

.. autosummary::
   :toctree: generated/

   DataFrame.dtypes
   DataFrame.ftypes
   DataFrame.get_dtype_counts
   DataFrame.get_ftype_counts
   DataFrame.select_dtypes
   DataFrame.values
   DataFrame.get_values
   DataFrame.axes
   DataFrame.ndim
   DataFrame.size
   DataFrame.shape
   DataFrame.memory_usage
   DataFrame.empty
   DataFrame.is_copy

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.astype
   DataFrame.convert_objects
   DataFrame.infer_objects
   DataFrame.copy
   DataFrame.isna
   DataFrame.notna
   DataFrame.bool

Indexing, iteration
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.head
   DataFrame.at
   DataFrame.iat
   DataFrame.loc
   DataFrame.iloc
   DataFrame.insert
   DataFrame.insert
   DataFrame.__iter__
   DataFrame.items
   DataFrame.keys
   DataFrame.iteritems
   DataFrame.iterrows
   DataFrame.itertuples
   DataFrame.lookup
   DataFrame.pop
   DataFrame.tail
   DataFrame.xs
   DataFrame.get
   DataFrame.isin
   DataFrame.where
   DataFrame.mask
   DataFrame.query

For more information on ``.at``, ``.iat``, ``.loc``, and
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
   DataFrame.dot
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
   DataFrame.combine_first

Function application, GroupBy & Window
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.apply
   DataFrame.applymap
   DataFrame.pipe
   DataFrame.agg
   DataFrame.aggregate
   DataFrame.transform
   DataFrame.groupby
   DataFrame.rolling
   DataFrame.expanding
   DataFrame.ewm

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
   DataFrame.compound
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
   DataFrame.kurtosis
   DataFrame.mad
   DataFrame.max
   DataFrame.mean
   DataFrame.median
   DataFrame.min
   DataFrame.mode
   DataFrame.pct_change
   DataFrame.prod
   DataFrame.product
   DataFrame.quantile
   DataFrame.rank
   DataFrame.round
   DataFrame.sem
   DataFrame.skew
   DataFrame.sum
   DataFrame.std
   DataFrame.var
   DataFrame.nunique

Reindexing / Selection / Label manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.add_prefix
   DataFrame.add_suffix
   DataFrame.align
   DataFrame.at_time
   DataFrame.between_time
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
   DataFrame.rename_axis
   DataFrame.reset_index
   DataFrame.sample
   DataFrame.select
   DataFrame.set_axis
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
   DataFrame.interpolate

Reshaping, sorting, transposing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.droplevel
   DataFrame.pivot
   DataFrame.pivot_table
   DataFrame.reorder_levels
   DataFrame.sort_values
   DataFrame.sort_index
   DataFrame.nlargest
   DataFrame.nsmallest
   DataFrame.swaplevel
   DataFrame.stack
   DataFrame.unstack
   DataFrame.swapaxes
   DataFrame.melt
   DataFrame.squeeze
   DataFrame.to_panel
   DataFrame.to_xarray
   DataFrame.T
   DataFrame.transpose

Combining / joining / merging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.append
   DataFrame.assign
   DataFrame.join
   DataFrame.merge
   DataFrame.update

Time series-related
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.asfreq
   DataFrame.asof
   DataFrame.shift
   DataFrame.slice_shift
   DataFrame.tshift
   DataFrame.first_valid_index
   DataFrame.last_valid_index
   DataFrame.resample
   DataFrame.to_period
   DataFrame.to_timestamp
   DataFrame.tz_convert
   DataFrame.tz_localize

.. _api.dataframe.plotting:

Plotting
~~~~~~~~

``DataFrame.plot`` is both a callable method and a namespace attribute for
specific plotting methods of the form ``DataFrame.plot.<kind>``.

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_callable.rst

   DataFrame.plot

.. autosummary::
   :toctree: generated/
   :template: autosummary/accessor_method.rst

   DataFrame.plot.area
   DataFrame.plot.bar
   DataFrame.plot.barh
   DataFrame.plot.box
   DataFrame.plot.density
   DataFrame.plot.hexbin
   DataFrame.plot.hist
   DataFrame.plot.kde
   DataFrame.plot.line
   DataFrame.plot.pie
   DataFrame.plot.scatter

.. autosummary::
   :toctree: generated/

   DataFrame.boxplot
   DataFrame.hist

Serialization / IO / Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.from_csv
   DataFrame.from_dict
   DataFrame.from_items
   DataFrame.from_records
   DataFrame.info
   DataFrame.to_parquet
   DataFrame.to_pickle
   DataFrame.to_csv
   DataFrame.to_hdf
   DataFrame.to_sql
   DataFrame.to_dict
   DataFrame.to_excel
   DataFrame.to_json
   DataFrame.to_html
   DataFrame.to_feather
   DataFrame.to_latex
   DataFrame.to_stata
   DataFrame.to_msgpack
   DataFrame.to_gbq
   DataFrame.to_records
   DataFrame.to_sparse
   DataFrame.to_dense
   DataFrame.to_string
   DataFrame.to_clipboard
   DataFrame.style

Sparse
~~~~~~
.. autosummary::
   :toctree: generated/

   SparseDataFrame.to_coo

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
   Panel.isna
   Panel.notna

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
   Panel.loc
   Panel.iloc
   Panel.__iter__
   Panel.iteritems
   Panel.pop
   Panel.xs
   Panel.major_xs
   Panel.minor_xs

For more information on ``.at``, ``.iat``, ``.loc``, and
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
   Panel.sample
   Panel.select
   Panel.take
   Panel.truncate


Missing data handling
~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.dropna

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
   Panel.to_sparse
   Panel.to_frame
   Panel.to_clipboard

.. _api.index:

Index
-----

**Many of these methods or variants thereof are available on the objects
that contain an index (Series/DataFrame) and those should most likely be
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
   Index.has_duplicates
   Index.hasnans
   Index.dtype
   Index.dtype_str
   Index.inferred_type
   Index.is_all_dates
   Index.shape
   Index.name
   Index.names
   Index.nbytes
   Index.ndim
   Index.size
   Index.empty
   Index.strides
   Index.itemsize
   Index.base
   Index.T
   Index.memory_usage

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
   :toctree: generated/

   Index.set_names
   Index.is_lexsorted_for_tuple
   Index.droplevel

Missing Values
~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.fillna
   Index.dropna
   Index.isna
   Index.notna

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.astype
   Index.item
   Index.map
   Index.ravel
   Index.tolist
   Index.to_native_types
   Index.to_series
   Index.to_frame
   Index.view

Sorting
~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.argsort
   Index.searchsorted
   Index.sort_values

Time-specific operations
~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.shift

Combining / joining / set operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.append
   Index.join
   Index.intersection
   Index.union
   Index.difference
   Index.symmetric_difference

Selecting
~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Index.asof
   Index.asof_locs
   Index.contains
   Index.get_duplicates
   Index.get_indexer
   Index.get_indexer_for
   Index.get_indexer_non_unique
   Index.get_level_values
   Index.get_loc
   Index.get_slice_bound
   Index.get_value
   Index.get_values
   Index.set_value
   Index.isin
   Index.slice_indexer
   Index.slice_locs

.. _api.numericindex:

Numeric Index
-------------

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_autosummary.rst

   RangeIndex
   Int64Index
   UInt64Index
   Float64Index

.. We need this autosummary so that the methods are generated.
.. Separate block, since they aren't classes.

.. autosummary::
   :toctree: generated/

   RangeIndex.from_range


.. _api.categoricalindex:

CategoricalIndex
----------------

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_autosummary.rst

   CategoricalIndex

Categorical Components
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

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
   CategoricalIndex.map

.. _api.intervalindex:

IntervalIndex
-------------

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_autosummary.rst

   IntervalIndex

IntervalIndex Components
~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   IntervalIndex.from_arrays
   IntervalIndex.from_tuples
   IntervalIndex.from_breaks
   IntervalIndex.contains
   IntervalIndex.left
   IntervalIndex.right
   IntervalIndex.mid
   IntervalIndex.closed
   IntervalIndex.length
   IntervalIndex.values
   IntervalIndex.is_non_overlapping_monotonic
   IntervalIndex.get_loc
   IntervalIndex.get_indexer
   IntervalIndex.set_closed


.. _api.multiindex:

MultiIndex
----------

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_autosummary.rst

   MultiIndex

.. autosummary::
   :toctree: generated/

   IndexSlice

MultiIndex Constructors
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   MultiIndex.from_arrays
   MultiIndex.from_tuples
   MultiIndex.from_product

MultiIndex Attributes
~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   MultiIndex.names
   MultiIndex.levels
   MultiIndex.labels
   MultiIndex.nlevels
   MultiIndex.levshape

MultiIndex Components
~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   MultiIndex.set_levels
   MultiIndex.set_labels
   MultiIndex.to_hierarchical
   MultiIndex.to_frame
   MultiIndex.is_lexsorted
   MultiIndex.sortlevel
   MultiIndex.droplevel
   MultiIndex.swaplevel
   MultiIndex.reorder_levels
   MultiIndex.remove_unused_levels
   MultiIndex.unique

MultiIndex Selecting
~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   MultiIndex.get_loc
   MultiIndex.get_indexer
   MultiIndex.get_level_values

.. _api.datetimeindex:

DatetimeIndex
-------------

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_autosummary.rst

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
   :toctree: generated/

   DatetimeIndex.indexer_at_time
   DatetimeIndex.indexer_between_time


Time-specific operations
~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DatetimeIndex.normalize
   DatetimeIndex.strftime
   DatetimeIndex.snap
   DatetimeIndex.tz_convert
   DatetimeIndex.tz_localize
   DatetimeIndex.is_dst
   DatetimeIndex.round
   DatetimeIndex.floor
   DatetimeIndex.ceil
   DatetimeIndex.month_name
   DatetimeIndex.day_name

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DatetimeIndex.to_period
   DatetimeIndex.to_perioddelta
   DatetimeIndex.to_pydatetime
   DatetimeIndex.to_series
   DatetimeIndex.to_frame

TimedeltaIndex
--------------

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_autosummary.rst

   TimedeltaIndex

Components
~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   TimedeltaIndex.days
   TimedeltaIndex.seconds
   TimedeltaIndex.microseconds
   TimedeltaIndex.nanoseconds
   TimedeltaIndex.components
   TimedeltaIndex.inferred_freq

Conversion
~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   TimedeltaIndex.to_pytimedelta
   TimedeltaIndex.to_series
   TimedeltaIndex.round
   TimedeltaIndex.floor
   TimedeltaIndex.ceil
   TimedeltaIndex.to_frame

.. currentmodule:: pandas

PeriodIndex
--------------

.. autosummary::
   :toctree: generated/
   :template: autosummary/class_without_autosummary.rst

   PeriodIndex

Attributes
~~~~~~~~~~
.. autosummary::
    :toctree: generated/

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
    :toctree: generated/

    PeriodIndex.asfreq
    PeriodIndex.strftime
    PeriodIndex.to_timestamp

Scalars
-------

Period
~~~~~~
.. autosummary::
    :toctree: generated/

    Period

Attributes
~~~~~~~~~~
.. autosummary::
    :toctree: generated/

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
    :toctree: generated/

    Period.asfreq
    Period.now
    Period.strftime
    Period.to_timestamp

Timestamp
~~~~~~~~~
.. autosummary::
    :toctree: generated/

    Timestamp

Properties
~~~~~~~~~~
.. autosummary::
    :toctree: generated/

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
    :toctree: generated/

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

Interval
~~~~~~~~
.. autosummary::
    :toctree: generated/

    Interval

Properties
~~~~~~~~~~
.. autosummary::
    :toctree: generated/

    Interval.closed
    Interval.closed_left
    Interval.closed_right
    Interval.left
    Interval.length
    Interval.mid
    Interval.open_left
    Interval.open_right
    Interval.right

Timedelta
~~~~~~~~~
.. autosummary::
    :toctree: generated/

    Timedelta

Properties
~~~~~~~~~~
.. autosummary::
    :toctree: generated/

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
    :toctree: generated/

    Timedelta.ceil
    Timedelta.floor
    Timedelta.isoformat
    Timedelta.round
    Timedelta.to_pytimedelta
    Timedelta.to_timedelta64
    Timedelta.total_seconds

.. _api.frequencies:

Frequencies
-----------

.. currentmodule:: pandas.tseries.frequencies

.. _api.offsets:

.. autosummary::
   :toctree: generated/

   to_offset


Window
------

.. currentmodule:: pandas.core.window

Rolling objects are returned by ``.rolling`` calls: :func:`pandas.DataFrame.rolling`, :func:`pandas.Series.rolling`, etc.
Expanding objects are returned by ``.expanding`` calls: :func:`pandas.DataFrame.expanding`, :func:`pandas.Series.expanding`, etc.
EWM objects are returned by ``.ewm`` calls: :func:`pandas.DataFrame.ewm`, :func:`pandas.Series.ewm`, etc.

Standard moving window functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: pandas.core.window

.. autosummary::
   :toctree: generated/

   Rolling.count
   Rolling.sum
   Rolling.mean
   Rolling.median
   Rolling.var
   Rolling.std
   Rolling.min
   Rolling.max
   Rolling.corr
   Rolling.cov
   Rolling.skew
   Rolling.kurt
   Rolling.apply
   Rolling.aggregate
   Rolling.quantile
   Window.mean
   Window.sum

.. _api.functions_expanding:

Standard expanding window functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: pandas.core.window

.. autosummary::
   :toctree: generated/

   Expanding.count
   Expanding.sum
   Expanding.mean
   Expanding.median
   Expanding.var
   Expanding.std
   Expanding.min
   Expanding.max
   Expanding.corr
   Expanding.cov
   Expanding.skew
   Expanding.kurt
   Expanding.apply
   Expanding.aggregate
   Expanding.quantile

Exponentially-weighted moving window functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: pandas.core.window

.. autosummary::
   :toctree: generated/

   EWM.mean
   EWM.std
   EWM.var
   EWM.corr
   EWM.cov

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
   :template: autosummary/class_without_autosummary.rst

   Grouper

.. currentmodule:: pandas.core.groupby

Function application
~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   GroupBy.apply
   GroupBy.aggregate
   GroupBy.transform
   GroupBy.pipe

Computations / Descriptive Stats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   GroupBy.all
   GroupBy.any
   GroupBy.bfill
   GroupBy.count
   GroupBy.cumcount
   GroupBy.ffill
   GroupBy.first
   GroupBy.head
   GroupBy.last
   GroupBy.max
   GroupBy.mean
   GroupBy.median
   GroupBy.min
   GroupBy.ngroup
   GroupBy.nth
   GroupBy.ohlc
   GroupBy.prod
   GroupBy.rank
   GroupBy.pct_change
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

   DataFrameGroupBy.agg
   DataFrameGroupBy.all
   DataFrameGroupBy.any
   DataFrameGroupBy.bfill
   DataFrameGroupBy.corr
   DataFrameGroupBy.count
   DataFrameGroupBy.cov
   DataFrameGroupBy.cummax
   DataFrameGroupBy.cummin
   DataFrameGroupBy.cumprod
   DataFrameGroupBy.cumsum
   DataFrameGroupBy.describe
   DataFrameGroupBy.diff
   DataFrameGroupBy.ffill
   DataFrameGroupBy.fillna
   DataFrameGroupBy.filter
   DataFrameGroupBy.hist
   DataFrameGroupBy.idxmax
   DataFrameGroupBy.idxmin
   DataFrameGroupBy.mad
   DataFrameGroupBy.pct_change
   DataFrameGroupBy.plot
   DataFrameGroupBy.quantile
   DataFrameGroupBy.rank
   DataFrameGroupBy.resample
   DataFrameGroupBy.shift
   DataFrameGroupBy.size
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
   SeriesGroupBy.is_monotonic_increasing
   SeriesGroupBy.is_monotonic_decreasing

The following methods are available only for ``DataFrameGroupBy`` objects.

.. autosummary::
   :toctree: generated/

   DataFrameGroupBy.corrwith
   DataFrameGroupBy.boxplot

Resampling
----------
.. currentmodule:: pandas.core.resample

Resampler objects are returned by resample calls: :func:`pandas.DataFrame.resample`, :func:`pandas.Series.resample`.

Indexing, iteration
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Resampler.__iter__
   Resampler.groups
   Resampler.indices
   Resampler.get_group

Function application
~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Resampler.apply
   Resampler.aggregate
   Resampler.transform
   Resampler.pipe

Upsampling
~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Resampler.ffill
   Resampler.backfill
   Resampler.bfill
   Resampler.pad
   Resampler.nearest
   Resampler.fillna
   Resampler.asfreq
   Resampler.interpolate

Computations / Descriptive Stats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Resampler.count
   Resampler.nunique
   Resampler.first
   Resampler.last
   Resampler.max
   Resampler.mean
   Resampler.median
   Resampler.min
   Resampler.ohlc
   Resampler.prod
   Resampler.size
   Resampler.sem
   Resampler.std
   Resampler.sum
   Resampler.var
   Resampler.quantile

Style
-----
.. currentmodule:: pandas.io.formats.style

``Styler`` objects are returned by :attr:`pandas.DataFrame.style`.

Styler Constructor
~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Styler
   Styler.from_custom_template


Styler Attributes
~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Styler.env
   Styler.template
   Styler.loader

Style Application
~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Styler.apply
   Styler.applymap
   Styler.where
   Styler.format
   Styler.set_precision
   Styler.set_table_styles
   Styler.set_table_attributes
   Styler.set_caption
   Styler.set_properties
   Styler.set_uuid
   Styler.clear

Builtin Styles
~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Styler.highlight_max
   Styler.highlight_min
   Styler.highlight_null
   Styler.background_gradient
   Styler.bar

Style Export and Import
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Styler.render
   Styler.export
   Styler.use
   Styler.to_excel

Plotting
--------

.. currentmodule:: pandas.plotting

The following functions are contained in the `pandas.plotting` module.

.. autosummary::
   :toctree: generated/

   andrews_curves
   bootstrap_plot
   deregister_matplotlib_converters
   lag_plot
   parallel_coordinates
   radviz
   register_matplotlib_converters
   scatter_matrix

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

Testing functions
~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   testing.assert_frame_equal
   testing.assert_series_equal
   testing.assert_index_equal


Exceptions and warnings
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   errors.DtypeWarning
   errors.EmptyDataError
   errors.OutOfBoundsDatetime
   errors.ParserError
   errors.ParserWarning
   errors.PerformanceWarning
   errors.UnsortedIndexError
   errors.UnsupportedFunctionCall


Data types related functionality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   api.types.union_categoricals
   api.types.infer_dtype
   api.types.pandas_dtype

Dtype introspection

.. autosummary::
   :toctree: generated/

    api.types.is_bool_dtype
    api.types.is_categorical_dtype
    api.types.is_complex_dtype
    api.types.is_datetime64_any_dtype
    api.types.is_datetime64_dtype
    api.types.is_datetime64_ns_dtype
    api.types.is_datetime64tz_dtype
    api.types.is_extension_type
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

.. autosummary::
   :toctree: generated/

    api.types.is_dict_like
    api.types.is_file_like
    api.types.is_list_like
    api.types.is_named_tuple
    api.types.is_iterator

Scalar introspection

.. autosummary::
   :toctree: generated/

    api.types.is_bool
    api.types.is_categorical
    api.types.is_complex
    api.types.is_datetimetz
    api.types.is_float
    api.types.is_hashable
    api.types.is_integer
    api.types.is_interval
    api.types.is_number
    api.types.is_period
    api.types.is_re
    api.types.is_re_compilable
    api.types.is_scalar

Extensions
----------

These are primarily intended for library authors looking to extend pandas
objects.

.. currentmodule:: pandas

.. autosummary::
   :toctree: generated/

   api.extensions.register_dataframe_accessor
   api.extensions.register_series_accessor
   api.extensions.register_index_accessor
   api.extensions.ExtensionDtype
   api.extensions.ExtensionArray

.. This is to prevent warnings in the doc build. We don't want to encourage
.. these methods.

.. toctree::
   :hidden:

   generated/pandas.DataFrame.blocks
   generated/pandas.DataFrame.as_matrix
   generated/pandas.DataFrame.ix
   generated/pandas.Index.asi8
   generated/pandas.Index.data
   generated/pandas.Index.flags
   generated/pandas.Index.holds_integer
   generated/pandas.Index.is_type_compatible
   generated/pandas.Index.nlevels
   generated/pandas.Index.sort
   generated/pandas.Panel.agg
   generated/pandas.Panel.aggregate
   generated/pandas.Panel.blocks
   generated/pandas.Panel.empty
   generated/pandas.Panel.is_copy
   generated/pandas.Panel.items
   generated/pandas.Panel.ix
   generated/pandas.Panel.major_axis
   generated/pandas.Panel.minor_axis
   generated/pandas.Series.asobject
   generated/pandas.Series.blocks
   generated/pandas.Series.from_array
   generated/pandas.Series.ix
   generated/pandas.Series.imag
   generated/pandas.Series.real
