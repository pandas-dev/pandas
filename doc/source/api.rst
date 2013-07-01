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

.. currentmodule:: pandas.io.pickle

.. autosummary::
   :toctree: generated/

   read_pickle

Flat File
~~~~~~~~~

.. currentmodule:: pandas.io.parsers

.. autosummary::
   :toctree: generated/

   read_table
   read_csv
   read_fwf
   read_clipboard

Excel
~~~~~

.. currentmodule:: pandas.io.excel

.. autosummary::
   :toctree: generated/

   read_excel
   ExcelFile.parse

JSON
~~~~

.. currentmodule:: pandas.io.json

.. autosummary::
   :toctree: generated/

   read_json

HTML
~~~~

.. currentmodule:: pandas.io.html

.. autosummary::
   :toctree: generated/

   read_html

HDFStore: PyTables (HDF5)
~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: pandas.io.pytables

.. autosummary::
   :toctree: generated/

   read_hdf
   HDFStore.put
   HDFStore.append
   HDFStore.get
   HDFStore.select

SQL
~~~
.. currentmodule:: pandas.io.sql

.. autosummary::
   :toctree: generated/

   read_sql
   read_frame
   write_frame


STATA
~~~~~

.. currentmodule:: pandas.io.stata

.. autosummary::
   :toctree: generated/

   read_stata
   StataReader.data
   StataReader.data_label
   StataReader.value_labels
   StataReader.variable_labels
   StataWriter.write_file


General functions
-----------------

Data manipulations
~~~~~~~~~~~~~~~~~~
.. currentmodule:: pandas.tools.pivot

.. autosummary::
   :toctree: generated/

   pivot_table

.. currentmodule:: pandas.tools.merge

.. autosummary::
   :toctree: generated/

   merge
   concat

Top-level missing data
~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: pandas.core.common

.. autosummary::
   :toctree: generated/

   isnull
   notnull

Top-level dealing with datetimes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. currentmodule:: pandas.tseries.tools

.. autosummary::
   :toctree: generated/

   to_datetime


Standard moving window functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Standard expanding window functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Exponentially-weighted moving window functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   ewma
   ewmstd
   ewmvar
   ewmcorr
   ewmcov

.. currentmodule:: pandas

.. _api.series:

Series
------

Attributes and underlying data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Axes**
  * **index**: axis labels

.. autosummary::
   :toctree: generated/

   Series.values
   Series.dtype
   Series.isnull
   Series.notnull

Conversion / Constructors
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Series.__init__
   Series.astype
   Series.copy

Indexing, iteration
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.get
   Series.ix
   Series.__iter__
   Series.iteritems

Binary operator functions
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.add
   Series.div
   Series.mul
   Series.sub
   Series.combine
   Series.combine_first
   Series.round

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
   Series.kurt
   Series.mad
   Series.max
   Series.mean
   Series.median
   Series.min
   Series.nunique
   Series.pct_change
   Series.prod
   Series.quantile
   Series.rank
   Series.skew
   Series.std
   Series.sum
   Series.unique
   Series.var
   Series.value_counts

Reindexing / Selection / Label manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.align
   Series.drop
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
   Series.weekday
   Series.resample
   Series.tz_convert
   Series.tz_localize

Plotting
~~~~~~~~
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
   Series.to_sparse
   Series.to_string
   Series.to_clipboard

.. _api.dataframe:

DataFrame
---------

Attributes and underlying data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**Axes**

  * **index**: row labels
  * **columns**: column labels

.. autosummary::
   :toctree: generated/

   DataFrame.as_matrix
   DataFrame.dtypes
   DataFrame.get_dtype_counts
   DataFrame.values
   DataFrame.axes
   DataFrame.ndim
   DataFrame.shape

Conversion / Constructors
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.__init__
   DataFrame.astype
   DataFrame.convert_objects
   DataFrame.copy

Indexing, iteration
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.head
   DataFrame.ix
   DataFrame.insert
   DataFrame.__iter__
   DataFrame.iteritems
   DataFrame.iterrows
   DataFrame.itertuples
   DataFrame.lookup
   DataFrame.pop
   DataFrame.tail
   DataFrame.xs

Binary operator functions
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.add
   DataFrame.div
   DataFrame.mul
   DataFrame.sub
   DataFrame.radd
   DataFrame.rdiv
   DataFrame.rmul
   DataFrame.rsub
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
   DataFrame.kurt
   DataFrame.mad
   DataFrame.max
   DataFrame.mean
   DataFrame.median
   DataFrame.min
   DataFrame.pct_change
   DataFrame.prod
   DataFrame.quantile
   DataFrame.rank
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

   DataFrame.delevel
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
   DataFrame.to_dict
   DataFrame.to_excel
   DataFrame.to_json
   DataFrame.to_html
   DataFrame.to_stata
   DataFrame.to_records
   DataFrame.to_sparse
   DataFrame.to_string
   DataFrame.to_clipboard

.. _api.panel:

Panel
------

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
   Panel.shape

Conversion / Constructors
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.__init__
   Panel.astype
   Panel.copy

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

   Panel.ix
   Panel.__iter__
   Panel.iteritems
   Panel.pop
   Panel.xs
   Panel.major_xs
   Panel.minor_xs

Binary operator functions
~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel.add
   Panel.div
   Panel.mul
   Panel.sub

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
   Panel.filter
   Panel.first
   Panel.last
   Panel.reindex
   Panel.reindex_axis
   Panel.reindex_like
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
   Panel.to_sparse
   Panel.to_frame
   Panel.to_clipboard

