.. currentmodule:: pandas
.. _api:

*************
API Reference
*************

.. _api.functions:

General functions
-----------------

Data manipulations
~~~~~~~~~~~~~~~~~~
.. currentmodule:: pandas

.. autosummary::
   :toctree: generated/

   pivot_table
   merge
   concat

Pickling
~~~~~~~~

.. currentmodule:: pandas.core.common

.. autosummary::
   :toctree: generated/

   load
   save

File IO
~~~~~~~

.. currentmodule:: pandas.io.parsers

.. autosummary::
   :toctree: generated/

   read_table
   read_csv
   ExcelFile.parse

HDFStore: PyTables (HDF5)
~~~~~~~~~~~~~~~~~~~~~~~~~
.. currentmodule:: pandas.io.pytables

.. autosummary::
   :toctree: generated/

   HDFStore.put
   HDFStore.get

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

   Series.autocorr
   Series.clip
   Series.clip_lower
   Series.clip_upper
   Series.corr
   Series.count
   Series.cumprod
   Series.cumsum
   Series.describe
   Series.diff
   Series.max
   Series.mean
   Series.median
   Series.min
   Series.prod
   Series.quantile
   Series.skew
   Series.std
   Series.sum
   Series.var
   Series.value_counts

Reindexing / Selection / Label manipulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.align
   Series.drop
   Series.reindex
   Series.reindex_like
   Series.rename
   Series.select
   Series.take
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
   Series.sort
   Series.sort_index
   Series.sortlevel
   Series.unstack

Combining / joining / merging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Series.append

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
   Series.load
   Series.save
   Series.to_csv
   Series.to_dict
   Series.to_sparse

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
   DataFrame.copy

Indexing, iteration
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.ix
   DataFrame.insert
   DataFrame.__iter__
   DataFrame.iteritems
   DataFrame.pop
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

   DataFrame.clip
   DataFrame.clip_lower
   DataFrame.clip_upper
   DataFrame.corr
   DataFrame.corrwith
   DataFrame.count
   DataFrame.cumprod
   DataFrame.cumsum
   DataFrame.describe
   DataFrame.diff
   DataFrame.mad
   DataFrame.max
   DataFrame.mean
   DataFrame.median
   DataFrame.min
   DataFrame.prod
   DataFrame.quantile
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
   DataFrame.filter
   DataFrame.reindex
   DataFrame.reindex_like
   DataFrame.rename
   DataFrame.select
   DataFrame.take
   DataFrame.truncate
   DataFrame.head
   DataFrame.tail

.. _api.dataframe.missing:

Missing data handling
~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.dropna
   DataFrame.fillna

Reshaping, sorting, transposing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.sort_index
   DataFrame.delevel
   DataFrame.pivot
   DataFrame.sortlevel
   DataFrame.swaplevel
   DataFrame.stack
   DataFrame.unstack
   DataFrame.T
   DataFrame.transpose

Combining / joining / merging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.join
   DataFrame.merge
   DataFrame.append

Time series-related
~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.asfreq
   DataFrame.shift
   DataFrame.first_valid_index
   DataFrame.last_valid_index

Plotting
~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.hist
   DataFrame.plot

Serialization / IO / Conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.from_csv
   DataFrame.from_records
   DataFrame.to_csv
   DataFrame.to_dict
   DataFrame.to_records
   DataFrame.to_sparse
   DataFrame.to_string
   DataFrame.save
   DataFrame.load
   DataFrame.info

.. _api.panel:

Panel
-----

.. _api.panel.stats:

Computations / Descriptive Stats
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

