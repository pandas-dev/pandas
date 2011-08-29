.. currentmodule:: pandas
.. _api:

***********
API Listing
***********

Series
------

DataFrame
---------

Attributes
~~~~~~~~~~
**Axes**

  * **index**:
  * **columns**:

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

Function application, GroupBy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.apply
   DataFrame.applymap
   DataFrame.groupby

Data / attributes
~~~~~~~~~~~~~~~~~
.. .. autosummary::
..    :toctree: generated/

..    DataFrame.as_matrix
..    DataFrame.values
..    DataFrame.axes
..    DataFrame.columns
..    DataFrame.ndim
..    DataFrame.index
..    DataFrame.shape

.. autosummary::
   :toctree: generated/

   DataFrame.as_matrix
   DataFrame.values
   DataFrame.axes
   DataFrame.ndim
   DataFrame.shape

Alignment / Selection / Reindexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.drop
   DataFrame.dropna
   DataFrame.fillna
   DataFrame.filter
   DataFrame.rename
   DataFrame.reindex
   DataFrame.reindex_like
   DataFrame.select
   DataFrame.take
   DataFrame.truncate
   DataFrame.head
   DataFrame.tail

Reshaping, sorting, transposing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.sort
   DataFrame.delevel
   DataFrame.pivot
   DataFrame.sortlevel
   DataFrame.stack
   DataFrame.unstack
   DataFrame.T
   DataFrame.transpose

Computations
~~~~~~~~~~~~
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
   DataFrame.product
   DataFrame.quantile
   DataFrame.skew
   DataFrame.sum
   DataFrame.std
   DataFrame.var

Combining / joining / merging
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   DataFrame.join
   DataFrame.append
   DataFrame.combineAdd
   DataFrame.combine_first
   DataFrame.combineMult

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

.. Deprecated (death row)
.. ~~~~~~~~~~~~~~~~~~~~~~
.. .. autosummary::
..    :toctree: generated/

..    DataFrame.asMatrix
..    DataFrame.cols
..    DataFrame.combineFirst
..    DataFrame.dropEmptyRows
..    DataFrame.dropIncompleteRows
..    DataFrame.getXS
..    DataFrame.merge
..    DataFrame.rows
..    DataFrame.fromRecords
..    DataFrame.fromcsv
..    DataFrame.tapply
..    DataFrame.tgroupby
..    DataFrame.toRecords
..    DataFrame.toCSV
..    DataFrame.toDataMatrix
..    DataFrame.toDict
..    DataFrame.fromRecords

WidePanel
---------

Input / Output
--------------

File IO
~~~~~~~

HDFStore: PyTables (HDF5)
~~~~~~~~~~~~~~~~~~~~~~~~~

Moving window statistics
------------------------

