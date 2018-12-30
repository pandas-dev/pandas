{{ header }}

.. _api.panel:

=====
Panel
=====
.. currentmodule:: pandas

Constructor
~~~~~~~~~~~
.. autosummary::
   :toctree: generated/

   Panel

Properties and underlying data
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
