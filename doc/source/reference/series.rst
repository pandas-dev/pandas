{{ header }}

.. _api.series:

======
Series
======
.. currentmodule:: pandas

Constructor
-----------
.. autosummary::
   :toctree: api/

   Series

Attributes
----------
**Axes**

.. autosummary::
   :toctree: api/

   Series.index
   Series.array
   Series.values
   Series.dtype
   Series.info
   Series.shape
   Series.nbytes
   Series.ndim
   Series.size
   Series.T
   Series.memory_usage
   Series.hasnans
   Series.empty
   Series.dtypes
   Series.name
   Series.flags
   Series.set_flags

Conversion
----------
.. autosummary::
   :toctree: api/

   Series.astype
   Series.convert_dtypes
   Series.infer_objects
   Series.copy
   Series.to_numpy
   Series.to_period
   Series.to_timestamp
   Series.to_list
   Series.__array__

Indexing, iteration
-------------------
.. autosummary::
   :toctree: api/

   Series.get
   Series.at
   Series.iat
   Series.loc
   Series.iloc
   Series.__iter__
   Series.items
   Series.keys
   Series.pop
   Series.item
   Series.xs

For more information on ``.at``, ``.iat``, ``.loc``, and
``.iloc``,  see the :ref:`indexing documentation <indexing>`.

Binary operator functions
-------------------------
.. autosummary::
   :toctree: api/

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

Function application, GroupBy & window
--------------------------------------
.. autosummary::
   :toctree: api/

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

Computations / descriptive stats
--------------------------------
.. autosummary::
   :toctree: api/

   Series.abs
   Series.all
   Series.any
   Series.autocorr
   Series.between
   Series.clip
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
   Series.is_monotonic_increasing
   Series.is_monotonic_decreasing
   Series.value_counts

Reindexing / selection / label manipulation
-------------------------------------------
.. autosummary::
   :toctree: api/

   Series.align
   Series.case_when
   Series.drop
   Series.droplevel
   Series.drop_duplicates
   Series.duplicated
   Series.equals
   Series.head
   Series.idxmax
   Series.idxmin
   Series.isin
   Series.reindex
   Series.reindex_like
   Series.rename
   Series.rename_axis
   Series.reset_index
   Series.sample
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
---------------------
.. autosummary::
   :toctree: api/

   Series.bfill
   Series.dropna
   Series.ffill
   Series.fillna
   Series.interpolate
   Series.isna
   Series.isnull
   Series.notna
   Series.notnull
   Series.replace

Reshaping, sorting
------------------
.. autosummary::
   :toctree: api/

   Series.argsort
   Series.argmin
   Series.argmax
   Series.reorder_levels
   Series.sort_values
   Series.sort_index
   Series.swaplevel
   Series.unstack
   Series.explode
   Series.searchsorted
   Series.repeat
   Series.squeeze

Combining / comparing / joining / merging
-----------------------------------------
.. autosummary::
   :toctree: api/

   Series.compare
   Series.update

Time Series-related
-------------------
.. autosummary::
   :toctree: api/

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

Accessors
---------

pandas provides dtype-specific methods under various accessors.
These are separate namespaces within :class:`Series` that only apply
to specific data types.

.. autosummary::
    :toctree: api/
    :nosignatures:
    :template: autosummary/accessor.rst

    Series.str
    Series.cat
    Series.dt
    Series.sparse
    DataFrame.sparse
    Index.str


=========================== =================================
Data Type                   Accessor
=========================== =================================
Datetime, Timedelta, Period :ref:`dt <api.series.dt>`
String                      :ref:`str <api.series.str>`
Categorical                 :ref:`cat <api.series.cat>`
Sparse                      :ref:`sparse <api.series.sparse>`
=========================== =================================

.. _api.series.dt:

Datetimelike properties
~~~~~~~~~~~~~~~~~~~~~~~

``Series.dt`` can be used to access the values of the series as
datetimelike and return several properties.
These can be accessed like ``Series.dt.<property>``.

Datetime properties
^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: api/
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
   Series.dt.dayofweek
   Series.dt.day_of_week
   Series.dt.weekday
   Series.dt.dayofyear
   Series.dt.day_of_year
   Series.dt.days_in_month
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
   Series.dt.unit

Datetime methods
^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_method.rst

   Series.dt.isocalendar
   Series.dt.to_period
   Series.dt.to_pydatetime
   Series.dt.tz_localize
   Series.dt.tz_convert
   Series.dt.normalize
   Series.dt.strftime
   Series.dt.round
   Series.dt.floor
   Series.dt.ceil
   Series.dt.month_name
   Series.dt.day_name
   Series.dt.as_unit

Period properties
^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_attribute.rst

   Series.dt.qyear
   Series.dt.start_time
   Series.dt.end_time

Timedelta properties
^^^^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_attribute.rst

   Series.dt.days
   Series.dt.seconds
   Series.dt.microseconds
   Series.dt.nanoseconds
   Series.dt.components
   Series.dt.unit

Timedelta methods
^^^^^^^^^^^^^^^^^

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_method.rst

   Series.dt.to_pytimedelta
   Series.dt.total_seconds
   Series.dt.as_unit


.. _api.series.str:

String handling
~~~~~~~~~~~~~~~

``Series.str`` can be used to access the values of the series as
strings and apply several methods to it. These can be accessed like
``Series.str.<function/property>``.

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_method.rst

   Series.str.capitalize
   Series.str.casefold
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
   Series.str.fullmatch
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
   Series.str.removeprefix
   Series.str.removesuffix
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

.. _api.series.cat:

Categorical accessor
~~~~~~~~~~~~~~~~~~~~

Categorical-dtype specific methods and attributes are available under
the ``Series.cat`` accessor.

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_attribute.rst

   Series.cat.categories
   Series.cat.ordered
   Series.cat.codes

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_method.rst

   Series.cat.rename_categories
   Series.cat.reorder_categories
   Series.cat.add_categories
   Series.cat.remove_categories
   Series.cat.remove_unused_categories
   Series.cat.set_categories
   Series.cat.as_ordered
   Series.cat.as_unordered


.. _api.series.sparse:

Sparse accessor
~~~~~~~~~~~~~~~

Sparse-dtype specific methods and attributes are provided under the
``Series.sparse`` accessor.

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_attribute.rst

   Series.sparse.npoints
   Series.sparse.density
   Series.sparse.fill_value
   Series.sparse.sp_values

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_method.rst

   Series.sparse.from_coo
   Series.sparse.to_coo


.. _api.series.list:

List accessor
~~~~~~~~~~~~~

Arrow list-dtype specific methods and attributes are provided under the
``Series.list`` accessor.

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_method.rst

   Series.list.flatten
   Series.list.len
   Series.list.__getitem__


.. _api.series.struct:

Struct accessor
~~~~~~~~~~~~~~~

Arrow struct-dtype specific methods and attributes are provided under the
``Series.struct`` accessor.

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_attribute.rst

   Series.struct.dtypes

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_method.rst

   Series.struct.field
   Series.struct.explode


.. _api.series.flags:

Flags
~~~~~

Flags refer to attributes of the pandas object. Properties of the dataset (like
the date is was recorded, the URL it was accessed from, etc.) should be stored
in :attr:`Series.attrs`.

.. autosummary::
   :toctree: api/

   Flags

.. _api.series.metadata:

Metadata
~~~~~~~~

:attr:`Series.attrs` is a dictionary for storing global metadata for this Series.

.. warning:: ``Series.attrs`` is considered experimental and may change without warning.

.. autosummary::
   :toctree: api/

   Series.attrs


Plotting
--------
``Series.plot`` is both a callable method and a namespace attribute for
specific plotting methods of the form ``Series.plot.<kind>``.

.. autosummary::
   :toctree: api/
   :template: autosummary/accessor_callable.rst

   Series.plot

.. autosummary::
   :toctree: api/
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
   :toctree: api/

   Series.hist

Serialization / IO / conversion
-------------------------------
.. autosummary::
   :toctree: api/

   Series.to_pickle
   Series.to_csv
   Series.to_dict
   Series.to_excel
   Series.to_frame
   Series.to_xarray
   Series.to_hdf
   Series.to_sql
   Series.to_json
   Series.to_string
   Series.to_clipboard
   Series.to_latex
   Series.to_markdown
