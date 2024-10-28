{{ header }}

.. _api.groupby:

=======
GroupBy
=======
.. currentmodule:: pandas.core.groupby

:class:`pandas.api.typing.DataFrameGroupBy` and :class:`pandas.api.typing.SeriesGroupBy`
instances are returned by groupby calls :func:`pandas.DataFrame.groupby` and
:func:`pandas.Series.groupby` respectively.

Indexing, iteration
-------------------
.. autosummary::
   :toctree: api/

   DataFrameGroupBy.__iter__
   SeriesGroupBy.__iter__
   DataFrameGroupBy.groups
   SeriesGroupBy.groups
   DataFrameGroupBy.indices
   SeriesGroupBy.indices
   DataFrameGroupBy.get_group
   SeriesGroupBy.get_group

.. currentmodule:: pandas

.. autosummary::
   :toctree: api/
   :template: autosummary/class_without_autosummary.rst

   Grouper

Function application helper
---------------------------
.. autosummary::
   :toctree: api/

   NamedAgg

.. currentmodule:: pandas.core.groupby

Function application
--------------------
.. autosummary::
   :toctree: api/

   SeriesGroupBy.apply
   DataFrameGroupBy.apply
   SeriesGroupBy.agg
   DataFrameGroupBy.agg
   SeriesGroupBy.aggregate
   DataFrameGroupBy.aggregate
   SeriesGroupBy.transform
   DataFrameGroupBy.transform
   SeriesGroupBy.pipe
   DataFrameGroupBy.pipe
   DataFrameGroupBy.filter
   SeriesGroupBy.filter

``DataFrameGroupBy`` computations / descriptive stats
-----------------------------------------------------
.. autosummary::
   :toctree: api/

   DataFrameGroupBy.all
   DataFrameGroupBy.any
   DataFrameGroupBy.bfill
   DataFrameGroupBy.corr
   DataFrameGroupBy.corrwith
   DataFrameGroupBy.count
   DataFrameGroupBy.cov
   DataFrameGroupBy.cumcount
   DataFrameGroupBy.cummax
   DataFrameGroupBy.cummin
   DataFrameGroupBy.cumprod
   DataFrameGroupBy.cumsum
   DataFrameGroupBy.describe
   DataFrameGroupBy.diff
   DataFrameGroupBy.ffill
   DataFrameGroupBy.first
   DataFrameGroupBy.head
   DataFrameGroupBy.idxmax
   DataFrameGroupBy.idxmin
   DataFrameGroupBy.last
   DataFrameGroupBy.max
   DataFrameGroupBy.mean
   DataFrameGroupBy.median
   DataFrameGroupBy.min
   DataFrameGroupBy.ngroup
   DataFrameGroupBy.nth
   DataFrameGroupBy.nunique
   DataFrameGroupBy.ohlc
   DataFrameGroupBy.pct_change
   DataFrameGroupBy.prod
   DataFrameGroupBy.quantile
   DataFrameGroupBy.rank
   DataFrameGroupBy.resample
   DataFrameGroupBy.rolling
   DataFrameGroupBy.sample
   DataFrameGroupBy.sem
   DataFrameGroupBy.shift
   DataFrameGroupBy.size
   DataFrameGroupBy.skew
   DataFrameGroupBy.std
   DataFrameGroupBy.sum
   DataFrameGroupBy.var
   DataFrameGroupBy.tail
   DataFrameGroupBy.take
   DataFrameGroupBy.value_counts

``SeriesGroupBy`` computations / descriptive stats
--------------------------------------------------
.. autosummary::
   :toctree: api/

   SeriesGroupBy.all
   SeriesGroupBy.any
   SeriesGroupBy.bfill
   SeriesGroupBy.corr
   SeriesGroupBy.count
   SeriesGroupBy.cov
   SeriesGroupBy.cumcount
   SeriesGroupBy.cummax
   SeriesGroupBy.cummin
   SeriesGroupBy.cumprod
   SeriesGroupBy.cumsum
   SeriesGroupBy.describe
   SeriesGroupBy.diff
   SeriesGroupBy.ffill
   SeriesGroupBy.first
   SeriesGroupBy.head
   SeriesGroupBy.last
   SeriesGroupBy.idxmax
   SeriesGroupBy.idxmin
   SeriesGroupBy.is_monotonic_increasing
   SeriesGroupBy.is_monotonic_decreasing
   SeriesGroupBy.max
   SeriesGroupBy.mean
   SeriesGroupBy.median
   SeriesGroupBy.min
   SeriesGroupBy.ngroup
   SeriesGroupBy.nlargest
   SeriesGroupBy.nsmallest
   SeriesGroupBy.nth
   SeriesGroupBy.nunique
   SeriesGroupBy.unique
   SeriesGroupBy.ohlc
   SeriesGroupBy.pct_change
   SeriesGroupBy.prod
   SeriesGroupBy.quantile
   SeriesGroupBy.rank
   SeriesGroupBy.resample
   SeriesGroupBy.rolling
   SeriesGroupBy.sample
   SeriesGroupBy.sem
   SeriesGroupBy.shift
   SeriesGroupBy.size
   SeriesGroupBy.skew
   SeriesGroupBy.std
   SeriesGroupBy.sum
   SeriesGroupBy.var
   SeriesGroupBy.tail
   SeriesGroupBy.take
   SeriesGroupBy.value_counts

Plotting and visualization
--------------------------
.. autosummary::
   :toctree: api/

   DataFrameGroupBy.boxplot
   DataFrameGroupBy.hist
   SeriesGroupBy.hist
   DataFrameGroupBy.plot
   SeriesGroupBy.plot
