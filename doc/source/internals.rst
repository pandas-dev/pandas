.. _internals:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   np.set_printoptions(precision=4, suppress=True)
   import pandas as pd
   pd.options.display.max_rows = 15

*********
Internals
*********

This section will provide a look into some of pandas internals.

Indexing
--------

In pandas there are a few objects implemented which can serve as valid
containers for the axis labels:

- ``Index``: the generic "ordered set" object, an ndarray of object dtype
  assuming nothing about its contents. The labels must be hashable (and
  likely immutable) and unique. Populates a dict of label to location in
  Cython to do ``O(1)`` lookups.
- ``Int64Index``: a version of ``Index`` highly optimized for 64-bit integer
  data, such as time stamps
- ``Float64Index``: a version of ``Index`` highly optimized for 64-bit float data
- ``MultiIndex``: the standard hierarchical index object
- ``DatetimeIndex``: An Index object with ``Timestamp`` boxed elements (impl are the int64 values)
- ``TimedeltaIndex``: An Index object with ``Timedelta`` boxed elements (impl are the in64 values)
- ``PeriodIndex``: An Index object with Period elements

There are functions that make the creation of a regular index easy:

- ``date_range``: fixed frequency date range generated from a time rule or
  DateOffset. An ndarray of Python datetime objects
- ``period_range``: fixed frequency date range generated from a time rule or
  DateOffset. An ndarray of ``Period`` objects, representing Timespans

The motivation for having an ``Index`` class in the first place was to enable
different implementations of indexing. This means that it's possible for you,
the user, to implement a custom ``Index`` subclass that may be better suited to
a particular application than the ones provided in pandas.

From an internal implementation point of view, the relevant methods that an
``Index`` must define are one or more of the following (depending on how
incompatible the new object internals are with the ``Index`` functions):

- ``get_loc``: returns an "indexer" (an integer, or in some cases a
  slice object) for a label
- ``slice_locs``: returns the "range" to slice between two labels
- ``get_indexer``: Computes the indexing vector for reindexing / data
  alignment purposes. See the source / docstrings for more on this
- ``get_indexer_non_unique``: Computes the indexing vector for reindexing / data
  alignment purposes when the index is non-unique. See the source / docstrings
  for more on this
- ``reindex``: Does any pre-conversion of the input index then calls
  ``get_indexer``
- ``union``, ``intersection``: computes the union or intersection of two
  Index objects
- ``insert``: Inserts a new label into an Index, yielding a new object
- ``delete``: Delete a label, yielding a new object
- ``drop``: Deletes a set of labels
- ``take``: Analogous to ndarray.take

MultiIndex
~~~~~~~~~~

Internally, the ``MultiIndex`` consists of a few things: the **levels**, the
integer **labels**, and the level **names**:

.. ipython:: python

   index = pd.MultiIndex.from_product([range(3), ['one', 'two']], names=['first', 'second'])
   index
   index.levels
   index.labels
   index.names

You can probably guess that the labels determine which unique element is
identified with that location at each layer of the index. It's important to
note that sortedness is determined **solely** from the integer labels and does
not check (or care) whether the levels themselves are sorted. Fortunately, the
constructors ``from_tuples`` and ``from_arrays`` ensure that this is true, but
if you compute the levels and labels yourself, please be careful.

.. _ref-subclassing-pandas:

Subclassing pandas Data Structures
----------------------------------

.. warning:: There are some easier alternatives before considering subclassing ``pandas`` data structures.

  1. Extensible method chains with :ref:`pipe <basics.pipe>`

  2. Use *composition*. See `here <http://en.wikipedia.org/wiki/Composition_over_inheritance>`_.

This section describes how to subclass ``pandas`` data structures to meet more specific needs. There are 2 points which need attention:

1. Override constructor properties.
2. Define original properties

.. note:: You can find a nice example in `geopandas <https://github.com/geopandas/geopandas>`_ project.

Override Constructor Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each data structure has constructor properties to specifying data constructors. By overriding these properties, you can retain defined-classes through ``pandas`` data manipulations.

There are 3 constructors to be defined:

- ``_constructor``: Used when a manipulation result has the same dimesions as the original.
- ``_constructor_sliced``: Used when a manipulation result has one lower dimension(s) as the original, such as ``DataFrame`` single columns slicing.
- ``_constructor_expanddim``: Used when a manipulation result has one higher dimension as the original, such as ``Series.to_frame()`` and ``DataFrame.to_panel()``.

Following table shows how ``pandas`` data structures define constructor properties by default.

===========================  ======================= =================== =======================
Property Attributes          ``Series``              ``DataFrame``       ``Panel``
===========================  ======================= =================== =======================
``_constructor``             ``Series``              ``DataFrame``       ``Panel``
``_constructor_sliced``      ``NotImplementedError`` ``Series``          ``DataFrame``
``_constructor_expanddim``   ``DataFrame``           ``Panel``           ``NotImplementedError``
===========================  ======================= =================== =======================

Below example shows how to define ``SubclassedSeries`` and ``SubclassedDataFrame`` overriding constructor properties.

.. code-block:: python

   class SubclassedSeries(Series):

       @property
       def _constructor(self):
           return SubclassedSeries

       @property
       def _constructor_expanddim(self):
           return SubclassedDataFrame

   class SubclassedDataFrame(DataFrame):

       @property
       def _constructor(self):
           return SubclassedDataFrame

       @property
       def _constructor_sliced(self):
           return SubclassedSeries

.. code-block:: python

   >>> s = SubclassedSeries([1, 2, 3])
   >>> type(s)
   <class '__main__.SubclassedSeries'>

   >>> to_framed = s.to_frame()
   >>> type(to_framed)
   <class '__main__.SubclassedDataFrame'>

   >>> df = SubclassedDataFrame({'A', [1, 2, 3], 'B': [4, 5, 6], 'C': [7, 8, 9]})
   >>> df
      A  B  C
   0  1  4  7
   1  2  5  8
   2  3  6  9

   >>> type(df)
   <class '__main__.SubclassedDataFrame'>

   >>> sliced1 = df[['A', 'B']]
   >>> sliced1
      A  B
   0  1  4
   1  2  5
   2  3  6
   >>> type(sliced1)
   <class '__main__.SubclassedDataFrame'>

   >>> sliced2 = df['A']
   >>> sliced2
   0    1
   1    2
   2    3
   Name: A, dtype: int64
   >>> type(sliced2)
   <class '__main__.SubclassedSeries'>

Define Original Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~

To let original data structures have additional properties, you should let ``pandas`` know what properties are added. ``pandas`` maps unknown properties to data names overriding ``__getattribute__``. Defining original properties can be done in one of 2 ways:

1. Define ``_internal_names`` and ``_internal_names_set`` for temporary properties which WILL NOT be passed to manipulation results.
2. Define ``_metadata`` for normal properties which will be passed to manipulation results.

Below is an example to define 2 original properties, "internal_cache" as a temporary property and "added_property" as a normal property

.. code-block:: python

   class SubclassedDataFrame2(DataFrame):

       # temporary properties
       _internal_names = pd.DataFrame._internal_names + ['internal_cache']
       _internal_names_set = set(_internal_names)

       # normal properties
       _metadata = ['added_property']

       @property
       def _constructor(self):
           return SubclassedDataFrame2

.. code-block:: python

   >>> df = SubclassedDataFrame2({'A', [1, 2, 3], 'B': [4, 5, 6], 'C': [7, 8, 9]})
   >>> df
      A  B  C
   0  1  4  7
   1  2  5  8
   2  3  6  9

   >>> df.internal_cache = 'cached'
   >>> df.added_property = 'property'

   >>> df.internal_cache
   cached
   >>> df.added_property
   property

   # properties defined in _internal_names is reset after manipulation
   >>> df[['A', 'B']].internal_cache
   AttributeError: 'SubclassedDataFrame2' object has no attribute 'internal_cache'

   # properties defined in _metadata are retained
   >>> df[['A', 'B']].added_property
   property
