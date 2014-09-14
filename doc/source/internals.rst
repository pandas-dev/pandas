.. _internals:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import random
   np.random.seed(123456)
   from pandas import *
   options.display.max_rows=15
   import pandas as pd
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)
   from pandas.compat import range, zip

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

These are range generates to make the creation of a regular index easy:

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

   index = MultiIndex.from_product([range(3), ['one', 'two']], names=['first', 'second'])
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


