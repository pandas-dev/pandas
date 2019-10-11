.. _internals:

{{ header }}

*********
Internals
*********

This section will provide a look into some of pandas internals. It's primarily
intended for developers of pandas itself.

Indexing
--------

In pandas there are a few objects implemented which can serve as valid
containers for the axis labels:

* ``Index``: the generic "ordered set" object, an ndarray of object dtype
  assuming nothing about its contents. The labels must be hashable (and
  likely immutable) and unique. Populates a dict of label to location in
  Cython to do ``O(1)`` lookups.
* ``Int64Index``: a version of ``Index`` highly optimized for 64-bit integer
  data, such as time stamps
* ``Float64Index``: a version of ``Index`` highly optimized for 64-bit float data
* ``MultiIndex``: the standard hierarchical index object
* ``DatetimeIndex``: An Index object with ``Timestamp`` boxed elements (impl are the int64 values)
* ``TimedeltaIndex``: An Index object with ``Timedelta`` boxed elements (impl are the in64 values)
* ``PeriodIndex``: An Index object with Period elements

There are functions that make the creation of a regular index easy:

* ``date_range``: fixed frequency date range generated from a time rule or
  DateOffset. An ndarray of Python datetime objects
* ``period_range``: fixed frequency date range generated from a time rule or
  DateOffset. An ndarray of ``Period`` objects, representing timespans

The motivation for having an ``Index`` class in the first place was to enable
different implementations of indexing. This means that it's possible for you,
the user, to implement a custom ``Index`` subclass that may be better suited to
a particular application than the ones provided in pandas.

From an internal implementation point of view, the relevant methods that an
``Index`` must define are one or more of the following (depending on how
incompatible the new object internals are with the ``Index`` functions):

* ``get_loc``: returns an "indexer" (an integer, or in some cases a
  slice object) for a label
* ``slice_locs``: returns the "range" to slice between two labels
* ``get_indexer``: Computes the indexing vector for reindexing / data
  alignment purposes. See the source / docstrings for more on this
* ``get_indexer_non_unique``: Computes the indexing vector for reindexing / data
  alignment purposes when the index is non-unique. See the source / docstrings
  for more on this
* ``reindex``: Does any pre-conversion of the input index then calls
  ``get_indexer``
* ``union``, ``intersection``: computes the union or intersection of two
  Index objects
* ``insert``: Inserts a new label into an Index, yielding a new object
* ``delete``: Delete a label, yielding a new object
* ``drop``: Deletes a set of labels
* ``take``: Analogous to ndarray.take

MultiIndex
~~~~~~~~~~

Internally, the ``MultiIndex`` consists of a few things: the **levels**, the
integer **codes** (until version 0.24 named *labels*), and the level **names**:

.. ipython:: python

   index = pd.MultiIndex.from_product([range(3), ['one', 'two']],
                                      names=['first', 'second'])
   index
   index.levels
   index.codes
   index.names

You can probably guess that the codes determine which unique element is
identified with that location at each layer of the index. It's important to
note that sortedness is determined **solely** from the integer codes and does
not check (or care) whether the levels themselves are sorted. Fortunately, the
constructors ``from_tuples`` and ``from_arrays`` ensure that this is true, but
if you compute the levels and codes yourself, please be careful.

Values
~~~~~~

Pandas extends NumPy's type system with custom types, like ``Categorical`` or
datetimes with a timezone, so we have multiple notions of "values". For 1-D
containers (``Index`` classes and ``Series``) we have the following convention:

* ``cls._ndarray_values`` is *always* a NumPy ``ndarray``. Ideally,
  ``_ndarray_values`` is cheap to compute. For example, for a ``Categorical``,
  this returns the codes, not the array of objects.
* ``cls._values`` refers is the "best possible" array. This could be an
  ``ndarray``, ``ExtensionArray``, or in ``Index`` subclass (note: we're in the
  process of removing the index subclasses here so that it's always an
  ``ndarray`` or ``ExtensionArray``).

So, for example, ``Series[category]._values`` is a ``Categorical``, while
``Series[category]._ndarray_values`` is the underlying codes.

.. _ref-subclassing-pandas:

Subclassing pandas data structures
----------------------------------

This section has been moved to :ref:`extending.subclassing-pandas`.
