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

* :class:`Index`: the generic "ordered set" object, an ndarray of object dtype
  assuming nothing about its contents. The labels must be hashable (and
  likely immutable) and unique. Populates a dict of label to location in
  Cython to do ``O(1)`` lookups.
* :class:`MultiIndex`: the standard hierarchical index object
* :class:`DatetimeIndex`: An Index object with :class:`Timestamp` boxed elements (impl are the int64 values)
* :class:`TimedeltaIndex`: An Index object with :class:`Timedelta` boxed elements (impl are the in64 values)
* :class:`PeriodIndex`: An Index object with Period elements

There are functions that make the creation of a regular index easy:

* :func:`date_range`: fixed frequency date range generated from a time rule or
  DateOffset. An ndarray of Python datetime objects
* :func:`period_range`: fixed frequency date range generated from a time rule or
  DateOffset. An ndarray of :class:`Period` objects, representing timespans

.. warning::

   Custom :class:`Index` subclasses are not supported, custom behavior should be implemented using the :class:`ExtensionArray` interface instead.

MultiIndex
~~~~~~~~~~

Internally, the :class:`MultiIndex` consists of a few things: the **levels**, the
integer **codes**, and the level **names**:

.. ipython:: python

   index = pd.MultiIndex.from_product(
       [range(3), ["one", "two"]], names=["first", "second"]
   )
   index
   index.levels
   index.codes
   index.names

You can probably guess that the codes determine which unique element is
identified with that location at each layer of the index. It's important to
note that sortedness is determined **solely** from the integer codes and does
not check (or care) whether the levels themselves are sorted. Fortunately, the
constructors :meth:`~MultiIndex.from_tuples` and :meth:`~MultiIndex.from_arrays` ensure
that this is true, but if you compute the levels and codes yourself, please be careful.

Values
~~~~~~

pandas extends NumPy's type system with custom types, like :class:`Categorical` or
datetimes with a timezone, so we have multiple notions of "values". For 1-D
containers (``Index`` classes and ``Series``) we have the following convention:

* ``cls._values`` refers is the "best possible" array. This could be an
  ``ndarray`` or ``ExtensionArray``.

So, for example, ``Series[category]._values`` is a ``Categorical``.

.. _ref-subclassing-pandas:

Subclassing pandas data structures
----------------------------------

This section has been moved to :ref:`extending.subclassing-pandas`.
