.. _series:

*******************
Series / TimeSeries
*******************

.. currentmodule:: pandas

:class:`Series` is a NumPy ndarray subclass which contains a vector
of labels corresponding to the data values. The labels, which will be
referred to everywhere as the **index**, can be any Python object;
common examples are datetimes or strings. The object is designed to
handle missing observations and support arithmetic between
different-sized Series by matching index values.

Because Series is an ndarray, it can be used interchangeably in
NumPy-based functions expecting one-dimensional ndarrays.

.. note::

    The basic concepts presented here apply to the higher dimensional
    data structures in pandas as well

Construction
------------

There are a number of ways to create Series objects. The most common
is to use the default constructor and pass two equal-length sequences:
one for the values, one for the index:

::

    >>> import numpy as np; from pandas import *

    >>> values = np.arange(5.)
    >>> labels = ['a', 'b', 'c', 'd', 'e']
    >>> s = Series(values, index=labels)
    >>> s
    a    0.0
    b    1.0
    c    2.0
    d    3.0
    e    4.0

We could also create this Series from a dict representing the data:

::

    >>> data = {'a': 0.0, 'b': 1.0, 'c': 2.0, 'd': 3.0, 'e': 4.0}
    >>> Series.fromDict(data)
    a    0.0
    b    1.0
    c    2.0
    d    3.0
    e    4.0


Any Series instance has the attribute **index** which is an Index
object containing the value labels:

::

    >>> s.index
    Index([a, b, c, d, e], dtype=object)

The index defines the *__contains__* behavior of the Series:

::

    >>> 'a' in s
    True

    >>> 'a' in s.index
    True

If an index contains all Python datetime objects, the created series
will be of type TimeSeries (so it is never necessary to explicitly use
the TimeSeries constructor):

::

    >>> dates
    [datetime.datetime(2009, 1, 1, 0, 0),
     datetime.datetime(2009, 1, 2, 0, 0),
     datetime.datetime(2009, 1, 5, 0, 0),
     datetime.datetime(2009, 1, 6, 0, 0),
     datetime.datetime(2009, 1, 7, 0, 0)]

    >>> ts = Series(values, index=dates)
    2009-01-01 00:00:00    0.0
    2009-01-02 00:00:00    1.0
    2009-01-05 00:00:00    2.0
    2009-01-06 00:00:00    3.0
    2009-01-07 00:00:00    4.0

    >>> type(ts)
    <class 'pandas.core.series.TimeSeries'>

Summary of constructors
~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Series
   Series.fromDict
   Series.fromValue

Indexing
--------

A Series operations similar to both an ndarray and a Python dict. For
example, values can be accessed either by integer location or by
index value:

::

    >>> s[1]
    1.0
    >>> s['b']
    1.0

If the index contains integers and there is ambiguity, the index will
be preferred.

For completeness of the dict-like interface, the **get** function is
provided for analogous behavior:

::

    >>> s.get('f')
    None

    >>> s.get('f', 0)
    0

Standard Python boolean indexing works as expected, as do slicing and
NumPy fancy indexing:

::

    >>> s[s > 2]
    d    3.0
    e    4.0

    >>> s[-3:]
    c    2.0
    d    3.0
    e    4.0

    >>> s[[4, 3, 1]]
    e    4.0
    d    3.0
    b    1.0

Observe that a new Index has been constructed corresponding to the
selected data.

Of course, the same behavior applies to *setting* values:

::

    >>> s[s > 2] = -1
    >>> print s
    a    0.0
    b    1.0
    c    2.0
    d    -1.0
    e    -1.0


Arithmetic, data alignment
--------------------------

Binary operations between Series objects ensure that two values being
combined have the same index value. This serves to prevent a lot of
the headache generally associated with heterogeneous data; the user is
**not** to ensure that all Series have the same index.

::

    >>> s + 2 * s
    a    0.0
    b    3.0
    c    6.0
    d    9.0
    e    12.0

    >>> s + s[2:]
    a    nan
    b    nan
    c    4.0
    d    6.0
    e    8.0

In this latter example, you can see that, since the **a** and **b**
values were missing in the second Series, the sum has NaN in those
locations. In general, pandas represents missing data as NaN (more on
this below).

Handling missing data and reindexing
------------------------------------

For all of the pandas data structures, we chose to represent missing
data as NaN. However, missing data could be represented in some other
forms (e.g. *None* values generated from DBNULL values in SQL
data). This problem is compounded by the fact that *numpy.isnan* is
only valid on float arrays. For this reason, pandas includes two
functions for testing for validity, **isnull** and **notnull**. These
functions are implemented in Cython and provide reasonably good
performance on object arrays. For numerical arrays, the performance
will be equivalent to *numpy.isfinite*.

::

    >>> s
    a    0.0
    b    1.0
    c    nan
    d    3.0
    e    4.0

    >>> isnull(s)
    a    False
    b    False
    c    True
    d    False
    e    False

    >>> isnull(None)
    True

These functions can be used, for example, to select only valid data
from a Series. Since this is such a common operation, a method
**valid** to do the same thing:

::

    >>> s[notnull(s)]
    a    0.0
    b    1.0
    d    3.0
    e    4.0

    >>> s.valid()
    a    0.0
    b    1.0
    d    3.0
    e    4.0

.. note::

    The choice of using NaN for missing data was one of practicality
    and ease-of-implementation. It differs from the MaskedArray
    approach of, for example, :mod:`scikits.timeseries`.

Filling, padding, and interpolating values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Series.fill
   Series.interpolate

Reindexing
~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Series.reindex
   Series.valid
   Series.merge
   Series.truncate

Iterating
---------

Series iterates by default over its values as though it were a regular
ndarray.

Otherwise, methods providing dict-like iteration are available:

.. autosummary::
   :toctree: generated/

   Series.keys
   Series.values
   Series.iteritems

Basic statistical functions
---------------------------

There are many built-in ndarray methods providing basic descriptive
statistics. Since these do not handle missing observations (which are
represented in our case as NaN), we've overrided these methods to do
the appropriate handling.

For example:

::

    >>> s
    a    0.0
    b    1.0
    c    nan
    d    3.0
    e    4.0

    >>> s.count()
    4

    >>> s.std()
    1.8257418583505536

    >>> s.cumsum()
    a    0.0
    b    1.0
    c    nan
    d    4.0
    e    8.0

Due to the way the numpy.{sum, mean, var, std} are implemented, they
can be used safely:

::

    >>> np.mean(s)
    2.0

Method summary
~~~~~~~~~~~~~~

.. autosummary::
   :toctree: generated/

   Series.count
   Series.sum
   Series.mean
   Series.min
   Series.max
   Series.std
   Series.var
   Series.skew
   Series.median
   Series.cumsum
   Series.cumprod

Additionally, some other useful methods not present in ndarray are
implemented:

.. autosummary::
   :toctree: generated/

   Series.corr
   Series.cap
   Series.floor

Sorting
-------

TODO

.. autosummary::
   :toctree: generated/

   Series.argsort
   Series.sort
   Series.order

TimeSeries-oriented methods
---------------------------

TODO

.. autosummary::
   :toctree: generated/

   Series.asfreq
   Series.shift
   Series.asOf
   Series.weekday

GroupBy functionality
---------------------

TODO

Plotting
--------

TODO

.. autosummary::
   :toctree: generated/

   Series.plot

Misc methods
------------

TODO

.. autosummary::
   :toctree: generated/

   Series.append
   Series.combineFunc
   Series.combineFirst
   Series.map
   Series.copy
   Series.toCSV
   Series.diff
