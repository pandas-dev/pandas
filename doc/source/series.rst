.. currentmodule:: pandas
.. _series:

.. _series.arithmetic:

Accessing underlying data
-------------------------

The data stored in a Series can be accessed via the **values**
attribute (which is a property returning a view of the Series as an
ndarray). The higher dimensional pandas data structures observe the
same interface for accessing the underlying data.

Handling missing data and reindexing
------------------------------------

For all of the pandas data structures, we chose to represent missing
data as NaN. However, missing data could be represented in some other
forms (e.g. *None* values generated from null values in SQL
data). This problem is compounded by the fact that *numpy.isnan* is
only valid on float arrays. For this reason, pandas includes two
functions for testing validity, **isnull** and **notnull**. These
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

.. _series.reindexing:

Reindexing
~~~~~~~~~~

Reindexing is one of the most important features of the Series and the
other pandas data structures. In essence it means: *conform data to a
specified index*.

Using our prior example, we can see the basic behavior:

::

    >>> s.reindex(['f', 'a', 'd', 'e'])
    f    nan
    a    0.0
    d    3.0
    e    4.0

As you can see, the new index order is as inputted, and values not
present in the Series appear as NaN.

For TimeSeries or other ordered Series, an additional argument can be
specified to perform forward- (referred to as "padding") or
back-filling:

::

    >>> ts
    2009-01-02 00:00:00    1.0
    2009-01-07 00:00:00    4.0

    >>> ts.reindex(dates, method='pad')
    2009-01-01 00:00:00    nan
    2009-01-02 00:00:00    1.0
    2009-01-05 00:00:00    1.0
    2009-01-06 00:00:00    1.0
    2009-01-07 00:00:00    4.0

    >>> ts.reindex(dates, method='backfill')
    2009-01-01 00:00:00    1.0
    2009-01-02 00:00:00    1.0
    2009-01-05 00:00:00    4.0
    2009-01-06 00:00:00    4.0
    2009-01-07 00:00:00    4.0

.. note::

    This filling logic assumes that the both the new and old Index
    objects have ordered values.

Two common reindexing methods are provided: **valid** (which we
already mentioned) and **truncate** (for selecting intervals of index
values).

::

    >>> ts
    2009-01-01 00:00:00    0.0
    2009-01-02 00:00:00    1.0
    2009-01-05 00:00:00    2.0
    2009-01-06 00:00:00    3.0
    2009-01-07 00:00:00    4.0

    >>> ts.truncate(before=datetime(2009, 1, 5), after=datetime(2009, 1, 6))
    2009-01-05 00:00:00    2.0
    2009-01-06 00:00:00    3.0

Since writing out datetimes interactively like that can be a bit
verbose, one can also pass a string date representation:

::

    >>> ts.truncate(after='1/5/2009')
    2009-01-01 00:00:00    0.0
    2009-01-02 00:00:00    1.0
    2009-01-05 00:00:00    2.0


.. autosummary::
   :toctree: generated/

   Series.reindex
   Series.valid
   Series.truncate

Filling, padding, and interpolating values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is often desirable to deal with missing values in some specific
way, especially for time series data. As seen above, reindexing can be
a useful approach, however we frequently will want to deal
specifically with missing data in a particular way.

The **fill** method provides two distinct behaviors: filling in NaNs
with a static value or alternately padding / backfilling as with
**reindex**:

::

    >>> ts
    2009-01-01 00:00:00	0.0
    2009-01-02 00:00:00	1.0
    2009-01-05 00:00:00	nan
    2009-01-06 00:00:00	nan
    2009-01-07 00:00:00	4.0

    >>> ts.fill(value=6)
    2009-01-01 00:00:00	0.0
    2009-01-02 00:00:00	1.0
    2009-01-05 00:00:00	6.0
    2009-01-06 00:00:00	6.0
    2009-01-07 00:00:00	4.0

    >>> ts.fill(method='pad')
    2009-01-01 00:00:00	0.0
    2009-01-02 00:00:00	1.0
    2009-01-05 00:00:00	1.0
    2009-01-06 00:00:00	1.0
    2009-01-07 00:00:00	4.0


In a similar vein, values can be linearly interpolated in either a
naive way or in a time-spaced way

::

    >>> ts.interpolate()
    2009-01-01 00:00:00	0.0
    2009-01-02 00:00:00	1.0
    2009-01-05 00:00:00	2.0
    2009-01-06 00:00:00	3.0
    2009-01-07 00:00:00	4.0

    >>> ts.interpolate(method='time')
    2009-01-01 00:00:00	0.0
    2009-01-02 00:00:00	1.0
    2009-01-05 00:00:00	2.8
    2009-01-06 00:00:00	3.4
    2009-01-07 00:00:00	4.0

.. autosummary::
   :toctree: generated/

   Series.reindex
   Series.fill
   Series.interpolate

Relabeling (renaming) Series index
----------------------------------

One might want to relabel the index based on a mapping or
function. For this the :func:`Series.rename` method is provided. It
can accept either a dict or a function:

::

    >>> s
	a    -0.544970223484
	b    -0.946388873158
	c    0.0360854957476
	d    -0.795018577574
	e    0.195977583894

	>>> s.rename(str.upper)
	A    -0.544970223484
	B    -0.946388873158
	C    0.0360854957476
	D    -0.795018577574
	E    0.195977583894

Iterating
---------

Series iterates by default over its values as though it were a regular
ndarray.

Otherwise, methods providing dict-like iteration are available:

::

    >>> for x in ts:
            print x
    0.0
    1.0
    2.0
    3.0
    4.0

    >>> for index, value in ts.iteritems():
            print index, value
    2009-01-01 00:00:00 0.0
    2009-01-02 00:00:00 1.0
    2009-01-05 00:00:00 2.0
    2009-01-06 00:00:00 3.0
    2009-01-07 00:00:00 4.0


.. autosummary::
   :toctree: generated/

   Series.values
   Series.iteritems

.. _series.statistics:

Basic statistical functions
---------------------------

There are many built-in ndarray methods providing basic descriptive
statistics. Since these do not handle missing observations (which are
represented in our case as NaN), we've overridden these methods to do
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

Merging Series based on key
---------------------------

You may be occasionally interested in joining data sets which are
keyed on different index values. This comes down to a simple mapping
problem in the one dimensional case and will be more interesting in
the 2- and 3-D cases, but the basic concept is the same:

::

    >>> s = Series(['six', 'seven', 'six', 'seven', 'six'],
                   index=['a', 'b', 'c', 'd', 'e'])
    >>> t = Series({'six' : 6., 'seven' : 7.})

    >>> s
    a	six
    b	seven
    c	six
    d	seven
    e	six

    >>> s.merge(t)
    a	6.0
    b	7.0
    c	6.0
    d	7.0
    e	6.0



.. autosummary::
   :toctree: generated/

   Series.merge

Sorting
-------

A number of methods for sorting Series data are provided:

::

    >>> s = Series(randn(5), index=['a', 'b', 'c', 'd', 'e'])
    >>> s
    a    -0.308339649397
    b    -0.447658314192
    c    -0.391847354208
    d    0.427084101354
    e    1.51816072219

    >>> s.order()
    b    -0.447658314192
    c    -0.391847354208
    a    -0.308339649397
    d    0.427084101354
    e    1.51816072219

    >>> s.argsort()
    a    1
    b    2
    c    0
    d    3
    e    4

    >>> s.sort()    # in-place sort
    >>> s
    b    -0.447658314192
    c    -0.391847354208
    a    -0.308339649397
    d    0.427084101354
    e    1.51816072219

:func:`Series.order` is intended to behave similarly to the R function
of the same name. In the presence of missing data it accepts an
optional argument specifying where to sort the NaN values (either the
end or the beginning). The default is to sort them to the end, which
is the new sorting behavior in NumPy >= 1.4.0:

::

    >>> s
    a    -2.21668112685
    b    -0.520791835078
    c    NaN
    d    -0.788775281233
    e    -0.779555719818

    >>> s.order()
    a    -2.21668112685
    d    -0.788775281233
    e    -0.779555719818
    b    -0.520791835078
    c    NaN

    >>> s.order(missingAtEnd=False)
    c    NaN
    a    -2.21668112685
    d    -0.788775281233
    e    -0.779555719818
    b    -0.520791835078


.. autosummary::
   :toctree: generated/

   Series.argsort
   Series.sort
   Series.order


TimeSeries-oriented methods
---------------------------

.. seealso::
    :ref:`Reindexing methods <series.reindexing>`;
    :ref:`DateRange and date offsets / time rules <datetools>`

.. note::

    While pandas does not force you to sort your dates, many of these
    methods may have unexpected or incorrect behavior in that case. In
    other words, *be careful*.

When working with time series data, a number of different
time-oriented operations may be useful. The first is **frequency
conversion**, which has similar options to :func:`Series.reindex`:

::

    >>> dr = DateRange('1/1/2010', periods=10,
                       offset=datetools.BMonthEnd())
    >>> ts = Series(np.arange(10.), index=dr)
    >>> ts
    2010-01-29 00:00:00    0.0
    2010-02-26 00:00:00    1.0
    2010-03-31 00:00:00    2.0
    2010-04-30 00:00:00    3.0
    2010-05-31 00:00:00    4.0
    2010-06-30 00:00:00    5.0
    2010-07-30 00:00:00    6.0
    2010-08-31 00:00:00    7.0
    2010-09-30 00:00:00    8.0
    2010-10-29 00:00:00    9.0

    >>> ts.asfreq('WEEKDAY', method='pad')
    2010-01-29 00:00:00    0.0
    2010-02-01 00:00:00    0.0
    2010-02-02 00:00:00    0.0
    2010-02-03 00:00:00    0.0
    2010-02-04 00:00:00    0.0
    <snip>
    2010-10-22 00:00:00    8.0
    2010-10-25 00:00:00    8.0
    2010-10-26 00:00:00    8.0
    2010-10-27 00:00:00    8.0
    2010-10-28 00:00:00    8.0
    2010-10-29 00:00:00    9.0

We often will also want to **shift** or *lag* a TimeSeries:

::

    >>> ts.shift(1)
    2010-01-29 00:00:00    NaN
    2010-02-26 00:00:00    0.0
    2010-03-31 00:00:00    1.0
    2010-04-30 00:00:00    2.0
    2010-05-31 00:00:00    3.0
    2010-06-30 00:00:00    4.0
    2010-07-30 00:00:00    5.0
    2010-08-31 00:00:00    6.0
    2010-09-30 00:00:00    7.0
    2010-10-29 00:00:00    8.0

    >>> ts.shift(5, offset=datetools.bday)
    2010-02-05 00:00:00    0.0
    2010-03-05 00:00:00    1.0
    2010-04-07 00:00:00    2.0
    2010-05-07 00:00:00    3.0
    2010-06-07 00:00:00    4.0
    2010-07-07 00:00:00    5.0
    2010-08-06 00:00:00    6.0
    2010-09-07 00:00:00    7.0
    2010-10-07 00:00:00    8.0
    2010-11-05 00:00:00    9.0

In the presence of missing data with sorted dates

A convenience method for selecting weekdays, similar to
:mod:`scikits.timeseries` is also provided:

::

    >>> dr = DateRange('1/1/2010', periods=10, offset=datetools.bday)
    >>> ts = Series(np.arange(10.), index=dr)
    >>> ts[ts.weekday == 2]
    2010-01-06 00:00:00    3.0
    2010-01-13 00:00:00    8.0
