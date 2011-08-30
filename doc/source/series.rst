.. currentmodule:: pandas
.. _series:

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
