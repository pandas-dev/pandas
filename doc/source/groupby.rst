.. currentmodule:: pandas
.. _groupby:

==================
GroupBy operations
==================

When we talk about *group by* operations, we are referring to performing
operations on subsets of a data structure determined by some group membership
criteria. For example, for a time series we might wish to group data by month or
year. Generally, there are two primary kinds of operations of interest:

* **Aggregation**: computing a single number from a group of data. Some
    examples:

 * Counting group sizes (using the Python :func:`len` function)
 * Computing group sums or means

* **Transformation**: computing new values for a group, leaving the data
    structure size unchanged. For example

 * Rescale a group by its maximum value
 * Standardize a group (subtract the mean, divide by the standard deviation)

In designing :mod:`pandas` we wished to provide a simple method for carrying out
the above kinds of operations. Sometimes procedures don't fit into the above two
categories of actions, so we provide additionally a simple way to iterate
through groups, which we illustrate below.

Grouping by dicts and functions
-------------------------------

As the axes of :mod:`pandas` data structure are uniquely labeled, the most
natural way to group data is by providing some assignment of index labels to
groups. To motivate, consider a simple :class:`Series` instance:

::

    >>> anim = Series(np.random.randn(6),
                      index=['dog', 'cat', 'lion', 'tiger', 'duck', 'buffalo'])
    >>> anim
    dog        -0.319350343654
    cat        -0.224467238308
    lion       -0.109077865448
    tiger      0.302085950154
    duck       0.573503465533
    buffalo    3.13284403296


Now, suppose we wanted to sort these data by "big" and "small" animals. We could
store this information in a dict

::

    >>> sizes = {'dog' : 'small', 'cat' : 'small',
    ...          'lion' : 'big', 'tiger' : 'big',
    ...          'duck' : 'small', 'buffalo' : 'big'}

Now, to group **anim**, we call **groupby**:

::

    >>> grouped = anim.groupby(sizes)
    >>> grouped
    <pandas.core.groupby.SeriesGroupBy object at 0x450f1b0>

    >>> grouped.groups
    <class 'pandas.core.groupby.GroupDict'>
     big   -> 3 values
     small -> 3 values

    >>> grouped.groups['big']
    ['lion', 'tiger', 'buffalo']

Using a dict is not the only way to compute group membership. The most general
way is to pass a *function* which will be **called on each axis label**. For
example, suppose we wanted to group the **anim** :class:`Series` by the length
of their names. We could then do:

::

    >>> anim.groupby(len).groups
    <class 'pandas.core.groupby.GroupDict'>
     3 -> 2 values
     4 -> 2 values
     5 -> 1 values
     7 -> 1 values

     >>> anim.groupby(len).groups[3]
     ['dog', 'cat']

As a perhaps more practical example, consider a :class:`TimeSeries`:

::

    >>> ts
    2000-01-03 00:00:00    0.56602563035
    2000-01-04 00:00:00    -1.18181160355
    2000-01-05 00:00:00    -1.24377783306
    2000-01-06 00:00:00    -1.12168338772
    2000-01-07 00:00:00    -0.296293517834
    2000-01-10 00:00:00    1.02113963739
    2000-01-11 00:00:00    0.249620178862
    2000-01-12 00:00:00    1.4956464977
    2000-01-13 00:00:00    -0.636055727754
    2000-01-14 00:00:00    1.54251281494
    2000-01-17 00:00:00    -1.39778010066
    2000-01-18 00:00:00    -1.27507586526
    2000-01-19 00:00:00    -0.729895375953
    2000-01-20 00:00:00    0.23422793374
    2000-01-21 00:00:00    -0.98459653972
    2000-01-24 00:00:00    -0.953545408705
    2000-01-25 00:00:00    0.488341257671
    2000-01-26 00:00:00    -1.06129410159
    2000-01-27 00:00:00    -0.475895459693
    2000-01-28 00:00:00    0.0227612634861

Suppose we wished to group the data by day of week. Since the index values are :class:`datetime` objects, we can do:

::

    >>> by_weekday = ts.groupby(lambda d: d.weekday())
    >>> by_weekday.groups
    <class 'pandas.core.groupby.GroupDict'>
     0 -> 4 values
     1 -> 4 values
     2 -> 4 values
     3 -> 4 values
     4 -> 4 values

Understanding the internals of :class:`GroupBy` objects is not essential;
understanding how to group the data, however, is.

Higher dimensional data structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The **groupby** methods for :class:`DataFrame`, :class:`WidePanel`, and
:class:`LongPanel` additionally have an **axis** argument which allows you to
select the axis to group on:

::

    >>> df
                           A              B              C              D
    2000-01-03 00:00:00    -0.5084        0.3694         -0.7599        -0.4266
    2000-01-04 00:00:00    0.1099         1.008          0.9618         0.9122
    2000-01-05 00:00:00    -0.885         0.1133         1.986          1.528
    2000-01-06 00:00:00    -0.03933       -0.2023        1.134          1.403
    2000-01-07 00:00:00    -0.9176        -0.4374        0.1027         -0.2443
    2000-01-10 00:00:00    -0.5625        -0.9719        0.6453         -1.826
    2000-01-11 00:00:00    -0.9899        0.3301         2.305          -2.087
    2000-01-12 00:00:00    0.8925         -1.198         -1.216         0.9092
    2000-01-13 00:00:00    0.3172         -0.8189        0.1865         -2.708
    2000-01-14 00:00:00    -0.1046        0.2199         0.06312        -2.157
    2000-01-17 00:00:00    0.2557         0.7938         2.275          0.8623
    2000-01-18 00:00:00    -0.6812        -0.7695        0.05992        0.1709
    2000-01-19 00:00:00      0.4017         1.032          -0.3528        -1.319
    2000-01-20 00:00:00    0.1693         -0.4478        1.238          -0.4015
    2000-01-21 00:00:00    -0.7986        -2.472         -0.1168        -0.1419

    >>> df.groupby({'A' : 0, 'B' : 0, 'C' : 1, 'D' : 1}, axis=1)
    >>> df.groupby(lambda d: d.weekday(), axis=0)

The axis names and numbers for **groupby** are as with **reindex** and other related
methods:

* **DataFrame**: *index* (0) and *columns* (1)
* **WidePanel**: *index* (0), *major* (1), and *minor* (2)
* **LongPanel**: **Not yet implemented**

Grouping by a `DataFrame` column
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes the grouping information is contained in a :class:`DataFrame` of
interest:

::

    >>> df
            country        value
    foo     US             1.715
    bar     US             -0.6018
    baz     JP             0.09215
    quux    JP             -0.9134

One can group by a column of the `DataFrame` and perform operations:

::

    >>> for name, group in df.groupby('country'):
    ...     print name, group['value'].mean()

    JP -0.410642609813
    US 0.556770179169

See more detail on iterating over a `GroupBy` object below.

Aggregation
-----------

Using a :class:`GroupBy` object, one can call its **aggregate** function with a
function producing a single value from an :class:`ndarray`.

::

    >>> anim.groupby(sizes).aggregate(len)
    big      3
    small    3

    >>> anim.groupby(sizes).aggregate(np.mean)
    big      -0.0740475497108
    small    -0.840956376614


Or, as in the :class:`DataFrame` above (note that **agg** is an alias for
**aggregate**):

::

    >>> df.groupby(lambda d: d.weekday()).agg(np.mean)
         A              B              C              D
    0    -0.2717        0.06377        0.7201         -0.4636
    1    -0.5204        0.1896         1.109          -0.3347
    2    0.1364         -0.01753       0.1391         0.3727
    3    0.149          -0.4897        0.853          -0.5686
    4    -0.6069        -0.8965        0.01633        -0.8476

Transformation
--------------



Iterating through groups
------------------------
