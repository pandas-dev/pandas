.. currentmodule:: pandas
.. _merging:

***************************
Merging / Joining data sets
***************************


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


Joining / merging DataFrames
----------------------------

The **join** method provides effectively SQL-like semantic for
combining related data sets. The basic join consists of two DataFrame
arguments and

::

    >>> df1
               A              B
    2000-01-03 00:00:00    -0.1174        -0.941
    2000-01-04 00:00:00    -0.6034        -0.008094
    2000-01-05 00:00:00    -0.3816        -0.9338
    2000-01-06 00:00:00    -0.3298        -0.9548
    2000-01-07 00:00:00    0.9576         0.4652
    2000-01-10 00:00:00    -0.7208        -1.131
    2000-01-11 00:00:00    1.568          0.8498
    2000-01-12 00:00:00    0.3717         -0.2323
    2000-01-13 00:00:00    -1.428         -1.997
    2000-01-14 00:00:00    -1.084         -0.271

    >>> df2
               C              D
    2000-01-03 00:00:00    0.2833         -0.1937
    2000-01-05 00:00:00    1.868          1.207
    2000-01-07 00:00:00    -0.8586        -0.7367
    2000-01-11 00:00:00    2.121          0.9104
    2000-01-13 00:00:00    0.7856         0.9063


    df1.join(df2)
               A              B              C              D
    2000-01-03 00:00:00    -0.1174        -0.941         0.2833         -0.1937
    2000-01-04 00:00:00    -0.6034        -0.008094      NaN            NaN
    2000-01-05 00:00:00    -0.3816        -0.9338        1.868          1.207
    2000-01-06 00:00:00    -0.3298        -0.9548        NaN            NaN
    2000-01-07 00:00:00    0.9576         0.4652         -0.8586        -0.7367
    2000-01-10 00:00:00    -0.7208        -1.131         NaN            NaN
    2000-01-11 00:00:00    1.568          0.8498         2.121          0.9104
    2000-01-12 00:00:00    0.3717         -0.2323        NaN            NaN
    2000-01-13 00:00:00    -1.428         -1.997         0.7856         0.9063
    2000-01-14 00:00:00    -1.084         -0.271         NaN            NaN

::

    >>> df1.join(df2, how='inner')
               A              B              C              D
    2000-01-03 00:00:00    -0.1174        -0.941         0.2833         -0.1937
    2000-01-05 00:00:00    -0.3816        -0.9338        1.868          1.207
    2000-01-07 00:00:00    0.9576         0.4652         -0.8586        -0.7367
    2000-01-11 00:00:00    1.568          0.8498         2.121          0.9104
    2000-01-13 00:00:00    -1.428         -1.997         0.7856         0.9063

The index (row labels) are the default key for joining, but a column
can also be used for a similar SQL-like join: It is also frequently
necessary to join (or *merge*) data sets based on some other key
mapping.

::

    >>> df2
               C              D              key
    2000-01-03 00:00:00    0.2833         -0.1937        0
    2000-01-05 00:00:00    1.868          1.207          1
    2000-01-07 00:00:00    -0.8586        -0.7367        0
    2000-01-11 00:00:00    2.121          0.9104         1
    2000-01-13 00:00:00    0.7856         0.9063         0

    >>> df3
     code
    0    foo
    1    bar

    >>> df2.join(df3, on='key')
               C              D              code           key
    2000-01-03 00:00:00    0.2833         -0.1937        foo            0
    2000-01-05 00:00:00    1.868          1.207          bar            1
    2000-01-07 00:00:00    -0.8586        -0.7367        foo            0
    2000-01-11 00:00:00    2.121          0.9104         bar            1
    2000-01-13 00:00:00    0.7856         0.9063         foo            0
