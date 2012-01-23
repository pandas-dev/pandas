.. currentmodule:: pandas
.. _gotchas:

.. ipython:: python
   :suppress:

   import numpy as np
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

*******************
Caveats and Gotchas
*******************

``NaN``, Integer ``NA`` values and ``NA`` type promotions
---------------------------------------------------------

Choice of ``NA`` representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For lack of ``NA`` (missing) support from the ground up in NumPy and Python in
general, we were given the difficult choice between either

- A *masked array* solution: an array of data and an array of boolean values
  indicating whether a value
- Using a special sentinel value, bit pattern, or set of sentinel values to
  denote ``NA`` across the dtypes


Support for integer ``NA``
~~~~~~~~~~~~~~~~~~~~~~~~~~

``NA`` type promotions
~~~~~~~~~~~~~~~~~~~~~~

Integer indexing
----------------

Label-based slicing conventions
-------------------------------

Non-monotonic indexes require exact matches
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Endpoints are inclusive
~~~~~~~~~~~~~~~~~~~~~~~

Compared with standard Python sequence slicing in which the slice endpoint is
not inclusive, label-based slicing in pandas **is inclusive**. The primary
reason for this is that it is often not possible to easily the "successor" or
next element after a particular label in an index. For example, consider the
following Series:

.. ipython:: python

   s = Series(randn(6), index=list('abcdef'))
   s

Suppose we wished to slice from ``c`` to ``e``, using integers this would be

.. ipython:: python

   s[2:5]

However, if you only had ``c`` and ``e``, determining the next element in the
index can be somewhat complicated. For example, the following does not work:

::

    s.ix['c':'e'+1]

A very common use case is to limit a time series to start and end at two
specific dates. To enable this, we made the design design to make label-based slicing include both endpoints:

.. ipython:: python

    s.ix['c':'e']

This is most definitely a "practicality beats purity" sort of thing, but it is
something to watch out for is you expect label-based slicing to behave exactly
in the way that standard Python integer slicing works.

Miscellaneous indexing gotchas
------------------------------

Reindex versus ix gotchas
~~~~~~~~~~~~~~~~~~~~~~~~~

Many users will find themselves using the ``ix`` indexing capabilities as a
concise means of selecting data from a pandas object:

.. ipython:: python

   df = DataFrame(randn(6, 4), columns=['one', 'two', 'three', 'four'],
                  index=list('abcdef'))
   df
   df.ix[['b', 'c', 'e']]

This is, of course, completely equivalent *in this case* to using th
``reindex`` method:

.. ipython:: python

   df.reindex(['b', 'c', 'e'])

Some might conclude that ``ix`` and ``reindex`` are 100% equivalent based on
this. This is indeed true **except in the case of integer indexing**. For
example, the above operation could alternately have been expressed as:

.. ipython:: python

   df.ix[[1, 2, 4]]

If you pass ``[1, 2, 4]`` to ``reindex`` you will get another thing entirely:

.. ipython:: python

   df.reindex([1, 2, 4])

So it's important to remember that ``reindex`` is **strict label indexing
only**. This can lead to some potentially surprising results in pathological
cases where an index contains, say, both integers and strings:

.. ipython:: python

   s = Series([1, 2, 3], index=['a', 0, 1])
   s
   s.ix[[0, 1]]
   s.reindex([0, 1])

Because the index in this case does not contain solely integers, ``ix`` falls
back on integer indexing. By contrast, ``reindex`` only looks for the values
passed in the index, thus finding the integers ``0`` and ``1``. While it would
be possible to insert some logic to check whether a passed sequence is all
contained in the index, that logic would exact a very high cost in large data
sets.
