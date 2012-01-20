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
