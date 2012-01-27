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

For many reasons we chose the latter. After years of production use it has
proven, at least in my opinion, to be the best decision given the state of
affairs in NumPy and Python in general. The special value ``NaN``
(Not-A-Number) is used everywhere as the ``NA`` value, and there are API
functions ``isnull`` and ``notnull`` which can be used across the dtypes to
detect NA values.

However, it comes with it a couple of trade-offs which I most certainly have
not ignored.

Support for integer ``NA``
~~~~~~~~~~~~~~~~~~~~~~~~~~

In the absence of high performance ``NA`` support being built into NumPy from
the ground up, the primary casualty is the ability to represent NAs in integer
arrays. For example:

.. ipython:: python

   s = Series([1, 2, 3, 4, 5], index=list('abcde'))
   s
   s.dtype

   s2 = s.reindex(['a', 'b', 'c', 'f', 'u'])
   s2
   s2.dtype

This trade-off is made largely for memory and performance reasons, and also so
that the resulting Series continues to be "numeric". One possibility is to use
``dtype=object`` arrays instead.

``NA`` type promotions
~~~~~~~~~~~~~~~~~~~~~~

When introducing NAs into an existing Series or DataFrame via ``reindex`` or
some other means, boolean and integer types will be promoted to a different
dtype in order to store the NAs. These are summarized by this table:

.. csv-table::
   :header: "Typeclass","Promotion dtype for storing NAs"
   :widths: 40,60

   ``floating``, no change
   ``object``, no change
   ``integer``, cast to ``float64``
   ``boolean``, cast to ``object``

While this may seem like a heavy trade-off, in practice I have found very few
cases where this is an issue in practice. Some explanation for the motivation
here in the next section.

Why not make NumPy like R?
~~~~~~~~~~~~~~~~~~~~~~~~~~

Many people have suggested that NumPy should simply emulate the ``NA`` support
present in the more domain-specific statistical programming langauge `R
<http://r-project.org>`__. Part of the reason is the NumPy type hierarchy:

.. csv-table::
   :header: "Typeclass","Dtypes"
   :widths: 30,70
   :delim: |

   ``numpy.floating`` | ``float16, float32, float64, float128``
   ``numpy.integer`` | ``int8, int16, int32, int64``
   ``numpy.unsignedinteger`` | ``uint8, uint16, uint32, uint64``
   ``numpy.object_`` | ``object_``
   ``numpy.bool_`` | ``bool_``
   ``numpy.character`` | ``string_, unicode_``

The R language, by contrast, only has a handful of built-in data types:
``integer``, ``numeric`` (floating-point), ``character``, and
``boolean``. ``NA`` types are implemented by reserving special bit patterns for
each type to be used as the missing value. While doing this with the full NumPy
type hierarchy would be possible, it would be a more substantial trade-off
(especially for the 8- and 16-bit data types) and implementation undertaking.

An alternate approach is that of using masked arrays. A masked array is an
array of data with an associated boolean *mask* denoting whether each value
should be considered ``NA`` or not. I am personally not in love with this
approach as I feel that overall it places a fairly heavy burden on the user and
the library implementer. Additionally, it exacts a fairly high performance cost
when working with numerical data compared with the simple approach of using
``NaN``. Thus, I have chosen the Pythonic "practicality beats purity" approach
and traded integer ``NA`` capability for a much simpler approach of using a
special value in float and object arrays to denote ``NA``, and promoting
integer arrays to floating when NAs must be introduced.

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
specific dates. To enable this, we made the design design to make label-based
slicing include both endpoints:

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
