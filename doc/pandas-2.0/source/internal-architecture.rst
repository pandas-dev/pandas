.. _internal-architecture:

.. ipython:: python
   :suppress:

   import numpy as np
   import pandas as pd
   np.set_printoptions(precision=4, suppress=True)
   pd.options.display.max_rows = 100

===============================
 Internal Architecture Changes
===============================

Logical types and Physical Storage Decoupling
=============================================

Removal of BlockManager / new DataFrame internals
=================================================

``pandas.Array`` and ``pandas.Table``
=====================================

Missing data consistency
========================

Once the physical memory representation has been effectively decoupled from the
user API, we can consider various approaches to implementing missing data in a
consistent way for every logical pandas data type.

To motivate this, let's look at some integer data:

.. ipython:: python

   s = pd.Series([1, 2, 3, 4, 5])
   s
   s.dtype
   s.values

If we assign a ``numpy.NaN``, see what happens:

.. ipython:: python

   s[2] = np.NaN
   s
   s.dtype
   s.values

The story for boolean data is similar:

.. ipython:: python

   s = pd.Series([True, False, True])
   s.dtype
   s[2] = np.NaN
   s.dtype
   s.values

This implicit behavior appears in many scenarios, such as:

* Loading data from any source: databases, CSV files, R data files, etc.
* Joins or reindexing operations introducing missing data
* Pivot / reshape operations
* Time series resampling
* Certain types of GroupBy operations

A proposed solution
~~~~~~~~~~~~~~~~~~~

My proposal for introducing missing data into any NumPy type outside of
floating point (which uses ``NaN`` for now) and Python object (which uses
``None`` or ``NaN`` interchangeably) is to **allocate and manage an internal
bitmap** (which the user never sees). This has numerous benefits:

* 1 byte of memory overhead for each 8 values
* Bitmaps can propagate their nulls in C through bitwise ``&`` or ``|``
  operations, which are inexpensive.
* Getting and setting bits on modern hardware is very CPU-inexpensive. For
  single-pass array operations (like groupbys) on very large arrays this may
  also result in better CPU cache utilization (fewer main-memory reads of the
  bitmap).
* Hardware and SIMD "popcount" intrinsics (which can operate on 64-128 bits at
  a time) can be used to count bits and skip null-handling on segments of data
  containing no nulls.

Notably, this is the way that PostgreSQL handles null values. For example, we
might have:

.. code-block::

   [0, 1, 2, NA, NA, 5, 6, NA]

        i: 7 6 5 4 3 2 1 0
   bitmap: 0 1 1 0 0 1 1 1

Here, the convention of 1 for "not null" (a la PostgreSQL) and
least-significant bit ordering (LSB "bit endianness") is being used.

Under the new regime, users could simply write:

.. code-block:: python

   s[2] = pandas.NA

and the data type would be unmodified. It may be necessary to write something
akin to:

.. code-block:: python

   s.to_numpy(dtype=np.float64, na_rep=np.nan)

and that would emulate the current behavior. Attempts to use ``__array__` (for
example: calling ``np.sqrt`` on the data) would result in an error since we
will likely want to refuse to make a guess as for what casting behavior the
user desires.

Tradeoffs
~~~~~~~~~

One potential downside of the bitmap approach is that missing data implemented
outside of NumPy's domain will need to be explicitly converted if it is needed
in another library that only knows about NumPy. I argue that this is better
than the current

Proper types for strings and some non-numeric data
==================================================

I believe that frequently-occurring data types, such as UTF8 strings, are
important enough to deserve a dedicated logical pandas data type. This will
enable us both to enforce tighter API semantics (i.e. attempts to assign a
non-string into string data will be a ``TypeError``) and improved performance
and memory use under the hood. I will devote an entire section to talking about
strings.

C++11/14 for lowest implementation tier
=======================================

3rd-party native API (i.e. Cython and C / C++)
==============================================
