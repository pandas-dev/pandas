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

Since this is the most important, but perhaps also most controversial, change
(in my opinion) to pandas, I'm going to go over it in great detail. I think the
hardest part of coming up with clear language and definitions for concepts so
that we can communicate effectively. For example the term "data type" is vague
and may mean different things to different people.

A motivating example
~~~~~~~~~~~~~~~~~~~~

Before digging too much into the technical details and problems/solutions,
let's look at some code examples. It is not unusual to find code like this in
pandas's internals:

.. code-block:: python

    def create_from_value(value, index, dtype):
        # return a new empty value suitable for the dtype

        if is_datetimetz(dtype):
            subarr = DatetimeIndex([value] * len(index), dtype=dtype)
        elif is_categorical_dtype(dtype):
            subarr = Categorical([value] * len(index))
        else:
            if not isinstance(dtype, (np.dtype, type(np.dtype))):
                dtype = dtype.dtype
            subarr = np.empty(len(index), dtype=dtype)
            subarr.fill(value)

or

.. code-block:: python

   if is_categorical_dtype(dtype):
       upcast_cls = 'category'
   elif is_datetimetz(dtype):
       upcast_cls = 'datetimetz'
   elif issubclass(dtype.type, np.bool_):
       upcast_cls = 'bool'
   elif issubclass(dtype.type, np.object_):
       upcast_cls = 'object'
   elif is_datetime64_dtype(dtype):
       upcast_cls = 'datetime'
   elif is_timedelta64_dtype(dtype):
       upcast_cls = 'timedelta'
   else:
       upcast_cls = 'float'

I've cherry-picked one of a number of places where this type of datatype-based
branching happens.

The primary reason for this complexity is that pandas is using both NumPy's
dtype objects (which describe *physical storage*) as well as its own custom
data type objects as a proxy for pandas's *semantic logical types*.

Let's step back for a second and come up with clear language to steer the
discussion.

Some definitions
~~~~~~~~~~~~~~~~

Here is my attempt at definitions of some of the key terms:

* **Metadata**: data that describes other data (such as its in-memory layout)

* **Semantics**: The meaning / abstract interpretation of something. We often
  discuss the semantics (meaning) of computer programs (i.e. what they do,
  fundamentally) without touching upon low level details like machine
  representation, programming languages, compilers, operating systems, etc.

* **Physical data (or storage) types**: these are metadata objects which
  provide a description of the precise structure of a piece of data in memory.

  * In NumPy, the ``numpy.dtype`` object (aka ``PyArray_Descr`` in the C API)
    is metadata describing a single cell / value in an array. Combined with the
    ``shape`` and ``strides`` attributes of the ``ndarray`` object, you have
    enough information to perform O(1) random access on any cell in an
    ``ndarray`` and to assign these values to a C type (or, in the case, of
    structured dtypes, assign to a packed C struct).

  * This may or may not include a physical representation of NULL or missing
    data (for example: nullable float64 might be a physical type indicating a
    normal float64 array along with a bitmap of null/not-null indicators).

* **Logical data type**: metadata which describes the semantic content of a
  single value in an array or other collection of values. Depending on the
  logical type, it may map 1-to-1 to a physical type or not at all. Here are
  some examples:

  * The ``double`` or ``float64`` type may be viewed both as a logical type as
    well as a physical type (a 1-to-1 correspondence).

  * pandas's ``category`` dtype contains its own auxiliary array of category
    values (for example, the distinct strings collected from a string
    array). Based on the number of categories, the category ``codes`` (which
    reference the categories array) are stored in the smallest possible integer
    physical type (from ``int8`` to ``int64``, depending whether the data type
    can accommodate the codes). For example, if there are 50 codes, the data is
    represented in ``int8`` storage. For 1000 codes, it would be ``int16``.

  * Another example: timestamps may be physically stored in ``int64``
    storage, and these values are interpreted in the context of a particular
    time unit or resolution (e.g. nanoseconds, milliseconds, seconds).

In general, new logical types may be formed either by placing new semantics on
top of a single physical data type or some composition of physical or logical
types. For example: you could have a categorical type (a logical construct
consisting of multiple arrays of data) whose categories are some other logical
type.

For historical reasons, **pandas never developed a clear semantic separation in
its user API between logical and physical data types**. Also, the addition of
new, pandas-only "synthetic" dtypes that are unknown to NumPy (like
categorical, datetimetz, etc.) has expanded this conflation considerably. If
you also consider pandas's custom missing / NULL data behavior, the addition of
ad hoc missing data semantics to a physical NumPy data type created, by the
definitions above, a logical data type (call it ``object[nullable]`` for an
object array) without ever explicitly saying so.

You might be thinking, "Good job, Wes. You really messed that up!" I'd be
inclined to agree with you now in retrospect, but back in 2011 pandas was not
the super popular project that it is today, and we were truly riding on NumPy's
coat tails. So the extent to which NumPy concepts and APIs were used explicitly
in pandas made the library easier to adopt. Now in 2016, this feels
anachronistic / outdated.

High-level logical type proposal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As we have been discussing periodically on the pandas-dev mailing list and
GitHub, I am proposing that we start to unravel our current mess by defining
pandas-specific metadata objects that model the current semantics / behavior of
the project. What does this mean, exactly?

* Each NumPy dtype object will map 1-to-1 to an equivalent ``pandas.DataType``
  object.
* Existing pandas "extension dtypes" (like ``CategoricalDtype`` and
  ``DatetimeTZDtype``), which have been designed to mimic ``numpy.dtype``, will
  become logical type subclasses of ``pandas.DataType`` like every other type
  in pandas.

Since pandas is about assisting with data manipulation and analysis, at some
point you must invoke functions that are specialized to the specific physical
memory representation of your data. For example, pandas has its own
implementations of ``ndarray.take`` that are used internally for arrays of
positive integers that may contain NULL / NA values (which are represented as
-1 -- search the codebase for implementations of ``take_1d``).

The major goals of introducing a logical type abstraction are the follows:

* Simplifying "dynamic dispatch": invoking the right functions or choosing the
  right code branches based on the data type.
* Enabling pandas to decouple both its internal semantics and physical storage
  from NumPy's metadata and APIs. Note that this is already happening with
  categorical types, since a particular instance ``CategoricalDtype`` may
  physically be stored in one of 4 NumPy data types.

Physical storage decoupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~

By separating pandas data from the presumption of using a particular physical
``numpy.dtype`` internally, we can:

* Begin to better protect users from NumPy data semantics (which are frequently
  different from pandas's!) leaking through to the pandas user API. This can
  enable us to address long-standing inconsistencies or "rough edges" in pandas
  that have persisted due to our tight semantic coupling to NumPy.

* We can consider adding new data structures to pandas, either custom to pandas
  or provided by 3rd-party libraries, that add new functionality alongside the
  existing code (presuming NumPy physical storage). As one concrete example,
  discussed in more detail below, we can enable missing data in integer pandas
  data by forming a composite data structure consisting of a NumPy array plus a
  bitmap marking the null / not-null values.

Note that neither of these points implies that we are trying to use NumPy
less. We already have large amounts of code that implement algorithms also
found in NumPy (see ``pandas.unique`` or the implementation of ``Series.sum``),
but taking into account pandas's missing data representation, etc. Internally,
we can use NumPy when its computational semantics match those we've chosen for
pandas, and elsewhere we can invoke pandas-specific code.

A major concern here based on these ideas is **preserving NumPy
interoperability**, so I'll examine this topic in some detail next.

Preserving NumPy interoperability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
* Getting and setting bits on modern hardware is CPU-inexpensive. For
  single-pass array operations (like groupbys) on large arrays this may also
  result in better CPU cache utilization (fewer main-memory reads of the
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

Memory accounting
=================

Proper types for strings and some non-numeric data
==================================================

I believe that frequently-occurring data types, such as UTF8 strings, are
important enough to deserve a dedicated logical pandas data type. This will
enable us both to enforce tighter API semantics (i.e. attempts to assign a
non-string into string data will be a ``TypeError``) and improved performance
and memory use under the hood. I will devote an entire section to talking about
strings.

In general, I would be supportive of making Python object (``numpy.object_``
dtype) arrays the solution only for mixed-type arrays and data types for which
pandas has no native handling.

Permitting "other" (non-NumPy) data structures
==============================================



C++11/14 for lowest implementation tier
=======================================

Currently, pandas architecturally is structured as follows:

* Pure Python implementation of internal data structure business logic
* Algorithms in Cython (more often) or C (less often) to accelerate
  computationally-intensive algorithms

While it's overall made pandas easier to develop and maintain internally
(perhaps increasingly less so over time!), this has had a number of drawbacks

Microperformance
~~~~~~~~~~~~~~~~

Microperformance (operations taking 1 microsecond to 1 millisecond) has
suffered considerably as pandas's internals have expanded to accommodate new
use cases. Fairly simple operations, from indexing to summary statistics, may
pass through multiple layers of scaffolding before hitting the lowest tier of
computations. Let's take for example:

.. ipython:: python

   s = pd.Series(np.random.randn(100))
   s.sum()

Profiling ``s.sum()`` with ``%prun`` in IPython, I am seeing 116 function
calls (pandas 0.18.1). Let's look at the microperformance:

.. code-block:: python

   In [14]: timeit s.sum()
   10000 loops, best of 3: 31.7 µs per loop

   In [15]: v = s.values

   In [16]: timeit v.sum()
   1000000 loops, best of 3: 1.07 µs per loop

While a slightly contrived example, the internal data structures and function
dispatch machinery add 30 microseconds of overhead. That may not be a
compelling number, but such a method called 1 million times has an additional
30 seconds of overhead. When you consider microperformance in the context of
custom ``groupby`` operations, for example, this may not be so unrealistic.

3rd-party native API (i.e. Cython and C / C++)
==============================================

Developers of 3rd-party projects (myself included) have often expressed a
desire to be able to inspect, construct, or otherwise manipulate pandas objects
(if even in a limited fashion) in compiled code (Cython, C, or C++).
