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
hardest part is coming up with clear language and definitions for concepts so
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

For historical reasons, **pandas never developed a clear or clean semantic
separation in its user API between logical and physical data types**. Also, the
addition of new, pandas-only "synthetic" dtypes that are unknown to NumPy (like
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
  categorical types, since a particular instance of ``CategoricalDtype`` may
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

* We can start to think about improved behavior around data ownership (like
  copy-on-write) which may yield many benefits. I will write a dedicated
  section about this.

Note that neither of these points implies that we are trying to use NumPy
less. We already have large amounts of code that implement algorithms similar
to those found in NumPy (e.g. ``pandas.unique`` or the implementation of
``Series.sum``), but taking into account pandas's missing data representation,
etc. Internally, we can use NumPy when its computational semantics match those
we've chosen for pandas, and elsewhere we can invoke pandas-specific code.

A major concern here based on these ideas is **preserving NumPy
interoperability**, so I'll examine this topic in some detail next.

Preserving NumPy interoperability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some of types of intended interoperability between NumPy and pandas are as
follows:

* Users can obtain the a ``numpy.ndarray`` (possibly a view depending on the
  internal block structure, more on this soon) in constant time and without
  copying the actual data. This has a couple other implications

  * Changes made to this array will be reflected in the source pandas object.
  * If you write C extension code (possibly in Cython) and respect pandas's
    missing data details, you can invoke certain kinds of fast custom code on
    pandas data (but it's somewhat inflexible -- see the latest discussion on
    adding a native code API to pandas).

* NumPy ufuncs (like ``np.sqrt`` or ``np.log``) can be invoked on
  pandas objects like Series and DataFrame

* ``numpy.asarray`` will always yield some array, even if it discards metadata
  or has to create a new array. For example ``asarray`` invoked on
  ``pandas.Categorical`` yields a reconstructed array (rather than either the
  categories or codes internal arrays)

* Many NumPy methods designed to work on subclasses (or duck-typed classes) of
  ``ndarray`` may be used. For example ``numpy.sum`` may be used on a Series
  even though it does not invoke NumPy's internal C sum algorithm. This means
  that a Series may be used as an interchangeable argument in a large set of
  functions that only know about NumPy arrays.

By and large, I think much of this can be preserved, but there will be some API
breakage.

If we add more composite data structures (Categorical can be thought of as
one existing composite data structure) to pandas or alternate non-NumPy data
structures, there will be cases where the semantic information in a Series
cannot be adequately represented in a NumPy array.

As one example, if we add pandas-only missing data support to integer and
boolean data (a long requested feature), calling ``np.asarray`` on such data
may not have well-defined behavior. As present, pandas is implicitly converting
these types to ``float64`` (see more below), which isn't too great. A decision
does not need to be made now, but the benefits of solving this long-standing
issue may merit breaking ``asarray`` as long as we provide an explicit way to
obtain the original casted ``float64`` NumPy array (with ``NaN`` for NULL/NA
values)

For pandas data that does not step outside NumPy's semantic realm, we can
continue to provide zero-copy views in many cases.

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
than the current implicit conversion which could yield data loss (for integers
falling outside the exact representable range for ``float64``).

Removal of BlockManager / new DataFrame internals
=================================================

Deep inside the belly pandas objects, there is a data structure called
``BlockManager`` which, at a high level, is responsible for managing the
physical arrays where the data inside a Series or DataFrame is looked
after (also Panel / PanelND structure, even though these are on their way to
deprecation).

While this data structure has served pandas well since its birth 5 years ago
(Summer 2011), it has a number of problems that make its removal and
replacement with something else an attractive option.

The goal of this section is to explain what the BlockManager is, why it exists
at all, and why we should consider removing it.

What is ``BlockManager`` and why does it exist?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reason that ``BlockManager`` exists at all goes back to some ancient pandas
history. Originally, the data in ``pandas.DataFrame`` was stored in a Python
``dict`` object. If you pull up pandas 0.1 or 0.2, you will see this.

Since the business logic of pandas's internals was originally implemented in
pure Python, as it is still is (but much larger / more complex), there was a
marked performance difference between column-oriented operations and
row-oriented operations. The reason for this is not really a memory layout
issue (NumPy users know about how contiguous memory access produces much better
performance) so much as a reliance on NumPy's two-dimensional array operations
for carrying out pandas's computations. So, to do anything row oriented on an
all-numeric DataFrame, pandas would concatenate all of the columns together
(using ``numpy.vstack`` or ``numpy.hstack``) then use array broadcasting or
methods like ``ndarray.sum`` (combined with ``np.isnan`` to mind missing data)
to carry out certain operations.

1. pandas's early users (i.e. AQR employees) beseeched me to address this
   performance issue. Thus ``DataMatrix`` was created, a roughly API-equivalent
   object whose internal storage was a 2D NumPy array, intended to be of a
   homogeneous type (e.g. ``numpy.float64``). The downside of this was that if
   you inserted a string column, everything would become ``numpy.object_``
   dtype. Users did not like that.

2. It had become apparent that the dichotomy between DataFrame and DataMatrix
   (and when to use each) was harming pandas's adoption and confusing users. So
   I set about creating a hybrid data structure that had "the best of both
   worlds".

3. The idea was that the BlockManager would track collections of NumPy arrays
   having the same dtype, particular as columns were inserted or removed
   (i.e. the *building* phase of the DataFrame's lifetime).

4. When you would invoke an operation that benefited from a single
   *consolidated* 2-dimensional ndarray of say ``float64`` dtype (for example:
   using ``reindex`` or performing a row-oriented operation), the BlockManager
   would glue together its accumulated pieces to create a single 2D ndarray of
   each data type. This is called **consolidation** in the codebase.

5. Since in practice, heterogeneous DataFrames had different types interspersed
   amongst their columns, the BlockManager maintains a mapping between the
   absolute column position and the relative position within the type-specific
   2D "block".

6. Over time, the BlockManager has been generalized for the 1 through N
   dimensional cases, not just the 2D case, so that even Series has a lean
   "SingleBlockManager" internally.

Drawbacks of BlockManager
~~~~~~~~~~~~~~~~~~~~~~~~~

While this data structure has enabled pandas to make it this far in life, it
has a number of drawbacks (not a complete list):

1. **Code complexity**: this has manifested in a number of ways (and probably
   others that I'm missing)

   * Making some of the most important algorithms in pandas fast, like joins
     and reshape operations, requires carefully constructing the precise block
     structure of the output DataFrame so that no further copying or
     consolidation will take place.

   * Adding new custom data types to DataFrame and not losing their metadata
     (e.g. time zones or categories) has had a sort of "fan out" effect
     touching numerous parts of the BlockManager internals.

2. **Loss of user visibility into memory use and memory layout**: With large
   data sets, some "naively" constructed DataFrame objects (e.g. from a dict of
   ndarrays) can produce a memory-doubling effect that may cause out-of-memory
   errors. Also, consolidated blocks can (depending on the version of pandas)
   result in columns having strided / non-contiguous data, resulting in
   degraded performance in column-oriented operations.

3. **Unavoidable consolidation**: Fairly common operations, like ``read_csv``,
   may require a consolidation step after completion, which for large data may
   result in performance or memory overhead (similar to the above bullet
   point).

4. **Microperformance issues / indexing slowness**: since a DataFrame can be a
   sort of many-layered onion, many common pandas operations may weave through
   dozens of different functions navigating the structure of the object and
   producing the appropriate output. I will talk more about microperformance
   later.

Replacing BlockManager without weakening pandas
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Our goal in replacing BlockManager would be to achieve:

* Substantially simpler code
* Easier extensibility with new logical types
* Performance on par (or better) the current implementation
* Better user control over memory use and layout
* Improved microperformance

I believe we can do this, but it's will require a significant inversion of the
internal code architecture to involve a more native code and less interpreted
Python. For example, it will be difficult or impossible to achieve comparable
performance in row-oriented operations (on consolidated DataFrame objects) with
pure Python code.

In the next section, I will start making my case for creating a "native core"
library where we can assemble the low level data structures, logical types, and
memory management for pandas. Additionally, we would want to port much of
pandas's helper Cython code to live inside this library and operate directly on
the internal data structures rather than being orchestrated from the Python
interpreter level.

Building "libpandas" in C++11/14 for lowest level implementation tier
=====================================================================

Currently, pandas architecturally is structured as follows:

* Pure Python implementation of internal data structure business logic
* Algorithms in Cython (more often) or C (less often) to accelerate
  computationally-intensive algorithms

While it's overall made pandas easier to develop and maintain internally
(perhaps increasingly less so over time!), this has had a number of drawbacks
as we've discussed. I mentioned microperformance above, so about that:

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

.. code-block:: text

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

C or C++ (C++11, to be specific)?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the risk of instigating a religious programming language debate, pandas's
use of Cython in many places is very C++-like:

* Generic programming through manual code generation (now using tempita)
  instead of templates
* Auxiliary types and data structures as ``cdef class`` extension types
* Relying on Python's reference counting for garbage collection and cleanup
  after exceptions are raised. The "blend C and Cython" style has aided
  developer productivity.

I argue that judicious and responsible use of modern C++ (and following a
reasonable style guide like `Google's guide
<http://google.github.io/styleguide/cppguide.html>`_, or some slight variation)
will enable us to:

* Simplify our existing Cython codebase by using templates (and very limited,
  template metaprogramming)

* Easier generic programming / inlining of data-type specific logic at compile
  time.

* Use RAII (exception-safe allocation) and smart pointers (``std::unique_ptr``
  and ``std::shared_ptr``) to simplify memory management

* Define performant C++ classes modeling the current internals, with various
  mechanisms for code reuse or type-specific dynamic dispatch (i.e. through
  template classes, CRTP, or simply virtual functions).

* Use C++11 standard library concurrency tools to more easily create concurrent
  / multithreaded implementations of common pandas algorithms.

By pushing down much of the business logic into C++ (with use of the Python and
NumPy C API where relevant), we'll be able to achieve macroperformance on par
or better than the current BlockManager-based implementation and handily better
microperformance in indexing and simple analytics.

``pandas.Array`` types
~~~~~~~~~~~~~~~~~~~~~~

My gut feeling is that we would want to create relatively simple container
classes having a common ``pandas::Array`` base type in C++, each of which
models a particular logical type. Each array type would have a corresponding
logical type implementation, in the vein of:

.. code-block:: c++

   class Array {
     // public API omitted
     private:
       std::shared_ptr<DataType> type_;
   }

   class CategoricalType : public DataType {
     // implementation

     private:
       std::shared_ptr<Array> categories_;
   };

   class CategoricalArray : public Array {
     public:
       std::shared_ptr<Array> codes() const;
       std::shared_ptr<Array> categories() const;
       // rest of implementation omitted
   };

An array containing a NumPy array will invoke ``Py_DECREF`` in its destructor,
so that after construction one can proceed largely with C++ programming
semantics without much need for manual memory management.

These Array types would be wrapped and exposed to pandas developers (probably
in Cython).

Index types
~~~~~~~~~~~

Like pandas's current code structure, Index types would be composed from the
Array types and some additional data structures (hash tables) for lookups and
other index operations. These can be similarly exposed to the world via Cython
(and wrapped in a convenient pandas.Index class).

``pandas.Table``
~~~~~~~~~~~~~~~~

My recommendation is to decommission the BlockManager in favor of a much
simpler low-level Table class, which operates more similarly to an R data.frame
(e.g. no row index). This would look something like

.. code-block:: c++

   class Table {
     public:
       std::shared_ptr<Array> GetColumn(int i);
       void SetColumn(int i, const std::shared_ptr<Array>& arr);

       // rest of public API omitted
     private:
       // Column index, possibly not necessary
       std::shared_ptr<Index> columns_;

       // List of arrays
       std::vector<std::shared_ptr<Array>> data_;
   };

Operators and dynamic dispatch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Under this proposed class structure, it may not make sense to add operations as
class methods. We could possibly do something like:

.. code-block:: c++

   #include "pandas/dispatch.h"

   // other includes omitted

   using ArrayRef = std::shared_ptr<Array>;

   template <typename U, typename V>
   inline ArrayRef TakeImpl(U, V) {
     // Implementation omitted
   }

   ArrayRef Take(ArrayRef values, ArrayRef indices) {
     return Dispatch<TakeImpl>(values, indices);
   }

Here, the Dispatch template would generate the matrix of logical type
combinations, some of which might throw a not implemented exception.

There's other approaches to dealing with runtime dispatch that don't feature
too much overhead.

Memory accounting
~~~~~~~~~~~~~~~~~

If pandas's internals are encapsulated in C++ classes inside the libpandas core
library, we could atomically track all memory allocations and deallocations to
produce a precise accounting of the number of bytes that pandas has currently
allocated (that are not opaque, so Python objects would only include their
``PyObject**`` array footprint).

Development toolchain
~~~~~~~~~~~~~~~~~~~~~

Introducing C++11 to pandas's development toolchain will add quite a bit of
complexity for developers, especially compared with pandas's current Cython and
C codebase which basically builds out of the box for most people. It would be
better for cross-platform support to use CMake than something else (distutils
doesn't have adequate support for C++).

Logical types for strings and possibly other non-numeric data
=============================================================

I believe that frequently-occurring data types, such as UTF8 strings, are
important enough to deserve a dedicated logical pandas data type. This will
enable us both to enforce tighter API semantics (i.e. attempts to assign a
non-string into string data will be a ``TypeError``) and improved performance
and memory use under the hood. I will devote an entire section to talking about
strings.

In general, I would be supportive of making Python object (``numpy.object_``
dtype) arrays the solution only for mixed-type arrays and data types for which
pandas has no native handling.

3rd-party native API (i.e. Cython and C / C++)
==============================================

Developers of 3rd-party projects (myself included) have often expressed a
desire to be able to inspect, construct, or otherwise manipulate pandas objects
(if even in a limited fashion) in compiled code (Cython, C, or C++).

Per the discussion of libpandas and a native core, I would propose the
following:

* Define public-facing ``.pxd`` files that allow developers to use ``cimport``
  and get access to pandas's internal extension types.
* Define factory function that enable fully formed Series and DataFrame objects
  to be constructed either by Cython API calls or potentially also C++
  libpandas API calls.
* Provide Cython APIs for 3rd-party developers to obtain pointers to access the
  underlying C++ objects contained in the wrapper Python objects
