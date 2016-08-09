.. _strings:

.. ipython:: python
   :suppress:

   import numpy as np
   import pandas as pd
   np.set_printoptions(precision=4, suppress=True)
   pd.options.display.max_rows = 100

==================================
 Enhanced string / UTF-8 handling
==================================

There are some things we can do to make pandas use less memory and perform
computations significantly faster on string data.

Current string problems
=======================

pandas offers support for columns containing strings (ASCII or Unicode) on a
somewhat ad hoc basis.

* Strings are stored in NumPy arrays of ``PyObject*`` / ``numpy.object_``
  dtype. This has several problems

  * Computations (e.g. ``groupby`` operations) typically utilize a code path
    for generic Python objects. For example comparisons or hashing goes through
    the ``PyObject_*`` C API functions. In addition to harming multithreading
    due to GIL contention (you must acquire the GIL to use these functions),
    these can also be significantly slower than algorithms that operate on
    ``const char*``, potentially taking advantage of hardware optimizations.

  * String arrays often feature many copies of or references to the same
    PyString. Thus, some algorithms may perform redundant computation. Some
    parts of pandas, like ``pandas.read_csv``, make an effort to deduplicate
    strings to free memory and accelerate computations (e.g. if you do ``x ==
    y``, and ``x`` and ``y`` are references to the same ``PyObject*``, Python
    skips comparing their internal data).

    * Note that this is somewhat mitigated by using ``pandas.Categorical``, but
      this is not the default storage mechanism. More on this below.

  * Using ``PyString`` objects and ``PyObject*`` NumPy storage adds non-trivial
    overhead (approximately 24 bytes per unique object, see `this exposition
    <http://www.gahcep.com/python-internals-pyobject/>`_ for a deeper drive) to
    each value.

Possible solution: new non-NumPy string memory layout
=====================================================

My proposed solution to the string conundrum is the following:

* Create a custom string array container type suitable for use in a
  ``pandas.Array``, and a ``pandas.string`` logical data type.
* Require that all strings be encoded as UTF-8.
* By default, represent all string arrays internally as dictionary-encoded
  a.k.a. categorical. Thus, we will typically only ever have 1 copy of any
  given string in an array.
* Store the actual string data in a packed UTF-8 buffer. I have seen this in a
  number of places, but notably it's the way that `Apache Arrow implements
  variable-length collections
  <https://github.com/apache/arrow/blob/master/format/Layout.md#list-type>`_.

Here is one possible C struct-like layout of this container:

.. code-block:: c++

   typedef struct {
     /* Category / dictionary indices into the string data */
     uint32_t* indices;

     /* The encoded string lengths */
     uint32_t* offsets;

     /* The packed UTF-8 data */
     const char* data;

     /* For nullness */
     uint8_t* bitmap;
   } string_array_t;

Here's an example of what the data would look like:

.. code-block:: text

   actual data : ['foo', 'bars', 'foo', null, 'bars']

   indices: [0, 1, 0, 0, 1]

                                    bitmap[0]
   bitmap (read right-to-left): 0 0 0 1 0 1 1 1 |

   offsets: [0, 3, 7]
   data: ['f', 'o', 'o', 'b', 'a', 'r', 's']

Some benefits of this approach include:

* Much better data locality for low-cardinality categorical data
* 8.125 bytes (8 bytes plus 1 bit) of memory overhead per value versus 24 bytes
  (the current)
* The data is already categorical: cast to ``category`` dtype can be perform
  very cheaply and without duplicating the underlying string memory buffer
* Computations like ``groupby`` on dictionary-encoded strings will be as
  performant as those on Categorical currently are.  performant

Some drawbacks

* This memory layout is best used as an immutable representation. Mutating
  slots here becomes more complex. Whether single value assignments or put /
  array-assignment may likely require constructing a new ``data`` buffer
  (either by ``realloc`` or some other copying mechanism). Without a compaction
  / "garbage collection" step on this buffer it will be possible to have "dead"
  memory inside it (for example, if you did ``arr[:] = 'a-new-string-value'``,
  all the existing values would be orphaned).

  * Some systems have addressed this issue by storing all string data in a
    "global string hash table". This is something we could explore, but it
    would add quite a bit of complexity to implement and may not be worthwhile
    at this time.

* Indexing into this data structure to obtain a single Python object will
  probably want to call ``PyUnicode_FromStringAndSize`` to construct a string
  (Python 3, therefore Unicode). This requires a memory allocation, whereas it
  currently only has to do a ``Py_INCREF``.

* Many of pandas's existing algorithms assuming Python objects would need to be
  specialized to take advantage of this new memory layout. This is both a pro
  and a con as it will most likely yield significantly better performance.

Concerns / problems
===================

Preserving code that assumes PyString objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Any alternate UTF-8 string in-memory representation should necessarily be able
to yield Python string objects using ``PyUnicode_FromStringAndSize``. Thus,
code like this could continue to work:

.. ipython:: python

   s = pd.Series(["como estás?"])
   s.map(lambda x: x.upper())

One trade-off is that creating the temporary Python strings is potentially
costly. This could be mitigated for Python ``str`` methods (optimized
array-oriented code path under the hood), but for arbitrary functions you would
have to pay.

Accommodating Non-UTF-8 data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some pandas users will have code that involves various non-UTF-8 Python string
types:

* Native unicode: Py_UCS1, Py_UCS2, Py_UCS4
* Non-UTF-8 PyBytes

.. ipython:: python

   s = pd.Series(["como estás?"])
   s
   s.str.encode('latin-1')
   s.str.encode('latin-1').str.decode('latin-1')

Such data could arise from reading a CSV file in a non-UTF-8 encoding, and you
did not indicate the encoding to ``pandas.read_csv``.

My proposed solution to this is to provide a ``binary`` logical type having the
same physical memory layout as UTF-8 strings, with only the metadata being
different. So you would have the following semantics:

* ``latin1_s = s.encode('latin-1')``: this yields a ``binary`` view and
  allocates new memory.
* ``utf8_s = s.encode('utf-8')``: this is a no-op, but yields a ``binary`` view.
* ``s2 = utf8_s.decode('utf-8')``: this requires using a Unicode codec to
  validate indicated codec.

Indexing and slicing
~~~~~~~~~~~~~~~~~~~~

Storing strings as UTF-8 bytes means that things like this become more
complicated:

.. ipython:: python

   s = pd.Series(["estás está estáis"])
   s.str[9]
   s.str[6:10]

Since UTF-8 is a variable length encoding, finding the logical character by
position will need to make use of the Python C API (expensive, requires
creating new Python objects) or a 3rd party library. We could make use of the
`ICU C++ Libraries <http://site.icu-project.org/>`_ to implement this.
