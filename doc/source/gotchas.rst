.. currentmodule:: pandas
.. _gotchas:

********************************
Frequently Asked Questions (FAQ)
********************************

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   np.set_printoptions(precision=4, suppress=True)
   import pandas as pd
   pd.options.display.max_rows = 15
   import matplotlib
   matplotlib.style.use('ggplot')
   import matplotlib.pyplot as plt
   plt.close('all')

.. _df-memory-usage:

DataFrame memory usage
----------------------
As of pandas version 0.15.0, the memory usage of a dataframe (including
the index) is shown when accessing the ``info`` method of a dataframe. A
configuration option, ``display.memory_usage`` (see :ref:`options`),
specifies if the dataframe's memory usage will be displayed when
invoking the ``df.info()`` method.

For example, the memory usage of the dataframe below is shown
when calling ``df.info()``:

.. ipython:: python

    dtypes = ['int64', 'float64', 'datetime64[ns]', 'timedelta64[ns]',
              'complex128', 'object', 'bool']
    n = 5000
    data = dict([ (t, np.random.randint(100, size=n).astype(t))
                    for t in dtypes])
    df = pd.DataFrame(data)
    df['categorical'] = df['object'].astype('category')

    df.info()

The ``+`` symbol indicates that the true memory usage could be higher, because
pandas does not count the memory used by values in columns with
``dtype=object``.

.. versionadded:: 0.17.1

Passing ``memory_usage='deep'`` will enable a more accurate memory usage report,
that accounts for the full usage of the contained objects. This is optional
as it can be expensive to do this deeper introspection.

.. ipython:: python

   df.info(memory_usage='deep')

By default the display option is set to ``True`` but can be explicitly
overridden by passing the ``memory_usage`` argument when invoking ``df.info()``.

The memory usage of each column can be found by calling the ``memory_usage``
method. This returns a Series with an index represented by column names
and memory usage of each column shown in bytes. For the dataframe above,
the memory usage of each column and the total memory usage of the
dataframe can be found with the memory_usage method:

.. ipython:: python

    df.memory_usage()

    # total memory usage of dataframe
    df.memory_usage().sum()

By default the memory usage of the dataframe's index is shown in the
returned Series, the memory usage of the index can be suppressed by passing
the ``index=False`` argument:

.. ipython:: python

    df.memory_usage(index=False)

The memory usage displayed by the ``info`` method utilizes the
``memory_usage`` method to determine the memory usage of a dataframe
while also formatting the output in human-readable units (base-2
representation; i.e., 1KB = 1024 bytes).

See also :ref:`Categorical Memory Usage <categorical.memory>`.

.. _gotchas.truth:

Using If/Truth Statements with pandas
-------------------------------------

pandas follows the numpy convention of raising an error when you try to convert something to a ``bool``.
This happens in a ``if`` or when using the boolean operations, ``and``, ``or``, or ``not``.  It is not clear
what the result of

.. code-block:: python

    >>> if pd.Series([False, True, False]):
         ...

should be. Should it be ``True`` because it's not zero-length? ``False`` because there are ``False`` values?
It is unclear, so instead, pandas raises a ``ValueError``:

.. code-block:: python

    >>> if pd.Series([False, True, False]):
        print("I was true")
    Traceback
        ...
    ValueError: The truth value of an array is ambiguous. Use a.empty, a.any() or a.all().


If you see that, you need to explicitly choose what you want to do with it (e.g., use `any()`, `all()` or `empty`).
or, you might want to compare if the pandas object is ``None``

.. code-block:: python

    >>> if pd.Series([False, True, False]) is not None:
           print("I was not None")
    >>> I was not None


or return if ``any`` value is ``True``.

.. code-block:: python

    >>> if pd.Series([False, True, False]).any():
           print("I am any")
    >>> I am any

To evaluate single-element pandas objects in a boolean context, use the method ``.bool()``:

.. ipython:: python

   pd.Series([True]).bool()
   pd.Series([False]).bool()
   pd.DataFrame([[True]]).bool()
   pd.DataFrame([[False]]).bool()

Bitwise boolean
~~~~~~~~~~~~~~~

Bitwise boolean operators like ``==`` and ``!=`` return a boolean ``Series``,
which is almost always what you want anyways.

.. code-block:: python

   >>> s = pd.Series(range(5))
   >>> s == 4
   0    False
   1    False
   2    False
   3    False
   4     True
   dtype: bool

See :ref:`boolean comparisons<basics.compare>` for more examples.

Using the ``in`` operator
~~~~~~~~~~~~~~~~~~~~~~~~~

Using the Python ``in`` operator on a Series tests for membership in the
index, not membership among the values.

.. ipython::

    s = pd.Series(range(5), index=list('abcde'))
    2 in s
    'b' in s

If this behavior is surprising, keep in mind that using ``in`` on a Python
dictionary tests keys, not values, and Series are dict-like.
To test for membership in the values, use the method :func:`~pandas.Series.isin`:

.. ipython::

    s.isin([2])
    s.isin([2]).any()

For DataFrames, likewise, ``in`` applies to the column axis,
testing for membership in the list of column names.

``NaN``, Integer ``NA`` values and ``NA`` type promotions
---------------------------------------------------------

Choice of ``NA`` representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For lack of ``NA`` (missing) support from the ground up in NumPy and Python in
general, we were given the difficult choice between either

- A *masked array* solution: an array of data and an array of boolean values
  indicating whether a value is there or is missing
- Using a special sentinel value, bit pattern, or set of sentinel values to
  denote ``NA`` across the dtypes

For many reasons we chose the latter. After years of production use it has
proven, at least in my opinion, to be the best decision given the state of
affairs in NumPy and Python in general. The special value ``NaN``
(Not-A-Number) is used everywhere as the ``NA`` value, and there are API
functions ``isna`` and ``notna`` which can be used across the dtypes to
detect NA values.

However, it comes with it a couple of trade-offs which I most certainly have
not ignored.

.. _gotchas.intna:

Support for integer ``NA``
~~~~~~~~~~~~~~~~~~~~~~~~~~

In the absence of high performance ``NA`` support being built into NumPy from
the ground up, the primary casualty is the ability to represent NAs in integer
arrays. For example:

.. ipython:: python

   s = pd.Series([1, 2, 3, 4, 5], index=list('abcde'))
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

While this may seem like a heavy trade-off, I have found very few cases where
this is an issue in practice i.e. storing values greater than 2**53. Some
explanation for the motivation is in the next section.

Why not make NumPy like R?
~~~~~~~~~~~~~~~~~~~~~~~~~~

Many people have suggested that NumPy should simply emulate the ``NA`` support
present in the more domain-specific statistical programming language `R
<https://r-project.org>`__. Part of the reason is the NumPy type hierarchy:

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


Differences with NumPy
----------------------
For Series and DataFrame objects, ``var`` normalizes by ``N-1`` to produce
unbiased estimates of the sample variance, while NumPy's ``var`` normalizes
by N, which measures the variance of the sample. Note that ``cov``
normalizes by ``N-1`` in both pandas and NumPy.


Thread-safety
-------------

As of pandas 0.11, pandas is not 100% thread safe. The known issues relate to
the ``DataFrame.copy`` method. If you are doing a lot of copying of DataFrame
objects shared among threads, we recommend holding locks inside the threads
where the data copying occurs.

See `this link <https://stackoverflow.com/questions/13592618/python-pandas-dataframe-thread-safe>`__
for more information.


Byte-Ordering Issues
--------------------
Occasionally you may have to deal with data that were created on a machine with
a different byte order than the one on which you are running Python. A common symptom of this issue is an error like

.. code-block:: python

    Traceback
        ...
    ValueError: Big-endian buffer not supported on little-endian compiler

To deal
with this issue you should convert the underlying NumPy array to the native
system byte order *before* passing it to Series/DataFrame/Panel constructors
using something similar to the following:

.. ipython:: python

   x = np.array(list(range(10)), '>i4') # big endian
   newx = x.byteswap().newbyteorder() # force native byteorder
   s = pd.Series(newx)

See `the NumPy documentation on byte order
<https://docs.scipy.org/doc/numpy/user/basics.byteswapping.html>`__ for more
details.
