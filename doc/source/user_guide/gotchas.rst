.. _gotchas:

{{ header }}

********************************
Frequently Asked Questions (FAQ)
********************************

.. _df-memory-usage:

DataFrame memory usage
----------------------
The memory usage of a :class:`DataFrame` (including the index) is shown when calling
the :meth:`~DataFrame.info`. A configuration option, ``display.memory_usage``
(see :ref:`the list of options <options.available>`), specifies if the
:class:`DataFrame` memory usage will be displayed when invoking the :meth:`~DataFrame.info`
method.

For example, the memory usage of the :class:`DataFrame` below is shown
when calling :meth:`~DataFrame.info`:

.. ipython:: python

    dtypes = [
        "int64",
        "float64",
        "datetime64[ns]",
        "timedelta64[ns]",
        "complex128",
        "object",
        "bool",
    ]
    n = 5000
    data = {t: np.random.randint(100, size=n).astype(t) for t in dtypes}
    df = pd.DataFrame(data)
    df["categorical"] = df["object"].astype("category")

    df.info()

The ``+`` symbol indicates that the true memory usage could be higher, because
pandas does not count the memory used by values in columns with
``dtype=object``.

Passing ``memory_usage='deep'`` will enable a more accurate memory usage report,
accounting for the full usage of the contained objects. This is optional
as it can be expensive to do this deeper introspection.

.. ipython:: python

   df.info(memory_usage="deep")

By default the display option is set to ``True`` but can be explicitly
overridden by passing the ``memory_usage`` argument when invoking :meth:`~DataFrame.info`.

The memory usage of each column can be found by calling the
:meth:`~DataFrame.memory_usage` method. This returns a :class:`Series` with an index
represented by column names and memory usage of each column shown in bytes. For
the :class:`DataFrame` above, the memory usage of each column and the total memory
usage can be found with the :meth:`~DataFrame.memory_usage` method:

.. ipython:: python

    df.memory_usage()

    # total memory usage of dataframe
    df.memory_usage().sum()

By default the memory usage of the :class:`DataFrame` index is shown in the
returned :class:`Series`, the memory usage of the index can be suppressed by passing
the ``index=False`` argument:

.. ipython:: python

    df.memory_usage(index=False)

The memory usage displayed by the :meth:`~DataFrame.info` method utilizes the
:meth:`~DataFrame.memory_usage` method to determine the memory usage of a
:class:`DataFrame` while also formatting the output in human-readable units (base-2
representation; i.e. 1KB = 1024 bytes).

See also :ref:`Categorical Memory Usage <categorical.memory>`.

.. _gotchas.truth:

Using if/truth statements with pandas
-------------------------------------

pandas follows the NumPy convention of raising an error when you try to convert
something to a ``bool``. This happens in an ``if``-statement or when using the
boolean operations: ``and``, ``or``, and ``not``. It is not clear what the result
of the following code should be:

.. code-block:: python

    >>> if pd.Series([False, True, False]):
    ...     pass

Should it be ``True`` because it's not zero-length, or ``False`` because there
are ``False`` values? It is unclear, so instead, pandas raises a ``ValueError``:

.. ipython:: python
    :okexcept:

    if pd.Series([False, True, False]):
        print("I was true")

You need to explicitly choose what you want to do with the :class:`DataFrame`, e.g.
use :meth:`~DataFrame.any`, :meth:`~DataFrame.all` or :meth:`~DataFrame.empty`.
Alternatively, you might want to compare if the pandas object is ``None``:

.. ipython:: python

    if pd.Series([False, True, False]) is not None:
        print("I was not None")


Below is how to check if any of the values are ``True``:

.. ipython:: python

    if pd.Series([False, True, False]).any():
        print("I am any")

Bitwise Boolean
~~~~~~~~~~~~~~~

Bitwise boolean operators like ``==`` and ``!=`` return a boolean :class:`Series`
which performs an element-wise comparison when compared to a scalar.

.. ipython:: python

   s = pd.Series(range(5))
   s == 4

See :ref:`boolean comparisons<basics.compare>` for more examples.

Using the ``in`` operator
~~~~~~~~~~~~~~~~~~~~~~~~~

Using the Python ``in`` operator on a :class:`Series` tests for membership in the
**index**, not membership among the values.

.. ipython:: python

    s = pd.Series(range(5), index=list("abcde"))
    2 in s
    'b' in s

If this behavior is surprising, keep in mind that using ``in`` on a Python
dictionary tests keys, not values, and :class:`Series` are dict-like.
To test for membership in the values, use the method :meth:`~pandas.Series.isin`:

.. ipython:: python

    s.isin([2])
    s.isin([2]).any()

For :class:`DataFrame`, likewise, ``in`` applies to the column axis,
testing for membership in the list of column names.

.. _gotchas.udf-mutation:

Mutating with User Defined Function (UDF) methods
-------------------------------------------------

This section applies to pandas methods that take a UDF. In particular, the methods
:meth:`DataFrame.apply`, :meth:`DataFrame.aggregate`, :meth:`DataFrame.transform`, and
:meth:`DataFrame.filter`.

It is a general rule in programming that one should not mutate a container
while it is being iterated over. Mutation will invalidate the iterator,
causing unexpected behavior. Consider the example:

.. ipython:: python

   values = [0, 1, 2, 3, 4, 5]
   n_removed = 0
   for k, value in enumerate(values):
       idx = k - n_removed
       if value % 2 == 1:
           del values[idx]
           n_removed += 1
       else:
           values[idx] = value + 1
   values

One probably would have expected that the result would be ``[1, 3, 5]``.
When using a pandas method that takes a UDF, internally pandas is often
iterating over the
:class:`DataFrame` or other pandas object. Therefore, if the UDF mutates (changes)
the :class:`DataFrame`, unexpected behavior can arise.

Here is a similar example with :meth:`DataFrame.apply`:

.. ipython:: python
   :okexcept:

   def f(s):
       s.pop("a")
       return s

   df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
   df.apply(f, axis="columns")

To resolve this issue, one can make a copy so that the mutation does
not apply to the container being iterated over.

.. ipython:: python

   values = [0, 1, 2, 3, 4, 5]
   n_removed = 0
   for k, value in enumerate(values.copy()):
       idx = k - n_removed
       if value % 2 == 1:
           del values[idx]
           n_removed += 1
       else:
           values[idx] = value + 1
   values

.. ipython:: python

   def f(s):
       s = s.copy()
       s.pop("a")
       return s

   df = pd.DataFrame({"a": [1, 2, 3], 'b': [4, 5, 6]})
   df.apply(f, axis="columns")

Missing value representation for NumPy types
--------------------------------------------

``np.nan`` as the ``NA`` representation for NumPy types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For lack of ``NA`` (missing) support from the ground up in NumPy and Python in
general, ``NA`` could have been represented with:

* A *masked array* solution: an array of data and an array of boolean values
  indicating whether a value is there or is missing.
* Using a special sentinel value, bit pattern, or set of sentinel values to
  denote ``NA`` across the dtypes.

The special value ``np.nan`` (Not-A-Number) was chosen as the ``NA`` value for NumPy types, and there are API
functions like :meth:`DataFrame.isna` and :meth:`DataFrame.notna` which can be used across the dtypes to
detect NA values. However, this choice has a downside of coercing missing integer data as float types as
shown in :ref:`gotchas.intna`.

``NA`` type promotions for NumPy types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When introducing NAs into an existing :class:`Series` or :class:`DataFrame` via
:meth:`~Series.reindex` or some other means, boolean and integer types will be
promoted to a different dtype in order to store the NAs. The promotions are
summarized in this table:

.. csv-table::
   :header: "Typeclass","Promotion dtype for storing NAs"
   :widths: 40,60

   ``floating``, no change
   ``object``, no change
   ``integer``, cast to ``float64``
   ``boolean``, cast to ``object``

.. _gotchas.intna:

Support for integer ``NA``
~~~~~~~~~~~~~~~~~~~~~~~~~~

In the absence of high performance ``NA`` support being built into NumPy from
the ground up, the primary casualty is the ability to represent NAs in integer
arrays. For example:

.. ipython:: python

   s = pd.Series([1, 2, 3, 4, 5], index=list("abcde"))
   s
   s.dtype

   s2 = s.reindex(["a", "b", "c", "f", "u"])
   s2
   s2.dtype

This trade-off is made largely for memory and performance reasons, and also so
that the resulting :class:`Series` continues to be "numeric".

If you need to represent integers with possibly missing values, use one of
the nullable-integer extension dtypes provided by pandas or pyarrow

* :class:`Int8Dtype`
* :class:`Int16Dtype`
* :class:`Int32Dtype`
* :class:`Int64Dtype`
* :class:`ArrowDtype`

.. ipython:: python

   s_int = pd.Series([1, 2, 3, 4, 5], index=list("abcde"), dtype=pd.Int64Dtype())
   s_int
   s_int.dtype

   s2_int = s_int.reindex(["a", "b", "c", "f", "u"])
   s2_int
   s2_int.dtype

   s_int_pa = pd.Series([1, 2, None], dtype="int64[pyarrow]")
   s_int_pa

See :ref:`integer_na` and :ref:`pyarrow` for more.

Why not make NumPy like R?
~~~~~~~~~~~~~~~~~~~~~~~~~~

Many people have suggested that NumPy should simply emulate the ``NA`` support
present in the more domain-specific statistical programming language `R
<https://www.r-project.org/>`__. Part of the reason is the
`NumPy type hierarchy <https://numpy.org/doc/stable/user/basics.types.html>`__.

The R language, by contrast, only has a handful of built-in data types:
``integer``, ``numeric`` (floating-point), ``character``, and
``boolean``. ``NA`` types are implemented by reserving special bit patterns for
each type to be used as the missing value. While doing this with the full NumPy
type hierarchy would be possible, it would be a more substantial trade-off
(especially for the 8- and 16-bit data types) and implementation undertaking.

However, R ``NA`` semantics are now available by using masked NumPy types such as :class:`Int64Dtype`
or PyArrow types (:class:`ArrowDtype`).


Differences with NumPy
----------------------
For :class:`Series` and :class:`DataFrame` objects, :meth:`~DataFrame.var` normalizes by
``N-1`` to produce `unbiased estimates of the population variance <https://en.wikipedia.org/wiki/Bias_of_an_estimator>`__, while NumPy's
:meth:`numpy.var` normalizes by N, which measures the variance of the sample. Note that
:meth:`~DataFrame.cov` normalizes by ``N-1`` in both pandas and NumPy.

.. _gotchas.thread-safety:

Thread-safety
-------------

pandas is not 100% thread safe. The known issues relate to
the :meth:`~DataFrame.copy` method. If you are doing a lot of copying of
:class:`DataFrame` objects shared among threads, we recommend holding locks inside
the threads where the data copying occurs.

See `this link <https://stackoverflow.com/questions/13592618/python-pandas-dataframe-thread-safe>`__
for more information.


Byte-ordering issues
--------------------
Occasionally you may have to deal with data that were created on a machine with
a different byte order than the one on which you are running Python. A common
symptom of this issue is an error like::

    Traceback
        ...
    ValueError: Big-endian buffer not supported on little-endian compiler

To deal
with this issue you should convert the underlying NumPy array to the native
system byte order *before* passing it to :class:`Series` or :class:`DataFrame`
constructors using something similar to the following:

.. ipython:: python

   x = np.array(list(range(10)), ">i4")  # big endian
   newx = x.byteswap().view(x.dtype.newbyteorder())  # force native byteorder
   s = pd.Series(newx)

See `the NumPy documentation on byte order
<https://numpy.org/doc/stable/user/byteswapping.html>`__ for more
details.
