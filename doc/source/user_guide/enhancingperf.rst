.. _enhancingperf:

{{ header }}

*********************
Enhancing performance
*********************

In this part of the tutorial, we will investigate how to speed up certain
functions operating on pandas :class:`DataFrame` using Cython, Numba and :func:`pandas.eval`.
Generally, using Cython and Numba can offer a larger speedup than using :func:`pandas.eval`
but will require a lot more code.

.. note::

   In addition to following the steps in this tutorial, users interested in enhancing
   performance are highly encouraged to install the
   :ref:`recommended dependencies<install.recommended_dependencies>` for pandas.
   These dependencies are often not installed by default, but will offer speed
   improvements if present.

.. _enhancingperf.cython:

Cython (writing C extensions for pandas)
----------------------------------------

For many use cases writing pandas in pure Python and NumPy is sufficient. In some
computationally heavy applications however, it can be possible to achieve sizable
speed-ups by offloading work to `cython <https://cython.org/>`__.

This tutorial assumes you have refactored as much as possible in Python, for example
by trying to remove for-loops and making use of NumPy vectorization. It's always worth
optimising in Python first.

This tutorial walks through a "typical" process of cythonizing a slow computation.
We use an `example from the Cython documentation <https://docs.cython.org/en/latest/src/quickstart/cythonize.html>`__
but in the context of pandas. Our final cythonized solution is around 100 times
faster than the pure Python solution.

.. _enhancingperf.pure:

Pure Python
~~~~~~~~~~~

We have a :class:`DataFrame` to which we want to apply a function row-wise.

.. ipython:: python

   df = pd.DataFrame(
       {
           "a": np.random.randn(1000),
           "b": np.random.randn(1000),
           "N": np.random.randint(100, 1000, (1000)),
           "x": "x",
       }
   )
   df

Here's the function in pure Python:

.. ipython:: python

   def f(x):
       return x * (x - 1)


   def integrate_f(a, b, N):
       s = 0
       dx = (b - a) / N
       for i in range(N):
           s += f(a + i * dx)
       return s * dx

We achieve our result by using :meth:`DataFrame.apply` (row-wise):

.. ipython:: python

   %timeit df.apply(lambda x: integrate_f(x["a"], x["b"], x["N"]), axis=1)

Let's take a look and see where the time is spent during this operation
using the `prun ipython magic function <https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-prun>`__:

.. ipython:: python

   # most time consuming 4 calls
   %prun -l 4 df.apply(lambda x: integrate_f(x["a"], x["b"], x["N"]), axis=1)  # noqa E999

By far the majority of time is spend inside either ``integrate_f`` or ``f``,
hence we'll concentrate our efforts cythonizing these two functions.

.. _enhancingperf.plain:

Plain Cython
~~~~~~~~~~~~

First we're going to need to import the Cython magic function to IPython:

.. ipython:: python
   :okwarning:

   %load_ext Cython


Now, let's simply copy our functions over to Cython:

.. ipython::

   In [2]: %%cython
      ...: def f_plain(x):
      ...:     return x * (x - 1)
      ...: def integrate_f_plain(a, b, N):
      ...:     s = 0
      ...:     dx = (b - a) / N
      ...:     for i in range(N):
      ...:         s += f_plain(a + i * dx)
      ...:     return s * dx
      ...:


.. ipython:: python

   %timeit df.apply(lambda x: integrate_f_plain(x["a"], x["b"], x["N"]), axis=1)

This has improved the performance compared to the pure Python approach by one-third.

.. _enhancingperf.type:

Declaring C types
~~~~~~~~~~~~~~~~~

We can annotate the function variables and return types as well as use ``cdef``
and ``cpdef`` to improve performance:

.. ipython::

   In [3]: %%cython
      ...: cdef double f_typed(double x) except? -2:
      ...:     return x * (x - 1)
      ...: cpdef double integrate_f_typed(double a, double b, int N):
      ...:     cdef int i
      ...:     cdef double s, dx
      ...:     s = 0
      ...:     dx = (b - a) / N
      ...:     for i in range(N):
      ...:         s += f_typed(a + i * dx)
      ...:     return s * dx
      ...:

.. ipython:: python

   %timeit df.apply(lambda x: integrate_f_typed(x["a"], x["b"], x["N"]), axis=1)

Annotating the functions with C types yields an over ten times performance improvement compared to
the original Python implementation.

.. _enhancingperf.ndarray:

Using ndarray
~~~~~~~~~~~~~

When re-profiling, time is spent creating a :class:`Series` from each row, and calling ``__getitem__`` from both
the index and the series (three times for each row). These Python function calls are expensive and
can be improved by passing an ``np.ndarray``.

.. ipython:: python

   %prun -l 4 df.apply(lambda x: integrate_f_typed(x["a"], x["b"], x["N"]), axis=1)

.. ipython::

   In [4]: %%cython
      ...: cimport numpy as np
      ...: import numpy as np
      ...: cdef double f_typed(double x) except? -2:
      ...:     return x * (x - 1)
      ...: cpdef double integrate_f_typed(double a, double b, int N):
      ...:     cdef int i
      ...:     cdef double s, dx
      ...:     s = 0
      ...:     dx = (b - a) / N
      ...:     for i in range(N):
      ...:         s += f_typed(a + i * dx)
      ...:     return s * dx
      ...: cpdef np.ndarray[double] apply_integrate_f(np.ndarray col_a, np.ndarray col_b,
      ...:                                            np.ndarray col_N):
      ...:     assert (col_a.dtype == np.float64
      ...:             and col_b.dtype == np.float64 and col_N.dtype == np.dtype(int))
      ...:     cdef Py_ssize_t i, n = len(col_N)
      ...:     assert (len(col_a) == len(col_b) == n)
      ...:     cdef np.ndarray[double] res = np.empty(n)
      ...:     for i in range(len(col_a)):
      ...:         res[i] = integrate_f_typed(col_a[i], col_b[i], col_N[i])
      ...:     return res
      ...:


This implementation creates an array of zeros and inserts the result
of ``integrate_f_typed`` applied over each row. Looping over an ``ndarray`` is faster
in Cython than looping over a :class:`Series` object.

Since ``apply_integrate_f`` is typed to accept an ``np.ndarray``, :meth:`Series.to_numpy`
calls are needed to utilize this function.

.. ipython:: python

   %timeit apply_integrate_f(df["a"].to_numpy(), df["b"].to_numpy(), df["N"].to_numpy())

Performance has improved from the prior implementation by almost ten times.

.. _enhancingperf.boundswrap:

Disabling compiler directives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The majority of the time is now spent in ``apply_integrate_f``. Disabling Cython's ``boundscheck``
and ``wraparound`` checks can yield more performance.

.. ipython:: python

   %prun -l 4 apply_integrate_f(df["a"].to_numpy(), df["b"].to_numpy(), df["N"].to_numpy())

.. ipython::

   In [5]: %%cython
      ...: cimport cython
      ...: cimport numpy as np
      ...: import numpy as np
      ...: cdef np.float64_t f_typed(np.float64_t x) except? -2:
      ...:     return x * (x - 1)
      ...: cpdef np.float64_t integrate_f_typed(np.float64_t a, np.float64_t b, np.int64_t N):
      ...:     cdef np.int64_t i
      ...:     cdef np.float64_t s = 0.0, dx
      ...:     dx = (b - a) / N
      ...:     for i in range(N):
      ...:         s += f_typed(a + i * dx)
      ...:     return s * dx
      ...: @cython.boundscheck(False)
      ...: @cython.wraparound(False)
      ...: cpdef np.ndarray[np.float64_t] apply_integrate_f_wrap(
      ...:     np.ndarray[np.float64_t] col_a,
      ...:     np.ndarray[np.float64_t] col_b,
      ...:     np.ndarray[np.int64_t] col_N
      ...: ):
      ...:     cdef np.int64_t i, n = len(col_N)
      ...:     assert len(col_a) == len(col_b) == n
      ...:     cdef np.ndarray[np.float64_t] res = np.empty(n, dtype=np.float64)
      ...:     for i in range(n):
      ...:         res[i] = integrate_f_typed(col_a[i], col_b[i], col_N[i])
      ...:     return res
      ...:

.. ipython:: python

   %timeit apply_integrate_f_wrap(df["a"].to_numpy(), df["b"].to_numpy(), df["N"].to_numpy())

However, a loop indexer ``i`` accessing an invalid location in an array would cause a segfault because memory access isn't checked.
For more about ``boundscheck`` and ``wraparound``, see the Cython docs on
`compiler directives <https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html#compiler-directives>`__.

.. _enhancingperf.numba:

Numba (JIT compilation)
-----------------------

An alternative to statically compiling Cython code is to use a dynamic just-in-time (JIT) compiler with `Numba <https://numba.pydata.org/>`__.

Numba allows you to write a pure Python function which can be JIT compiled to native machine instructions, similar in performance to C, C++ and Fortran,
by decorating your function with ``@jit``.

Numba works by generating optimized machine code using the LLVM compiler infrastructure at import time, runtime, or statically (using the included pycc tool).
Numba supports compilation of Python to run on either CPU or GPU hardware and is designed to integrate with the Python scientific software stack.

.. note::

    The ``@jit`` compilation will add overhead to the runtime of the function, so performance benefits may not be realized especially when using small data sets.
    Consider `caching <https://numba.readthedocs.io/en/stable/developer/caching.html>`__ your function to avoid compilation overhead each time your function is run.

Numba can be used in 2 ways with pandas:

#. Specify the ``engine="numba"`` keyword in select pandas methods
#. Define your own Python function decorated with ``@jit`` and pass the underlying NumPy array of :class:`Series` or :class:`DataFrame` (using :meth:`Series.to_numpy`) into the function

pandas Numba Engine
~~~~~~~~~~~~~~~~~~~

If Numba is installed, one can specify ``engine="numba"`` in select pandas methods to execute the method using Numba.
Methods that support ``engine="numba"`` will also have an ``engine_kwargs`` keyword that accepts a dictionary that allows one to specify
``"nogil"``, ``"nopython"`` and ``"parallel"`` keys with boolean values to pass into the ``@jit`` decorator.
If ``engine_kwargs`` is not specified, it defaults to ``{"nogil": False, "nopython": True, "parallel": False}`` unless otherwise specified.

.. note::

   In terms of performance, **the first time a function is run using the Numba engine will be slow**
   as Numba will have some function compilation overhead. However, the JIT compiled functions are cached,
   and subsequent calls will be fast. In general, the Numba engine is performant with
   a larger amount of data points (e.g. 1+ million).

   .. code-block:: ipython

      In [1]: data = pd.Series(range(1_000_000))  # noqa: E225

      In [2]: roll = data.rolling(10)

      In [3]: def f(x):
         ...:     return np.sum(x) + 5
      # Run the first time, compilation time will affect performance
      In [4]: %timeit -r 1 -n 1 roll.apply(f, engine='numba', raw=True)
      1.23 s ± 0 ns per loop (mean ± std. dev. of 1 run, 1 loop each)
      # Function is cached and performance will improve
      In [5]: %timeit roll.apply(f, engine='numba', raw=True)
      188 ms ± 1.93 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)

      In [6]: %timeit roll.apply(f, engine='cython', raw=True)
      3.92 s ± 59 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

If your compute hardware contains multiple CPUs, the largest performance gain can be realized by setting ``parallel`` to ``True``
to leverage more than 1 CPU. Internally, pandas leverages numba to parallelize computations over the columns of a :class:`DataFrame`;
therefore, this performance benefit is only beneficial for a :class:`DataFrame` with a large number of columns.

.. code-block:: ipython

   In [1]: import numba

   In [2]: numba.set_num_threads(1)

   In [3]: df = pd.DataFrame(np.random.randn(10_000, 100))

   In [4]: roll = df.rolling(100)

   In [5]: %timeit roll.mean(engine="numba", engine_kwargs={"parallel": True})
   347 ms ± 26 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

   In [6]: numba.set_num_threads(2)

   In [7]: %timeit roll.mean(engine="numba", engine_kwargs={"parallel": True})
   201 ms ± 2.97 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

Custom Function Examples
~~~~~~~~~~~~~~~~~~~~~~~~

A custom Python function decorated with ``@jit`` can be used with pandas objects by passing their NumPy array
representations with :meth:`Series.to_numpy`.

.. code-block:: python

   import numba


   @numba.jit
   def f_plain(x):
       return x * (x - 1)


   @numba.jit
   def integrate_f_numba(a, b, N):
       s = 0
       dx = (b - a) / N
       for i in range(N):
           s += f_plain(a + i * dx)
       return s * dx


   @numba.jit
   def apply_integrate_f_numba(col_a, col_b, col_N):
       n = len(col_N)
       result = np.empty(n, dtype="float64")
       assert len(col_a) == len(col_b) == n
       for i in range(n):
           result[i] = integrate_f_numba(col_a[i], col_b[i], col_N[i])
       return result


   def compute_numba(df):
       result = apply_integrate_f_numba(
           df["a"].to_numpy(), df["b"].to_numpy(), df["N"].to_numpy()
       )
       return pd.Series(result, index=df.index, name="result")


.. code-block:: ipython

   In [4]: %timeit compute_numba(df)
   1000 loops, best of 3: 798 us per loop

In this example, using Numba was faster than Cython.

Numba can also be used to write vectorized functions that do not require the user to explicitly
loop over the observations of a vector; a vectorized function will be applied to each row automatically.
Consider the following example of doubling each observation:

.. code-block:: python

   import numba


   def double_every_value_nonumba(x):
       return x * 2


   @numba.vectorize
   def double_every_value_withnumba(x):  # noqa E501
       return x * 2

.. code-block:: ipython

   # Custom function without numba
   In [5]: %timeit df["col1_doubled"] = df["a"].apply(double_every_value_nonumba)  # noqa E501
   1000 loops, best of 3: 797 us per loop

   # Standard implementation (faster than a custom function)
   In [6]: %timeit df["col1_doubled"] = df["a"] * 2
   1000 loops, best of 3: 233 us per loop

   # Custom function with numba
   In [7]: %timeit df["col1_doubled"] = double_every_value_withnumba(df["a"].to_numpy())
   1000 loops, best of 3: 145 us per loop

Caveats
~~~~~~~

Numba is best at accelerating functions that apply numerical functions to NumPy
arrays. If you try to ``@jit`` a function that contains unsupported `Python <https://numba.readthedocs.io/en/stable/reference/pysupported.html>`__
or `NumPy <https://numba.readthedocs.io/en/stable/reference/numpysupported.html>`__
code, compilation will revert `object mode <https://numba.readthedocs.io/en/stable/glossary.html#term-object-mode>`__ which
will mostly likely not speed up your function. If you would
prefer that Numba throw an error if it cannot compile a function in a way that
speeds up your code, pass Numba the argument
``nopython=True`` (e.g.  ``@jit(nopython=True)``). For more on
troubleshooting Numba modes, see the `Numba troubleshooting page
<https://numba.pydata.org/numba-doc/latest/user/troubleshoot.html#the-compiled-code-is-too-slow>`__.

Using ``parallel=True`` (e.g. ``@jit(parallel=True)``) may result in a ``SIGABRT`` if the threading layer leads to unsafe
behavior. You can first `specify a safe threading layer <https://numba.readthedocs.io/en/stable/user/threading-layer.html#selecting-a-threading-layer-for-safe-parallel-execution>`__
before running a JIT function with ``parallel=True``.

Generally if the you encounter a segfault (``SIGSEGV``) while using Numba, please report the issue
to the `Numba issue tracker. <https://github.com/numba/numba/issues/new/choose>`__

.. _enhancingperf.eval:

Expression evaluation via :func:`~pandas.eval`
----------------------------------------------

The top-level function :func:`pandas.eval` implements performant expression evaluation of
:class:`~pandas.Series` and :class:`~pandas.DataFrame`. Expression evaluation allows operations
to be expressed as strings and can potentially provide a performance improvement
by evaluate arithmetic and boolean expression all at once for large :class:`~pandas.DataFrame`.

.. note::

   You should not use :func:`~pandas.eval` for simple
   expressions or for expressions involving small DataFrames. In fact,
   :func:`~pandas.eval` is many orders of magnitude slower for
   smaller expressions or objects than plain Python. A good rule of thumb is
   to only use :func:`~pandas.eval` when you have a
   :class:`.DataFrame` with more than 10,000 rows.

Supported syntax
~~~~~~~~~~~~~~~~

These operations are supported by :func:`pandas.eval`:

* Arithmetic operations except for the left shift (``<<``) and right shift
  (``>>``) operators, e.g., ``df + 2 * pi / s ** 4 % 42 - the_golden_ratio``
* Comparison operations, including chained comparisons, e.g., ``2 < df < df2``
* Boolean operations, e.g., ``df < df2 and df3 < df4 or not df_bool``
* ``list`` and ``tuple`` literals, e.g., ``[1, 2]`` or ``(1, 2)``
* Attribute access, e.g., ``df.a``
* Subscript expressions, e.g., ``df[0]``
* Simple variable evaluation, e.g., ``pd.eval("df")`` (this is not very useful)
* Math functions: ``sin``, ``cos``, ``exp``, ``log``, ``expm1``, ``log1p``,
  ``sqrt``, ``sinh``, ``cosh``, ``tanh``, ``arcsin``, ``arccos``, ``arctan``, ``arccosh``,
  ``arcsinh``, ``arctanh``, ``abs``, ``arctan2`` and ``log10``.

The following Python syntax is **not** allowed:

* Expressions

    * Function calls other than math functions.
    * ``is``/``is not`` operations
    * ``if`` expressions
    * ``lambda`` expressions
    * ``list``/``set``/``dict`` comprehensions
    * Literal ``dict`` and ``set`` expressions
    * ``yield`` expressions
    * Generator expressions
    * Boolean expressions consisting of only scalar values

* Statements

    * Neither `simple <https://docs.python.org/3/reference/simple_stmts.html>`__
      or `compound <https://docs.python.org/3/reference/compound_stmts.html>`__
      statements are allowed. This includes ``for``, ``while``, and
      ``if``.

Local variables
~~~~~~~~~~~~~~~

You must *explicitly reference* any local variable that you want to use in an
expression by placing the ``@`` character in front of the name. This mechanism is
the same for both :meth:`DataFrame.query` and :meth:`DataFrame.eval`. For example,

.. ipython:: python

   df = pd.DataFrame(np.random.randn(5, 2), columns=list("ab"))
   newcol = np.random.randn(len(df))
   df.eval("b + @newcol")
   df.query("b < @newcol")

If you don't prefix the local variable with ``@``, pandas will raise an
exception telling you the variable is undefined.

When using :meth:`DataFrame.eval` and :meth:`DataFrame.query`, this allows you
to have a local variable and a :class:`~pandas.DataFrame` column with the same
name in an expression.


.. ipython:: python

   a = np.random.randn()
   df.query("@a < a")
   df.loc[a < df["a"]]  # same as the previous expression

.. warning::

   :func:`pandas.eval` will raise an exception if you cannot use the ``@`` prefix because it
   isn't defined in that context.

   .. ipython:: python
      :okexcept:

      a, b = 1, 2
      pd.eval("@a + b")

   In this case, you should simply refer to the variables like you would in
   standard Python.

   .. ipython:: python

      pd.eval("a + b")


:func:`pandas.eval` parsers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two different expression syntax parsers.

The default ``'pandas'`` parser allows a more intuitive syntax for expressing
query-like operations (comparisons, conjunctions and disjunctions). In
particular, the precedence of the ``&`` and ``|`` operators is made equal to
the precedence of the corresponding boolean operations ``and`` and ``or``.

For example, the above conjunction can be written without parentheses.
Alternatively, you can use the ``'python'`` parser to enforce strict Python
semantics.

.. ipython:: python

   nrows, ncols = 20000, 100
   df1, df2, df3, df4 = [pd.DataFrame(np.random.randn(nrows, ncols)) for _ in range(4)]

   expr = "(df1 > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)"
   x = pd.eval(expr, parser="python")
   expr_no_parens = "df1 > 0 & df2 > 0 & df3 > 0 & df4 > 0"
   y = pd.eval(expr_no_parens, parser="pandas")
   np.all(x == y)


The same expression can be "anded" together with the word :keyword:`and` as
well:

.. ipython:: python

   expr = "(df1 > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)"
   x = pd.eval(expr, parser="python")
   expr_with_ands = "df1 > 0 and df2 > 0 and df3 > 0 and df4 > 0"
   y = pd.eval(expr_with_ands, parser="pandas")
   np.all(x == y)

The :keyword:`and` and :keyword:`or` operators here have the same precedence that they would
in Python.


:func:`pandas.eval` engines
~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two different expression engines.

The ``'numexpr'`` engine is the more performant engine that can yield performance improvements
compared to standard Python syntax for large :class:`DataFrame`. This engine requires the
optional dependency ``numexpr`` to be installed.

The ``'python'`` engine is generally *not* useful except for testing
other evaluation engines against it. You will achieve **no** performance
benefits using :func:`~pandas.eval` with ``engine='python'`` and may
incur a performance hit.

.. ipython:: python

   %timeit df1 + df2 + df3 + df4

.. ipython:: python

   %timeit pd.eval("df1 + df2 + df3 + df4", engine="python")


The :meth:`DataFrame.eval` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to the top level :func:`pandas.eval` function you can also
evaluate an expression in the "context" of a :class:`~pandas.DataFrame`.

.. ipython:: python
   :suppress:

   try:
       del a
   except NameError:
       pass

   try:
       del b
   except NameError:
       pass

.. ipython:: python

   df = pd.DataFrame(np.random.randn(5, 2), columns=["a", "b"])
   df.eval("a + b")

Any expression that is a valid :func:`pandas.eval` expression is also a valid
:meth:`DataFrame.eval` expression, with the added benefit that you don't have to
prefix the name of the :class:`~pandas.DataFrame` to the column(s) you're
interested in evaluating.

In addition, you can perform assignment of columns within an expression.
This allows for *formulaic evaluation*. The assignment target can be a
new column name or an existing column name, and it must be a valid Python
identifier.

.. ipython:: python

   df = pd.DataFrame(dict(a=range(5), b=range(5, 10)))
   df = df.eval("c = a + b")
   df = df.eval("d = a + b + c")
   df = df.eval("a = 1")
   df

A copy of the :class:`DataFrame` with the
new or modified columns is returned, and the original frame is unchanged.

.. ipython:: python

   df
   df.eval("e = a - c")
   df

Multiple column assignments can be performed by using a multi-line string.

.. ipython:: python

   df.eval(
       """
   c = a + b
   d = a + b + c
   a = 1""",
   )

The equivalent in standard Python would be

.. ipython:: python

   df = pd.DataFrame(dict(a=range(5), b=range(5, 10)))
   df["c"] = df["a"] + df["b"]
   df["d"] = df["a"] + df["b"] + df["c"]
   df["a"] = 1
   df


:func:`~pandas.eval` performance comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`pandas.eval` works well with expressions containing large arrays.

.. ipython:: python

   nrows, ncols = 20000, 100
   df1, df2, df3, df4 = [pd.DataFrame(np.random.randn(nrows, ncols)) for _ in range(4)]


:class:`DataFrame` arithmetic:

.. ipython:: python

   %timeit df1 + df2 + df3 + df4

.. ipython:: python

   %timeit pd.eval("df1 + df2 + df3 + df4")


:class:`DataFrame` comparison:

.. ipython:: python

   %timeit (df1 > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)

.. ipython:: python

   %timeit pd.eval("(df1 > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)")


:class:`DataFrame` arithmetic with unaligned axes.

.. ipython:: python

   s = pd.Series(np.random.randn(50))
   %timeit df1 + df2 + df3 + df4 + s

.. ipython:: python

   %timeit pd.eval("df1 + df2 + df3 + df4 + s")

.. note::

   Operations such as

   .. code-block:: python

      1 and 2  # would parse to 1 & 2, but should evaluate to 2
      3 or 4  # would parse to 3 | 4, but should evaluate to 3
      ~1  # this is okay, but slower when using eval

   should be performed in Python. An exception will be raised if you try to
   perform any boolean/bitwise operations with scalar operands that are not
   of type ``bool`` or ``np.bool_``.

Here is a plot showing the running time of
:func:`pandas.eval` as function of the size of the frame involved in the
computation. The two lines are two different engines.

..
    The eval-perf.png figure below was generated with /doc/scripts/eval_performance.py

.. image:: ../_static/eval-perf.png

You will only see the performance benefits of using the ``numexpr`` engine with :func:`pandas.eval` if your :class:`~pandas.DataFrame`
has more than approximately 100,000 rows.

This plot was created using a :class:`DataFrame` with 3 columns each containing
floating point values generated using ``numpy.random.randn()``.

Expression evaluation limitations with ``numexpr``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Expressions that would result in an object dtype or involve datetime operations
because of ``NaT`` must be evaluated in Python space, but part of an expression
can still be evaluated with ``numexpr``. For example:

.. ipython:: python

   df = pd.DataFrame(
       {"strings": np.repeat(list("cba"), 3), "nums": np.repeat(range(3), 3)}
   )
   df
   df.query("strings == 'a' and nums == 1")

The numeric part of the comparison (``nums == 1``) will be evaluated by
``numexpr`` and the object part of the comparison (``"strings == 'a'``) will
be evaluated by Python.
