.. _enhancingperf:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   np.set_printoptions(precision=4, suppress=True)
   import pandas as pd
   pd.options.display.max_rows=15

   import os
   import csv


*********************
Enhancing Performance
*********************

.. _enhancingperf.cython:

Cython (Writing C extensions for pandas)
----------------------------------------

For many use cases writing pandas in pure python and numpy is sufficient. In some
computationally heavy applications however, it can be possible to achieve sizeable
speed-ups by offloading work to `cython <http://cython.org/>`__.

This tutorial assumes you have refactored as much as possible in python, for example
trying to remove for loops and making use of numpy vectorization, it's always worth
optimising in python first.

This tutorial walks through a "typical" process of cythonizing a slow computation.
We use an `example from the cython documentation <http://docs.cython.org/src/quickstart/cythonize.html>`__
but in the context of pandas. Our final cythonized solution is around 100 times
faster than the pure python.

.. _enhancingperf.pure:

Pure python
~~~~~~~~~~~

We have a DataFrame to which we want to apply a function row-wise.

.. ipython:: python

   df = pd.DataFrame({'a': np.random.randn(1000),
                      'b': np.random.randn(1000),
                      'N': np.random.randint(100, 1000, (1000)),
                      'x': 'x'})
   df

Here's the function in pure python:

.. ipython:: python

   def f(x):
       return x * (x - 1)

   def integrate_f(a, b, N):
       s = 0
       dx = (b - a) / N
       for i in range(N):
           s += f(a + i * dx)
       return s * dx

We achieve our result by using ``apply`` (row-wise):

.. code-block:: ipython

   In [7]: %timeit df.apply(lambda x: integrate_f(x['a'], x['b'], x['N']), axis=1)
   10 loops, best of 3: 174 ms per loop

But clearly this isn't fast enough for us. Let's take a look and see where the
time is spent during this operation (limited to the most time consuming
four calls) using the `prun ipython magic function <http://ipython.org/ipython-doc/stable/api/generated/IPython.core.magics.execution.html#IPython.core.magics.execution.ExecutionMagics.prun>`__:

.. ipython:: python

   %prun -l 4 df.apply(lambda x: integrate_f(x['a'], x['b'], x['N']), axis=1)

By far the majority of time is spend inside either ``integrate_f`` or ``f``,
hence we'll concentrate our efforts cythonizing these two functions.

.. note::

  In python 2 replacing the ``range`` with its generator counterpart (``xrange``)
  would mean the ``range`` line would vanish. In python 3 ``range`` is already a generator.

.. _enhancingperf.plain:

Plain cython
~~~~~~~~~~~~

First we're going to need to import the cython magic function to ipython (for
cython versions  < 0.21 you can use ``%load_ext cythonmagic``):

.. ipython:: python
   :okwarning:

   %load_ext Cython


Now, let's simply copy our functions over to cython as is (the suffix
is here to distinguish between function versions):

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

.. note::

  If you're having trouble pasting the above into your ipython, you may need
  to be using bleeding edge ipython for paste to play well with cell magics.


.. code-block:: ipython

   In [4]: %timeit df.apply(lambda x: integrate_f_plain(x['a'], x['b'], x['N']), axis=1)
   10 loops, best of 3: 85.5 ms per loop

Already this has shaved a third off, not too bad for a simple copy and paste.

.. _enhancingperf.type:

Adding type
~~~~~~~~~~~

We get another huge improvement simply by providing type information:

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

.. code-block:: ipython

   In [4]: %timeit df.apply(lambda x: integrate_f_typed(x['a'], x['b'], x['N']), axis=1)
   10 loops, best of 3: 20.3 ms per loop

Now, we're talking! It's now over ten times faster than the original python
implementation, and we haven't *really* modified the code. Let's have another
look at what's eating up time:

.. ipython:: python

   %prun -l 4 df.apply(lambda x: integrate_f_typed(x['a'], x['b'], x['N']), axis=1)

.. _enhancingperf.ndarray:

Using ndarray
~~~~~~~~~~~~~

It's calling series... a lot! It's creating a Series from each row, and get-ting from both
the index and the series (three times for each row). Function calls are expensive
in python, so maybe we could minimise these by cythonizing the apply part.

.. note::

  We are now passing ndarrays into the cython function, fortunately cython plays
  very nicely with numpy.

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
      ...: cpdef np.ndarray[double] apply_integrate_f(np.ndarray col_a, np.ndarray col_b, np.ndarray col_N):
      ...:     assert (col_a.dtype == np.float and col_b.dtype == np.float and col_N.dtype == np.int)
      ...:     cdef Py_ssize_t i, n = len(col_N)
      ...:     assert (len(col_a) == len(col_b) == n)
      ...:     cdef np.ndarray[double] res = np.empty(n)
      ...:     for i in range(len(col_a)):
      ...:         res[i] = integrate_f_typed(col_a[i], col_b[i], col_N[i])
      ...:     return res
      ...:


The implementation is simple, it creates an array of zeros and loops over
the rows, applying our ``integrate_f_typed``, and putting this in the zeros array.


.. warning::

   In 0.13.0 since ``Series`` has internaly been refactored to no longer sub-class ``ndarray``
   but instead subclass ``NDFrame``, you can **not pass** a ``Series`` directly as a ``ndarray`` typed parameter
   to a cython function. Instead pass the actual ``ndarray`` using the ``.values`` attribute of the Series.

   Prior to 0.13.0

   .. code-block:: python

        apply_integrate_f(df['a'], df['b'], df['N'])

   Use ``.values`` to get the underlying ``ndarray``

   .. code-block:: python

        apply_integrate_f(df['a'].values, df['b'].values, df['N'].values)

.. note::

    Loops like this would be *extremely* slow in python, but in Cython looping
    over numpy arrays is *fast*.

.. code-block:: ipython

   In [4]: %timeit apply_integrate_f(df['a'].values, df['b'].values, df['N'].values)
   1000 loops, best of 3: 1.25 ms per loop

We've gotten another big improvement. Let's check again where the time is spent:

.. ipython:: python

   %prun -l 4 apply_integrate_f(df['a'].values, df['b'].values, df['N'].values)

As one might expect, the majority of the time is now spent in ``apply_integrate_f``,
so if we wanted to make anymore efficiencies we must continue to concentrate our
efforts here.

.. _enhancingperf.boundswrap:

More advanced techniques
~~~~~~~~~~~~~~~~~~~~~~~~

There is still hope for improvement. Here's an example of using some more
advanced cython techniques:

.. ipython::

   In [5]: %%cython
      ...: cimport cython
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
      ...: @cython.boundscheck(False)
      ...: @cython.wraparound(False)
      ...: cpdef np.ndarray[double] apply_integrate_f_wrap(np.ndarray[double] col_a, np.ndarray[double] col_b, np.ndarray[int] col_N):
      ...:     cdef int i, n = len(col_N)
      ...:     assert len(col_a) == len(col_b) == n
      ...:     cdef np.ndarray[double] res = np.empty(n)
      ...:     for i in range(n):
      ...:         res[i] = integrate_f_typed(col_a[i], col_b[i], col_N[i])
      ...:     return res
      ...:

.. code-block:: ipython

   In [4]: %timeit apply_integrate_f_wrap(df['a'].values, df['b'].values, df['N'].values)
   1000 loops, best of 3: 987 us per loop

Even faster, with the caveat that a bug in our cython code (an off-by-one error,
for example) might cause a segfault because memory access isn't checked.


.. _enhancingperf.numba:

Using numba
-----------

A recent alternative to statically compiling cython code, is to use a *dynamic jit-compiler*, ``numba``.

Numba gives you the power to speed up your applications with high performance functions written directly in Python. With a few annotations, array-oriented and math-heavy Python code can be just-in-time compiled to native machine instructions, similar in performance to C, C++ and Fortran, without having to switch languages or Python interpreters.

Numba works by generating optimized machine code using the LLVM compiler infrastructure at import time, runtime, or statically (using the included pycc tool). Numba supports compilation of Python to run on either CPU or GPU hardware, and is designed to integrate with the Python scientific software stack.

.. note::

    You will need to install ``numba``. This is easy with ``conda``, by using: ``conda install numba``, see :ref:`installing using miniconda<install.miniconda>`.

.. note::

    As of ``numba`` version 0.20, pandas objects cannot be passed directly to numba-compiled functions. Instead, one must pass the ``numpy`` array underlying the ``pandas`` object to the numba-compiled function as demonstrated below.

Jit
~~~

Using ``numba`` to just-in-time compile your code. We simply take the plain python code from above and annotate with the ``@jit`` decorator.

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
       result = np.empty(n, dtype='float64')
       assert len(col_a) == len(col_b) == n
       for i in range(n):
          result[i] = integrate_f_numba(col_a[i], col_b[i], col_N[i])
       return result

    def compute_numba(df):
       result = apply_integrate_f_numba(df['a'].values, df['b'].values, df['N'].values)
       return pd.Series(result, index=df.index, name='result')

Note that we directly pass ``numpy`` arrays to the numba function. ``compute_numba`` is just a wrapper that provides a nicer interface by passing/returning pandas objects.

.. code-block:: ipython

    In [4]: %timeit compute_numba(df)
    1000 loops, best of 3: 798 us per loop

Vectorize
~~~~~~~~~

``numba`` can also be used to write vectorized functions that do not require the user to explicitly
loop over the observations of a vector; a vectorized function will be applied to each row automatically.
Consider the following toy example of doubling each observation:

.. code-block:: python

    import numba

    def double_every_value_nonumba(x):
        return x*2

    @numba.vectorize
    def double_every_value_withnumba(x):
        return x*2


    # Custom function without numba
    In [5]: %timeit df['col1_doubled'] = df.a.apply(double_every_value_nonumba)
    1000 loops, best of 3: 797 us per loop

    # Standard implementation (faster than a custom function)
    In [6]: %timeit df['col1_doubled'] = df.a*2
    1000 loops, best of 3: 233 us per loop

    # Custom function with numba
    In [7]: %timeit df['col1_doubled'] = double_every_value_withnumba(df.a.values)
    1000 loops, best of 3: 145 us per loop

Caveats
~~~~~~~

.. note::

    ``numba`` will execute on any function, but can only accelerate certain classes of functions.

``numba`` is best at accelerating functions that apply numerical functions to numpy arrays. When passed a function that only uses operations it knows how to accelerate, it will execute in ``nopython`` mode.

If ``numba`` is passed a function that includes something it doesn't know how to work with -- a category that currently includes sets, lists, dictionaries, or string functions -- it will revert to ``object mode``. In ``object mode``, numba will execute but your code will not speed up significantly. If you would prefer that ``numba`` throw an error if it cannot compile a function in a way that speeds up your code, pass numba the argument ``nopython=True`` (e.g.  ``@numba.jit(nopython=True)``). For more on troubleshooting ``numba`` modes, see the `numba troubleshooting page <http://numba.pydata.org/numba-doc/0.20.0/user/troubleshoot.html#the-compiled-code-is-too-slow>`__.

Read more in the `numba docs <http://numba.pydata.org/>`__.

.. _enhancingperf.eval:

Expression Evaluation via :func:`~pandas.eval` (Experimental)
-------------------------------------------------------------

.. versionadded:: 0.13

The top-level function :func:`pandas.eval` implements expression evaluation of
:class:`~pandas.Series` and :class:`~pandas.DataFrame` objects.

.. note::

   To benefit from using :func:`~pandas.eval` you need to
   install ``numexpr``. See the :ref:`recommended dependencies section
   <install.recommended_dependencies>` for more details.

The point of using :func:`~pandas.eval` for expression evaluation rather than
plain Python is two-fold: 1) large :class:`~pandas.DataFrame` objects are
evaluated more efficiently and 2) large arithmetic and boolean expressions are
evaluated all at once by the underlying engine (by default ``numexpr`` is used
for evaluation).

.. note::

   You should not use :func:`~pandas.eval` for simple
   expressions or for expressions involving small DataFrames. In fact,
   :func:`~pandas.eval` is many orders of magnitude slower for
   smaller expressions/objects than plain ol' Python. A good rule of thumb is
   to only use :func:`~pandas.eval` when you have a
   :class:`~pandas.core.frame.DataFrame` with more than 10,000 rows.


:func:`~pandas.eval` supports all arithmetic expressions supported by the
engine in addition to some extensions available only in pandas.

.. note::

   The larger the frame and the larger the expression the more speedup you will
   see from using :func:`~pandas.eval`.

Supported Syntax
~~~~~~~~~~~~~~~~

These operations are supported by :func:`pandas.eval`:

- Arithmetic operations except for the left shift (``<<``) and right shift
  (``>>``) operators, e.g., ``df + 2 * pi / s ** 4 % 42 - the_golden_ratio``
- Comparison operations, including chained comparisons, e.g., ``2 < df < df2``
- Boolean operations, e.g., ``df < df2 and df3 < df4 or not df_bool``
- ``list`` and ``tuple`` literals, e.g., ``[1, 2]`` or ``(1, 2)``
- Attribute access, e.g., ``df.a``
- Subscript expressions, e.g., ``df[0]``
- Simple variable evaluation, e.g., ``pd.eval('df')`` (this is not very useful)
- Math functions, `sin`, `cos`, `exp`, `log`, `expm1`, `log1p`,
  `sqrt`, `sinh`, `cosh`, `tanh`, `arcsin`, `arccos`, `arctan`, `arccosh`,
  `arcsinh`, `arctanh`, `abs` and `arctan2`.

This Python syntax is **not** allowed:

* Expressions

  - Function calls other than math functions.
  - ``is``/``is not`` operations
  - ``if`` expressions
  - ``lambda`` expressions
  - ``list``/``set``/``dict`` comprehensions
  - Literal ``dict`` and ``set`` expressions
  - ``yield`` expressions
  - Generator expressions
  - Boolean expressions consisting of only scalar values

* Statements

  - Neither `simple <http://docs.python.org/2/reference/simple_stmts.html>`__
    nor `compound <http://docs.python.org/2/reference/compound_stmts.html>`__
    statements are allowed. This includes things like ``for``, ``while``, and
    ``if``.



:func:`~pandas.eval` Examples
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`pandas.eval` works well with expressions containing large arrays.

First let's create a few decent-sized arrays to play with:

.. ipython:: python

   nrows, ncols = 20000, 100
   df1, df2, df3, df4 = [pd.DataFrame(np.random.randn(nrows, ncols)) for _ in range(4)]


Now let's compare adding them together using plain ol' Python versus
:func:`~pandas.eval`:

.. ipython:: python

   %timeit df1 + df2 + df3 + df4

.. ipython:: python

   %timeit pd.eval('df1 + df2 + df3 + df4')


Now let's do the same thing but with comparisons:

.. ipython:: python

   %timeit (df1 > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)

.. ipython:: python

   %timeit pd.eval('(df1 > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)')


:func:`~pandas.eval` also works with unaligned pandas objects:

.. ipython:: python

   s = pd.Series(np.random.randn(50))
   %timeit df1 + df2 + df3 + df4 + s

.. ipython:: python

   %timeit pd.eval('df1 + df2 + df3 + df4 + s')

.. note::

   Operations such as

      .. code-block:: python

         1 and 2  # would parse to 1 & 2, but should evaluate to 2
         3 or 4  # would parse to 3 | 4, but should evaluate to 3
         ~1  # this is okay, but slower when using eval

   should be performed in Python. An exception will be raised if you try to
   perform any boolean/bitwise operations with scalar operands that are not
   of type ``bool`` or ``np.bool_``. Again, you should perform these kinds of
   operations in plain Python.

The ``DataFrame.eval`` method (Experimental)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.13

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

   df = pd.DataFrame(np.random.randn(5, 2), columns=['a', 'b'])
   df.eval('a + b')

Any expression that is a valid :func:`pandas.eval` expression is also a valid
:meth:`DataFrame.eval` expression, with the added benefit that you don't have to
prefix the name of the :class:`~pandas.DataFrame` to the column(s) you're
interested in evaluating.

In addition, you can perform assignment of columns within an expression.
This allows for *formulaic evaluation*.  The assignment target can be a
new column name or an existing column name, and it must be a valid Python
identifier.

.. versionadded:: 0.18.0

The ``inplace`` keyword determines whether this assignment will performed
on the original ``DataFrame`` or return a copy with the new column.

.. warning::

   For backwards compatability, ``inplace`` defaults to ``True`` if not
   specified. This will change in a future version of pandas - if your
   code depends on an inplace assignment you should update to explicitly
   set ``inplace=True``

.. ipython:: python

   df = pd.DataFrame(dict(a=range(5), b=range(5, 10)))
   df.eval('c = a + b', inplace=True)
   df.eval('d = a + b + c', inplace=True)
   df.eval('a = 1', inplace=True)
   df

When ``inplace`` is set to ``False``, a copy of the ``DataFrame`` with the
new or modified columns is returned and the original frame is unchanged.

.. ipython:: python

   df
   df.eval('e = a - c', inplace=False)
   df

.. versionadded:: 0.18.0

As a convenience, multiple assignments can be performed by using a
multi-line string.

.. ipython:: python

   df.eval("""
   c = a + b
   d = a + b + c
   a = 1""", inplace=False)

The equivalent in standard Python would be

.. ipython:: python

   df = pd.DataFrame(dict(a=range(5), b=range(5, 10)))
   df['c'] = df.a + df.b
   df['d'] = df.a + df.b + df.c
   df['a'] = 1
   df

.. versionadded:: 0.18.0

The ``query`` method gained the ``inplace`` keyword which determines
whether the query modifies the original frame.

.. ipython:: python

   df = pd.DataFrame(dict(a=range(5), b=range(5, 10)))
   df.query('a > 2')
   df.query('a > 2', inplace=True)
   df

.. warning::

   Unlike with ``eval``, the default value for ``inplace`` for ``query``
   is ``False``.  This is consistent with prior versions of pandas.

Local Variables
~~~~~~~~~~~~~~~

In pandas version 0.14 the local variable API has changed. In pandas 0.13.x,
you could refer to local variables the same way you would in standard Python.
For example,

.. code-block:: python

   df = pd.DataFrame(np.random.randn(5, 2), columns=['a', 'b'])
   newcol = np.random.randn(len(df))
   df.eval('b + newcol')

   UndefinedVariableError: name 'newcol' is not defined

As you can see from the exception generated, this syntax is no longer allowed.
You must *explicitly reference* any local variable that you want to use in an
expression by placing the ``@`` character in front of the name. For example,

.. ipython:: python

   df = pd.DataFrame(np.random.randn(5, 2), columns=list('ab'))
   newcol = np.random.randn(len(df))
   df.eval('b + @newcol')
   df.query('b < @newcol')

If you don't prefix the local variable with ``@``, pandas will raise an
exception telling you the variable is undefined.

When using :meth:`DataFrame.eval` and :meth:`DataFrame.query`, this allows you
to have a local variable and a :class:`~pandas.DataFrame` column with the same
name in an expression.


.. ipython:: python

   a = np.random.randn()
   df.query('@a < a')
   df.loc[a < df.a]  # same as the previous expression

With :func:`pandas.eval` you cannot use the ``@`` prefix *at all*, because it
isn't defined in that context. ``pandas`` will let you know this if you try to
use ``@`` in a top-level call to :func:`pandas.eval`. For example,

.. ipython:: python
   :okexcept:

   a, b = 1, 2
   pd.eval('@a + b')

In this case, you should simply refer to the variables like you would in
standard Python.

.. ipython:: python

   pd.eval('a + b')


:func:`pandas.eval` Parsers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two different parsers and two different engines you can use as
the backend.

The default ``'pandas'`` parser allows a more intuitive syntax for expressing
query-like operations (comparisons, conjunctions and disjunctions). In
particular, the precedence of the ``&`` and ``|`` operators is made equal to
the precedence of the corresponding boolean operations ``and`` and ``or``.

For example, the above conjunction can be written without parentheses.
Alternatively, you can use the ``'python'`` parser to enforce strict Python
semantics.

.. ipython:: python

   expr = '(df1 > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)'
   x = pd.eval(expr, parser='python')
   expr_no_parens = 'df1 > 0 & df2 > 0 & df3 > 0 & df4 > 0'
   y = pd.eval(expr_no_parens, parser='pandas')
   np.all(x == y)


The same expression can be "anded" together with the word :keyword:`and` as
well:

.. ipython:: python

   expr = '(df1 > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)'
   x = pd.eval(expr, parser='python')
   expr_with_ands = 'df1 > 0 and df2 > 0 and df3 > 0 and df4 > 0'
   y = pd.eval(expr_with_ands, parser='pandas')
   np.all(x == y)


The ``and`` and ``or`` operators here have the same precedence that they would
in vanilla Python.


:func:`pandas.eval` Backends
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There's also the option to make :func:`~pandas.eval` operate identical to plain
ol' Python.

.. note::

   Using the ``'python'`` engine is generally *not* useful, except for testing
   other evaluation engines against it. You will achieve **no** performance
   benefits using :func:`~pandas.eval` with ``engine='python'`` and in fact may
   incur a performance hit.

You can see this by using :func:`pandas.eval` with the ``'python'`` engine. It
is a bit slower (not by much) than evaluating the same expression in Python

.. ipython:: python

   %timeit df1 + df2 + df3 + df4

.. ipython:: python

   %timeit pd.eval('df1 + df2 + df3 + df4', engine='python')


:func:`pandas.eval` Performance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`~pandas.eval` is intended to speed up certain kinds of operations. In
particular, those operations involving complex expressions with large
:class:`~pandas.DataFrame`/:class:`~pandas.Series` objects should see a
significant performance benefit.  Here is a plot showing the running time of
:func:`pandas.eval` as function of the size of the frame involved in the
computation. The two lines are two different engines.


.. image:: _static/eval-perf.png


.. note::

   Operations with smallish objects (around 15k-20k rows) are faster using
   plain Python:

       .. image:: _static/eval-perf-small.png


This plot was created using a ``DataFrame`` with 3 columns each containing
floating point values generated using ``numpy.random.randn()``.

Technical Minutia Regarding Expression Evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Expressions that would result in an object dtype or involve datetime operations
(because of ``NaT``) must be evaluated in Python space. The main reason for
this behavior is to maintain backwards compatibility with versions of numpy <
1.7. In those versions of ``numpy`` a call to ``ndarray.astype(str)`` will
truncate any strings that are more than 60 characters in length. Second, we
can't pass ``object`` arrays to ``numexpr`` thus string comparisons must be
evaluated in Python space.

The upshot is that this *only* applies to object-dtype'd expressions. So, if
you have an expression--for example

.. ipython:: python

   df = pd.DataFrame({'strings': np.repeat(list('cba'), 3),
                      'nums': np.repeat(range(3), 3)})
   df
   df.query('strings == "a" and nums == 1')

the numeric part of the comparison (``nums == 1``) will be evaluated by
``numexpr``.

In general, :meth:`DataFrame.query`/:func:`pandas.eval` will
evaluate the subexpressions that *can* be evaluated by ``numexpr`` and those
that must be evaluated in Python space transparently to the user. This is done
by inferring the result type of an expression from its arguments and operators.
