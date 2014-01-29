.. _enhancingperf:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import os
   import csv
   from pandas import DataFrame
   import pandas as pd
   pd.options.display.max_rows=15

   import numpy as np
   np.random.seed(123456)
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)


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

   df = DataFrame({'a': randn(1000), 'b': randn(1000),'N': randint(100, 1000, (1000)), 'x': 'x'})
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

We achieve our result by by using ``apply`` (row-wise):

.. ipython:: python

   %timeit df.apply(lambda x: integrate_f(x['a'], x['b'], x['N']), axis=1)

But clearly this isn't fast enough for us. Let's take a look and see where the
time is spent during this operation (limited to the most time consuming
four calls) using the `prun ipython magic function <http://ipython.org/ipython-doc/stable/api/generated/IPython.core.magics.execution.html#IPython.core.magics.execution.ExecutionMagics.prun>`__:

.. ipython:: python

   %prun -l 4 df.apply(lambda x: integrate_f(x['a'], x['b'], x['N']), axis=1)

By far the majority of time is spend inside either ``integrate_f`` or ``f``,
hence we'll concentrate our efforts cythonizing these two functions.

.. note::

  In python 2 replacing the ``range`` with its generator counterpart (``xrange``)
  would mean the ``range`` line would vanish. In python 3 range is already a generator.

.. _enhancingperf.plain:

Plain cython
~~~~~~~~~~~~

First we're going to need to import the cython magic function to ipython:

.. ipython:: python

   %load_ext cythonmagic


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


.. ipython:: python

   %timeit df.apply(lambda x: integrate_f_plain(x['a'], x['b'], x['N']), axis=1)

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

.. ipython:: python

   %timeit df.apply(lambda x: integrate_f_typed(x['a'], x['b'], x['N']), axis=1)

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

.. ipython:: python

   %timeit apply_integrate_f(df['a'].values, df['b'].values, df['N'].values)

We've gotten another big improvement. Let's check again where the time is spent:

.. ipython:: python

   %prun -l 4 apply_integrate_f(df['a'].values, df['b'].values, df['N'].values)

As one might expect, the majority of the time is now spent in ``apply_integrate_f``,
so if we wanted to make anymore efficiencies we must continue to concentrate our
efforts here.

.. _enhancingperf.boundswrap:

More advanced techniques
~~~~~~~~~~~~~~~~~~~~~~~~

There is still scope for improvement, here's an example of using some more
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
      ...: cpdef np.ndarray[double] apply_integrate_f_wrap(np.ndarray[double] col_a, np.ndarray[double] col_b, np.ndarray[Py_ssize_t] col_N):
      ...:     cdef Py_ssize_t i, n = len(col_N)
      ...:     assert len(col_a) == len(col_b) == n
      ...:     cdef np.ndarray[double] res = np.empty(n)
      ...:     for i in range(n):
      ...:         res[i] = integrate_f_typed(col_a[i], col_b[i], col_N[i])
      ...:     return res
      ...:

.. ipython:: python

   %timeit apply_integrate_f_wrap(df['a'].values, df['b'].values, df['N'].values)

Even faster, with the caveat that a bug in our cython code (an off-by-one error,
for example) might cause a segfault because memory access isn't checked.


Further topics
~~~~~~~~~~~~~~

- Loading C modules into cython.

Read more in the `cython docs <http://docs.cython.org/>`__.

.. _enhancingperf.eval:

Expression Evaluation via :func:`~pandas.eval` (Experimental)
-------------------------------------------------------------

.. versionadded:: 0.13

The top-level function :func:`~pandas.eval` implements expression evaluation of
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

These operations are supported by :func:`~pandas.eval`:

- Arithmetic operations except for the left shift (``<<``) and right shift
  (``>>``) operators, e.g., ``df + 2 * pi / s ** 4 % 42 - the_golden_ratio``
- Comparison operations, e.g., ``2 < df < df2``
- Boolean operations, e.g., ``df < df2 and df3 < df4 or not df_bool``
- ``list`` and ``tuple`` literals, e.g., ``[1, 2]`` or ``(1, 2)``
- Attribute access, e.g., ``df.a``
- Subscript expressions, e.g., ``df[0]``
- Simple variable evaluation, e.g., ``pd.eval('df')`` (this is not very useful)

This Python syntax is **not** allowed:

* Expressions

  - Function calls
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

:func:`~pandas.eval` works wonders for expressions containing large arrays

First let's create 4 decent-sized arrays to play with:

.. ipython:: python

   import pandas as pd
   from pandas import DataFrame, Series
   from numpy.random import randn
   import numpy as np
   nrows, ncols = 20000, 100
   df1, df2, df3, df4 = [DataFrame(randn(nrows, ncols)) for _ in xrange(4)]


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

   s = Series(randn(50))
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

In addition to the top level :func:`~pandas.eval` function you can also
evaluate an expression in the "context" of a ``DataFrame``.

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

   df = DataFrame(randn(5, 2), columns=['a', 'b'])
   df.eval('a + b')

Any expression that is a valid :func:`~pandas.eval` expression is also a valid
``DataFrame.eval`` expression, with the added benefit that *you don't have to
prefix the name of the* ``DataFrame`` *to the column(s) you're interested in
evaluating*.

In addition, you can perform assignment of columns within an expression.
This allows for *formulaic evaluation*. Only a single assignment is permitted.
The assignment target can be a new column name or an existing column name, and
it must be a valid Python identifier.

.. ipython:: python

   df = DataFrame(dict(a=range(5), b=range(5, 10)))
   df.eval('c = a + b')
   df.eval('d = a + b + c')
   df.eval('a = 1')
   df

Local Variables
~~~~~~~~~~~~~~~

You can refer to local variables the same way you would in vanilla Python

.. ipython:: python

   df = DataFrame(randn(5, 2), columns=['a', 'b'])
   newcol = randn(len(df))
   df.eval('b + newcol')

.. note::

   The one exception is when you have a local (or global) with the same name as
   a column in the ``DataFrame``

    .. code-block:: python

       df = DataFrame(randn(5, 2), columns=['a', 'b'])
       a = randn(len(df))
       df.eval('a + b')
       NameResolutionError: resolvers and locals overlap on names ['a']


   To deal with these conflicts, a special syntax exists for referring
   variables with the same name as a column

    .. ipython:: python
       :suppress:

       a = randn(len(df))

    .. ipython:: python

       df.eval('@a + b')

   The same is true for :meth:`~pandas.DataFrame.query`

    .. ipython:: python

       df.query('@a < b')

    .. ipython:: python
       :suppress:

       del a


:func:`~pandas.eval` Parsers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two different parsers and and two different engines you can use as
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


:func:`~pandas.eval` Backends
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There's also the option to make :func:`~pandas.eval` operate identical to plain
ol' Python.

.. note::

   Using the ``'python'`` engine is generally *not* useful, except for testing
   other :func:`~pandas.eval` engines against it. You will acheive **no**
   performance benefits using :func:`~pandas.eval` with ``engine='python'``.

You can see this by using :func:`~pandas.eval` with the ``'python'`` engine is
actually a bit slower (not by much) than evaluating the same expression in
Python:

.. ipython:: python

   %timeit df1 + df2 + df3 + df4

.. ipython:: python

   %timeit pd.eval('df1 + df2 + df3 + df4', engine='python')


:func:`~pandas.eval` Performance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`~pandas.eval` is intended to speed up certain kinds of operations. In
particular, those operations involving complex expressions with large
``DataFrame``/``Series`` objects should see a significant performance benefit.
Here is a plot showing the running time of :func:`~pandas.eval` as function of
the size of the frame involved in the computation. The two lines are two
different engines.


.. image:: _static/eval-perf.png


.. note::

   Operations with smallish objects (around 15k-20k rows) are faster using
   plain Python:

       .. image:: _static/eval-perf-small.png


This plot was created using a ``DataFrame`` with 3 columns each containing
floating point values generated using ``numpy.random.randn()``.

Technical Minutia
~~~~~~~~~~~~~~~~~
- Expressions that would result in an object dtype (including simple
  variable evaluation) have to be evaluated in Python space. The main reason
  for this behavior is to maintain backwards compatbility with versions of
  numpy < 1.7. In those versions of ``numpy`` a call to ``ndarray.astype(str)``
  will truncate any strings that are more than 60 characters in length. Second,
  we can't pass ``object`` arrays to ``numexpr`` thus string comparisons must
  be evaluated in Python space.
- The upshot is that this *only* applies to object-dtype'd expressions. So,
  if you have an expression--for example--that's a string comparison
  ``and``-ed together with another boolean expression that's from a numeric
  comparison, the numeric comparison will be evaluated by ``numexpr``. In fact,
  in general, :func:`~pandas.query`/:func:`~pandas.eval` will "pick out" the
  subexpressions that are ``eval``-able by ``numexpr`` and those that must be
  evaluated in Python space transparently to the user.
