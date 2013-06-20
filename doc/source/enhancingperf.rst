.. _enhancingperf:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import os
   import csv
   from pandas import DataFrame
   import pandas as pd

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
speed-ups by offloading work to `cython <http://cython.org/>`_.

This tutorial assumes you have refactored as much as possible in python, for example
trying to remove for loops and making use of numpy vectorization, it's always worth
optimising in python first.

This tutorial walks through a "typical" process of cythonizing a slow computation.
We use an `example from the cython documentation <http://docs.cython.org/src/quickstart/cythonize.html>`_ 
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
four calls) using the `prun ipython magic function <http://ipython.org/ipython-doc/stable/api/generated/IPython.core.magics.execution.html#IPython.core.magics.execution.ExecutionMagics.prun>`_:

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


.. note::

    Loop like this would be *extremely* slow in python, but in cython looping over
    numpy arrays is *fast*.

.. ipython:: python

   %timeit apply_integrate_f(df['a'], df['b'], df['N'])

We've gone another three times faster! Let's check again where the time is spent:

.. ipython:: python

   %prun -l 4 apply_integrate_f(df['a'], df['b'], df['N'])

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

   %timeit apply_integrate_f_wrap(df['a'], df['b'], df['N'])

This shaves another third off!

Further topics
~~~~~~~~~~~~~~

- Loading C modules into cython.

Read more in the `cython docs <http://docs.cython.org/>`_.