.. _cython:

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


****************************************
Cython (Writing C extensions for pandas)
****************************************

For many use cases writing pandas in pure python and numpy is sufficient. In some computationally heavy applications however, it can be possible to achieve sizeable speed-ups by offloading work to `cython <http://cython.org/>`_.

- Say something about this being tutorial for "advanced" users?

.. note::

    The first thing to do here is to see if we can refactor in python, removing for loops (TODO add some waffle, and maybe trivial example, maybe even just using a for loop rather than apply in this example) a way which could make use of numpy...


This tutorial walksthrough a "typical" process of cythonizing a slow computation, we use an `example from the cython documentation <http://docs.cython.org/src/quickstart/cythonize.html>`_ in the context of pandas:

We have a function, ``integrate_f``, which we want to apply row-wise across a DataFrame, ``df``:

.. ipython:: python

   df = DataFrame({'x': 'x', 'a': randn(1000), 'b': randn(1000),'N': randint(100, 1000, (1000))})
   df

.. ipython:: python

   def f(x):
       return x * (x - 1)

   def integrate_f(a, b, N):
       s = 0
       dx = (b - a) / N
       for i in range(N):
           s += f(a + i * dx)
       return s * dx

In pure pandas we might achieve this using a row-wise ``apply``:

.. ipython:: python
   
   %timeit df.apply(lambda x: integrate_f(x['a'], x['b'], x['N']), axis=1)

Clearly this isn't fast enough for us, so let's take a look and see where the time is spent performing this operation (limited to the most time consuming four calls) using the `prun ipython magic function <http://ipython.org/ipython-doc/stable/api/generated/IPython.core.magics.execution.html#IPython.core.magics.execution.ExecutionMagics.prun>`_:

.. ipython:: python

   %prun -l 4 df.apply(lambda x: integrate_f(x['a'], x['b'], x['N']), axis=1)

By far the majority of time is spend inside either ``integrate_f`` or ``f``, hence we concentrate our efforts cythonizing these two functions.

.. note::
 
    In python 2 replacing the ``range`` with its generator counterpart (``xrange``) would mean the ``range`` line would vanish. In python 3 range is already a generator.

First, let's simply just copy our function over to cython as is (here the ``_plain`` suffix stands for "plain cython", allowing us to distinguish between our cython functions):

.. ipython:: python

   %load_ext cythonmagic

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

   %timeit df.apply(lambda x: integrate_f_plain(x['a'], x['b'], x['N']), axis=1)


We're already shaved a third off, not too bad for a simple copy and paste. We'll get another huge improvement simply by providing type information:

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

Now, we're talking! Already we're over ten times faster than the original python version, and we haven't *really* modified the code. Let's go back and have another look at what's eating up time now:

.. ipython:: python

   %prun -l 4 df.apply(lambda x: integrate_f_typed(x['a'], x['b'], x['N']), axis=1)

It's calling series and frames... a lot, in fact they're getting called for every row in the DataFrame. Function calls are expensive in python, so maybe we should cythonize the apply part and see if we can minimise these.

We are now passing ndarrays into the cython function, fortunately cython plays very nicely with numpy. TODO mention the ``Py_ssize_t``.

.. ipython:: 

   In [4]: %%cython
      ...: cimport numpy as np
      ...: import numpy as np
      ...: cdef double f_typed(double x) except? -2:
      ...:     return x**2-x
      ...: cpdef double integrate_f_typed(double a, double b, int N):
      ...:     cdef int i
      ...:     cdef double s, dx
      ...:     s = 0
      ...:     dx = (b-a)/N
      ...:     for i in range(N):
      ...:         s += f_typed(a+i*dx)
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


We create an array of zeros and loop over the rows, applying our ``integrate_f_typed`` function to fill it up. It's worth mentioning here that although a loop like this would be extremely slow in python (TODO: "as we saw" considerably slower than the apply?) while looping over a numpy array in cython is *fast*.

.. ipython:: python

   %timeit apply_integrate_f(df['a'], df['b'], df['N'])

We've gone another three times faster! Let's check again where the time is spent:

.. ipython:: python

   %prun -l 4 apply_integrate_f(df['a'], df['b'], df['N'])

As on might expect, the majority of the time is now spent in ``apply_integrate_f``, so if we wanted to make anymore efficiencies we must continue to concentrate our efforts here...

TODO explain decorators, and why they make it so fast!

.. ipython::

   In [5]: %%cython
      ...: cimport cython
      ...: cimport numpy as np
      ...: import numpy as np
      ...: cdef double f_typed(double x) except? -2:
      ...:     return x**2-x
      ...: cpdef double integrate_f_typed(double a, double b, int N):
      ...:     cdef int i
      ...:     cdef double s, dx
      ...:     s = 0
      ...:     dx = (b-a)/N
      ...:     for i in range(N):
      ...:         s += f_typed(a+i*dx)
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

Again we've shaved another third off, so let's have a look at where the time is spent:

.. ipython:: python

   %prun -l 4 apply_integrate_f_wrap(df['a'], df['b'], df['N'])

We can see that now all the time appears to be spent in ``apply_integrate_f_wrap`` and not much anywhere else. It would make sense to continue looking here for efficiencies...

TODO more? Have a 2D ndarray example?

Using cython has made our calculation around 100 times faster than the original python only version, and yet we're left with something which doesn't look too dissimilar.

TODO some warning that you don't need to cythonize every function (!)

Further topics:

- One can also load in functions from other C modules you've already written.
- More??

Read more in the `cython docs <http://docs.cython.org/>`_.