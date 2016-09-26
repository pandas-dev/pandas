.. currentmodule:: pandas
.. _faq:

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

Byte-Ordering Issues
--------------------
Occasionally you may have to deal with data that were created on a machine with
a different byte order than the one on which you are running Python. To deal
with this issue you should convert the underlying NumPy array to the native
system byte order *before* passing it to Series/DataFrame/Panel constructors
using something similar to the following:

.. ipython:: python

   x = np.array(list(range(10)), '>i4') # big endian
   newx = x.byteswap().newbyteorder() # force native byteorder
   s = pd.Series(newx)

See `the NumPy documentation on byte order
<http://docs.scipy.org/doc/numpy/user/basics.byteswapping.html>`__ for more
details.


Visualizing Data in Qt applications
-----------------------------------

There is no support for such visualization in pandas. However, the external
package `pandas-qt <https://github.com/datalyze-solutions/pandas-qt>`_ does
provide this functionality.
