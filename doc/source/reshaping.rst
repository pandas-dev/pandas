.. currentmodule:: pandas
.. _reshaping:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   import pandas.util.testing as tm
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

***************************
Pivoting and reshaping data
***************************

.. note::

   Since some of the functionality documented in this section is very new, the
   user should keep an eye on any changes to the API or behavior which may
   occur by the next release.

Reshaping by pivoting DataFrame objects
---------------------------------------

.. ipython::
   :suppress:

   In [1]: import pandas.util.testing as tm; tm.N = 3

   In [2]: def unpivot(frame):
      ...:         N, K = frame.shape
      ...:         data = {'value' : frame.values.ravel('F'),
      ...:                 'variable' : np.asarray(frame.columns).repeat(N),
      ...:                 'date' : np.tile(np.asarray(frame.index), K)}
      ...:         columns = ['date', 'variable', 'value']
      ...:         return DataFrame(data, columns=columns)
      ...:

   In [3]: df = unpivot(tm.makeTimeDataFrame())

Data is often stored in CSV files or databases in so-called "stacked" or
"record" format:

.. ipython:: python

   df


For the curious here is how the above DataFrame was created:

.. code-block:: python

   import pandas.util.testing as tm; tm.N = 3
   def unpivot(frame):
       N, K = frame.shape
       data = {'value' : frame.values.ravel('F'),
               'variable' : np.asarray(frame.columns).repeat(N),
               'date' : np.tile(np.asarray(frame.index), K)}
       return DataFrame(data, columns=['date', 'variable', 'value'])
   df = unpivot(tm.makeTimeDataFrame())

To select out everything for variable ``A`` we could do:

.. ipython:: python

   df[df['variable'] == 'A']

But if we wished to do time series operations between variables, this will
hardly do at all. This is really just a representation of a DataFrame whose
``columns`` are formed from the unique ``variable`` values and ``index`` from
the ``date`` values. To reshape the data into this form, use the ``pivot``
function:

.. ipython:: python

   df.pivot(index='date', columns='variable', values='value')

If the ``values`` argument is omitted, the resulting "pivoted" DataFrame will
have :ref:`hierarchical columns <indexing.hierarchical>` with the top level
being the set of value columns:

.. ipython:: python

   df['value2'] = df['value'] * 2
   pivoted = df.pivot('date', 'variable')
   pivoted

You of course can then select subsets from the pivoted DataFrame:

.. ipython:: python

   pivoted['value2']

Note that this returns a view on the underlying data in the case where the data
are homogeneously-typed.

Reshaping by stacking and unstacking
------------------------------------

Closely related to the ``pivot`` function are the related ``stack`` and
``unstack`` functions currently available on Series and DataFrame. These
functions are designed to tie together with ``MultiIndex`` objects (see the
section on :ref:`hierarchical indexing <indexing.hierarchical>`). Here are
essentially what these functions do:

  - ``stack``: collapse level in ``axis=1`` to produce new object whose index
    has the collapsed columns as its lowest level
  - ``unstack``: inverse operation from ``stack``; "pivot" index level to
    produce reshaped DataFrame

Actually very hard to explain in words; the clearest way is by example.
