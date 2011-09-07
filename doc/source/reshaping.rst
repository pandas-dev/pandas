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

Actually very hard to explain in words; the clearest way is by example. Let's
take a prior example data set from the hierarchical indexing section:

.. ipython:: python

   tuples = zip(*[['bar', 'bar', 'baz', 'baz',
                   'foo', 'foo', 'qux', 'qux'],
                  ['one', 'two', 'one', 'two',
                   'one', 'two', 'one', 'two']])
   index = MultiIndex.from_tuples(tuples)
   df = DataFrame(randn(8, 2), index=index, columns=['A', 'B'])
   df2 = df[:4]
   df2

The ``stack`` function "compresses" a level in the DataFrame's columns to
produce either:

  - A Series, in the case of a simple column Index
  - A DataFrame, in the case of a ``MultiIndex`` in the columns

If the columns have a ``MultiIndex``, you can choose which level to stack. The
stacked level becomes the new lowest level in a ``MultiIndex`` on the columns:

.. ipython:: python

   stacked = df2.stack()
   stacked

With a "stacked" DataFrame or Series (having a ``MultiIndex`` as the
``index``), the inverse operation of ``stack`` is ``unstack``, which by default
unstacks the **last level**:

.. ipython:: python

   stacked.unstack()
   stacked.unstack(1)
   stacked.unstack(0)

These functions are very intelligent about handling missing data and do not
expect each subgroup within the hierarchical index to have the same set of
labels. They also can handle the index being unsorted (but you can make it
sorted by calling ``sortlevel``, of course). Here is a more complex example:

.. ipython:: python

   columns = MultiIndex.from_tuples([('A', 'cat'), ('B', 'dog'),
                                     ('B', 'cat'), ('A', 'dog')])
   df = DataFrame(randn(8, 4), index=index, columns=columns)
   df2 = df.ix[[0, 1, 2, 4, 5, 7]]
   df2

As mentioned above, ``stack`` can be called with a ``level`` argument to select
which level in the columns to stack:

.. ipython:: python

   df2.stack(1)
   df2.stack(0)

Unstacking when the columns are a ``MultiIndex`` is also careful about doing
the right thing:

.. ipython:: python

   df[:3].unstack(0)
   df2.unstack(1)

Combining with stats and GroupBy
--------------------------------

It should be no shock that combining ``pivot`` / ``stack`` / ``unstack`` with
GroupBy and the basic Series and DataFrame statistical functions can produce
some very expressive and fast data manipulations.

.. ipython:: python

   df
   df.stack().mean(1).unstack()

   # same result, another way
   df.groupby(level=1, axis=1).mean()

   df.stack().groupby(level=1).mean()

   df.mean().unstack(0)
