.. currentmodule:: pandas
.. _merging:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

***************************
Merging / Joining data sets
***************************

Appending disjoint objects
--------------------------

Series and DataFrame have an ``append`` method which will glue together objects
each of whose ``index`` (Series labels or DataFrame rows) is mutually
exclusive.

.. ipython:: python

   s = Series(randn(10), index=np.arange(10))
   s1 = s[:5]
   s2 = s[-5:]
   s1.append(s2)

In the case of DataFrame, the indexes must be disjoint but the columns do not need to be:

.. ipython:: python

   df = DataFrame(randn(6, 4), index=DateRange('1/1/2000', periods=6),
                  columns=['A', 'B', 'C', 'D'])
   df1 = df.ix[:3]
   df2 = df.ix[3:, :3]
   df1
   df2
   df1.append(df2)

Joining / merging DataFrames
----------------------------

The ``join`` method on DataFrame combines objects with **disjoint columns** by
lining up the rows based on some logic. The default alignment matches ``index``
on ``index``, with a ``how`` argument with the following options:

  - ``how='left'``: use the calling DataFrame's index
  - ``how='right'``: use the index of the DataFrame passed to ``join``
  - ``how='inner'``: intersect the indexes
  - ``how='outer'``: take the union of the indexes

Here are some examples:

.. ipython:: python

   df1 = df.ix[1:, ['A', 'B']]
   df2 = df.ix[:5, ['C', 'D']]
   df1
   df2
   df1.join(df2) # defaults to how='left'
   df1.join(df2, how='outer')
   df1.join(df2, how='inner')

Joining on a key
~~~~~~~~~~~~~~~~

``join`` takes an optional ``on`` argument which should be a column name in the
calling DataFrame which will be used to "align" the passed DataFrame. This is
best illustrated by example:

.. ipython:: python

   df['key'] = ['foo', 'bar'] * 3
   to_join = DataFrame(randn(2, 2), index=['bar', 'foo'],
                       columns=['j1', 'j2'])
   df
   to_join
   df.join(to_join, on='key')

Merging ordered records
~~~~~~~~~~~~~~~~~~~~~~~

This has not been implemented yet

Joining multiple DataFrame objects at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This has not been implemented yet, but is due to be implemented soon.
