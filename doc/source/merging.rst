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

Appending DataFrame objects
---------------------------

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

Appending record-array like DataFrames
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For DataFrames which don't have a meaningful index, you may wish to append them
and ignore the fact that they may have overlapping indexes:

.. ipython:: python

   df1 = DataFrame(randn(6, 4), columns=['A', 'B', 'C', 'D'])
   df2 = DataFrame(randn(3, 4), columns=['A', 'B', 'C', 'D'])

   df1
   df2

.. _merging.ignore_index:

To do this, use the ``ignore_index`` argument:

.. ipython:: python

   df1.append(df2, ignore_index=True)


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
calling DataFrame which will be used to "align" the passed DataFrame. The
joining currently aligns the calling DataFrame's column (or columns) on the
passed DataFrame's index. This is best illustrated by example:

.. ipython:: python

   df['key'] = ['foo', 'bar'] * 3
   to_join = DataFrame(randn(2, 2), index=['bar', 'foo'],
                       columns=['j1', 'j2'])
   df
   to_join
   df.join(to_join, on='key')

.. _merging.multikey_join:

To join on multiple keys, the passed DataFrame must have a ``MultiIndex``:

.. ipython:: python

   index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                              ['one', 'two', 'three']],
                      labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                              [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                      names=['first', 'second'])
   to_join = DataFrame(np.random.randn(10, 3), index=index,
                       columns=['j_one', 'j_two', 'j_three'])

   # a little relevant example with NAs
   key1 = ['bar', 'bar', 'bar', 'foo', 'foo', 'baz', 'baz', 'qux',
           'qux', 'snap']
   key2 = ['two', 'one', 'three', 'one', 'two', 'one', 'two', 'two',
           'three', 'one']

   data = np.random.randn(len(key1))
   data = DataFrame({'key1' : key1, 'key2' : key2,
                     'data' : data})
   data
   to_join


.. ipython:: python

   data.join(to_join, on=['key1', 'key2'])

.. _merging.df_inner_join:

This is by default a "many-to-one" or "VLOOKUP"-style left join operation. An
inner join is also supported:

.. ipython:: python

   data.join(to_join, on=['key1', 'key2'], how='inner')

This drops any rows where there was no match.

Merging ordered records
~~~~~~~~~~~~~~~~~~~~~~~

This has not been implemented yet

Joining multiple DataFrame objects at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This has not been implemented yet, but is due to be implemented soon.
