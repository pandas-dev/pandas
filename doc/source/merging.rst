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

pandas has full-featured, **high performance** in-memory join operations
idiomatically very similar to relational databases like SQL. These methods
perform significantly better (in some cases well over an order of magnitude
better) than other open source implementations (like ``base::merge.data.frame``
in R). The reason for this is careful algorithmic design and internal layout of
the data in DataFrame.

.. _merging.join:

Database-style DataFrame joining/merging
----------------------------------------

pandas provides a single function, ``merge``, as the entry point for all
standard database join operations between DataFrame objects:

::

    merge(left, right, how='left', on=None, left_on=None, right_on=None,
          left_index=False, right_index=False, sort=True,
          suffixes=('.x', '.y'), copy=True)

Here's a description of what each argument is for:

  - ``left``: A DataFrame object
  - ``right``: Another DataFrame object
  - ``on``: Columns (names) to join on. Must be found in both the left and
    right DataFrame objects. If not passed and ``left_index`` and
    ``right_index`` are ``False``, the intersectino of the columns in the
    DataFrames will be inferred to be the join keys
  - ``left_on``: Columns from the left DataFrame to use as keys. Can either be
    column names or arrays with length equal to the length of the DataFrame
  - ``right_on``: Columns from the left DataFrame to use as keys. Can either be
    column names or arrays with length equal to the length of the DataFrame
  - ``left_index``: If ``True``, use the index (row labels) from the left
    DataFrame as its join key(s). In the case of a DataFrame with a MultiIndex
    (hierarchical), the number of levels must match the number of join keys
    from the right DataFrame
  - ``right_index``: Same usage as ``left_index`` for the right DataFrame
  - ``how``: One of ``'left'``, ``'right'``, ``'outer'``, ``'inner'``. Defaults
    to ``inner``. See below for more detailed description of each method
  - ``sort``: Sort the result DataFrame by the join keys in lexicographical
    order. Defaults to ``True``, setting to ``False`` will improve performance
    substantially in many cases
  - ``suffixes``: A tuple of string suffixes to apply to overlapping
    columns. Defaults to ``('.x', '.y')``.
  - ``copy``: Always copy data (default ``True``) from the passed DataFrame
    objects, even when reindexing is not necessary. Cannot be avoided in many
    cases but may improve performance / memory usage. The cases where copying
    can be avoided are somewhat pathological but this option is provided
    nonetheless.

``merge`` is a function in the pandas namespace, and it is also available as a
DataFrame instance method, with the calling DataFrame being implicitly
considered the left object in the join.

The related ``DataFrame.join`` method, uses ``merge`` internally for the
index-on-index and index-on-column(s) joins, but *joins on indexes* by default
rather than trying to join on common columns (the default behavior for
``merge``). If you are joining on index, you may wish to use ``DataFrame.join``
to save yourself some typing.

Brief primer on merge methods (relational algebra)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Experienced users of relational databases like SQL will be familiar with the
terminology used to describe join operations between two SQL-table like
structures (DataFrame objects). There are several cases to consider which are
very important to understand:

  - **one-to-one** joins: for example when joining two DataFrame objects on
    their indexes (which must contain unique values)
  - **many-to-one** joins: for example when joining an index (unique) to one or
    more columns in a DataFrame
  - **many-to-many** joins: joining columns on columns.

.. note::

   When joining columns on columns (potentially a many-to-many join), any
   indexes on the passed DataFrame objects **will be discarded**.


It is worth spending some time understanding the result of the **many-to-many**
join case. In SQL / standard relational algebra, if a key combination appears
more than once in both tables, the resulting table will have the **Cartesian
product** of the associated data. Here is a very basic example with one unique
key combination:

.. ipython:: python

   left = DataFrame({'key': ['foo', 'foo'], 'lval': [1, 2]})
   right = DataFrame({'key': ['foo', 'foo'], 'rval': [4, 5]})
   left
   right
   merge(left, right, on='key')

Here is a more complicated example with multiple join keys:

.. ipython:: python

   left = DataFrame({'key1': ['foo', 'foo', 'bar'],
                     'key2': ['one', 'two', 'one'],
                     'lval': [1, 2, 3]})
   right = DataFrame({'key1': ['foo', 'foo', 'bar', 'bar'],
                      'key2': ['one', 'one', 'one', 'two'],
                      'rval': [4, 5, 6, 7]})
   merge(left, right, how='outer')
   merge(left, right, how='inner')

The ``how`` argument to ``merge`` specifies how to determine which keys are to
be included in the resulting table. If a key combination **does not appear** in
either the left or right tables, the values in the joined table will be
``NA``. Here is a summary of the ``how`` options and their SQL equivalent names:

.. csv-table::
    :header: "Merge method", "SQL Join Name", "Description"
    :widths: 20, 20, 60

    ``left``, ``LEFT OUTER JOIN``, Use keys from left frame only
    ``right``, ``RIGHT OUTER JOIN``, Use keys from right frame only
    ``outer``, ``FULL OUTER JOIN``, Use union of keys from both frames
    ``inner``, ``INNER JOIN``, Use intersection of keys from both frames

Note that if using the index from either the left or right DataFrame (or both)
using the ``left_index`` / ``right_index`` options, the join operation is no
longer a many-to-many join by construction, as the index values are necessarily
unique. There will be some examples of this below.

.. _merging.join.index:

Joining on index
~~~~~~~~~~~~~~~~

``DataFrame.join`` is a convenient method for combining the columns of two
potentially differently-indexed DataFrames into a single result DataFrame. Here
is a very basic example:

.. ipython:: python

   df = DataFrame(np.random.randn(8, 4), columns=['A','B','C','D'])
   df1 = df.ix[1:, ['A', 'B']]
   df2 = df.ix[:5, ['C', 'D']]
   df1
   df2
   df1.join(df2)
   df1.join(df2, how='outer')
   df1.join(df2, how='inner')

The data alignment here is on the indexes (row labels). This same behavior can
be achieved using ``merge`` plus additional arguments instructing it to use the
indexes:

.. ipython:: python

   merge(df1, df2, left_index=True, right_index=True, how='outer')

Joining key columns on an index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``join`` takes an optional ``on`` argument which may be a column or multiple
column names, which specifies that the passed DataFrame is to be aligned on
that column in the DataFrame. These two function calls are completely
equivalent:

::

    left.join(right, on=key_or_keys)
    merge(left, right, left_on=key_or_keys, right_index=True,
          how='left', sort=False)

Obviously you can choose whichever form you find more convenient. For
many-to-one joins (where one of the DataFrame's is already indexed by the join
key), using ``join`` may be more convenient. Here is a simple example:

.. ipython:: python

   df['key'] = ['foo', 'bar'] * 4
   to_join = DataFrame(randn(2, 2), index=['bar', 'foo'],
                       columns=['j1', 'j2'])
   df
   to_join
   df.join(to_join, on='key')
   merge(df, to_join, left_on='key', right_index=True,
         how='left', sort=False)

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

Now this can be joined by passing the two key column names:

.. ipython:: python

   data.join(to_join, on=['key1', 'key2'])

.. _merging.df_inner_join:

The default for ``DataFrame.join`` is to perform a left join (essentially a
"VLOOKUP" operation, for Excel users), which uses only the keys found in the
calling DataFrame. Other join types, for example inner join, can be just as
easily performed:

.. ipython:: python

   data.join(to_join, on=['key1', 'key2'], how='inner')

As you can see, this drops any rows where there was no match.

Overlapping value columns
~~~~~~~~~~~~~~~~~~~~~~~~~

The merge ``suffixes`` argument takes a tuple of list of strings to append to
overlapping column names in the input DataFrames to disambiguate the result
columns:

.. ipython:: python

   left = DataFrame({'key': ['foo', 'foo'], 'value': [1, 2]})
   right = DataFrame({'key': ['foo', 'foo'], 'value': [4, 5]})
   merge(left, right, on='key', suffixes=['_left', '_right'])

``DataFrame.join`` has ``lsuffix`` and ``rsuffix`` arguments which behave
similarly.

Joining multiple DataFrame objects at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This has not been implemented yet, but is due to be implemented soon.

.. _merging.append:

Appending DataFrame objects (row-wise)
--------------------------------------

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

.. _merging.append.row:

Appending single rows to a DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While not especially efficient (since a new object must be created), you can
append a row to a DataFrame by passing a Series to ``append``, which returns a
new DataFrame as above:

.. ipython:: python

   df = DataFrame(np.random.randn(8, 4), columns=['A','B','C','D'])
   df
   s = df.xs(3)
   df.append(s, ignore_index=True)
