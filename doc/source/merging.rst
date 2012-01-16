.. currentmodule:: pandas
.. _merging:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

****************************
Merge, join, and concatenate
****************************

pandas provides various facilities for easily combining together Series,
DataFrame, and Panel objects with various kinds of set logic for the indexes
and relational algebra functionality in the case of join / merge-type
operations.

.. _merging.concat:

Concatenating objects
---------------------

The ``concat`` function (in the main pandas namespace) does all of the heavy
lifting of performing concatenation operations along an axis while performing
optional set logic (union or intersection) of the indexes (if any) on the other
axes. Note that I say "if any" because there is only a single possible axis of
concatenation for Series.

Before diving into all of the details of ``concat`` and what it can do, here is
a simple example:

.. ipython:: python

   df = DataFrame(np.random.randn(10, 4))
   df

   # break it into pieces
   pieces = [df[:3], df[3:7], df[7:]]

   concatenated = concat(pieces)
   concatenated

Like its sibling function on ndarrays, ``numpy.concatenate``, ``pandas.concat``
takes a list or dict of homogeneously-typed objects and concatenates them with
some configurable handling of "what to do with the other axes":

::

    concat(objs, axis=0, join='outer', join_axes=None, ignore_index=False,
           keys=None, levels=None, names=None, verify_integrity=False)

- ``objs``: list or dict of Series, DataFrame, or Panel objects. If a dict is
  passed, the sorted keys will be used as the `keys` argument, unless it is
  passed, in which case the values will be selected (see below)
- ``axis``: {0, 1, ...}, default 0. The axis to concatenate along
- ``join``: {'inner', 'outer'}, default 'outer'. How to handle indexes on
  other axis(es). Outer for union and inner for intersection
- ``join_axes``: list of Index objects. Specific indexes to use for the other
  n - 1 axes instead of performing inner/outer set logic
- ``keys``: sequence, default None. Construct hierarchical index using the
  passed keys as the outermost level If multiple levels passed, should
  contain tuples.
- ``levels`` : list of sequences, default None. If keys passed, specific
  levels to use for the resulting MultiIndex. Otherwise they will be inferred
  from the keys
- ``names``: list, default None. Names for the levels in the resulting
  hierarchical index
- ``verify_integrity``: boolean, default False. Check whether the new
  concatenated axis contains duplicates. This can be very expensive relative
  to the actual data concatenation
- ``ignore_index`` : boolean, default False. If True, do not use the index
  values on the concatenation axis. The resulting axis will be labeled 0, ...,
  n - 1. This is useful if you are concatenating objects where the
  concatenation axis does not have meaningful indexing information.

Without a little bit of context and example many of these arguments don't make
much sense. Let's take the above example. Suppose we wanted to associate
specific keys with each of the pieces of the chopped up DataFrame. We can do
this using the ``keys`` argument:

.. ipython:: python

   concatenated = concat(pieces, keys=['first', 'second', 'third'])
   concatenated

As you can see (if you've read the rest of the documentation), the resulting
object's index has a :ref:`hierarchical index <indexing.hierarchical>`. This
means that we can now do stuff like select out each chunk by key:

.. ipython:: python

   concatenated.ix['second']

It's not a stretch to see how this can be very useful. More detail on this
functionality below.

Set logic on the other axes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When gluing together multiple DataFrames (or Panels or...), for example, you
have a choice of how to handle the other axes (other than the one being
concatenated). This can be done in three ways:

- Take the (sorted) union of them all, ``join='outer'``. This is the default
  option as it results in zero information loss.
- Take the intersection, ``join='inner'``.
- Use a specific index (in the case of DataFrame) or indexes (in the case of
  Panel or future higher dimensional objects), i.e. the ``join_axes`` argument

Here is a example of each of these methods. First, the default ``join='outer'``
behavior:

.. ipython:: python

   from pandas.util.testing import rands
   df = DataFrame(np.random.randn(10, 4), columns=['a', 'b', 'c', 'd'],
                  index=[rands(5) for _ in xrange(10)])
   df

   concat([df.ix[:7, ['a', 'b']], df.ix[2:-2, ['c']],
           df.ix[-7:, ['d']]], axis=1)

Note that the row indexes have been unioned and sorted. Here is the same thing
with ``join='inner'``:

.. ipython:: python

   concat([df.ix[:7, ['a', 'b']], df.ix[2:-2, ['c']],
           df.ix[-7:, ['d']]], axis=1, join='inner')

Lastly, suppose we just wanted to reuse the *exact index* from the original
DataFrame:

.. ipython:: python

   concat([df.ix[:7, ['a', 'b']], df.ix[2:-2, ['c']],
           df.ix[-7:, ['d']]], axis=1, join_axes=[df.index])

.. _merging.concatenation:

Concatenating using ``append``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A useful shortcut to ``concat`` are the ``append`` instance methods on Series
and DataFrame. These methods actually predated ``concat``. They concatenate
along ``axis=0``, namely the index:

.. ipython:: python

   s = Series(randn(10), index=np.arange(10))
   s1 = s[:5] # note we're slicing with labels here, so 5 is included
   s2 = s[6:]
   s1.append(s2)

In the case of DataFrame, the indexes must be disjoint but the columns do not
need to be:

.. ipython:: python

   df = DataFrame(randn(6, 4), index=DateRange('1/1/2000', periods=6),
                  columns=['A', 'B', 'C', 'D'])
   df1 = df.ix[:3]
   df2 = df.ix[3:, :3]
   df1
   df2
   df1.append(df2)

``append`` may take multiple objects to concatenate:

.. ipython:: python

   df1 = df.ix[:2]
   df2 = df.ix[2:4]
   df3 = df.ix[4:]
   df1.append([df2,df3])

.. note::

   Unlike `list.append` method, which appends to the original list and
   returns nothing, ``append`` here **does not** modify ``df1`` and
   returns its copy with ``df2`` appended.

.. _merging.ignore_index:

Ignoring indexes on the concatenation axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For DataFrames which don't have a meaningful index, you may wish to append them
and ignore the fact that they may have overlapping indexes:

.. ipython:: python

   df1 = DataFrame(randn(6, 4), columns=['A', 'B', 'C', 'D'])
   df2 = DataFrame(randn(3, 4), columns=['A', 'B', 'C', 'D'])

   df1
   df2

To do this, use the ``ignore_index`` argument:

.. ipython:: python

   concat([df1, df2], ignore_index=True)

This is also a valid argument to ``DataFrame.append``:

.. ipython:: python

   df1.append(df2, ignore_index=True)


More concatenating with group keys
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let's consider a variation on the first example presented:

.. ipython:: python

   df = DataFrame(np.random.randn(10, 4))
   df

   # break it into pieces
   pieces = [df.ix[:, [0, 1]], df.ix[:, [2]], df.ix[:, [3]]]

   result = concat(pieces, axis=1, keys=['one', 'two', 'three'])
   result

You can also pass a dict to ``concat`` in which case the dict keys will be used
for the ``keys`` argument (unless other keys are specified):

.. ipython:: python

   pieces = {'one': df.ix[:, [0, 1]],
             'two': df.ix[:, [2]],
             'three': df.ix[:, [3]]}
   concat(pieces, axis=1)
   concat(pieces, keys=['three', 'two'])

The MultiIndex created has levels that are constructed from the passed keys and
the columns of the DataFrame pieces:

.. ipython:: python

   result.columns.levels

If you wish to specify other levels (as will occasionally be the case), you can
do so using the ``levels`` argument:

.. ipython:: python

   result = concat(pieces, axis=1, keys=['one', 'two', 'three'],
                   levels=[['three', 'two', 'one', 'zero']],
                   names=['group_key'])
   result
   result.columns.levels

Yes, this is fairly esoteric, but is actually necessary for implementing things
like GroupBy where the order of a categorical variable is meaningful.

.. _merging.append.row:

Appending rows to a DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While not especially efficient (since a new object must be created), you can
append a single row to a DataFrame by passing a Series or dict to ``append``,
which returns a new DataFrame as above.

.. ipython:: python

   df = DataFrame(np.random.randn(8, 4), columns=['A','B','C','D'])
   df
   s = df.xs(3)
   df.append(s, ignore_index=True)

You should use ``ignore_index`` with this method to instruct DataFrame to
discard its index. If you wish to preserve the index, you should construct an
appropriately-indexed DataFrame and append or concatenate those objects.

You can also pass a list of dicts or Series:

.. ipython:: python

   df = DataFrame(np.random.randn(5, 4),
                  columns=['foo', 'bar', 'baz', 'qux'])
   dicts = [{'foo': 1, 'bar': 2, 'baz': 3, 'peekaboo': 4},
            {'foo': 5, 'bar': 6, 'baz': 7, 'peekaboo': 8}]
   result = df.append(dicts, ignore_index=True)
   result

.. _merging.join:

Database-style DataFrame joining/merging
----------------------------------------

pandas has full-featured, **high performance** in-memory join operations
idiomatically very similar to relational databases like SQL. These methods
perform significantly better (in some cases well over an order of magnitude
better) than other open source implementations (like ``base::merge.data.frame``
in R). The reason for this is careful algorithmic design and internal layout of
the data in DataFrame.

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

.. _merging.multiple_join:

Joining multiple DataFrame or Panel objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A list or tuple of DataFrames can also be passed to ``DataFrame.join`` to join
them together on their indexes. The same is true for ``Panel.join``.

.. ipython:: python

   df1 = df.ix[:, ['A', 'B']]
   df2 = df.ix[:, ['C', 'D']]
   df3 = df.ix[:, ['key']]
   df1
   df1.join([df2, df3])

