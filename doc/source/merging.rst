.. currentmodule:: pandas
.. _merging:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   import pandas as pd
   pd.options.display.max_rows=15
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

   import matplotlib.pyplot as plt
   plt.close('all')
   import pandas.util._doctools as doctools
   p = doctools.TablePlotter()


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

   df1 = pd.DataFrame({'A': ['A0', 'A1', 'A2', 'A3'],
                       'B': ['B0', 'B1', 'B2', 'B3'],
                       'C': ['C0', 'C1', 'C2', 'C3'],
                       'D': ['D0', 'D1', 'D2', 'D3']},
                       index=[0, 1, 2, 3])

   df2 = pd.DataFrame({'A': ['A4', 'A5', 'A6', 'A7'],
                       'B': ['B4', 'B5', 'B6', 'B7'],
                       'C': ['C4', 'C5', 'C6', 'C7'],
                       'D': ['D4', 'D5', 'D6', 'D7']},
                        index=[4, 5, 6, 7])

   df3 = pd.DataFrame({'A': ['A8', 'A9', 'A10', 'A11'],
                       'B': ['B8', 'B9', 'B10', 'B11'],
                       'C': ['C8', 'C9', 'C10', 'C11'],
                       'D': ['D8', 'D9', 'D10', 'D11']},
                       index=[8, 9, 10, 11])

   frames = [df1, df2, df3]
   result = pd.concat(frames)

.. ipython:: python
   :suppress:

   @savefig merging_concat_basic.png
   p.plot(frames, result,
          labels=['df1', 'df2', 'df3'], vertical=True);
   plt.close('all');

Like its sibling function on ndarrays, ``numpy.concatenate``, ``pandas.concat``
takes a list or dict of homogeneously-typed objects and concatenates them with
some configurable handling of "what to do with the other axes":

::

    pd.concat(objs, axis=0, join='outer', join_axes=None, ignore_index=False,
              keys=None, levels=None, names=None, verify_integrity=False,
              copy=True)

- ``objs`` : a sequence or mapping of Series, DataFrame, or Panel objects. If a
  dict is passed, the sorted keys will be used as the `keys` argument, unless
  it is passed, in which case the values will be selected (see below). Any None
  objects will be dropped silently unless they are all None in which case a
  ValueError will be raised.
- ``axis`` : {0, 1, ...}, default 0. The axis to concatenate along.
- ``join`` : {'inner', 'outer'}, default 'outer'. How to handle indexes on
  other axis(es). Outer for union and inner for intersection.
- ``ignore_index`` : boolean, default False. If True, do not use the index
  values on the concatenation axis. The resulting axis will be labeled 0, ...,
  n - 1. This is useful if you are concatenating objects where the
  concatenation axis does not have meaningful indexing information. Note
  the index values on the other axes are still respected in the join.
- ``join_axes`` : list of Index objects. Specific indexes to use for the other
  n - 1 axes instead of performing inner/outer set logic.
- ``keys`` : sequence, default None. Construct hierarchical index using the
  passed keys as the outermost level. If multiple levels passed, should
  contain tuples.
- ``levels`` : list of sequences, default None. Specific levels (unique values)
  to use for constructing a MultiIndex. Otherwise they will be inferred from the
  keys.
- ``names`` : list, default None. Names for the levels in the resulting
  hierarchical index.
- ``verify_integrity`` : boolean, default False. Check whether the new
  concatenated axis contains duplicates. This can be very expensive relative
  to the actual data concatenation.
- ``copy`` : boolean, default True. If False, do not copy data unnecessarily.

Without a little bit of context and example many of these arguments don't make
much sense. Let's take the above example. Suppose we wanted to associate
specific keys with each of the pieces of the chopped up DataFrame. We can do
this using the ``keys`` argument:

.. ipython:: python

   result = pd.concat(frames, keys=['x', 'y', 'z'])

.. ipython:: python
   :suppress:

   @savefig merging_concat_keys.png
   p.plot(frames, result,
          labels=['df1', 'df2', 'df3'], vertical=True)
   plt.close('all');

As you can see (if you've read the rest of the documentation), the resulting
object's index has a :ref:`hierarchical index <advanced.hierarchical>`. This
means that we can now do stuff like select out each chunk by key:

.. ipython:: python

   result.loc['y']

It's not a stretch to see how this can be very useful. More detail on this
functionality below.

.. note::
   It is worth noting however, that ``concat`` (and therefore ``append``) makes
   a full copy of the data, and that constantly reusing this function can
   create a significant performance hit. If you need to use the operation over
   several datasets, use a list comprehension.

::

   frames = [ process_your_file(f) for f in files ]
   result = pd.concat(frames)


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

   df4 = pd.DataFrame({'B': ['B2', 'B3', 'B6', 'B7'],
                    'D': ['D2', 'D3', 'D6', 'D7'],
                    'F': ['F2', 'F3', 'F6', 'F7']},
                   index=[2, 3, 6, 7])
   result = pd.concat([df1, df4], axis=1)


.. ipython:: python
   :suppress:

   @savefig merging_concat_axis1.png
   p.plot([df1, df4], result,
          labels=['df1', 'df4'], vertical=False);
   plt.close('all');

Note that the row indexes have been unioned and sorted. Here is the same thing
with ``join='inner'``:

.. ipython:: python

   result = pd.concat([df1, df4], axis=1, join='inner')

.. ipython:: python
   :suppress:

   @savefig merging_concat_axis1_inner.png
   p.plot([df1, df4], result,
          labels=['df1', 'df4'], vertical=False);
   plt.close('all');

Lastly, suppose we just wanted to reuse the *exact index* from the original
DataFrame:

.. ipython:: python

   result = pd.concat([df1, df4], axis=1, join_axes=[df1.index])

.. ipython:: python
   :suppress:

   @savefig merging_concat_axis1_join_axes.png
   p.plot([df1, df4], result,
          labels=['df1', 'df4'], vertical=False);
   plt.close('all');

.. _merging.concatenation:

Concatenating using ``append``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A useful shortcut to ``concat`` are the ``append`` instance methods on Series
and DataFrame. These methods actually predated ``concat``. They concatenate
along ``axis=0``, namely the index:

.. ipython:: python

   result = df1.append(df2)

.. ipython:: python
   :suppress:

   @savefig merging_append1.png
   p.plot([df1, df2], result,
          labels=['df1', 'df2'], vertical=True);
   plt.close('all');

In the case of DataFrame, the indexes must be disjoint but the columns do not
need to be:

.. ipython:: python

   result = df1.append(df4)

.. ipython:: python
   :suppress:

   @savefig merging_append2.png
   p.plot([df1, df4], result,
          labels=['df1', 'df4'], vertical=True);
   plt.close('all');

``append`` may take multiple objects to concatenate:

.. ipython:: python

   result = df1.append([df2, df3])

.. ipython:: python
   :suppress:

   @savefig merging_append3.png
   p.plot([df1, df2, df3], result,
          labels=['df1', 'df2', 'df3'], vertical=True);
   plt.close('all');

.. note::

   Unlike `list.append` method, which appends to the original list and
   returns nothing, ``append`` here **does not** modify ``df1`` and
   returns its copy with ``df2`` appended.

.. _merging.ignore_index:

Ignoring indexes on the concatenation axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
For DataFrames which don't have a meaningful index, you may wish to append them
and ignore the fact that they may have overlapping indexes:

To do this, use the ``ignore_index`` argument:

.. ipython:: python

   result = pd.concat([df1, df4], ignore_index=True)

.. ipython:: python
   :suppress:

   @savefig merging_concat_ignore_index.png
   p.plot([df1, df4], result,
          labels=['df1', 'df4'], vertical=True);
   plt.close('all');

This is also a valid argument to ``DataFrame.append``:

.. ipython:: python

   result = df1.append(df4, ignore_index=True)

.. ipython:: python
   :suppress:

   @savefig merging_append_ignore_index.png
   p.plot([df1, df4], result,
          labels=['df1', 'df4'], vertical=True);
   plt.close('all');

.. _merging.mixed_ndims:

Concatenating with mixed ndims
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can concatenate a mix of Series and DataFrames. The
Series will be transformed to DataFrames with the column name as
the name of the Series.

.. ipython:: python

   s1 = pd.Series(['X0', 'X1', 'X2', 'X3'], name='X')
   result = pd.concat([df1, s1], axis=1)

.. ipython:: python
   :suppress:

   @savefig merging_concat_mixed_ndim.png
   p.plot([df1, s1], result,
          labels=['df1', 's1'], vertical=False);
   plt.close('all');

If unnamed Series are passed they will be numbered consecutively.

.. ipython:: python

   s2 = pd.Series(['_0', '_1', '_2', '_3'])
   result = pd.concat([df1, s2, s2, s2], axis=1)

.. ipython:: python
   :suppress:

   @savefig merging_concat_unnamed_series.png
   p.plot([df1, s2], result,
          labels=['df1', 's2'], vertical=False);
   plt.close('all');

Passing ``ignore_index=True`` will drop all name references.

.. ipython:: python

   result = pd.concat([df1, s1], axis=1, ignore_index=True)

.. ipython:: python
   :suppress:

   @savefig merging_concat_series_ignore_index.png
   p.plot([df1, s1], result,
          labels=['df1', 's1'], vertical=False);
   plt.close('all');

More concatenating with group keys
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A fairly common use of the ``keys`` argument is to override the column names when creating a new DataFrame based on existing Series.
Notice how the default behaviour consists on letting the resulting DataFrame inherits the parent Series' name, when these existed.

.. ipython:: python

   s3 = pd.Series([0, 1, 2, 3], name='foo')
   s4 = pd.Series([0, 1, 2, 3])
   s5 = pd.Series([0, 1, 4, 5])

   pd.concat([s3, s4, s5], axis=1)

Through the ``keys`` argument we can override the existing column names.

.. ipython:: python

   pd.concat([s3, s4, s5], axis=1, keys=['red','blue','yellow'])

Let's consider now a variation on the very first example presented:

.. ipython:: python

   result = pd.concat(frames, keys=['x', 'y', 'z'])

.. ipython:: python
   :suppress:

   @savefig merging_concat_group_keys2.png
   p.plot(frames, result,
          labels=['df1', 'df2', 'df3'], vertical=True);
   plt.close('all');

You can also pass a dict to ``concat`` in which case the dict keys will be used
for the ``keys`` argument (unless other keys are specified):

.. ipython:: python

   pieces = {'x': df1, 'y': df2, 'z': df3}
   result = pd.concat(pieces)

.. ipython:: python
   :suppress:

   @savefig merging_concat_dict.png
   p.plot([df1, df2, df3], result,
          labels=['df1', 'df2', 'df3'], vertical=True);
   plt.close('all');

.. ipython:: python

   result = pd.concat(pieces, keys=['z', 'y'])

.. ipython:: python
   :suppress:

   @savefig merging_concat_dict_keys.png
   p.plot([df1, df2, df3], result,
          labels=['df1', 'df2', 'df3'], vertical=True);
   plt.close('all');

The MultiIndex created has levels that are constructed from the passed keys and
the index of the DataFrame pieces:

.. ipython:: python

   result.index.levels

If you wish to specify other levels (as will occasionally be the case), you can
do so using the ``levels`` argument:

.. ipython:: python

   result = pd.concat(pieces, keys=['x', 'y', 'z'],
                   levels=[['z', 'y', 'x', 'w']],
                   names=['group_key'])

.. ipython:: python
   :suppress:

   @savefig merging_concat_dict_keys_names.png
   p.plot([df1, df2, df3], result,
          labels=['df1', 'df2', 'df3'], vertical=True);
   plt.close('all');

.. ipython:: python

   result.index.levels

Yes, this is fairly esoteric, but is actually necessary for implementing things
like GroupBy where the order of a categorical variable is meaningful.

.. _merging.append.row:

Appending rows to a DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While not especially efficient (since a new object must be created), you can
append a single row to a DataFrame by passing a Series or dict to ``append``,
which returns a new DataFrame as above.

.. ipython:: python

   s2 = pd.Series(['X0', 'X1', 'X2', 'X3'], index=['A', 'B', 'C', 'D'])
   result = df1.append(s2, ignore_index=True)

.. ipython:: python
   :suppress:

   @savefig merging_append_series_as_row.png
   p.plot([df1, s2], result,
          labels=['df1', 's2'], vertical=True);
   plt.close('all');

You should use ``ignore_index`` with this method to instruct DataFrame to
discard its index. If you wish to preserve the index, you should construct an
appropriately-indexed DataFrame and append or concatenate those objects.

You can also pass a list of dicts or Series:

.. ipython:: python

   dicts = [{'A': 1, 'B': 2, 'C': 3, 'X': 4},
            {'A': 5, 'B': 6, 'C': 7, 'Y': 8}]
   result = df1.append(dicts, ignore_index=True)

.. ipython:: python
   :suppress:

   @savefig merging_append_dits.png
   p.plot([df1, pd.DataFrame(dicts)], result,
          labels=['df1', 'dicts'], vertical=True);
   plt.close('all');

.. _merging.join:

Database-style DataFrame joining/merging
----------------------------------------

pandas has full-featured, **high performance** in-memory join operations
idiomatically very similar to relational databases like SQL. These methods
perform significantly better (in some cases well over an order of magnitude
better) than other open source implementations (like ``base::merge.data.frame``
in R). The reason for this is careful algorithmic design and internal layout of
the data in DataFrame.

See the :ref:`cookbook<cookbook.merge>` for some advanced strategies.

Users who are familiar with SQL but new to pandas might be interested in a
:ref:`comparison with SQL<compare_with_sql.join>`.

pandas provides a single function, ``merge``, as the entry point for all
standard database join operations between DataFrame objects:

::

    pd.merge(left, right, how='inner', on=None, left_on=None, right_on=None,
             left_index=False, right_index=False, sort=True,
             suffixes=('_x', '_y'), copy=True, indicator=False,
             validate=None)

- ``left``: A DataFrame object
- ``right``: Another DataFrame object
- ``on``: Column or index level names to join on. Must be found in both the left
  and right DataFrame objects. If not passed and ``left_index`` and
  ``right_index`` are ``False``, the intersection of the columns in the
  DataFrames will be inferred to be the join keys
- ``left_on``: Columns or index levels from the left DataFrame to use as
  keys. Can either be column names, index level names, or arrays with length
  equal to the length of the DataFrame
- ``right_on``: Columns or index levels from the right DataFrame to use as
  keys. Can either be column names, index level names, or arrays with length
  equal to the length of the DataFrame
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
  columns. Defaults to ``('_x', '_y')``.
- ``copy``: Always copy data (default ``True``) from the passed DataFrame
  objects, even when reindexing is not necessary. Cannot be avoided in many
  cases but may improve performance / memory usage. The cases where copying
  can be avoided are somewhat pathological but this option is provided
  nonetheless.
- ``indicator``: Add a column to the output DataFrame called ``_merge``
  with information on the source of each row. ``_merge`` is Categorical-type
  and takes on a value of ``left_only`` for observations whose merge key
  only appears in ``'left'`` DataFrame, ``right_only`` for observations whose
  merge key only appears in ``'right'`` DataFrame, and ``both`` if the
  observation's merge key is found in both.

- ``validate`` : string, default None.
  If specified, checks if merge is of specified type.

  * "one_to_one" or "1:1": checks if merge keys are unique in both
    left and right datasets.
  * "one_to_many" or "1:m": checks if merge keys are unique in left
    dataset.
  * "many_to_one" or "m:1": checks if merge keys are unique in right
    dataset.
  * "many_to_many" or "m:m": allowed, but does not result in checks.

  .. versionadded:: 0.21.0

.. note::

   Support for specifying index levels as the ``on``, ``left_on``, and
   ``right_on`` parameters was added in version 0.22.0.

The return type will be the same as ``left``. If ``left`` is a ``DataFrame``
and ``right`` is a subclass of DataFrame, the return type will still be
``DataFrame``.

``merge`` is a function in the pandas namespace, and it is also available as a
DataFrame instance method, with the calling DataFrame being implicitly
considered the left object in the join.

The related ``DataFrame.join`` method, uses ``merge`` internally for the
index-on-index (by default) and column(s)-on-index join. If you are joining on
index only, you may wish to use ``DataFrame.join`` to save yourself some typing.

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

   left = pd.DataFrame({'key': ['K0', 'K1', 'K2', 'K3'],
                        'A': ['A0', 'A1', 'A2', 'A3'],
                        'B': ['B0', 'B1', 'B2', 'B3']})

   right = pd.DataFrame({'key': ['K0', 'K1', 'K2', 'K3'],
                         'C': ['C0', 'C1', 'C2', 'C3'],
                         'D': ['D0', 'D1', 'D2', 'D3']})
   result = pd.merge(left, right, on='key')

.. ipython:: python
   :suppress:

   @savefig merging_merge_on_key.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

Here is a more complicated example with multiple join keys:

.. ipython:: python

   left = pd.DataFrame({'key1': ['K0', 'K0', 'K1', 'K2'],
                        'key2': ['K0', 'K1', 'K0', 'K1'],
                        'A': ['A0', 'A1', 'A2', 'A3'],
                        'B': ['B0', 'B1', 'B2', 'B3']})

   right = pd.DataFrame({'key1': ['K0', 'K1', 'K1', 'K2'],
                         'key2': ['K0', 'K0', 'K0', 'K0'],
                         'C': ['C0', 'C1', 'C2', 'C3'],
                         'D': ['D0', 'D1', 'D2', 'D3']})

   result = pd.merge(left, right, on=['key1', 'key2'])

.. ipython:: python
   :suppress:

   @savefig merging_merge_on_key_multiple.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

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

.. ipython:: python

   result = pd.merge(left, right, how='left', on=['key1', 'key2'])

.. ipython:: python
   :suppress:

   @savefig merging_merge_on_key_left.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. ipython:: python

   result = pd.merge(left, right, how='right', on=['key1', 'key2'])

.. ipython:: python
   :suppress:

   @savefig merging_merge_on_key_right.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);

.. ipython:: python

   result = pd.merge(left, right, how='outer', on=['key1', 'key2'])

.. ipython:: python
   :suppress:

   @savefig merging_merge_on_key_outer.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. ipython:: python

   result = pd.merge(left, right, how='inner', on=['key1', 'key2'])

.. ipython:: python
   :suppress:

   @savefig merging_merge_on_key_inner.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

Here is another example with duplicate join keys in DataFrames:

.. ipython:: python

   left = pd.DataFrame({'A' : [1,2], 'B' : [2, 2]})

   right = pd.DataFrame({'A' : [4,5,6], 'B': [2,2,2]})

   result = pd.merge(left, right, on='B', how='outer')

.. ipython:: python
   :suppress:

   @savefig merging_merge_on_key_dup.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');


.. warning::

  Joining / merging on duplicate keys can cause a returned frame that is the multiplication of the row dimensions, which may result in memory overflow. It is the user' s responsibility to manage duplicate values in keys before joining large DataFrames.

.. _merging.validation:

Checking for duplicate keys
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.21.0

Users can use the ``validate`` argument to automatically check whether there are unexpected duplicates in their merge keys. Key uniqueness is checked before merge operations and so should protect against memory overflows. Checking key uniqueness is also a good way to ensure user data structures are as expected. 

In the following example, there are duplicate values of ``B`` in the right DataFrame. As this is not a one-to-one merge -- as specified in the ``validate`` argument -- an exception will be raised.


.. ipython:: python

  left = pd.DataFrame({'A' : [1,2], 'B' : [1, 2]})
  right = pd.DataFrame({'A' : [4,5,6], 'B': [2, 2, 2]})

.. code-block:: ipython

  In [53]: result = pd.merge(left, right, on='B', how='outer', validate="one_to_one")
  ...
  MergeError: Merge keys are not unique in right dataset; not a one-to-one merge    

If the user is aware of the duplicates in the right `DataFrame` but wants to ensure there are no duplicates in the left DataFrame, one can use the `validate='one_to_many'` argument instead, which will not raise an exception. 

.. ipython:: python

   pd.merge(left, right, on='B', how='outer', validate="one_to_many")


.. _merging.indicator:

The merge indicator
~~~~~~~~~~~~~~~~~~~

``merge`` accepts the argument ``indicator``. If ``True``, a Categorical-type column called ``_merge`` will be added to the output object that takes on values:

  ===================================   ================
  Observation Origin                    ``_merge`` value
  ===================================   ================
  Merge key only in ``'left'`` frame    ``left_only``
  Merge key only in ``'right'`` frame   ``right_only``
  Merge key in both frames              ``both``
  ===================================   ================

.. ipython:: python

   df1 = pd.DataFrame({'col1': [0, 1], 'col_left':['a', 'b']})
   df2 = pd.DataFrame({'col1': [1, 2, 2],'col_right':[2, 2, 2]})
   pd.merge(df1, df2, on='col1', how='outer', indicator=True)

The ``indicator`` argument will also accept string arguments, in which case the indicator function will use the value of the passed string as the name for the indicator column.

.. ipython:: python

   pd.merge(df1, df2, on='col1', how='outer', indicator='indicator_column')


.. _merging.dtypes:

Merge Dtypes
~~~~~~~~~~~~

.. versionadded:: 0.19.0

Merging will preserve the dtype of the join keys.

.. ipython:: python

   left = pd.DataFrame({'key': [1], 'v1': [10]})
   left
   right = pd.DataFrame({'key': [1, 2], 'v1': [20, 30]})
   right

We are able to preserve the join keys

.. ipython:: python

   pd.merge(left, right, how='outer')
   pd.merge(left, right, how='outer').dtypes

Of course if you have missing values that are introduced, then the
resulting dtype will be upcast.

.. ipython:: python

   pd.merge(left, right, how='outer', on='key')
   pd.merge(left, right, how='outer', on='key').dtypes

.. versionadded:: 0.20.0

Merging will preserve ``category`` dtypes of the mergands. See also the section on :ref:`categoricals <categorical.merge>`

The left frame.

.. ipython:: python

   from pandas.api.types import CategoricalDtype

   X = pd.Series(np.random.choice(['foo', 'bar'], size=(10,)))
   X = X.astype(CategoricalDtype(categories=['foo', 'bar']))

   left = pd.DataFrame({'X': X,
                        'Y': np.random.choice(['one', 'two', 'three'], size=(10,))})
   left
   left.dtypes

The right frame.

.. ipython:: python

   right = pd.DataFrame({
        'X': pd.Series(['foo', 'bar'],
                       dtype=CategoricalDtype(['foo', 'bar'])),
        'Z': [1, 2]
   })
   right
   right.dtypes

The merged result

.. ipython:: python

   result = pd.merge(left, right, how='outer')
   result
   result.dtypes

.. note::

   The category dtypes must be *exactly* the same, meaning the same categories and the ordered attribute.
   Otherwise the result will coerce to ``object`` dtype.

.. note::

   Merging on ``category`` dtypes that are the same can be quite performant compared to ``object`` dtype merging.

.. _merging.join.index:

Joining on index
~~~~~~~~~~~~~~~~

``DataFrame.join`` is a convenient method for combining the columns of two
potentially differently-indexed DataFrames into a single result DataFrame. Here
is a very basic example:

.. ipython:: python

   left = pd.DataFrame({'A': ['A0', 'A1', 'A2'],
                        'B': ['B0', 'B1', 'B2']},
                        index=['K0', 'K1', 'K2'])

   right = pd.DataFrame({'C': ['C0', 'C2', 'C3'],
                         'D': ['D0', 'D2', 'D3']},
                         index=['K0', 'K2', 'K3'])

   result = left.join(right)

.. ipython:: python
   :suppress:

   @savefig merging_join.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. ipython:: python

   result = left.join(right, how='outer')

.. ipython:: python
   :suppress:

   @savefig merging_join_outer.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. ipython:: python

   result = left.join(right, how='inner')

.. ipython:: python
   :suppress:

   @savefig merging_join_inner.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

The data alignment here is on the indexes (row labels). This same behavior can
be achieved using ``merge`` plus additional arguments instructing it to use the
indexes:

.. ipython:: python

   result = pd.merge(left, right, left_index=True, right_index=True, how='outer')

.. ipython:: python
   :suppress:

   @savefig merging_merge_index_outer.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. ipython:: python

   result = pd.merge(left, right, left_index=True, right_index=True, how='inner');

.. ipython:: python
   :suppress:

   @savefig merging_merge_index_inner.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

Joining key columns on an index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``join`` takes an optional ``on`` argument which may be a column or multiple
column names, which specifies that the passed DataFrame is to be aligned on
that column in the DataFrame. These two function calls are completely
equivalent:

::

    left.join(right, on=key_or_keys)
    pd.merge(left, right, left_on=key_or_keys, right_index=True,
          how='left', sort=False)

Obviously you can choose whichever form you find more convenient. For
many-to-one joins (where one of the DataFrame's is already indexed by the join
key), using ``join`` may be more convenient. Here is a simple example:

.. ipython:: python

   left = pd.DataFrame({'A': ['A0', 'A1', 'A2', 'A3'],
                        'B': ['B0', 'B1', 'B2', 'B3'],
                        'key': ['K0', 'K1', 'K0', 'K1']})

   right = pd.DataFrame({'C': ['C0', 'C1'],
                         'D': ['D0', 'D1']},
                         index=['K0', 'K1'])

   result = left.join(right, on='key')

.. ipython:: python
   :suppress:

   @savefig merging_join_key_columns.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. ipython:: python

   result = pd.merge(left, right, left_on='key', right_index=True,
                     how='left', sort=False);

.. ipython:: python
   :suppress:

   @savefig merging_merge_key_columns.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. _merging.multikey_join:

To join on multiple keys, the passed DataFrame must have a ``MultiIndex``:

.. ipython:: python

   left = pd.DataFrame({'A': ['A0', 'A1', 'A2', 'A3'],
                        'B': ['B0', 'B1', 'B2', 'B3'],
                        'key1': ['K0', 'K0', 'K1', 'K2'],
                        'key2': ['K0', 'K1', 'K0', 'K1']})

   index = pd.MultiIndex.from_tuples([('K0', 'K0'), ('K1', 'K0'),
                                     ('K2', 'K0'), ('K2', 'K1')])
   right = pd.DataFrame({'C': ['C0', 'C1', 'C2', 'C3'],
                      'D': ['D0', 'D1', 'D2', 'D3']},
                     index=index)

Now this can be joined by passing the two key column names:

.. ipython:: python

   result = left.join(right, on=['key1', 'key2'])

.. ipython:: python
   :suppress:

   @savefig merging_join_multikeys.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. _merging.df_inner_join:

The default for ``DataFrame.join`` is to perform a left join (essentially a
"VLOOKUP" operation, for Excel users), which uses only the keys found in the
calling DataFrame. Other join types, for example inner join, can be just as
easily performed:

.. ipython:: python

   result = left.join(right, on=['key1', 'key2'], how='inner')

.. ipython:: python
   :suppress:

   @savefig merging_join_multikeys_inner.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

As you can see, this drops any rows where there was no match.

.. _merging.join_on_mi:

Joining a single Index to a Multi-index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can join a singly-indexed ``DataFrame`` with a level of a multi-indexed ``DataFrame``.
The level will match on the name of the index of the singly-indexed frame against
a level name of the multi-indexed frame.

..  ipython:: python

   left = pd.DataFrame({'A': ['A0', 'A1', 'A2'],
                        'B': ['B0', 'B1', 'B2']},
                        index=pd.Index(['K0', 'K1', 'K2'], name='key'))

   index = pd.MultiIndex.from_tuples([('K0', 'Y0'), ('K1', 'Y1'),
                                     ('K2', 'Y2'), ('K2', 'Y3')],
                                      names=['key', 'Y'])
   right = pd.DataFrame({'C': ['C0', 'C1', 'C2', 'C3'],
                         'D': ['D0', 'D1', 'D2', 'D3']},
                         index=index)

   result = left.join(right, how='inner')

.. ipython:: python
   :suppress:

   @savefig merging_join_multiindex_inner.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

This is equivalent but less verbose and more memory efficient / faster than this.

..  ipython:: python

    result = pd.merge(left.reset_index(), right.reset_index(),
          on=['key'], how='inner').set_index(['key','Y'])

.. ipython:: python
   :suppress:

   @savefig merging_merge_multiindex_alternative.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

Joining with two multi-indexes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is not Implemented via ``join`` at-the-moment, however it can be done using the following.

.. ipython:: python

   index = pd.MultiIndex.from_tuples([('K0', 'X0'), ('K0', 'X1'),
                                      ('K1', 'X2')],
                                       names=['key', 'X'])
   left = pd.DataFrame({'A': ['A0', 'A1', 'A2'],
                        'B': ['B0', 'B1', 'B2']},
                         index=index)

   result = pd.merge(left.reset_index(), right.reset_index(),
                     on=['key'], how='inner').set_index(['key','X','Y'])

.. ipython:: python
   :suppress:

   @savefig merging_merge_two_multiindex.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. _merging.merge_on_columns_and_levels:

Merging on a combination of columns and index levels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.22

Strings passed as the ``on``, ``left_on``, and ``right_on`` parameters
may refer to either column names or index level names.  This enables merging
``DataFrame`` instances on a combination of index levels and columns without
resetting indexes.

.. ipython:: python

   left_index = pd.Index(['K0', 'K0', 'K1', 'K2'], name='key1')

   left = pd.DataFrame({'A': ['A0', 'A1', 'A2', 'A3'],
                        'B': ['B0', 'B1', 'B2', 'B3'],
                        'key2': ['K0', 'K1', 'K0', 'K1']},
                       index=left_index)

   right_index = pd.Index(['K0', 'K1', 'K2', 'K2'], name='key1')

   right = pd.DataFrame({'C': ['C0', 'C1', 'C2', 'C3'],
                         'D': ['D0', 'D1', 'D2', 'D3'],
                         'key2': ['K0', 'K0', 'K0', 'K1']},
                        index=right_index)

   result = left.merge(right, on=['key1', 'key2'])

.. ipython:: python
   :suppress:

   @savefig merge_on_index_and_column.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. note::

   When DataFrames are merged on a string that matches an index level in both
   frames, the index level is preserved as an index level in the resulting
   DataFrame.

.. note::

   If a string matches both a column name and an index level name, then a
   warning is issued and the column takes precedence. This will result in an
   ambiguity error in a future version.

Overlapping value columns
~~~~~~~~~~~~~~~~~~~~~~~~~

The merge ``suffixes`` argument takes a tuple of list of strings to append to
overlapping column names in the input DataFrames to disambiguate the result
columns:

.. ipython:: python

   left = pd.DataFrame({'k': ['K0', 'K1', 'K2'], 'v': [1, 2, 3]})
   right = pd.DataFrame({'k': ['K0', 'K0', 'K3'], 'v': [4, 5, 6]})

   result = pd.merge(left, right, on='k')

.. ipython:: python
   :suppress:

   @savefig merging_merge_overlapped.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. ipython:: python

   result = pd.merge(left, right, on='k', suffixes=['_l', '_r'])

.. ipython:: python
   :suppress:

   @savefig merging_merge_overlapped_suffix.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

``DataFrame.join`` has ``lsuffix`` and ``rsuffix`` arguments which behave
similarly.

.. ipython:: python

   left = left.set_index('k')
   right = right.set_index('k')
   result = left.join(right, lsuffix='_l', rsuffix='_r')

.. ipython:: python
   :suppress:

   @savefig merging_merge_overlapped_multi_suffix.png
   p.plot([left, right], result,
          labels=['left', 'right'], vertical=False);
   plt.close('all');

.. _merging.multiple_join:

Joining multiple DataFrame or Panel objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A list or tuple of DataFrames can also be passed to ``DataFrame.join`` to join
them together on their indexes. The same is true for ``Panel.join``.

.. ipython:: python

   right2 = pd.DataFrame({'v': [7, 8, 9]}, index=['K1', 'K1', 'K2'])
   result = left.join([right, right2])

.. ipython:: python
   :suppress:

   @savefig merging_join_multi_df.png
   p.plot([left, right, right2], result,
          labels=['left', 'right', 'right2'], vertical=False);
   plt.close('all');

.. _merging.combine_first.update:

Merging together values within Series or DataFrame columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another fairly common situation is to have two like-indexed (or similarly
indexed) Series or DataFrame objects and wanting to "patch" values in one
object from values for matching indices in the other. Here is an example:

.. ipython:: python

   df1 = pd.DataFrame([[np.nan, 3., 5.], [-4.6, np.nan, np.nan],
                      [np.nan, 7., np.nan]])
   df2 = pd.DataFrame([[-42.6, np.nan, -8.2], [-5., 1.6, 4]],
                      index=[1, 2])

For this, use the ``combine_first`` method:

.. ipython:: python

   result = df1.combine_first(df2)

.. ipython:: python
   :suppress:

   @savefig merging_combine_first.png
   p.plot([df1, df2], result,
          labels=['df1', 'df2'], vertical=False);
   plt.close('all');

Note that this method only takes values from the right DataFrame if they are
missing in the left DataFrame. A related method, ``update``, alters non-NA
values inplace:

.. ipython:: python
   :suppress:

   df1_copy = df1.copy()

.. ipython:: python

   df1.update(df2)

.. ipython:: python
   :suppress:

   @savefig merging_update.png
   p.plot([df1_copy, df2], df1,
          labels=['df1', 'df2'], vertical=False);
   plt.close('all');

.. _merging.time_series:

Timeseries friendly merging
---------------------------

.. _merging.merge_ordered:

Merging Ordered Data
~~~~~~~~~~~~~~~~~~~~

A :func:`merge_ordered` function allows combining time series and other
ordered data. In particular it has an optional ``fill_method`` keyword to
fill/interpolate missing data:

.. ipython:: python

   left = pd.DataFrame({'k': ['K0', 'K1', 'K1', 'K2'],
                        'lv': [1, 2, 3, 4],
                        's': ['a', 'b', 'c', 'd']})

   right = pd.DataFrame({'k': ['K1', 'K2', 'K4'],
                         'rv': [1, 2, 3]})

   pd.merge_ordered(left, right, fill_method='ffill', left_by='s')

.. _merging.merge_asof:

Merging AsOf
~~~~~~~~~~~~

.. versionadded:: 0.19.0

A :func:`merge_asof` is similar to an ordered left-join except that we match on nearest key rather than equal keys. For each row in the ``left`` DataFrame, we select the last row in the ``right`` DataFrame whose ``on`` key is less than the left's key. Both DataFrames must be sorted by the key.

Optionally an asof merge can perform a group-wise merge. This matches the ``by`` key equally,
in addition to the nearest match on the ``on`` key.

For example; we might have ``trades`` and ``quotes`` and we want to ``asof`` merge them.

.. ipython:: python

   trades = pd.DataFrame({
       'time': pd.to_datetime(['20160525 13:30:00.023',
                               '20160525 13:30:00.038',
                               '20160525 13:30:00.048',
                               '20160525 13:30:00.048',
                               '20160525 13:30:00.048']),
       'ticker': ['MSFT', 'MSFT',
                  'GOOG', 'GOOG', 'AAPL'],
       'price': [51.95, 51.95,
                 720.77, 720.92, 98.00],
       'quantity': [75, 155,
                    100, 100, 100]},
       columns=['time', 'ticker', 'price', 'quantity'])

   quotes = pd.DataFrame({
       'time': pd.to_datetime(['20160525 13:30:00.023',
                               '20160525 13:30:00.023',
                               '20160525 13:30:00.030',
                               '20160525 13:30:00.041',
                               '20160525 13:30:00.048',
                               '20160525 13:30:00.049',
                               '20160525 13:30:00.072',
                               '20160525 13:30:00.075']),
       'ticker': ['GOOG', 'MSFT', 'MSFT',
                  'MSFT', 'GOOG', 'AAPL', 'GOOG',
                  'MSFT'],
       'bid': [720.50, 51.95, 51.97, 51.99,
               720.50, 97.99, 720.50, 52.01],
       'ask': [720.93, 51.96, 51.98, 52.00,
               720.93, 98.01, 720.88, 52.03]},
       columns=['time', 'ticker', 'bid', 'ask'])

.. ipython:: python

   trades
   quotes

By default we are taking the asof of the quotes.

.. ipython:: python

   pd.merge_asof(trades, quotes,
                 on='time',
                 by='ticker')

We only asof within ``2ms`` between the quote time and the trade time.

.. ipython:: python

   pd.merge_asof(trades, quotes,
                 on='time',
                 by='ticker',
                 tolerance=pd.Timedelta('2ms'))

We only asof within ``10ms`` between the quote time and the trade time and we exclude exact matches on time.
Note that though we exclude the exact matches (of the quotes), prior quotes DO propagate to that point
in time.

.. ipython:: python

   pd.merge_asof(trades, quotes,
                 on='time',
                 by='ticker',
                 tolerance=pd.Timedelta('10ms'),
                 allow_exact_matches=False)
