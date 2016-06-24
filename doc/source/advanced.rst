.. _advanced:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import pandas as pd
   np.random.seed(123456)
   np.set_printoptions(precision=4, suppress=True)
   pd.options.display.max_rows=15

******************************
MultiIndex / Advanced Indexing
******************************

This section covers indexing with a ``MultiIndex`` and more advanced indexing features.

See the :ref:`Indexing and Selecting Data <indexing>` for general indexing documentation.

.. warning::

   Whether a copy or a reference is returned for a setting operation, may
   depend on the context.  This is sometimes called ``chained assignment`` and
   should be avoided.  See :ref:`Returning a View versus Copy
   <indexing.view_versus_copy>`

.. warning::

   In 0.15.0 ``Index`` has internally been refactored to no longer sub-class ``ndarray``
   but instead subclass ``PandasObject``, similarly to the rest of the pandas objects. This should be
   a transparent change with only very limited API implications (See the :ref:`Internal Refactoring <whatsnew_0150.refactoring>`)

See the :ref:`cookbook<cookbook.selection>` for some advanced strategies

.. _advanced.hierarchical:

Hierarchical indexing (MultiIndex)
----------------------------------

Hierarchical / Multi-level indexing is very exciting as it opens the door to some
quite sophisticated data analysis and manipulation, especially for working with
higher dimensional data. In essence, it enables you to store and manipulate
data with an arbitrary number of dimensions in lower dimensional data
structures like Series (1d) and DataFrame (2d).

In this section, we will show what exactly we mean by "hierarchical" indexing
and how it integrates with the all of the pandas indexing functionality
described above and in prior sections. Later, when discussing :ref:`group by
<groupby>` and :ref:`pivoting and reshaping data <reshaping>`, we'll show
non-trivial applications to illustrate how it aids in structuring data for
analysis.

See the :ref:`cookbook<cookbook.multi_index>` for some advanced strategies

Creating a MultiIndex (hierarchical index) object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``MultiIndex`` object is the hierarchical analogue of the standard
``Index`` object which typically stores the axis labels in pandas objects. You
can think of ``MultiIndex`` an array of tuples where each tuple is unique. A
``MultiIndex`` can be created from a list of arrays (using
``MultiIndex.from_arrays``), an array of tuples (using
``MultiIndex.from_tuples``), or a crossed set of iterables (using
``MultiIndex.from_product``).  The ``Index`` constructor will attempt to return
a ``MultiIndex`` when it is passed a list of tuples.  The following examples
demo different ways to initialize MultiIndexes.


.. ipython:: python

   arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
             ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
   tuples = list(zip(*arrays))
   tuples

   index = pd.MultiIndex.from_tuples(tuples, names=['first', 'second'])
   index

   s = pd.Series(np.random.randn(8), index=index)
   s

When you want every pairing of the elements in two iterables, it can be easier
to use the ``MultiIndex.from_product`` function:

.. ipython:: python

   iterables = [['bar', 'baz', 'foo', 'qux'], ['one', 'two']]
   pd.MultiIndex.from_product(iterables, names=['first', 'second'])

As a convenience, you can pass a list of arrays directly into Series or
DataFrame to construct a MultiIndex automatically:

.. ipython:: python

   arrays = [np.array(['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux']),
             np.array(['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two'])]
   s = pd.Series(np.random.randn(8), index=arrays)
   s
   df = pd.DataFrame(np.random.randn(8, 4), index=arrays)
   df

All of the ``MultiIndex`` constructors accept a ``names`` argument which stores
string names for the levels themselves. If no names are provided, ``None`` will
be assigned:

.. ipython:: python

   df.index.names

This index can back any axis of a pandas object, and the number of **levels**
of the index is up to you:

.. ipython:: python

   df = pd.DataFrame(np.random.randn(3, 8), index=['A', 'B', 'C'], columns=index)
   df
   pd.DataFrame(np.random.randn(6, 6), index=index[:6], columns=index[:6])

We've "sparsified" the higher levels of the indexes to make the console output a
bit easier on the eyes.

It's worth keeping in mind that there's nothing preventing you from using
tuples as atomic labels on an axis:

.. ipython:: python

   pd.Series(np.random.randn(8), index=tuples)

The reason that the ``MultiIndex`` matters is that it can allow you to do
grouping, selection, and reshaping operations as we will describe below and in
subsequent areas of the documentation. As you will see in later sections, you
can find yourself working with hierarchically-indexed data without creating a
``MultiIndex`` explicitly yourself. However, when loading data from a file, you
may wish to generate your own ``MultiIndex`` when preparing the data set.

Note that how the index is displayed by be controlled using the
``multi_sparse`` option in ``pandas.set_printoptions``:

.. ipython:: python

   pd.set_option('display.multi_sparse', False)
   df
   pd.set_option('display.multi_sparse', True)

.. _advanced.get_level_values:

Reconstructing the level labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The method ``get_level_values`` will return a vector of the labels for each
location at a particular level:

.. ipython:: python

   index.get_level_values(0)
   index.get_level_values('second')

Basic indexing on axis with MultiIndex
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the important features of hierarchical indexing is that you can select
data by a "partial" label identifying a subgroup in the data. **Partial**
selection "drops" levels of the hierarchical index in the result in a
completely analogous way to selecting a column in a regular DataFrame:

.. ipython:: python

   df['bar']
   df['bar', 'one']
   df['bar']['one']
   s['qux']

See :ref:`Cross-section with hierarchical index <advanced.xs>` for how to select
on a deeper level.

.. note::

   The repr of a ``MultiIndex`` shows ALL the defined levels of an index, even
   if the they are not actually used. When slicing an index, you may notice this.
   For example:

   .. ipython:: python

      # original multi-index
      df.columns

      # sliced
      df[['foo','qux']].columns

   This is done to avoid a recomputation of the levels in order to make slicing
   highly performant. If you want to see the actual used levels.

   .. ipython:: python

      df[['foo','qux']].columns.values

      # for a specific level
      df[['foo','qux']].columns.get_level_values(0)

   To reconstruct the multiindex with only the used levels

   .. ipython:: python

      pd.MultiIndex.from_tuples(df[['foo','qux']].columns.values)

Data alignment and using ``reindex``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Operations between differently-indexed objects having ``MultiIndex`` on the
axes will work as you expect; data alignment will work the same as an Index of
tuples:

.. ipython:: python

   s + s[:-2]
   s + s[::2]

``reindex`` can be called with another ``MultiIndex`` or even a list or array
of tuples:

.. ipython:: python

   s.reindex(index[:3])
   s.reindex([('foo', 'two'), ('bar', 'one'), ('qux', 'one'), ('baz', 'one')])

.. _advanced.advanced_hierarchical:

Advanced indexing with hierarchical index
-----------------------------------------

Syntactically integrating ``MultiIndex`` in advanced indexing with ``.loc/.ix`` is a
bit challenging, but we've made every effort to do so. for example the
following works as you would expect:

.. ipython:: python

   df = df.T
   df
   df.loc['bar']
   df.loc['bar', 'two']

"Partial" slicing also works quite nicely.

.. ipython:: python

   df.loc['baz':'foo']

You can slice with a 'range' of values, by providing a slice of tuples.

.. ipython:: python

   df.loc[('baz', 'two'):('qux', 'one')]
   df.loc[('baz', 'two'):'foo']

Passing a list of labels or tuples works similar to reindexing:

.. ipython:: python

   df.ix[[('bar', 'two'), ('qux', 'one')]]

.. _advanced.mi_slicers:

Using slicers
~~~~~~~~~~~~~

.. versionadded:: 0.14.0

In 0.14.0 we added a new way to slice multi-indexed objects.
You can slice a multi-index by providing multiple indexers.

You can provide any of the selectors as if you are indexing by label, see :ref:`Selection by Label <indexing.label>`,
including slices, lists of labels, labels, and boolean indexers.

You can use ``slice(None)`` to select all the contents of *that* level. You do not need to specify all the
*deeper* levels, they will be implied as ``slice(None)``.

As usual, **both sides** of the slicers are included as this is label indexing.

.. warning::

   You should specify all axes in the ``.loc`` specifier, meaning the indexer for the **index** and
   for the **columns**. There are some ambiguous cases where the passed indexer could be mis-interpreted
   as indexing *both* axes, rather than into say the MuliIndex for the rows.

   You should do this:

   .. code-block:: python

      df.loc[(slice('A1','A3'),.....),:]

   rather than this:

   .. code-block:: python

      df.loc[(slice('A1','A3'),.....)]

.. ipython:: python

   def mklbl(prefix,n):
       return ["%s%s" % (prefix,i)  for i in range(n)]

   miindex = pd.MultiIndex.from_product([mklbl('A',4),
                                         mklbl('B',2),
                                         mklbl('C',4),
                                         mklbl('D',2)])
   micolumns = pd.MultiIndex.from_tuples([('a','foo'),('a','bar'),
                                          ('b','foo'),('b','bah')],
                                         names=['lvl0', 'lvl1'])
   dfmi = pd.DataFrame(np.arange(len(miindex)*len(micolumns)).reshape((len(miindex),len(micolumns))),
                       index=miindex,
                       columns=micolumns).sort_index().sort_index(axis=1)
   dfmi

Basic multi-index slicing using slices, lists, and labels.

.. ipython:: python

   dfmi.loc[(slice('A1','A3'),slice(None), ['C1','C3']),:]

You can use a ``pd.IndexSlice`` to have a more natural syntax using ``:`` rather than using ``slice(None)``

.. ipython:: python

   idx = pd.IndexSlice
   dfmi.loc[idx[:,:,['C1','C3']],idx[:,'foo']]

It is possible to perform quite complicated selections using this method on multiple
axes at the same time.

.. ipython:: python

   dfmi.loc['A1',(slice(None),'foo')]
   dfmi.loc[idx[:,:,['C1','C3']],idx[:,'foo']]

Using a boolean indexer you can provide selection related to the *values*.

.. ipython:: python

   mask = dfmi[('a','foo')]>200
   dfmi.loc[idx[mask,:,['C1','C3']],idx[:,'foo']]

You can also specify the ``axis`` argument to ``.loc`` to interpret the passed
slicers on a single axis.

.. ipython:: python

   dfmi.loc(axis=0)[:,:,['C1','C3']]

Furthermore you can *set* the values using these methods

.. ipython:: python

   df2 = dfmi.copy()
   df2.loc(axis=0)[:,:,['C1','C3']] = -10
   df2

You can use a right-hand-side of an alignable object as well.

.. ipython:: python

   df2 = dfmi.copy()
   df2.loc[idx[:,:,['C1','C3']],:] = df2*1000
   df2

.. _advanced.xs:

Cross-section
~~~~~~~~~~~~~

The ``xs`` method of ``DataFrame`` additionally takes a level argument to make
selecting data at a particular level of a MultiIndex easier.

.. ipython:: python

   df
   df.xs('one', level='second')

.. ipython:: python

   # using the slicers (new in 0.14.0)
   df.loc[(slice(None),'one'),:]

You can also select on the columns with :meth:`~pandas.MultiIndex.xs`, by
providing the axis argument

.. ipython:: python

   df = df.T
   df.xs('one', level='second', axis=1)

.. ipython:: python

   # using the slicers (new in 0.14.0)
   df.loc[:,(slice(None),'one')]

:meth:`~pandas.MultiIndex.xs` also allows selection with multiple keys

.. ipython:: python

   df.xs(('one', 'bar'), level=('second', 'first'), axis=1)

.. ipython:: python

   # using the slicers (new in 0.14.0)
   df.loc[:,('bar','one')]

.. versionadded:: 0.13.0

You can pass ``drop_level=False`` to :meth:`~pandas.MultiIndex.xs` to retain
the level that was selected

.. ipython:: python

   df.xs('one', level='second', axis=1, drop_level=False)

versus the result with ``drop_level=True`` (the default value)

.. ipython:: python

   df.xs('one', level='second', axis=1, drop_level=True)

.. ipython:: python
   :suppress:

   df = df.T

.. _advanced.advanced_reindex:

Advanced reindexing and alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The parameter ``level`` has been added to the ``reindex`` and ``align`` methods
of pandas objects. This is useful to broadcast values across a level. For
instance:

.. ipython:: python

   midx = pd.MultiIndex(levels=[['zero', 'one'], ['x','y']],
                        labels=[[1,1,0,0],[1,0,1,0]])
   df = pd.DataFrame(np.random.randn(4,2), index=midx)
   df
   df2 = df.mean(level=0)
   df2
   df2.reindex(df.index, level=0)

   # aligning
   df_aligned, df2_aligned = df.align(df2, level=0)
   df_aligned
   df2_aligned


Swapping levels with :meth:`~pandas.MultiIndex.swaplevel`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``swaplevel`` function can switch the order of two levels:

.. ipython:: python

   df[:5]
   df[:5].swaplevel(0, 1, axis=0)

.. _advanced.reorderlevels:

Reordering levels with :meth:`~pandas.MultiIndex.reorder_levels`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``reorder_levels`` function generalizes the ``swaplevel`` function,
allowing you to permute the hierarchical index levels in one step:

.. ipython:: python

   df[:5].reorder_levels([1,0], axis=0)

Sorting a :class:`~pandas.MultiIndex`
-------------------------------------

For MultiIndex-ed objects to be indexed & sliced effectively, they need
to be sorted. As with any index, you can use ``sort_index``.

.. ipython:: python

   import random; random.shuffle(tuples)
   s = pd.Series(np.random.randn(8), index=pd.MultiIndex.from_tuples(tuples))
   s
   s.sort_index()
   s.sort_index(level=0)
   s.sort_index(level=1)

.. _advanced.sortlevel_byname:

You may also pass a level name to ``sort_index`` if the MultiIndex levels
are named.

.. ipython:: python

   s.index.set_names(['L1', 'L2'], inplace=True)
   s.sort_index(level='L1')
   s.sort_index(level='L2')

On higher dimensional objects, you can sort any of the other axes by level if
they have a MultiIndex:

.. ipython:: python

   df.T.sort_index(level=1, axis=1)

Indexing will work even if the data are not sorted, but will be rather
inefficient (and show a ``PerformanceWarning``). It will also
return a copy of the data rather than a view:

.. ipython:: python

   dfm = pd.DataFrame({'jim': [0, 0, 1, 1],
                       'joe': ['x', 'x', 'z', 'y'],
                       'jolie': np.random.rand(4)})
   dfm = dfm.set_index(['jim', 'joe'])
   dfm

.. code-block:: ipython

   In [4]: dfm.loc[(1, 'z')]
   PerformanceWarning: indexing past lexsort depth may impact performance.

   Out[4]:
              jolie
   jim joe
   1   z    0.64094

Furthermore if you try to index something that is not fully lexsorted, this can raise:

.. code-block:: ipython

    In [5]: dfm.loc[(0,'y'):(1, 'z')]
    KeyError: 'Key length (2) was greater than MultiIndex lexsort depth (1)'

The ``is_lexsorted()`` method on an ``Index`` show if the index is sorted, and the ``lexsort_depth`` property returns the sort depth:

.. ipython:: python

   dfm.index.is_lexsorted()
   dfm.index.lexsort_depth

.. ipython:: python

   dfm = dfm.sort_index()
   dfm
   dfm.index.is_lexsorted()
   dfm.index.lexsort_depth

And now selection works as expected.

.. ipython:: python

   dfm.loc[(0,'y'):(1, 'z')]

Take Methods
------------

.. _advanced.take:

Similar to numpy ndarrays, pandas Index, Series, and DataFrame also provides
the ``take`` method that retrieves elements along a given axis at the given
indices. The given indices must be either a list or an ndarray of integer
index positions. ``take`` will also accept negative integers as relative positions to the end of the object.

.. ipython:: python

   index = pd.Index(np.random.randint(0, 1000, 10))
   index

   positions = [0, 9, 3]

   index[positions]
   index.take(positions)

   ser = pd.Series(np.random.randn(10))

   ser.iloc[positions]
   ser.take(positions)

For DataFrames, the given indices should be a 1d list or ndarray that specifies
row or column positions.

.. ipython:: python

   frm = pd.DataFrame(np.random.randn(5, 3))

   frm.take([1, 4, 3])

   frm.take([0, 2], axis=1)

It is important to note that the ``take`` method on pandas objects are not
intended to work on boolean indices and may return unexpected results.

.. ipython:: python

   arr = np.random.randn(10)
   arr.take([False, False, True, True])
   arr[[0, 1]]

   ser = pd.Series(np.random.randn(10))
   ser.take([False, False, True, True])
   ser.ix[[0, 1]]

Finally, as a small note on performance, because the ``take`` method handles
a narrower range of inputs, it can offer performance that is a good deal
faster than fancy indexing.

.. ipython::

   arr = np.random.randn(10000, 5)
   indexer = np.arange(10000)
   random.shuffle(indexer)

   timeit arr[indexer]
   timeit arr.take(indexer, axis=0)

   ser = pd.Series(arr[:, 0])
   timeit ser.ix[indexer]
   timeit ser.take(indexer)

.. _indexing.index_types:

Index Types
-----------

We have discussed ``MultiIndex`` in the previous sections pretty extensively. ``DatetimeIndex`` and ``PeriodIndex``
are shown :ref:`here <timeseries.overview>`. ``TimedeltaIndex`` are :ref:`here <timedeltas.timedeltas>`.

In the following sub-sections we will highlite some other index types.

.. _indexing.categoricalindex:

CategoricalIndex
~~~~~~~~~~~~~~~~

.. versionadded:: 0.16.1

We introduce a ``CategoricalIndex``, a new type of index object that is useful for supporting
indexing with duplicates. This is a container around a ``Categorical`` (introduced in v0.15.0)
and allows efficient indexing and storage of an index with a large number of duplicated elements. Prior to 0.16.1,
setting the index of a ``DataFrame/Series`` with a ``category`` dtype would convert this to regular object-based ``Index``.

.. ipython:: python

   df = pd.DataFrame({'A': np.arange(6),
                      'B': list('aabbca')})
   df['B'] = df['B'].astype('category', categories=list('cab'))
   df
   df.dtypes
   df.B.cat.categories

Setting the index, will create create a ``CategoricalIndex``

.. ipython:: python

   df2 = df.set_index('B')
   df2.index

Indexing with ``__getitem__/.iloc/.loc/.ix`` works similarly to an ``Index`` with duplicates.
The indexers MUST be in the category or the operation will raise.

.. ipython:: python

   df2.loc['a']

These PRESERVE the ``CategoricalIndex``

.. ipython:: python

   df2.loc['a'].index

Sorting will order by the order of the categories

.. ipython:: python

   df2.sort_index()

Groupby operations on the index will preserve the index nature as well

.. ipython:: python

   df2.groupby(level=0).sum()
   df2.groupby(level=0).sum().index

Reindexing operations, will return a resulting index based on the type of the passed
indexer, meaning that passing a list will return a plain-old-``Index``; indexing with
a ``Categorical`` will return a ``CategoricalIndex``, indexed according to the categories
of the PASSED ``Categorical`` dtype. This allows one to arbitrarly index these even with
values NOT in the categories, similarly to how you can reindex ANY pandas index.

.. ipython :: python

   df2.reindex(['a','e'])
   df2.reindex(['a','e']).index
   df2.reindex(pd.Categorical(['a','e'],categories=list('abcde')))
   df2.reindex(pd.Categorical(['a','e'],categories=list('abcde'))).index

.. warning::

   Reshaping and Comparison operations on a ``CategoricalIndex`` must have the same categories
   or a ``TypeError`` will be raised.

   .. code-block:: python

      In [9]: df3 = pd.DataFrame({'A' : np.arange(6),
                                  'B' : pd.Series(list('aabbca')).astype('category')})

      In [11]: df3 = df3.set_index('B')

      In [11]: df3.index
      Out[11]: CategoricalIndex([u'a', u'a', u'b', u'b', u'c', u'a'], categories=[u'a', u'b', u'c'], ordered=False, name=u'B', dtype='category')

      In [12]: pd.concat([df2, df3]
      TypeError: categories must match existing categories when appending

.. _indexing.rangeindex:

Int64Index and RangeIndex
~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   Indexing on an integer-based Index with floats has been clarified in 0.18.0, for a summary of the changes, see :ref:`here <whatsnew_0180.float_indexers>`.

``Int64Index`` is a fundamental basic index in *pandas*. This is an Immutable array implementing an ordered, sliceable set.
Prior to 0.18.0, the ``Int64Index`` would provide the default index for all ``NDFrame`` objects.

``RangeIndex`` is a sub-class of ``Int64Index`` added in version 0.18.0, now providing the default index for all ``NDFrame`` objects.
``RangeIndex`` is an optimized version of ``Int64Index`` that can represent a monotonic ordered set. These are analagous to python `range types <https://docs.python.org/3/library/stdtypes.html#typesseq-range>`__.

.. _indexing.float64index:

Float64Index
~~~~~~~~~~~~

.. note::

   As of 0.14.0, ``Float64Index`` is backed by a native ``float64`` dtype
   array. Prior to 0.14.0, ``Float64Index`` was backed by an ``object`` dtype
   array. Using a ``float64`` dtype in the backend speeds up arithmetic
   operations by about 30x and boolean indexing operations on the
   ``Float64Index`` itself are about 2x as fast.

.. versionadded:: 0.13.0

By default a ``Float64Index`` will be automatically created when passing floating, or mixed-integer-floating values in index creation.
This enables a pure label-based slicing paradigm that makes ``[],ix,loc`` for scalar indexing and slicing work exactly the
same.

.. ipython:: python

   indexf = pd.Index([1.5, 2, 3, 4.5, 5])
   indexf
   sf = pd.Series(range(5), index=indexf)
   sf

Scalar selection for ``[],.ix,.loc`` will always be label based. An integer will match an equal float index (e.g. ``3`` is equivalent to ``3.0``)

.. ipython:: python

   sf[3]
   sf[3.0]
   sf.ix[3]
   sf.ix[3.0]
   sf.loc[3]
   sf.loc[3.0]

The only positional indexing is via ``iloc``

.. ipython:: python

   sf.iloc[3]

A scalar index that is not found will raise ``KeyError``

Slicing is ALWAYS on the values of the index, for ``[],ix,loc`` and ALWAYS positional with ``iloc``

.. ipython:: python

   sf[2:4]
   sf.ix[2:4]
   sf.loc[2:4]
   sf.iloc[2:4]

In float indexes, slicing using floats is allowed

.. ipython:: python

   sf[2.1:4.6]
   sf.loc[2.1:4.6]

In non-float indexes, slicing using floats will raise a ``TypeError``

.. code-block:: ipython

   In [1]: pd.Series(range(5))[3.5]
   TypeError: the label [3.5] is not a proper indexer for this index type (Int64Index)

   In [1]: pd.Series(range(5))[3.5:4.5]
   TypeError: the slice start [3.5] is not a proper indexer for this index type (Int64Index)

.. warning::

   Using a scalar float indexer for ``.iloc`` has been removed in 0.18.0, so the following will raise a ``TypeError``

   .. code-block:: ipython

      In [3]: pd.Series(range(5)).iloc[3.0]
      TypeError: cannot do positional indexing on <class 'pandas.indexes.range.RangeIndex'> with these indexers [3.0] of <type 'float'>

   Further the treatment of ``.ix`` with a float indexer on a non-float index, will be label based, and thus coerce the index.

   .. ipython:: python

      s2 = pd.Series([1, 2, 3], index=list('abc'))
      s2
      s2.ix[1.0] = 10
      s2

Here is a typical use-case for using this type of indexing. Imagine that you have a somewhat
irregular timedelta-like indexing scheme, but the data is recorded as floats. This could for
example be millisecond offsets.

.. ipython:: python

   dfir = pd.concat([pd.DataFrame(np.random.randn(5,2),
                                  index=np.arange(5) * 250.0,
                                  columns=list('AB')),
                     pd.DataFrame(np.random.randn(6,2),
                                  index=np.arange(4,10) * 250.1,
                                  columns=list('AB'))])
   dfir

Selection operations then will always work on a value basis, for all selection operators.

.. ipython:: python

   dfir[0:1000.4]
   dfir.loc[0:1001,'A']
   dfir.loc[1000.4]

You could then easily pick out the first 1 second (1000 ms) of data then.

.. ipython:: python

   dfir[0:1000]

Of course if you need integer based selection, then use ``iloc``

.. ipython:: python

   dfir.iloc[0:5]
