.. _advanced:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import random
   np.random.seed(123456)
   from pandas import *
   options.display.max_rows=15
   import pandas as pd
   randn = np.random.randn
   randint = np.random.randint
   np.set_printoptions(precision=4, suppress=True)
   from pandas.compat import range, zip

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

   index = MultiIndex.from_tuples(tuples, names=['first', 'second'])
   index

   s = Series(randn(8), index=index)
   s

When you want every pairing of the elements in two iterables, it can be easier
to use the ``MultiIndex.from_product`` function:

.. ipython:: python

   iterables = [['bar', 'baz', 'foo', 'qux'], ['one', 'two']]
   MultiIndex.from_product(iterables, names=['first', 'second'])

As a convenience, you can pass a list of arrays directly into Series or
DataFrame to construct a MultiIndex automatically:

.. ipython:: python

   arrays = [np.array(['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux']),
             np.array(['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two'])]
   s = Series(randn(8), index=arrays)
   s
   df = DataFrame(randn(8, 4), index=arrays)
   df

All of the ``MultiIndex`` constructors accept a ``names`` argument which stores
string names for the levels themselves. If no names are provided, ``None`` will
be assigned:

.. ipython:: python

   df.index.names

This index can back any axis of a pandas object, and the number of **levels**
of the index is up to you:

.. ipython:: python

   df = DataFrame(randn(3, 8), index=['A', 'B', 'C'], columns=index)
   df
   DataFrame(randn(6, 6), index=index[:6], columns=index[:6])

We've "sparsified" the higher levels of the indexes to make the console output a
bit easier on the eyes.

It's worth keeping in mind that there's nothing preventing you from using
tuples as atomic labels on an axis:

.. ipython:: python

   Series(randn(8), index=tuples)

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
   for the **columns**. Their are some ambiguous cases where the passed indexer could be mis-interpreted
   as indexing *both* axes, rather than into say the MuliIndex for the rows.

   You should do this:

   .. code-block:: python

      df.loc[(slice('A1','A3'),.....),:]

   rather than this:

   .. code-block:: python

      df.loc[(slice('A1','A3'),.....)]

.. warning::

   You will need to make sure that the selection axes are fully lexsorted!

.. ipython:: python

   def mklbl(prefix,n):
       return ["%s%s" % (prefix,i)  for i in range(n)]

   miindex = MultiIndex.from_product([mklbl('A',4),
                                      mklbl('B',2),
                                      mklbl('C',4),
                                      mklbl('D',2)])
   micolumns = MultiIndex.from_tuples([('a','foo'),('a','bar'),
                                       ('b','foo'),('b','bah')],
                                        names=['lvl0', 'lvl1'])
   dfmi = DataFrame(np.arange(len(miindex)*len(micolumns)).reshape((len(miindex),len(micolumns))),
                    index=miindex,
                    columns=micolumns).sortlevel().sortlevel(axis=1)
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

   midx = MultiIndex(levels=[['zero', 'one'], ['x','y']],
                     labels=[[1,1,0,0],[1,0,1,0]])
   df = DataFrame(randn(4,2), index=midx)
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

The need for sortedness with :class:`~pandas.MultiIndex`
--------------------------------------------------------

**Caveat emptor**: the present implementation of ``MultiIndex`` requires that
the labels be sorted for some of the slicing / indexing routines to work
correctly. You can think about breaking the axis into unique groups, where at
the hierarchical level of interest, each distinct group shares a label, but no
two have the same label. However, the ``MultiIndex`` does not enforce this:
**you are responsible for ensuring that things are properly sorted**. There is
an important new method ``sortlevel`` to sort an axis within a ``MultiIndex``
so that its labels are grouped and sorted by the original ordering of the
associated factor at that level. Note that this does not necessarily mean the
labels will be sorted lexicographically!

.. ipython:: python

   import random; random.shuffle(tuples)
   s = Series(randn(8), index=MultiIndex.from_tuples(tuples))
   s
   s.sortlevel(0)
   s.sortlevel(1)

.. _advanced.sortlevel_byname:

Note, you may also pass a level name to ``sortlevel`` if the MultiIndex levels
are named.

.. ipython:: python

   s.index.set_names(['L1', 'L2'], inplace=True)
   s.sortlevel(level='L1')
   s.sortlevel(level='L2')

Some indexing will work even if the data are not sorted, but will be rather
inefficient and will also return a copy of the data rather than a view:

.. ipython:: python

   s['qux']
   s.sortlevel(1)['qux']

On higher dimensional objects, you can sort any of the other axes by level if
they have a MultiIndex:

.. ipython:: python

   df.T.sortlevel(1, axis=1)

The ``MultiIndex`` object has code to **explicity check the sort depth**. Thus,
if you try to index at a depth at which the index is not sorted, it will raise
an exception. Here is a concrete example to illustrate this:

.. ipython:: python

   tuples = [('a', 'a'), ('a', 'b'), ('b', 'a'), ('b', 'b')]
   idx = MultiIndex.from_tuples(tuples)
   idx.lexsort_depth

   reordered = idx[[1, 0, 3, 2]]
   reordered.lexsort_depth

   s = Series(randn(4), index=reordered)
   s.ix['a':'a']

However:

::

   >>> s.ix[('a', 'b'):('b', 'a')]
   Traceback (most recent call last)
        ...
   KeyError: Key length (3) was greater than MultiIndex lexsort depth (2)


Take Methods
------------

.. _advanced.take:

Similar to numpy ndarrays, pandas Index, Series, and DataFrame also provides
the ``take`` method that retrieves elements along a given axis at the given
indices. The given indices must be either a list or an ndarray of integer
index positions. ``take`` will also accept negative integers as relative positions to the end of the object.

.. ipython:: python

   index = Index(randint(0, 1000, 10))
   index

   positions = [0, 9, 3]

   index[positions]
   index.take(positions)

   ser = Series(randn(10))

   ser.iloc[positions]
   ser.take(positions)

For DataFrames, the given indices should be a 1d list or ndarray that specifies
row or column positions.

.. ipython:: python

   frm = DataFrame(randn(5, 3))

   frm.take([1, 4, 3])

   frm.take([0, 2], axis=1)

It is important to note that the ``take`` method on pandas objects are not
intended to work on boolean indices and may return unexpected results.

.. ipython:: python

   arr = randn(10)
   arr.take([False, False, True, True])
   arr[[0, 1]]

   ser = Series(randn(10))
   ser.take([False, False, True, True])
   ser.ix[[0, 1]]

Finally, as a small note on performance, because the ``take`` method handles
a narrower range of inputs, it can offer performance that is a good deal
faster than fancy indexing.

.. ipython::

   arr = randn(10000, 5)
   indexer = np.arange(10000)
   random.shuffle(indexer)

   timeit arr[indexer]
   timeit arr.take(indexer, axis=0)

   ser = Series(arr[:, 0])
   timeit ser.ix[indexer]
   timeit ser.take(indexer)

.. _indexing.float64index:

Float64Index
------------

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

   indexf = Index([1.5, 2, 3, 4.5, 5])
   indexf
   sf = Series(range(5),index=indexf)
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

.. code-block:: python

   In [1]: Series(range(5))[3.5]
   TypeError: the label [3.5] is not a proper indexer for this index type (Int64Index)

   In [1]: Series(range(5))[3.5:4.5]
   TypeError: the slice start [3.5] is not a proper indexer for this index type (Int64Index)

Using a scalar float indexer will be deprecated in a future version, but is allowed for now.

.. code-block:: python

   In [3]: Series(range(5))[3.0]
   Out[3]: 3

Here is a typical use-case for using this type of indexing. Imagine that you have a somewhat
irregular timedelta-like indexing scheme, but the data is recorded as floats. This could for
example be millisecond offsets.

.. ipython:: python

   dfir = concat([DataFrame(randn(5,2),
                     index=np.arange(5) * 250.0,
                     columns=list('AB')),
                  DataFrame(randn(6,2),
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

