.. _indexing:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

***************************
Indexing and selecting data
***************************

The axis labeling information in pandas objects serves many purposes:

  - Identifies data (i.e. provides *metadata*) using known indicators,
    important for for analysis, visualization, and interactive console display
  - Enables automatic and explicit data alignment
  - Allows intuitive getting and setting of subsets of the data set

In this section / chapter, we will focus on the final point: namely, how to
slice, dice, and generally get and set subsets of pandas objects. The primary
focus will be on Series and DataFrame as they have received more development
attention in this area. Expect more work to be invested higher-dimensional data
structures (including Panel) in the future, especially in label-based advanced
indexing.

.. _indexing.basics:

Basics
------

As mentioned when introducing the data structures in the :ref:`last section
<basics>`, the primary function of indexing with ``[]`` (a.k.a. ``__getitem__``
for those familiar with implementing class behavior in Python) is selecting out
lower-dimensional slices. Thus,

  - **Series**: ``series[label]`` returns a scalar value
  - **DataFrame**: ``frame[colname]`` returns a Series corresponding to the
    passed column name
  - **Panel**: ``panel[itemname]`` returns a DataFrame corresponding to the
    passed item name

Here we construct a simple time series data set to use for illustrating the
indexing functionality:

.. ipython:: python

   dates = np.asarray(DateRange('1/1/2000', periods=8))
   df = DataFrame(randn(8, 4), index=dates, columns=['A', 'B', 'C', 'D'])
   df
   panel = Panel({'one' : df, 'two' : df - df.mean()})
   panel

.. note::

   None of the indexing functionality is time series specific unless
   specifically stated.

Thus, as per above, we have the most basic indexing using ``[]``:

.. ipython:: python

   s = df['A']
   s[dates[5]]
   panel['two']


.. _indexing.basics.get_value:

Fast scalar value getting and setting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since indexing with ``[]`` must handle a lot of cases (single-label access,
slicing, boolean indexing, etc.), it has a bit of overhead in order to figure
out what you're asking for. If you only want to access a scalar value, the
fastest way is to use the ``get_value`` method, which is implemented on all of
the data structures:

.. ipython:: python

   s.get_value(dates[5])
   df.get_value(dates[5], 'A')

There is an analogous ``set_value`` method which has the additional capability
of enlarging an object. This method *always* returns a reference to the object
it modified, which in the fast of enlargement, will be a **new object**:

.. ipython:: python

   df.set_value(dates[5], 'E', 7)

Additional Column Access
~~~~~~~~~~~~~~~~~~~~~~~~

.. _indexing.columns.multiple:

.. _indexing.df_cols:

You may access a column on a dataframe directly as an attribute:

.. ipython:: python

   df.A

If you are using the IPython environment, you may also use tab-completion to
see the accessible columns of a DataFrame.

You can pass a list of columns to ``[]`` to select columns in that order:
If a column is not contained in the DataFrame, an exception will be
raised. Multiple columns can also be set in this manner:

.. ipython:: python

   df
   df[['B', 'A']] = df[['A', 'B']]
   df

You may find this useful for applying a transform (in-place) to a subset of the
columns.

Data slices on other axes
~~~~~~~~~~~~~~~~~~~~~~~~~

It's certainly possible to retrieve data slices along the other axes of a
DataFrame or Panel. We tend to refer to these slices as
*cross-sections*. DataFrame has the ``xs`` function for retrieving rows as
Series and Panel has the analogous ``major_xs`` and ``minor_xs`` functions for
retrieving slices as DataFrames for a given ``major_axis`` or ``minor_axis``
label, respectively.

.. ipython:: python

   date = dates[5]
   df.xs(date)
   panel.major_xs(date)
   panel.minor_xs('A')


Slicing ranges
~~~~~~~~~~~~~~

The most robust and consistent way of slicing ranges along arbitrary axes is
described in the :ref:`Advanced indexing <indexing.advanced>` section detailing
the ``.ix`` method. For now, we explain the semantics of slicing using the
``[]`` operator.

With Series, the syntax works exactly as with an ndarray, returning a slice of
the values and the corresponding labels:

.. ipython:: python

   s[:5]
   s[::2]
   s[::-1]

Note that setting works as well:

.. ipython:: python

   s2 = s.copy()
   s2[:5] = 0
   s2

With DataFrame, slicing inside of ``[]`` **slices the rows**. This is provided
largely as a convenience since it is such a common operation.

.. ipython:: python

   df[:3]
   df[::-1]

Boolean indexing
~~~~~~~~~~~~~~~~

.. _indexing.boolean:

Using a boolean vector to index a Series works exactly as in a numpy ndarray:

.. ipython:: python

   s[s > 0]
   s[(s < 0) & (s > -0.5)]

You may select rows from a DataFrame using a boolean vector the same length as
the DataFrame's index (for example, something derived from one of the columns
of the DataFrame):

.. ipython:: python

   df[df['A'] > 0]

Consider the ``isin`` method of Series, which returns a boolean vector that is
true wherever the Series elements exist in the passed list. This allows you to
select out rows where one or more columns have values you want:

.. ipython:: python

   df2 = DataFrame({'a' : ['one', 'one', 'two', 'three', 'two', 'one', 'six'],
                    'b' : ['x', 'y', 'y', 'x', 'y', 'x', 'x'],
                    'c' : np.random.randn(7)})
   df2[df2['a'].isin(['one', 'two'])]

Note, with the :ref:`advanced indexing <indexing.advanced>` ``ix`` method, you
may select along more than one axis using boolean vectors combined with other
indexing expressions.

Indexing a DataFrame with a boolean DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may wish to set values on a DataFrame based on some boolean criteria
derived from itself or another DataFrame or set of DataFrames. This can be done
intuitively like so:

.. ipython:: python

   df2 = df.copy()
   df2 < 0
   df2[df2 < 0] = 0
   df2

Note that such an operation requires that the boolean DataFrame is indexed
exactly the same.


Take Methods
~~~~~~~~~~~~

.. _indexing.take:

TODO: Fill Me In

Duplicate Data
~~~~~~~~~~~~~~

.. _indexing.duplicate:

If you want to indentify and remove duplicate rows in a DataFrame,  there are
two methods that will help: ``duplicated`` and ``drop_duplicates``. Each
takes as an argument the columns to use to identify duplicated rows.

``duplicated`` returns a boolean vector whose length is the number of rows, and
which indicates whether a row is duplicated.

``drop_duplicates`` removes duplicate rows.

By default, the first observed row of a duplicate set is considered unique, but
each method has a ``take_last`` parameter that indicates the last observed row
should be taken instead.

.. ipython:: python

   df2 = DataFrame({'a' : ['one', 'one', 'two', 'three', 'two', 'one', 'six'],
                    'b' : ['x', 'y', 'y', 'x', 'y', 'x', 'x'],
                    'c' : np.random.randn(7)})
   df2.duplicated(['a','b'])
   df2.drop_duplicates(['a','b'])
   df2.drop_duplicates(['a','b'], take_last=True)

.. _indexing.dictionarylike:

Dictionary-like ``get`` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each of Series, DataFrame, and Panel have a ``get`` method which can return a
default value.

.. ipython:: python

   s = Series([1,2,3], index=['a','b','c'])
   s.get('a')               # equivalent to s['a']
   s.get('x', default=-1)


.. _indexing.advanced:

Advanced indexing with labels
-----------------------------

We have avoided excessively overloading the ``[]`` / ``__getitem__`` operator
to keep the basic functionality of the pandas objects straightforward and
simple. However, there are often times when you may wish get a subset (or
analogously set a subset) of the data in a way that is not straightforward
using the combination of ``reindex`` and ``[]``. Complicated setting operations
are actually quite difficult because ``reindex`` usually returns a copy.

By *advanced* indexing we are referring to a special ``.ix`` attribute on
pandas objects which enable you to do getting/setting operations on a
DataFrame, for example, with matrix/ndarray-like semantics. Thus you can
combine the following kinds of indexing:

  - An integer or single label, e.g. ``5`` or ``'a'``
  - A list or array of labels ``['a', 'b', 'c']`` or integers ``[4, 3, 0]``
  - A slice object with ints ``1:7`` or labels ``'a':'f'``
  - A boolean array

We'll illustrate all of these methods. First, note that this provides a concise
way of reindexing on multiple axes at once:

.. ipython:: python

   subindex = dates[[3,4,5]]
   df.reindex(index=subindex, columns=['C', 'B'])
   df.ix[subindex, ['C', 'B']]

Assignment / setting values is possible when using ``ix``:

.. ipython:: python

   df2 = df.copy()
   df2.ix[subindex, ['C', 'B']] = 0
   df2

Indexing with an array of integers can also be done:

.. ipython:: python

   df.ix[[4,3,1]]
   df.ix[dates[[4,3,1]]]

**Slicing** has standard Python semantics for integer slices:

.. ipython:: python

   df.ix[1:7, :2]

Slicing with labels is semantically slightly different because the slice start
and stop are **inclusive** in the label-based case:

.. ipython:: python

   date1, date2 = dates[[2, 4]]
   print date1, date2
   df.ix[date1:date2]
   df['A'].ix[date1:date2]

Getting and setting rows in a DataFrame, especially by their location, is much
easier:

.. ipython:: python

   df2 = df[:5].copy()
   df2.ix[3]
   df2.ix[3] = np.arange(len(df2.columns))
   df2

Column or row selection can be combined as you would expect with arrays of
labels or even boolean vectors:

.. ipython:: python

   df.ix[df['A'] > 0, 'B']
   df.ix[date1:date2, 'B']
   df.ix[date1, 'B']

Slicing with labels is closely related to the ``truncate`` method which does
precisely ``.ix[start:stop]`` but returns a copy (for legacy reasons).

Returning a view versus a copy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The rules about when a view on the data is returned are entirely dependent on
NumPy. Whenever an array of labels or a boolean vector are involved in the
indexing operation, the result will be a copy. With single label / scalar
indexing and slicing, e.g. ``df.ix[3:6]`` or ``df.ix[:, 'A']``, a view will be
returned.

The ``select`` method
~~~~~~~~~~~~~~~~~~~~~

Another way to extract slices from an object is with the ``select`` method of
Series, DataFrame, and Panel. This method should be used only when there is no
more direct way.  ``select`` takes a function which operates on labels along
``axis`` and returns a boolean.  For instance:

.. ipython:: python

   df.select(lambda x: x == 'A', axis=1)

The ``lookup`` method
~~~~~~~~~~~~~~~~~~~~~

Sometimes you want to extract a set of values given a sequence of row labels
and column labels, and the ``lookup`` method allows for this and returns a
numpy array.  For instance,

.. ipython:: python

  dflookup = DataFrame(np.random.rand(20,4), columns = ['A','B','C','D'])
  dflookup.lookup(xrange(0,10,2), ['B','C','A','B','D'])


Advanced indexing with integer labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Label-based indexing with integer axis labels is a thorny topic. It has been
discussed heavily on mailing lists and among various members of the scientific
Python community. In pandas, our general viewpoint is that labels matter more
than integer locations. Therefore, advanced indexing with ``.ix`` will always
attempt label-based indexing, before falling back on integer-based indexing.

Setting values in mixed-type DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _indexing.mixed_type_setting:

Setting values on a mixed-type DataFrame or Panel is supported when using scalar
values, though setting arbitrary vectors is not yet supported:

.. ipython:: python

   df2 = df[:4]
   df2['foo'] = 'bar'
   print df2
   df2.ix[2] = np.nan
   print df2
   print df2.dtypes

.. _indexing.class:

Index objects
-------------

The pandas Index class and its subclasses can be viewed as implementing an
*ordered set* in addition to providing the support infrastructure necessary for
lookups, data alignment, and reindexing. The easiest way to create one directly
is to pass a list or other sequence to ``Index``:

.. ipython:: python

   index = Index(['e', 'd', 'a', 'b'])
   index
   'd' in index

You can also pass a ``name`` to be stored in the index:


.. ipython:: python

   index = Index(['e', 'd', 'a', 'b'], name='something')
   index.name

Starting with pandas 0.5, the name, if set, will be shown in the console
display:

.. ipython:: python

   index = Index(range(5), name='rows')
   columns = Index(['A', 'B', 'C'], name='cols')
   df = DataFrame(np.random.randn(5, 3), index=index, columns=columns)
   df
   df['A']


Set operations on Index objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _indexing.set_ops:

The three main operations are ``union (|)``, ``intersection (&)``, and ``diff
(-)``. These can be directly called as instance methods or used via overloaded
operators:

.. ipython:: python

   a = Index(['c', 'b', 'a'])
   b = Index(['c', 'e', 'd'])
   a.union(b)
   a | b
   a & b
   a - b

``isin`` method of Index objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One additional operation is the ``isin`` method that works analogously to the
``Series.isin`` method found :ref:`here <indexing.boolean>`.

.. _indexing.hierarchical:

Hierarchical indexing (MultiIndex)
----------------------------------

Hierarchical indexing (also referred to as "multi-level" indexing) is brand new
in the pandas 0.4 release. It is very exciting as it opens the door to some
quite sophisticated data analysis and manipulation, especially for working with
higher dimensional data. In essence, it enables you to effectively store and
manipulate arbitrarily high dimension data in a 2-dimensional tabular structure
(DataFrame), for example. It is not limited to DataFrame

In this section, we will show what exactly we mean by "hierarchical" indexing
and how it integrates with the all of the pandas indexing functionality
described above and in prior sections. Later, when discussing :ref:`group by
<groupby>` and :ref:`pivoting and reshaping data <reshaping>`, we'll show
non-trivial applications to illustrate how it aids in structuring data for
analysis.

.. note::

   Given that hierarchical indexing is so new to the library, it is definitely
   "bleeding-edge" functionality but is certainly suitable for production. But,
   there may inevitably be some minor API changes as more use cases are explored
   and any weaknesses in the design / implementation are identified. pandas aims
   to be "eminently usable" so any feedback about new functionality like this is
   extremely helpful.

Creating a MultiIndex (hierarchical index) object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``MultiIndex`` object is the hierarchical analogue of the standard
``Index`` object which typically stores the axis labels in pandas objects. You
can think of ``MultiIndex`` an array of tuples where each tuple is unique. A
``MultiIndex`` can be created from a list of arrays (using
``MultiIndex.from_arrays``) or an array of tuples (using
``MultiIndex.from_tuples``).

.. ipython:: python

   arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
             ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
   tuples = zip(*arrays)
   tuples
   index = MultiIndex.from_tuples(tuples, names=['first', 'second'])
   s = Series(randn(8), index=index)
   s

All of the ``MultiIndex`` constructors accept a ``names`` argument which stores
string names for the levels themselves. If no names are provided, some
arbitrary ones will be assigned:

.. ipython:: python

   index.names

This index can back any axis of a pandas object, and the number of **levels**
of the index is up to you:

.. ipython:: python

   df = DataFrame(randn(3, 8), index=['A', 'B', 'C'], columns=index)
   df
   DataFrame(randn(6, 6), index=index[:6], columns=index[:6])

We've "sparsified" the higher levels of the indexes to make the console output a
bit easier on the eyes.

It's worth keeping in mind that there's nothing preventing you from using tuples
as atomic labels on an axis:

.. ipython:: python

   Series(randn(8), index=tuples)

The reason that the ``MultiIndex`` matters is that it can allow you to do
grouping, selection, and reshaping operations as we will describe below and in
subsequent areas of the documentation. As you will see in later sections, you
can find yourself working with hierarchically-indexed data without creating a
``MultiIndex`` explicitly yourself. However, when loading data from a file, you
may wish to generate your own ``MultiIndex`` when preparing the data set.

Reconstructing the level labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _indexing.get_level_values:

The method ``get_level_values`` will return a vector of the labels for each
location at a particular level:

.. ipython:: python

   index.get_level_values(0)
   index.get_level_values('second')


Basic indexing on axis with MultiIndex
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One of the important features of hierarchical indexing is that you can select
data by a "partial" label identifying a subgroup in the data. **Partial**
selection "drops" levels of the hierarchical index in the result in a completely
analogous way to selecting a column in a regular DataFrame:

.. ipython:: python

   df['bar']
   df['bar', 'one']
   df['bar']['one']
   s['qux']

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

Advanced indexing with hierarchical index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Syntactically integrating ``MultiIndex`` in advanced indexing with ``.ix`` is a
bit challenging, but we've made every effort to do so. for example the
following works as you would expect:

.. ipython:: python

   df = df.T
   df
   df.ix['bar']
   df.ix['bar', 'two']

"Partial" slicing also works quite nicely:

.. ipython:: python

   df.ix['baz':'foo']
   df.ix[('baz', 'two'):('qux', 'one')]
   df.ix[('baz', 'two'):'foo']

Passing a list of labels or tuples works similar to reindexing:

.. ipython:: python

   df.ix[[('bar', 'two'), ('qux', 'one')]]

The following does not work, and it's not clear if it should or not:

::

   >>> df.ix[['bar', 'qux']]

The code for implementing ``.ix`` makes every attempt to "do the right thing"
but as you use it you may uncover corner cases or unintuitive behavior. If you
do find something like this, do not hesitate to report the issue or ask on the
mailing list.

.. _indexing.xs:

Cross-section with hierarchical index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``xs`` method of ``DataFrame`` additionally takes a level argument to make
selecting data at a particular level of a MultiIndex easier.

.. ipython:: python

    df.xs('one', level='second')

.. _indexing.advanced_reindex:

Advanced reindexing and alignment with hierarchical index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The parameter ``level`` has been added to the ``reindex`` and ``align`` methods
of pandas objects. This is useful to broadcast values across a level. For
instance:

.. ipython:: python

   midx = MultiIndex(levels=[['zero', 'one'], ['x','y']],
                     labels=[[1,1,0,0],[1,0,1,0]])
   df = DataFrame(randn(4,2), index=midx)
   print df
   df2 = df.mean(level=0)
   print df2
   print df2.reindex(df.index, level=0)
   df_aligned, df2_aligned = df.align(df2, level=0)
   print df_aligned
   print df2_aligned


The need for sortedness
~~~~~~~~~~~~~~~~~~~~~~~

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

.. _indexing.sortlevel_byname:

Note, you may also pass a level name to ``sortlevel`` if the MultiIndex levels
are named.

.. ipython:: python

   s.index.names = ['L1', 'L2']
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
   Exception: MultiIndex lexsort depth 1, key was length 2

Swapping levels with ``swaplevel``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``swaplevel`` function can switch the order of two levels:

.. ipython:: python

   df[:5]
   df[:5].swaplevel(0, 1, axis=0)

.. _indexing.reorderlevels:

Reordering levels with ``reorder_levels``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``reorder_levels`` function generalizes the ``swaplevel`` function,
allowing you to permute the hierarchical index levels in one step:

.. ipython:: python

   df[:5].reorder_levels([1,0], axis=0)


Some gory internal details
~~~~~~~~~~~~~~~~~~~~~~~~~~

Internally, the ``MultiIndex`` consists of a few things: the **levels**, the
integer **labels**, and the level **names**:

.. ipython:: python

   index
   index.levels
   index.labels
   index.names

You can probably guess that the labels determine which unique element is
identified with that location at each layer of the index. It's important to
note that sortedness is determined **solely** from the integer labels and does
not check (or care) whether the levels themselves are sorted. Fortunately, the
constructors ``from_tuples`` and ``from_arrays`` ensure that this is true, but
if you compute the levels and labels yourself, please be careful.

Adding an index to an existing DataFrame
----------------------------------------

Occasionally you will load or create a data set into a DataFrame and want to
add an index after you've already done so. There are a couple of different
ways.

Add an index using DataFrame columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _indexing.set_index:

DataFrame has a ``set_index`` method which takes a column name (for a regular
``Index``) or a list of column names (for a ``MultiIndex``), to create a new,
indexed DataFrame:

.. ipython:: python
   :suppress:

   data = DataFrame({'a' : ['bar', 'bar', 'foo', 'foo'],
                     'b' : ['one', 'two', 'one', 'two'],
                     'c' : ['z', 'y', 'x', 'w'],
                     'd' : [1., 2., 3, 4]})

.. ipython:: python

   data
   indexed1 = data.set_index('c')
   indexed1
   indexed2 = data.set_index(['a', 'b'])
   indexed2

Other options in ``set_index`` allow you not drop the index columns or to add
the index in-place (without creating a new object):

.. ipython:: python

   data.set_index('c', drop=False)
   df = data.set_index(['a', 'b'], inplace=True)
   data

Remove / reset the index,  ``reset_index``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As a convenience, there is a new function on DataFrame called ``reset_index``
which transfers the index values into the DataFrame's columns and sets a simple
integer index. This is the inverse operation to ``set_index``

.. ipython:: python

   df
   df.reset_index()

The output is more similar to a SQL table or a record array. The names for the
columns derived from the index are the ones stored in the ``names`` attribute.

``reset_index`` takes an optional parameter ``drop`` which if true simply
discards the index, instead of putting index values in the DataFrame's columns.

.. note::

   The ``reset_index`` method used to be called ``delevel`` which is now deprecated.

Adding an ad hoc index
~~~~~~~~~~~~~~~~~~~~~~

If you create an index yourself, you can just assign it to the ``index`` field:

.. code-block:: python

   df.index = index

Indexing internal details
-------------------------

.. note::

    The following is largely relevant for those actually working on the pandas
    codebase. And the source code is still the best place to look at the
    specifics of how things are implemented.

In pandas there are a few objects implemented which can serve as valid
containers for the axis labels:

  - ``Index``: the generic "ordered set" object, an ndarray of object dtype
    assuming nothing about its contents. The labels must be hashable (and
    likely immutable) and unique. Populates a dict of label to location in
    Cython to do :math:`O(1)` lookups.
  - ``Int64Index``: a version of ``Index`` highly optimized for 64-bit integer
    data, such as time stamps
  - ``MultiIndex``: the standard hierarchical index object
  - ``DateRange``: fixed frequency date range generated from a time rule or
    DateOffset. An ndarray of Python datetime objects

The motivation for having an ``Index`` class in the first place was to enable
different implementations of indexing. This means that it's possible for you,
the user, to implement a custom ``Index`` subclass that may be better suited to
a particular application than the ones provided in pandas. For example, we plan
to add a more efficient datetime index which leverages the new
``numpy.datetime64`` dtype in the relatively near future.

From an internal implementation point of view, the relevant methods that an
``Index`` must define are one or more of the following (depending on how
incompatible the new object internals are with the ``Index`` functions):

  - ``get_loc``: returns an "indexer" (an integer, or in some cases a
    slice object) for a label
  - ``slice_locs``: returns the "range" to slice between two labels
  - ``get_indexer``: Computes the indexing vector for reindexing / data
    alignment purposes. See the source / docstrings for more on this
  - ``reindex``: Does any pre-conversion of the input index then calls
    ``get_indexer``
  - ``union``, ``intersection``: computes the union or intersection of two
    Index objects
  - ``insert``: Inserts a new label into an Index, yielding a new object
  - ``delete``: Delete a label, yielding a new object
  - ``drop``: Deletes a set of labels
  - ``take``: Analogous to ndarray.take
