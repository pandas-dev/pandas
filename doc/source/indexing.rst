.. _indexing:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(12345)
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

In this section / chapter, we will focus on the latter set of functionality,
namely how to slice, dice, and generally get and set subsets of pandas
objects. The primary focus will be on Series and DataFrame as they have
received more development attention in this area. More work will be invested in
WidePanel and future higher-dimensional data structures in the future,
especially in label-based advanced indexing.

.. _indexing.basics:

Basics
------

As mentioned when introducing the data structures in the :ref:`last section
<basics>`, the primary function of indexing with ``[]``
(a.k.a. ``__getitem__``) for those familiar with implementing class behavior in
Python) is selecting out lower-dimensional slices. Thus,

  - **Series**: ``series[label]`` returns a scalar value
  - **DataFrame**: ``frame[colname]`` returns a Series corresponding to the
    passed column name
  - **WidePanel**: ``panel[itemname]`` returns a DataFrame corresponding to the
    passed item name

Here we construct a simple time series data set to use for illustrating the
indexing functionality:

.. ipython:: python

   dates = np.asarray(DateRange('1/1/2000', periods=8))
   df = DataFrame(randn(8, 4), index=dates, columns=['A', 'B', 'C', 'D'])
   df
   panel = WidePanel({'one' : df, 'two' : df - df.mean()})
   panel

.. note::

   None of the indexing functionality is time series specific unless
   specifically stated.

Thus, as per above, we have the most basic indexing using ``[]``:

.. ipython:: python

   s = df['A']
   s[dates[5]]
   panel['two']

Data slices on other axes
~~~~~~~~~~~~~~~~~~~~~~~~~

It's certainly possible to retrieve data slices along the other axes of a
DataFrame or WidePanel. We tend to refer to these slices as
*cross-sections*. DataFrame has the ``xs`` function for retrieving rows as
Series and WidePanel has the analogous ``major_xs`` and ``minor_xs`` functions
for retrieving slices as DataFrames for a given ``major_axis`` or
``minor_axis`` label, respectively.

.. ipython:: python

   date = dates[5]
   df.xs(date)
   panel.major_xs(date)
   panel.minor_xs('A')

.. note::

   See :ref:`advanced indexing <indexing.advanced>` below for an alternate and
   more concise way of doing the same thing.

Slicing ranges
~~~~~~~~~~~~~~

:ref:`Advanced indexing <indexing.advanced>` detailed below is the most robust
and consistent way of slicing integer ranges, e.g. ``obj[5:10]``, across all of
the data structures and their axes. On Series, this syntax works exactly as
expected as with an ndarray, returning a slice of the values and the
corresponding labels:

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

Using a boolean vector to index a Series works exactly like an ndarray:

.. ipython:: python

   s[s > 0]
   s[(s < 0) & (s > -0.5)]

Again as a convenience, selecting rows from a DataFrame using a boolean vector
the same length as the DataFrame's index (for example, something derived from
one of the columns of the DataFrame) is supported:

.. ipython:: python

   df[df['A'] > 0]

With the advanced indexing capabilities discussed later, you are able to do
boolean indexing in any of axes or combine a boolean vector with an indexing
expression on one of the other axes

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

Slicing ranges
~~~~~~~~~~~~~~

Similar to Python lists and ndarrays, for convenience DataFrame
supports slicing:

.. ipython:: python

    df[:2]
    df[::-1]
    df[-3:].T

Boolean indexing
~~~~~~~~~~~~~~~~

As another indexing convenience, it is possible to use boolean
indexing to select rows of a DataFrame:

.. ipython:: python

    df[df['A'] > 0.5]

As we will see later on, the same operation could be accomplished by
reindexing. However, the syntax would be more verbose; hence, the inclusion of
this indexing method.

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

We'll illustrate all of these methods on the above time series data set.

Advanced indexing with integer labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _indexing.hierarchical:

Hierarchical indexing (MultiIndex)
----------------------------------

Advanced indexing with hierarchical index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Indexing internal details
-------------------------

.. note::

    The following is largely relevant for those actually working on the pandas
    codebase. And the source code is still the best place to look at the
    specifics of how things are implemented.
