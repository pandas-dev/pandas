.. _10min:

{{ header }}

********************
10 Minutes to pandas
********************

This is a short introduction to pandas, geared mainly for new users.
You can see more complex recipes in the :ref:`Cookbook<cookbook>`.

Customarily, we import as follows:

.. ipython:: python

   import numpy as np
   import pandas as pd

Object Creation
---------------

See the :ref:`Data Structure Intro section <dsintro>`.

Creating a :class:`Series` by passing a list of values, letting pandas create
a default integer index:

.. ipython:: python

   s = pd.Series([1, 3, 5, np.nan, 6, 8])
   s

Creating a :class:`DataFrame` by passing a NumPy array, with a datetime index
and labeled columns:

.. ipython:: python

   dates = pd.date_range('20130101', periods=6)
   dates
   df = pd.DataFrame(np.random.randn(6, 4), index=dates, columns=list('ABCD'))
   df

Creating a ``DataFrame`` by passing a dict of objects that can be converted to series-like.

.. ipython:: python

   df2 = pd.DataFrame({'A': 1.,
                       'B': pd.Timestamp('20130102'),
                       'C': pd.Series(1, index=list(range(4)), dtype='float32'),
                       'D': np.array([3] * 4, dtype='int32'),
                       'E': pd.Categorical(["test", "train", "test", "train"]),
                       'F': 'foo'})
   df2

The columns of the resulting ``DataFrame`` have different
:ref:`dtypes <basics.dtypes>`.

.. ipython:: python

   df2.dtypes

If you're using IPython, tab completion for column names (as well as public
attributes) is automatically enabled. Here's a subset of the attributes that
will be completed:

.. ipython::

   @verbatim
   In [1]: df2.<TAB>  # noqa: E225, E999
   df2.A                  df2.bool
   df2.abs                df2.boxplot
   df2.add                df2.C
   df2.add_prefix         df2.clip
   df2.add_suffix         df2.clip_lower
   df2.align              df2.clip_upper
   df2.all                df2.columns
   df2.any                df2.combine
   df2.append             df2.combine_first
   df2.apply              df2.compound
   df2.applymap           df2.consolidate
   df2.D

As you can see, the columns ``A``, ``B``, ``C``, and ``D`` are automatically
tab completed. ``E`` is there as well; the rest of the attributes have been
truncated for brevity.

Viewing Data
------------

See the :ref:`Basics section <basics>`.

Here is how to view the top and bottom rows of the frame:

.. ipython:: python

   df.head()
   df.tail(3)

Display the index, columns:

.. ipython:: python

   df.index
   df.columns

:meth:`DataFrame.to_numpy` gives a NumPy representation of the underlying data.
Note that this can be an expensive operation when your :class:`DataFrame` has
columns with different data types, which comes down to a fundamental difference
between pandas and NumPy: **NumPy arrays have one dtype for the entire array,
while pandas DataFrames have one dtype per column**. When you call
:meth:`DataFrame.to_numpy`, pandas will find the NumPy dtype that can hold *all*
of the dtypes in the DataFrame. This may end up being ``object``, which requires
casting every value to a Python object.

For ``df``, our :class:`DataFrame` of all floating-point values,
:meth:`DataFrame.to_numpy` is fast and doesn't require copying data.

.. ipython:: python

   df.to_numpy()

For ``df2``, the :class:`DataFrame` with multiple dtypes,
:meth:`DataFrame.to_numpy` is relatively expensive.

.. ipython:: python

   df2.to_numpy()

.. note::

   :meth:`DataFrame.to_numpy` does *not* include the index or column
   labels in the output.

:func:`~DataFrame.describe` shows a quick statistic summary of your data:

.. ipython:: python

   df.describe()

Transposing your data:

.. ipython:: python

   df.T

Sorting by an axis:

.. ipython:: python

   df.sort_index(axis=1, ascending=False)

Sorting by values:

.. ipython:: python

   df.sort_values(by='B')

Selection
---------

.. note::

   While standard Python / Numpy expressions for selecting and setting are
   intuitive and come in handy for interactive work, for production code, we
   recommend the optimized pandas data access methods, ``.at``, ``.iat``,
   ``.loc`` and ``.iloc``.

See the indexing documentation :ref:`Indexing and Selecting Data <indexing>` and :ref:`MultiIndex / Advanced Indexing <advanced>`.

Getting
~~~~~~~

Selecting a single column, which yields a ``Series``,
equivalent to ``df.A``:

.. ipython:: python

   df['A']

Selecting via ``[]``, which slices the rows.

.. ipython:: python

   df[0:3]
   df['20130102':'20130104']

Selection by Label
~~~~~~~~~~~~~~~~~~

See more in :ref:`Selection by Label <indexing.label>`.

For getting a cross section using a label:

.. ipython:: python

   df.loc[dates[0]]

Selecting on a multi-axis by label:

.. ipython:: python

   df.loc[:, ['A', 'B']]

Showing label slicing, both endpoints are *included*:

.. ipython:: python

   df.loc['20130102':'20130104', ['A', 'B']]

Reduction in the dimensions of the returned object:

.. ipython:: python

   df.loc['20130102', ['A', 'B']]

For getting a scalar value:

.. ipython:: python

   df.loc[dates[0], 'A']

For getting fast access to a scalar (equivalent to the prior method):

.. ipython:: python

   df.at[dates[0], 'A']

Selection by Position
~~~~~~~~~~~~~~~~~~~~~

See more in :ref:`Selection by Position <indexing.integer>`.

Select via the position of the passed integers:

.. ipython:: python

   df.iloc[3]

By integer slices, acting similar to numpy/python:

.. ipython:: python

   df.iloc[3:5, 0:2]

By lists of integer position locations, similar to the numpy/python style:

.. ipython:: python

   df.iloc[[1, 2, 4], [0, 2]]

For slicing rows explicitly:

.. ipython:: python

   df.iloc[1:3, :]

For slicing columns explicitly:

.. ipython:: python

   df.iloc[:, 1:3]

For getting a value explicitly:

.. ipython:: python

   df.iloc[1, 1]

For getting fast access to a scalar (equivalent to the prior method):

.. ipython:: python

   df.iat[1, 1]

Boolean Indexing
~~~~~~~~~~~~~~~~

Using a single column's values to select data.

.. ipython:: python

   df[df.A > 0]

Selecting values from a DataFrame where a boolean condition is met.

.. ipython:: python

   df[df > 0]

Using the :func:`~Series.isin` method for filtering:

.. ipython:: python

   df2 = df.copy()
   df2['E'] = ['one', 'one', 'two', 'three', 'four', 'three']
   df2
   df2[df2['E'].isin(['two', 'four'])]

Setting
~~~~~~~

Setting a new column automatically aligns the data
by the indexes.

.. ipython:: python

   s1 = pd.Series([1, 2, 3, 4, 5, 6], index=pd.date_range('20130102', periods=6))
   s1
   df['F'] = s1

Setting values by label:

.. ipython:: python

   df.at[dates[0], 'A'] = 0

Setting values by position:

.. ipython:: python

   df.iat[0, 1] = 0

Setting by assigning with a NumPy array:

.. ipython:: python

   df.loc[:, 'D'] = np.array([5] * len(df))

The result of the prior setting operations.

.. ipython:: python

   df

A ``where`` operation with setting.

.. ipython:: python

   df2 = df.copy()
   df2[df2 > 0] = -df2
   df2


Missing Data
------------

pandas primarily uses the value ``np.nan`` to represent missing data. It is by
default not included in computations. See the :ref:`Missing Data section
<missing_data>`.

Reindexing allows you to change/add/delete the index on a specified axis. This
returns a copy of the data.

.. ipython:: python

   df1 = df.reindex(index=dates[0:4], columns=list(df.columns) + ['E'])
   df1.loc[dates[0]:dates[1], 'E'] = 1
   df1

To drop any rows that have missing data.

.. ipython:: python

   df1.dropna(how='any')

Filling missing data.

.. ipython:: python

   df1.fillna(value=5)

To get the boolean mask where values are ``nan``.

.. ipython:: python

   pd.isna(df1)


Operations
----------

See the :ref:`Basic section on Binary Ops <basics.binop>`.

Stats
~~~~~

Operations in general *exclude* missing data.

Performing a descriptive statistic:

.. ipython:: python

   df.mean()

Same operation on the other axis:

.. ipython:: python

   df.mean(1)

Operating with objects that have different dimensionality and need alignment.
In addition, pandas automatically broadcasts along the specified dimension.

.. ipython:: python

   s = pd.Series([1, 3, 5, np.nan, 6, 8], index=dates).shift(2)
   s
   df.sub(s, axis='index')


Apply
~~~~~

Applying functions to the data:

.. ipython:: python

   df.apply(np.cumsum)
   df.apply(lambda x: x.max() - x.min())

Histogramming
~~~~~~~~~~~~~

See more at :ref:`Histogramming and Discretization <basics.discretization>`.

.. ipython:: python

   s = pd.Series(np.random.randint(0, 7, size=10))
   s
   s.value_counts()

String Methods
~~~~~~~~~~~~~~

Series is equipped with a set of string processing methods in the `str`
attribute that make it easy to operate on each element of the array, as in the
code snippet below. Note that pattern-matching in `str` generally uses `regular
expressions <https://docs.python.org/3/library/re.html>`__ by default (and in
some cases always uses them). See more at :ref:`Vectorized String Methods
<text.string_methods>`.

.. ipython:: python

   s = pd.Series(['A', 'B', 'C', 'Aaba', 'Baca', np.nan, 'CABA', 'dog', 'cat'])
   s.str.lower()

Merge
-----

Concat
~~~~~~

pandas provides various facilities for easily combining together Series and
DataFrame objects with various kinds of set logic for the indexes
and relational algebra functionality in the case of join / merge-type
operations.

See the :ref:`Merging section <merging>`.

Concatenating pandas objects together with :func:`concat`:

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 4))
   df

   # break it into pieces
   pieces = [df[:3], df[3:7], df[7:]]

   pd.concat(pieces)

Join
~~~~

SQL style merges. See the :ref:`Database style joining <merging.join>` section.

.. ipython:: python

   left = pd.DataFrame({'key': ['foo', 'foo'], 'lval': [1, 2]})
   right = pd.DataFrame({'key': ['foo', 'foo'], 'rval': [4, 5]})
   left
   right
   pd.merge(left, right, on='key')

Another example that can be given is:

.. ipython:: python

   left = pd.DataFrame({'key': ['foo', 'bar'], 'lval': [1, 2]})
   right = pd.DataFrame({'key': ['foo', 'bar'], 'rval': [4, 5]})
   left
   right
   pd.merge(left, right, on='key')


Append
~~~~~~

Append rows to a dataframe. See the :ref:`Appending <merging.concatenation>`
section.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(8, 4), columns=['A', 'B', 'C', 'D'])
   df
   s = df.iloc[3]
   df.append(s, ignore_index=True)


Grouping
--------

By "group by" we are referring to a process involving one or more of the
following steps:

 - **Splitting** the data into groups based on some criteria
 - **Applying** a function to each group independently
 - **Combining** the results into a data structure

See the :ref:`Grouping section <groupby>`.

.. ipython:: python

   df = pd.DataFrame({'A': ['foo', 'bar', 'foo', 'bar',
                            'foo', 'bar', 'foo', 'foo'],
                      'B': ['one', 'one', 'two', 'three',
                            'two', 'two', 'one', 'three'],
                      'C': np.random.randn(8),
                      'D': np.random.randn(8)})
   df

Grouping and then applying the :meth:`~DataFrame.sum` function to the resulting
groups.

.. ipython:: python

   df.groupby('A').sum()

Grouping by multiple columns forms a hierarchical index, and again we can
apply the ``sum`` function.

.. ipython:: python

   df.groupby(['A', 'B']).sum()

Reshaping
---------

See the sections on :ref:`Hierarchical Indexing <advanced.hierarchical>` and
:ref:`Reshaping <reshaping.stacking>`.

Stack
~~~~~

.. ipython:: python

   tuples = list(zip(*[['bar', 'bar', 'baz', 'baz',
                        'foo', 'foo', 'qux', 'qux'],
                       ['one', 'two', 'one', 'two',
                        'one', 'two', 'one', 'two']]))
   index = pd.MultiIndex.from_tuples(tuples, names=['first', 'second'])
   df = pd.DataFrame(np.random.randn(8, 2), index=index, columns=['A', 'B'])
   df2 = df[:4]
   df2

The :meth:`~DataFrame.stack` method "compresses" a level in the DataFrame's
columns.

.. ipython:: python

   stacked = df2.stack()
   stacked

With a "stacked" DataFrame or Series (having a ``MultiIndex`` as the
``index``), the inverse operation of :meth:`~DataFrame.stack` is
:meth:`~DataFrame.unstack`, which by default unstacks the **last level**:

.. ipython:: python

   stacked.unstack()
   stacked.unstack(1)
   stacked.unstack(0)

Pivot Tables
~~~~~~~~~~~~
See the section on :ref:`Pivot Tables <reshaping.pivot>`.

.. ipython:: python

   df = pd.DataFrame({'A': ['one', 'one', 'two', 'three'] * 3,
                      'B': ['A', 'B', 'C'] * 4,
                      'C': ['foo', 'foo', 'foo', 'bar', 'bar', 'bar'] * 2,
                      'D': np.random.randn(12),
                      'E': np.random.randn(12)})
   df

We can produce pivot tables from this data very easily:

.. ipython:: python

   pd.pivot_table(df, values='D', index=['A', 'B'], columns=['C'])


Time Series
-----------

pandas has simple, powerful, and efficient functionality for performing
resampling operations during frequency conversion (e.g., converting secondly
data into 5-minutely data). This is extremely common in, but not limited to,
financial applications. See the :ref:`Time Series section <timeseries>`.

.. ipython:: python

   rng = pd.date_range('1/1/2012', periods=100, freq='S')
   ts = pd.Series(np.random.randint(0, 500, len(rng)), index=rng)
   ts.resample('5Min').sum()

Time zone representation:

.. ipython:: python

   rng = pd.date_range('3/6/2012 00:00', periods=5, freq='D')
   ts = pd.Series(np.random.randn(len(rng)), rng)
   ts
   ts_utc = ts.tz_localize('UTC')
   ts_utc

Converting to another time zone:

.. ipython:: python

   ts_utc.tz_convert('US/Eastern')

Converting between time span representations:

.. ipython:: python

   rng = pd.date_range('1/1/2012', periods=5, freq='M')
   ts = pd.Series(np.random.randn(len(rng)), index=rng)
   ts
   ps = ts.to_period()
   ps
   ps.to_timestamp()

Converting between period and timestamp enables some convenient arithmetic
functions to be used. In the following example, we convert a quarterly
frequency with year ending in November to 9am of the end of the month following
the quarter end:

.. ipython:: python

   prng = pd.period_range('1990Q1', '2000Q4', freq='Q-NOV')
   ts = pd.Series(np.random.randn(len(prng)), prng)
   ts.index = (prng.asfreq('M', 'e') + 1).asfreq('H', 's') + 9
   ts.head()

Categoricals
------------

pandas can include categorical data in a ``DataFrame``. For full docs, see the
:ref:`categorical introduction <categorical>` and the :ref:`API documentation <api.arrays.categorical>`.

.. ipython:: python

    df = pd.DataFrame({"id": [1, 2, 3, 4, 5, 6],
                       "raw_grade": ['a', 'b', 'b', 'a', 'a', 'e']})

Convert the raw grades to a categorical data type.

.. ipython:: python

    df["grade"] = df["raw_grade"].astype("category")
    df["grade"]

Rename the categories to more meaningful names (assigning to
``Series.cat.categories`` is inplace!).

.. ipython:: python

    df["grade"].cat.categories = ["very good", "good", "very bad"]

Reorder the categories and simultaneously add the missing categories (methods under ``Series
.cat`` return a new ``Series`` by default).

.. ipython:: python

    df["grade"] = df["grade"].cat.set_categories(["very bad", "bad", "medium",
                                                  "good", "very good"])
    df["grade"]

Sorting is per order in the categories, not lexical order.

.. ipython:: python

    df.sort_values(by="grade")

Grouping by a categorical column also shows empty categories.

.. ipython:: python

    df.groupby("grade").size()


Plotting
--------

See the :ref:`Plotting <visualization>` docs.

.. ipython:: python
   :suppress:

   import matplotlib.pyplot as plt
   plt.close('all')

.. ipython:: python

   ts = pd.Series(np.random.randn(1000),
                  index=pd.date_range('1/1/2000', periods=1000))
   ts = ts.cumsum()

   @savefig series_plot_basic.png
   ts.plot()

On a DataFrame, the :meth:`~DataFrame.plot` method is a convenience to plot all
of the columns with labels:

.. ipython:: python

   df = pd.DataFrame(np.random.randn(1000, 4), index=ts.index,
                     columns=['A', 'B', 'C', 'D'])
   df = df.cumsum()

   plt.figure()
   df.plot()
   @savefig frame_plot_basic.png
   plt.legend(loc='best')

Getting Data In/Out
-------------------

CSV
~~~

:ref:`Writing to a csv file. <io.store_in_csv>`

.. ipython:: python

   df.to_csv('foo.csv')

:ref:`Reading from a csv file. <io.read_csv_table>`

.. ipython:: python

   pd.read_csv('foo.csv')

.. ipython:: python
   :suppress:

   import os
   os.remove('foo.csv')

HDF5
~~~~

Reading and writing to :ref:`HDFStores <io.hdf5>`.

Writing to a HDF5 Store.

.. ipython:: python

   df.to_hdf('foo.h5', 'df')

Reading from a HDF5 Store.

.. ipython:: python

   pd.read_hdf('foo.h5', 'df')

.. ipython:: python
   :suppress:

   os.remove('foo.h5')

Excel
~~~~~

Reading and writing to :ref:`MS Excel <io.excel>`.

Writing to an excel file.

.. ipython:: python

   df.to_excel('foo.xlsx', sheet_name='Sheet1')

Reading from an excel file.

.. ipython:: python

   pd.read_excel('foo.xlsx', 'Sheet1', index_col=None, na_values=['NA'])

.. ipython:: python
   :suppress:

   os.remove('foo.xlsx')

Gotchas
-------

If you are attempting to perform an operation you might see an exception like:

.. code-block:: python

    >>> if pd.Series([False, True, False]):
    ...     print("I was true")
    Traceback
        ...
    ValueError: The truth value of an array is ambiguous. Use a.empty, a.any() or a.all().

See :ref:`Comparisons<basics.compare>` for an explanation and what to do.

See :ref:`Gotchas<gotchas>` as well.
