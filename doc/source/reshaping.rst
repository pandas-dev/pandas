.. currentmodule:: pandas
.. _reshaping:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   options.display.max_rows=15
   from pandas.core.reshape import *
   import pandas.util.testing as tm
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   from pandas.tools.tile import *
   from pandas.compat import zip

**************************
Reshaping and Pivot Tables
**************************

Reshaping by pivoting DataFrame objects
---------------------------------------

.. ipython::
   :suppress:

   In [1]: import pandas.util.testing as tm; tm.N = 3

   In [2]: def unpivot(frame):
      ...:         N, K = frame.shape
      ...:         data = {'value' : frame.values.ravel('F'),
      ...:                 'variable' : np.asarray(frame.columns).repeat(N),
      ...:                 'date' : np.tile(np.asarray(frame.index), K)}
      ...:         columns = ['date', 'variable', 'value']
      ...:         return DataFrame(data, columns=columns)
      ...:

   In [3]: df = unpivot(tm.makeTimeDataFrame())

Data is often stored in CSV files or databases in so-called "stacked" or
"record" format:

.. ipython:: python

   df


For the curious here is how the above DataFrame was created:

.. code-block:: python

   import pandas.util.testing as tm; tm.N = 3
   def unpivot(frame):
       N, K = frame.shape
       data = {'value' : frame.values.ravel('F'),
               'variable' : np.asarray(frame.columns).repeat(N),
               'date' : np.tile(np.asarray(frame.index), K)}
       return DataFrame(data, columns=['date', 'variable', 'value'])
   df = unpivot(tm.makeTimeDataFrame())

To select out everything for variable ``A`` we could do:

.. ipython:: python

   df[df['variable'] == 'A']

But suppose we wish to do time series operations with the variables. A better
representation would be where the ``columns`` are the unique variables and an
``index`` of dates identifies individual observations. To reshape the data into
this form, use the ``pivot`` function:

.. ipython:: python

   df.pivot(index='date', columns='variable', values='value')

If the ``values`` argument is omitted, and the input DataFrame has more than
one column of values which are not used as column or index inputs to ``pivot``,
then the resulting "pivoted" DataFrame will have :ref:`hierarchical columns
<indexing.hierarchical>` whose topmost level indicates the respective value
column:

.. ipython:: python

   df['value2'] = df['value'] * 2
   pivoted = df.pivot('date', 'variable')
   pivoted

You of course can then select subsets from the pivoted DataFrame:

.. ipython:: python

   pivoted['value2']

Note that this returns a view on the underlying data in the case where the data
are homogeneously-typed.

.. _reshaping.stacking:

Reshaping by stacking and unstacking
------------------------------------

Closely related to the ``pivot`` function are the related ``stack`` and
``unstack`` functions currently available on Series and DataFrame. These
functions are designed to work together with ``MultiIndex`` objects (see the
section on :ref:`hierarchical indexing <indexing.hierarchical>`). Here are
essentially what these functions do:

  - ``stack``: "pivot" a level of the (possibly hierarchical) column labels,
    returning a DataFrame with an index with a new inner-most level of row
    labels.
  - ``unstack``: inverse operation from ``stack``: "pivot" a level of the
    (possibly hierarchical) row index to the column axis, producing a reshaped
    DataFrame with a new inner-most level of column labels.

The clearest way to explain is by example. Let's take a prior example data set
from the hierarchical indexing section:

.. ipython:: python

   tuples = list(zip(*[['bar', 'bar', 'baz', 'baz',
                   'foo', 'foo', 'qux', 'qux'],
                  ['one', 'two', 'one', 'two',
                   'one', 'two', 'one', 'two']]))
   index = MultiIndex.from_tuples(tuples, names=['first', 'second'])
   df = DataFrame(randn(8, 2), index=index, columns=['A', 'B'])
   df2 = df[:4]
   df2

The ``stack`` function "compresses" a level in the DataFrame's columns to
produce either:

  - A Series, in the case of a simple column Index
  - A DataFrame, in the case of a ``MultiIndex`` in the columns

If the columns have a ``MultiIndex``, you can choose which level to stack. The
stacked level becomes the new lowest level in a ``MultiIndex`` on the columns:

.. ipython:: python

   stacked = df2.stack()
   stacked

With a "stacked" DataFrame or Series (having a ``MultiIndex`` as the
``index``), the inverse operation of ``stack`` is ``unstack``, which by default
unstacks the **last level**:

.. ipython:: python

   stacked.unstack()
   stacked.unstack(1)
   stacked.unstack(0)

.. _reshaping.unstack_by_name:

If the indexes have names, you can use the level names instead of specifying
the level numbers:

.. ipython:: python

   stacked.unstack('second')

You may also stack or unstack more than one level at a time by passing a list
of levels, in which case the end result is as if each level in the list were
processed individually.

These functions are intelligent about handling missing data and do not expect
each subgroup within the hierarchical index to have the same set of labels.
They also can handle the index being unsorted (but you can make it sorted by
calling ``sortlevel``, of course). Here is a more complex example:

.. ipython:: python

   columns = MultiIndex.from_tuples([('A', 'cat'), ('B', 'dog'),
                                     ('B', 'cat'), ('A', 'dog')],
                                    names=['exp', 'animal'])
   df = DataFrame(randn(8, 4), index=index, columns=columns)
   df2 = df.ix[[0, 1, 2, 4, 5, 7]]
   df2

As mentioned above, ``stack`` can be called with a ``level`` argument to select
which level in the columns to stack:

.. ipython:: python

   df2.stack('exp')
   df2.stack('animal')

Unstacking when the columns are a ``MultiIndex`` is also careful about doing
the right thing:

.. ipython:: python

   df[:3].unstack(0)
   df2.unstack(1)

.. _reshaping.melt:

Reshaping by Melt
-----------------

The :func:`~pandas.melt` function is useful to massage a
DataFrame into a format where one or more columns are identifier variables,
while all other columns, considered measured variables, are "unpivoted" to the
row axis, leaving just two non-identifier columns, "variable" and "value". The
names of those columns can be customized by supplying the ``var_name`` and
``value_name`` parameters.

For instance,

.. ipython:: python

   cheese = DataFrame({'first' : ['John', 'Mary'],
                       'last' : ['Doe', 'Bo'],
                       'height' : [5.5, 6.0],
                       'weight' : [130, 150]})
   cheese
   melt(cheese, id_vars=['first', 'last'])
   melt(cheese, id_vars=['first', 'last'], var_name='quantity')

Another way to transform is to use the ``wide_to_long`` panel data convenience function.

.. ipython:: python

  dft = pd.DataFrame({"A1970" : {0 : "a", 1 : "b", 2 : "c"},
                      "A1980" : {0 : "d", 1 : "e", 2 : "f"},
                      "B1970" : {0 : 2.5, 1 : 1.2, 2 : .7},
                      "B1980" : {0 : 3.2, 1 : 1.3, 2 : .1},
                      "X"     : dict(zip(range(3), np.random.randn(3)))
                     })
  dft["id"] = dft.index
  dft
  pd.wide_to_long(dft, ["A", "B"], i="id", j="year")

Combining with stats and GroupBy
--------------------------------

It should be no shock that combining ``pivot`` / ``stack`` / ``unstack`` with
GroupBy and the basic Series and DataFrame statistical functions can produce
some very expressive and fast data manipulations.

.. ipython:: python

   df
   df.stack().mean(1).unstack()

   # same result, another way
   df.groupby(level=1, axis=1).mean()

   df.stack().groupby(level=1).mean()

   df.mean().unstack(0)


Pivot tables and cross-tabulations
----------------------------------

.. _reshaping.pivot:

The function ``pandas.pivot_table`` can be used to create spreadsheet-style pivot
tables. See the :ref:`cookbook<cookbook.pivot>` for some advanced strategies

It takes a number of arguments

- ``data``: A DataFrame object
- ``values``: a column or a list of columns to aggregate
- ``index``: a column, Grouper, array which has the same length as data, or list of them.
  Keys to group by on the pivot table index. If an array is passed, it is being used as the same manner as column values.
- ``columns``: a column, Grouper, array which has the same length as data, or list of them.
  Keys to group by on the pivot table column. If an array is passed, it is being used as the same manner as column values.
- ``aggfunc``: function to use for aggregation, defaulting to ``numpy.mean``

Consider a data set like this:

.. ipython:: python

   import datetime
   df = DataFrame({'A' : ['one', 'one', 'two', 'three'] * 6,
                   'B' : ['A', 'B', 'C'] * 8,
                   'C' : ['foo', 'foo', 'foo', 'bar', 'bar', 'bar'] * 4,
                   'D' : np.random.randn(24),
                   'E' : np.random.randn(24),
                   'F' : [datetime.datetime(2013, i, 1) for i in range(1, 13)] +
                         [datetime.datetime(2013, i, 15) for i in range(1, 13)]})
   df

We can produce pivot tables from this data very easily:

.. ipython:: python

   pivot_table(df, values='D', index=['A', 'B'], columns=['C'])
   pivot_table(df, values='D', index=['B'], columns=['A', 'C'], aggfunc=np.sum)
   pivot_table(df, values=['D','E'], index=['B'], columns=['A', 'C'], aggfunc=np.sum)

The result object is a DataFrame having potentially hierarchical indexes on the
rows and columns. If the ``values`` column name is not given, the pivot table
will include all of the data that can be aggregated in an additional level of
hierarchy in the columns:

.. ipython:: python

   pivot_table(df, index=['A', 'B'], columns=['C'])

Also, you can use ``Grouper`` for ``index`` and ``columns`` keywords. For detail of ``Grouper``, see :ref:`Grouping with a Grouper specification <groupby.specify>`.

.. ipython:: python

   pivot_table(df, values='D', index=Grouper(freq='M', key='F'), columns='C')

You can render a nice output of the table omitting the missing values by
calling ``to_string`` if you wish:

.. ipython:: python

   table = pivot_table(df, index=['A', 'B'], columns=['C'])
   print(table.to_string(na_rep=''))

Note that ``pivot_table`` is also available as an instance method on DataFrame.

Cross tabulations
~~~~~~~~~~~~~~~~~

Use the ``crosstab`` function to compute a cross-tabulation of two (or more)
factors. By default ``crosstab`` computes a frequency table of the factors
unless an array of values and an aggregation function are passed.

It takes a number of arguments

- ``index``: array-like, values to group by in the rows
- ``columns``: array-like, values to group by in the columns
- ``values``: array-like, optional, array of values to aggregate according to
  the factors
- ``aggfunc``: function, optional, If no values array is passed, computes a
  frequency table
- ``rownames``: sequence, default None, must match number of row arrays passed
- ``colnames``: sequence, default None, if passed, must match number of column
  arrays passed
- ``margins``: boolean, default False, Add row/column margins (subtotals)

Any Series passed will have their name attributes used unless row or column
names for the cross-tabulation are specified

For example:

.. ipython:: python

    foo, bar, dull, shiny, one, two = 'foo', 'bar', 'dull', 'shiny', 'one', 'two'
    a = np.array([foo, foo, bar, bar, foo, foo], dtype=object)
    b = np.array([one, one, two, one, two, one], dtype=object)
    c = np.array([dull, dull, shiny, dull, dull, shiny], dtype=object)
    crosstab(a, [b, c], rownames=['a'], colnames=['b', 'c'])

.. _reshaping.pivot.margins:

Adding margins (partial aggregates)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you pass ``margins=True`` to ``pivot_table``, special ``All`` columns and
rows will be added with partial group aggregates across the categories on the
rows and columns:

.. ipython:: python

   df.pivot_table(index=['A', 'B'], columns='C', margins=True, aggfunc=np.std)

.. _reshaping.tile:

Tiling
------

.. _reshaping.tile.cut:

The ``cut`` function computes groupings for the values of the input array and
is often used to transform continuous variables to discrete or categorical
variables:

.. ipython:: python

   ages = np.array([10, 15, 13, 12, 23, 25, 28, 59, 60])


   cut(ages, bins=3)

If the ``bins`` keyword is an integer, then equal-width bins are formed.
Alternatively we can specify custom bin-edges:

.. ipython:: python

   cut(ages, bins=[0, 18, 35, 70])


.. _reshaping.dummies:

Computing indicator / dummy variables
-------------------------------------

To convert a categorical variable into a "dummy" or "indicator" DataFrame, for example
a column in a DataFrame (a Series) which has ``k`` distinct values, can derive a DataFrame
containing ``k`` columns of 1s and 0s:

.. ipython:: python

   df = DataFrame({'key': list('bbacab'), 'data1': range(6)})


   get_dummies(df['key'])

Sometimes it's useful to prefix the column names, for example when merging the result
with the original DataFrame:

.. ipython:: python

   dummies = get_dummies(df['key'], prefix='key')
   dummies


   df[['data1']].join(dummies)

This function is often used along with discretization functions like ``cut``:

.. ipython:: python

   values = randn(10)
   values


   bins = [0, 0.2, 0.4, 0.6, 0.8, 1]


   get_dummies(cut(values, bins))

See also :func:`Series.str.get_dummies <pandas.core.strings.StringMethods.get_dummies>`.

Factorizing values
------------------

To encode 1-d values as an enumerated type use ``factorize``:

.. ipython:: python

   x = pd.Series(['A', 'A', np.nan, 'B', 3.14, np.inf])
   x
   labels, uniques = pd.factorize(x)
   labels
   uniques

Note that ``factorize`` is similar to ``numpy.unique``, but differs in its
handling of NaN:

.. note::
   The following ``numpy.unique`` will fail under Python 3 with a ``TypeError``
   because of an ordering bug. See also
   `Here <https://github.com/numpy/numpy/issues/641>`__

.. ipython:: python

   pd.factorize(x, sort=True)
   np.unique(x, return_inverse=True)[::-1]
