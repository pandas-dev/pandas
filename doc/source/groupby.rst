.. currentmodule:: pandas
.. _groupby:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   np.set_printoptions(precision=4, suppress=True)
   import pandas as pd
   pd.options.display.max_rows = 15
   import matplotlib
   # matplotlib.style.use('default')
   import matplotlib.pyplot as plt
   plt.close('all')
   from collections import OrderedDict

*****************************
Group By: split-apply-combine
*****************************

By "group by" we are referring to a process involving one or more of the following
steps:

* **Splitting** the data into groups based on some criteria.
* **Applying** a function to each group independently.
* **Combining** the results into a data structure.

Out of these, the split step is the most straightforward. In fact, in many
situations we may wish to split the data set into groups and do something with
those groups. In the apply step, we might wish to one of the
following:

* **Aggregation**: compute a summary statistic (or statistics) for each
  group. Some examples:

    * Compute group sums or means.
    * Compute group sizes / counts.

* **Transformation**: perform some group-specific computations and return a
  like-indexed object. Some examples:

    * Standardize data (zscore) within a group.
    * Filling NAs within groups with a value derived from each group.

* **Filtration**: discard some groups, according to a group-wise computation
  that evaluates True or False. Some examples:

    * Discard data that belongs to groups with only a few members.
    * Filter out data based on the group sum or mean.

* Some combination of the above: GroupBy will examine the results of the apply
  step and try to return a sensibly combined result if it doesn't fit into
  either of the above two categories.

Since the set of object instance methods on pandas data structures are generally
rich and expressive, we often simply want to invoke, say, a DataFrame function
on each group. The name GroupBy should be quite familiar to those who have used
a SQL-based tool (or ``itertools``), in which you can write code like:

.. code-block:: sql

   SELECT Column1, Column2, mean(Column3), sum(Column4)
   FROM SomeTable
   GROUP BY Column1, Column2

We aim to make operations like this natural and easy to express using
pandas. We'll address each area of GroupBy functionality then provide some
non-trivial examples / use cases.

See the :ref:`cookbook<cookbook.grouping>` for some advanced strategies.

.. _groupby.split:

Splitting an object into groups
-------------------------------

pandas objects can be split on any of their axes. The abstract definition of
grouping is to provide a mapping of labels to group names. To create a GroupBy
object (more on what the GroupBy object is later), you may do the following:

.. code-block:: ipython

   # default is axis=0
   >>> grouped = obj.groupby(key)
   >>> grouped = obj.groupby(key, axis=1)
   >>> grouped = obj.groupby([key1, key2])

The mapping can be specified many different ways:

* A Python function, to be called on each of the axis labels.
* A list or NumPy array of the same length as the selected axis.
* A dict or ``Series``, providing a ``label -> group name`` mapping.
* For ``DataFrame`` objects, a string indicating a column to be used to group.
  Of course ``df.groupby('A')`` is just syntactic sugar for
  ``df.groupby(df['A'])``, but it makes life simpler.
* For ``DataFrame`` objects, a string indicating an index level to be used to
  group.
* A list of any of the above things.

Collectively we refer to the grouping objects as the **keys**. For example,
consider the following ``DataFrame``:

.. note::

   .. versionadded:: 0.20

   A string passed to ``groupby`` may refer to either a column or an index level.
   If a string matches both a column name and an index level name then a warning is
   issued and the column takes precedence. This will result in an ambiguity error
   in a future version.

.. ipython:: python

   df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
                             'foo', 'bar', 'foo', 'foo'],
                      'B' : ['one', 'one', 'two', 'three',
                             'two', 'two', 'one', 'three'],
                      'C' : np.random.randn(8),
                      'D' : np.random.randn(8)})
   df

On a DataFrame, we obtain a GroupBy object by calling :meth:`~DataFrame.groupby`.
We could naturally group by either the ``A`` or ``B`` columns, or both:

.. ipython:: python

   grouped = df.groupby('A')
   grouped = df.groupby(['A', 'B'])

These will split the DataFrame on its index (rows). We could also split by the
columns:

.. ipython::

    In [4]: def get_letter_type(letter):
       ...:     if letter.lower() in 'aeiou':
       ...:         return 'vowel'
       ...:     else:
       ...:         return 'consonant'
       ...:

    In [5]: grouped = df.groupby(get_letter_type, axis=1)

pandas :class:`~pandas.Index` objects support duplicate values. If a
non-unique index is used as the group key in a groupby operation, all values
for the same index value will be considered to be in one group and thus the
output of aggregation functions will only contain unique index values:

.. ipython:: python

   lst = [1, 2, 3, 1, 2, 3]

   s = pd.Series([1, 2, 3, 10, 20, 30], lst)

   grouped = s.groupby(level=0)

   grouped.first()

   grouped.last()

   grouped.sum()

Note that **no splitting occurs** until it's needed. Creating the GroupBy object
only verifies that you've passed a valid mapping.

.. note::

   Many kinds of complicated data manipulations can be expressed in terms of
   GroupBy operations (though can't be guaranteed to be the most
   efficient). You can get quite creative with the label mapping functions.

.. _groupby.sorting:

GroupBy sorting
~~~~~~~~~~~~~~~~~~~~~~~~~

By default the group keys are sorted during the ``groupby`` operation. You may however pass ``sort=False`` for potential speedups:

.. ipython:: python

   df2 = pd.DataFrame({'X' : ['B', 'B', 'A', 'A'], 'Y' : [1, 2, 3, 4]})
   df2.groupby(['X']).sum()
   df2.groupby(['X'], sort=False).sum()


Note that ``groupby`` will preserve the order in which *observations* are sorted *within* each group.
For example, the groups created by ``groupby()`` below are in the order they appeared in the original ``DataFrame``:

.. ipython:: python

   df3 = pd.DataFrame({'X' : ['A', 'B', 'A', 'B'], 'Y' : [1, 4, 3, 2]})
   df3.groupby(['X']).get_group('A')

   df3.groupby(['X']).get_group('B')



.. _groupby.attributes:

GroupBy object attributes
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``groups`` attribute is a dict whose keys are the computed unique groups
and corresponding values being the axis labels belonging to each group. In the
above example we have:

.. ipython:: python

   df.groupby('A').groups
   df.groupby(get_letter_type, axis=1).groups

Calling the standard Python ``len`` function on the GroupBy object just returns
the length of the ``groups`` dict, so it is largely just a convenience:

.. ipython:: python

   grouped = df.groupby(['A', 'B'])
   grouped.groups
   len(grouped)


.. _groupby.tabcompletion:

``GroupBy`` will tab complete column names (and other attributes):

.. ipython:: python
   :suppress:

   n = 10
   weight = np.random.normal(166, 20, size=n)
   height = np.random.normal(60, 10, size=n)
   time = pd.date_range('1/1/2000', periods=n)
   gender = np.random.choice(['male', 'female'], size=n)
   df = pd.DataFrame({'height': height, 'weight': weight,
                      'gender': gender}, index=time)

.. ipython:: python

   df
   gb = df.groupby('gender')


.. ipython::

   @verbatim
   In [1]: gb.<TAB>
   gb.agg        gb.boxplot    gb.cummin     gb.describe   gb.filter     gb.get_group  gb.height     gb.last       gb.median     gb.ngroups    gb.plot       gb.rank       gb.std        gb.transform
   gb.aggregate  gb.count      gb.cumprod    gb.dtype      gb.first      gb.groups     gb.hist       gb.max        gb.min        gb.nth        gb.prod       gb.resample   gb.sum        gb.var
   gb.apply      gb.cummax     gb.cumsum     gb.fillna     gb.gender     gb.head       gb.indices    gb.mean       gb.name       gb.ohlc       gb.quantile   gb.size       gb.tail       gb.weight

.. _groupby.multiindex:

GroupBy with MultiIndex
~~~~~~~~~~~~~~~~~~~~~~~

With :ref:`hierarchically-indexed data <advanced.hierarchical>`, it's quite
natural to group by one of the levels of the hierarchy.

Let's create a Series with a two-level ``MultiIndex``.

.. ipython:: python


   arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
             ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
   index = pd.MultiIndex.from_arrays(arrays, names=['first', 'second'])
   s = pd.Series(np.random.randn(8), index=index)
   s

We can then group by one of the levels in ``s``.

.. ipython:: python

   grouped = s.groupby(level=0)
   grouped.sum()

If the MultiIndex has names specified, these can be passed instead of the level
number:

.. ipython:: python

   s.groupby(level='second').sum()

The aggregation functions such as ``sum`` will take the level parameter
directly. Additionally, the resulting index will be named according to the
chosen level:

.. ipython:: python

   s.sum(level='second')

Grouping with multiple levels is supported.

.. ipython:: python
   :suppress:

   arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
             ['doo', 'doo', 'bee', 'bee', 'bop', 'bop', 'bop', 'bop'],
             ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
   tuples = list(zip(*arrays))
   index = pd.MultiIndex.from_tuples(tuples, names=['first', 'second', 'third'])
   s = pd.Series(np.random.randn(8), index=index)

.. ipython:: python

   s
   s.groupby(level=['first', 'second']).sum()

.. versionadded:: 0.20

Index level names may be supplied as keys.

.. ipython:: python

   s.groupby(['first', 'second']).sum()

More on the ``sum`` function and aggregation later.

Grouping DataFrame with Index Levels and Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A DataFrame may be grouped by a combination of columns and index levels by
specifying the column names as strings and the index levels as ``pd.Grouper``
objects.

.. ipython:: python

   arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
             ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]

   index = pd.MultiIndex.from_arrays(arrays, names=['first', 'second'])

   df = pd.DataFrame({'A': [1, 1, 1, 1, 2, 2, 3, 3],
                      'B': np.arange(8)},
                     index=index)

   df

The following example groups ``df`` by the ``second`` index level and
the ``A`` column.

.. ipython:: python

   df.groupby([pd.Grouper(level=1), 'A']).sum()

Index levels may also be specified by name.

.. ipython:: python

   df.groupby([pd.Grouper(level='second'), 'A']).sum()

.. versionadded:: 0.20

Index level names may be specified as keys directly to ``groupby``.

.. ipython:: python

   df.groupby(['second', 'A']).sum()

DataFrame column selection in GroupBy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have created the GroupBy object from a DataFrame, you might want to do
something different for each of the columns. Thus, using ``[]`` similar to
getting a column from a DataFrame, you can do:

.. ipython:: python
   :suppress:

   df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
                             'foo', 'bar', 'foo', 'foo'],
                      'B' : ['one', 'one', 'two', 'three',
                             'two', 'two', 'one', 'three'],
                      'C' : np.random.randn(8),
                      'D' : np.random.randn(8)})

.. ipython:: python

   grouped = df.groupby(['A'])
   grouped_C = grouped['C']
   grouped_D = grouped['D']

This is mainly syntactic sugar for the alternative and much more verbose:

.. ipython:: python

   df['C'].groupby(df['A'])

Additionally this method avoids recomputing the internal grouping information
derived from the passed key.

.. _groupby.iterating:

Iterating through groups
------------------------

With the GroupBy object in hand, iterating through the grouped data is very
natural and functions similarly to :py:func:`itertools.groupby`:

.. ipython::

   In [4]: grouped = df.groupby('A')

   In [5]: for name, group in grouped:
      ...:        print(name)
      ...:        print(group)
      ...:

In the case of grouping by multiple keys, the group name will be a tuple:

.. ipython::

   In [5]: for name, group in df.groupby(['A', 'B']):
      ...:        print(name)
      ...:        print(group)
      ...:

It's standard Python-fu but remember you can unpack the tuple in the for loop
statement if you wish: ``for (k1, k2), group in grouped:``.

Selecting a group
-----------------

A single group can be selected using
:meth:`~pandas.core.groupby.DataFrameGroupBy.get_group`:

.. ipython:: python

   grouped.get_group('bar')

Or for an object grouped on multiple columns:

.. ipython:: python

   df.groupby(['A', 'B']).get_group(('bar', 'one'))

.. _groupby.aggregate:

Aggregation
-----------

Once the GroupBy object has been created, several methods are available to
perform a computation on the grouped data. These operations are similar to the
:ref:`aggregating API <basics.aggregate>`, :ref:`window functions API <stats.aggregate>`,
and :ref:`resample API <timeseries.aggregate>`.

An obvious one is aggregation via the
:meth:`~pandas.core.groupby.DataFrameGroupBy.aggregate` or equivalently
:meth:`~pandas.core.groupby.DataFrameGroupBy.agg` method:

.. ipython:: python

   grouped = df.groupby('A')
   grouped.aggregate(np.sum)

   grouped = df.groupby(['A', 'B'])
   grouped.aggregate(np.sum)

As you can see, the result of the aggregation will have the group names as the
new index along the grouped axis. In the case of multiple keys, the result is a
:ref:`MultiIndex <advanced.hierarchical>` by default, though this can be
changed by using the ``as_index`` option:

.. ipython:: python

   grouped = df.groupby(['A', 'B'], as_index=False)
   grouped.aggregate(np.sum)

   df.groupby('A', as_index=False).sum()

Note that you could use the ``reset_index`` DataFrame function to achieve the
same result as the column names are stored in the resulting ``MultiIndex``:

.. ipython:: python

   df.groupby(['A', 'B']).sum().reset_index()

Another simple aggregation example is to compute the size of each group.
This is included in GroupBy as the ``size`` method. It returns a Series whose
index are the group names and whose values are the sizes of each group.

.. ipython:: python

   grouped.size()

.. ipython:: python

   grouped.describe()

.. note::

   Aggregation functions **will not** return the groups that you are aggregating over
   if they are named *columns*, when ``as_index=True``, the default. The grouped columns will
   be the **indices** of the returned object.

   Passing ``as_index=False`` **will** return the groups that you are aggregating over, if they are
   named *columns*.

Aggregating functions are the ones that reduce the dimension of the returned objects.
Some common aggregating functions are tabulated below:

.. csv-table::
    :header: "Function", "Description"
    :widths: 20, 80
    :delim: ;

	:meth:`~pd.core.groupby.DataFrameGroupBy.mean`;Compute mean of groups
	:meth:`~pd.core.groupby.DataFrameGroupBy.sum`;Compute sum of group values
	:meth:`~pd.core.groupby.DataFrameGroupBy.size`;Compute group sizes
	:meth:`~pd.core.groupby.DataFrameGroupBy.count`;Compute count of group
	:meth:`~pd.core.groupby.DataFrameGroupBy.std`;Standard deviation of groups
	:meth:`~pd.core.groupby.DataFrameGroupBy.var`;Compute variance of groups
	:meth:`~pd.core.groupby.DataFrameGroupBy.sem`;Standard error of the mean of groups
	:meth:`~pd.core.groupby.DataFrameGroupBy.describe`;Generates descriptive statistics
	:meth:`~pd.core.groupby.DataFrameGroupBy.first`;Compute first of group values
	:meth:`~pd.core.groupby.DataFrameGroupBy.last`;Compute last of group values
	:meth:`~pd.core.groupby.DataFrameGroupBy.nth`;Take nth value, or a subset if n is a list
	:meth:`~pd.core.groupby.DataFrameGroupBy.min`;Compute min of group values
	:meth:`~pd.core.groupby.DataFrameGroupBy.max`;Compute max of group values


The aggregating functions above will exclude NA values. Any function which
reduces a :class:`Series` to a scalar value is an aggregation function and will work,
a trivial example is ``df.groupby('A').agg(lambda ser: 1)``. Note that
:meth:`~pd.core.groupby.DataFrameGroupBy.nth` can act as a reducer *or* a
filter, see :ref:`here <groupby.nth>`.

.. _groupby.aggregate.multifunc:

Applying multiple functions at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With grouped ``Series`` you can also pass a list or dict of functions to do
aggregation with, outputting a DataFrame:

.. ipython:: python

   grouped = df.groupby('A')
   grouped['C'].agg([np.sum, np.mean, np.std])

On a grouped ``DataFrame``, you can pass a list of functions to apply to each
column, which produces an aggregated result with a hierarchical index:

.. ipython:: python

   grouped.agg([np.sum, np.mean, np.std])


The resulting aggregations are named for the functions themselves. If you
need to rename, then you can add in a chained operation for a ``Series`` like this:

.. ipython:: python

   (grouped['C'].agg([np.sum, np.mean, np.std])
                .rename(columns={'sum': 'foo',
                                 'mean': 'bar',
                                 'std': 'baz'})
   )

For a grouped ``DataFrame``, you can rename in a similar manner:

.. ipython:: python

   (grouped.agg([np.sum, np.mean, np.std])
           .rename(columns={'sum': 'foo',
                            'mean': 'bar',
                            'std': 'baz'})
    )


Applying different functions to DataFrame columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By passing a dict to ``aggregate`` you can apply a different aggregation to the
columns of a DataFrame:

.. ipython:: python

   grouped.agg({'C' : np.sum,
                'D' : lambda x: np.std(x, ddof=1)})

The function names can also be strings. In order for a string to be valid it
must be either implemented on GroupBy or available via :ref:`dispatching
<groupby.dispatch>`:

.. ipython:: python

   grouped.agg({'C' : 'sum', 'D' : 'std'})

.. note::

    If you pass a dict to ``aggregate``, the ordering of the output columns is
    non-deterministic. If you want to be sure the output columns will be in a specific
    order, you can use an ``OrderedDict``.  Compare the output of the following two commands:

.. ipython:: python

   grouped.agg({'D': 'std', 'C': 'mean'})
   grouped.agg(OrderedDict([('D', 'std'), ('C', 'mean')]))

.. _groupby.aggregate.cython:

Cython-optimized aggregation functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some common aggregations, currently only ``sum``, ``mean``, ``std``, and ``sem``, have
optimized Cython implementations:

.. ipython:: python

   df.groupby('A').sum()
   df.groupby(['A', 'B']).mean()

Of course ``sum`` and ``mean`` are implemented on pandas objects, so the above
code would work even without the special versions via dispatching (see below).

.. _groupby.transform:

Transformation
--------------

The ``transform`` method returns an object that is indexed the same (same size)
as the one being grouped. The transform function must:

* Return a result that is either the same size as the group chunk or
  broadcastable to the size of the group chunk (e.g., a scalar,
  ``grouped.transform(lambda x: x.iloc[-1])``).
* Operate column-by-column on the group chunk.  The transform is applied to
  the first group chunk using chunk.apply.
* Not perform in-place operations on the group chunk. Group chunks should
  be treated as immutable, and changes to a group chunk may produce unexpected
  results. For example, when using ``fillna``, ``inplace`` must be ``False``
  (``grouped.transform(lambda x: x.fillna(inplace=False))``).
* (Optionally) operates on the entire group chunk. If this is supported, a
  fast path is used starting from the *second* chunk.

For example, suppose we wished to standardize the data within each group:

.. ipython:: python

   index = pd.date_range('10/1/1999', periods=1100)
   ts = pd.Series(np.random.normal(0.5, 2, 1100), index)
   ts = ts.rolling(window=100,min_periods=100).mean().dropna()

   ts.head()
   ts.tail()
   key = lambda x: x.year
   zscore = lambda x: (x - x.mean()) / x.std()
   transformed = ts.groupby(key).transform(zscore)

We would expect the result to now have mean 0 and standard deviation 1 within
each group, which we can easily check:

.. ipython:: python

   # Original Data
   grouped = ts.groupby(key)
   grouped.mean()
   grouped.std()

   # Transformed Data
   grouped_trans = transformed.groupby(key)
   grouped_trans.mean()
   grouped_trans.std()

We can also visually compare the original and transformed data sets.

.. ipython:: python

   compare = pd.DataFrame({'Original': ts, 'Transformed': transformed})

   @savefig groupby_transform_plot.png
   compare.plot()

Transformation functions that have lower dimension outputs are broadcast to
match the shape of the input array.

.. ipython:: python

   data_range = lambda x: x.max() - x.min()
   ts.groupby(key).transform(data_range)

Alternatively, the built-in methods could be used to produce the same outputs.

.. ipython:: python

   ts.groupby(key).transform('max') - ts.groupby(key).transform('min')

Another common data transform is to replace missing data with the group mean.

.. ipython:: python
   :suppress:

   cols = ['A', 'B', 'C']
   values = np.random.randn(1000, 3)
   values[np.random.randint(0, 1000, 100), 0] = np.nan
   values[np.random.randint(0, 1000, 50), 1] = np.nan
   values[np.random.randint(0, 1000, 200), 2] = np.nan
   data_df = pd.DataFrame(values, columns=cols)

.. ipython:: python

   data_df

   countries = np.array(['US', 'UK', 'GR', 'JP'])
   key = countries[np.random.randint(0, 4, 1000)]

   grouped = data_df.groupby(key)

   # Non-NA count in each group
   grouped.count()

   f = lambda x: x.fillna(x.mean())

   transformed = grouped.transform(f)

We can verify that the group means have not changed in the transformed data
and that the transformed data contains no NAs.

.. ipython:: python

   grouped_trans = transformed.groupby(key)

   grouped.mean() # original group means
   grouped_trans.mean() # transformation did not change group means

   grouped.count() # original has some missing data points
   grouped_trans.count() # counts after transformation
   grouped_trans.size() # Verify non-NA count equals group size

.. note::

   Some functions will automatically transform the input when applied to a
   GroupBy object, but returning an object of the same shape as the original.
   Passing ``as_index=False`` will not affect these transformation methods.

   For example: ``fillna, ffill, bfill, shift.``.

   .. ipython:: python

      grouped.ffill()


.. _groupby.transform.window_resample:

New syntax to window and resample operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. versionadded:: 0.18.1

Working with the resample, expanding or rolling operations on the groupby
level used to require the application of helper functions. However,
now it is possible to use ``resample()``, ``expanding()`` and
``rolling()`` as methods on groupbys.

The example below will apply the ``rolling()`` method on the samples of
the column B based on the groups of column A.

.. ipython:: python

   df_re = pd.DataFrame({'A': [1] * 10 + [5] * 10,
                         'B': np.arange(20)})
   df_re

   df_re.groupby('A').rolling(4).B.mean()


The ``expanding()`` method will accumulate a given operation
(``sum()`` in the example) for all the members of each particular
group.

.. ipython:: python

   df_re.groupby('A').expanding().sum()


Suppose you want to use the ``resample()`` method to get a daily
frequency in each group of your dataframe and wish to complete the
missing values with the ``ffill()`` method.

.. ipython:: python

   df_re = pd.DataFrame({'date': pd.date_range(start='2016-01-01',
                                 periods=4,
                         freq='W'),
                        'group': [1, 1, 2, 2],
                        'val': [5, 6, 7, 8]}).set_index('date')
   df_re

   df_re.groupby('group').resample('1D').ffill()

.. _groupby.filter:

Filtration
----------

The ``filter`` method returns a subset of the original object. Suppose we
want to take only elements that belong to groups with a group sum greater
than 2.

.. ipython:: python

   sf = pd.Series([1, 1, 2, 3, 3, 3])
   sf.groupby(sf).filter(lambda x: x.sum() > 2)

The argument of ``filter`` must be a function that, applied to the group as a
whole, returns ``True`` or ``False``.

Another useful operation is filtering out elements that belong to groups
with only a couple members.

.. ipython:: python

   dff = pd.DataFrame({'A': np.arange(8), 'B': list('aabbbbcc')})
   dff.groupby('B').filter(lambda x: len(x) > 2)

Alternatively, instead of dropping the offending groups, we can return a
like-indexed objects where the groups that do not pass the filter are filled
with NaNs.

.. ipython:: python

   dff.groupby('B').filter(lambda x: len(x) > 2, dropna=False)

For DataFrames with multiple columns, filters should explicitly specify a column as the filter criterion.

.. ipython:: python

   dff['C'] = np.arange(8)
   dff.groupby('B').filter(lambda x: len(x['C']) > 2)

.. note::

   Some functions when applied to a groupby object will act as a **filter** on the input, returning
   a reduced shape of the original (and potentially eliminating groups), but with the index unchanged.
   Passing ``as_index=False`` will not affect these transformation methods.

   For example: ``head, tail``.

   .. ipython:: python

      dff.groupby('B').head(2)


.. _groupby.dispatch:

Dispatching to instance methods
-------------------------------

When doing an aggregation or transformation, you might just want to call an
instance method on each data group. This is pretty easy to do by passing lambda
functions:

.. ipython:: python

   grouped = df.groupby('A')
   grouped.agg(lambda x: x.std())

But, it's rather verbose and can be untidy if you need to pass additional
arguments. Using a bit of metaprogramming cleverness, GroupBy now has the
ability to "dispatch" method calls to the groups:

.. ipython:: python

   grouped.std()

What is actually happening here is that a function wrapper is being
generated. When invoked, it takes any passed arguments and invokes the function
with any arguments on each group (in the above example, the ``std``
function). The results are then combined together much in the style of ``agg``
and ``transform`` (it actually uses ``apply`` to infer the gluing, documented
next). This enables some operations to be carried out rather succinctly:

.. ipython:: python

   tsdf = pd.DataFrame(np.random.randn(1000, 3),
                       index=pd.date_range('1/1/2000', periods=1000),
                       columns=['A', 'B', 'C'])
   tsdf.iloc[::2] = np.nan
   grouped = tsdf.groupby(lambda x: x.year)
   grouped.fillna(method='pad')

In this example, we chopped the collection of time series into yearly chunks
then independently called :ref:`fillna <missing_data.fillna>` on the
groups.

The ``nlargest`` and ``nsmallest`` methods work on ``Series`` style groupbys:

.. ipython:: python

   s = pd.Series([9, 8, 7, 5, 19, 1, 4.2, 3.3])
   g = pd.Series(list('abababab'))
   gb = s.groupby(g)
   gb.nlargest(3)
   gb.nsmallest(3)

.. _groupby.apply:

Flexible ``apply``
------------------

Some operations on the grouped data might not fit into either the aggregate or
transform categories. Or, you may simply want GroupBy to infer how to combine
the results. For these, use the ``apply`` function, which can be substituted
for both ``aggregate`` and ``transform`` in many standard use cases. However,
``apply`` can handle some exceptional use cases, for example:

.. ipython:: python

   df
   grouped = df.groupby('A')

   # could also just call .describe()
   grouped['C'].apply(lambda x: x.describe())

The dimension of the returned result can also change:

.. ipython::

    In [8]: grouped = df.groupby('A')['C']

    In [10]: def f(group):
       ....:     return pd.DataFrame({'original' : group,
       ....:                          'demeaned' : group - group.mean()})
       ....:

    In [11]: grouped.apply(f)

``apply`` on a Series can operate on a returned value from the applied function,
that is itself a series, and possibly upcast the result to a DataFrame:

.. ipython:: python

    def f(x):
      return pd.Series([ x, x**2 ], index = ['x', 'x^2'])
    s = pd.Series(np.random.rand(5))
    s
    s.apply(f)


.. note::

   ``apply`` can act as a reducer, transformer, *or* filter function, depending on exactly what is passed to it.
   So depending on the path taken, and exactly what you are grouping. Thus the grouped columns(s) may be included in
   the output as well as set the indices.

.. warning::

    In the current implementation apply calls func twice on the
    first group to decide whether it can take a fast or slow code
    path. This can lead to unexpected behavior if func has
    side-effects, as they will take effect twice for the first
    group.

    .. ipython:: python

        d = pd.DataFrame({"a":["x", "y"], "b":[1,2]})
        def identity(df):
            print(df)
            return df

        d.groupby("a").apply(identity)


Other useful features
---------------------

Automatic exclusion of "nuisance" columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Again consider the example DataFrame we've been looking at:

.. ipython:: python

   df

Suppose we wish to compute the standard deviation grouped by the ``A``
column. There is a slight problem, namely that we don't care about the data in
column ``B``. We refer to this as a "nuisance" column. If the passed
aggregation function can't be applied to some columns, the troublesome columns
will be (silently) dropped. Thus, this does not pose any problems:

.. ipython:: python

   df.groupby('A').std()

Note that ``df.groupby('A').colname.std().`` is more efficient than
``df.groupby('A').std().colname``, so if the result of an aggregation function
is only interesting over one column (here ``colname``), it may be filtered
*before* applying the aggregation function.

.. _groupby.observed:

Handling of (un)observed Categorical values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using a ``Categorical`` grouper (as a single grouper, or as part of multiple groupers), the ``observed`` keyword
controls whether to return a cartesian product of all possible groupers values (``observed=False``) or only those
that are observed groupers (``observed=True``).

Show all values:

.. ipython:: python

   pd.Series([1, 1, 1]).groupby(pd.Categorical(['a', 'a', 'a'], categories=['a', 'b']), observed=False).count()

Show only the observed values:

.. ipython:: python

   pd.Series([1, 1, 1]).groupby(pd.Categorical(['a', 'a', 'a'], categories=['a', 'b']), observed=True).count()

The returned dtype of the grouped will *always* include *all* of the categories that were grouped.

.. ipython:: python

   s = pd.Series([1, 1, 1]).groupby(pd.Categorical(['a', 'a', 'a'], categories=['a', 'b']), observed=False).count()
   s.index.dtype

.. _groupby.missing:

NA and NaT group handling
~~~~~~~~~~~~~~~~~~~~~~~~~

If there are any NaN or NaT values in the grouping key, these will be
automatically excluded. In other words, there will never be an "NA group" or
"NaT group". This was not the case in older versions of pandas, but users were
generally discarding the NA group anyway (and supporting it was an
implementation headache).

Grouping with ordered factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Categorical variables represented as instance of pandas's ``Categorical`` class
can be used as group keys. If so, the order of the levels will be preserved:

.. ipython:: python

   data = pd.Series(np.random.randn(100))

   factor = pd.qcut(data, [0, .25, .5, .75, 1.])

   data.groupby(factor).mean()

.. _groupby.specify:

Grouping with a Grouper specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may need to specify a bit more data to properly group. You can
use the ``pd.Grouper`` to provide this local control.

.. ipython:: python

   import datetime

   df = pd.DataFrame({
            'Branch' : 'A A A A A A A B'.split(),
            'Buyer': 'Carl Mark Carl Carl Joe Joe Joe Carl'.split(),
            'Quantity': [1,3,5,1,8,1,9,3],
            'Date' : [
                datetime.datetime(2013,1,1,13,0),
                datetime.datetime(2013,1,1,13,5),
                datetime.datetime(2013,10,1,20,0),
                datetime.datetime(2013,10,2,10,0),
                datetime.datetime(2013,10,1,20,0),
                datetime.datetime(2013,10,2,10,0),
                datetime.datetime(2013,12,2,12,0),
                datetime.datetime(2013,12,2,14,0),
                ]
            })

   df

Groupby a specific column with the desired frequency. This is like resampling.

.. ipython:: python

   df.groupby([pd.Grouper(freq='1M',key='Date'),'Buyer']).sum()

You have an ambiguous specification in that you have a named index and a column
that could be potential groupers.

.. ipython:: python

   df = df.set_index('Date')
   df['Date'] = df.index + pd.offsets.MonthEnd(2)
   df.groupby([pd.Grouper(freq='6M',key='Date'),'Buyer']).sum()

   df.groupby([pd.Grouper(freq='6M',level='Date'),'Buyer']).sum()


Taking the first rows of each group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Just like for a DataFrame or Series you can call head and tail on a groupby:

.. ipython:: python

   df = pd.DataFrame([[1, 2], [1, 4], [5, 6]], columns=['A', 'B'])
   df

   g = df.groupby('A')
   g.head(1)

   g.tail(1)

This shows the first or last n rows from each group.

.. _groupby.nth:

Taking the nth row of each group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To select from a DataFrame or Series the nth item, use
:meth:`~pd.core.groupby.DataFrameGroupBy.nth`. This is a reduction method, and
will return a single row (or no row) per group if you pass an int for n:

.. ipython:: python

   df = pd.DataFrame([[1, np.nan], [1, 4], [5, 6]], columns=['A', 'B'])
   g = df.groupby('A')

   g.nth(0)
   g.nth(-1)
   g.nth(1)

If you want to select the nth not-null item, use the ``dropna`` kwarg. For a DataFrame this should be either ``'any'`` or ``'all'`` just like you would pass to dropna:

.. ipython:: python

   # nth(0) is the same as g.first()
   g.nth(0, dropna='any')
   g.first()

   # nth(-1) is the same as g.last()
   g.nth(-1, dropna='any')  # NaNs denote group exhausted when using dropna
   g.last()

   g.B.nth(0, dropna='all')

As with other methods, passing ``as_index=False``, will achieve a filtration, which returns the grouped row.

.. ipython:: python

   df = pd.DataFrame([[1, np.nan], [1, 4], [5, 6]], columns=['A', 'B'])
   g = df.groupby('A',as_index=False)

   g.nth(0)
   g.nth(-1)

You can also select multiple rows from each group by specifying multiple nth values as a list of ints.

.. ipython:: python

   business_dates = pd.date_range(start='4/1/2014', end='6/30/2014', freq='B')
   df = pd.DataFrame(1, index=business_dates, columns=['a', 'b'])
   # get the first, 4th, and last date index for each month
   df.groupby([df.index.year, df.index.month]).nth([0, 3, -1])

Enumerate group items
~~~~~~~~~~~~~~~~~~~~~

To see the order in which each row appears within its group, use the
``cumcount`` method:

.. ipython:: python

   dfg = pd.DataFrame(list('aaabba'), columns=['A'])
   dfg

   dfg.groupby('A').cumcount()

   dfg.groupby('A').cumcount(ascending=False)

.. _groupby.ngroup:

Enumerate groups
~~~~~~~~~~~~~~~~

.. versionadded:: 0.20.2

To see the ordering of the groups (as opposed to the order of rows
within a group given by ``cumcount``) you can use
:meth:`~pandas.core.groupby.DataFrameGroupBy.ngroup`.



Note that the numbers given to the groups match the order in which the
groups would be seen when iterating over the groupby object, not the
order they are first observed.

.. ipython:: python

   dfg = pd.DataFrame(list('aaabba'), columns=['A'])
   dfg

   dfg.groupby('A').ngroup()

   dfg.groupby('A').ngroup(ascending=False)

Plotting
~~~~~~~~

Groupby also works with some plotting methods.  For example, suppose we
suspect that some features in a DataFrame may differ by group, in this case,
the values in column 1 where the group is "B" are 3 higher on average.

.. ipython:: python

   np.random.seed(1234)
   df = pd.DataFrame(np.random.randn(50, 2))
   df['g'] = np.random.choice(['A', 'B'], size=50)
   df.loc[df['g'] == 'B', 1] += 3

We can easily visualize this with a boxplot:

.. ipython:: python
   :okwarning:

   @savefig groupby_boxplot.png
   df.groupby('g').boxplot()

The result of calling ``boxplot`` is a dictionary whose keys are the values
of our grouping column ``g`` ("A" and "B"). The values of the resulting dictionary
can be controlled by the ``return_type`` keyword of ``boxplot``.
See the :ref:`visualization documentation<visualization.box>` for more.

.. warning::

  For historical reasons, ``df.groupby("g").boxplot()`` is not equivalent
  to ``df.boxplot(by="g")``. See :ref:`here<visualization.box.return>` for
  an explanation.

.. _groupby.pipe:

Piping function calls
~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.21.0

Similar to the functionality provided by ``DataFrame`` and ``Series``, functions
that take ``GroupBy`` objects can be chained together using a ``pipe`` method to
allow for a cleaner, more readable syntax. To read about ``.pipe`` in general terms,
see :ref:`here <basics.pipe>`.

Combining ``.groupby`` and ``.pipe`` is often useful when you need to reuse
GroupBy objects.

As an example, imagine having a DataFrame with columns for stores, products,
revenue and quantity sold. We'd like to do a groupwise calculation of *prices*
(i.e. revenue/quantity) per store and per product. We could do this in a
multi-step operation, but expressing it in terms of piping can make the
code more readable. First we set the data:

.. ipython:: python

   import numpy as np
   n = 1000
   df = pd.DataFrame({'Store': np.random.choice(['Store_1', 'Store_2'], n),
                      'Product': np.random.choice(['Product_1',
                                                   'Product_2'], n),
                      'Revenue': (np.random.random(n)*50+10).round(2),
                      'Quantity': np.random.randint(1, 10, size=n)})
   df.head(2)

Now, to find prices per store/product, we can simply do:

.. ipython:: python

   (df.groupby(['Store', 'Product'])
      .pipe(lambda grp: grp.Revenue.sum()/grp.Quantity.sum())
      .unstack().round(2))

Piping can also be expressive when you want to deliver a grouped object to some
arbitrary function, for example:

.. code-block:: python

   (df.groupby(['Store', 'Product']).pipe(report_func)

where ``report_func`` takes a GroupBy object and creates a report
from that.

Examples
--------

Regrouping by factor
~~~~~~~~~~~~~~~~~~~~

Regroup columns of a DataFrame according to their sum, and sum the aggregated ones.

.. ipython:: python

   df = pd.DataFrame({'a':[1,0,0], 'b':[0,1,0], 'c':[1,0,0], 'd':[2,3,4]})
   df
   df.groupby(df.sum(), axis=1).sum()

.. _groupby.multicolumn_factorization:

Multi-column factorization
~~~~~~~~~~~~~~~~~~~~~~~~~~

By using :meth:`~pandas.core.groupby.DataFrameGroupBy.ngroup`, we can extract
information about the groups in a way similar to :func:`factorize` (as described
further in the :ref:`reshaping API <reshaping.factorize>`) but which applies
naturally to multiple columns of mixed type and different
sources. This can be useful as an intermediate categorical-like step
in processing, when the relationships between the group rows are more
important than their content, or as input to an algorithm which only
accepts the integer encoding. (For more information about support in
pandas for full categorical data, see the :ref:`Categorical
introduction <categorical>` and the
:ref:`API documentation <api.categorical>`.)

.. ipython:: python

    dfg = pd.DataFrame({"A": [1, 1, 2, 3, 2], "B": list("aaaba")})

    dfg

    dfg.groupby(["A", "B"]).ngroup()

    dfg.groupby(["A", [0, 0, 0, 1, 1]]).ngroup()

Groupby by Indexer to 'resample' data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Resampling produces new hypothetical samples (resamples) from already existing observed data or from a model that generates data. These new samples are similar to the pre-existing samples.

In order to resample to work on indices that are non-datetimelike, the following procedure can be utilized.

In the following examples, **df.index // 5** returns a binary array which is used to determine what gets selected for the groupby operation.

.. note:: The below example shows how we can downsample by consolidation of samples into fewer samples. Here by using **df.index // 5**, we are aggregating the samples in bins. By applying **std()** function, we aggregate the information contained in many samples into a small subset of values which is their standard deviation thereby reducing the number of samples.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10,2))
   df
   df.index // 5
   df.groupby(df.index // 5).std()

Returning a Series to propagate names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Group DataFrame columns, compute a set of metrics and return a named Series.
The Series name is used as the name for the column index. This is especially
useful in conjunction with reshaping operations such as stacking in which the
column index name will be used as the name of the inserted column:

.. ipython:: python

   df = pd.DataFrame({
            'a':  [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2],
            'b':  [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1],
            'c':  [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            'd':  [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1],
            })

   def compute_metrics(x):
       result = {'b_sum': x['b'].sum(), 'c_mean': x['c'].mean()}
       return pd.Series(result, name='metrics')

   result = df.groupby('a').apply(compute_metrics)

   result

   result.stack()
