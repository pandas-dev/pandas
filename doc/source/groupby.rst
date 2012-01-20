.. currentmodule:: pandas
.. _groupby:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

*****************************
Group By: split-apply-combine
*****************************

By "group by" we are refer to a process involving one or more of the following
steps

 - **Splitting** the data into groups based on some criteria
 - **Applying** a function to each group independently
 - **Combining** the results into a data structure

Of these, the split step is the most straightforward. In fact, in many
situations you may wish to split the data set into groups and do something with
those groups yourself. In the apply step, we might wish to one of the
following:

 - **Aggregation**: computing a summary statistic (or statistics) about each
   group. Some examples:

    - Compute group sums or means
    - Compute group sizes / counts

 - **Transformation**: perform some group-specific computations and return a
   like-indexed. Some examples:

    - Standardizing data (zscore) within group
    - Filling NAs within groups with a value derived from each group

 - Some combination of the above: GroupBy will examine the results of the apply
   step and try to return a sensibly combined result if it doesn't fit into
   either of the above two categories

Since the set of object instance method on pandas data structures are generally
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

.. _groupby.split:

Splitting an object into groups
-------------------------------

pandas objects can be split on any of their axes. The abstract definition of
grouping is to provide a mapping of labels to group names. To create a GroupBy
object (more on what the GroupBy object is later), you do the following:

.. code-block:: ipython

   # default is axis=0
   >>> grouped = obj.groupby(key)
   >>> grouped = obj.groupby(key, axis=1)
   >>> grouped = obj.groupby([key1, key2])

The mapping can be specified many different ways:

  - A Python function, to be called on each of the axis labels
  - A list or NumPy array of the same length as the selected axis
  - A dict or Series, providing a ``label -> group name`` mapping
  - For DataFrame objects, a string indicating a column to be used to group. Of
    course ``df.groupby('A')`` is just syntactic sugar for
    ``df.groupby(df['A'])``, but it makes life simpler
  - A list of any of the above things

Collectively we refer to the grouping objects as the **keys**. For example,
consider the following DataFrame:

.. ipython:: python

   df = DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
                          'foo', 'bar', 'foo', 'foo'],
                   'B' : ['one', 'one', 'two', 'three',
                          'two', 'two', 'one', 'three'],
                   'C' : randn(8), 'D' : randn(8)})
   df

We could naturally group by either the ``A`` or ``B`` columns or both:

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

Note that **no splitting occurs** until it's needed. Creating the GroupBy object
only verifies that you've passed a valid mapping.

.. note::

   Many kinds of complicated data manipulations can be expressed in terms of
   GroupBy operations (though can't be guaranteed to be the most
   efficient). You can get quite creative with the label mapping functions.

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

By default the group keys are sorted during the groupby operation. You may
however pass ``sort``=``False`` for potential speedups:

.. ipython:: python

   df2 = DataFrame({'X' : ['B', 'B', 'A', 'A'], 'Y' : [1, 2, 3, 4]})
   df2.groupby(['X'], sort=True).sum()
   df2.groupby(['X'], sort=False).sum()

.. _groupby.multiindex:

GroupBy with MultiIndex
~~~~~~~~~~~~~~~~~~~~~~~

With :ref:`hierarchically-indexed data <indexing.hierarchical>`, it's quite
natural to group by one of the levels of the hierarchy.

.. ipython:: python
   :suppress:

   arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
             ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
   tuples = zip(*arrays)
   tuples
   index = MultiIndex.from_tuples(tuples, names=['first', 'second'])
   s = Series(randn(8), index=index)

.. ipython:: python

   s
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

Also as of v0.6, grouping with multiple levels is supported.

.. ipython:: python
   :suppress:

   arrays = [['bar', 'bar', 'baz', 'baz', 'foo', 'foo', 'qux', 'qux'],
             ['doo', 'doo', 'bee', 'bee', 'bop', 'bop', 'bop', 'bop'],
             ['one', 'two', 'one', 'two', 'one', 'two', 'one', 'two']]
   tuples = zip(*arrays)
   index = MultiIndex.from_tuples(tuples, names=['first', 'second', 'third'])
   s = Series(randn(8), index=index)

.. ipython:: python

   s
   s.groupby(level=['first','second']).sum()

More on the ``sum`` function and aggregation later.

DataFrame column selection in GroupBy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have created the GroupBy object from a DataFrame, for example, you
might want to do something different for each of the columns. Thus, using
``[]`` similar to getting a column from a DataFrame, you can do:

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
natural and functions similarly to ``itertools.groupby``:

.. ipython::

   In [4]: grouped = df.groupby('A')

   In [5]: for name, group in grouped:
      ...:        print name
      ...:        print group
      ...:

In the case of grouping by multiple keys, the group name will be a tuple:

.. ipython::

   In [5]: for name, group in df.groupby(['A', 'B']):
      ...:        print name
      ...:        print group
      ...:

It's standard Python-fu but remember you can unpack the tuple in the for loop
statement if you wish: ``for (k1, k2), group in grouped:``.

.. _groupby.aggregate:

Aggregation
-----------

Once the GroupBy object has been created, several methods are available to
perform a computation on the grouped data. An obvious one is aggregation via
the ``aggregate`` or equivalently ``agg`` method:

.. ipython:: python

   grouped = df.groupby('A')
   grouped.aggregate(np.sum)

   grouped = df.groupby(['A', 'B'])
   grouped.aggregate(np.sum)

As you can see, the result of the aggregation will have the group names as the
new index along the grouped axis. In the case of multiple keys, the result is a
:ref:`MultiIndex <indexing.hierarchical>` by default, though this can be
changed by using the ``as_index`` option:

.. ipython:: python

   grouped = df.groupby(['A', 'B'], as_index=False)
   grouped.aggregate(np.sum)

   df.groupby('A', as_index=False).sum()

Note that you could use the ``reset_index`` DataFrame function to achieve the
same result as the column names are stored in the resulting ``MultiIndex``:

.. ipython:: python

   df.groupby(['A', 'B']).sum().reset_index()

.. _groupby.aggregate.multifunc:

Applying multiple functions at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With grouped Series you can also pass a list or dict of functions to do
aggregation with, outputting a DataFrame:

.. ipython:: python

   grouped = df.groupby('A')
   grouped['C'].agg([np.sum, np.mean, np.std])

If a dict is passed, the keys will be used to name the columns. Otherwise the
function's name (stored in the function object) will be used.

.. ipython:: python

   grouped['D'].agg({'result1' : np.sum,
                     'result2' : np.mean})

On a grouped DataFrame, you can pass a list of functions to apply to each
column, which produces an aggregated result with a hierarchical index:

.. ipython:: python

   grouped.agg([np.sum, np.mean, np.std])

Passing a dict of functions has different behavior by default, see the next
section.

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

.. _groupby.aggregate.cython:

Cython-optimized aggregation functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some common aggregations, currently only ``sum``, ``mean``, and ``std``, have
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
as the one being grouped. Thus, the passed transform function should return a
result that is the same size as the group chunk. For example, suppose we wished
to standardize a data set within a group:

.. ipython:: python

   tsdf = DataFrame(randn(1000, 3),
                    index=DateRange('1/1/2000', periods=1000),
                    columns=['A', 'B', 'C'])
   tsdf

   zscore = lambda x: (x - x.mean()) / x.std()
   transformed = tsdf.groupby(lambda x: x.year).transform(zscore)

We would expect the result to now have mean 0 and standard deviation 1 within
each group, which we can easily check:

.. ipython:: python

   grouped = transformed.groupby(lambda x: x.year)

   # OK, close enough to zero
   grouped.mean()
   grouped.std()

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

   tsdf.ix[::2] = np.nan
   grouped = tsdf.groupby(lambda x: x.year)
   grouped.fillna(method='pad')

In this example, we chopped the collection of time series into yearly chunks
then independently called :ref:`fillna <missing_data.fillna>` on the
groups.

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
       ....:     return DataFrame({'original' : group,
       ....:                       'demeaned' : group - group.mean()})
       ....:

    In [11]: grouped.apply(f)


Other useful features
---------------------

Automatic exclusion of "nuisance" columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Again consider the example DataFrame we've been looking at:

.. ipython:: python

   df

Supposed we wished to compute the standard deviation grouped by the ``A``
column. There is a slight problem, namely that we don't care about the data in
column ``B``. We refer to this as a "nuisance" column. If the passed
aggregation function can't be applied to some columns, the troublesome columns
will be (silently) dropped. Thus, this does not pose any problems:

.. ipython:: python

   df.groupby('A').std()

NA group handling
~~~~~~~~~~~~~~~~~

If there are any NaN values in the grouping key, these will be automatically
excluded. So there will never be an "NA group". This was not the case in older
versions of pandas, but users were generally discarding the NA group anyway
(and supporting it was an implementation headache).
