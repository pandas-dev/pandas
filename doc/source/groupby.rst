.. currentmodule:: pandas
.. _groupby:

.. ipython:: python
   :suppress:

   import numpy as np
   np.random.seed(123456)
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

****************************
GroupBy: split-apply-combine
****************************

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
   like-indexed.
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
   index = MultiIndex.from_tuples(tuples)
   s = Series(randn(8), index=index)

.. ipython:: python

   s
   grouped = s.groupby(level=0)
   grouped.sum()

More on the ``sum`` function and aggregation later. Grouping with multiple
levels (as opposed to a single level) is not yet supported, though implementing
it is not difficult.

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
:ref:`MultiIndex <indexing.hierarchical>` by default.

Cython-optimized aggregation functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some common aggregations, currently only ``sum`` and ``mean``, have optimized
Cython implementations:

.. ipython:: python

   df.groupby('A').sum()
   df.groupby(['A', 'B']).mean()

Of course ``sum`` and ``mean`` are implemented on pandas objects, so the above
code would work even without the special versions via the "dispatching" feature
described below.

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
   grouped.mean()
   grouped.std()

Flexible ``apply``
------------------

Other useful features
---------------------

Invoking instance methods on groups, "dispatching"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Grouping single columns of a DataFrame
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Automatic exclusion of "nuisance" columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
