.. currentmodule:: pandas
.. _basics:

.. ipython:: python
   :suppress:

   import numpy as np
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

*****************************
Essential basic functionality
*****************************

Here we discuss a lot of the essential functionality common to the pandas data
structures. Here's how to create some of the objects used in the examples from
the previous section:

.. ipython:: python

   index = DateRange('1/1/2000', periods=8)
   s = Series(randn(5), index=['a', 'b', 'c', 'd', 'e'])
   df = DataFrame(randn(8, 3), index=index,
                  columns=['A', 'B', 'C'])
   wp = Panel(randn(2, 5, 4), items=['Item1', 'Item2'],
              major_axis=DateRange('1/1/2000', periods=5),
              minor_axis=['A', 'B', 'C', 'D'])

.. _basics.head_tail:

Head and Tail
-------------

To view a small sample of a Series or DataFrame object, use the ``head`` and
``tail`` methods. The default number of elements to display is five, but you
may pass a custom number.

.. ipython:: python

   long_series = Series(randn(1000))
   long_series.head()
   long_series.tail(3)

.. _basics.attrs:

Attributes and the raw ndarray(s)
---------------------------------

pandas objects have a number of attributes enabling you to access the metadata

  * **shape**: gives the axis dimensions of the object, consistent with ndarray
  * Axis labels

    * **Series**: *index* (only axis)
    * **DataFrame**: *index* (rows) and *columns*
    * **Panel**: *items*, *major_axis*, and *minor_axis*

Note, **these attributes can be safely assigned to**!

.. ipython:: python

   df[:2]
   df.columns = [x.lower() for x in df.columns]
   df

To get the actual data inside a data structure, one need only access the
**values** property:

.. ipython:: python

    s.values
    df.values
    wp.values

If a DataFrame or Panel contains homogeneously-typed data, the ndarray can
actually be modified in-place, and the changes will be reflected in the data
structure. For heterogeneous data (e.g. some of the DataFrame's columns are not
all the same dtype), this will not be the case. The values attribute itself,
unlike the axis labels, cannot be assigned to.

.. note::

    When working with heterogeneous data, the dtype of the resulting ndarray
    will be chosen to accommodate all of the data involved. For example, if
    strings are involved, the result will be of object dtype. If there are only
    floats and integers, the resulting array will be of float dtype.

.. _basics.binop:

Flexible binary operations
--------------------------

With binary operations between pandas data structures, there are two key points
of interest:

  * Broadcasting behavior between higher- (e.g. DataFrame) and
    lower-dimensional (e.g. Series) objects.
  * Missing data in computations

We will demonstrate how to manage these issues independently, though they can
be handled simultaneously.

Matching / broadcasting behavior
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

DataFrame has the methods **add, sub, mul, div** and related functions **radd,
rsub, ...** for carrying out binary operations. For broadcasting behavior,
Series input is of primary interest. Using these functions, you can use to
either match on the *index* or *columns* via the **axis** keyword:

.. ipython:: python
   :suppress:

   d = {'one' : Series(randn(3), index=['a', 'b', 'c']),
        'two' : Series(randn(4), index=['a', 'b', 'c', 'd']),
        'three' : Series(randn(3), index=['b', 'c', 'd'])}
   df = DataFrame(d)

.. ipython:: python

   df
   row = df.ix[1]
   column = df['two']

   df.sub(row, axis='columns')
   df.sub(row, axis=1)

   df.sub(column, axis='index')
   df.sub(column, axis=0)

With Panel, describing the matching behavior is a bit more difficult, so
the arithmetic methods instead (and perhaps confusingly?) give you the option
to specify the *broadcast axis*. For example, suppose we wished to demean the
data over a particular axis. This can be accomplished by taking the mean over
an axis and broadcasting over the same axis:

.. ipython:: python

   major_mean = wp.mean(axis='major')
   major_mean
   wp.sub(major_mean, axis='major')

And similarly for axis="items" and axis="minor".

.. note::

   I could be convinced to make the **axis** argument in the DataFrame methods
   match the broadcasting behavior of Panel. Though it would require a
   transition period so users can change their code...

Missing data / operations with fill values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Series and DataFrame (though not yet in Panel), the arithmetic functions
have the option of inputting a *fill_value*, namely a value to substitute when
at most one of the values at a location are missing. For example, when adding
two DataFrame objects, you may wish to treat NaN as 0 unless both DataFrames
are missing that value, in which case the result will be NaN (you can later
replace NaN with some other value using ``fillna`` if you wish).

.. ipython:: python
   :suppress:

   df2 = df.copy()
   df2['three']['a'] = 1.

.. ipython:: python

   df
   df2
   df + df2
   df.add(df2, fill_value=0)

Combining overlapping data sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A problem occasionally arising is the combination of two similar data sets
where values in one are preferred over the other. An example would be two data
series representing a particular economic indicator where one is considered to
be of "higher quality". However, the lower quality series might extend further
back in history or have more complete data coverage. As such, we would like to
combine two DataFrame objects where missing values in one DataFrame are
conditionally filled with like-labeled values from the other DataFrame. The
function implementing this operation is ``combine_first``, which we illustrate:

.. ipython:: python

   df1 = DataFrame({'A' : [1., np.nan, 3., 5., np.nan],
                    'B' : [np.nan, 2., 3., np.nan, 6.]})
   df2 = DataFrame({'A' : [5., 2., 4., np.nan, 3., 7.],
                    'B' : [np.nan, np.nan, 3., 4., 6., 8.]})
   df1
   df2
   df1.combine_first(df2)

General DataFrame Combine
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``combine_first`` method above calls the more general DataFrame method
``combine``. This method takes another DataFrame and a combiner function,
aligns the input DataFrame and then passes the combiner function pairs of
Series (ie, columns whose names are the same).

So, for instance, to reproduce ``combine_first`` as above:

.. ipython:: python

   combiner = lambda x, y: np.where(isnull(x), y, x)
   df1.combine(df2, combiner)

.. _basics.stats:

Descriptive statistics
----------------------

A large number of methods for computing descriptive statistics and other related
operations on :ref:`Series <api.series.stats>`, :ref:`DataFrame
<api.dataframe.stats>`, and :ref:`Panel <api.panel.stats>`. Most of these
are aggregations (hence producing a lower-dimensional result) like **sum**,
**mean**, and **quantile**, but some of them, like **cumsum** and **cumprod**,
produce an object of the same size. Generally speaking, these methods take an
**axis** argument, just like *ndarray.{sum, std, ...}*, but the axis can be
specified by name or integer:

  - **Series**: no axis argument needed
  - **DataFrame**: "index" (axis=0, default), "columns" (axis=1)
  - **Panel**: "items" (axis=0), "major" (axis=1, default), "minor"
    (axis=2)

For example:

.. ipython:: python

   df
   df.mean(0)
   df.mean(1)

All such methods have a ``skipna`` option signaling whether to exclude missing
data (``True`` by default):

.. ipython:: python

   df.sum(0, skipna=False)
   df.sum(axis=1, skipna=True)

Combined with the broadcasting / arithmetic behavior, one can describe various
statistical procedures, like standardization (rendering data zero mean and
standard deviation 1), very concisely:

.. ipython:: python

   ts_stand = (df - df.mean()) / df.std()
   ts_stand.std()
   xs_stand = df.sub(df.mean(1), axis=0).div(df.std(1), axis=0)
   xs_stand.std(1)

Note that methods like **cumsum** and **cumprod** preserve the location of NA
values:

.. ipython:: python

   df.cumsum()

Here is a quick reference summary table of common functions. Each also takes an
optional ``level`` parameter which applies only if the object has a
:ref:`hierarchical index<indexing.hierarchical>`.

.. csv-table::
    :header: "Function", "Description"
    :widths: 20, 80

    ``count``, Number of non-null observations
    ``sum``, Sum of values
    ``mean``, Mean of values
    ``mad``, Mean absolute deviation
    ``median``, Arithmetic median of values
    ``min``, Minimum
    ``max``, Maximum
    ``abs``, Absolute Value
    ``prod``, Product of values
    ``std``, Unbiased standard deviation
    ``var``, Unbiased variance
    ``skew``, Unbiased skewness (3rd moment)
    ``kurt``, Unbiased kurtosis (4th moment)
    ``quantile``, Sample quantile (value at %)
    ``cumsum``, Cumulative sum
    ``cumprod``, Cumulative product
    ``cummax``, Cumulative maximum
    ``cummin``, Cumulative minimum

Note that by chance some NumPy methods, like ``mean``, ``std``, and ``sum``,
will exclude NAs on Series input by default:

.. ipython:: python

   np.mean(df['one'])
   np.mean(df['one'].values)

``Series`` also has a method ``nunique`` which will return the number of unique
non-null values:

.. ipython:: python

   series = Series(randn(500))
   series[20:500] = np.nan
   series[10:20]  = 5
   series.nunique()


Summarizing data: describe
~~~~~~~~~~~~~~~~~~~~~~~~~~

There is a convenient ``describe`` function which computes a variety of summary
statistics about a Series or the columns of a DataFrame (excluding NAs of
course):

.. ipython:: python

    series = Series(randn(1000))
    series[::2] = np.nan
    series.describe()
    frame = DataFrame(randn(1000, 5), columns=['a', 'b', 'c', 'd', 'e'])
    frame.ix[::2] = np.nan
    frame.describe()

.. _basics.describe:

For a non-numerical Series object, `describe` will give a simple summary of the
number of unique values and most frequently occurring values:


.. ipython:: python

   s = Series(['a', 'a', 'b', 'b', 'a', 'a', np.nan, 'c', 'd', 'a'])
   s.describe()

There also is a utility function, ``value_range`` which takes a DataFrame and
returns a series with the minimum/maximum values in the DataFrame.

.. _basics.idxmin:

Index of Min/Max Values
~~~~~~~~~~~~~~~~~~~~~~~

The ``idxmin`` and ``idxmax`` functions on Series and DataFrame compute the
index labels with the minimum and maximum corresponding values:

.. ipython:: python

   s1 = Series(randn(5))
   s1
   s1.idxmin(), s1.idxmax()

   df1 = DataFrame(randn(5,3), columns=['A','B','C'])
   df1
   df1.idxmin(axis=0)
   df1.idxmax(axis=1)

.. _basics.apply:

Function application
--------------------

Arbitrary functions can be applied along the axes of a DataFrame or Panel
using the ``apply`` method, which, like the descriptive statistics methods,
take an optional ``axis`` argument:

.. ipython:: python

   df.apply(np.mean)
   df.apply(np.mean, axis=1)
   df.apply(lambda x: x.max() - x.min())
   df.apply(np.cumsum)
   df.apply(np.exp)

Depending on the return type of the function passed to ``apply``, the result
will either be of lower dimension or the same dimension.

``apply`` combined with some cleverness can be used to answer many questions
about a data set. For example, suppose we wanted to extract the date where the
maximum value for each column occurred:


.. ipython:: python

   tsdf = DataFrame(randn(1000, 3), columns=['A', 'B', 'C'],
                    index=DateRange('1/1/2000', periods=1000))
   tsdf.apply(lambda x: x.index[x.dropna().argmax()])

You may also pass additional arguments and keyword arguments to the ``apply``
method. For instance, consider the following function you would like to apply:

.. code-block:: python

   def subtract_and_divide(x, sub, divide=1):
       return (x - sub) / divide

You may then apply this function as follows:

.. code-block:: python

   df.apply(subtract_and_divide, args=(5,), divide=3)

Another useful feature is the ability to pass Series methods to carry out some
Series operation on each column or row:

.. ipython:: python
   :suppress:

   tsdf = DataFrame(randn(10, 3), columns=['A', 'B', 'C'],
                    index=DateRange('1/1/2000', periods=10))
   tsdf.values[3:7] = np.nan

.. ipython:: python

   tsdf
   tsdf.apply(Series.interpolate)

Finally, ``apply`` takes an argument ``raw`` which is False by default, which
converts each row or column into a Series before applying the function. When
set to True, the passed function will instead receive an ndarray object, which
has positive performance implications if you do not need the indexing
functionality.

.. seealso::

   The section on :ref:`GroupBy <groupby>` demonstrates related, flexible
   functionality for grouping by some criterion, applying, and combining the
   results into a Series, DataFrame, etc.

Applying elementwise Python functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since not all functions can be vectorized (accept NumPy arrays and return
another array or value), the methods ``applymap`` on DataFrame and analogously
``map`` on Series accept any Python function taking a single value and
returning a single value. For example:

.. ipython:: python

   f = lambda x: len(str(x))
   df['one'].map(f)
   df.applymap(f)

``Series.map`` has an additional feature which is that it can be used to easily
"link" or "map" values defined by a secondary series. This is closely related
to :ref:`merging/joining functionality <merging>`:


.. ipython:: python

   s = Series(['six', 'seven', 'six', 'seven', 'six'],
              index=['a', 'b', 'c', 'd', 'e'])
   t = Series({'six' : 6., 'seven' : 7.})
   s
   s.map(t)

.. _basics.reindexing:

Reindexing and altering labels
------------------------------

``reindex`` is the fundamental data alignment method in pandas. It is used to
implement nearly all other features relying on label-alignment
functionality. To *reindex* means to conform the data to match a given set of
labels along a particular axis. This accomplishes several things:

  * Reorders the existing data to match a new set of labels
  * Inserts missing value (NA) markers in label locations where no data for
    that label existed
  * If specified, **fill** data for missing labels using logic (highly relevant
    to working with time series data)

Here is a simple example:

.. ipython:: python

   s = Series(randn(5), index=['a', 'b', 'c', 'd', 'e'])
   s
   s.reindex(['e', 'b', 'f', 'd'])

Here, the ``f`` label was not contained in the Series and hence appears as
``NaN`` in the result.

With a DataFrame, you can simultaneously reindex the index and columns:

.. ipython:: python

   df
   df.reindex(index=['c', 'f', 'b'], columns=['three', 'two', 'one'])

For convenience, you may utilize the ``reindex_axis`` method, which takes the
labels and a keyword ``axis`` paramater.

Note that the ``Index`` objects containing the actual axis labels can be
**shared** between objects. So if we have a Series and a DataFrame, the
following can be done:

.. ipython:: python

   rs = s.reindex(df.index)
   rs
   rs.index is df.index

This means that the reindexed Series's index is the same Python object as the
DataFrame's index.


.. seealso::

   :ref:`Advanced indexing <indexing.advanced>` is an even more concise way of
   doing reindexing.

.. note::

    When writing performance-sensitive code, there is a good reason to spend
    some time becoming a reindexing ninja: **many operations are faster on
    pre-aligned data**. Adding two unaligned DataFrames internally triggers a
    reindexing step. For exploratory analysis you will hardly notice the
    difference (because ``reindex`` has been heavily optimized), but when CPU
    cycles matter sprinking a few explicit ``reindex`` calls here and there can
    have an impact.

.. _basics.reindex_like:

Reindexing to align with another object
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may wish to take an object and reindex its axes to be labeled the same as
another object. While the syntax for this is straightforward albeit verbose, it
is a common enough operation that the ``reindex_like`` method is available to
make this simpler:

.. ipython:: python
   :suppress:

   df2 = df.reindex(['a', 'b', 'c'], columns=['one', 'two'])
   df2 = df2 - df2.mean()


.. ipython:: python

   df
   df2
   df.reindex_like(df2)

Reindexing with ``reindex_axis``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _basics.align:

Aligning objects with each other with ``align``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``align`` method is the fastest way to simultaneously align two objects. It
supports a ``join`` argument (related to :ref:`joining and merging <merging>`):

  - ``join='outer'``: take the union of the indexes
  - ``join='left'``: use the calling object's index
  - ``join='right'``: use the passed object's index
  - ``join='inner'``: intersect the indexes

It returns a tuple with both of the reindexed Series:

.. ipython:: python

   s = Series(randn(5), index=['a', 'b', 'c', 'd', 'e'])
   s1 = s[:4]
   s2 = s[1:]
   s1.align(s2)
   s1.align(s2, join='inner')
   s1.align(s2, join='left')

.. _basics.df_join:

For DataFrames, the join method will be applied to both the index and the
columns by default:

.. ipython:: python

   df.align(df2, join='inner')

You can also pass an ``axis`` option to only align on the specified axis:

.. ipython:: python

   df.align(df2, join='inner', axis=0)

.. _basics.align.frame.series:

If you pass a Series to ``DataFrame.align``, you can choose to align both
objects either on the DataFrame's index or columns using the ``axis`` argument:

.. ipython:: python

   df.align(df2.ix[0], axis=1)

.. _basics.reindex_fill:

Filling while reindexing
~~~~~~~~~~~~~~~~~~~~~~~~

``reindex`` takes an optional parameter ``method`` which is a filling method
chosen from the following table:

.. csv-table::
    :header: "Method", "Action"
    :widths: 30, 50

    pad / ffill, Fill values forward
    bfill / backfill, Fill values backward

Other fill methods could be added, of course, but these are the two most
commonly used for time series data. In a way they only make sense for time
series or otherwise ordered data, but you may have an application on non-time
series data where this sort of "interpolation" logic is the correct thing to
do. More sophisticated interpolation of missing values would be an obvious
extension.

We illustrate these fill methods on a simple TimeSeries:

.. ipython:: python

   rng = DateRange('1/3/2000', periods=8)
   ts = Series(randn(8), index=rng)
   ts2 = ts[[0, 3, 6]]
   ts
   ts2

   ts2.reindex(ts.index)
   ts2.reindex(ts.index, method='ffill')
   ts2.reindex(ts.index, method='bfill')

Note the same result could have been achieved using :ref:`fillna
<missing_data.fillna>`:

.. ipython:: python

   ts2.reindex(ts.index).fillna(method='ffill')

Note these methods generally assume that the indexes are **sorted**. They may
be modified in the future to be a bit more flexible but as time series data is
ordered most of the time anyway, this has not been a major priority.

.. _basics.drop:

Dropping labels from an axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A method closely related to ``reindex`` is the ``drop`` function. It removes a
set of labels from an axis:

.. ipython:: python

   df
   df.drop(['a', 'd'], axis=0)
   df.drop(['one'], axis=1)

Note that the following also works, but a bit less obvious / clean:

.. ipython:: python

   df.reindex(df.index - ['a', 'd'])

.. _basics.rename:

Renaming / mapping labels
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``rename`` method allows you to relabel an axis based on some mapping (a
dict or Series) or an arbitrary function.

.. ipython:: python

   s
   s.rename(str.upper)

If you pass a function, it must return a value when called with any of the
labels (and must produce a set of unique values). But if you pass a dict or
Series, it need only contain a subset of the labels as keys:

.. ipython:: python

   df.rename(columns={'one' : 'foo', 'two' : 'bar'},
             index={'a' : 'apple', 'b' : 'banana', 'd' : 'durian'})

The ``rename`` method also provides a ``copy`` named parameter that is by
default ``True`` and copies the underlying data. Pass ``copy=False`` to rename
the data in place.

.. _basics.rename_axis:

The Panel class has an a related ``rename_axis`` class which can rename any of
its three axes.

Iteration
---------

Considering the pandas as somewhat dict-like structure, basic iteration
produces the "keys" of the objects, namely:

  * **Series**: the index label
  * **DataFrame**: the column labels
  * **Panel**: the item labels

Thus, for example:

.. ipython::

   In [0]: for col in df:
      ...:     print col
      ...:

iteritems
~~~~~~~~~

Consistent with the dict-like interface, **iteritems** iterates through
key-value pairs:

  * **Series**: (index, scalar value) pairs
  * **DataFrame**: (column, Series) pairs
  * **Panel**: (item, DataFrame) pairs

For example:

.. ipython::

   In [0]: for item, frame in wp.iteritems():
      ...:     print item
      ...:     print frame
      ...:


.. _basics.iterrows:

iterrows
~~~~~~~~

New in v0.7 is the ability to iterate efficiently through rows of a
DataFrame. It returns an iterator yielding each index value along with a Series
containing the data in each row:

.. ipython::

   In [0]: for row_index, row in df2.iterrows():
      ...:     print '%s\n%s' % (row_index, row)
      ...:


For instance, a contrived way to transpose the dataframe would be:

.. ipython:: python

   df2 = DataFrame({'x': [1, 2, 3], 'y': [4, 5, 6]})
   print df2
   print df2.T

   df2_t = DataFrame(dict((idx,values) for idx, values in df2.iterrows()))
   print df2_t

.. _basics.sorting:

Sorting by index and value
--------------------------

There are two obvious kinds of sorting that you may be interested in: sorting
by label and sorting by actual values. The primary method for sorting axis
labels (indexes) across data structures is the ``sort_index`` method.

.. ipython:: python

   unsorted_df = df.reindex(index=['a', 'd', 'c', 'b'],
                            columns=['three', 'two', 'one'])
   unsorted_df.sort_index()
   unsorted_df.sort_index(ascending=False)
   unsorted_df.sort_index(axis=1)

``DataFrame.sort_index`` can accept an optional ``by`` argument for ``axis=0``
which will use an arbitrary vector or a column name of the DataFrame to
determine the sort order:

.. ipython:: python

   df.sort_index(by='two')

The ``by`` argument can take a list of column names, e.g.:

.. ipython:: python

   df = DataFrame({'one':[2,1,1,1],'two':[1,3,2,4],'three':[5,4,3,2]})
   df[['one', 'two', 'three']].sort_index(by=['one','two'])

Series has the method ``order`` (analogous to `R's order function
<http://stat.ethz.ch/R-manual/R-patched/library/base/html/order.html>`__) which
sorts by value, with special treatment of NA values via the ``na_last``
argument:

.. ipython:: python

   s[2] = np.nan
   s.order()
   s.order(na_last=False)

Some other sorting notes / nuances:

  * ``Series.sort`` sorts a Series by value in-place. This is to provide
    compatibility with NumPy methods which expect the ``ndarray.sort``
    behavior.
  * ``DataFrame.sort`` takes a ``column`` argument instead of ``by``. This
    method will likely be deprecated in a future release in favor of just using
    ``sort_index``.

.. _basics.cast:

Copying, type casting
---------------------

The ``copy`` method on pandas objects copies the underlying data (though not
the axis indexes, since they are immutable) and returns a new object. Note that
**it is seldom necessary to copy objects**. For example, there are only a
handful of ways to alter a DataFrame *in-place*:

  * Inserting, deleting, or modifying a column
  * Assigning to the ``index`` or ``columns`` attributes
  * For homogeneous data, directly modifying the values via the ``values``
    attribute or advanced indexing

To be clear, no pandas methods have the side effect of modifying your data;
almost all methods return new objects, leaving the original object
untouched. If data is modified, it is because you did so explicitly.

Data can be explicitly cast to a NumPy dtype by using the ``astype`` method or
alternately passing the ``dtype`` keyword argument to the object constructor.

.. ipython:: python

   df = DataFrame(np.arange(12).reshape((4, 3)))
   df[0].dtype
   df.astype(float)[0].dtype
   df = DataFrame(np.arange(12).reshape((4, 3)), dtype=float)
   df[0].dtype

.. _basics.cast.infer:

Inferring better types for object columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``convert_objects`` DataFrame method will attempt to convert
``dtype=object`` columns to a better NumPy dtype. Occasionally (after
transposing multiple times, for example), a mixed-type DataFrame will end up
with everything as ``dtype=object``. This method attempts to fix that:

.. ipython:: python

   df = DataFrame(randn(6, 3), columns=['a', 'b', 'c'])
   df['d'] = 'foo'
   df
   df = df.T.T
   df.dtypes
   converted = df.convert_objects()
   converted.dtypes

.. _basics.serialize:

Pickling and serialization
--------------------------

All pandas objects are equipped with ``save`` methods which use Python's
``cPickle`` module to save data structures to disk using the pickle format.

.. ipython:: python

   df
   df.save('foo.pickle')

The ``load`` function in the ``pandas`` namespace can be used to load any
pickled pandas object (or any other pickled object) from file:


.. ipython:: python

   load('foo.pickle')

There is also a ``save`` function which takes any object as its first argument:

.. ipython:: python

   save(df, 'foo.pickle')
   load('foo.pickle')

.. ipython:: python
   :suppress:

   import os
   os.remove('foo.pickle')

Console Output Formatting
-------------------------

.. _basics.console_output:

Use the ``set_eng_float_format`` function in the ``pandas.core.common`` module
to alter the floating-point formatting of pandas objects to produce a particular
format.

For instance:

.. ipython:: python

   set_eng_float_format(accuracy=3, use_eng_prefix=True)
   df['a']/1.e3
   df['a']/1.e6

.. ipython:: python
   :suppress:

   reset_printoptions()


The ``set_printoptions`` function has a number of options for controlling how
floating point numbers are formatted (using hte ``precision`` argument) in the
console and . The ``max_rows`` and ``max_columns`` control how many rows and
columns of DataFrame objects are shown by default. If ``max_columns`` is set to
0 (the default, in fact), the library will attempt to fit the DataFrame's
string representation into the current terminal width, and defaulting to the
summary view otherwise.
