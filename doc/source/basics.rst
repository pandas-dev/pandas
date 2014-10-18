.. currentmodule:: pandas
.. _basics:

.. ipython:: python
   :suppress:

   import numpy as np
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   from pandas.compat import lrange
   options.display.max_rows=15

==============================
 Essential Basic Functionality
==============================

Here we discuss a lot of the essential functionality common to the pandas data
structures. Here's how to create some of the objects used in the examples from
the previous section:

.. ipython:: python

   index = date_range('1/1/2000', periods=8)
   s = Series(randn(5), index=['a', 'b', 'c', 'd', 'e'])
   df = DataFrame(randn(8, 3), index=index,
                  columns=['A', 'B', 'C'])
   wp = Panel(randn(2, 5, 4), items=['Item1', 'Item2'],
              major_axis=date_range('1/1/2000', periods=5),
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

.. _basics.accelerate:

Accelerated operations
----------------------

pandas has support for accelerating certain types of binary numerical and boolean operations using
the ``numexpr`` library (starting in 0.11.0) and the ``bottleneck`` libraries.

These libraries are especially useful when dealing with large data sets, and provide large
speedups. ``numexpr`` uses smart chunking, caching, and multiple cores. ``bottleneck`` is
a set of specialized cython routines that are especially fast when dealing with arrays that have
``nans``.

Here is a sample (using 100 column x 100,000 row ``DataFrames``):

.. csv-table::
    :header: "Operation", "0.11.0 (ms)", "Prior Version (ms)", "Ratio to Prior"
    :widths: 25, 25, 25, 25
    :delim: ;

    ``df1 > df2``; 13.32; 125.35;  0.1063
    ``df1 * df2``; 21.71;  36.63;  0.5928
    ``df1 + df2``; 22.04;  36.50;  0.6039

You are highly encouraged to install both libraries. See the section
:ref:`Recommended Dependencies <install.recommended_dependencies>` for more installation info.

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

   df = DataFrame({'one' : Series(randn(3), index=['a', 'b', 'c']),
                   'two' : Series(randn(4), index=['a', 'b', 'c', 'd']),
                   'three' : Series(randn(3), index=['b', 'c', 'd'])})
   df
   row = df.ix[1]
   column = df['two']

   df.sub(row, axis='columns')
   df.sub(row, axis=1)

   df.sub(column, axis='index')
   df.sub(column, axis=0)

.. ipython:: python
   :suppress:

   df_orig = df

Furthermore you can align a level of a multi-indexed DataFrame with a Series.

.. ipython:: python

   dfmi = df.copy()
   dfmi.index = MultiIndex.from_tuples([(1,'a'),(1,'b'),(1,'c'),(2,'a')],
                                       names=['first','second'])
   dfmi.sub(column, axis=0, level='second')

With Panel, describing the matching behavior is a bit more difficult, so
the arithmetic methods instead (and perhaps confusingly?) give you the option
to specify the *broadcast axis*. For example, suppose we wished to demean the
data over a particular axis. This can be accomplished by taking the mean over
an axis and broadcasting over the same axis:

.. ipython:: python

   major_mean = wp.mean(axis='major')
   major_mean
   wp.sub(major_mean, axis='major')

And similarly for ``axis="items"`` and ``axis="minor"``.

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

.. _basics.compare:

Flexible Comparisons
~~~~~~~~~~~~~~~~~~~~

Starting in v0.8, pandas introduced binary comparison methods eq, ne, lt, gt,
le, and ge to Series and DataFrame whose behavior is analogous to the binary
arithmetic operations described above:

.. ipython:: python

   df.gt(df2)
   df2.ne(df)

These operations produce a pandas object the same type as the left-hand-side input
that if of dtype ``bool``. These ``boolean`` objects can be used in indexing operations,
see :ref:`here<indexing.boolean>`

.. _basics.reductions:

Boolean Reductions
~~~~~~~~~~~~~~~~~~

You can apply the reductions: ``empty``, ``any()``, ``all()``, and ``bool()`` to provide a
way to summarize a boolean result.

.. ipython:: python

   (df>0).all()
   (df>0).any()

You can reduce to a final boolean value.

.. ipython:: python

   (df>0).any().any()

You can test if a pandas object is empty, via the ``empty`` property.

.. ipython:: python

   df.empty
   DataFrame(columns=list('ABC')).empty

To evaluate single-element pandas objects in a boolean context, use the method ``.bool()``:

.. ipython:: python

   Series([True]).bool()
   Series([False]).bool()
   DataFrame([[True]]).bool()
   DataFrame([[False]]).bool()

.. warning::

   You might be tempted to do the following:

   .. code-block:: python

       >>>if df:
            ...

   Or

   .. code-block:: python

       >>> df and df2

   These both will raise as you are trying to compare multiple values.

   .. code-block:: python

       ValueError: The truth value of an array is ambiguous. Use a.empty, a.any() or a.all().

See :ref:`gotchas<gotchas.truth>` for a more detailed discussion.

.. _basics.equals:

Comparing if objects are equivalent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often you may find there is more than one way to compute the same
result.  As a simple example, consider ``df+df`` and ``df*2``. To test
that these two computations produce the same result, given the tools
shown above, you might imagine using ``(df+df == df*2).all()``. But in
fact, this expression is False:

.. ipython:: python

   df+df == df*2
   (df+df == df*2).all()

Notice that the boolean DataFrame ``df+df == df*2`` contains some False values!
That is because NaNs do not compare as equals:

.. ipython:: python

   np.nan == np.nan

So, as of v0.13.1, NDFrames (such as Series, DataFrames, and Panels)
have an ``equals`` method for testing equality, with NaNs in corresponding
locations treated as equal.

.. ipython:: python

   (df+df).equals(df*2)

Note that the Series or DataFrame index needs to be in the same order for
equality to be True:

.. ipython:: python

   df1 = DataFrame({'col':['foo', 0, np.nan]})
   df2 = DataFrame({'col':[np.nan, 0, 'foo']}, index=[2,1,0])
   df1.equals(df2)
   df1.equals(df2.sort())


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
Series (i.e., columns whose names are the same).

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
:ref:`hierarchical index<advanced.hierarchical>`.

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
    ``mode``, Mode
    ``abs``, Absolute Value
    ``prod``, Product of values
    ``std``, Unbiased standard deviation
    ``var``, Unbiased variance
    ``sem``, Unbiased standard error of the mean
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

.. _basics.describe:

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

You can select specific percentiles to include in the output:

.. ipython:: python

    series.describe(percentiles=[.05, .25, .75, .95])

By default, the median is always included.

For a non-numerical Series object, `describe` will give a simple summary of the
number of unique values and most frequently occurring values:


.. ipython:: python

   s = Series(['a', 'a', 'b', 'b', 'a', 'a', np.nan, 'c', 'd', 'a'])
   s.describe()

Note that on a mixed-type DataFrame object, `describe` will restrict the summary to
include only numerical columns or, if none are, only categorical columns:

.. ipython:: python

    frame = DataFrame({'a': ['Yes', 'Yes', 'No', 'No'], 'b': range(4)})
    frame.describe()

This behaviour can be controlled by providing a list of types as ``include``/``exclude``
arguments. The special value ``all`` can also be used:

.. ipython:: python

    frame.describe(include=['object'])
    frame.describe(include=['number'])
    frame.describe(include='all')

That feature relies on :ref:`select_dtypes <basics.selectdtypes>`. Refer to there for details about accepted inputs.

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

When there are multiple rows (or columns) matching the minimum or maximum
value, ``idxmin`` and ``idxmax`` return the first matching index:

.. ipython:: python

   df3 = DataFrame([2, 1, 1, 3, np.nan], columns=['A'], index=list('edcba'))
   df3
   df3['A'].idxmin()

.. note::

   ``idxmin`` and ``idxmax`` are called ``argmin`` and ``argmax`` in NumPy.

.. _basics.discretization:

Value counts (histogramming) / Mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``value_counts`` Series method and top-level function computes a histogram
of a 1D array of values. It can also be used as a function on regular arrays:

.. ipython:: python

   data = np.random.randint(0, 7, size=50)
   data
   s = Series(data)
   s.value_counts()
   value_counts(data)

Similarly, you can get the most frequently occurring value(s) (the mode) of the values in a Series or DataFrame:

.. ipython:: python

    s5 = Series([1, 1, 3, 3, 3, 5, 5, 7, 7, 7])
    s5.mode()
    df5 = DataFrame({"A": np.random.randint(0, 7, size=50),
                     "B": np.random.randint(-10, 15, size=50)})
    df5.mode()


Discretization and quantiling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Continuous values can be discretized using the ``cut`` (bins based on values)
and ``qcut`` (bins based on sample quantiles) functions:

.. ipython:: python

   arr = np.random.randn(20)
   factor = cut(arr, 4)
   factor

   factor = cut(arr, [-5, -1, 0, 1, 5])
   factor

``qcut`` computes sample quantiles. For example, we could slice up some
normally distributed data into equal-size quartiles like so:

.. ipython:: python

   arr = np.random.randn(30)
   factor = qcut(arr, [0, .25, .5, .75, 1])
   factor
   value_counts(factor)

We can also pass infinite values to define the bins:

.. ipython:: python

   arr = np.random.randn(20)
   factor = cut(arr, [-np.inf, 0, np.inf])
   factor

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
                    index=date_range('1/1/2000', periods=1000))
   tsdf.apply(lambda x: x.idxmax())

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
                    index=date_range('1/1/2000', periods=10))
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
   :suppress:

   df4 = df_orig.copy()

.. ipython:: python

   df4
   f = lambda x: len(str(x))
   df4['one'].map(f)
   df4.applymap(f)

``Series.map`` has an additional feature which is that it can be used to easily
"link" or "map" values defined by a secondary series. This is closely related
to :ref:`merging/joining functionality <merging>`:


.. ipython:: python

   s = Series(['six', 'seven', 'six', 'seven', 'six'],
              index=['a', 'b', 'c', 'd', 'e'])
   t = Series({'six' : 6., 'seven' : 7.})
   s
   s.map(t)


.. _basics.apply_panel:

Applying with a Panel
~~~~~~~~~~~~~~~~~~~~~

Applying with a ``Panel`` will pass a ``Series`` to the applied function. If the applied
function returns a ``Series``, the result of the application will be a ``Panel``. If the applied function
reduces to a scalar, the result of the application will be a ``DataFrame``.

.. note::

   Prior to 0.13.1 ``apply`` on a ``Panel`` would only work on ``ufuncs`` (e.g. ``np.sum/np.max``).

.. ipython:: python

   import pandas.util.testing as tm
   panel = tm.makePanel(5)
   panel
   panel['ItemA']

A transformational apply.

.. ipython:: python

   result = panel.apply(lambda x: x*2, axis='items')
   result
   result['ItemA']

A reduction operation.

.. ipython:: python

   panel.apply(lambda x: x.dtype, axis='items')

A similar reduction type operation

.. ipython:: python

   panel.apply(lambda x: x.sum(), axis='major_axis')

This last reduction is equivalent to

.. ipython:: python

   panel.sum('major_axis')

A transformation operation that returns a ``Panel``, but is computing
the z-score across the ``major_axis``.

.. ipython:: python

   result = panel.apply(
              lambda x: (x-x.mean())/x.std(),
              axis='major_axis')
   result
   result['ItemA']

Apply can also accept multiple axes in the ``axis`` argument. This will pass a
``DataFrame`` of the cross-section to the applied function.

.. ipython:: python

   f = lambda x: ((x.T-x.mean(1))/x.std(1)).T

   result = panel.apply(f, axis = ['items','major_axis'])
   result
   result.loc[:,:,'ItemA']

This is equivalent to the following

.. ipython:: python

   result = Panel(dict([ (ax,f(panel.loc[:,:,ax]))
                           for ax in panel.minor_axis ]))
   result
   result.loc[:,:,'ItemA']

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
labels and a keyword ``axis`` parameter.

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

   :ref:`MultiIndex / Advanced Indexing <advanced>` is an even more concise way of
   doing reindexing.

.. note::

    When writing performance-sensitive code, there is a good reason to spend
    some time becoming a reindexing ninja: **many operations are faster on
    pre-aligned data**. Adding two unaligned DataFrames internally triggers a
    reindexing step. For exploratory analysis you will hardly notice the
    difference (because ``reindex`` has been heavily optimized), but when CPU
    cycles matter sprinkling a few explicit ``reindex`` calls here and there can
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
   df3 = df2 - df2.mean()


.. ipython:: python

   df2
   df3
   df.reindex_like(df2)

Reindexing with ``reindex_axis``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _basics.align:

Aligning objects with each other with ``align``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``align`` method is the fastest way to simultaneously align two objects. It
supports a ``join`` argument (related to :ref:`joining and merging <merging>`):

  - ``join='outer'``: take the union of the indexes (default)
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

   rng = date_range('1/3/2000', periods=8)
   ts = Series(randn(8), index=rng)
   ts2 = ts[[0, 3, 6]]
   ts
   ts2

   ts2.reindex(ts.index)
   ts2.reindex(ts.index, method='ffill')
   ts2.reindex(ts.index, method='bfill')

Note these methods require that the indexes are **order increasing**.

Note the same result could have been achieved using :ref:`fillna
<missing_data.fillna>`:

.. ipython:: python

   ts2.reindex(ts.index).fillna(method='ffill')

Note that ``reindex`` will raise a ValueError if the index is not
monotonic. ``fillna`` will not make any checks on the order of the index.

.. _basics.drop:

Dropping labels from an axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A method closely related to ``reindex`` is the ``drop`` function. It removes a
set of labels from an axis:

.. ipython:: python

   df
   df.drop(['a', 'd'], axis=0)
   df.drop(['one'], axis=1)

Note that the following also works, but is a bit less obvious / clean:

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

The ``rename`` method also provides an ``inplace`` named parameter that is by
default ``False`` and copies the underlying data. Pass ``inplace=True`` to
rename the data in place.

.. _basics.rename_axis:

The Panel class has a related ``rename_axis`` class which can rename any of
its three axes.

Iteration
---------

Because Series is array-like, basic iteration produces the values. Other data
structures follow the dict-like convention of iterating over the "keys" of the
objects. In short:

  * **Series**: values
  * **DataFrame**: column labels
  * **Panel**: item labels

Thus, for example:

.. ipython::

   In [0]: for col in df:
      ...:     print(col)
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
      ...:     print(item)
      ...:     print(frame)
      ...:


.. _basics.iterrows:

iterrows
~~~~~~~~

New in v0.7 is the ability to iterate efficiently through rows of a
DataFrame. It returns an iterator yielding each index value along with a Series
containing the data in each row:

.. ipython::

   In [0]: for row_index, row in df2.iterrows():
      ...:     print('%s\n%s' % (row_index, row))
      ...:

For instance, a contrived way to transpose the DataFrame would be:

.. ipython:: python

   df2 = DataFrame({'x': [1, 2, 3], 'y': [4, 5, 6]})
   print(df2)
   print(df2.T)

   df2_t = DataFrame(dict((idx,values) for idx, values in df2.iterrows()))
   print(df2_t)

.. note::

   ``iterrows`` does **not** preserve dtypes across the rows (dtypes are
   preserved across columns for DataFrames). For example,

    .. ipython:: python

      df_iter = DataFrame([[1, 1.0]], columns=['x', 'y'])
      row = next(df_iter.iterrows())[1]
      print(row['x'].dtype)
      print(df_iter['x'].dtype)

itertuples
~~~~~~~~~~

This method will return an iterator yielding a tuple for each row in the
DataFrame. The first element of the tuple will be the row's corresponding index
value, while the remaining values are the row values proper.

For instance,

.. ipython:: python

   for r in df2.itertuples():
       print(r)

.. _basics.dt_accessors:

.dt accessor
~~~~~~~~~~~~

``Series`` has an accessor to succinctly return datetime like properties for the *values* of the Series, if its a datetime/period like Series.
This will return a Series, indexed like the existing Series.

.. ipython:: python

   # datetime
   s = Series(date_range('20130101 09:10:12',periods=4))
   s
   s.dt.hour
   s.dt.second
   s.dt.day

This enables nice expressions like this:

.. ipython:: python

   s[s.dt.day==2]

You can easily produces tz aware transformations:

.. ipython:: python

   stz = s.dt.tz_localize('US/Eastern')
   stz
   stz.dt.tz

You can also chain these types of operations:

.. ipython:: python

   s.dt.tz_localize('UTC').dt.tz_convert('US/Eastern')

The ``.dt`` accessor works for period and timedelta dtypes.

.. ipython:: python

   # period
   s = Series(period_range('20130101',periods=4,freq='D'))
   s
   s.dt.year
   s.dt.day

.. ipython:: python

   # timedelta
   s = Series(timedelta_range('1 day 00:00:05',periods=4,freq='s'))
   s
   s.dt.days
   s.dt.seconds
   s.dt.components

.. note::

   ``Series.dt`` will raise a ``TypeError`` if you access with a non-datetimelike values

Vectorized string methods
-------------------------

Series is equipped with a set of string processing methods that make it easy to
operate on each element of the array. Perhaps most importantly, these methods
exclude missing/NA values automatically. These are accessed via the Series's
``str`` attribute and generally have names matching the equivalent (scalar)
built-in string methods. For example:

 .. ipython:: python

  s = Series(['A', 'B', 'C', 'Aaba', 'Baca', np.nan, 'CABA', 'dog', 'cat'])
  s.str.lower()

Powerful pattern-matching methods are provided as well, but note that
pattern-matching generally uses `regular expressions
<https://docs.python.org/2/library/re.html>`__ by default (and in some cases
always uses them).

Please see :ref:`Vectorized String Methods <text.string_methods>` for a complete
description.

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

   df1 = DataFrame({'one':[2,1,1,1],'two':[1,3,2,4],'three':[5,4,3,2]})
   df1.sort_index(by='two')

The ``by`` argument can take a list of column names, e.g.:

.. ipython:: python

   df1[['one', 'two', 'three']].sort_index(by=['one','two'])

Series has the method ``order`` (analogous to `R's order function
<http://stat.ethz.ch/R-manual/R-patched/library/base/html/order.html>`__) which
sorts by value, with special treatment of NA values via the ``na_position``
argument:

.. ipython:: python

   s[2] = np.nan
   s.order()
   s.order(na_position='first')

.. note::

   ``Series.sort`` sorts a Series by value in-place. This is to provide
   compatibility with NumPy methods which expect the ``ndarray.sort``
   behavior. ``Series.order`` returns a copy of the sorted data.

Series has the ``searchsorted`` method, which works similar to
``np.ndarray.searchsorted``.

.. ipython:: python

   ser = Series([1, 2, 3])
   ser.searchsorted([0, 3])
   ser.searchsorted([0, 4])
   ser.searchsorted([1, 3], side='right')
   ser.searchsorted([1, 3], side='left')
   ser = Series([3, 1, 2])
   ser.searchsorted([0, 3], sorter=np.argsort(ser))

.. _basics.nsorted:

smallest / largest values
~~~~~~~~~~~~~~~~~~~~~~~~~

.. versionadded:: 0.14.0

``Series`` has the ``nsmallest`` and ``nlargest`` methods which return the
smallest or largest :math:`n` values. For a large ``Series`` this can be much
faster than sorting the entire Series and calling ``head(n)`` on the result.

.. ipython:: python

   s = Series(np.random.permutation(10))
   s
   s.order()
   s.nsmallest(3)
   s.nlargest(3)


.. _basics.multi-index_sorting:

Sorting by a multi-index column
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You must be explicit about sorting when the column is a multi-index, and fully specify
all levels to ``by``.

.. ipython:: python

   df1.columns = MultiIndex.from_tuples([('a','one'),('a','two'),('b','three')])
   df1.sort_index(by=('a','two'))


Copying
-------

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

.. _basics.dtypes:

dtypes
------

The main types stored in pandas objects are ``float``, ``int``, ``bool``, ``datetime64[ns]``, ``timedelta[ns]``,
and ``object``. In addition these dtypes have item sizes, e.g. ``int64`` and ``int32``. A convenient ``dtypes``
attribute for DataFrames returns a Series with the data type of each column.

.. ipython:: python

   dft = DataFrame(dict( A = np.random.rand(3),
                         B = 1,
                         C = 'foo',
                         D = Timestamp('20010102'),
                         E = Series([1.0]*3).astype('float32'),
			 F = False,
			 G = Series([1]*3,dtype='int8')))
   dft
   dft.dtypes

On a ``Series`` use the ``dtype`` method.

.. ipython:: python

   dft['A'].dtype

If a pandas object contains data multiple dtypes *IN A SINGLE COLUMN*, the dtype of the
column will be chosen to accommodate all of the data types (``object`` is the most
general).

.. ipython:: python

   # these ints are coerced to floats
   Series([1, 2, 3, 4, 5, 6.])

   # string data forces an ``object`` dtype
   Series([1, 2, 3, 6., 'foo'])

The method ``get_dtype_counts`` will return the number of columns of
each type in a ``DataFrame``:

.. ipython:: python

   dft.get_dtype_counts()

Numeric dtypes will propagate and can coexist in DataFrames (starting in v0.11.0).
If a dtype is passed (either directly via the ``dtype`` keyword, a passed ``ndarray``,
or a passed ``Series``, then it will be preserved in DataFrame operations. Furthermore,
different numeric dtypes will **NOT** be combined. The following example will give you a taste.

.. ipython:: python

   df1 = DataFrame(randn(8, 1), columns = ['A'], dtype = 'float32')
   df1
   df1.dtypes
   df2 = DataFrame(dict( A = Series(randn(8),dtype='float16'),
                         B = Series(randn(8)),
                         C = Series(np.array(randn(8),dtype='uint8')) ))
   df2
   df2.dtypes

defaults
~~~~~~~~

By default integer types are ``int64`` and float types are ``float64``,
*REGARDLESS* of platform (32-bit or 64-bit). The following will all result in ``int64`` dtypes.

.. ipython:: python

   DataFrame([1, 2], columns=['a']).dtypes
   DataFrame({'a': [1, 2]}).dtypes
   DataFrame({'a': 1 }, index=list(range(2))).dtypes

Numpy, however will choose *platform-dependent* types when creating arrays.
The following **WILL** result in ``int32`` on 32-bit platform.

.. ipython:: python

   frame = DataFrame(np.array([1, 2]))


upcasting
~~~~~~~~~

Types can potentially be *upcasted* when combined with other types, meaning they are promoted
from the current type (say ``int`` to ``float``)

.. ipython:: python

   df3 = df1.reindex_like(df2).fillna(value=0.0) + df2
   df3
   df3.dtypes

The ``values`` attribute on a DataFrame return the *lower-common-denominator* of the dtypes, meaning
the dtype that can accommodate **ALL** of the types in the resulting homogeneous dtyped numpy array. This can
force some *upcasting*.

.. ipython:: python

   df3.values.dtype

astype
~~~~~~

.. _basics.cast:

You can use the ``astype`` method to explicitly convert dtypes from one to another. These will by default return a copy,
even if the dtype was unchanged (pass ``copy=False`` to change this behavior). In addition, they will raise an
exception if the astype operation is invalid.

Upcasting is always according to the **numpy** rules. If two different dtypes are involved in an operation,
then the more *general* one will be used as the result of the operation.

.. ipython:: python

   df3
   df3.dtypes

   # conversion of dtypes
   df3.astype('float32').dtypes

object conversion
~~~~~~~~~~~~~~~~~

``convert_objects`` is a method to try to force conversion of types from the ``object`` dtype to other types.
To force conversion of specific types that are *number like*, e.g. could be a string that represents a number,
pass ``convert_numeric=True``. This will force strings and numbers alike to be numbers if possible, otherwise
they will be set to ``np.nan``.

.. ipython:: python

   df3['D'] = '1.'
   df3['E'] = '1'
   df3.convert_objects(convert_numeric=True).dtypes

   # same, but specific dtype conversion
   df3['D'] = df3['D'].astype('float16')
   df3['E'] = df3['E'].astype('int32')
   df3.dtypes

To force conversion to ``datetime64[ns]``, pass ``convert_dates='coerce'``.
This will convert any datetime-like object to dates, forcing other values to ``NaT``.
This might be useful if you are reading in data which is mostly dates,
but occasionally has non-dates intermixed and you want to represent as missing.

.. ipython:: python

   s = Series([datetime(2001,1,1,0,0),
              'foo', 1.0, 1, Timestamp('20010104'),
              '20010105'],dtype='O')
   s
   s.convert_objects(convert_dates='coerce')

In addition, ``convert_objects`` will attempt the *soft* conversion of any *object* dtypes, meaning that if all
the objects in a Series are of the same type, the Series will have that dtype.

gotchas
~~~~~~~

Performing selection operations on ``integer`` type data can easily upcast the data to ``floating``.
The dtype of the input data will be preserved in cases where ``nans`` are not introduced (starting in 0.11.0)
See also :ref:`integer na gotchas <gotchas.intna>`

.. ipython:: python

   dfi = df3.astype('int32')
   dfi['E'] = 1
   dfi
   dfi.dtypes

   casted = dfi[dfi>0]
   casted
   casted.dtypes

While float dtypes are unchanged.

.. ipython:: python

   dfa = df3.copy()
   dfa['A'] = dfa['A'].astype('float32')
   dfa.dtypes

   casted = dfa[df2>0]
   casted
   casted.dtypes

Selecting columns based on ``dtype``
------------------------------------

.. _basics.selectdtypes:

.. versionadded:: 0.14.1

The :meth:`~pandas.DataFrame.select_dtypes` method implements subsetting of columns
based on their ``dtype``.

First, let's create a :class:`~pandas.DataFrame` with a slew of different
dtypes:

.. ipython:: python

   df = DataFrame({'string': list('abc'),
                   'int64': list(range(1, 4)),
                   'uint8': np.arange(3, 6).astype('u1'),
                   'float64': np.arange(4.0, 7.0),
                   'bool1': [True, False, True],
                   'bool2': [False, True, False],
                   'dates': pd.date_range('now', periods=3).values,
                   'category': pd.Categorical(list("ABC"))})
   df['tdeltas'] = df.dates.diff()
   df['uint64'] = np.arange(3, 6).astype('u8')
   df['other_dates'] = pd.date_range('20130101', periods=3).values
   df


``select_dtypes`` has two parameters ``include`` and ``exclude`` that allow you to
say "give me the columns WITH these dtypes" (``include``) and/or "give the
columns WITHOUT these dtypes" (``exclude``).

For example, to select ``bool`` columns

.. ipython:: python

   df.select_dtypes(include=[bool])

You can also pass the name of a dtype in the `numpy dtype hierarchy
<http://docs.scipy.org/doc/numpy/reference/arrays.scalars.html>`__:

.. ipython:: python

   df.select_dtypes(include=['bool'])

:meth:`~pandas.DataFrame.select_dtypes` also works with generic dtypes as well.

For example, to select all numeric and boolean columns while excluding unsigned
integers

.. ipython:: python

   df.select_dtypes(include=['number', 'bool'], exclude=['unsignedinteger'])

To select string columns you must use the ``object`` dtype:

.. ipython:: python

   df.select_dtypes(include=['object'])

To see all the child dtypes of a generic ``dtype`` like ``numpy.number`` you
can define a function that returns a tree of child dtypes:

.. ipython:: python

   def subdtypes(dtype):
       subs = dtype.__subclasses__()
       if not subs:
           return dtype
       return [dtype, [subdtypes(dt) for dt in subs]]

All numpy dtypes are subclasses of ``numpy.generic``:

.. ipython:: python

    subdtypes(np.generic)

.. note::

    Pandas also defines an additional ``category`` dtype, which is not integrated into the normal
    numpy hierarchy and wont show up with the above function.

.. note::

   The ``include`` and ``exclude`` parameters must be non-string sequences.
