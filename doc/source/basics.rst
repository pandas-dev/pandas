.. currentmodule:: pandas
.. _basics:

*********************
Data structure basics
*********************

We'll start with a quick, non-comprehensive overview of the fundamental data
structures in pandas to get you started. The fundamental behavior about data
types, indexing, and axis labeling / alignment apply across all of the
objects. To get started, import numpy and load pandas into your namespace:

.. ipython:: python
   :suppress:

   import numpy as np
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)

.. ipython:: python

   import numpy as np
   # will use a lot in examples
   randn = np.random.randn
   from pandas import *

Here is a basic tenet to keep in mind: **data alignment is intrinsic**. Link
between labels and data will not be broken unless done so explicitly by you.

We'll give a brief intro to the data structures, then consider all of the broad
categories of functionality and methods in separate sections.

.. _basics.series:

Series
------

:class:`Series` is a one-dimensional labeled array (technically a subclass of
ndarray) capable of holding any data type (integers, strings, floating point
numbers, Python objects, etc.). The axis labels are collectively referred to as
the **index**. The basic method to create a Series is to call:

::

    >>> s = Series(data, index=index)

Here, ``data`` can be many different things:

 - a Python dict
 - an ndarray
 - a scalar value (like 5)

The passed **index** is a list of axis labels. Thus, this separates into a few
cases depending on what **data is**:

**From ndarray**

If ``data`` is an ndarray, **index** must be the same length as **data**. If no
index is passed, one will be created having values ``[0, ..., len(data) - 1]``.

.. ipython:: python

   s = Series(randn(5), index=['a', 'b', 'c', 'd', 'e'])
   s
   s.index

   Series(randn(5))

.. note::

    The values in the index must be unique. If they are not, an exception will
    **not** be raised immediately, but attempting any operation involving the
    index will later result in an exception. In other words, the Index object
    containing the labels "lazily" checks whether the values are unique. The
    reason for being lazy is nearly all performance-based (there are many
    instances in computations, like parts of GroupBy, where the index is not
    used).

**From dict**

If ``data`` is a dict, if **index** is passed the values in data corresponding
to the labels in the index will be pulled out. Otherwise, an index will be
constructed from the sorted keys of the dict, if possible.

.. ipython:: python

   d = {'a' : 0., 'b' : 1., 'c' : 2.}
   Series(d)
   Series(d, index=['b', 'c', 'd', 'a'])

.. note::

    NaN (not a number) is the standard missing data marker used in pandas

**From scalar value** If ``data`` is a scalar value, an index must be
provided. The value will be repeated to match the length of **index**

.. ipython:: python

   Series(5., index=['a', 'b', 'c', 'd', 'e'])

Series is ndarray-like
~~~~~~~~~~~~~~~~~~~~~~

As a subclass of ndarray, Series is a valid argument to most NumPy functions
and behaves similarly to a NumPy array. However, things like slicing also slice
the index.

.. ipython :: python

    s[0]
    s[:3]
    s[s > s.median()]
    s[[4, 3, 1]]
    np.exp(s)

Series is dict-like
~~~~~~~~~~~~~~~~~~~

A Series is alike a fixed-size dict in that you can get and set values by index
label:

.. ipython :: python

    s['a']
    s['e'] = 12.
    s
    'e' in s
    'f' in s

If a label is not contained, an exception

.. code-block:: python

    >>> s['f']
    KeyError: 'f'

    >>> s.get('f')
    nan

Vectorized operations and label alignment with Series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When doing data analysis, as with raw NumPy arrays looping through Series
value-by-value is usually not necessary. Series can be also be passed into most
NumPy methods expecting an ndarray.


.. ipython:: python

    s + s
    s * 2
    np.exp(s)

A key difference between Series and ndarray is that operations between Series
automatically align the data based on label. Thus, you can write computations
without giving consideration to whether the Series involved have the same
labels.

.. ipython:: python

    s[1:] + s[:-1]

The result of an operation between unaligned Series will have the **union** of
the indexes involved. If a label is not found in one Series or the other, the
result will be marked as missing (NaN). Being able to write code without doing
any explicit data alignment grants immense freedom and flexibility in
interactive data analysis and research. The integrated data alignment features
of the pandas data structures set pandas apart from the majority of related
tools for working with labeled data.

.. note::

    In general, we chose to make the default result of operations between
    differently indexed objects yield the **union** of the indexes in order to
    avoid loss of information. Having an index label, though the data is
    missing, is typically important information as part of a computation. You
    of course have the option of dropping labels with missing data via the
    **dropna** function.

.. _basics.dataframe:

DataFrame
---------

**DataFrame** is a 2-dimensional labeled data structure with columns of
potentially different types. You can think of it like a spreadsheet or SQL
table, or a dict of Series objects. It is generally the most commonly used
pandas object. Like Series, DataFrame accepts many different kinds of input:

 - Dict of 1D ndarrays, lists, dicts, or Series
 - 2-D numpy.ndarray
 - `Structured or record
   <http://docs.scipy.org/doc/numpy/user/basics.rec.html>`__ ndarray
 - Another DataFrame

Along with the data, you can optionally pass **index** (row labels) and
**columns** (column labels) arguments. If you pass an index and / or columns,
pyou are guaranteeing the index and / or columns of the resulting
DataFrame. Thus, a dict of Series plus a specific index will discard all data
not matching up to the passed index.

If axis labels are not passed, they will be constructed from the input data
based on common sense rules.

**From dict of Series or dicts**

the result **index** will be the **union** of the indexes of the various
Series. If there are any nested dicts, these will be first converted to
Series. If no columns are passed, the columns will be the sorted list of dict
keys.

.. ipython:: python

    d = {'one' : Series([1., 2., 3.], index=['a', 'b', 'c']),
         'two' : Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd'])}
    df = DataFrame(d)
    df

    DataFrame(d, index=['d', 'b', 'a'])
    DataFrame(d, index=['d', 'b', 'a'], columns=['two', 'three'])

The row and column labels can be accessed respectively by accessing the
**index** and **columns** attributes:

.. ipython:: python

   df.index
   df.columns

**From dict of ndarrays / lists**

The ndarrays must all be the same length. If an index is passed, it must
clearly also be the same length as the arrays. If no index is passed, the
result will be ``range(n)``, where ``n`` is the array length.

.. ipython:: python

   d = {'one' : [1., 2., 3., 4.],
        'two' : [4., 3., 2., 1.]}
   DataFrame(d)
   DataFrame(d, index=['a', 'b', 'c', 'd'])

**From structured or record array**

This case is handled identically to a dict of arrays.

.. ipython:: python

   data = np.zeros((2,),dtype=[('A', 'i4'),('B', 'f4'),('C', 'a10')])
   data[:] = [(1,2.,'Hello'),(2,3.,"World")]

   DataFrame(data)
   DataFrame(data, index=['first', 'second'])
   DataFrame(data, columns=['C', 'A', 'B'])

.. note::

    DataFrame is not intended to work exactly like a 2-dimensional NumPy
    ndarray.

Column selection, addition, deletion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can treat a DataFrame semantically like a dict of like-indexed Series
objects. Getting, setting, and deleting columns works with the same syntax as
the analogous dict operations:

.. ipython:: python

   df['one']
   df['three'] = df['one'] * df['two']
   df['flag'] = df['one'] > 2
   df

Columns can be deleted or popped like with a dict:

.. ipython:: python

   del df['two']
   three = df.pop('three')
   df

When inserting a scalar value, it will naturally be propagated to fill the
column:

.. ipython:: python

   df['foo'] = 'bar'
   df

When inserting a Series that does not have the same index as the DataFrame, it
will be conformed to the DataFrame's index:

.. ipython:: python

   df['one_trunc'] = df['one'][:2]
   df

You can insert raw ndarrays but their length must match the length of the
DataFrame's index.

By default, columns get inserted at the end. The ``insert`` function is
available to insert at a particular location in the columns:

.. ipython:: python

   df.insert(1, 'bar', df['one'])
   df

Indexing / Selection
~~~~~~~~~~~~~~~~~~~~
The basics of indexing are as follows:

.. csv-table::
    :header: "Operation", "Syntax", "Result"
    :widths: 30, 20, 10

    Select column, ``df[col]``, Series
    Select row by label, ``df.xs(label)`` or ``df.ix[label]``, Series
    Select row by location (int), ``df.ix[loc]``, Series
    Slice rows, ``df[5:10]``, DataFrame
    Select rows by boolean vector, ``df[bool_vec]``, DataFrame

Row selection, for example, returns a Series whose index is the columns of the
DataFrame:

.. ipython:: python

   df.xs('b')
   df.ix[2]

Note if a DataFrame contains columns of multiple dtypes, the dtype of the row
will be chosen to accommodate all of the data types (dtype=object is the most
general).

For a more exhaustive treatment of more sophisticated label-based indexing and
slicing, see the :ref:`section on indexing <indexing>`. We will address the
fundamentals of reindexing / conforming to new sets of lables in the
:ref:`section on reindexing <basics.reindexing>`.

Data alignment and arithmetic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data alignment between DataFrame objects automatically align on **both the
columns and the index (row labels)**. Again, the resulting object will have the
union of the column and row labels.

.. ipython:: python

    df = DataFrame(randn(10, 4), columns=['A', 'B', 'C', 'D'])
    df2 = DataFrame(randn(7, 3), columns=['A', 'B', 'C'])
    df + df2

When doing an operation between DataFrame and Series, the default behavior is
to align the Series **index** on the DataFrame **columns**, thus `broadcasting
<http://docs.scipy.org/doc/numpy/user/basics.broadcasting.html>`__
row-wise. For example:

.. ipython:: python

   df - df.ix[0]

In the special case of working with time series data, if the Series is a
TimeSeries (which it will be automatically if the index contains datetime
objects), and the DataFrame index also contains dates, the broadcasting will be
column-wise:

.. ipython:: python

   index = DateRange('1/1/2000', periods=8)
   df = DataFrame(randn(8, 3), index=index,
                  columns=['A', 'B', 'C'])
   df
   type(df['A'])
   df - df['A']

Technical purity aside, this case is so common in practice that supporting the
special case is preferable to the alternative of forcing the user to transpose
and do column-based alignment like so:

.. ipython:: python

   (df.T - df['A']).T

For explicit control over the matching and broadcasting behavior, see the
section on :ref:`flexible binary operations <basics.binop>`.

Operations with scalars are just as you would expect:

.. ipython:: python

   df * 5 + 2
   1 / df
   df ** 4

Transposing
~~~~~~~~~~~

To transpose, access the ``T`` attribute (also the ``transpose`` function),
similar to an ndarray:

.. ipython:: python

   # only show the first 5 rows
   df[:5].T

DataFrame interoperability with NumPy functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Elementwise NumPy ufuncs (log, exp, sqrt, ...) and various other NumPy functions
can be used with no issues on DataFrame, assuming the data within are numeric:

.. ipython:: python

   np.exp(df)
   np.asarray(df)

DataFrame is not intended to be a drop-in replacement for ndarray as its
indexing semantics are quite different in places from a matrix.

Console display
~~~~~~~~~~~~~~~

For very large DataFrame objects, only a summary will be printed to the console
(here I am reading a CSV version of the **baseball** dataset from the **plyr**
R package):

.. ipython:: python

   baseball = read_csv('baseball.csv')
   baseball

However, using ``to_string`` will display any DataFrame in tabular form, though
it won't always fit the console width:

.. ipython:: python

   baseball.ix[-20:, :12].to_string()

.. _basics.panel:

WidePanel
---------

WidePanel is a somewhat less-used, but still important container for
3-dimensional data. The term `panel data
<http://en.wikipedia.org/wiki/Panel_data>`__ is derived from econometrics and
is partially responsible for the name pandas: pan(el)-da(ta)-s. The names for
the 3 axes are intended to give some semantic meaning to describing operations
involving panel data and, in particular, econometric analysis of panel
data. However, for the strict purposes of slicing and dicing a collection of
DataFrame objects, you may find the axis names slightly arbitrary:

  - **items**: axis 0, each item corresponds to a DataFrame contained inside
  - **major_axis**: axis 1, it is the **index** (rows) of each of the
    DataFrames
  - **minor_axis**: axis 2, it is the **columns** of each of the DataFrames

.. note::

    The "wide" in **WidePanel** name comes from the notion of "long" and "wide"
    formats of grouped data. The R `reshape function
    <http://stat.ethz.ch/R-manual/R-patched/library/stats/html/reshape.html>`__
    has some more to say about these.

Construction of WidePanels works about like you would expect:

**3D ndarray with optional axis labels**

.. ipython:: python

   wp = WidePanel(randn(2, 5, 4), items=['Item1', 'Item2'],
                  major_axis=DateRange('1/1/2000', periods=5),
                  minor_axis=['A', 'B', 'C', 'D'])
   wp


**dict of DataFrame objects**

.. ipython:: python

   data = {'Item1' : DataFrame(randn(4, 3)),
           'Item2' : DataFrame(randn(4, 2))}
   WidePanel(data)

Note that the values in the dict need only be **convertible to
DataFrame**. Thus, they can be any of the other valid inputs to DataFrame as
per above.

.. note::

   Unfortunately WidePanel, being less commonly used than Series and DataFrame,
   has been slightly neglected feature-wise. A number of methods and options
   available in DataFrame are not available in WidePanel. This will get worked
   on, of course, in future releases. And faster if you join me in working on
   the codebase.

Item selection / addition / deletion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to DataFrame functioning as a dict of Series, WidePanel is like a dict
of DataFrames:

.. ipython:: python

   wp['Item1']
   wp['Item3'] = wp['Item1'] / wp['Item2']

The API for insertion and deletion is the same as for DataFrame.

Indexing / Selection
~~~~~~~~~~~~~~~~~~~~

As of this writing, indexing with WidePanel is a bit more restrictive than in
DataFrame. Notably, :ref:`fancy indexing <indexing>` via the **ix** property
has not yet been integrated in WidePanel. This will be done, however, in a
future release.

.. csv-table::
    :header: "Operation", "Syntax", "Result"
    :widths: 30, 20, 10

    Select item, ``wp[item]``, DataFrame
    Get slice at major_axis label, ``wp.major_xs(val)``, DataFrame
    Get slice at minor_axis label, ``wp.minor_xs(val)``, DataFrame

For example, using the earlier example data, we could do:

.. ipython:: python

    wp['Item1']
    wp.major_xs(wp.major_axis[2])
    wp.minor_axis
    wp.minor_xs('C')

.. _basics.attrs:

Attributes and the raw ndarray(s)
---------------------------------

pandas objects have a number of attributes enabling you to access the metadata

  * **shape**: gives the axis dimensions of the object, consistent with ndarray
  * Axis labels

    * **Series**: *index* (only axis)
    * **DataFrame**: *index* (rows) and *columns*
    * **WidePanel**: *items*, *major_axis*, and *minor_axis*

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

If a DataFrame or WidePanel contains homogeneously-typed data, the ndarray can
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

With binary operations between pandas data structures, we have a couple items
of interest:

  * How to describe broadcasting behavior between higher- (e.g. DataFrame) and
    lower-dimensional (e.g. Series) objects.
  * Behavior of missing data in computations

We will demonstrate the currently-available functions to illustrate these
issues independently, though they can be performed simultaneously.

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

With WidePanel, describing the matching behavior is a bit more difficult, so
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
   match the broadcasting behavior of WidePanel. Though it would require a
   transition period so users can change their code...

Missing data / operations with fill values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Series and DataFrame (though not yet in WidePanel), the arithmetic functions
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

.. _basics.stats:

Descriptive statistics
----------------------

A large number of methods for computing descriptive statistics and other related
operations on :ref:`Series <api.series.stats>`, :ref:`DataFrame
<api.dataframe.stats>`, and :ref:`WidePanel <api.panel.stats>`. Most of these
are aggregations (hence producing a lower-dimensional result) like **sum**,
**mean**, and **quantile**, but some of them, like **cumsum** and **cumprod**,
produce an object of the same size. Generally speaking, these methods take an
**axis** argument, just like *ndarray.{sum, std, ...}*, but the axis can be
specified by name or integer:

  - **Series**: no axis argument needed
  - **DataFrame**: "index" (axis=0, default), "columns" (axis=1)
  - **WidePanel**: "items" (axis=0), "major" (axis=1, default), "minor"
    (axis=2)

For example:

.. ipython:: python

   df
   df.mean(0)
   df.mean(1)

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

Here is a quick reference summary table of common functions

.. csv-table::
    :header: "Function", "Description"
    :widths: 20, 80

    ``count``, Number of non-null observations
    ``sum``, Sum of values
    ``mean``, Mean of values
    ``median``, Arithmetic median of values
    ``min``, Minimum
    ``max``, Maximum
    ``prod``, Product of values
    ``std``, Unbiased standard deviation
    ``var``, Unbiased variance
    ``skew``, Unbiased skewness (3rd moment)
    ``kurt``, Unbiased kurtosis (4th moment)
    ``quantile``, Sample quantile (value at %)
    ``cumsum``, Cumulative sum
    ``cumprod``, Cumulative product

Summarizing data: describe
~~~~~~~~~~~~~~~~~~~~~~~~~~

For floating point data, there is a convenient ``describe`` function which
computes a variety of summary statistics about a Series or the columns of a
DataFrame (excluding NAs of course):

.. ipython:: python

    series = Series(randn(1000))
    series[::2] = np.nan
    series.describe()
    frame = DataFrame(randn(1000, 5))
    frame.ix[::2] = np.nan
    frame.describe()

Correlations between objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several handy methods for computing correlations are provided. The only
behavior available at the moment is to compute "pairwise complete
observations". In computing correlations in the presence of missing data, one
must be careful internally to compute the standard deviation of each Series
over the labels with valid data in both objects.

.. ipython:: python

   # Series with Series
   frame[0].corr(frame[0])

   # Pairwise correlation of DataFrame columns
   frame.corr()

A related method ``corrwith`` is implemented on DataFrame to compute the
correlation between like-labeled Series contained in different DataFrame
objects.

.. ipython:: python

   index = ['a', 'b', 'c', 'd', 'e']
   columns = ['one', 'two', 'three', 'four']
   df1 = DataFrame(randn(5, 4), index=index, columns=columns)
   df2 = DataFrame(randn(4, 4), index=index[:4], columns=columns)
   df1.corrwith(df2)
   df2.corrwith(df1, axis=1)

.. _basics.apply:

Function application
--------------------

Arbitrary functions can be applied along the axes of a DataFrame or WidePanel
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

   :ref:`Fancy indexing <indexing.fancy>` is an even more concise way of doing
   reindexing.

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
another object. While the syntax for this is straightforwad albeit verbose, it
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

Iteration
---------

Considering the pandas as somewhat dict-like structure, basic iteration
produces the "keys" of the objects, namely:

  * **Series**: the index label
  * **DataFrame**: the column labels
  * **WidePanel**: the item labels

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
  * **WidePanel**: (item, DataFrame) pairs

For example:

.. ipython::

   In [0]: for item, frame in wp.iteritems():
      ...:     print item
      ...:     print frame
      ...:

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

.. _basics.serialize:

Pickling and serialization
--------------------------

