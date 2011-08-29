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
   np.set_printoptions(precision=4, suppress=True)

.. ipython:: python

   import numpy as np
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

Here, **data** can be many different things:

 - a Python dict
 - an ndarray
 - a scalar value (like 5)

The passed **index** is a list of axis labels. Thus, this separates into a few
cases depending on what **data is**:

**From ndarray**

If **data** is an ndarray, **index** must be the same length as **data**. If no
index is passed, one will be created having values ``[0, ..., len(data) - 1]``.

.. ipython:: python

   s = Series(np.random.randn(5), index=['a', 'b', 'c', 'd', 'e'])
   s
   s.index

   Series(np.random.randn(5))

**From dict**

If **data** is a dict, if **index** is passed the values in data corresponding
to the labels in the index will be pulled out. Otherwise, an index will be
constructed from the sorted keys of the dict, if possible.

.. ipython:: python

   d = {'a' : 0., 'b' : 1., 'c' : 2.}
   Series(d)
   Series(d, index=['b', 'c', 'd', 'a'])

.. note::

    NaN (not a number) is the standard missing data marker used in pandas

**From scalar value** If **data** is a scalar value, an index must be provided. The value
will be repeated to match the length of **index**

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

By default, columns get inserted at the end. The **insert**
function is available to insert at a particular location in the columns:

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

    df = DataFrame(np.random.randn(10, 4), columns=['A', 'B', 'C', 'D'])
    df2 = DataFrame(np.random.randn(7, 3), columns=['A', 'B', 'C'])
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
   df = DataFrame(np.random.randn(8, 3), index=index,
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

Console display
~~~~~~~~~~~~~~~

For very large DataFrame objects, only a summary will be printed to the console
(here I am reading a CSV version of the **baseball** dataset from the **plyr**
R package):

.. ipython:: python

   baseball = read_csv('baseball.csv')
   baseball

However, using **to_string** will display any DataFrame in tabular form, though
it won't always fit the he console width:

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

   wp = WidePanel(np.random.randn(2, 5, 4), items=['Item1', 'Item2'],
                  major_axis=DateRange('1/1/2000', periods=5),
                  minor_axis=['A', 'B', 'C', 'D'])
   wp


**dict of DataFrame objects**

.. ipython:: python

   data = {'Item1' : DataFrame(np.random.randn(4, 3)),
           'Item2' : DataFrame(np.random.randn(4, 2))}
   WidePanel(data)

Note that the values in the dict need only be **convertible to
DataFrame**. Thus, they can be any of the other valid inputs to DataFrame as
per above.

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

.. _basics.apply:

Function application and basic stats
------------------------------------

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

   d = {'one' : Series(np.random.randn(3), index=['a', 'b', 'c']),
        'two' : Series(np.random.randn(4), index=['a', 'b', 'c', 'd']),
        'three' : Series(np.random.randn(3), index=['b', 'c', 'd'])}
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
to specify the *broadcast axis*.

Missing data / operations with fill values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


.. _basics.reindexing:

Reindexing
----------



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

Dropping labels from an axis
----------------------------

.. _basics.cast:

Copying, type casting
---------------------

.. _basics.serialize:

Pickling and serialization
--------------------------

