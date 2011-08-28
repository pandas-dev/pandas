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

The passed **index** is a list of axis labels. Thus, this separates into three
cases depending on what **data is**:

**Case 1:** If **data** is an ndarray, **index** must be the same length as
**data**. If no index is passed, one will be created having values ``[0, ...,
len(data) - 1]``.

.. ipython:: python

   s = Series(np.random.randn(5), index=['a', 'b', 'c', 'd', 'e'])
   s
   s.index

   Series(np.random.randn(5))

**Case 2:** If **data** is a dict, if **index** is passed the values in data
corresponding to the labels in the index will be pulled out. Otherwise, an
index will be constructed from the sorted keys of the dict, if possible.

.. ipython:: python

   d = {'a' : 0., 'b' : 1., 'c' : 2.}
   Series(d)
   Series(d, index=['b', 'c', 'd', 'a'])

.. note::

    NaN (not a number) is the standard missing data marker used in pandas

**Case 3:** If **data** is a scalar value, an index must be provided. The value
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

 - Dict of 1D ndarrays, lists, or Series
 - 2-D numpy.ndarray
 - `Structured or record
   <http://docs.scipy.org/doc/numpy/user/basics.rec.html>`__ ndarray
 - Another DataFrame

Along with the data, you can optionally pass **index** (row labels) and
**columns** (column labels) arguments. If you pass an index and / or columns,
you are guaranteeing the index and / or columns of the resulting
DataFrame. Thus, a dict of Series plus a specific index will discard all data
not matching up to the passed index.

If axis labels are not passed, they will be constructed from the input data
based on common sense rules. Here are the various cases:

**Case 1, dict of Series**: the result **index** will be the **union** of the
indexes of the various Series. If no columns are passed, the columns will be
the sorted list of dict keys.

.. ipython:: python

    d = {'one' : Series([1., 2., 3.], index=['a', 'b', 'c']),
         'two' : Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd'])}
    df = DataFrame(d)
    df

    DataFrame(d, index=['d', 'b', 'a'])
    DataFrame(d, index=['d', 'b', 'a'], columns=['two', 'three'])

**Case 2, dict of ndarrays / lists**: The ndarrays must all be the same
length. If an index is passed, it must clearly also be the same length as the
arrays. If no index is passed, the result will be ``range(n)``, where ``n`` is
the array length.

.. ipython:: python

   d = {'one' : [1., 2., 3., 4.],
        'two' : [4., 3., 2., 1.]}
   DataFrame(d)
   DataFrame(d, index=['a', 'b', 'c', 'd'])

**Case 3, structured or record array**: This case is handled identically to a
dict of arrays.

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
   del df['two']
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
slicing, see the :ref:`section on indexing <indexing>`.

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
DataFrame objects, the axis names are slightly arbitrary:

  - **items**: axis 0, each item corresponds to a DataFrame contained inside
  - **major_axis**: axis 1, it is the **index** (rows) of each of the
    DataFrames
  - **minor_axis**: axis 2, it is the **columns** of each of the DataFrames

.. note::

    The "wide" in **WidePanel** name comes from the notion of "long" and "wide"
    formats of grouped data. The R `reshape function
    <http://stat.ethz.ch/R-manual/R-patched/library/stats/html/reshape.html>`__
    has some more to say about these.


.. _basics.attrs:

Attributes and the raw ndarray(s)
---------------------------------

.. _basics.binop:

Flexible binary operations
--------------------------

Combining with fill values
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _basics.reindexing:

Reindexing
----------



Dropping labels from an axis
----------------------------

.. _basics.apply:

Function application and basic stats
------------------------------------

.. _basics.cast:

Copying, type casting
---------------------

.. _basics.serialize:

Pickling and serialization
--------------------------

