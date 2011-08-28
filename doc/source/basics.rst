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

Series
------

:class:`Series` is a one-dimensional labeled array. The axis labels are
collectively referred to as the **index**. The basic method to create a Series
is to call:

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

   Series(np.random.randn(5))

**Case 2:** If **data** is a dict, if **index** is passed the values in data**
**corresponding to the labels in the index will be pulled out. Otherwise, an
**index will be constructed from the sorted keys of the dict, if possible.

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

If a label is not contained, an exception

.. code-block:: python

    >>> s['f']
    KeyError: 'f'

    >>> s.get('f')
    nan



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
    DataFrame(d)
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

Column selection, addition, deletion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

For a more exhaustive treatment of more sophisticated label-based indexing and
slicing, see the `section on indexing <indexing>`__.

WidePanel
---------

WidePanel is a less-used, but still important container for 3-dimensional
data. The term `panel data <http://en.wikipedia.org/wiki/Panel_data>`__ is
derived from econometrics and is partially responsible for the name pandas:
pan(el)-da(ta)-s.

Binary operations between objects
---------------------------------

Alignment and reindexing
------------------------

Deleting labels from an axis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Function application and basic stats
------------------------------------

Copying, type casting
---------------------

Pickling and serialization
--------------------------


