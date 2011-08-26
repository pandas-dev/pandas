.. currentmodule:: pandas
.. _basics:

*********************
Data structure basics
*********************

This is a quick, non-comprehensive overview of the fundamental data structures
in pandas to get you started. The fundamental behavior about data types,
indexing, and axis labeling / alignment apply across all of the objects. To get
started, load pandas into your namespace:

.. ipython:: python

   from pandas import *


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

   Series(np.random.randn(5), index=['a', 'b', 'c', 'd', 'e'])
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

DataFrame
---------

**DataFrame** is a 2-dimensional labeled data structure with columns of
potentially different types. You can think of it like a spreadsheet or SQL
table, or preferably a dict of Series objects. It is generally the most
commonly used pandas object.

WidePanel
---------
