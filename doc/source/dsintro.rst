.. currentmodule:: pandas
.. _dsintro:

************************
Intro to Data Structures
************************

We'll start with a quick, non-comprehensive overview of the fundamental data
structures in pandas to get you started. The fundamental behavior about data
types, indexing, and axis labeling / alignment apply across all of the
objects. To get started, import numpy and load pandas into your namespace:

.. ipython:: python
   :suppress:

   import numpy as np
   from pandas import *
   randn = np.random.randn
   np.set_printoptions(precision=4, suppress=True)
   set_option('display.precision', 4, 'display.max_columns', 8)

.. ipython:: python

   import numpy as np
   # will use a lot in examples
   randn = np.random.randn
   from pandas import *

Here is a basic tenet to keep in mind: **data alignment is intrinsic**. The link
between labels and data will not be broken unless done so explicitly by you.

We'll give a brief intro to the data structures, then consider all of the broad
categories of functionality and methods in separate sections.

When using pandas, we recommend the following import convention:

.. code-block:: python

   import pandas as pd


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

    Starting in v0.8.0, pandas supports non-unique index values. If an operation
    that does not support duplicate index values is attempted, an exception
    will be raised at that time. The reason for being lazy is nearly all performance-based
    (there are many instances in computations, like parts of GroupBy, where the index
    is not used).

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

We will address array-based indexing in a separate :ref:`section <indexing>`.

Series is dict-like
~~~~~~~~~~~~~~~~~~~

A Series is like a fixed-size dict in that you can get and set values by index
label:

.. ipython :: python

    s['a']
    s['e'] = 12.
    s
    'e' in s
    'f' in s

If a label is not contained, an exception is raised:

.. code-block:: python

    >>> s['f']
    KeyError: 'f'

Using the ``get`` method, a missing label will return None or specified default:

.. ipython:: python

   s.get('f')

   s.get('f', np.nan)

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

Name attribute
~~~~~~~~~~~~~~

.. _dsintro.name_attribute:

Series can also have a ``name`` attribute:

.. ipython:: python

   s = Series(np.random.randn(5), name='something')
   s
   s.name

The Series ``name`` will be assigned automatically in many cases, in particular
when taking 1D slices of DataFrame as you will see below.

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
 - A ``Series``
 - Another ``DataFrame``

Along with the data, you can optionally pass **index** (row labels) and
**columns** (column labels) arguments. If you pass an index and / or columns,
you are guaranteeing the index and / or columns of the resulting
DataFrame. Thus, a dict of Series plus a specific index will discard all data
not matching up to the passed index.

If axis labels are not passed, they will be constructed from the input data
based on common sense rules.

From dict of Series or dicts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The result **index** will be the **union** of the indexes of the various
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

.. note::

   When a particular set of columns is passed along with a dict of data, the
   passed columns override the keys in the dict.

.. ipython:: python

   df.index
   df.columns

From dict of ndarrays / lists
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ndarrays must all be the same length. If an index is passed, it must
clearly also be the same length as the arrays. If no index is passed, the
result will be ``range(n)``, where ``n`` is the array length.

.. ipython:: python

   d = {'one' : [1., 2., 3., 4.],
        'two' : [4., 3., 2., 1.]}
   DataFrame(d)
   DataFrame(d, index=['a', 'b', 'c', 'd'])

From structured or record array
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. _basics.dataframe.from_list_of_dicts:

From a list of dicts
~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

   data2 = [{'a': 1, 'b': 2}, {'a': 5, 'b': 10, 'c': 20}]
   DataFrame(data2)
   DataFrame(data2, index=['first', 'second'])
   DataFrame(data2, columns=['a', 'b'])

.. _basics.dataframe.from_series:

From a Series
~~~~~~~~~~~~~

The result will be a DataFrame with the same index as the input Series, and
with one column whose name is the original name of the Series (only if no other
column name provided).

**Missing Data**

Much more will be said on this topic in the :ref:`Missing data <missing_data>`
section. To construct a DataFrame with missing data, use ``np.nan`` for those
values which are missing. Alternatively, you may pass a ``numpy.MaskedArray``
as the data argument to the DataFrame constructor, and its masked entries will
be considered missing.

Alternate Constructors
~~~~~~~~~~~~~~~~~~~~~~

.. _basics.dataframe.from_dict:

**DataFrame.from_dict**

``DataFrame.from_dict`` takes a dict of dicts or a dict of array-like sequences
and returns a DataFrame. It operates like the ``DataFrame`` constructor except
for the ``orient`` parameter which is ``'columns'`` by default, but which can be
set to ``'index'`` in order to use the dict keys as row labels.

.. _basics.dataframe.from_records:

**DataFrame.from_records**

``DataFrame.from_records`` takes a list of tuples or an ndarray with structured
dtype. Works analogously to the normal ``DataFrame`` constructor, except that
index maybe be a specific field of the structured dtype to use as the index.
For example:

.. ipython:: python

   data
   DataFrame.from_records(data, index='C')

.. _basics.dataframe.from_items:

**DataFrame.from_items**

``DataFrame.from_items`` works analogously to the form of the ``dict``
constructor that takes a sequence of ``(key, value)`` pairs, where the keys are
column (or row, in the case of ``orient='index'``) names, and the value are the
column values (or row values). This can be useful for constructing a DataFrame
with the columns in a particular order without having to pass an explicit list
of columns:

.. ipython:: python

   DataFrame.from_items([('A', [1, 2, 3]), ('B', [4, 5, 6])])

If you pass ``orient='index'``, the keys will be the row labels. But in this
case you must also pass the desired column names:

.. ipython:: python

   DataFrame.from_items([('A', [1, 2, 3]), ('B', [4, 5, 6])],
                        orient='index', columns=['one', 'two', 'three'])

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
    Select row by label, ``df.loc[label]``, Series
    Select row by integer location, ``df.iloc[loc]``, Series
    Slice rows, ``df[5:10]``, DataFrame
    Select rows by boolean vector, ``df[bool_vec]``, DataFrame

Row selection, for example, returns a Series whose index is the columns of the
DataFrame:

.. ipython:: python

   df.loc['b']
   df.iloc[2]

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

   df - df.iloc[0]

In the special case of working with time series data, if the Series is a
TimeSeries (which it will be automatically if the index contains datetime
objects), and the DataFrame index also contains dates, the broadcasting will be
column-wise:

.. ipython:: python

   index = date_range('1/1/2000', periods=8)
   df = DataFrame(randn(8, 3), index=index, columns=list('ABC'))
   df
   type(df['A'])
   df - df['A']

.. warning::

   .. code-block:: python

      df - df['A']

   is now deprecated and will be removed in a future release. The preferred way
   to replicate this behavior is

   .. code-block:: python

      df.sub(df['A'], axis=0)

For explicit control over the matching and broadcasting behavior, see the
section on :ref:`flexible binary operations <basics.binop>`.

Operations with scalars are just as you would expect:

.. ipython:: python

   df * 5 + 2
   1 / df
   df ** 4

.. _dsintro.boolean:

Boolean operators work as well:

.. ipython:: python

   df1 = DataFrame({'a' : [1, 0, 1], 'b' : [0, 1, 1] }, dtype=bool)
   df2 = DataFrame({'a' : [0, 1, 1], 'b' : [1, 1, 0] }, dtype=bool)
   df1 & df2
   df1 | df2
   df1 ^ df2
   -df1

Transposing
~~~~~~~~~~~

To transpose, access the ``T`` attribute (also the ``transpose`` function),
similar to an ndarray:

.. ipython:: python

   # only show the first 5 rows
   df[:5].T

DataFrame interoperability with NumPy functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _dsintro.numpy_interop:

Elementwise NumPy ufuncs (log, exp, sqrt, ...) and various other NumPy functions
can be used with no issues on DataFrame, assuming the data within are numeric:

.. ipython:: python

   np.exp(df)
   np.asarray(df)

The dot method on DataFrame implements matrix multiplication:

.. ipython:: python

   df.T.dot(df)

Similarly, the dot method on Series implements dot product:

.. ipython:: python

   s1 = Series(np.arange(5,10))
   s1.dot(s1)

DataFrame is not intended to be a drop-in replacement for ndarray as its
indexing semantics are quite different in places from a matrix.

Console display
~~~~~~~~~~~~~~~

For very large DataFrame objects, only a summary will be printed to the console
(here I am reading a CSV version of the **baseball** dataset from the **plyr**
R package):

.. ipython:: python
   :suppress:

   # force a summary to be printed
   pd.set_option('display.max_rows', 5)

.. ipython:: python

   baseball = read_csv('data/baseball.csv')
   print baseball

.. ipython:: python
   :suppress:

   # restore GlobalPrintConfig
   pd.reset_option('^display\.')

However, using ``to_string`` will return a string representation of the
DataFrame in tabular form, though it won't always fit the console width:

.. ipython:: python

   print baseball.iloc[-20:, :12].to_string()

New since 0.10.0, wide DataFrames will now be printed across multiple rows by
default:

.. ipython:: python

   DataFrame(randn(3, 12))

You can change how much to print on a single row by setting the ``line_width``
option:

.. ipython:: python

   set_option('line_width', 40) # default is 80

   DataFrame(randn(3, 12))

.. ipython:: python
   :suppress:

   reset_option('line_width')

You can also disable this feature via the ``expand_frame_repr`` option:

.. ipython:: python

   set_option('expand_frame_repr', False)

   DataFrame(randn(3, 12))

.. ipython:: python
   :suppress:

   reset_option('expand_frame_repr')


DataFrame column attribute access and IPython completion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a DataFrame column label is a valid Python variable name, the column can be
accessed like attributes:

.. ipython:: python

   df = DataFrame({'foo1' : np.random.randn(5),
                   'foo2' : np.random.randn(5)})
   df
   df.foo1

The columns are also connected to the `IPython <http://ipython.org>`__
completion mechanism so they can be tab-completed:

.. code-block:: ipython

    In [5]: df.fo<TAB>
    df.foo1  df.foo2

.. _basics.panel:

Panel
-----

Panel is a somewhat less-used, but still important container for 3-dimensional
data. The term `panel data <http://en.wikipedia.org/wiki/Panel_data>`__ is
derived from econometrics and is partially responsible for the name pandas:
pan(el)-da(ta)-s. The names for the 3 axes are intended to give some semantic
meaning to describing operations involving panel data and, in particular,
econometric analysis of panel data. However, for the strict purposes of slicing
and dicing a collection of DataFrame objects, you may find the axis names
slightly arbitrary:

  - **items**: axis 0, each item corresponds to a DataFrame contained inside
  - **major_axis**: axis 1, it is the **index** (rows) of each of the
    DataFrames
  - **minor_axis**: axis 2, it is the **columns** of each of the DataFrames

Construction of Panels works about like you would expect:

From 3D ndarray with optional axis labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

   wp = Panel(randn(2, 5, 4), items=['Item1', 'Item2'],
              major_axis=date_range('1/1/2000', periods=5),
              minor_axis=['A', 'B', 'C', 'D'])
   wp


From dict of DataFrame objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

   data = {'Item1' : DataFrame(randn(4, 3)),
           'Item2' : DataFrame(randn(4, 2))}
   Panel(data)

Note that the values in the dict need only be **convertible to
DataFrame**. Thus, they can be any of the other valid inputs to DataFrame as
per above.

One helpful factory method is ``Panel.from_dict``, which takes a
dictionary of DataFrames as above, and the following named parameters:

.. csv-table::
   :header: "Parameter", "Default", "Description"
   :widths: 10, 10, 40

   intersect, ``False``, drops elements whose indices do not align
   orient, ``items``, use ``minor`` to use DataFrames' columns as panel items

For example, compare to the construction above:

.. ipython:: python

   Panel.from_dict(data, orient='minor')

Orient is especially useful for mixed-type DataFrames. If you pass a dict of
DataFrame objects with mixed-type columns, all of the data will get upcasted to
``dtype=object`` unless you pass ``orient='minor'``:

.. ipython:: python

   df = DataFrame({'a': ['foo', 'bar', 'baz'],
                   'b': np.random.randn(3)})
   df
   data = {'item1': df, 'item2': df}
   panel = Panel.from_dict(data, orient='minor')
   panel['a']
   panel['b']
   panel['b'].dtypes

.. note::

   Unfortunately Panel, being less commonly used than Series and DataFrame,
   has been slightly neglected feature-wise. A number of methods and options
   available in DataFrame are not available in Panel. This will get worked
   on, of course, in future releases. And faster if you join me in working on
   the codebase.

.. _dsintro.to_panel:

From DataFrame using ``to_panel`` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This method was introduced in v0.7 to replace ``LongPanel.to_long``, and converts
a DataFrame with a two-level index to a Panel.

.. ipython:: python

   midx = MultiIndex(levels=[['one', 'two'], ['x','y']], labels=[[1,1,0,0],[1,0,1,0]])
   df = DataFrame({'A' : [1, 2, 3, 4], 'B': [5, 6, 7, 8]}, index=midx)
   df.to_panel()

.. _dsintro.panel_item_selection:

Item selection / addition / deletion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to DataFrame functioning as a dict of Series, Panel is like a dict
of DataFrames:

.. ipython:: python

   wp['Item1']
   wp['Item3'] = wp['Item1'] / wp['Item2']

The API for insertion and deletion is the same as for DataFrame. And as with
DataFrame, if the item is a valid python identifier, you can access it as an
attribute and tab-complete it in IPython.

Transposing
~~~~~~~~~~~

A Panel can be rearranged using its ``transpose`` method (which does not make a
copy by default unless the data are heterogeneous):

.. ipython:: python

   wp.transpose(2, 0, 1)

Indexing / Selection
~~~~~~~~~~~~~~~~~~~~

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

Squeezing
~~~~~~~~~

Another way to change the dimensionality of an object is to ``squeeze`` a 1-len object, similar to ``wp['Item1']``

.. ipython:: python

   wp.reindex(items=['Item1']).squeeze()
   wp.reindex(items=['Item1'],minor=['B']).squeeze()


Conversion to DataFrame
~~~~~~~~~~~~~~~~~~~~~~~

A Panel can be represented in 2D form as a hierarchically indexed
DataFrame. See the section :ref:`hierarchical indexing <indexing.hierarchical>`
for more on this. To convert a Panel to a DataFrame, use the ``to_frame``
method:

.. ipython:: python

   panel = Panel(np.random.randn(3, 5, 4), items=['one', 'two', 'three'],
                 major_axis=date_range('1/1/2000', periods=5),
                 minor_axis=['a', 'b', 'c', 'd'])
   panel.to_frame()


.. _dsintro.panel4d:

Panel4D (Experimental)
----------------------

``Panel4D`` is a 4-Dimensional named container very much like a ``Panel``, but
having 4 named dimensions. It is intended as a test bed for more N-Dimensional named
containers.

  - **labels**: axis 0, each item corresponds to a Panel contained inside
  - **items**: axis 1, each item corresponds to a DataFrame contained inside
  - **major_axis**: axis 2, it is the **index** (rows) of each of the
    DataFrames
  - **minor_axis**: axis 3, it is the **columns** of each of the DataFrames


``Panel4D`` is a sub-class of ``Panel``, so most methods that work on Panels are
applicable to Panel4D. The following methods are disabled:

  - ``join , to_frame , to_excel , to_sparse , groupby``

Construction of Panel4D works in a very similar manner to a ``Panel``

From 4D ndarray with optional axis labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

   p4d = Panel4D(randn(2, 2, 5, 4),
              labels=['Label1','Label2'],
              items=['Item1', 'Item2'],
              major_axis=date_range('1/1/2000', periods=5),
              minor_axis=['A', 'B', 'C', 'D'])
   p4d


From dict of Panel objects
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

   data = { 'Label1' : Panel({ 'Item1' : DataFrame(randn(4, 3)) }),
            'Label2' : Panel({ 'Item2' : DataFrame(randn(4, 2)) }) }
   Panel4D(data)

Note that the values in the dict need only be **convertible to Panels**.
Thus, they can be any of the other valid inputs to Panel as per above.

Slicing
~~~~~~~

Slicing works in a similar manner to a Panel. ``[]`` slices the first dimension.
``.ix`` allows you to slice abitrarily and get back lower dimensional objects

.. ipython:: python

   p4d['Label1']

4D -> Panel

.. ipython:: python

   p4d.ix[:,:,:,'A']

4D -> DataFrame

.. ipython:: python

   p4d.ix[:,:,0,'A']

4D -> Series

.. ipython:: python

   p4d.ix[:,0,0,'A']

Transposing
~~~~~~~~~~~

A Panel4D can be rearranged using its ``transpose`` method (which does not make a
copy by default unless the data are heterogeneous):

.. ipython:: python

   p4d.transpose(3, 2, 1, 0)

PanelND (Experimental)
----------------------

PanelND is a module with a set of factory functions to enable a user to construct N-dimensional named
containers like Panel4D, with a custom set of axis labels. Thus a domain-specific container can easily be
created.

The following creates a Panel5D. A new panel type object must be sliceable into a lower dimensional object.
Here we slice to a Panel4D.

.. ipython:: python

    from pandas.core import panelnd
    Panel5D = panelnd.create_nd_panel_factory(
        klass_name   = 'Panel5D',
        axis_orders  = [ 'cool', 'labels','items','major_axis','minor_axis'],
        axis_slices  = { 'labels' : 'labels', 'items' : 'items',
	                 'major_axis' : 'major_axis', 'minor_axis' : 'minor_axis' },
        slicer       = Panel4D,
        axis_aliases = { 'major' : 'major_axis', 'minor' : 'minor_axis' },
        stat_axis    = 2)

    p5d = Panel5D(dict(C1 = p4d))
    p5d

    # print a slice of our 5D
    p5d.ix['C1',:,:,0:3,:]

    # transpose it
    p5d.transpose(1,2,3,4,0)

    # look at the shape & dim
    p5d.shape
    p5d.ndim
