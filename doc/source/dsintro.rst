.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   np.set_printoptions(precision=4, suppress=True)
   import pandas as pd
   pd.set_option('display.precision', 4, 'display.max_columns', 8)
   pd.options.display.max_rows = 15

   import matplotlib
   # matplotlib.style.use('default')
   import matplotlib.pyplot as plt
   plt.close('all')

.. _dsintro:

************************
Intro to Data Structures
************************

We'll start with a quick, non-comprehensive overview of the fundamental data
structures in pandas to get you started. The fundamental behavior about data
types, indexing, and axis labeling / alignment apply across all of the
objects. To get started, import NumPy and load pandas into your namespace:

.. ipython:: python

   import numpy as np
   import pandas as pd

Here is a basic tenet to keep in mind: **data alignment is intrinsic**. The link
between labels and data will not be broken unless done so explicitly by you.

We'll give a brief intro to the data structures, then consider all of the broad
categories of functionality and methods in separate sections.

.. _basics.series:

Series
------

:class:`Series` is a one-dimensional labeled array capable of holding any data
type (integers, strings, floating point numbers, Python objects, etc.). The axis
labels are collectively referred to as the **index**. The basic method to create a Series is to call:

::

    >>> s = pd.Series(data, index=index)

Here, ``data`` can be many different things:

* a Python dict
* an ndarray
* a scalar value (like 5)

The passed **index** is a list of axis labels. Thus, this separates into a few
cases depending on what **data is**:

**From ndarray**

If ``data`` is an ndarray, **index** must be the same length as **data**. If no
index is passed, one will be created having values ``[0, ..., len(data) - 1]``.

.. ipython:: python

   s = pd.Series(np.random.randn(5), index=['a', 'b', 'c', 'd', 'e'])
   s
   s.index

   pd.Series(np.random.randn(5))

.. note::

    pandas supports non-unique index values. If an operation
    that does not support duplicate index values is attempted, an exception
    will be raised at that time. The reason for being lazy is nearly all performance-based
    (there are many instances in computations, like parts of GroupBy, where the index
    is not used).

**From dict**

Series can be instantiated from dicts:

.. ipython:: python

   d = {'b' : 1, 'a' : 0, 'c' : 2}
   pd.Series(d)

.. note::

   When the data is a dict, and an index is not passed, the ``Series`` index
   will be ordered by the dict's insertion order, if you're using Python
   version >= 3.6 and Pandas version >= 0.23.

   If you're using Python < 3.6 or Pandas < 0.23, and an index is not passed,
   the ``Series`` index will be the lexically ordered list of dict keys.

In the example above, if you were on a Python version lower than 3.6 or a
Pandas version lower than 0.23, the ``Series`` would be ordered by the lexical
order of the dict keys (i.e. ``['a', 'b', 'c']`` rather than ``['b', 'a', 'c']``).

If an index is passed, the values in data corresponding to the labels in the
index will be pulled out.

.. ipython:: python

   d = {'a' : 0., 'b' : 1., 'c' : 2.}
   pd.Series(d)
   pd.Series(d, index=['b', 'c', 'd', 'a'])

.. note::

    NaN (not a number) is the standard missing data marker used in pandas.

**From scalar value**

If ``data`` is a scalar value, an index must be
provided. The value will be repeated to match the length of **index**.

.. ipython:: python

   pd.Series(5., index=['a', 'b', 'c', 'd', 'e'])

Series is ndarray-like
~~~~~~~~~~~~~~~~~~~~~~

``Series`` acts very similarly to a ``ndarray``, and is a valid argument to most NumPy functions.
However, operations such as slicing will also slice the index.

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

See also the :ref:`section on attribute access<indexing.attribute_access>`.

Vectorized operations and label alignment with Series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When working with raw NumPy arrays, looping through value-by-value is usually
not necessary. The same is true when working with Series in pandas.
Series can also be passed into most NumPy methods expecting an ndarray.

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
result will be marked as missing ``NaN``. Being able to write code without doing
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

   s = pd.Series(np.random.randn(5), name='something')
   s
   s.name

The Series ``name`` will be assigned automatically in many cases, in particular
when taking 1D slices of DataFrame as you will see below.

.. versionadded:: 0.18.0

You can rename a Series with the :meth:`pandas.Series.rename` method.

.. ipython:: python

   s2 = s.rename("different")
   s2.name

Note that ``s`` and ``s2`` refer to different objects.

.. _basics.dataframe:

DataFrame
---------

**DataFrame** is a 2-dimensional labeled data structure with columns of
potentially different types. You can think of it like a spreadsheet or SQL
table, or a dict of Series objects. It is generally the most commonly used
pandas object. Like Series, DataFrame accepts many different kinds of input:

* Dict of 1D ndarrays, lists, dicts, or Series
* 2-D numpy.ndarray
* `Structured or record
  <http://docs.scipy.org/doc/numpy/user/basics.rec.html>`__ ndarray
* A ``Series``
* Another ``DataFrame``

Along with the data, you can optionally pass **index** (row labels) and
**columns** (column labels) arguments. If you pass an index and / or columns,
you are guaranteeing the index and / or columns of the resulting
DataFrame. Thus, a dict of Series plus a specific index will discard all data
not matching up to the passed index.

If axis labels are not passed, they will be constructed from the input data
based on common sense rules.

.. note::

   When the data is a dict, and ``columns`` is not specified, the ``DataFrame``
   columns will be ordered by the dict's insertion order, if you are using
   Python version >= 3.6 and Pandas >= 0.23.

   If you are using Python < 3.6 or Pandas < 0.23, and ``columns`` is not
   specified, the ``DataFrame`` columns will be the lexically ordered list of dict
   keys.

From dict of Series or dicts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The resulting **index** will be the **union** of the indexes of the various
Series. If there are any nested dicts, these will first be converted to
Series. If no columns are passed, the columns will be the ordered list of dict
keys.

.. ipython:: python

    d = {'one' : pd.Series([1., 2., 3.], index=['a', 'b', 'c']),
         'two' : pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd'])}
    df = pd.DataFrame(d)
    df

    pd.DataFrame(d, index=['d', 'b', 'a'])
    pd.DataFrame(d, index=['d', 'b', 'a'], columns=['two', 'three'])

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
   pd.DataFrame(d)
   pd.DataFrame(d, index=['a', 'b', 'c', 'd'])

From structured or record array
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This case is handled identically to a dict of arrays.

.. ipython:: python

   data = np.zeros((2,), dtype=[('A', 'i4'),('B', 'f4'),('C', 'a10')])
   data[:] = [(1,2.,'Hello'), (2,3.,"World")]

   pd.DataFrame(data)
   pd.DataFrame(data, index=['first', 'second'])
   pd.DataFrame(data, columns=['C', 'A', 'B'])

.. note::

    DataFrame is not intended to work exactly like a 2-dimensional NumPy
    ndarray.

.. _basics.dataframe.from_list_of_dicts:

From a list of dicts
~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

   data2 = [{'a': 1, 'b': 2}, {'a': 5, 'b': 10, 'c': 20}]
   pd.DataFrame(data2)
   pd.DataFrame(data2, index=['first', 'second'])
   pd.DataFrame(data2, columns=['a', 'b'])

.. _basics.dataframe.from_dict_of_tuples:

From a dict of tuples
~~~~~~~~~~~~~~~~~~~~~

You can automatically create a MultiIndexed frame by passing a tuples
dictionary.

.. ipython:: python

   pd.DataFrame({('a', 'b'): {('A', 'B'): 1, ('A', 'C'): 2},
                 ('a', 'a'): {('A', 'C'): 3, ('A', 'B'): 4},
                 ('a', 'c'): {('A', 'B'): 5, ('A', 'C'): 6},
                 ('b', 'a'): {('A', 'C'): 7, ('A', 'B'): 8},
                 ('b', 'b'): {('A', 'D'): 9, ('A', 'B'): 10}})

.. _basics.dataframe.from_series:

From a Pandas Series
~~~~~~~~~~~~~~~~~~~~

The result will be a DataFrame with the same index as the input Series, and
with one column whose name is the original name of the Series (only if no other
column name provided).

.. ipython:: python

   df = pd.DataFrame()
   df
   x1 = pd.Series([1,2,3,4], index=[1,2,3,4])
   x1
   df['A'] = x1
   df
   x2 = pd.Series([5,5], index=[1, 3])
   x2
   df['B'] = x2
   df
   lst = [9, 9, 9, 9]
   x3 = pd.Series(lst, df.index)
   x3
   x4 = df['C']
   df['D'] = pd.Series(df['C']+1, df.index)
   df

From a Numpy Series
~~~~~~~~~~~~~~~~~~~

Numpy series is list or list of lists with no index.

.. ipython:: python

    x1 = np.random.randint(0, 10, (5))
    x1
    x2 = np.random.randint(-1.0, 10, (5))
    x2
    df1 = pd.DataFrame(x1, columns=['A'])
    df1
    df2 = pd.DataFrame([[x[0], x[1]] for x in zip(x1, x2)], columns=['B','C'])
    df2
    df3 = pd.DataFrame(np.random.randint(100, 999, (5, 2)), columns=['A', 'B'])
    df3
    x1 = np.random.rand(5, 2)
    x1
    df3['A'] = x1
    df3
    df3[['A', 'B']] = x1
    df3

**Missing Data**

Much more will be said on this topic in the :ref:`Missing data <missing_data>`
section. To construct a DataFrame with missing data, we use ``np.nan`` to
represent missing values. Alternatively, you may pass a ``numpy.MaskedArray``
as the data argument to the DataFrame constructor, and its masked entries will
be considered missing.

.. ipython:: python

   df = pd.DataFrame()
   df
   x1 = pd.Series([1,2,3,4], index=[1,2,3,4])
   x1
   df['A'] = x1
   df
   x2 = pd.Series([1, np.nan, np.nan, 4], df.index)
   x2
   
Alternate Constructors
~~~~~~~~~~~~~~~~~~~~~~

.. _basics.dataframe.from_dict:

**DataFrame.from_dict**

``DataFrame.from_dict`` takes a dict of dicts or a dict of array-like sequences
and returns a DataFrame. It operates like the ``DataFrame`` constructor except
for the ``orient`` parameter which is ``'columns'`` by default, but which can be
set to ``'index'`` in order to use the dict keys as row labels.


.. ipython:: python

   pd.DataFrame.from_dict(dict([('A', [1, 2, 3]), ('B', [4, 5, 6])]))

If you pass ``orient='index'``, the keys will be the row labels. In this
case, you can also pass the desired column names:

.. ipython:: python

   pd.DataFrame.from_dict(dict([('A', [1, 2, 3]), ('B', [4, 5, 6])]),
                          orient='index', columns=['one', 'two', 'three'])

.. _basics.dataframe.from_records:

**DataFrame.from_records**

``DataFrame.from_records`` takes a list of tuples or an ndarray with structured
dtype. It works analogously to the normal ``DataFrame`` constructor, except that
the resulting DataFrame index may be a specific field of the structured
dtype. For example:

.. ipython:: python

   data
   pd.DataFrame.from_records(data, index='C')


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

.. _dsintro.chained_assignment:

Assigning New Columns in Method Chains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inspired by `dplyr's
<http://cran.rstudio.com/web/packages/dplyr/vignettes/introduction.html#mutate>`__
``mutate`` verb, DataFrame has an :meth:`~pandas.DataFrame.assign`
method that allows you to easily create new columns that are potentially
derived from existing columns.

.. ipython:: python

   iris = pd.read_csv('data/iris.data')
   iris.head()

   (iris.assign(sepal_ratio = iris['SepalWidth'] / iris['SepalLength'])
        .head())

In the example above, we inserted a precomputed value. We can also pass in
a function of one argument to be evaluated on the DataFrame being assigned to.

.. ipython:: python

   iris.assign(sepal_ratio = lambda x: (x['SepalWidth'] /
                                        x['SepalLength'])).head()

``assign`` **always** returns a copy of the data, leaving the original
DataFrame untouched.

Passing a callable, as opposed to an actual value to be inserted, is
useful when you don't have a reference to the DataFrame at hand. This is
common when using ``assign`` in a chain of operations. For example,
we can limit the DataFrame to just those observations with a Sepal Length
greater than 5, calculate the ratio, and plot:

.. ipython:: python

   @savefig basics_assign.png
   (iris.query('SepalLength > 5')
        .assign(SepalRatio = lambda x: x.SepalWidth / x.SepalLength,
                PetalRatio = lambda x: x.PetalWidth / x.PetalLength)
        .plot(kind='scatter', x='SepalRatio', y='PetalRatio'))

Since a function is passed in, the function is computed on the DataFrame
being assigned to. Importantly, this is the DataFrame that's been filtered
to those rows with sepal length greater than 5. The filtering happens first,
and then the ratio calculations. This is an example where we didn't
have a reference to the *filtered* DataFrame available.

The function signature for ``assign`` is simply ``**kwargs``. The keys
are the column names for the new fields, and the values are either a value
to be inserted (for example, a ``Series`` or NumPy array), or a function
of one argument to be called on the ``DataFrame``. A *copy* of the original
DataFrame is returned, with the new values inserted.

.. versionchanged:: 0.23.0

Starting with Python 3.6 the order of ``**kwargs`` is preserved. This allows
for *dependent* assignment, where an expression later in ``**kwargs`` can refer
to a column created earlier in the same :meth:`~DataFrame.assign`.

.. ipython:: python

   dfa = pd.DataFrame({"A": [1, 2, 3],
                       "B": [4, 5, 6]})
   dfa.assign(C=lambda x: x['A'] + x['B'],
              D=lambda x: x['A'] + x['C'])

In the second expression, ``x['C']`` will refer to the newly created column,
that's equal to ``dfa['A'] + dfa['B']``.

To write code compatible with all versions of Python, split the assignment in two.

.. ipython:: python

   dependent = pd.DataFrame({"A": [1, 1, 1]})
   (dependent.assign(A=lambda x: x['A'] + 1)
             .assign(B=lambda x: x['A'] + 2))

.. warning::

   Dependent assignment maybe subtly change the behavior of your code between
   Python 3.6 and older versions of Python.

   If you wish write code that supports versions of python before and after 3.6,
   you'll need to take care when passing ``assign`` expressions that

   * Updating an existing column
   * Referring to the newly updated column in the same ``assign``

   For example, we'll update column "A" and then refer to it when creating "B".

   .. code-block:: python

      >>> dependent = pd.DataFrame({"A": [1, 1, 1]})
      >>> dependent.assign(A=lambda x: x["A"] + 1,
                           B=lambda x: x["A"] + 2)

   For Python 3.5 and earlier the expression creating ``B`` refers to the
   "old" value of ``A``, ``[1, 1, 1]``. The output is then

   .. code-block:: python

         A  B
      0  2  3
      1  2  3
      2  2  3

   For Python 3.6 and later, the expression creating ``A`` refers to the
   "new" value of ``A``, ``[2, 2, 2]``, which results in

   .. code-block:: python

         A  B
      0  2  4
      1  2  4
      2  2  4



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

For a more exhaustive treatment of sophisticated label-based indexing and
slicing, see the :ref:`section on indexing <indexing>`. We will address the
fundamentals of reindexing / conforming to new sets of labels in the
:ref:`section on reindexing <basics.reindexing>`.

Data alignment and arithmetic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data alignment between DataFrame objects automatically align on **both the
columns and the index (row labels)**. Again, the resulting object will have the
union of the column and row labels.

.. ipython:: python

    df = pd.DataFrame(np.random.randn(10, 4), columns=['A', 'B', 'C', 'D'])
    df2 = pd.DataFrame(np.random.randn(7, 3), columns=['A', 'B', 'C'])
    df + df2

When doing an operation between DataFrame and Series, the default behavior is
to align the Series **index** on the DataFrame **columns**, thus `broadcasting
<http://docs.scipy.org/doc/numpy/user/basics.broadcasting.html>`__
row-wise. For example:

.. ipython:: python

   df - df.iloc[0]

In the special case of working with time series data, and the DataFrame index
also contains dates, the broadcasting will be column-wise:

.. ipython:: python
   :okwarning:

   index = pd.date_range('1/1/2000', periods=8)
   df = pd.DataFrame(np.random.randn(8, 3), index=index, columns=list('ABC'))
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

   df1 = pd.DataFrame({'a' : [1, 0, 1], 'b' : [0, 1, 1] }, dtype=bool)
   df2 = pd.DataFrame({'a' : [0, 1, 1], 'b' : [1, 1, 0] }, dtype=bool)
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

   s1 = pd.Series(np.arange(5,10))
   s1.dot(s1)

DataFrame is not intended to be a drop-in replacement for ndarray as its
indexing semantics are quite different in places from a matrix.

Console display
~~~~~~~~~~~~~~~

Very large DataFrames will be truncated to display them in the console.
You can also get a summary using :meth:`~pandas.DataFrame.info`.
(Here I am reading a CSV version of the **baseball** dataset from the **plyr**
R package):

.. ipython:: python
   :suppress:

   # force a summary to be printed
   pd.set_option('display.max_rows', 5)

.. ipython:: python

   baseball = pd.read_csv('data/baseball.csv')
   print(baseball)
   baseball.info()

.. ipython:: python
   :suppress:
   :okwarning:

   # restore GlobalPrintConfig
   pd.reset_option('^display\.')

However, using ``to_string`` will return a string representation of the
DataFrame in tabular form, though it won't always fit the console width:

.. ipython:: python

   print(baseball.iloc[-20:, :12].to_string())

Wide DataFrames will be printed across multiple rows by
default:

.. ipython:: python

   pd.DataFrame(np.random.randn(3, 12))

You can change how much to print on a single row by setting the ``display.width``
option:

.. ipython:: python

   pd.set_option('display.width', 40) # default is 80

   pd.DataFrame(np.random.randn(3, 12))

You can adjust the max width of the individual columns by setting ``display.max_colwidth``

.. ipython:: python

   datafile={'filename': ['filename_01','filename_02'],
             'path': ["media/user_name/storage/folder_01/filename_01",
                      "media/user_name/storage/folder_02/filename_02"]}

   pd.set_option('display.max_colwidth',30)
   pd.DataFrame(datafile)

   pd.set_option('display.max_colwidth',100)
   pd.DataFrame(datafile)

.. ipython:: python
   :suppress:

   pd.reset_option('display.width')
   pd.reset_option('display.max_colwidth')

You can also disable this feature via the ``expand_frame_repr`` option.
This will print the table in one block.

DataFrame column attribute access and IPython completion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a DataFrame column label is a valid Python variable name, the column can be
accessed like an attribute:

.. ipython:: python

   df = pd.DataFrame({'foo1' : np.random.randn(5),
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

.. warning::

    In 0.20.0, ``Panel`` is deprecated and will be removed in
    a future version. See the section :ref:`Deprecate Panel <dsintro.deprecate_panel>`.

Panel is a somewhat less-used, but still important container for 3-dimensional
data. The term `panel data <http://en.wikipedia.org/wiki/Panel_data>`__ is
derived from econometrics and is partially responsible for the name pandas:
pan(el)-da(ta)-s. The names for the 3 axes are intended to give some semantic
meaning to describing operations involving panel data and, in particular,
econometric analysis of panel data. However, for the strict purposes of slicing
and dicing a collection of DataFrame objects, you may find the axis names
slightly arbitrary:

* **items**: axis 0, each item corresponds to a DataFrame contained inside
* **major_axis**: axis 1, it is the **index** (rows) of each of the
  DataFrames
* **minor_axis**: axis 2, it is the **columns** of each of the DataFrames

Construction of Panels works about like you would expect:

From 3D ndarray with optional axis labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python
   :okwarning:

   wp = pd.Panel(np.random.randn(2, 5, 4), items=['Item1', 'Item2'],
                 major_axis=pd.date_range('1/1/2000', periods=5),
                 minor_axis=['A', 'B', 'C', 'D'])
   wp


From dict of DataFrame objects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. ipython:: python
   :okwarning:

   data = {'Item1' : pd.DataFrame(np.random.randn(4, 3)),
           'Item2' : pd.DataFrame(np.random.randn(4, 2))}
   pd.Panel(data)

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
   :okwarning:

   pd.Panel.from_dict(data, orient='minor')

Orient is especially useful for mixed-type DataFrames. If you pass a dict of
DataFrame objects with mixed-type columns, all of the data will get upcasted to
``dtype=object`` unless you pass ``orient='minor'``:

.. ipython:: python
   :okwarning:

   df = pd.DataFrame({'a': ['foo', 'bar', 'baz'],
                      'b': np.random.randn(3)})
   df
   data = {'item1': df, 'item2': df}
   panel = pd.Panel.from_dict(data, orient='minor')
   panel['a']
   panel['b']
   panel['b'].dtypes

.. note::

   Panel, being less commonly used than Series and DataFrame,
   has been slightly neglected feature-wise. A number of methods and options
   available in DataFrame are not available in Panel.

.. _dsintro.to_panel:

From DataFrame using ``to_panel`` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``to_panel`` converts a DataFrame with a two-level index to a Panel.

.. ipython:: python
   :okwarning:

   midx = pd.MultiIndex(levels=[['one', 'two'], ['x','y']], labels=[[1,1,0,0],[1,0,1,0]])
   df = pd.DataFrame({'A' : [1, 2, 3, 4], 'B': [5, 6, 7, 8]}, index=midx)
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
DataFrame, if the item is a valid Python identifier, you can access it as an
attribute and tab-complete it in IPython.

Transposing
~~~~~~~~~~~

A Panel can be rearranged using its ``transpose`` method (which does not make a
copy by default unless the data are heterogeneous):

.. ipython:: python
   :okwarning:

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

Another way to change the dimensionality of an object is to ``squeeze`` a 1-len
object, similar to ``wp['Item1']``.

.. ipython:: python
   :okwarning:

   wp.reindex(items=['Item1']).squeeze()
   wp.reindex(items=['Item1'], minor=['B']).squeeze()


Conversion to DataFrame
~~~~~~~~~~~~~~~~~~~~~~~

A Panel can be represented in 2D form as a hierarchically indexed
DataFrame. See the section :ref:`hierarchical indexing <advanced.hierarchical>`
for more on this. To convert a Panel to a DataFrame, use the ``to_frame``
method:

.. ipython:: python
   :okwarning:

   panel = pd.Panel(np.random.randn(3, 5, 4), items=['one', 'two', 'three'],
                    major_axis=pd.date_range('1/1/2000', periods=5),
                    minor_axis=['a', 'b', 'c', 'd'])
   panel.to_frame()


.. _dsintro.deprecate_panel:

Deprecate Panel
---------------

Over the last few years, pandas has increased in both breadth and depth, with new features,
datatype support, and manipulation routines. As a result, supporting efficient indexing and functional
routines for ``Series``, ``DataFrame`` and ``Panel`` has contributed to an increasingly fragmented and
difficult-to-understand code base.

The 3-D structure of a ``Panel`` is much less common for many types of data analysis,
than the 1-D of the ``Series`` or the 2-D of the ``DataFrame``. Going forward it makes sense for
pandas to focus on these areas exclusively.

Oftentimes, one can simply use a MultiIndex ``DataFrame`` for easily working with higher dimensional data.

In addition, the ``xarray`` package was built from the ground up, specifically in order to
support the multi-dimensional analysis that is one of ``Panel`` s main use cases.
`Here is a link to the xarray panel-transition documentation <http://xarray.pydata.org/en/stable/pandas.html#panel-transition>`__.

.. ipython:: python
   :okwarning:

   p = tm.makePanel()
   p

Convert to a MultiIndex DataFrame.

.. ipython:: python
   :okwarning:

   p.to_frame()

Alternatively, one can convert to an xarray ``DataArray``.

.. ipython:: python
   :okwarning:

   p.to_xarray()

You can see the full-documentation for the `xarray package <http://xarray.pydata.org/en/stable/>`__.
