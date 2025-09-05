.. _dsintro:

{{ header }}

************************
Intro to data structures
************************

We'll start with a quick, non-comprehensive overview of the fundamental data
structures in pandas to get you started. The fundamental behavior about data
types, indexing, axis labeling, and alignment apply across all of the
objects. To get started, import NumPy and load pandas into your namespace:

.. ipython:: python

   import numpy as np
   import pandas as pd

Fundamentally, **data alignment is intrinsic**. The link
between labels and data will not be broken unless done so explicitly by you.

We'll give a brief intro to the data structures, then consider all of the broad
categories of functionality and methods in separate sections.

.. _basics.series:

Series
------

:class:`Series` is a one-dimensional labeled array capable of holding any data
type (integers, strings, floating point numbers, Python objects, etc.). The axis
labels are collectively referred to as the **index**. The basic method to create a :class:`Series` is to call:

.. code-block:: python

   s = pd.Series(data, index=index)

Here, ``data`` can be many different things:

* a Python dict
* an ndarray
* a scalar value (like 5)

The passed **index** is a list of axis labels. The constructor's behavior
depends on **data**'s type:

**From ndarray**

If ``data`` is an ndarray, **index** must be the same length as **data**. If no
index is passed, one will be created having values ``[0, ..., len(data) - 1]``.

.. ipython:: python

   s = pd.Series(np.random.randn(5), index=["a", "b", "c", "d", "e"])
   s
   s.index

   pd.Series(np.random.randn(5))

.. note::

    pandas supports non-unique index values. If an operation
    that does not support duplicate index values is attempted, an exception
    will be raised at that time.

**From dict**

:class:`Series` can be instantiated from dicts:

.. ipython:: python

   d = {"b": 1, "a": 0, "c": 2}
   pd.Series(d)

If an index is passed, the values in data corresponding to the labels in the
index will be pulled out.

.. ipython:: python

   d = {"a": 0.0, "b": 1.0, "c": 2.0}
   pd.Series(d)
   pd.Series(d, index=["b", "c", "d", "a"])

.. note::

    NaN (not a number) is the standard missing data marker used in pandas.

**From scalar value**

If ``data`` is a scalar value, the value will be repeated to match
the length of **index**.  If the **index** is not provided, it defaults
to ``RangeIndex(1)``.

.. ipython:: python

   pd.Series(5.0, index=["a", "b", "c", "d", "e"])

Series is ndarray-like
~~~~~~~~~~~~~~~~~~~~~~

:class:`Series` acts very similarly to a :class:`numpy.ndarray` and is a valid argument to most NumPy functions.
However, operations such as slicing will also slice the index.

.. ipython:: python

    s.iloc[0]
    s.iloc[:3]
    s[s > s.median()]
    s.iloc[[4, 3, 1]]
    np.exp(s)

.. note::

   We will address array-based indexing like ``s.iloc[[4, 3, 1]]``
   in the :ref:`section on indexing <indexing>`.

Like a NumPy array, a pandas :class:`Series` has a single :attr:`~Series.dtype`.

.. ipython:: python

   s.dtype

This is often a NumPy dtype. However, pandas and 3rd-party libraries
extend NumPy's type system in a few places, in which case the dtype would
be an :class:`~pandas.api.extensions.ExtensionDtype`. Some examples within
pandas are :ref:`categorical` and :ref:`integer_na`. See :ref:`basics.dtypes`
for more.

If you need the actual array backing a :class:`Series`, use :attr:`Series.array`.

.. ipython:: python

   s.array

Accessing the array can be useful when you need to do some operation without the
index (to disable :ref:`automatic alignment <dsintro.alignment>`, for example).

:attr:`Series.array` will always be an :class:`~pandas.api.extensions.ExtensionArray`.
Briefly, an ExtensionArray is a thin wrapper around one or more *concrete* arrays like a
:class:`numpy.ndarray`. pandas knows how to take an :class:`~pandas.api.extensions.ExtensionArray` and
store it in a :class:`Series` or a column of a :class:`DataFrame`.
See :ref:`basics.dtypes` for more.

While :class:`Series` is ndarray-like, if you need an *actual* ndarray, then use
:meth:`Series.to_numpy`.

.. ipython:: python

   s.to_numpy()

Even if the :class:`Series` is backed by a :class:`~pandas.api.extensions.ExtensionArray`,
:meth:`Series.to_numpy` will return a NumPy ndarray.

Series is dict-like
~~~~~~~~~~~~~~~~~~~

A :class:`Series` is also like a fixed-size dict in that you can get and set values by index
label:

.. ipython:: python

    s["a"]
    s["e"] = 12.0
    s
    "e" in s
    "f" in s

If a label is not contained in the index, an exception is raised:

.. ipython:: python
    :okexcept:

    s["f"]

Using the :meth:`Series.get` method, a missing label will return None or specified default:

.. ipython:: python

   s.get("f")

   s.get("f", np.nan)

These labels can also be accessed by :ref:`attribute<indexing.attribute_access>`.

Vectorized operations and label alignment with Series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When working with raw NumPy arrays, looping through value-by-value is usually
not necessary. The same is true when working with :class:`Series` in pandas.
:class:`Series` can also be passed into most NumPy methods expecting an ndarray.

.. ipython:: python

    s + s
    s * 2
    np.exp(s)

A key difference between :class:`Series` and ndarray is that operations between :class:`Series`
automatically align the data based on label. Thus, you can write computations
without giving consideration to whether the :class:`Series` involved have the same
labels.

.. ipython:: python

    s.iloc[1:] + s.iloc[:-1]

The result of an operation between unaligned :class:`Series` will have the **union** of
the indexes involved. If a label is not found in one :class:`Series` or the other, the
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

:class:`Series` also has a ``name`` attribute:

.. ipython:: python

   s = pd.Series(np.random.randn(5), name="something")
   s
   s.name

The :class:`Series` ``name`` can be assigned automatically in many cases, in particular,
when selecting a single column from a :class:`DataFrame`, the ``name`` will be assigned
the column label.

You can rename a :class:`Series` with the :meth:`pandas.Series.rename` method.

.. ipython:: python

   s2 = s.rename("different")
   s2.name

Note that ``s`` and ``s2`` refer to different objects.

.. _basics.dataframe:

DataFrame
---------

:class:`DataFrame` is a 2-dimensional labeled data structure with columns of
potentially different types. You can think of it like a spreadsheet or SQL
table, or a dict of Series objects. It is generally the most commonly used
pandas object. Like Series, DataFrame accepts many different kinds of input:

* Dict of 1D ndarrays, lists, dicts, or :class:`Series`
* 2-D numpy.ndarray
* `Structured or record
  <https://numpy.org/doc/stable/user/basics.rec.html>`__ ndarray
* A :class:`Series`
* Another :class:`DataFrame`

Along with the data, you can optionally pass **index** (row labels) and
**columns** (column labels) arguments. If you pass an index and / or columns,
you are guaranteeing the index and / or columns of the resulting
DataFrame. Thus, a dict of Series plus a specific index will discard all data
not matching up to the passed index.

If axis labels are not passed, they will be constructed from the input data
based on common sense rules.

From dict of Series or dicts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The resulting **index** will be the **union** of the indexes of the various
Series. If there are any nested dicts, these will first be converted to
Series. If no columns are passed, the columns will be the ordered list of dict
keys.

.. ipython:: python

    d = {
        "one": pd.Series([1.0, 2.0, 3.0], index=["a", "b", "c"]),
        "two": pd.Series([1.0, 2.0, 3.0, 4.0], index=["a", "b", "c", "d"]),
    }
    df = pd.DataFrame(d)
    df

    pd.DataFrame(d, index=["d", "b", "a"])
    pd.DataFrame(d, index=["d", "b", "a"], columns=["two", "three"])

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

All ndarrays must share the same length. If an index is passed, it must
also be the same length as the arrays. If no index is passed, the
result will be ``range(n)``, where ``n`` is the array length.

.. ipython:: python

   d = {"one": [1.0, 2.0, 3.0, 4.0], "two": [4.0, 3.0, 2.0, 1.0]}
   pd.DataFrame(d)
   pd.DataFrame(d, index=["a", "b", "c", "d"])

From structured or record array
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This case is handled identically to a dict of arrays.

.. ipython:: python

   data = np.zeros((2,), dtype=[("A", "i4"), ("B", "f4"), ("C", "S10")])
   data[:] = [(1, 2.0, "Hello"), (2, 3.0, "World")]

   pd.DataFrame(data)
   pd.DataFrame(data, index=["first", "second"])
   pd.DataFrame(data, columns=["C", "A", "B"])

.. note::

    DataFrame is not intended to work exactly like a 2-dimensional NumPy
    ndarray.

.. _basics.dataframe.from_list_of_dicts:

From a list of dicts
~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

   data2 = [{"a": 1, "b": 2}, {"a": 5, "b": 10, "c": 20}]
   pd.DataFrame(data2)
   pd.DataFrame(data2, index=["first", "second"])
   pd.DataFrame(data2, columns=["a", "b"])

.. _basics.dataframe.from_dict_of_tuples:

From a dict of tuples
~~~~~~~~~~~~~~~~~~~~~

You can automatically create a MultiIndexed frame by passing a tuples
dictionary.

.. ipython:: python

   pd.DataFrame(
       {
           ("a", "b"): {("A", "B"): 1, ("A", "C"): 2},
           ("a", "a"): {("A", "C"): 3, ("A", "B"): 4},
           ("a", "c"): {("A", "B"): 5, ("A", "C"): 6},
           ("b", "a"): {("A", "C"): 7, ("A", "B"): 8},
           ("b", "b"): {("A", "D"): 9, ("A", "B"): 10},
       }
   )

.. _basics.dataframe.from_series:

From a Series
~~~~~~~~~~~~~

The result will be a DataFrame with the same index as the input Series, and
with one column whose name is the original name of the Series (only if no other
column name provided).

.. ipython:: python

   ser = pd.Series(range(3), index=list("abc"), name="ser")
   pd.DataFrame(ser)

.. _basics.dataframe.from_list_namedtuples:

From a list of namedtuples
~~~~~~~~~~~~~~~~~~~~~~~~~~

The field names of the first ``namedtuple`` in the list determine the columns
of the :class:`DataFrame`. The remaining namedtuples (or tuples) are simply unpacked
and their values are fed into the rows of the :class:`DataFrame`. If any of those
tuples is shorter than the first ``namedtuple`` then the later columns in the
corresponding row are marked as missing values. If any are longer than the
first ``namedtuple``, a ``ValueError`` is raised.

.. ipython:: python

    from collections import namedtuple

    Point = namedtuple("Point", "x y")

    pd.DataFrame([Point(0, 0), Point(0, 3), (2, 3)])

    Point3D = namedtuple("Point3D", "x y z")

    pd.DataFrame([Point3D(0, 0, 0), Point3D(0, 3, 5), Point(2, 3)])


.. _basics.dataframe.from_list_dataclasses:

From a list of dataclasses
~~~~~~~~~~~~~~~~~~~~~~~~~~

Data Classes as introduced in `PEP557 <https://www.python.org/dev/peps/pep-0557>`__,
can be passed into the DataFrame constructor.
Passing a list of dataclasses is equivalent to passing a list of dictionaries.

Please be aware, that all values in the list should be dataclasses, mixing
types in the list would result in a ``TypeError``.

.. ipython:: python

    from dataclasses import make_dataclass

    Point = make_dataclass("Point", [("x", int), ("y", int)])

    pd.DataFrame([Point(0, 0), Point(0, 3), Point(2, 3)])

**Missing data**

To construct a DataFrame with missing data, we use ``np.nan`` to
represent missing values. Alternatively, you may pass a ``numpy.MaskedArray``
as the data argument to the DataFrame constructor, and its masked entries will
be considered missing. See :ref:`Missing data <missing_data>` for more.

Alternate constructors
~~~~~~~~~~~~~~~~~~~~~~

.. _basics.dataframe.from_dict:

**DataFrame.from_dict**

:meth:`DataFrame.from_dict` takes a dict of dicts or a dict of array-like sequences
and returns a DataFrame. It operates like the :class:`DataFrame` constructor except
for the ``orient`` parameter which is ``'columns'`` by default, but which can be
set to ``'index'`` in order to use the dict keys as row labels.


.. ipython:: python

   pd.DataFrame.from_dict(dict([("A", [1, 2, 3]), ("B", [4, 5, 6])]))

If you pass ``orient='index'``, the keys will be the row labels. In this
case, you can also pass the desired column names:

.. ipython:: python

   pd.DataFrame.from_dict(
       dict([("A", [1, 2, 3]), ("B", [4, 5, 6])]),
       orient="index",
       columns=["one", "two", "three"],
   )

.. _basics.dataframe.from_records:

**DataFrame.from_records**

:meth:`DataFrame.from_records` takes a list of tuples or an ndarray with structured
dtype. It works analogously to the normal :class:`DataFrame` constructor, except that
the resulting DataFrame index may be a specific field of the structured
dtype.

.. ipython:: python

   data
   pd.DataFrame.from_records(data, index="C")

.. _basics.dataframe.sel_add_del:

Column selection, addition, deletion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can treat a :class:`DataFrame` semantically like a dict of like-indexed :class:`Series`
objects. Getting, setting, and deleting columns works with the same syntax as
the analogous dict operations:

.. ipython:: python

   df["one"]
   df["three"] = df["one"] * df["two"]
   df["flag"] = df["one"] > 2
   df

Columns can be deleted or popped like with a dict:

.. ipython:: python

   del df["two"]
   three = df.pop("three")
   df

When inserting a scalar value, it will naturally be propagated to fill the
column:

.. ipython:: python

   df["foo"] = "bar"
   df

When inserting a :class:`Series` that does not have the same index as the :class:`DataFrame`, it
will be conformed to the DataFrame's index:

.. ipython:: python

   df["one_trunc"] = df["one"][:2]
   df

You can insert raw ndarrays but their length must match the length of the
DataFrame's index.

By default, columns get inserted at the end. :meth:`DataFrame.insert`
inserts at a particular location in the columns:

.. ipython:: python

   df.insert(1, "bar", df["one"])
   df

.. _dsintro.chained_assignment:

Assigning new columns in method chains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inspired by `dplyr's
<https://dplyr.tidyverse.org/reference/mutate.html>`__
``mutate`` verb, DataFrame has an :meth:`~pandas.DataFrame.assign`
method that allows you to easily create new columns that are potentially
derived from existing columns.

.. ipython:: python

   iris = pd.read_csv("data/iris.data")
   iris.head()
   iris.assign(sepal_ratio=iris["SepalWidth"] / iris["SepalLength"]).head()

In the example above, we inserted a precomputed value. We can also pass in
a function of one argument to be evaluated on the DataFrame being assigned to.

.. ipython:: python

   iris.assign(sepal_ratio=lambda x: (x["SepalWidth"] / x["SepalLength"])).head()

or, using :meth:`pandas.col`:

.. ipython:: python

   iris.assign(sepal_ratio=pd.col("SepalWidth") / pd.col("SepalLength")).head()

:meth:`~pandas.DataFrame.assign` **always** returns a copy of the data, leaving the original
DataFrame untouched.

Passing a callable, as opposed to an actual value to be inserted, is
useful when you don't have a reference to the DataFrame at hand. This is
common when using :meth:`~pandas.DataFrame.assign` in a chain of operations. For example,
we can limit the DataFrame to just those observations with a Sepal Length
greater than 5, calculate the ratio, and plot:

.. ipython:: python

   @savefig basics_assign.png
   (
       iris.query("SepalLength > 5")
       .assign(
           SepalRatio=lambda x: x.SepalWidth / x.SepalLength,
           PetalRatio=lambda x: x.PetalWidth / x.PetalLength,
       )
       .plot(kind="scatter", x="SepalRatio", y="PetalRatio")
   )

Since a function is passed in, the function is computed on the DataFrame
being assigned to. Importantly, this is the DataFrame that's been filtered
to those rows with sepal length greater than 5. The filtering happens first,
and then the ratio calculations. This is an example where we didn't
have a reference to the *filtered* DataFrame available.

The function signature for :meth:`~pandas.DataFrame.assign` is simply ``**kwargs``. The keys
are the column names for the new fields, and the values are either a value
to be inserted (for example, a :class:`Series` or NumPy array), or a function
of one argument to be called on the :class:`DataFrame`. A *copy* of the original
:class:`DataFrame` is returned, with the new values inserted.

The order of ``**kwargs`` is preserved. This allows
for *dependent* assignment, where an expression later in ``**kwargs`` can refer
to a column created earlier in the same :meth:`~DataFrame.assign`.

.. ipython:: python

   dfa = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
   dfa.assign(C=lambda x: x["A"] + x["B"], D=lambda x: x["A"] + x["C"])

In the second expression, ``x['C']`` will refer to the newly created column,
that's equal to ``dfa['A'] + dfa['B']``.


Indexing / selection
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

Row selection, for example, returns a :class:`Series` whose index is the columns of the
:class:`DataFrame`:

.. ipython:: python

   df.loc["b"]
   df.iloc[2]

For a more exhaustive treatment of sophisticated label-based indexing and
slicing, see the :ref:`section on indexing <indexing>`. We will address the
fundamentals of reindexing / conforming to new sets of labels in the
:ref:`section on reindexing <basics.reindexing>`.

.. _dsintro.alignment:

Data alignment and arithmetic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Data alignment between :class:`DataFrame` objects automatically align on **both the
columns and the index (row labels)**. Again, the resulting object will have the
union of the column and row labels.

.. ipython:: python

    df = pd.DataFrame(np.random.randn(10, 4), columns=["A", "B", "C", "D"])
    df2 = pd.DataFrame(np.random.randn(7, 3), columns=["A", "B", "C"])
    df + df2

When doing an operation between :class:`DataFrame` and :class:`Series`, the default behavior is
to align the :class:`Series` **index** on the :class:`DataFrame` **columns**, thus `broadcasting
<https://numpy.org/doc/stable/user/basics.broadcasting.html>`__
row-wise. For example:

.. ipython:: python

   df - df.iloc[0]

For explicit control over the matching and broadcasting behavior, see the
section on :ref:`flexible binary operations <basics.binop>`.

Arithmetic operations with scalars operate element-wise:

.. ipython:: python

   df * 5 + 2
   1 / df
   df ** 4

.. _dsintro.boolean:

Boolean operators operate element-wise as well:

.. ipython:: python

   df1 = pd.DataFrame({"a": [1, 0, 1], "b": [0, 1, 1]}, dtype=bool)
   df2 = pd.DataFrame({"a": [0, 1, 1], "b": [1, 1, 0]}, dtype=bool)
   df1 & df2
   df1 | df2
   df1 ^ df2
   -df1

Transposing
~~~~~~~~~~~

To transpose, access the ``T`` attribute or :meth:`DataFrame.transpose`,
similar to an ndarray:

.. ipython:: python

   # only show the first 5 rows
   df[:5].T

.. _dsintro.numpy_interop:

DataFrame interoperability with NumPy functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most NumPy functions can be called directly on :class:`Series` and :class:`DataFrame`.

.. ipython:: python

   np.exp(df)
   np.asarray(df)

:class:`DataFrame` is not intended to be a drop-in replacement for ndarray as its
indexing semantics and data model are quite different in places from an n-dimensional
array.

:class:`Series` implements ``__array_ufunc__``, which allows it to work with NumPy's
`universal functions <https://numpy.org/doc/stable/reference/ufuncs.html>`_.

The ufunc is applied to the underlying array in a :class:`Series`.

.. ipython:: python

   ser = pd.Series([1, 2, 3, 4])
   np.exp(ser)

When multiple :class:`Series` are passed to a ufunc, they are aligned before
performing the operation.

Like other parts of the library, pandas will automatically align labeled inputs
as part of a ufunc with multiple inputs. For example, using :meth:`numpy.remainder`
on two :class:`Series` with differently ordered labels will align before the operation.

.. ipython:: python

   ser1 = pd.Series([1, 2, 3], index=["a", "b", "c"])
   ser2 = pd.Series([1, 3, 5], index=["b", "a", "c"])
   ser1
   ser2
   np.remainder(ser1, ser2)

As usual, the union of the two indices is taken, and non-overlapping values are filled
with missing values.

.. ipython:: python

   ser3 = pd.Series([2, 4, 6], index=["b", "c", "d"])
   ser3
   np.remainder(ser1, ser3)

When a binary ufunc is applied to a :class:`Series` and :class:`Index`, the :class:`Series`
implementation takes precedence and a :class:`Series` is returned.

.. ipython:: python

   ser = pd.Series([1, 2, 3])
   idx = pd.Index([4, 5, 6])

   np.maximum(ser, idx)

NumPy ufuncs are safe to apply to :class:`Series` backed by non-ndarray arrays,
for example :class:`arrays.SparseArray` (see :ref:`sparse.calculation`). If possible,
the ufunc is applied without converting the underlying data to an ndarray.

Console display
~~~~~~~~~~~~~~~

A very large :class:`DataFrame` will be truncated to display them in the console.
You can also get a summary using :meth:`~pandas.DataFrame.info`.
(The **baseball** dataset is from the **plyr** R package):

.. ipython:: python
   :suppress:

   # force a summary to be printed
   pd.set_option("display.max_rows", 5)

.. ipython:: python

   baseball = pd.read_csv("data/baseball.csv")
   print(baseball)
   baseball.info()

.. ipython:: python
   :suppress:
   :okwarning:

   # restore GlobalPrintConfig
   pd.reset_option(r"^display\.")

However, using :meth:`DataFrame.to_string` will return a string representation of the
:class:`DataFrame` in tabular form, though it won't always fit the console width:

.. ipython:: python

   print(baseball.iloc[-20:, :12].to_string())

Wide DataFrames will be printed across multiple rows by
default:

.. ipython:: python

   pd.DataFrame(np.random.randn(3, 12))

You can change how much to print on a single row by setting the ``display.width``
option:

.. ipython:: python

   pd.set_option("display.width", 40)  # default is 80

   pd.DataFrame(np.random.randn(3, 12))

You can adjust the max width of the individual columns by setting ``display.max_colwidth``

.. ipython:: python

   datafile = {
       "filename": ["filename_01", "filename_02"],
       "path": [
           "media/user_name/storage/folder_01/filename_01",
           "media/user_name/storage/folder_02/filename_02",
       ],
   }

   pd.set_option("display.max_colwidth", 30)
   pd.DataFrame(datafile)

   pd.set_option("display.max_colwidth", 100)
   pd.DataFrame(datafile)

.. ipython:: python
   :suppress:

   pd.reset_option("display.width")
   pd.reset_option("display.max_colwidth")

You can also disable this feature via the ``expand_frame_repr`` option.
This will print the table in one block.

DataFrame column attribute access and IPython completion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a :class:`DataFrame` column label is a valid Python variable name, the column can be
accessed like an attribute:

.. ipython:: python

   df = pd.DataFrame({"foo1": np.random.randn(5), "foo2": np.random.randn(5)})
   df
   df.foo1

The columns are also connected to the `IPython <https://ipython.org>`__
completion mechanism so they can be tab-completed:

.. code-block:: ipython

    In [5]: df.foo<TAB>  # noqa: E225, E999
    df.foo1  df.foo2
