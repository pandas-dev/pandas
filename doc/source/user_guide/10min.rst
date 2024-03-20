.. _10min:

{{ header }}

********************
10 minutes to pandas
********************

This is a short introduction to pandas, geared mainly for new users.
You can see more complex recipes in the :ref:`Cookbook<cookbook>`.

Customarily, we import as follows:

.. ipython:: python

   import numpy as np
   import pandas as pd

Basic data structures in pandas
-------------------------------

pandas provides two types of classes for handling data:

1. :class:`Series`: a one-dimensional labeled array holding data of any type
    such as integers, strings, Python objects etc.
2. :class:`DataFrame`: a two-dimensional data structure that holds data like
   a two-dimension array or a table with rows and columns.

Object creation
---------------

See the :ref:`Intro to data structures section <dsintro>`.

Creating a :class:`Series` by passing a list of values, letting pandas create
a default :class:`RangeIndex`.

.. ipython:: python

   s = pd.Series([1, 3, 5, np.nan, 6, 8])
   s

Creating a :class:`DataFrame` by passing a NumPy array with a datetime index using :func:`date_range`
and labeled columns:

.. ipython:: python

   dates = pd.date_range("20130101", periods=6)
   dates
   df = pd.DataFrame(np.random.randn(6, 4), index=dates, columns=list("ABCD"))
   df

Creating a :class:`DataFrame` by passing a dictionary of objects where the keys are the column
labels and the values are the column values.

.. ipython:: python

   df2 = pd.DataFrame(
       {
           "A": 1.0,
           "B": pd.Timestamp("20130102"),
           "C": pd.Series(1, index=list(range(4)), dtype="float32"),
           "D": np.array([3] * 4, dtype="int32"),
           "E": pd.Categorical(["test", "train", "test", "train"]),
           "F": "foo",
       }
   )
   df2

The columns of the resulting :class:`DataFrame` have different
:ref:`dtypes <basics.dtypes>`:

.. ipython:: python

   df2.dtypes

If you're using IPython, tab completion for column names (as well as public
attributes) is automatically enabled. Here's a subset of the attributes that
will be completed:

.. ipython::

   @verbatim
   In [1]: df2.<TAB>  # noqa: E225, E999
   df2.A                  df2.bool
   df2.abs                df2.boxplot
   df2.add                df2.C
   df2.add_prefix         df2.clip
   df2.add_suffix         df2.columns
   df2.align              df2.copy
   df2.all                df2.count
   df2.any                df2.combine
   df2.append             df2.D
   df2.apply              df2.describe
   df2.B                  df2.duplicated
   df2.diff

As you can see, the columns ``A``, ``B``, ``C``, and ``D`` are automatically
tab completed. ``E`` and ``F`` are there as well; the rest of the attributes have been
truncated for brevity.

Viewing data
------------

See the :ref:`Essentially basics functionality section <basics>`.

Use :meth:`DataFrame.head` and :meth:`DataFrame.tail` to view the top and bottom rows of the frame
respectively:

.. ipython:: python

   df.head()
   df.tail(3)

Display the :attr:`DataFrame.index` or :attr:`DataFrame.columns`:

.. ipython:: python

   df.index
   df.columns

Return a NumPy representation of the underlying data with :meth:`DataFrame.to_numpy`
without the index or column labels:

.. ipython:: python

   df.to_numpy()

.. note::

   **NumPy arrays have one dtype for the entire array while pandas DataFrames
   have one dtype per column**. When you call :meth:`DataFrame.to_numpy`, pandas will
   find the NumPy dtype that can hold *all* of the dtypes in the DataFrame.
   If the common data type is ``object``, :meth:`DataFrame.to_numpy` will require
   copying data.

   .. ipython:: python

      df2.dtypes
      df2.to_numpy()

:func:`~DataFrame.describe` shows a quick statistic summary of your data:

.. ipython:: python

   df.describe()

Transposing your data:

.. ipython:: python

   df.T

:meth:`DataFrame.sort_index` sorts by an axis:

.. ipython:: python

   df.sort_index(axis=1, ascending=False)

:meth:`DataFrame.sort_values` sorts by values:

.. ipython:: python

   df.sort_values(by="B")

Selection
---------

.. note::

   While standard Python / NumPy expressions for selecting and setting are
   intuitive and come in handy for interactive work, for production code, we
   recommend the optimized pandas data access methods, :meth:`DataFrame.at`, :meth:`DataFrame.iat`,
   :meth:`DataFrame.loc` and :meth:`DataFrame.iloc`.

See the indexing documentation :ref:`Indexing and Selecting Data <indexing>` and :ref:`MultiIndex / Advanced Indexing <advanced>`.

Getitem (``[]``)
~~~~~~~~~~~~~~~~

For a :class:`DataFrame`, passing a single label selects a columns and
yields a :class:`Series` equivalent to ``df.A``:

.. ipython:: python

   df["A"]

For a :class:`DataFrame`, passing a slice ``:`` selects matching rows:

.. ipython:: python

   df[0:3]
   df["20130102":"20130104"]

Selection by label
~~~~~~~~~~~~~~~~~~

See more in :ref:`Selection by Label <indexing.label>` using :meth:`DataFrame.loc` or :meth:`DataFrame.at`.

Selecting a row matching a label:

.. ipython:: python

   df.loc[dates[0]]

Selecting all rows (``:``) with a select column labels:

.. ipython:: python

   df.loc[:, ["A", "B"]]

For label slicing, both endpoints are *included*:

.. ipython:: python

   df.loc["20130102":"20130104", ["A", "B"]]

Selecting a single row and column label returns a scalar:

.. ipython:: python

   df.loc[dates[0], "A"]

For getting fast access to a scalar (equivalent to the prior method):

.. ipython:: python

   df.at[dates[0], "A"]

Selection by position
~~~~~~~~~~~~~~~~~~~~~

See more in :ref:`Selection by Position <indexing.integer>` using :meth:`DataFrame.iloc` or :meth:`DataFrame.iat`.

Select via the position of the passed integers:

.. ipython:: python

   df.iloc[3]

Integer slices acts similar to NumPy/Python:

.. ipython:: python

   df.iloc[3:5, 0:2]

Lists of integer position locations:

.. ipython:: python

   df.iloc[[1, 2, 4], [0, 2]]

For slicing rows explicitly:

.. ipython:: python

   df.iloc[1:3, :]

For slicing columns explicitly:

.. ipython:: python

   df.iloc[:, 1:3]

For getting a value explicitly:

.. ipython:: python

   df.iloc[1, 1]

For getting fast access to a scalar (equivalent to the prior method):

.. ipython:: python

   df.iat[1, 1]

Boolean indexing
~~~~~~~~~~~~~~~~

Select rows where ``df.A`` is greater than ``0``.

.. ipython:: python

   df[df["A"] > 0]

Selecting values from a :class:`DataFrame` where a boolean condition is met:

.. ipython:: python

   df[df > 0]

Using :func:`~Series.isin` method for filtering:

.. ipython:: python

   df2 = df.copy()
   df2["E"] = ["one", "one", "two", "three", "four", "three"]
   df2
   df2[df2["E"].isin(["two", "four"])]

Setting
~~~~~~~

Setting a new column automatically aligns the data by the indexes:

.. ipython:: python

   s1 = pd.Series([1, 2, 3, 4, 5, 6], index=pd.date_range("20130102", periods=6))
   s1
   df["F"] = s1

Setting values by label:

.. ipython:: python

   df.at[dates[0], "A"] = 0

Setting values by position:

.. ipython:: python

   df.iat[0, 1] = 0

Setting by assigning with a NumPy array:

.. ipython:: python
   :okwarning:

   df.loc[:, "D"] = np.array([5] * len(df))

The result of the prior setting operations:

.. ipython:: python

   df

A ``where`` operation with setting:

.. ipython:: python

   df2 = df.copy()
   df2[df2 > 0] = -df2
   df2


Missing data
------------

For NumPy data types, ``np.nan`` represents missing data. It is by
default not included in computations. See the :ref:`Missing Data section
<missing_data>`.

Reindexing allows you to change/add/delete the index on a specified axis. This
returns a copy of the data:

.. ipython:: python

   df1 = df.reindex(index=dates[0:4], columns=list(df.columns) + ["E"])
   df1.loc[dates[0] : dates[1], "E"] = 1
   df1

:meth:`DataFrame.dropna` drops any rows that have missing data:

.. ipython:: python

   df1.dropna(how="any")

:meth:`DataFrame.fillna` fills missing data:

.. ipython:: python

   df1.fillna(value=5)

:func:`isna` gets the boolean mask where values are ``nan``:

.. ipython:: python

   pd.isna(df1)


Operations
----------

See the :ref:`Basic section on Binary Ops <basics.binop>`.

Stats
~~~~~

Operations in general *exclude* missing data.

Calculate the mean value for each column:

.. ipython:: python

   df.mean()

Calculate the mean value for each row:

.. ipython:: python

   df.mean(axis=1)

Operating with another :class:`Series` or :class:`DataFrame` with a different index or column
will align the result with the union of the index or column labels. In addition, pandas
automatically broadcasts along the specified dimension and will fill unaligned labels with ``np.nan``.

.. ipython:: python

   s = pd.Series([1, 3, 5, np.nan, 6, 8], index=dates).shift(2)
   s
   df.sub(s, axis="index")


User defined functions
~~~~~~~~~~~~~~~~~~~~~~

:meth:`DataFrame.agg` and :meth:`DataFrame.transform` applies a user defined function
that reduces or broadcasts its result respectively.

.. ipython:: python

   df.agg(lambda x: np.mean(x) * 5.6)
   df.transform(lambda x: x * 101.2)

Value Counts
~~~~~~~~~~~~~

See more at :ref:`Histogramming and Discretization <basics.discretization>`.

.. ipython:: python

   s = pd.Series(np.random.randint(0, 7, size=10))
   s
   s.value_counts()

String Methods
~~~~~~~~~~~~~~

:class:`Series` is equipped with a set of string processing methods in the ``str``
attribute that make it easy to operate on each element of the array, as in the
code snippet below. See more at :ref:`Vectorized String Methods
<text.string_methods>`.

.. ipython:: python

   s = pd.Series(["A", "B", "C", "Aaba", "Baca", np.nan, "CABA", "dog", "cat"])
   s.str.lower()

Merge
-----

Concat
~~~~~~

pandas provides various facilities for easily combining together :class:`Series` and
:class:`DataFrame` objects with various kinds of set logic for the indexes
and relational algebra functionality in the case of join / merge-type
operations.

See the :ref:`Merging section <merging>`.

Concatenating pandas objects together row-wise with :func:`concat`:

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 4))
   df

   # break it into pieces
   pieces = [df[:3], df[3:7], df[7:]]

   pd.concat(pieces)

.. note::

   Adding a column to a :class:`DataFrame` is relatively fast. However, adding
   a row requires a copy, and may be expensive. We recommend passing a
   pre-built list of records to the :class:`DataFrame` constructor instead
   of building a :class:`DataFrame` by iteratively appending records to it.

Join
~~~~

:func:`merge` enables SQL style join types along specific columns. See the :ref:`Database style joining <merging.join>` section.

.. ipython:: python

   left = pd.DataFrame({"key": ["foo", "foo"], "lval": [1, 2]})
   right = pd.DataFrame({"key": ["foo", "foo"], "rval": [4, 5]})
   left
   right
   pd.merge(left, right, on="key")

:func:`merge` on unique keys:

.. ipython:: python

   left = pd.DataFrame({"key": ["foo", "bar"], "lval": [1, 2]})
   right = pd.DataFrame({"key": ["foo", "bar"], "rval": [4, 5]})
   left
   right
   pd.merge(left, right, on="key")

Grouping
--------

By "group by" we are referring to a process involving one or more of the
following steps:

* **Splitting** the data into groups based on some criteria
* **Applying** a function to each group independently
* **Combining** the results into a data structure

See the :ref:`Grouping section <groupby>`.

.. ipython:: python

   df = pd.DataFrame(
       {
           "A": ["foo", "bar", "foo", "bar", "foo", "bar", "foo", "foo"],
           "B": ["one", "one", "two", "three", "two", "two", "one", "three"],
           "C": np.random.randn(8),
           "D": np.random.randn(8),
       }
   )
   df

Grouping by a column label, selecting column labels, and then applying the
:meth:`.DataFrameGroupBy.sum` function to the resulting
groups:

.. ipython:: python

   df.groupby("A")[["C", "D"]].sum()

Grouping by multiple columns label forms :class:`MultiIndex`.

.. ipython:: python

   df.groupby(["A", "B"]).sum()

Reshaping
---------

See the sections on :ref:`Hierarchical Indexing <advanced.hierarchical>` and
:ref:`Reshaping <reshaping.stacking>`.

Stack
~~~~~

.. ipython:: python

   arrays = [
      ["bar", "bar", "baz", "baz", "foo", "foo", "qux", "qux"],
      ["one", "two", "one", "two", "one", "two", "one", "two"],
   ]
   index = pd.MultiIndex.from_arrays(arrays, names=["first", "second"])
   df = pd.DataFrame(np.random.randn(8, 2), index=index, columns=["A", "B"])
   df2 = df[:4]
   df2

The :meth:`~DataFrame.stack` method "compresses" a level in the DataFrame's
columns:

.. ipython:: python

   stacked = df2.stack()
   stacked

With a "stacked" DataFrame or Series (having a :class:`MultiIndex` as the
``index``), the inverse operation of :meth:`~DataFrame.stack` is
:meth:`~DataFrame.unstack`, which by default unstacks the **last level**:

.. ipython:: python

   stacked.unstack()
   stacked.unstack(1)
   stacked.unstack(0)

Pivot tables
~~~~~~~~~~~~
See the section on :ref:`Pivot Tables <reshaping.pivot>`.

.. ipython:: python

   df = pd.DataFrame(
       {
           "A": ["one", "one", "two", "three"] * 3,
           "B": ["A", "B", "C"] * 4,
           "C": ["foo", "foo", "foo", "bar", "bar", "bar"] * 2,
           "D": np.random.randn(12),
           "E": np.random.randn(12),
       }
   )
   df

:func:`pivot_table` pivots a :class:`DataFrame` specifying the ``values``, ``index`` and ``columns``

.. ipython:: python

   pd.pivot_table(df, values="D", index=["A", "B"], columns=["C"])


Time series
-----------

pandas has simple, powerful, and efficient functionality for performing
resampling operations during frequency conversion (e.g., converting secondly
data into 5-minutely data). This is extremely common in, but not limited to,
financial applications. See the :ref:`Time Series section <timeseries>`.

.. ipython:: python

   rng = pd.date_range("1/1/2012", periods=100, freq="s")
   ts = pd.Series(np.random.randint(0, 500, len(rng)), index=rng)
   ts.resample("5Min").sum()

:meth:`Series.tz_localize` localizes a time series to a time zone:

.. ipython:: python

   rng = pd.date_range("3/6/2012 00:00", periods=5, freq="D")
   ts = pd.Series(np.random.randn(len(rng)), rng)
   ts
   ts_utc = ts.tz_localize("UTC")
   ts_utc

:meth:`Series.tz_convert` converts a timezones aware time series to another time zone:

.. ipython:: python

   ts_utc.tz_convert("US/Eastern")

Adding a non-fixed duration (:class:`~pandas.tseries.offsets.BusinessDay`) to a time series:

.. ipython:: python

   rng
   rng + pd.offsets.BusinessDay(5)

Categoricals
------------

pandas can include categorical data in a :class:`DataFrame`. For full docs, see the
:ref:`categorical introduction <categorical>` and the :ref:`API documentation <api.arrays.categorical>`.

.. ipython:: python

    df = pd.DataFrame(
        {"id": [1, 2, 3, 4, 5, 6], "raw_grade": ["a", "b", "b", "a", "a", "e"]}
    )

Converting the raw grades to a categorical data type:

.. ipython:: python

    df["grade"] = df["raw_grade"].astype("category")
    df["grade"]

Rename the categories to more meaningful names:

.. ipython:: python

    new_categories = ["very good", "good", "very bad"]
    df["grade"] = df["grade"].cat.rename_categories(new_categories)

Reorder the categories and simultaneously add the missing categories (methods under :meth:`Series.cat` return a new :class:`Series` by default):

.. ipython:: python

    df["grade"] = df["grade"].cat.set_categories(
        ["very bad", "bad", "medium", "good", "very good"]
    )
    df["grade"]

Sorting is per order in the categories, not lexical order:

.. ipython:: python

    df.sort_values(by="grade")

Grouping by a categorical column with ``observed=False`` also shows empty categories:

.. ipython:: python

    df.groupby("grade", observed=False).size()


Plotting
--------

See the :ref:`Plotting <visualization>` docs.

We use the standard convention for referencing the matplotlib API:

.. ipython:: python

   import matplotlib.pyplot as plt

   plt.close("all")

The ``plt.close`` method is used to `close <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.close.html>`__ a figure window:

.. ipython:: python

   ts = pd.Series(np.random.randn(1000), index=pd.date_range("1/1/2000", periods=1000))
   ts = ts.cumsum()

   @savefig series_plot_basic.png
   ts.plot();

.. note::

   When using Jupyter, the plot will appear using :meth:`~Series.plot`.  Otherwise use
   `matplotlib.pyplot.show <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.show.html>`__ to show it or
   `matplotlib.pyplot.savefig <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html>`__ to write it to a file.

:meth:`~DataFrame.plot` plots all columns:

.. ipython:: python

   df = pd.DataFrame(
       np.random.randn(1000, 4), index=ts.index, columns=["A", "B", "C", "D"]
   )

   df = df.cumsum()

   plt.figure();
   df.plot();
   @savefig frame_plot_basic.png
   plt.legend(loc='best');

Importing and exporting data
----------------------------

See the :ref:`IO Tools <io>` section.

CSV
~~~

:ref:`Writing to a csv file: <io.store_in_csv>` using :meth:`DataFrame.to_csv`

.. ipython:: python

   df = pd.DataFrame(np.random.randint(0, 5, (10, 5)))
   df.to_csv("foo.csv")

:ref:`Reading from a csv file: <io.read_csv_table>` using :func:`read_csv`

.. ipython:: python

   pd.read_csv("foo.csv")

.. ipython:: python
   :suppress:

   import os

   os.remove("foo.csv")

Parquet
~~~~~~~

Writing to a Parquet file:

.. ipython:: python

   df.to_parquet("foo.parquet")

Reading from a Parquet file Store using :func:`read_parquet`:

.. ipython:: python

   pd.read_parquet("foo.parquet")

.. ipython:: python
   :suppress:

   os.remove("foo.parquet")

Excel
~~~~~

Reading and writing to :ref:`Excel <io.excel>`.

Writing to an excel file using :meth:`DataFrame.to_excel`:

.. ipython:: python

   df.to_excel("foo.xlsx", sheet_name="Sheet1")

Reading from an excel file using :func:`read_excel`:

.. ipython:: python

   pd.read_excel("foo.xlsx", "Sheet1", index_col=None, na_values=["NA"])

.. ipython:: python
   :suppress:

   os.remove("foo.xlsx")

Gotchas
-------

If you are attempting to perform a boolean operation on a :class:`Series` or :class:`DataFrame`
you might see an exception like:

.. ipython:: python
   :okexcept:

    if pd.Series([False, True, False]):
        print("I was true")

See :ref:`Comparisons<basics.compare>` and :ref:`Gotchas<gotchas>` for an explanation and what to do.
