.. _reshaping:

{{ header }}

**************************
Reshaping and pivot tables
**************************

.. _reshaping.reshaping:


pandas provides methods for manipulating a :class:`Series` and :class:`DataFrame` to alter the
representation of the data for further data processing or data summarization.

* :func:`~pandas.pivot` and :func:`~pandas.pivot_table`: Group unique values within one or more discrete categories.
* :meth:`~DataFrame.stack` and :meth:`~DataFrame.unstack`: Pivot a column or row level to the opposite axis respectively.
* :func:`~pandas.melt` and :func:`~pandas.wide_to_long`: Unpivot a wide :class:`DataFrame` to a long format.
* :func:`~pandas.get_dummies` and :func:`~pandas.from_dummies`: Conversions with indicator variables.
* :meth:`~Series.explode`: Convert a column of list-like values to individual rows.
* :func:`~pandas.crosstab`: Calculate a cross-tabulation of multiple 1 dimensional factor arrays.
* :func:`~pandas.cut`: Transform continuous variables to discrete, categorical values
* :func:`~pandas.factorize`: Encode 1 dimensional variables into integer labels.


:func:`~pandas.pivot` and :func:`~pandas.pivot_table`
-----------------------------------------------------

.. image:: ../_static/reshaping_pivot.png

:func:`~pandas.pivot`
~~~~~~~~~~~~~~~~~~~~~

Data is often stored in so-called "stacked" or "record" format. In a "record" or "wide" format,
typically there is one row for each subject. In the "stacked" or "long" format there are
multiple rows for each subject where applicable.

.. ipython:: python

   data = {
      "value": range(12),
      "variable": ["A"] * 3 + ["B"] * 3 + ["C"] * 3 + ["D"] * 3,
      "date": pd.to_datetime(["2020-01-03", "2020-01-04", "2020-01-05"] * 4)
   }
   df = pd.DataFrame(data)

To perform time series operations with each unique variable, a better
representation would be where the ``columns`` are the unique variables and an
``index`` of dates identifies individual observations. To reshape the data into
this form, we use the :meth:`DataFrame.pivot` method (also implemented as a
top level function :func:`~pandas.pivot`):

.. ipython:: python

   pivoted = df.pivot(index="date", columns="variable", values="value")
   pivoted

If the ``values`` argument is omitted, and the input :class:`DataFrame` has more than
one column of values which are not used as column or index inputs to :meth:`~DataFrame.pivot`,
then the resulting "pivoted" :class:`DataFrame` will have :ref:`hierarchical columns
<advanced.hierarchical>` whose topmost level indicates the respective value
column:

.. ipython:: python

   df["value2"] = df["value"] * 2
   pivoted = df.pivot(index="date", columns="variable")
   pivoted

You can then select subsets from the pivoted :class:`DataFrame`:

.. ipython:: python

   pivoted["value2"]

Note that this returns a view on the underlying data in the case where the data
are homogeneously-typed.

.. note::

   :func:`~pandas.pivot` can only handle unique rows specified by ``index`` and ``columns``.
   If you data contains duplicates, use :func:`~pandas.pivot_table`.


.. _reshaping.pivot:

:func:`~pandas.pivot_table`
~~~~~~~~~~~~~~~~~~~~~~~~~~~

While :meth:`~DataFrame.pivot` provides general purpose pivoting with various
data types, pandas also provides :func:`~pandas.pivot_table` or :meth:`~DataFrame.pivot_table`
for pivoting with aggregation of numeric data.

The function :func:`~pandas.pivot_table` can be used to create spreadsheet-style
pivot tables. See the :ref:`cookbook<cookbook.pivot>` for some advanced
strategies.

.. ipython:: python

   import datetime

   df = pd.DataFrame(
       {
           "A": ["one", "one", "two", "three"] * 6,
           "B": ["A", "B", "C"] * 8,
           "C": ["foo", "foo", "foo", "bar", "bar", "bar"] * 4,
           "D": np.random.randn(24),
           "E": np.random.randn(24),
           "F": [datetime.datetime(2013, i, 1) for i in range(1, 13)]
           + [datetime.datetime(2013, i, 15) for i in range(1, 13)],
       }
   )
   df
   pd.pivot_table(df, values="D", index=["A", "B"], columns=["C"])
   pd.pivot_table(
       df, values=["D", "E"],
       index=["B"],
       columns=["A", "C"],
       aggfunc="sum",
   )
   pd.pivot_table(
       df, values="E",
       index=["B", "C"],
       columns=["A"],
       aggfunc=["sum", "mean"],
   )

The result is a :class:`DataFrame` potentially having a :class:`MultiIndex` on the
index or column. If the ``values`` column name is not given, the pivot table
will include all of the data in an additional level of hierarchy in the columns:

.. ipython:: python

   pd.pivot_table(df[["A", "B", "C", "D", "E"]], index=["A", "B"], columns=["C"])

Also, you can use :class:`Grouper` for ``index`` and ``columns`` keywords. For detail of :class:`Grouper`, see :ref:`Grouping with a Grouper specification <groupby.specify>`.

.. ipython:: python

   pd.pivot_table(df, values="D", index=pd.Grouper(freq="ME", key="F"), columns="C")

.. _reshaping.pivot.margins:

Adding margins
^^^^^^^^^^^^^^

Passing ``margins=True`` to :meth:`~DataFrame.pivot_table` will add a row and column with an
``All`` label with partial group aggregates across the categories on the
rows and columns:

.. ipython:: python

   table = df.pivot_table(
       index=["A", "B"],
       columns="C",
       values=["D", "E"],
       margins=True,
       aggfunc="std"
   )
   table

Additionally, you can call :meth:`DataFrame.stack` to display a pivoted DataFrame
as having a multi-level index:

.. ipython:: python

    table.stack()

.. _reshaping.stacking:

:meth:`~DataFrame.stack` and :meth:`~DataFrame.unstack`
-------------------------------------------------------

.. image:: ../_static/reshaping_stack.png

Closely related to the :meth:`~DataFrame.pivot` method are the related
:meth:`~DataFrame.stack` and :meth:`~DataFrame.unstack` methods available on
:class:`Series` and :class:`DataFrame`. These methods are designed to work together with
:class:`MultiIndex` objects (see the section on :ref:`hierarchical indexing
<advanced.hierarchical>`).

* :meth:`~DataFrame.stack`: "pivot" a level of the (possibly hierarchical) column labels,
  returning a :class:`DataFrame` with an index with a new inner-most level of row
  labels.
* :meth:`~DataFrame.unstack`: (inverse operation of :meth:`~DataFrame.stack`) "pivot" a level of the
  (possibly hierarchical) row index to the column axis, producing a reshaped
  :class:`DataFrame` with a new inner-most level of column labels.

.. image:: ../_static/reshaping_unstack.png

.. ipython:: python

   tuples = [
      ["bar", "bar", "baz", "baz", "foo", "foo", "qux", "qux"],
      ["one", "two", "one", "two", "one", "two", "one", "two"],
   ]
   index = pd.MultiIndex.from_arrays(tuples, names=["first", "second"])
   df = pd.DataFrame(np.random.randn(8, 2), index=index, columns=["A", "B"])
   df2 = df[:4]
   df2

The :meth:`~DataFrame.stack` function "compresses" a level in the :class:`DataFrame` columns to
produce either:

* A :class:`Series`, in the case of a :class:`Index` in the columns.
* A :class:`DataFrame`, in the case of a :class:`MultiIndex` in the columns.

If the columns have a :class:`MultiIndex`, you can choose which level to stack. The
stacked level becomes the new lowest level in a :class:`MultiIndex` on the columns:

.. ipython:: python

   stacked = df2.stack()
   stacked

With a "stacked" :class:`DataFrame` or :class:`Series` (having a :class:`MultiIndex` as the
``index``), the inverse operation of :meth:`~DataFrame.stack` is :meth:`~DataFrame.unstack`, which by default
unstacks the **last level**:

.. ipython:: python

   stacked.unstack()
   stacked.unstack(1)
   stacked.unstack(0)

.. _reshaping.unstack_by_name:

.. image:: ../_static/reshaping_unstack_1.png

If the indexes have names, you can use the level names instead of specifying
the level numbers:

.. ipython:: python

   stacked.unstack("second")


.. image:: ../_static/reshaping_unstack_0.png

Notice that the :meth:`~DataFrame.stack` and :meth:`~DataFrame.unstack` methods implicitly sort the index
levels involved. Hence a call to :meth:`~DataFrame.stack` and then :meth:`~DataFrame.unstack`, or vice versa,
will result in a **sorted** copy of the original :class:`DataFrame` or :class:`Series`:

.. ipython:: python

   index = pd.MultiIndex.from_product([[2, 1], ["a", "b"]])
   df = pd.DataFrame(np.random.randn(4), index=index, columns=["A"])
   df
   all(df.unstack().stack() == df.sort_index())

.. _reshaping.stack_multiple:

Multiple levels
~~~~~~~~~~~~~~~

You may also stack or unstack more than one level at a time by passing a list
of levels, in which case the end result is as if each level in the list were
processed individually.

.. ipython:: python

    columns = pd.MultiIndex.from_tuples(
        [
            ("A", "cat", "long"),
            ("B", "cat", "long"),
            ("A", "dog", "short"),
            ("B", "dog", "short"),
        ],
        names=["exp", "animal", "hair_length"],
    )
    df = pd.DataFrame(np.random.randn(4, 4), columns=columns)
    df

    df.stack(level=["animal", "hair_length"])

The list of levels can contain either level names or level numbers but
not a mixture of the two.

.. ipython:: python

    # df.stack(level=['animal', 'hair_length'])
    # from above is equivalent to:
    df.stack(level=[1, 2])

Missing data
~~~~~~~~~~~~

Unstacking can result in missing values if subgroups do not have the same
set of labels. By default, missing values will be replaced with the default
fill value for that data type.

.. ipython:: python

   columns = pd.MultiIndex.from_tuples(
       [
           ("A", "cat"),
           ("B", "dog"),
           ("B", "cat"),
           ("A", "dog"),
       ],
       names=["exp", "animal"],
   )
   index = pd.MultiIndex.from_product(
       [("bar", "baz", "foo", "qux"), ("one", "two")], names=["first", "second"]
   )
   df = pd.DataFrame(np.random.randn(8, 4), index=index, columns=columns)
   df3 = df.iloc[[0, 1, 4, 7], [1, 2]]
   df3
   df3.unstack()

The missing value can be filled with a specific value with the ``fill_value`` argument.

.. ipython:: python

   df3.unstack(fill_value=-1e9)

.. _reshaping.melt:

:func:`~pandas.melt` and :func:`~pandas.wide_to_long`
-----------------------------------------------------

.. image:: ../_static/reshaping_melt.png

The top-level :func:`~pandas.melt` function and the corresponding :meth:`DataFrame.melt`
are useful to massage a :class:`DataFrame` into a format where one or more columns
are *identifier variables*, while all other columns, considered *measured
variables*, are "unpivoted" to the row axis, leaving just two non-identifier
columns, "variable" and "value". The names of those columns can be customized
by supplying the ``var_name`` and ``value_name`` parameters.

.. ipython:: python

   cheese = pd.DataFrame(
       {
           "first": ["John", "Mary"],
           "last": ["Doe", "Bo"],
           "height": [5.5, 6.0],
           "weight": [130, 150],
       }
   )
   cheese
   cheese.melt(id_vars=["first", "last"])
   cheese.melt(id_vars=["first", "last"], var_name="quantity")

When transforming a DataFrame using :func:`~pandas.melt`, the index will be ignored.
The original index values can be kept by setting the ``ignore_index=False`` parameter to ``False`` (default is ``True``).
``ignore_index=False`` will however duplicate index values.

.. ipython:: python

   index = pd.MultiIndex.from_tuples([("person", "A"), ("person", "B")])
   cheese = pd.DataFrame(
       {
           "first": ["John", "Mary"],
           "last": ["Doe", "Bo"],
           "height": [5.5, 6.0],
           "weight": [130, 150],
       },
       index=index,
   )
   cheese
   cheese.melt(id_vars=["first", "last"])
   cheese.melt(id_vars=["first", "last"], ignore_index=False)

:func:`~pandas.wide_to_long` is similar to :func:`~pandas.melt` with more customization for
column matching.

.. ipython:: python

  dft = pd.DataFrame(
      {
          "A1970": {0: "a", 1: "b", 2: "c"},
          "A1980": {0: "d", 1: "e", 2: "f"},
          "B1970": {0: 2.5, 1: 1.2, 2: 0.7},
          "B1980": {0: 3.2, 1: 1.3, 2: 0.1},
          "X": dict(zip(range(3), np.random.randn(3))),
      }
  )
  dft["id"] = dft.index
  dft
  pd.wide_to_long(dft, ["A", "B"], i="id", j="year")

.. _reshaping.dummies:

:func:`~pandas.get_dummies` and :func:`~pandas.from_dummies`
------------------------------------------------------------

To convert categorical variables of a :class:`Series` into a "dummy" or "indicator",
:func:`~pandas.get_dummies` creates a new :class:`DataFrame` with columns of the unique
variables and the values representing the presence of those variables per row.

.. ipython:: python

   df = pd.DataFrame({"key": list("bbacab"), "data1": range(6)})

   pd.get_dummies(df["key"])
   df["key"].str.get_dummies()

``prefix`` adds a prefix to the the column names which is useful for merging the result
with the original :class:`DataFrame`:

.. ipython:: python

   dummies = pd.get_dummies(df["key"], prefix="key")
   dummies

   df[["data1"]].join(dummies)

This function is often used along with discretization functions like :func:`~pandas.cut`:

.. ipython:: python

   values = np.random.randn(10)
   values

   bins = [0, 0.2, 0.4, 0.6, 0.8, 1]

   pd.get_dummies(pd.cut(values, bins))


:func:`get_dummies` also accepts a :class:`DataFrame`. By default, ``object``, ``string``,
or ``categorical`` type columns are encoded as dummy variables with other columns unaltered.

.. ipython:: python

    df = pd.DataFrame({"A": ["a", "b", "a"], "B": ["c", "c", "b"], "C": [1, 2, 3]})
    pd.get_dummies(df)

Specifying the ``columns`` keyword will encode a column of any type.

.. ipython:: python

    pd.get_dummies(df, columns=["A"])

As with the :class:`Series` version, you can pass values for the ``prefix`` and
``prefix_sep``. By default the column name is used as the prefix and ``_`` as
the prefix separator. You can specify ``prefix`` and ``prefix_sep`` in 3 ways:

* string: Use the same value for ``prefix`` or ``prefix_sep`` for each column
  to be encoded.
* list: Must be the same length as the number of columns being encoded.
* dict: Mapping column name to prefix.

.. ipython:: python

    simple = pd.get_dummies(df, prefix="new_prefix")
    simple
    from_list = pd.get_dummies(df, prefix=["from_A", "from_B"])
    from_list
    from_dict = pd.get_dummies(df, prefix={"B": "from_B", "A": "from_A"})
    from_dict

To avoid collinearity when feeding the result to statistical models,
specify ``drop_first=True``.

.. ipython:: python

    s = pd.Series(list("abcaa"))

    pd.get_dummies(s)

    pd.get_dummies(s, drop_first=True)

When a column contains only one level, it will be omitted in the result.

.. ipython:: python

    df = pd.DataFrame({"A": list("aaaaa"), "B": list("ababc")})

    pd.get_dummies(df)

    pd.get_dummies(df, drop_first=True)

The values can be cast to a different type using the ``dtype`` argument.

.. ipython:: python

    df = pd.DataFrame({"A": list("abc"), "B": [1.1, 2.2, 3.3]})

    pd.get_dummies(df, dtype=np.float32).dtypes

.. versionadded:: 1.5.0

:func:`~pandas.from_dummies` converts the output of :func:`~pandas.get_dummies` back into
a :class:`Series` of categorical values from indicator values.

.. ipython:: python

   df = pd.DataFrame({"prefix_a": [0, 1, 0], "prefix_b": [1, 0, 1]})
   df

   pd.from_dummies(df, sep="_")

Dummy coded data only requires ``k - 1`` categories to be included, in this case
the last category is the default category. The default category can be modified with
``default_category``.

.. ipython:: python

   df = pd.DataFrame({"prefix_a": [0, 1, 0]})
   df

   pd.from_dummies(df, sep="_", default_category="b")

.. _reshaping.explode:

:meth:`~Series.explode`
-----------------------

For a :class:`DataFrame` column with nested, list-like values, :meth:`~Series.explode` will transform
each list-like value to a separate row. The resulting :class:`Index` will be duplicated corresponding
to the index label from the original row:

.. ipython:: python

   keys = ["panda1", "panda2", "panda3"]
   values = [["eats", "shoots"], ["shoots", "leaves"], ["eats", "leaves"]]
   df = pd.DataFrame({"keys": keys, "values": values})
   df
   df["values"].explode()

:class:`DataFrame.explode` can also explode the column in the :class:`DataFrame`.

.. ipython:: python

   df.explode("values")

:meth:`Series.explode` will replace empty lists with a missing value indicator and preserve scalar entries.

.. ipython:: python

   s = pd.Series([[1, 2, 3], "foo", [], ["a", "b"]])
   s
   s.explode()

A comma-separated string value can be split into individual values in a list and then exploded to a new row.

.. ipython:: python

    df = pd.DataFrame([{"var1": "a,b,c", "var2": 1}, {"var1": "d,e,f", "var2": 2}])
    df.assign(var1=df.var1.str.split(",")).explode("var1")

.. _reshaping.crosstabulations:

:func:`~pandas.crosstab`
------------------------

Use :func:`~pandas.crosstab` to compute a cross-tabulation of two (or more)
factors. By default :func:`~pandas.crosstab` computes a frequency table of the factors
unless an array of values and an aggregation function are passed.

Any :class:`Series` passed will have their name attributes used unless row or column
names for the cross-tabulation are specified

.. ipython:: python

    a = np.array(["foo", "foo", "bar", "bar", "foo", "foo"], dtype=object)
    b = np.array(["one", "one", "two", "one", "two", "one"], dtype=object)
    c = np.array(["dull", "dull", "shiny", "dull", "dull", "shiny"], dtype=object)
    pd.crosstab(a, [b, c], rownames=["a"], colnames=["b", "c"])


If :func:`~pandas.crosstab` receives only two :class:`Series`, it will provide a frequency table.

.. ipython:: python

    df = pd.DataFrame(
        {"A": [1, 2, 2, 2, 2], "B": [3, 3, 4, 4, 4], "C": [1, 1, np.nan, 1, 1]}
    )
    df

    pd.crosstab(df["A"], df["B"])

:func:`~pandas.crosstab` can also summarize to :class:`Categorical` data.

.. ipython:: python

    foo = pd.Categorical(["a", "b"], categories=["a", "b", "c"])
    bar = pd.Categorical(["d", "e"], categories=["d", "e", "f"])
    pd.crosstab(foo, bar)

For :class:`Categorical` data, to include **all** of data categories even if the actual data does
not contain any instances of a particular category, use ``dropna=False``.

.. ipython:: python

    pd.crosstab(foo, bar, dropna=False)

Normalization
~~~~~~~~~~~~~

Frequency tables can also be normalized to show percentages rather than counts
using the ``normalize`` argument:

.. ipython:: python

   pd.crosstab(df["A"], df["B"], normalize=True)

``normalize`` can also normalize values within each row or within each column:

.. ipython:: python

   pd.crosstab(df["A"], df["B"], normalize="columns")

:func:`~pandas.crosstab` can also accept a third :class:`Series` and an aggregation function
(``aggfunc``) that will be applied to the values of the third :class:`Series` within
each group defined by the first two :class:`Series`:

.. ipython:: python

   pd.crosstab(df["A"], df["B"], values=df["C"], aggfunc="sum")

Adding margins
~~~~~~~~~~~~~~

``margins=True`` will add a row and column with an ``All`` label with partial group aggregates
across the categories on the rows and columns:

.. ipython:: python

   pd.crosstab(
       df["A"], df["B"], values=df["C"], aggfunc="sum", normalize=True, margins=True
   )

.. _reshaping.tile:
.. _reshaping.tile.cut:

:func:`~pandas.cut`
-------------------

The :func:`~pandas.cut` function computes groupings for the values of the input
array and is often used to transform continuous variables to discrete or
categorical variables:


An integer ``bins`` will form equal-width bins.

.. ipython:: python

   ages = np.array([10, 15, 13, 12, 23, 25, 28, 59, 60])

   pd.cut(ages, bins=3)

A list of ordered bin edges will assign an interval for each variable.

.. ipython:: python

   pd.cut(ages, bins=[0, 18, 35, 70])

If the ``bins`` keyword is an :class:`IntervalIndex`, then these will be
used to bin the passed data.

.. ipython:: python

   pd.cut(ages, bins=pd.IntervalIndex.from_breaks([0, 40, 70]))

.. _reshaping.factorize:

:func:`~pandas.factorize`
-------------------------

:func:`~pandas.factorize` encodes 1 dimensional values into integer labels. Missing values
are encoded as ``-1``.

.. ipython:: python

   x = pd.Series(["A", "A", np.nan, "B", 3.14, np.inf])
   x
   labels, uniques = pd.factorize(x)
   labels
   uniques

:class:`Categorical` will similarly encode 1 dimensional values for further
categorical operations

.. ipython:: python

   pd.Categorical(x)
