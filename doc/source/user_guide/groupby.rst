.. _groupby:

{{ header }}

*****************************
Group by: split-apply-combine
*****************************

By "group by" we are referring to a process involving one or more of the following
steps:

* **Splitting** the data into groups based on some criteria.
* **Applying** a function to each group independently.
* **Combining** the results into a data structure.

Out of these, the split step is the most straightforward. In the apply step, we
might wish to do one of the following:

* **Aggregation**: compute a summary statistic (or statistics) for each
  group. Some examples:

    * Compute group sums or means.
    * Compute group sizes / counts.

* **Transformation**: perform some group-specific computations and return a
  like-indexed object. Some examples:

    * Standardize data (zscore) within a group.
    * Filling NAs within groups with a value derived from each group.

* **Filtration**: discard some groups, according to a group-wise computation
  that evaluates to True or False. Some examples:

    * Discard data that belong to groups with only a few members.
    * Filter out data based on the group sum or mean.

Many of these operations are defined on GroupBy objects. These operations are similar
to those of the :ref:`aggregating API <basics.aggregate>`,
:ref:`window API <window.overview>`, and :ref:`resample API <timeseries.aggregate>`.

It is possible that a given operation does not fall into one of these categories or
is some combination of them. In such a case, it may be possible to compute the
operation using GroupBy's ``apply`` method. This method will examine the results of the
apply step and try to sensibly combine them into a single result if it doesn't fit into either
of the above three categories.

.. note::

   An operation that is split into multiple steps using built-in GroupBy operations
   will be more efficient than using the ``apply`` method with a user-defined Python
   function.


The name GroupBy should be quite familiar to those who have used
a SQL-based tool (or ``itertools``), in which you can write code like:

.. code-block:: sql

   SELECT Column1, Column2, mean(Column3), sum(Column4)
   FROM SomeTable
   GROUP BY Column1, Column2

We aim to make operations like this natural and easy to express using
pandas. We'll address each area of GroupBy functionality, then provide some
non-trivial examples / use cases.

See the :ref:`cookbook<cookbook.grouping>` for some advanced strategies.

.. _groupby.split:

Splitting an object into groups
-------------------------------

The abstract definition of grouping is to provide a mapping of labels to
group names. To create a GroupBy object (more on what the GroupBy object is
later), you may do the following:

.. ipython:: python

    speeds = pd.DataFrame(
        [
            ("bird", "Falconiformes", 389.0),
            ("bird", "Psittaciformes", 24.0),
            ("mammal", "Carnivora", 80.2),
            ("mammal", "Primates", np.nan),
            ("mammal", "Carnivora", 58),
        ],
        index=["falcon", "parrot", "lion", "monkey", "leopard"],
        columns=("class", "order", "max_speed"),
    )
    speeds

    grouped = speeds.groupby("class")
    grouped = speeds.groupby(["class", "order"])

The mapping can be specified many different ways:

* A Python function, to be called on each of the index labels.
* A list or NumPy array of the same length as the index.
* A dict or ``Series``, providing a ``label -> group name`` mapping.
* For ``DataFrame`` objects, a string indicating either a column name or
  an index level name to be used to group.
* A list of any of the above things.

Collectively we refer to the grouping objects as the **keys**. For example,
consider the following ``DataFrame``:

.. note::

   A string passed to ``groupby`` may refer to either a column or an index level.
   If a string matches both a column name and an index level name, a
   ``ValueError`` will be raised.

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

On a DataFrame, we obtain a GroupBy object by calling :meth:`~DataFrame.groupby`.
This method returns a ``pandas.api.typing.DataFrameGroupBy`` instance.
We could naturally group by either the ``A`` or ``B`` columns, or both:

.. ipython:: python

   grouped = df.groupby("A")
   grouped = df.groupby("B")
   grouped = df.groupby(["A", "B"])

.. note::

   ``df.groupby('A')`` is just syntactic sugar for ``df.groupby(df['A'])``.

If we also have a MultiIndex on columns ``A`` and ``B``, we can group by all
the columns except the one we specify:

.. ipython:: python

   df2 = df.set_index(["A", "B"])
   grouped = df2.groupby(level=df2.index.names.difference(["B"]))
   grouped.sum()

The above GroupBy will split the DataFrame on its index (rows). To split by columns, first do
a transpose:

.. ipython::

    In [4]: def get_letter_type(letter):
       ...:     if letter.lower() in 'aeiou':
       ...:         return 'vowel'
       ...:     else:
       ...:         return 'consonant'
       ...:

    In [5]: grouped = df.T.groupby(get_letter_type)

pandas :class:`~pandas.Index` objects support duplicate values. If a
non-unique index is used as the group key in a groupby operation, all values
for the same index value will be considered to be in one group and thus the
output of aggregation functions will only contain unique index values:

.. ipython:: python

   index = [1, 2, 3, 1, 2, 3]

   s = pd.Series([1, 2, 3, 10, 20, 30], index=index)

   s

   grouped = s.groupby(level=0)

   grouped.first()

   grouped.last()

   grouped.sum()

Note that **no splitting occurs** until it's needed. Creating the GroupBy object
only verifies that you've passed a valid mapping.

.. note::

   Many kinds of complicated data manipulations can be expressed in terms of
   GroupBy operations (though it can't be guaranteed to be the most efficient implementation).
   You can get quite creative with the label mapping functions.

.. _groupby.sorting:

GroupBy sorting
~~~~~~~~~~~~~~~~~~~~~~~~~

By default the group keys are sorted during the ``groupby`` operation. You may however pass ``sort=False`` for potential speedups. With ``sort=False`` the order among group-keys follows the order of appearance of the keys in the original dataframe:

.. ipython:: python

   df2 = pd.DataFrame({"X": ["B", "B", "A", "A"], "Y": [1, 2, 3, 4]})
   df2.groupby(["X"]).sum()
   df2.groupby(["X"], sort=False).sum()


Note that ``groupby`` will preserve the order in which *observations* are sorted *within* each group.
For example, the groups created by ``groupby()`` below are in the order they appeared in the original ``DataFrame``:

.. ipython:: python

   df3 = pd.DataFrame({"X": ["A", "B", "A", "B"], "Y": [1, 4, 3, 2]})
   df3.groupby("X").get_group("A")

   df3.groupby(["X"]).get_group(("B",))


.. _groupby.dropna:

GroupBy dropna
^^^^^^^^^^^^^^

By default ``NA`` values are excluded from group keys during the ``groupby`` operation. However,
in case you want to include ``NA`` values in group keys, you could pass ``dropna=False`` to achieve it.

.. ipython:: python

    df_list = [[1, 2, 3], [1, None, 4], [2, 1, 3], [1, 2, 2]]
    df_dropna = pd.DataFrame(df_list, columns=["a", "b", "c"])

    df_dropna

.. ipython:: python

    # Default ``dropna`` is set to True, which will exclude NaNs in keys
    df_dropna.groupby(by=["b"], dropna=True).sum()

    # In order to allow NaN in keys, set ``dropna`` to False
    df_dropna.groupby(by=["b"], dropna=False).sum()

The default setting of ``dropna`` argument is ``True`` which means ``NA`` are not included in group keys.


.. _groupby.attributes:

GroupBy object attributes
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``groups`` attribute is a dictionary whose keys are the computed unique groups
and corresponding values are the axis labels belonging to each group. In the
above example we have:

.. ipython:: python

   df.groupby("A").groups
   df.T.groupby(get_letter_type).groups

Calling the standard Python ``len`` function on the GroupBy object returns
the number of groups, which is the same as the length of the ``groups`` dictionary:

.. ipython:: python

   grouped = df.groupby(["A", "B"])
   grouped.groups
   len(grouped)


.. _groupby.tabcompletion:

``GroupBy`` will tab complete column names, GroupBy operations, and other attributes:

.. ipython:: python

   n = 10
   weight = np.random.normal(166, 20, size=n)
   height = np.random.normal(60, 10, size=n)
   time = pd.date_range("1/1/2000", periods=n)
   gender = np.random.choice(["male", "female"], size=n)
   df = pd.DataFrame(
       {"height": height, "weight": weight, "gender": gender}, index=time
   )
   df
   gb = df.groupby("gender")


.. ipython::

   @verbatim
   In [1]: gb.<TAB>  # noqa: E225, E999
   gb.agg        gb.boxplot    gb.cummin     gb.describe   gb.filter     gb.get_group  gb.height     gb.last       gb.median     gb.ngroups    gb.plot       gb.rank       gb.std        gb.transform
   gb.aggregate  gb.count      gb.cumprod    gb.dtype      gb.first      gb.groups     gb.hist       gb.max        gb.min        gb.nth        gb.prod       gb.resample   gb.sum        gb.var
   gb.apply      gb.cummax     gb.cumsum     gb.fillna     gb.gender     gb.head       gb.indices    gb.mean       gb.name       gb.ohlc       gb.quantile   gb.size       gb.tail       gb.weight

.. _groupby.multiindex:

GroupBy with MultiIndex
~~~~~~~~~~~~~~~~~~~~~~~

With :ref:`hierarchically-indexed data <advanced.hierarchical>`, it's quite
natural to group by one of the levels of the hierarchy.

Let's create a Series with a two-level ``MultiIndex``.

.. ipython:: python


   arrays = [
       ["bar", "bar", "baz", "baz", "foo", "foo", "qux", "qux"],
       ["one", "two", "one", "two", "one", "two", "one", "two"],
   ]
   index = pd.MultiIndex.from_arrays(arrays, names=["first", "second"])
   s = pd.Series(np.random.randn(8), index=index)
   s

We can then group by one of the levels in ``s``.

.. ipython:: python

   grouped = s.groupby(level=0)
   grouped.sum()

If the MultiIndex has names specified, these can be passed instead of the level
number:

.. ipython:: python

   s.groupby(level="second").sum()

Grouping with multiple levels is supported.

.. ipython:: python

   arrays = [
       ["bar", "bar", "baz", "baz", "foo", "foo", "qux", "qux"],
       ["doo", "doo", "bee", "bee", "bop", "bop", "bop", "bop"],
       ["one", "two", "one", "two", "one", "two", "one", "two"],
   ]
   index = pd.MultiIndex.from_arrays(arrays, names=["first", "second", "third"])
   s = pd.Series(np.random.randn(8), index=index)
   s
   s.groupby(level=["first", "second"]).sum()

Index level names may be supplied as keys.

.. ipython:: python

   s.groupby(["first", "second"]).sum()

More on the ``sum`` function and aggregation later.

Grouping DataFrame with Index levels and columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A DataFrame may be grouped by a combination of columns and index levels. You
can specify both column and index names, or use a :class:`Grouper`.

Let's first create a DataFrame with a MultiIndex:

.. ipython:: python

   arrays = [
       ["bar", "bar", "baz", "baz", "foo", "foo", "qux", "qux"],
       ["one", "two", "one", "two", "one", "two", "one", "two"],
   ]

   index = pd.MultiIndex.from_arrays(arrays, names=["first", "second"])

   df = pd.DataFrame({"A": [1, 1, 1, 1, 2, 2, 3, 3], "B": np.arange(8)}, index=index)

   df

Then we group ``df`` by the ``second`` index level and the ``A`` column.

.. ipython:: python

   df.groupby([pd.Grouper(level=1), "A"]).sum()

Index levels may also be specified by name.

.. ipython:: python

   df.groupby([pd.Grouper(level="second"), "A"]).sum()

Index level names may be specified as keys directly to ``groupby``.

.. ipython:: python

   df.groupby(["second", "A"]).sum()

DataFrame column selection in GroupBy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once you have created the GroupBy object from a DataFrame, you might want to do
something different for each of the columns. Thus, by using ``[]`` on the GroupBy
object in a similar way as the one used to get a column from a DataFrame, you can do:

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

   grouped = df.groupby(["A"])
   grouped_C = grouped["C"]
   grouped_D = grouped["D"]

This is mainly syntactic sugar for the alternative, which is much more verbose:

.. ipython:: python

   df["C"].groupby(df["A"])

Additionally, this method avoids recomputing the internal grouping information
derived from the passed key.

You can also include the grouping columns if you want to operate on them.

.. ipython:: python

   grouped[["A", "B"]].sum()

.. _groupby.iterating-label:

Iterating through groups
------------------------

With the GroupBy object in hand, iterating through the grouped data is very
natural and functions similarly to :py:func:`itertools.groupby`:

.. ipython::

   In [4]: grouped = df.groupby('A')

   In [5]: for name, group in grouped:
      ...:     print(name)
      ...:     print(group)
      ...:

In the case of grouping by multiple keys, the group name will be a tuple:

.. ipython::

   In [5]: for name, group in df.groupby(['A', 'B']):
      ...:     print(name)
      ...:     print(group)
      ...:

See :ref:`timeseries.iterating-label`.

Selecting a group
-----------------

A single group can be selected using
:meth:`.DataFrameGroupBy.get_group`:

.. ipython:: python

   grouped.get_group("bar")

Or for an object grouped on multiple columns:

.. ipython:: python

   df.groupby(["A", "B"]).get_group(("bar", "one"))

.. _groupby.aggregate:

Aggregation
-----------

An aggregation is a GroupBy operation that reduces the dimension of the grouping
object. The result of an aggregation is, or at least is treated as,
a scalar value for each column in a group. For example, producing the sum of each
column in a group of values.

.. ipython:: python

   animals = pd.DataFrame(
       {
           "kind": ["cat", "dog", "cat", "dog"],
           "height": [9.1, 6.0, 9.5, 34.0],
           "weight": [7.9, 7.5, 9.9, 198.0],
       }
   )
   animals
   animals.groupby("kind").sum()

In the result, the keys of the groups appear in the index by default. They can be
instead included in the columns by passing ``as_index=False``.

.. ipython:: python

   animals.groupby("kind", as_index=False).sum()

.. _groupby.aggregate.builtin:

Built-in aggregation methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Many common aggregations are built-in to GroupBy objects as methods. Of the methods
listed below, those with a ``*`` do *not* have an efficient, GroupBy-specific, implementation.

.. csv-table::
    :header: "Method", "Description"
    :widths: 20, 80
    :delim: ;

        :meth:`~.DataFrameGroupBy.any`;Compute whether any of the values in the groups are truthy
        :meth:`~.DataFrameGroupBy.all`;Compute whether all of the values in the groups are truthy
        :meth:`~.DataFrameGroupBy.count`;Compute the number of non-NA values in the groups
        :meth:`~.DataFrameGroupBy.cov` * ;Compute the covariance of the groups
        :meth:`~.DataFrameGroupBy.first`;Compute the first occurring value in each group
        :meth:`~.DataFrameGroupBy.idxmax`;Compute the index of the maximum value in each group
        :meth:`~.DataFrameGroupBy.idxmin`;Compute the index of the minimum value in each group
        :meth:`~.DataFrameGroupBy.last`;Compute the last occurring value in each group
        :meth:`~.DataFrameGroupBy.max`;Compute the maximum value in each group
        :meth:`~.DataFrameGroupBy.mean`;Compute the mean of each group
        :meth:`~.DataFrameGroupBy.median`;Compute the median of each group
        :meth:`~.DataFrameGroupBy.min`;Compute the minimum value in each group
        :meth:`~.DataFrameGroupBy.nunique`;Compute the number of unique values in each group
        :meth:`~.DataFrameGroupBy.prod`;Compute the product of the values in each group
        :meth:`~.DataFrameGroupBy.quantile`;Compute a given quantile of the values in each group
        :meth:`~.DataFrameGroupBy.sem`;Compute the standard error of the mean of the values in each group
        :meth:`~.DataFrameGroupBy.size`;Compute the number of values in each group
        :meth:`~.DataFrameGroupBy.skew` *;Compute the skew of the values in each group
        :meth:`~.DataFrameGroupBy.std`;Compute the standard deviation of the values in each group
        :meth:`~.DataFrameGroupBy.sum`;Compute the sum of the values in each group
        :meth:`~.DataFrameGroupBy.var`;Compute the variance of the values in each group

Some examples:

.. ipython:: python

   df.groupby("A")[["C", "D"]].max()
   df.groupby(["A", "B"]).mean()

Another aggregation example is to compute the size of each group.
This is included in GroupBy as the ``size`` method. It returns a Series whose
index consists of the group names and the values are the sizes of each group.

.. ipython:: python

   grouped = df.groupby(["A", "B"])
   grouped.size()

While the :meth:`.DataFrameGroupBy.describe` method is not itself a reducer, it
can be used to conveniently produce a collection of summary statistics about each of
the groups.

.. ipython:: python

   grouped.describe()

Another aggregation example is to compute the number of unique values of each group.
This is similar to the :meth:`.DataFrameGroupBy.value_counts` function, except that it only counts the
number of unique values.

.. ipython:: python

   ll = [['foo', 1], ['foo', 2], ['foo', 2], ['bar', 1], ['bar', 1]]
   df4 = pd.DataFrame(ll, columns=["A", "B"])
   df4
   df4.groupby("A")["B"].nunique()

.. note::

   Aggregation functions **will not** return the groups that you are aggregating over
   as named *columns* when ``as_index=True``, the default. The grouped columns will
   be the **indices** of the returned object.

   Passing ``as_index=False`` **will** return the groups that you are aggregating over as
   named columns, regardless if they are named **indices** or *columns* in the inputs.


.. _groupby.aggregate.agg:

The :meth:`~.DataFrameGroupBy.aggregate` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::
    The :meth:`~.DataFrameGroupBy.aggregate` method can accept many different types of
    inputs. This section details using string aliases for various GroupBy methods; other
    inputs are detailed in the sections below.

Any reduction method that pandas implements can be passed as a string to
:meth:`~.DataFrameGroupBy.aggregate`. Users are encouraged to use the shorthand,
``agg``. It will operate as if the corresponding method was called.

.. ipython:: python

   grouped = df.groupby("A")
   grouped[["C", "D"]].aggregate("sum")

   grouped = df.groupby(["A", "B"])
   grouped.agg("sum")

The result of the aggregation will have the group names as the
new index. In the case of multiple keys, the result is a
:ref:`MultiIndex <advanced.hierarchical>` by default. As mentioned above, this can be
changed by using the ``as_index`` option:

.. ipython:: python

   grouped = df.groupby(["A", "B"], as_index=False)
   grouped.agg("sum")

   df.groupby("A", as_index=False)[["C", "D"]].agg("sum")

Note that you could use the :meth:`DataFrame.reset_index` DataFrame function to achieve
the same result as the column names are stored in the resulting ``MultiIndex``, although
this will make an extra copy.

.. ipython:: python

   df.groupby(["A", "B"]).agg("sum").reset_index()

.. _groupby.aggregate.udf:

Aggregation with User-Defined Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users can also provide their own User-Defined Functions (UDFs) for custom aggregations.

.. warning::

    When aggregating with a UDF, the UDF should not mutate the
    provided ``Series``. See :ref:`gotchas.udf-mutation` for more information.

.. note::

    Aggregating with a UDF is often less performant than using
    the pandas built-in methods on GroupBy. Consider breaking up a complex operation
    into a chain of operations that utilize the built-in methods.

.. ipython:: python

   animals
   animals.groupby("kind")[["height"]].agg(lambda x: set(x))

The resulting dtype will reflect that of the aggregating function. If the results from different groups have
different dtypes, then a common dtype will be determined in the same way as ``DataFrame`` construction.

.. ipython:: python

   animals.groupby("kind")[["height"]].agg(lambda x: x.astype(int).sum())

.. _groupby.aggregate.multifunc:

Applying multiple functions at once
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On a grouped ``Series``, you can pass a list or dict of functions to
:meth:`SeriesGroupBy.agg`, outputting a DataFrame:

.. ipython:: python

   grouped = df.groupby("A")
   grouped["C"].agg(["sum", "mean", "std"])

On a grouped ``DataFrame``, you can pass a list of functions to
:meth:`DataFrameGroupBy.agg` to aggregate each
column, which produces an aggregated result with a hierarchical column index:

.. ipython:: python

   grouped[["C", "D"]].agg(["sum", "mean", "std"])


The resulting aggregations are named after the functions themselves. If you
need to rename, then you can add in a chained operation for a ``Series`` like this:

.. ipython:: python

   (
       grouped["C"]
       .agg(["sum", "mean", "std"])
       .rename(columns={"sum": "foo", "mean": "bar", "std": "baz"})
   )

For a grouped ``DataFrame``, you can rename in a similar manner:

.. ipython:: python

   (
       grouped[["C", "D"]].agg(["sum", "mean", "std"]).rename(
           columns={"sum": "foo", "mean": "bar", "std": "baz"}
       )
   )

.. note::

   In general, the output column names should be unique, but pandas will allow
   you apply to the same function (or two functions with the same name) to the same
   column.

   .. ipython:: python

      grouped["C"].agg(["sum", "sum"])


   pandas also allows you to provide multiple lambdas. In this case, pandas
   will mangle the name of the (nameless) lambda functions, appending ``_<i>``
   to each subsequent lambda.

   .. ipython:: python

      grouped["C"].agg([lambda x: x.max() - x.min(), lambda x: x.median() - x.mean()])


.. _groupby.aggregate.named:

Named aggregation
~~~~~~~~~~~~~~~~~

To support column-specific aggregation *with control over the output column names*, pandas
accepts the special syntax in :meth:`.DataFrameGroupBy.agg` and :meth:`.SeriesGroupBy.agg`, known as "named aggregation", where

- The keywords are the *output* column names
- The values are tuples whose first element is the column to select
  and the second element is the aggregation to apply to that column. pandas
  provides the :class:`NamedAgg` namedtuple with the fields ``['column', 'aggfunc']``
  to make it clearer what the arguments are. As usual, the aggregation can
  be a callable or a string alias.

.. ipython:: python

   animals

   animals.groupby("kind").agg(
       min_height=pd.NamedAgg(column="height", aggfunc="min"),
       max_height=pd.NamedAgg(column="height", aggfunc="max"),
       average_weight=pd.NamedAgg(column="weight", aggfunc="mean"),
   )


:class:`NamedAgg` is just a ``namedtuple``. Plain tuples are allowed as well.

.. ipython:: python

   animals.groupby("kind").agg(
       min_height=("height", "min"),
       max_height=("height", "max"),
       average_weight=("weight", "mean"),
   )


If the column names you want are not valid Python keywords, construct a dictionary
and unpack the keyword arguments

.. ipython:: python

   animals.groupby("kind").agg(
       **{
           "total weight": pd.NamedAgg(column="weight", aggfunc="sum")
       }
   )

When using named aggregation, additional keyword arguments are not passed through
to the aggregation functions; only pairs
of ``(column, aggfunc)`` should be passed as ``**kwargs``. If your aggregation functions
require additional arguments, apply them partially with :meth:`functools.partial`.

Named aggregation is also valid for Series groupby aggregations. In this case there's
no column selection, so the values are just the functions.

.. ipython:: python

   animals.groupby("kind").height.agg(
       min_height="min",
       max_height="max",
   )

Applying different functions to DataFrame columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By passing a dict to ``aggregate`` you can apply a different aggregation to the
columns of a DataFrame:

.. ipython:: python

   grouped.agg({"C": "sum", "D": lambda x: np.std(x, ddof=1)})

The function names can also be strings. In order for a string to be valid it
must be implemented on GroupBy:

.. ipython:: python

   grouped.agg({"C": "sum", "D": "std"})

.. _groupby.transform:

Transformation
--------------

A transformation is a GroupBy operation whose result is indexed the same
as the one being grouped. Common examples include :meth:`~.DataFrameGroupBy.cumsum` and
:meth:`~.DataFrameGroupBy.diff`.

.. ipython:: python

    speeds
    grouped = speeds.groupby("class")["max_speed"]
    grouped.cumsum()
    grouped.diff()

Unlike aggregations, the groupings that are used to split
the original object are not included in the result.

.. note::

    Since transformations do not include the groupings that are used to split the result,
    the arguments ``as_index`` and ``sort`` in :meth:`DataFrame.groupby` and
    :meth:`Series.groupby` have no effect.

A common use of a transformation is to add the result back into the original DataFrame.

.. ipython:: python

    result = speeds.copy()
    result["cumsum"] = grouped.cumsum()
    result["diff"] = grouped.diff()
    result

Built-in transformation methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following methods on GroupBy act as transformations.

.. csv-table::
    :header: "Method", "Description"
    :widths: 20, 80
    :delim: ;

        :meth:`~.DataFrameGroupBy.bfill`;Back fill NA values within each group
        :meth:`~.DataFrameGroupBy.cumcount`;Compute the cumulative count within each group
        :meth:`~.DataFrameGroupBy.cummax`;Compute the cumulative max within each group
        :meth:`~.DataFrameGroupBy.cummin`;Compute the cumulative min within each group
        :meth:`~.DataFrameGroupBy.cumprod`;Compute the cumulative product within each group
        :meth:`~.DataFrameGroupBy.cumsum`;Compute the cumulative sum within each group
        :meth:`~.DataFrameGroupBy.diff`;Compute the difference between adjacent values within each group
        :meth:`~.DataFrameGroupBy.ffill`;Forward fill NA values within each group
        :meth:`~.DataFrameGroupBy.pct_change`;Compute the percent change between adjacent values within each group
        :meth:`~.DataFrameGroupBy.rank`;Compute the rank of each value within each group
        :meth:`~.DataFrameGroupBy.shift`;Shift values up or down within each group

In addition, passing any built-in aggregation method as a string to
:meth:`~.DataFrameGroupBy.transform` (see the next section) will broadcast the result
across the group, producing a transformed result. If the aggregation method has an efficient
implementation, this will be performant as well.

.. _groupby.transformation.transform:

The :meth:`~.DataFrameGroupBy.transform` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Similar to the :ref:`aggregation method <groupby.aggregate.agg>`, the
:meth:`~.DataFrameGroupBy.transform` method can accept string aliases to the built-in
transformation methods in the previous section. It can *also* accept string aliases to
the built-in aggregation methods. When an aggregation method is provided, the result
will be broadcast across the group.

.. ipython:: python

    speeds
    grouped = speeds.groupby("class")[["max_speed"]]
    grouped.transform("cumsum")
    grouped.transform("sum")

In addition to string aliases, the :meth:`~.DataFrameGroupBy.transform` method can
also accept User-Defined Functions (UDFs). The UDF must:

* Return a result that is either the same size as the group chunk or
  broadcastable to the size of the group chunk (e.g., a scalar,
  ``grouped.transform(lambda x: x.iloc[-1])``).
* Operate column-by-column on the group chunk.  The transform is applied to
  the first group chunk using chunk.apply.
* Not perform in-place operations on the group chunk. Group chunks should
  be treated as immutable, and changes to a group chunk may produce unexpected
  results. See :ref:`gotchas.udf-mutation` for more information.
* (Optionally) operates on all columns of the entire group chunk at once. If this is
  supported, a fast path is used starting from the *second* chunk.

.. note::

    Transforming by supplying ``transform`` with a UDF is
    often less performant than using the built-in methods on GroupBy.
    Consider breaking up a complex operation into a chain of operations that utilize
    the built-in methods.

    All of the examples in this section can be made more performant by calling
    built-in methods instead of using UDFs.
    See :ref:`below for examples <groupby_efficient_transforms>`.

.. versionchanged:: 2.0.0

    When using ``.transform`` on a grouped DataFrame and the transformation function
    returns a DataFrame, pandas now aligns the result's index
    with the input's index. You can call ``.to_numpy()`` within the transformation
    function to avoid alignment.

Similar to :ref:`groupby.aggregate.agg`, the resulting dtype will reflect that of the
transformation function. If the results from different groups have different dtypes, then
a common dtype will be determined in the same way as ``DataFrame`` construction.

Suppose we wish to standardize the data within each group:

.. ipython:: python

   index = pd.date_range("10/1/1999", periods=1100)
   ts = pd.Series(np.random.normal(0.5, 2, 1100), index)
   ts = ts.rolling(window=100, min_periods=100).mean().dropna()

   ts.head()
   ts.tail()

   transformed = ts.groupby(lambda x: x.year).transform(
       lambda x: (x - x.mean()) / x.std()
   )


We would expect the result to now have mean 0 and standard deviation 1 within
each group (up to floating-point error), which we can easily check:

.. ipython:: python

   # Original Data
   grouped = ts.groupby(lambda x: x.year)
   grouped.mean()
   grouped.std()

   # Transformed Data
   grouped_trans = transformed.groupby(lambda x: x.year)
   grouped_trans.mean()
   grouped_trans.std()

We can also visually compare the original and transformed data sets.

.. ipython:: python

   compare = pd.DataFrame({"Original": ts, "Transformed": transformed})

   @savefig groupby_transform_plot.png
   compare.plot()

Transformation functions that have lower dimension outputs are broadcast to
match the shape of the input array.

.. ipython:: python

   ts.groupby(lambda x: x.year).transform(lambda x: x.max() - x.min())

Another common data transform is to replace missing data with the group mean.

.. ipython:: python

   cols = ["A", "B", "C"]
   values = np.random.randn(1000, 3)
   values[np.random.randint(0, 1000, 100), 0] = np.nan
   values[np.random.randint(0, 1000, 50), 1] = np.nan
   values[np.random.randint(0, 1000, 200), 2] = np.nan
   data_df = pd.DataFrame(values, columns=cols)
   data_df

   countries = np.array(["US", "UK", "GR", "JP"])
   key = countries[np.random.randint(0, 4, 1000)]

   grouped = data_df.groupby(key)

   # Non-NA count in each group
   grouped.count()

   transformed = grouped.transform(lambda x: x.fillna(x.mean()))

We can verify that the group means have not changed in the transformed data,
and that the transformed data contains no NAs.

.. ipython:: python

   grouped_trans = transformed.groupby(key)

   grouped.mean()  # original group means
   grouped_trans.mean()  # transformation did not change group means

   grouped.count()  # original has some missing data points
   grouped_trans.count()  # counts after transformation
   grouped_trans.size()  # Verify non-NA count equals group size

.. _groupby_efficient_transforms:

As mentioned in the note above, each of the examples in this section can be computed
more efficiently using built-in methods. In the code below, the inefficient way
using a UDF is commented out and the faster alternative appears below.

.. ipython:: python

    # result = ts.groupby(lambda x: x.year).transform(
    #     lambda x: (x - x.mean()) / x.std()
    # )
    grouped = ts.groupby(lambda x: x.year)
    result = (ts - grouped.transform("mean")) / grouped.transform("std")

    # result = ts.groupby(lambda x: x.year).transform(lambda x: x.max() - x.min())
    grouped = ts.groupby(lambda x: x.year)
    result = grouped.transform("max") - grouped.transform("min")

    # grouped = data_df.groupby(key)
    # result = grouped.transform(lambda x: x.fillna(x.mean()))
    grouped = data_df.groupby(key)
    result = data_df.fillna(grouped.transform("mean"))

.. _groupby.transform.window_resample:

Window and resample operations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to use ``resample()``, ``expanding()`` and
``rolling()`` as methods on groupbys.

The example below will apply the ``rolling()`` method on the samples of
the column B, based on the groups of column A.

.. ipython:: python

   df_re = pd.DataFrame({"A": [1] * 10 + [5] * 10, "B": np.arange(20)})
   df_re

   df_re.groupby("A").rolling(4).B.mean()


The ``expanding()`` method will accumulate a given operation
(``sum()`` in the example) for all the members of each particular
group.

.. ipython:: python

   df_re.groupby("A").expanding().sum()


Suppose you want to use the ``resample()`` method to get a daily
frequency in each group of your dataframe, and wish to complete the
missing values with the ``ffill()`` method.

.. ipython:: python

   df_re = pd.DataFrame(
       {
           "date": pd.date_range(start="2016-01-01", periods=4, freq="W"),
           "group": [1, 1, 2, 2],
           "val": [5, 6, 7, 8],
       }
   ).set_index("date")
   df_re

   df_re.groupby("group").resample("1D", include_groups=False).ffill()

.. _groupby.filter:

Filtration
----------

A filtration is a GroupBy operation that subsets the original grouping object. It
may either filter out entire groups, part of groups, or both. Filtrations return
a filtered version of the calling object, including the grouping columns when provided.
In the following example, ``class`` is included in the result.

.. ipython:: python

    speeds
    speeds.groupby("class").nth(1)

.. note::

    Unlike aggregations, filtrations do not add the group keys to the index of the
    result. Because of this, passing ``as_index=False`` or ``sort=True`` will not
    affect these methods.

Filtrations will respect subsetting the columns of the GroupBy object.

.. ipython:: python

    speeds.groupby("class")[["order", "max_speed"]].nth(1)

Built-in filtrations
~~~~~~~~~~~~~~~~~~~~

The following methods on GroupBy act as filtrations. All these methods have an
efficient, GroupBy-specific, implementation.

.. csv-table::
    :header: "Method", "Description"
    :widths: 20, 80
    :delim: ;

        :meth:`~.DataFrameGroupBy.head`;Select the top row(s) of each group
        :meth:`~.DataFrameGroupBy.nth`;Select the nth row(s) of each group
        :meth:`~.DataFrameGroupBy.tail`;Select the bottom row(s) of each group

Users can also use transformations along with Boolean indexing to construct complex
filtrations within groups. For example, suppose we are given groups of products and
their volumes, and we wish to subset the data to only the largest products capturing no
more than 90% of the total volume within each group.

.. ipython:: python

    product_volumes = pd.DataFrame(
        {
            "group": list("xxxxyyy"),
            "product": list("abcdefg"),
            "volume": [10, 30, 20, 15, 40, 10, 20],
        }
    )
    product_volumes

    # Sort by volume to select the largest products first
    product_volumes = product_volumes.sort_values("volume", ascending=False)
    grouped = product_volumes.groupby("group")["volume"]
    cumpct = grouped.cumsum() / grouped.transform("sum")
    cumpct
    significant_products = product_volumes[cumpct <= 0.9]
    significant_products.sort_values(["group", "product"])

The :class:`~DataFrameGroupBy.filter` method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    Filtering by supplying ``filter`` with a User-Defined Function (UDF) is
    often less performant than using the built-in methods on GroupBy.
    Consider breaking up a complex operation into a chain of operations that utilize
    the built-in methods.

The ``filter`` method takes a User-Defined Function (UDF) that, when applied to
an entire group, returns either ``True`` or ``False``. The result of the ``filter``
method is then the subset of groups for which the UDF returned ``True``.

Suppose we want to take only elements that belong to groups with a group sum greater
than 2.

.. ipython:: python

   sf = pd.Series([1, 1, 2, 3, 3, 3])
   sf.groupby(sf).filter(lambda x: x.sum() > 2)

Another useful operation is filtering out elements that belong to groups
with only a couple members.

.. ipython:: python

   dff = pd.DataFrame({"A": np.arange(8), "B": list("aabbbbcc")})
   dff.groupby("B").filter(lambda x: len(x) > 2)

Alternatively, instead of dropping the offending groups, we can return a
like-indexed objects where the groups that do not pass the filter are filled
with NaNs.

.. ipython:: python

   dff.groupby("B").filter(lambda x: len(x) > 2, dropna=False)

For DataFrames with multiple columns, filters should explicitly specify a column as the filter criterion.

.. ipython:: python

   dff["C"] = np.arange(8)
   dff.groupby("B").filter(lambda x: len(x["C"]) > 2)

.. _groupby.apply:

Flexible ``apply``
------------------

Some operations on the grouped data might not fit into the aggregation,
transformation, or filtration categories. For these, you can use the ``apply``
function.

.. warning::

   ``apply`` has to try to infer from the result whether it should act as a reducer,
   transformer, *or* filter, depending on exactly what is passed to it. Thus the
   grouped column(s) may be included in the output or not. While
   it tries to intelligently guess how to behave, it can sometimes guess wrong.

.. note::

   All of the examples in this section can be more reliably, and more efficiently,
   computed using other pandas functionality.

.. ipython:: python

   df
   grouped = df.groupby("A")

   # could also just call .describe()
   grouped["C"].apply(lambda x: x.describe())

The dimension of the returned result can also change:

.. ipython:: python

    grouped = df.groupby('A')['C']

    def f(group):
        return pd.DataFrame({'original': group,
                             'demeaned': group - group.mean()})

    grouped.apply(f)

``apply`` on a Series can operate on a returned value from the applied function
that is itself a series, and possibly upcast the result to a DataFrame:

.. ipython:: python

    def f(x):
        return pd.Series([x, x ** 2], index=["x", "x^2"])


    s = pd.Series(np.random.rand(5))
    s
    s.apply(f)

Similar to :ref:`groupby.aggregate.agg`, the resulting dtype will reflect that of the
apply function. If the results from different groups have different dtypes, then
a common dtype will be determined in the same way as ``DataFrame`` construction.

Control grouped column(s) placement with ``group_keys``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To control whether the grouped column(s) are included in the indices, you can use
the argument ``group_keys`` which defaults to ``True``. Compare

.. ipython:: python

    df.groupby("A", group_keys=True).apply(lambda x: x, include_groups=False)

with

.. ipython:: python

    df.groupby("A", group_keys=False).apply(lambda x: x, include_groups=False)


Numba Accelerated Routines
--------------------------

.. versionadded:: 1.1

If `Numba <https://numba.pydata.org/>`__ is installed as an optional dependency, the ``transform`` and
``aggregate`` methods support ``engine='numba'`` and ``engine_kwargs`` arguments.
See :ref:`enhancing performance with Numba <enhancingperf.numba>` for general usage of the arguments
and performance considerations.

The function signature must start with ``values, index`` **exactly** as the data belonging to each group
will be passed into ``values``, and the group index will be passed into ``index``.

.. warning::

   When using ``engine='numba'``, there will be no "fall back" behavior internally. The group
   data and group index will be passed as NumPy arrays to the JITed user defined function, and no
   alternative execution attempts will be tried.

Other useful features
---------------------

Exclusion of non-numeric columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Again consider the example DataFrame we've been looking at:

.. ipython:: python

   df

Suppose we wish to compute the standard deviation grouped by the ``A``
column. There is a slight problem, namely that we don't care about the data in
column ``B`` because it is not numeric. You can avoid non-numeric columns by
specifying ``numeric_only=True``:

.. ipython:: python

   df.groupby("A").std(numeric_only=True)

Note that ``df.groupby('A').colname.std().`` is more efficient than
``df.groupby('A').std().colname``. So if the result of an aggregation function
is only needed over one column (here ``colname``), it may be filtered
*before* applying the aggregation function.

.. ipython:: python

    from decimal import Decimal

    df_dec = pd.DataFrame(
        {
            "id": [1, 2, 1, 2],
            "int_column": [1, 2, 3, 4],
            "dec_column": [
                Decimal("0.50"),
                Decimal("0.15"),
                Decimal("0.25"),
                Decimal("0.40"),
            ],
        }
    )
    df_dec.groupby(["id"])[["dec_column"]].sum()

.. _groupby.observed:

Handling of (un)observed Categorical values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using a ``Categorical`` grouper (as a single grouper, or as part of multiple groupers), the ``observed`` keyword
controls whether to return a cartesian product of all possible groupers values (``observed=False``) or only those
that are observed groupers (``observed=True``).

Show all values:

.. ipython:: python

   pd.Series([1, 1, 1]).groupby(
       pd.Categorical(["a", "a", "a"], categories=["a", "b"]), observed=False
   ).count()

Show only the observed values:

.. ipython:: python

   pd.Series([1, 1, 1]).groupby(
       pd.Categorical(["a", "a", "a"], categories=["a", "b"]), observed=True
   ).count()

The returned dtype of the grouped will *always* include *all* of the categories that were grouped.

.. ipython:: python

   s = (
       pd.Series([1, 1, 1])
       .groupby(pd.Categorical(["a", "a", "a"], categories=["a", "b"]), observed=True)
       .count()
   )
   s.index.dtype

.. _groupby.missing:

NA group handling
~~~~~~~~~~~~~~~~~

By ``NA``, we are referring to any ``NA`` values, including
:class:`NA`, ``NaN``, ``NaT``, and ``None``. If there are any ``NA`` values in the
grouping key, by default these will be excluded. In other words, any
"``NA`` group" will be dropped. You can include NA groups by specifying ``dropna=False``.

.. ipython:: python

   df = pd.DataFrame({"key": [1.0, 1.0, np.nan, 2.0, np.nan], "A": [1, 2, 3, 4, 5]})
   df

   df.groupby("key", dropna=True).sum()

   df.groupby("key", dropna=False).sum()

Grouping with ordered factors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Categorical variables represented as instances of pandas's ``Categorical`` class
can be used as group keys. If so, the order of the levels will be preserved. When
``observed=False`` and ``sort=False``, any unobserved categories will be at the
end of the result in order.

.. ipython:: python

    days = pd.Categorical(
        values=["Wed", "Mon", "Thu", "Mon", "Wed", "Sat"],
        categories=["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"],
    )
    data = pd.DataFrame(
       {
           "day": days,
           "workers": [3, 4, 1, 4, 2, 2],
       }
    )
    data

    data.groupby("day", observed=False, sort=True).sum()

    data.groupby("day", observed=False, sort=False).sum()

.. _groupby.specify:

Grouping with a grouper specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may need to specify a bit more data to properly group. You can
use the ``pd.Grouper`` to provide this local control.

.. ipython:: python

   import datetime

   df = pd.DataFrame(
       {
           "Branch": "A A A A A A A B".split(),
           "Buyer": "Carl Mark Carl Carl Joe Joe Joe Carl".split(),
           "Quantity": [1, 3, 5, 1, 8, 1, 9, 3],
           "Date": [
               datetime.datetime(2013, 1, 1, 13, 0),
               datetime.datetime(2013, 1, 1, 13, 5),
               datetime.datetime(2013, 10, 1, 20, 0),
               datetime.datetime(2013, 10, 2, 10, 0),
               datetime.datetime(2013, 10, 1, 20, 0),
               datetime.datetime(2013, 10, 2, 10, 0),
               datetime.datetime(2013, 12, 2, 12, 0),
               datetime.datetime(2013, 12, 2, 14, 0),
           ],
       }
   )

   df

Groupby a specific column with the desired frequency. This is like resampling.

.. ipython:: python

   df.groupby([pd.Grouper(freq="1ME", key="Date"), "Buyer"])[["Quantity"]].sum()

When ``freq`` is specified, the object returned by ``pd.Grouper`` will be an
instance of ``pandas.api.typing.TimeGrouper``. When there is a column and index
with the same name, you can use ``key`` to group by the column and ``level``
to group by the index.

.. ipython:: python

   df = df.set_index("Date")
   df["Date"] = df.index + pd.offsets.MonthEnd(2)
   df.groupby([pd.Grouper(freq="6ME", key="Date"), "Buyer"])[["Quantity"]].sum()

   df.groupby([pd.Grouper(freq="6ME", level="Date"), "Buyer"])[["Quantity"]].sum()


Taking the first rows of each group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Just like for a DataFrame or Series you can call head and tail on a groupby:

.. ipython:: python

   df = pd.DataFrame([[1, 2], [1, 4], [5, 6]], columns=["A", "B"])
   df

   g = df.groupby("A")
   g.head(1)

   g.tail(1)

This shows the first or last n rows from each group.

.. _groupby.nth:

Taking the nth row of each group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To select the nth item from each group, use :meth:`.DataFrameGroupBy.nth` or
:meth:`.SeriesGroupBy.nth`. Arguments supplied can be any integer, lists of integers,
slices, or lists of slices; see below for examples. When the nth element of a group
does not exist an error is *not* raised; instead no corresponding rows are returned.

In general this operation acts as a filtration. In certain cases it will also return
one row per group, making it also a reduction. However because in general it can
return zero or multiple rows per group, pandas treats it as a filtration in all cases.

.. ipython:: python

   df = pd.DataFrame([[1, np.nan], [1, 4], [5, 6]], columns=["A", "B"])
   g = df.groupby("A")

   g.nth(0)
   g.nth(-1)
   g.nth(1)

If the nth element of a group does not exist, then no corresponding row is included
in the result. In particular, if the specified ``n`` is larger than any group, the
result will be an empty DataFrame.

.. ipython:: python

   g.nth(5)

If you want to select the nth not-null item, use the ``dropna`` kwarg. For a DataFrame this should be either ``'any'`` or ``'all'`` just like you would pass to dropna:

.. ipython:: python

   # nth(0) is the same as g.first()
   g.nth(0, dropna="any")
   g.first()

   # nth(-1) is the same as g.last()
   g.nth(-1, dropna="any")
   g.last()

   g.B.nth(0, dropna="all")

You can also select multiple rows from each group by specifying multiple nth values as a list of ints.

.. ipython:: python

   business_dates = pd.date_range(start="4/1/2014", end="6/30/2014", freq="B")
   df = pd.DataFrame(1, index=business_dates, columns=["a", "b"])
   # get the first, 4th, and last date index for each month
   df.groupby([df.index.year, df.index.month]).nth([0, 3, -1])

You may also use slices or lists of slices.

.. ipython:: python

   df.groupby([df.index.year, df.index.month]).nth[1:]
   df.groupby([df.index.year, df.index.month]).nth[1:, :-1]

Enumerate group items
~~~~~~~~~~~~~~~~~~~~~

To see the order in which each row appears within its group, use the
``cumcount`` method:

.. ipython:: python

   dfg = pd.DataFrame(list("aaabba"), columns=["A"])
   dfg

   dfg.groupby("A").cumcount()

   dfg.groupby("A").cumcount(ascending=False)

.. _groupby.ngroup:

Enumerate groups
~~~~~~~~~~~~~~~~

To see the ordering of the groups (as opposed to the order of rows
within a group given by ``cumcount``) you can use
:meth:`.DataFrameGroupBy.ngroup`.



Note that the numbers given to the groups match the order in which the
groups would be seen when iterating over the groupby object, not the
order they are first observed.

.. ipython:: python

   dfg = pd.DataFrame(list("aaabba"), columns=["A"])
   dfg

   dfg.groupby("A").ngroup()

   dfg.groupby("A").ngroup(ascending=False)

Plotting
~~~~~~~~

Groupby also works with some plotting methods.  In this case, suppose we
suspect that the values in column 1 are 3 times higher on average in group "B".


.. ipython:: python

   np.random.seed(1234)
   df = pd.DataFrame(np.random.randn(50, 2))
   df["g"] = np.random.choice(["A", "B"], size=50)
   df.loc[df["g"] == "B", 1] += 3

We can easily visualize this with a boxplot:

.. ipython:: python
   :okwarning:

   @savefig groupby_boxplot.png
   df.groupby("g").boxplot()

The result of calling ``boxplot`` is a dictionary whose keys are the values
of our grouping column ``g`` ("A" and "B"). The values of the resulting dictionary
can be controlled by the ``return_type`` keyword of ``boxplot``.
See the :ref:`visualization documentation<visualization.box>` for more.

.. warning::

  For historical reasons, ``df.groupby("g").boxplot()`` is not equivalent
  to ``df.boxplot(by="g")``. See :ref:`here<visualization.box.return>` for
  an explanation.

.. _groupby.pipe:

Piping function calls
~~~~~~~~~~~~~~~~~~~~~

Similar to the functionality provided by ``DataFrame`` and ``Series``, functions
that take ``GroupBy`` objects can be chained together using a ``pipe`` method to
allow for a cleaner, more readable syntax. To read about ``.pipe`` in general terms,
see :ref:`here <basics.pipe>`.

Combining ``.groupby`` and ``.pipe`` is often useful when you need to reuse
GroupBy objects.

As an example, imagine having a DataFrame with columns for stores, products,
revenue and quantity sold. We'd like to do a groupwise calculation of *prices*
(i.e. revenue/quantity) per store and per product. We could do this in a
multi-step operation, but expressing it in terms of piping can make the
code more readable. First we set the data:

.. ipython:: python

   n = 1000
   df = pd.DataFrame(
       {
           "Store": np.random.choice(["Store_1", "Store_2"], n),
           "Product": np.random.choice(["Product_1", "Product_2"], n),
           "Revenue": (np.random.random(n) * 50 + 10).round(2),
           "Quantity": np.random.randint(1, 10, size=n),
       }
   )
   df.head(2)

We now find the prices per store/product.

.. ipython:: python

   (
       df.groupby(["Store", "Product"])
       .pipe(lambda grp: grp.Revenue.sum() / grp.Quantity.sum())
       .unstack()
       .round(2)
   )

Piping can also be expressive when you want to deliver a grouped object to some
arbitrary function, for example:

.. ipython:: python

   def mean(groupby):
       return groupby.mean()


   df.groupby(["Store", "Product"]).pipe(mean)

Here ``mean`` takes a GroupBy object and finds the mean of the Revenue and Quantity
columns respectively for each Store-Product combination. The ``mean`` function can
be any function that takes in a GroupBy object; the ``.pipe`` will pass the GroupBy
object as a parameter into the function you specify.

Examples
--------

.. _groupby.multicolumn_factorization:

Multi-column factorization
~~~~~~~~~~~~~~~~~~~~~~~~~~

By using :meth:`.DataFrameGroupBy.ngroup`, we can extract
information about the groups in a way similar to :func:`factorize` (as described
further in the :ref:`reshaping API <reshaping.factorize>`) but which applies
naturally to multiple columns of mixed type and different
sources. This can be useful as an intermediate categorical-like step
in processing, when the relationships between the group rows are more
important than their content, or as input to an algorithm which only
accepts the integer encoding. (For more information about support in
pandas for full categorical data, see the :ref:`Categorical
introduction <categorical>` and the
:ref:`API documentation <api.arrays.categorical>`.)

.. ipython:: python

    dfg = pd.DataFrame({"A": [1, 1, 2, 3, 2], "B": list("aaaba")})

    dfg

    dfg.groupby(["A", "B"]).ngroup()

    dfg.groupby(["A", [0, 0, 0, 1, 1]]).ngroup()

Groupby by indexer to 'resample' data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Resampling produces new hypothetical samples (resamples) from already existing observed data or from a model that generates data. These new samples are similar to the pre-existing samples.

In order for resample to work on indices that are non-datetimelike, the following procedure can be utilized.

In the following examples, **df.index // 5** returns an integer array which is used to determine what gets selected for the groupby operation.

.. note::

   The example below shows how we can downsample by consolidation of samples into fewer ones.
   Here by using **df.index // 5**, we are aggregating the samples in bins. By applying **std()**
   function, we aggregate the information contained in many samples into a small subset of values
   which is their standard deviation thereby reducing the number of samples.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 2))
   df
   df.index // 5
   df.groupby(df.index // 5).std()

Returning a Series to propagate names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Group DataFrame columns, compute a set of metrics and return a named Series.
The Series name is used as the name for the column index. This is especially
useful in conjunction with reshaping operations such as stacking, in which the
column index name will be used as the name of the inserted column:

.. ipython:: python

   df = pd.DataFrame(
       {
           "a": [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2],
           "b": [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1],
           "c": [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
           "d": [0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1],
       }
   )

   def compute_metrics(x):
       result = {"b_sum": x["b"].sum(), "c_mean": x["c"].mean()}
       return pd.Series(result, name="metrics")

   result = df.groupby("a").apply(compute_metrics, include_groups=False)

   result

   result.stack(future_stack=True)
