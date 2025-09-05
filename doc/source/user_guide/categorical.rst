.. _categorical:

{{ header }}

****************
Categorical data
****************

This is an introduction to pandas categorical data type, including a short comparison
with R's ``factor``.

``Categoricals`` are a pandas data type corresponding to categorical variables in
statistics. A categorical variable takes on a limited, and usually fixed,
number of possible values (``categories``; ``levels`` in R). Examples are gender,
social class, blood type, country affiliation, observation time or rating via
Likert scales.

In contrast to statistical categorical variables, categorical data might have an order (e.g.
'strongly agree' vs 'agree' or 'first observation' vs. 'second observation'), but numerical
operations (additions, divisions, ...) are not possible.

All values of categorical data are either in ``categories`` or ``np.nan``. Order is defined by
the order of ``categories``, not lexical order of the values. Internally, the data structure
consists of a ``categories`` array and an integer array of ``codes`` which point to the real value in
the ``categories`` array.

The categorical data type is useful in the following cases:

* A string variable consisting of only a few different values. Converting such a string
  variable to a categorical variable will save some memory, see :ref:`here <categorical.memory>`.
* The lexical order of a variable is not the same as the logical order ("one", "two", "three").
  By converting to a categorical and specifying an order on the categories, sorting and
  min/max will use the logical order instead of the lexical order, see :ref:`here <categorical.sort>`.
* As a signal to other Python libraries that this column should be treated as a categorical
  variable (e.g. to use suitable statistical methods or plot types).

See also the :ref:`API docs on categoricals<api.arrays.categorical>`.

.. _categorical.objectcreation:

Object creation
---------------

Series creation
~~~~~~~~~~~~~~~

Categorical ``Series`` or columns in a ``DataFrame`` can be created in several ways:

By specifying ``dtype="category"`` when constructing a ``Series``:

.. ipython:: python

    s = pd.Series(["a", "b", "c", "a"], dtype="category")
    s

By converting an existing ``Series`` or column to a ``category`` dtype:

.. ipython:: python

    df = pd.DataFrame({"A": ["a", "b", "c", "a"]})
    df["B"] = df["A"].astype("category")
    df

By using special functions, such as :func:`~pandas.cut`, which groups data into
discrete bins. See the :ref:`example on tiling <reshaping.tile.cut>` in the docs.

.. ipython:: python

    df = pd.DataFrame({"value": np.random.randint(0, 100, 20)})
    labels = ["{0} - {1}".format(i, i + 9) for i in range(0, 100, 10)]

    df["group"] = pd.cut(df.value, range(0, 105, 10), right=False, labels=labels)
    df.head(10)

By passing a :class:`pandas.Categorical` object to a ``Series`` or assigning it to a ``DataFrame``.

.. ipython:: python

    raw_cat = pd.Categorical(
        ["a", "b", "c", "a"], categories=["b", "c", "d"], ordered=False
    )
    s = pd.Series(raw_cat)
    s
    df = pd.DataFrame({"A": ["a", "b", "c", "a"]})
    df["B"] = raw_cat
    df

Categorical data has a specific ``category`` :ref:`dtype <basics.dtypes>`:

.. ipython:: python

    df.dtypes

DataFrame creation
~~~~~~~~~~~~~~~~~~

Similar to the previous section where a single column was converted to categorical, all columns in a
``DataFrame`` can be batch converted to categorical either during or after construction.

This can be done during construction by specifying ``dtype="category"`` in the ``DataFrame`` constructor:

.. ipython:: python

    df = pd.DataFrame({"A": list("abca"), "B": list("bccd")}, dtype="category")
    df.dtypes

Note that the categories present in each column differ; the conversion is done column by column, so
only labels present in a given column are categories:

.. ipython:: python

    df["A"]
    df["B"]


Analogously, all columns in an existing ``DataFrame`` can be batch converted using :meth:`DataFrame.astype`:

.. ipython:: python

    df = pd.DataFrame({"A": list("abca"), "B": list("bccd")})
    df_cat = df.astype("category")
    df_cat.dtypes

This conversion is likewise done column by column:

.. ipython:: python

    df_cat["A"]
    df_cat["B"]


Controlling behavior
~~~~~~~~~~~~~~~~~~~~

In the examples above where we passed ``dtype='category'``, we used the default
behavior:

1. Categories are inferred from the data.
2. Categories are unordered.

To control those behaviors, instead of passing ``'category'``, use an instance
of :class:`~pandas.api.types.CategoricalDtype`.

.. ipython:: python

    from pandas.api.types import CategoricalDtype

    s = pd.Series(["a", "b", "c", "a"])
    cat_type = CategoricalDtype(categories=["b", "c", "d"], ordered=True)
    s_cat = s.astype(cat_type)
    s_cat

Similarly, a ``CategoricalDtype`` can be used with a ``DataFrame`` to ensure that categories
are consistent among all columns.

.. ipython:: python

    from pandas.api.types import CategoricalDtype

    df = pd.DataFrame({"A": list("abca"), "B": list("bccd")})
    cat_type = CategoricalDtype(categories=list("abcd"), ordered=True)
    df_cat = df.astype(cat_type)
    df_cat["A"]
    df_cat["B"]

.. note::

    To perform table-wise conversion, where all labels in the entire ``DataFrame`` are used as
    categories for each column, the ``categories`` parameter can be determined programmatically by
    ``categories = pd.unique(df.to_numpy().ravel())``.

If you already have ``codes`` and ``categories``, you can use the
:func:`~pandas.Categorical.from_codes` constructor to save the factorize step
during normal constructor mode:

.. ipython:: python

    splitter = np.random.choice([0, 1], 5, p=[0.5, 0.5])
    s = pd.Series(pd.Categorical.from_codes(splitter, categories=["train", "test"]))


Regaining original data
~~~~~~~~~~~~~~~~~~~~~~~

To get back to the original ``Series`` or NumPy array, use
``Series.astype(original_dtype)`` or ``np.asarray(categorical)``:

.. ipython:: python

    s = pd.Series(["a", "b", "c", "a"])
    s
    s2 = s.astype("category")
    s2
    s2.astype(str)
    np.asarray(s2)

.. note::

    In contrast to R's ``factor`` function, categorical data is not converting input values to
    strings; categories will end up the same data type as the original values.

.. note::

    In contrast to R's ``factor`` function, there is currently no way to assign/change labels at
    creation time. Use ``categories`` to change the categories after creation time.

.. _categorical.categoricaldtype:

CategoricalDtype
----------------

A categorical's type is fully described by

1. ``categories``: a sequence of unique values and no missing values
2. ``ordered``: a boolean

This information can be stored in a :class:`~pandas.api.types.CategoricalDtype`.
The ``categories`` argument is optional, which implies that the actual categories
should be inferred from whatever is present in the data when the
:class:`pandas.Categorical` is created. The categories are assumed to be unordered
by default.

.. ipython:: python

   from pandas.api.types import CategoricalDtype

   CategoricalDtype(["a", "b", "c"])
   CategoricalDtype(["a", "b", "c"], ordered=True)
   CategoricalDtype()

A :class:`~pandas.api.types.CategoricalDtype` can be used in any place pandas
expects a ``dtype``. For example :func:`pandas.read_csv`,
:func:`pandas.DataFrame.astype`, or in the ``Series`` constructor.

.. note::

    As a convenience, you can use the string ``'category'`` in place of a
    :class:`~pandas.api.types.CategoricalDtype` when you want the default behavior of
    the categories being unordered, and equal to the set values present in the
    array. In other words, ``dtype='category'`` is equivalent to
    ``dtype=CategoricalDtype()``.

Equality semantics
~~~~~~~~~~~~~~~~~~

Two instances of :class:`~pandas.api.types.CategoricalDtype` compare equal
whenever they have the same categories and order. When comparing two
unordered categoricals, the order of the ``categories`` is not considered. Note
that categories with different dtypes are not the same.

.. ipython:: python

   c1 = CategoricalDtype(["a", "b", "c"], ordered=False)

   # Equal, since order is not considered when ordered=False
   c1 == CategoricalDtype(["b", "c", "a"], ordered=False)

   # Unequal, since the second CategoricalDtype is ordered
   c1 == CategoricalDtype(["a", "b", "c"], ordered=True)

All instances of ``CategoricalDtype`` compare equal to the string ``'category'``.

.. ipython:: python

   c1 == "category"

Notice that the ``categories_dtype`` should be considered, especially when comparing with
two empty ``CategoricalDtype`` instances.

.. ipython:: python

    c2 = pd.Categorical(np.array([], dtype=object))
    c3 = pd.Categorical(np.array([], dtype=float))

    c2.dtype == c3.dtype

Description
-----------

Using :meth:`~DataFrame.describe` on categorical data will produce similar
output to a ``Series`` or ``DataFrame`` of type ``string``.

.. ipython:: python

    cat = pd.Categorical(["a", "c", "c", np.nan], categories=["b", "a", "c"])
    df = pd.DataFrame({"cat": cat, "s": ["a", "c", "c", np.nan]})
    df.describe()
    df["cat"].describe()

.. _categorical.cat:

Working with categories
-----------------------

Categorical data has a ``categories`` and a ``ordered`` property, which list their
possible values and whether the ordering matters or not. These properties are
exposed as ``s.cat.categories`` and ``s.cat.ordered``. If you don't manually
specify categories and ordering, they are inferred from the passed arguments.

.. ipython:: python

    s = pd.Series(["a", "b", "c", "a"], dtype="category")
    s.cat.categories
    s.cat.ordered

It's also possible to pass in the categories in a specific order:

.. ipython:: python

    s = pd.Series(pd.Categorical(["a", "b", "c", "a"], categories=["c", "b", "a"]))
    s.cat.categories
    s.cat.ordered

.. note::

    New categorical data are **not** automatically ordered. You must explicitly
    pass ``ordered=True`` to indicate an ordered ``Categorical``.


.. note::

    The result of :meth:`~Series.unique` is not always the same as ``Series.cat.categories``,
    because ``Series.unique()`` has a couple of guarantees, namely that it returns categories
    in the order of appearance, and it only includes values that are actually present.

    .. ipython:: python

         s = pd.Series(list("babc")).astype(CategoricalDtype(list("abcd")))
         s

         # categories
         s.cat.categories

         # uniques
         s.unique()

Renaming categories
~~~~~~~~~~~~~~~~~~~

Renaming categories is done by using the
:meth:`~pandas.Categorical.rename_categories` method:


.. ipython:: python

    s = pd.Series(["a", "b", "c", "a"], dtype="category")
    s
    new_categories = ["Group %s" % g for g in s.cat.categories]
    s = s.cat.rename_categories(new_categories)
    s
    # You can also pass a dict-like object to map the renaming
    s = s.cat.rename_categories({1: "x", 2: "y", 3: "z"})
    s

.. note::

    In contrast to R's ``factor``, categorical data can have categories of other types than string.

Categories must be unique or a ``ValueError`` is raised:

.. ipython:: python

    try:
        s = s.cat.rename_categories([1, 1, 1])
    except ValueError as e:
        print("ValueError:", str(e))

Categories must also not be ``NaN`` or a ``ValueError`` is raised:

.. ipython:: python

    try:
        s = s.cat.rename_categories([1, 2, np.nan])
    except ValueError as e:
        print("ValueError:", str(e))

Appending new categories
~~~~~~~~~~~~~~~~~~~~~~~~

Appending categories can be done by using the
:meth:`~pandas.Categorical.add_categories` method:

.. ipython:: python

    s = s.cat.add_categories([4])
    s.cat.categories
    s

Removing categories
~~~~~~~~~~~~~~~~~~~

Removing categories can be done by using the
:meth:`~pandas.Categorical.remove_categories` method. Values which are removed
are replaced by ``np.nan``.:

.. ipython:: python

    s = s.cat.remove_categories([4])
    s

Removing unused categories
~~~~~~~~~~~~~~~~~~~~~~~~~~

Removing unused categories can also be done:

.. ipython:: python

    s = pd.Series(pd.Categorical(["a", "b", "a"], categories=["a", "b", "c", "d"]))
    s
    s.cat.remove_unused_categories()

Setting categories
~~~~~~~~~~~~~~~~~~

If you want to do remove and add new categories in one step (which has some
speed advantage), or simply set the categories to a predefined scale,
use :meth:`~pandas.Categorical.set_categories`.


.. ipython:: python

    s = pd.Series(["one", "two", "four", "-"], dtype="category")
    s
    s = s.cat.set_categories(["one", "two", "three", "four"])
    s

.. note::
    Be aware that :func:`Categorical.set_categories` cannot know whether some category is omitted
    intentionally or because it is misspelled or (under Python3) due to a type difference (e.g.,
    NumPy S1 dtype and Python strings). This can result in surprising behaviour!

Sorting and order
-----------------

.. _categorical.sort:

If categorical data is ordered (``s.cat.ordered == True``), then the order of the categories has a
meaning and certain operations are possible. If the categorical is unordered, ``.min()/.max()`` will raise a ``TypeError``.

.. ipython:: python

    s = pd.Series(pd.Categorical(["a", "b", "c", "a"], ordered=False))
    s = s.sort_values()
    s = pd.Series(["a", "b", "c", "a"]).astype(CategoricalDtype(ordered=True))
    s = s.sort_values()
    s
    s.min(), s.max()

You can set categorical data to be ordered by using ``as_ordered()`` or unordered by using ``as_unordered()``. These will by
default return a *new* object.

.. ipython:: python

    s.cat.as_ordered()
    s.cat.as_unordered()

Sorting will use the order defined by categories, not any lexical order present on the data type.
This is even true for strings and numeric data:

.. ipython:: python

    s = pd.Series([1, 2, 3, 1], dtype="category")
    s = s.cat.set_categories([2, 3, 1], ordered=True)
    s
    s = s.sort_values()
    s
    s.min(), s.max()


Reordering
~~~~~~~~~~

Reordering the categories is possible via the :meth:`Categorical.reorder_categories` and
the :meth:`Categorical.set_categories` methods. For :meth:`Categorical.reorder_categories`, all
old categories must be included in the new categories and no new categories are allowed. This will
necessarily make the sort order the same as the categories order.

.. ipython:: python

    s = pd.Series([1, 2, 3, 1], dtype="category")
    s = s.cat.reorder_categories([2, 3, 1], ordered=True)
    s
    s = s.sort_values()
    s
    s.min(), s.max()

.. note::

    Note the difference between assigning new categories and reordering the categories: the first
    renames categories and therefore the individual values in the ``Series``, but if the first
    position was sorted last, the renamed value will still be sorted last. Reordering means that the
    way values are sorted is different afterwards, but not that individual values in the
    ``Series`` are changed.

.. note::

    If the ``Categorical`` is not ordered, :meth:`Series.min` and :meth:`Series.max` will raise
    ``TypeError``. Numeric operations like ``+``, ``-``, ``*``, ``/`` and operations based on them
    (e.g. :meth:`Series.median`, which would need to compute the mean between two values if the length
    of an array is even) do not work and raise a ``TypeError``.

Multi column sorting
~~~~~~~~~~~~~~~~~~~~

A categorical dtyped column will participate in a multi-column sort in a similar manner to other columns.
The ordering of the categorical is determined by the ``categories`` of that column.

.. ipython:: python

   dfs = pd.DataFrame(
       {
           "A": pd.Categorical(
               list("bbeebbaa"),
               categories=["e", "a", "b"],
               ordered=True,
           ),
           "B": [1, 2, 1, 2, 2, 1, 2, 1],
       }
   )
   dfs.sort_values(by=["A", "B"])

Reordering the ``categories`` changes a future sort.

.. ipython:: python

   dfs["A"] = dfs["A"].cat.reorder_categories(["a", "b", "e"])
   dfs.sort_values(by=["A", "B"])

Comparisons
-----------

Comparing categorical data with other objects is possible in three cases:

* Comparing equality (``==`` and ``!=``) to a list-like object (list, Series, array,
  ...) of the same length as the categorical data.
* All comparisons (``==``, ``!=``, ``>``, ``>=``, ``<``, and ``<=``) of categorical data to
  another categorical Series, when ``ordered==True`` and the ``categories`` are the same.
* All comparisons of a categorical data to a scalar.

All other comparisons, especially "non-equality" comparisons of two categoricals with different
categories or a categorical with any list-like object, will raise a ``TypeError``.

.. note::

    Any "non-equality" comparisons of categorical data with a ``Series``, ``np.array``, ``list`` or
    categorical data with different categories or ordering will raise a ``TypeError`` because custom
    categories ordering could be interpreted in two ways: one with taking into account the
    ordering and one without.

.. ipython:: python

    cat = pd.Series([1, 2, 3]).astype(CategoricalDtype([3, 2, 1], ordered=True))
    cat_base = pd.Series([2, 2, 2]).astype(CategoricalDtype([3, 2, 1], ordered=True))
    cat_base2 = pd.Series([2, 2, 2]).astype(CategoricalDtype(ordered=True))

    cat
    cat_base
    cat_base2

Comparing to a categorical with the same categories and ordering or to a scalar works:

.. ipython:: python

    cat > cat_base
    cat > 2

Equality comparisons work with any list-like object of same length and scalars:

.. ipython:: python

    cat == cat_base
    cat == np.array([1, 2, 3])
    cat == 2

This doesn't work because the categories are not the same:

.. ipython:: python

    try:
        cat > cat_base2
    except TypeError as e:
        print("TypeError:", str(e))

If you want to do a "non-equality" comparison of a categorical series with a list-like object
which is not categorical data, you need to be explicit and convert the categorical data back to
the original values:

.. ipython:: python

    base = np.array([1, 2, 3])

    try:
        cat > base
    except TypeError as e:
        print("TypeError:", str(e))

    np.asarray(cat) > base

When you compare two unordered categoricals with the same categories, the order is not considered:

.. ipython:: python

   c1 = pd.Categorical(["a", "b"], categories=["a", "b"], ordered=False)
   c2 = pd.Categorical(["a", "b"], categories=["b", "a"], ordered=False)
   c1 == c2

Operations
----------

Apart from :meth:`Series.min`, :meth:`Series.max` and :meth:`Series.mode`, the
following operations are possible with categorical data:

``Series`` methods like :meth:`Series.value_counts` will use all categories,
even if some categories are not present in the data:

.. ipython:: python

    s = pd.Series(pd.Categorical(["a", "b", "c", "c"], categories=["c", "a", "b", "d"]))
    s.value_counts()

``DataFrame`` methods like :meth:`DataFrame.sum` also show "unused" categories when ``observed=False``.

.. ipython:: python

    columns = pd.Categorical(
        ["One", "One", "Two"], categories=["One", "Two", "Three"], ordered=True
    )
    df = pd.DataFrame(
        data=[[1, 2, 3], [4, 5, 6]],
        columns=pd.MultiIndex.from_arrays([["A", "B", "B"], columns]),
    ).T
    df.groupby(level=1, observed=False).sum()

Groupby will also show "unused" categories when ``observed=False``:

.. ipython:: python

    cats = pd.Categorical(
        ["a", "b", "b", "b", "c", "c", "c"], categories=["a", "b", "c", "d"]
    )
    df = pd.DataFrame({"cats": cats, "values": [1, 2, 2, 2, 3, 4, 5]})
    df.groupby("cats", observed=False).mean()

    cats2 = pd.Categorical(["a", "a", "b", "b"], categories=["a", "b", "c"])
    df2 = pd.DataFrame(
        {
            "cats": cats2,
            "B": ["c", "d", "c", "d"],
            "values": [1, 2, 3, 4],
        }
    )
    df2.groupby(["cats", "B"], observed=False).mean()


Pivot tables:

.. ipython:: python

    raw_cat = pd.Categorical(["a", "a", "b", "b"], categories=["a", "b", "c"])
    df = pd.DataFrame({"A": raw_cat, "B": ["c", "d", "c", "d"], "values": [1, 2, 3, 4]})
    pd.pivot_table(df, values="values", index=["A", "B"], observed=False)

Data munging
------------

The optimized pandas data access methods  ``.loc``, ``.iloc``, ``.at``, and ``.iat``,
work as normal. The only difference is the return type (for getting) and
that only values already in ``categories`` can be assigned.

Getting
~~~~~~~

If the slicing operation returns either a ``DataFrame`` or a column of type
``Series``, the ``category`` dtype is preserved.

.. ipython:: python

    idx = pd.Index(["h", "i", "j", "k", "l", "m", "n"])
    cats = pd.Series(["a", "b", "b", "b", "c", "c", "c"], dtype="category", index=idx)
    values = [1, 2, 2, 2, 3, 4, 5]
    df = pd.DataFrame({"cats": cats, "values": values}, index=idx)
    df.iloc[2:4, :]
    df.iloc[2:4, :].dtypes
    df.loc["h":"j", "cats"]
    df[df["cats"] == "b"]

An example where the category type is not preserved is if you take one single
row: the resulting ``Series`` is of dtype ``object``:

.. ipython:: python

    # get the complete "h" row as a Series
    df.loc["h", :]

Returning a single item from categorical data will also return the value, not a categorical
of length "1".

.. ipython:: python

    df.iat[0, 0]
    df["cats"] = df["cats"].cat.rename_categories(["x", "y", "z"])
    df.at["h", "cats"]  # returns a string

.. note::
    The is in contrast to R's ``factor`` function, where ``factor(c(1,2,3))[1]``
    returns a single value ``factor``.

To get a single value ``Series`` of type ``category``, you pass in a list with
a single value:

.. ipython:: python

    df.loc[["h"], "cats"]

String and datetime accessors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The accessors  ``.dt`` and ``.str`` will work if the ``s.cat.categories`` are of
an appropriate type:


.. ipython:: python

    str_s = pd.Series(list("aabb"))
    str_cat = str_s.astype("category")
    str_cat
    str_cat.str.contains("a")

    date_s = pd.Series(pd.date_range("1/1/2015", periods=5))
    date_cat = date_s.astype("category")
    date_cat
    date_cat.dt.day

.. note::

    The returned ``Series`` (or ``DataFrame``) is of the same type as if you used the
    ``.str.<method>`` / ``.dt.<method>`` on a ``Series`` of that type (and not of
    type ``category``!).

That means, that the returned values from methods and properties on the accessors of a
``Series`` and the returned values from methods and properties on the accessors of this
``Series`` transformed to one of type ``category`` will be equal:

.. ipython:: python

    ret_s = str_s.str.contains("a")
    ret_cat = str_cat.str.contains("a")
    ret_s.dtype == ret_cat.dtype
    ret_s == ret_cat

.. note::

    The work is done on the ``categories`` and then a new ``Series`` is constructed. This has
    some performance implication if you have a ``Series`` of type string, where lots of elements
    are repeated (i.e. the number of unique elements in the ``Series`` is a lot smaller than the
    length of the ``Series``). In this case it can be faster to convert the original ``Series``
    to one of type ``category`` and use ``.str.<method>`` or ``.dt.<property>`` on that.

Setting
~~~~~~~

Setting values in a categorical column (or ``Series``) works as long as the
value is included in the ``categories``:

.. ipython:: python

    idx = pd.Index(["h", "i", "j", "k", "l", "m", "n"])
    cats = pd.Categorical(["a", "a", "a", "a", "a", "a", "a"], categories=["a", "b"])
    values = [1, 1, 1, 1, 1, 1, 1]
    df = pd.DataFrame({"cats": cats, "values": values}, index=idx)

    df.iloc[2:4, :] = [["b", 2], ["b", 2]]
    df
    try:
        df.iloc[2:4, :] = [["c", 3], ["c", 3]]
    except TypeError as e:
        print("TypeError:", str(e))

Setting values by assigning categorical data will also check that the ``categories`` match:

.. ipython:: python

    df.loc["j":"k", "cats"] = pd.Categorical(["a", "a"], categories=["a", "b"])
    df
    try:
        df.loc["j":"k", "cats"] = pd.Categorical(["b", "b"], categories=["a", "b", "c"])
    except TypeError as e:
        print("TypeError:", str(e))

Assigning a ``Categorical`` to parts of a column of other types will use the values:

.. ipython:: python
    :okwarning:

    df = pd.DataFrame({"a": [1, 1, 1, 1, 1], "b": ["a", "a", "a", "a", "a"]})
    df.loc[1:2, "a"] = pd.Categorical([2, 2], categories=[2, 3])
    df.loc[2:3, "b"] = pd.Categorical(["b", "b"], categories=["a", "b"])
    df
    df.dtypes

.. _categorical.merge:
.. _categorical.concat:

Merging / concatenation
~~~~~~~~~~~~~~~~~~~~~~~

By default, combining ``Series`` or ``DataFrames`` which contain the same
categories results in ``category`` dtype, otherwise results will depend on the
dtype of the underlying categories. Merges that result in non-categorical
dtypes will likely have higher memory usage. Use ``.astype`` or
``union_categoricals`` to ensure ``category`` results.

.. ipython:: python

   from pandas.api.types import union_categoricals

   # same categories
   s1 = pd.Series(["a", "b"], dtype="category")
   s2 = pd.Series(["a", "b", "a"], dtype="category")
   pd.concat([s1, s2])

   # different categories
   s3 = pd.Series(["b", "c"], dtype="category")
   pd.concat([s1, s3])

   # Output dtype is inferred based on categories values
   int_cats = pd.Series([1, 2], dtype="category")
   float_cats = pd.Series([3.0, 4.0], dtype="category")
   pd.concat([int_cats, float_cats])

   pd.concat([s1, s3]).astype("category")
   union_categoricals([s1.array, s3.array])

The following table summarizes the results of merging ``Categoricals``:

+-------------------+------------------------+----------------------+-----------------------------+
| arg1              | arg2                   |      identical       | result                      |
+===================+========================+======================+=============================+
| category          | category               | True                 | category                    |
+-------------------+------------------------+----------------------+-----------------------------+
| category (object) | category (object)      | False                | object (dtype is inferred)  |
+-------------------+------------------------+----------------------+-----------------------------+
| category (int)    | category (float)       | False                | float (dtype is inferred)   |
+-------------------+------------------------+----------------------+-----------------------------+

.. _categorical.union:

Unioning
~~~~~~~~

If you want to combine categoricals that do not necessarily have the same
categories, the :func:`~pandas.api.types.union_categoricals` function will
combine a list-like of categoricals. The new categories will be the union of
the categories being combined.

.. ipython:: python

    from pandas.api.types import union_categoricals

    a = pd.Categorical(["b", "c"])
    b = pd.Categorical(["a", "b"])
    union_categoricals([a, b])

By default, the resulting categories will be ordered as
they appear in the data. If you want the categories to
be lexsorted, use ``sort_categories=True`` argument.

.. ipython:: python

    union_categoricals([a, b], sort_categories=True)

``union_categoricals`` also works with the "easy" case of combining two
categoricals of the same categories and order information
(e.g. what you could also ``append`` for).

.. ipython:: python

    a = pd.Categorical(["a", "b"], ordered=True)
    b = pd.Categorical(["a", "b", "a"], ordered=True)
    union_categoricals([a, b])

The below raises ``TypeError`` because the categories are ordered and not identical.

.. ipython:: python
   :okexcept:

   a = pd.Categorical(["a", "b"], ordered=True)
   b = pd.Categorical(["a", "b", "c"], ordered=True)
   union_categoricals([a, b])

Ordered categoricals with different categories or orderings can be combined by
using the ``ignore_ordered=True`` argument.

.. ipython:: python

    a = pd.Categorical(["a", "b", "c"], ordered=True)
    b = pd.Categorical(["c", "b", "a"], ordered=True)
    union_categoricals([a, b], ignore_order=True)

:func:`~pandas.api.types.union_categoricals` also works with a
``CategoricalIndex``, or ``Series`` containing categorical data, but note that
the resulting array will always be a plain ``Categorical``:

.. ipython:: python

    a = pd.Series(["b", "c"], dtype="category")
    b = pd.Series(["a", "b"], dtype="category")
    union_categoricals([a, b])

.. note::

   ``union_categoricals`` may recode the integer codes for categories
   when combining categoricals.  This is likely what you want,
   but if you are relying on the exact numbering of the categories, be
   aware.

   .. ipython:: python

      c1 = pd.Categorical(["b", "c"])
      c2 = pd.Categorical(["a", "b"])

      c1
      # "b" is coded to 0
      c1.codes

      c2
      # "b" is coded to 1
      c2.codes

      c = union_categoricals([c1, c2])
      c
      # "b" is coded to 0 throughout, same as c1, different from c2
      c.codes


Getting data in/out
-------------------

You can write data that contains ``category`` dtypes to a ``HDFStore``.
See :ref:`here <io.hdf5-categorical>` for an example and caveats.

It is also possible to write data to and reading data from *Stata* format files.
See :ref:`here <io.stata-categorical>` for an example and caveats.

Writing to a CSV file will convert the data, effectively removing any information about the
categorical (categories and ordering). So if you read back the CSV file you have to convert the
relevant columns back to ``category`` and assign the right categories and categories ordering.

.. ipython:: python

    import io

    s = pd.Series(pd.Categorical(["a", "b", "b", "a", "a", "d"]))
    # rename the categories
    s = s.cat.rename_categories(["very good", "good", "bad"])
    # reorder the categories and add missing categories
    s = s.cat.set_categories(["very bad", "bad", "medium", "good", "very good"])
    df = pd.DataFrame({"cats": s, "vals": [1, 2, 3, 4, 5, 6]})
    csv = io.StringIO()
    df.to_csv(csv)
    df2 = pd.read_csv(io.StringIO(csv.getvalue()))
    df2.dtypes
    df2["cats"]
    # Redo the category
    df2["cats"] = df2["cats"].astype("category")
    df2["cats"] = df2["cats"].cat.set_categories(
        ["very bad", "bad", "medium", "good", "very good"]
    )
    df2.dtypes
    df2["cats"]

The same holds for writing to a SQL database with ``to_sql``.

Missing data
------------

pandas primarily uses the value ``np.nan`` to represent missing data. It is by
default not included in computations. See the :ref:`Missing Data section
<missing_data>`.

Missing values should **not** be included in the Categorical's ``categories``,
only in the ``values``.
Instead, it is understood that NaN is different, and is always a possibility.
When working with the Categorical's ``codes``, missing values will always have
a code of ``-1``.

.. ipython:: python

    s = pd.Series(["a", "b", np.nan, "a"], dtype="category")
    # only two categories
    s
    s.cat.codes


Methods for working with missing data, e.g. :meth:`~Series.isna`, :meth:`~Series.fillna`,
:meth:`~Series.dropna`, all work normally:

.. ipython:: python

    s = pd.Series(["a", "b", np.nan], dtype="category")
    s
    pd.isna(s)
    s.fillna("a")

Differences to R's ``factor``
-----------------------------

The following differences to R's factor functions can be observed:

* R's ``levels`` are named ``categories``.
* R's ``levels`` are always of type string, while ``categories`` in pandas can be of any dtype.
* It's not possible to specify labels at creation time. Use ``s.cat.rename_categories(new_labels)``
  afterwards.
* In contrast to R's ``factor`` function, using categorical data as the sole input to create a
  new categorical series will *not* remove unused categories but create a new categorical series
  which is equal to the passed in one!
* R allows for missing values to be included in its ``levels`` (pandas' ``categories``). pandas
  does not allow ``NaN`` categories, but missing values can still be in the ``values``.


Gotchas
-------

.. _categorical.rfactor:

Memory usage
~~~~~~~~~~~~

.. _categorical.memory:

The memory usage of a ``Categorical`` is proportional to the number of categories plus the length of the data. In contrast,
an ``object`` dtype is a constant times the length of the data.

.. ipython:: python

   s = pd.Series(["foo", "bar"] * 1000)

   # object dtype
   s.nbytes

   # category dtype
   s.astype("category").nbytes

.. note::

   If the number of categories approaches the length of the data, the ``Categorical`` will use nearly the same or
   more memory than an equivalent ``object`` dtype representation.

   .. ipython:: python

      s = pd.Series(["foo%04d" % i for i in range(2000)])

      # object dtype
      s.nbytes

      # category dtype
      s.astype("category").nbytes


``Categorical`` is not a ``numpy`` array
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently, categorical data and the underlying ``Categorical`` is implemented as a Python
object and not as a low-level NumPy array dtype. This leads to some problems.

NumPy itself doesn't know about the new ``dtype``:

.. ipython:: python

    try:
        np.dtype("category")
    except TypeError as e:
        print("TypeError:", str(e))

    dtype = pd.Categorical(["a"]).dtype
    try:
        np.dtype(dtype)
    except TypeError as e:
        print("TypeError:", str(e))

Dtype comparisons work:

.. ipython:: python

    dtype == np.str_
    np.str_ == dtype

To check if a Series contains Categorical data, use ``hasattr(s, 'cat')``:

.. ipython:: python

    hasattr(pd.Series(["a"], dtype="category"), "cat")
    hasattr(pd.Series(["a"]), "cat")

Using NumPy functions on a ``Series`` of type ``category`` should not work as ``Categoricals``
are not numeric data (even in the case that ``.categories`` is numeric).

.. ipython:: python

    s = pd.Series(pd.Categorical([1, 2, 3, 4]))
    try:
        np.sum(s)
        # same with np.log(s),...
    except TypeError as e:
        print("TypeError:", str(e))

.. note::
    If such a function works, please file a bug at https://github.com/pandas-dev/pandas!

dtype in apply
~~~~~~~~~~~~~~

pandas currently does not preserve the dtype in apply functions: If you apply along rows you get
a ``Series`` of ``object`` ``dtype`` (same as getting a row -> getting one element will return a
basic type) and applying along columns will also convert to object. ``NaN`` values are unaffected.
You can use ``fillna`` to handle missing values before applying a function.

.. ipython:: python

    df = pd.DataFrame(
        {
            "a": [1, 2, 3, 4],
            "b": ["a", "b", "c", "d"],
            "cats": pd.Categorical([1, 2, 3, 2]),
        }
    )
    df.apply(lambda row: type(row["cats"]), axis=1)
    df.apply(lambda col: col.dtype, axis=0)

Categorical index
~~~~~~~~~~~~~~~~~

``CategoricalIndex`` is a type of index that is useful for supporting
indexing with duplicates. This is a container around a ``Categorical``
and allows efficient indexing and storage of an index with a large number of duplicated elements.
See the :ref:`advanced indexing docs <advanced.categoricalindex>` for a more detailed
explanation.

Setting the index will create a ``CategoricalIndex``:

.. ipython:: python

    cats = pd.Categorical([1, 2, 3, 4], categories=[4, 2, 3, 1])
    strings = ["a", "b", "c", "d"]
    values = [4, 2, 3, 1]
    df = pd.DataFrame({"strings": strings, "values": values}, index=cats)
    df.index
    # This now sorts by the categories order
    df.sort_index()

Side effects
~~~~~~~~~~~~

Constructing a ``Series`` from a ``Categorical`` will not copy the input
``Categorical``. This means that changes to the ``Series`` will in most cases
change the original ``Categorical``:

.. ipython:: python

    cat = pd.Categorical([1, 2, 3, 10], categories=[1, 2, 3, 4, 10])
    s = pd.Series(cat, name="cat")
    cat
    s.iloc[0:2] = 10
    cat

Use ``copy=True`` to prevent such a behaviour or simply don't reuse ``Categoricals``:

.. ipython:: python

    cat = pd.Categorical([1, 2, 3, 10], categories=[1, 2, 3, 4, 10])
    s = pd.Series(cat, name="cat", copy=True)
    cat
    s.iloc[0:2] = 10
    cat

.. note::

    This also happens in some cases when you supply a NumPy array instead of a ``Categorical``:
    using an int array (e.g. ``np.array([1,2,3,4])``) will exhibit the same behavior, while using
    a string array (e.g. ``np.array(["a","b","c","a"])``) will not.
