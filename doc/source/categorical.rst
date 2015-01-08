.. _categorical:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import random
   import os
   np.random.seed(123456)
   from pandas import options
   from pandas import *
   import pandas as pd
   np.set_printoptions(precision=4, suppress=True)
   options.display.mpl_style='default'
   options.display.max_rows=15


****************
Categorical Data
****************

.. versionadded:: 0.15

.. note::
    While there was in `pandas.Categorical` in earlier versions, the ability to use
    categorical data in `Series` and `DataFrame` is new.


This is a introduction to pandas categorical data type, including a short comparison
with R's ``factor``.

`Categoricals` are a pandas data type, which correspond to categorical variables in
statistics: a variable, which can take on only a limited, and usually fixed,
number of possible values (`categories`; `levels` in R). Examples are gender, social class,
blood types, country affiliations, observation time or ratings via Likert scales.

In contrast to statistical categorical variables, categorical data might have an order (e.g.
'strongly agree' vs 'agree' or 'first observation' vs. 'second observation'), but numerical
operations (additions, divisions, ...) are not possible.

All values of categorical data are either in `categories` or `np.nan`. Order is defined by
the order of `categories`, not lexical order of the values. Internally, the data structure
consists of a `categories` array and an integer array of `codes` which point to the real value in
the `categories` array.

The categorical data type is useful in the following cases:

* A string variable consisting of only a few different values. Converting such a string
  variable to a categorical variable will save some memory, see :ref:`here <categorical.memory>`.
* The lexical order of a variable is not the same as the logical order ("one", "two", "three").
  By converting to a categorical and specifying an order on the categories, sorting and
  min/max will use the logical order instead of the lexical order, see :ref:`here <categorical.sort>`.
* As a signal to other python libraries that this column should be treated as a categorical
  variable (e.g. to use suitable statistical methods or plot types).

See also the :ref:`API docs on categoricals<api.categorical>`.

Object Creation
---------------

Categorical `Series` or columns in a `DataFrame` can be created in several ways:

By specifying ``dtype="category"`` when constructing a `Series`:

.. ipython:: python

    s = Series(["a","b","c","a"], dtype="category")
    s

By converting an existing `Series` or column to a ``category`` dtype:

.. ipython:: python

    df = DataFrame({"A":["a","b","c","a"]})
    df["B"] = df["A"].astype('category')
    df

By using some special functions:

.. ipython:: python

    df = DataFrame({'value': np.random.randint(0, 100, 20)})
    labels = [ "{0} - {1}".format(i, i + 9) for i in range(0, 100, 10) ]

    df['group'] = pd.cut(df.value, range(0, 105, 10), right=False, labels=labels)
    df.head(10)

See :ref:`documentation <reshaping.tile.cut>` for :func:`~pandas.cut`.

By passing a :class:`pandas.Categorical` object to a `Series` or assigning it to a `DataFrame`.
This is the only possibility to specify differently ordered categories (or no order at all) at
creation time and the only reason to use :class:`pandas.Categorical` directly:

.. ipython:: python

    raw_cat = Categorical(["a","b","c","a"], categories=["b","c","d"],
                             ordered=False)
    s = Series(raw_cat)
    s
    df = DataFrame({"A":["a","b","c","a"]})
    df["B"] = raw_cat
    df

Categorical data has a specific ``category`` :ref:`dtype <basics.dtypes>`:

.. ipython:: python

    df.dtypes

.. note::

    In contrast to R's `factor` function, categorical data is not converting input values to
    strings and categories will end up the same data type as the original values.

.. note::

    In contrast to R's `factor` function, there is currently no way to assign/change labels at
    creation time. Use `categories` to change the categories after creation time.

To get back to the original Series or `numpy` array, use ``Series.astype(original_dtype)`` or
``np.asarray(categorical)``:

.. ipython:: python

    s = Series(["a","b","c","a"])
    s
    s2 = s.astype('category')
    s2
    s3 = s2.astype('string')
    s3
    np.asarray(s2)

If you have already `codes` and `categories`, you can use the :func:`~pandas.Categorical.from_codes`
constructor to save the factorize step during normal constructor mode:

.. ipython:: python

    splitter = np.random.choice([0,1], 5, p=[0.5,0.5])
    s = Series(Categorical.from_codes(splitter, categories=["train", "test"]))

Description
-----------

Using ``.describe()`` on categorical data will produce similar output to a `Series` or
`DataFrame` of type ``string``.

.. ipython:: python

    cat = Categorical(["a","c","c",np.nan], categories=["b","a","c",np.nan] )
    df = DataFrame({"cat":cat, "s":["a","c","c",np.nan]})
    df.describe()
    df["cat"].describe()

Working with categories
-----------------------

Categorical data has a `categories` and a `ordered` property, which list their possible values and
whether the ordering matters or not. These properties are exposed as ``s.cat.categories`` and
``s.cat.ordered``. If you don't manually specify categories and ordering, they are inferred from the
passed in values.

.. ipython:: python

    s = Series(["a","b","c","a"], dtype="category")
    s.cat.categories
    s.cat.ordered

It's also possible to pass in the categories in a specific order:

.. ipython:: python

    s = Series(Categorical(["a","b","c","a"], categories=["c","b","a"]))
    s.cat.categories
    s.cat.ordered

.. note::
    New categorical data is automatically ordered if the passed in values are sortable or a
    `categories` argument is supplied. This is a difference to R's `factors`, which are unordered
    unless explicitly told to be ordered (``ordered=TRUE``). You can of course overwrite that by
    passing in an explicit ``ordered=False``.


Renaming categories
~~~~~~~~~~~~~~~~~~~

Renaming categories is done by assigning new values to the ``Series.cat.categories`` property or
by using the :func:`Categorical.rename_categories` method:

.. ipython:: python

    s = Series(["a","b","c","a"], dtype="category")
    s
    s.cat.categories = ["Group %s" % g for g in s.cat.categories]
    s
    s.cat.rename_categories([1,2,3])

.. note::

    In contrast to R's `factor`, categorical data can have categories of other types than string.

.. note::

    Be aware that assigning new categories is an inplace operations, while most other operation
    under ``Series.cat`` per default return a new Series of dtype `category`.

Categories must be unique or a `ValueError` is raised:

.. ipython:: python

    try:
        s.cat.categories = [1,1,1]
    except ValueError as e:
        print("ValueError: " + str(e))

Appending new categories
~~~~~~~~~~~~~~~~~~~~~~~~

Appending categories can be done by using the :func:`Categorical.add_categories` method:

.. ipython:: python

    s = s.cat.add_categories([4])
    s.cat.categories
    s

Removing categories
~~~~~~~~~~~~~~~~~~~

Removing categories can be done by using the :func:`Categorical.remove_categories` method. Values
which are removed are replaced by ``np.nan``.:

.. ipython:: python

    s = s.cat.remove_categories([4])
    s

Removing unused categories
~~~~~~~~~~~~~~~~~~~~~~~~~~

Removing unused categories can also be done:

.. ipython:: python

    s = Series(Categorical(["a","b","a"], categories=["a","b","c","d"]))
    s
    s.cat.remove_unused_categories()

Setting categories
~~~~~~~~~~~~~~~~~~

If you want to do remove and add new categories in one step (which has some speed advantage),
or simply set the categories to a predefined scale, use :func:`Categorical.set_categories`.

.. ipython:: python

    s = Series(["one","two","four", "-"], dtype="category")
    s
    s = s.cat.set_categories(["one","two","three","four"])
    s

.. note::
    Be aware that :func:`Categorical.set_categories` cannot know whether some category is omitted
    intentionally or because it is misspelled or (under Python3) due to a type difference (e.g.,
    numpys S1 dtype and python strings). This can result in surprising behaviour!

Sorting and Order
-----------------

.. _categorical.sort:

If categorical data is ordered (``s.cat.ordered == True``), then the order of the categories has a
meaning and certain operations are possible. If the categorical is unordered, a `TypeError` is
raised.

.. ipython:: python

    s = Series(Categorical(["a","b","c","a"], ordered=False))
    try:
        s.sort()
    except TypeError as e:
        print("TypeError: " + str(e))
    s = Series(["a","b","c","a"], dtype="category") # ordered per default!
    s.sort()
    s
    s.min(), s.max()

Sorting will use the order defined by categories, not any lexical order present on the data type.
This is even true for strings and numeric data:

.. ipython:: python

    s = Series([1,2,3,1], dtype="category")
    s.cat.categories = [2,3,1]
    s
    s.sort()
    s
    s.min(), s.max()


Reordering
~~~~~~~~~~

Reordering the categories is possible via the :func:`Categorical.reorder_categories` and
the :func:`Categorical.set_categories` methods. For :func:`Categorical.reorder_categories`, all
old categories must be included in the new categories and no new categories are allowed. This will
necessarily make the sort order the same as the categories order.

.. ipython:: python

    s = Series([1,2,3,1], dtype="category")
    s = s.cat.reorder_categories([2,3,1])
    s
    s.sort()
    s
    s.min(), s.max()

.. note::

    Note the difference between assigning new categories and reordering the categories: the first
    renames categories and therefore the individual values in the `Series`, but if the first
    position was sorted last, the renamed value will still be sorted last. Reordering means that the
    way values are sorted is different afterwards, but not that individual values in the
    `Series` are changed.

.. note::

    If the `Categorical` is not ordered, ``Series.min()`` and ``Series.max()`` will raise
    `TypeError`. Numeric operations like ``+``, ``-``, ``*``, ``/`` and operations based on them
    (e.g.``Series.median()``, which would need to compute the mean between two values if the length
    of an array is even) do not work and raise a `TypeError`.

Multi Column Sorting
~~~~~~~~~~~~~~~~~~~~

A categorical dtyped column will partcipate in a multi-column sort in a similar manner to other columns.
The ordering of the categorical is determined by the ``categories`` of that columns.

.. ipython:: python

   dfs = DataFrame({'A' : Categorical(list('bbeebbaa'),categories=['e','a','b']),
                    'B' : [1,2,1,2,2,1,2,1] })
   dfs.sort(['A','B'])

Reordering the ``categories``, changes a future sort.

.. ipython:: python

   dfs['A'] = dfs['A'].cat.reorder_categories(['a','b','e'])
   dfs.sort(['A','B'])

Comparisons
-----------

Comparing categorical data with other objects is possible in three cases:

 * comparing equality (``==`` and ``!=``) to a list-like object (list, Series, array,
   ...) of the same length as the categorical data.
 * all comparisons (``==``, ``!=``, ``>``, ``>=``, ``<``, and ``<=``) of categorical data to
   another categorical Series, when ``ordered==True`` and the `categories` are the same.
 * all comparisons of a categorical data to a scalar.

All other comparisons, especially "non-equality" comparisons of two categoricals with different
categories or a categorical with any list-like object, will raise a TypeError.

.. note::

    Any "non-equality" comparisons of categorical data with a `Series`, `np.array`, `list` or
    categorical data with different categories or ordering will raise an `TypeError` because custom
    categories ordering could be interpreted in two ways: one with taking in account the
    ordering and one without.

.. ipython:: python

    cat = Series(Categorical([1,2,3], categories=[3,2,1]))
    cat_base = Series(Categorical([2,2,2], categories=[3,2,1]))
    cat_base2 = Series(Categorical([2,2,2]))

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
    cat == np.array([1,2,3])
    cat == 2

This doesn't work because the categories are not the same:

.. ipython:: python

    try:
        cat > cat_base2
    except TypeError as e:
         print("TypeError: " + str(e))

If you want to do a "non-equality" comparison of a categorical series with a list-like object
which is not categorical data, you need to be explicit and convert the categorical data back to
the original values:

.. ipython:: python

    base = np.array([1,2,3])

    try:
        cat > base
    except TypeError as e:
         print("TypeError: " + str(e))

    np.asarray(cat) > base

Operations
----------

Apart from ``Series.min()``, ``Series.max()`` and ``Series.mode()``, the following operations are
possible with categorical data:

`Series` methods like `Series.value_counts()` will use all categories, even if some categories are not
present in the data:

.. ipython:: python

    s = Series(Categorical(["a","b","c","c"], categories=["c","a","b","d"]))
    s.value_counts()

Groupby will also show "unused" categories:

.. ipython:: python

    cats = Categorical(["a","b","b","b","c","c","c"], categories=["a","b","c","d"])
    df = DataFrame({"cats":cats,"values":[1,2,2,2,3,4,5]})
    df.groupby("cats").mean()

    cats2 = Categorical(["a","a","b","b"], categories=["a","b","c"])
    df2 = DataFrame({"cats":cats2,"B":["c","d","c","d"], "values":[1,2,3,4]})
    df2.groupby(["cats","B"]).mean()


Pivot tables:

.. ipython:: python

    raw_cat = Categorical(["a","a","b","b"], categories=["a","b","c"])
    df = DataFrame({"A":raw_cat,"B":["c","d","c","d"], "values":[1,2,3,4]})
    pd.pivot_table(df, values='values', index=['A', 'B'])

Data munging
------------

The optimized pandas data access methods  ``.loc``, ``.iloc``, ``.ix`` ``.at``, and ``.iat``,
work as normal, the only difference is the return type (for getting) and
that only values already in `categories` can be assigned.

Getting
~~~~~~~

If the slicing operation returns either a `DataFrame` or a column of type `Series`,
the ``category`` dtype is preserved.

.. ipython:: python

    idx = Index(["h","i","j","k","l","m","n",])
    cats = Series(["a","b","b","b","c","c","c"], dtype="category", index=idx)
    values= [1,2,2,2,3,4,5]
    df = DataFrame({"cats":cats,"values":values}, index=idx)
    df.iloc[2:4,:]
    df.iloc[2:4,:].dtypes
    df.loc["h":"j","cats"]
    df.ix["h":"j",0:1]
    df[df["cats"] == "b"]

An example where the category type is not preserved is if you take one single row: the
resulting `Series` is of dtype ``object``:

.. ipython:: python

    # get the complete "h" row as a Series
    df.loc["h", :]

Returning a single item from categorical data will also return the value, not a categorical
of length "1".

.. ipython:: python

    df.iat[0,0]
    df["cats"].cat.categories = ["x","y","z"]
    df.at["h","cats"] # returns a string

.. note::
    This is a difference to R's `factor` function, where ``factor(c(1,2,3))[1]``
    returns a single value `factor`.

To get a single value `Series` of type ``category`` pass in a list with a single value:

.. ipython:: python

    df.loc[["h"],"cats"]

Setting
~~~~~~~

Setting values in a categorical column (or `Series`) works as long as the value is included in the
`categories`:

.. ipython:: python

    idx = Index(["h","i","j","k","l","m","n"])
    cats = Categorical(["a","a","a","a","a","a","a"], categories=["a","b"])
    values = [1,1,1,1,1,1,1]
    df = DataFrame({"cats":cats,"values":values}, index=idx)

    df.iloc[2:4,:] = [["b",2],["b",2]]
    df
    try:
        df.iloc[2:4,:] = [["c",3],["c",3]]
    except ValueError as e:
        print("ValueError: " + str(e))

Setting values by assigning categorical data will also check that the `categories` match:

.. ipython:: python

    df.loc["j":"k","cats"] = Categorical(["a","a"], categories=["a","b"])
    df
    try:
        df.loc["j":"k","cats"] = Categorical(["b","b"], categories=["a","b","c"])
    except ValueError as e:
        print("ValueError: " + str(e))

Assigning a `Categorical` to parts of a column of other types will use the values:

.. ipython:: python

    df = DataFrame({"a":[1,1,1,1,1], "b":["a","a","a","a","a"]})
    df.loc[1:2,"a"] = Categorical(["b","b"], categories=["a","b"])
    df.loc[2:3,"b"] = Categorical(["b","b"], categories=["a","b"])
    df
    df.dtypes


Merging
~~~~~~~

You can concat two `DataFrames` containing categorical data together,
but the categories of these categoricals need to be the same:

.. ipython:: python

    cat = Series(["a","b"], dtype="category")
    vals = [1,2]
    df = DataFrame({"cats":cat, "vals":vals})
    res = pd.concat([df,df])
    res
    res.dtypes

In this case the categories are not the same and so an error is raised:

.. ipython:: python

    df_different = df.copy()
    df_different["cats"].cat.categories = ["c","d"]
    try:
        pd.concat([df,df_different])
    except ValueError as e:
        print("ValueError: " + str(e))

The same applies to ``df.append(df_different)``.

Getting Data In/Out
-------------------

.. versionadded:: 0.15.2

Writing data (`Series`, `Frames`) to a HDF store that contains a ``category`` dtype was implemented
in 0.15.2. See :ref:`here <io.hdf5-categorical>` for an example and caveats.

Writing data to and reading data from *Stata* format files was implemented in
0.15.2. See :ref:`here <io.stata-categorical>` for an example and caveats.

Writing to a CSV file will convert the data, effectively removing any information about the
categorical (categories and ordering). So if you read back the CSV file you have to convert the
relevant columns back to `category` and assign the right categories and categories ordering.

.. ipython:: python
    :suppress:

    from pandas.compat import StringIO

.. ipython:: python

    s = Series(Categorical(['a', 'b', 'b', 'a', 'a', 'd']))
    # rename the categories
    s.cat.categories = ["very good", "good", "bad"]
    # reorder the categories and add missing categories
    s = s.cat.set_categories(["very bad", "bad", "medium", "good", "very good"])
    df = DataFrame({"cats":s, "vals":[1,2,3,4,5,6]})
    csv = StringIO()
    df.to_csv(csv)
    df2 = pd.read_csv(StringIO(csv.getvalue()))
    df2.dtypes
    df2["cats"]
    # Redo the category
    df2["cats"] = df2["cats"].astype("category")
    df2["cats"].cat.set_categories(["very bad", "bad", "medium", "good", "very good"],
                                   inplace=True)
    df2.dtypes
    df2["cats"]

The same holds for writing to a SQL database with ``to_sql``.

Missing Data
------------

pandas primarily uses the value `np.nan` to represent missing data. It is by
default not included in computations. See the :ref:`Missing Data section
<missing_data>`

There are two ways a `np.nan` can be represented in categorical data: either the value is not
available ("missing value") or `np.nan` is a valid category.

.. ipython:: python

    s = Series(["a","b",np.nan,"a"], dtype="category")
    # only two categories
    s
    s2 = Series(["a","b","c","a"], dtype="category")
    s2.cat.categories = [1,2,np.nan]
    # three categories, np.nan included
    s2

.. note::
    As integer `Series` can't include NaN, the categories were converted to `object`.

.. note::
    Missing value methods like ``isnull`` and ``fillna`` will take both missing values as well as
    `np.nan` categories into account:

.. ipython:: python

    c = Series(["a","b",np.nan], dtype="category")
    c.cat.set_categories(["a","b",np.nan], inplace=True)
    # will be inserted as a NA category:
    c[0] = np.nan
    s = Series(c)
    s
    pd.isnull(s)
    s.fillna("a")

Differences to R's `factor`
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following differences to R's factor functions can be observed:

* R's `levels` are named `categories`
* R's `levels` are always of type string, while `categories` in pandas can be of any dtype.
* New categorical data is automatically ordered if the passed in values are sortable or a
  `categories` argument is supplied. This is a difference to R's `factors`, which are unordered
  unless explicitly told to be ordered (``ordered=TRUE``).
* It's not possible to specify labels at creation time. Use ``s.cat.rename_categories(new_labels)``
  afterwards.
* In contrast to R's `factor` function, using categorical data as the sole input to create a
  new categorical series will *not* remove unused categories but create a new categorical series
  which is equal to the passed in one!

Gotchas
-------

.. _categorical.rfactor:

Memory Usage
~~~~~~~~~~~~

.. _categorical.memory:

The memory usage of a ``Categorical`` is proportional to the number of categories times the length of the data. In contrast,
an ``object`` dtype is a constant times the length of the data.

.. ipython:: python

   s = Series(['foo','bar']*1000)

   # object dtype
   s.nbytes

   # category dtype
   s.astype('category').nbytes

.. note::

   If the number of categories approaches the length of the data, the ``Categorical`` will use nearly (or more) memory than an
   equivalent ``object`` dtype representation.

   .. ipython:: python

      s = Series(['foo%04d' % i for i in range(2000)])

      # object dtype
      s.nbytes

      # category dtype
      s.astype('category').nbytes


Old style constructor usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In earlier versions than pandas 0.15, a `Categorical` could be constructed by passing in precomputed
`codes` (called then `labels`) instead of values with categories. The `codes` were interpreted as
pointers to the categories with `-1` as `NaN`. This type of constructor useage is replaced by
the special constructor :func:`Categorical.from_codes`.

Unfortunately, in some special cases, using code which assumes the old style constructor usage
will work with the current pandas version, resulting in subtle bugs:

.. code-block:: python

    >>> cat = Categorical([1,2], [1,2,3])
    >>> # old version
    >>> cat.get_values()
    array([2, 3], dtype=int64)
    >>> # new version
    >>> cat.get_values()
    array([1, 2], dtype=int64)

.. warning::
    If you used `Categoricals` with older versions of pandas, please audit your code before
    upgrading and change your code to use the :func:`~pandas.Categorical.from_codes`
    constructor.

`Categorical` is not a `numpy` array
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently, categorical data and the underlying `Categorical` is implemented as a python
object and not as a low-level `numpy` array dtype. This leads to some problems.

`numpy` itself doesn't know about the new `dtype`:

.. ipython:: python

    try:
        np.dtype("category")
    except TypeError as e:
        print("TypeError: " + str(e))

    dtype = Categorical(["a"]).dtype
    try:
        np.dtype(dtype)
    except TypeError as e:
         print("TypeError: " + str(e))

Dtype comparisons work:

.. ipython:: python

    dtype == np.str_
    np.str_ == dtype

Using `numpy` functions on a `Series` of type ``category`` should not work as `Categoricals`
are not numeric data (even in the case that ``.categories`` is numeric).

.. ipython:: python

    s = Series(Categorical([1,2,3,4]))
    try:
        np.sum(s)
        #same with np.log(s),..
    except TypeError as e:
         print("TypeError: " + str(e))

.. note::
    If such a function works, please file a bug at https://github.com/pydata/pandas!

dtype in apply
~~~~~~~~~~~~~~

Pandas currently does not preserve the dtype in apply functions: If you apply along rows you get
a `Series` of ``object`` `dtype` (same as getting a row -> getting one element will return a
basic type) and applying along columns will also convert to object.

.. ipython:: python

    df = DataFrame({"a":[1,2,3,4],
                    "b":["a","b","c","d"],
                    "cats":Categorical([1,2,3,2])})
    df.apply(lambda row: type(row["cats"]), axis=1)
    df.apply(lambda col: col.dtype, axis=0)

No Categorical Index
~~~~~~~~~~~~~~~~~~~~

There is currently no index of type ``category``, so setting the index to categorical column will
convert the categorical data to a "normal" dtype first and therefore remove any custom
ordering of the categories:

.. ipython:: python

    cats = Categorical([1,2,3,4], categories=[4,2,3,1])
    strings = ["a","b","c","d"]
    values = [4,2,3,1]
    df = DataFrame({"strings":strings, "values":values}, index=cats)
    df.index
    # This should sort by categories but does not as there is no CategoricalIndex!
    df.sort_index()

.. note::
    This could change if a `CategoricalIndex` is implemented (see
    https://github.com/pydata/pandas/issues/7629)


Side Effects
~~~~~~~~~~~~

Constructing a `Series` from a `Categorical` will not copy the input `Categorical`. This
means that changes to the `Series` will in most cases change the original `Categorical`:

.. ipython:: python

    cat = Categorical([1,2,3,10], categories=[1,2,3,4,10])
    s = Series(cat, name="cat")
    cat
    s.iloc[0:2] = 10
    cat
    df = DataFrame(s)
    df["cat"].cat.categories = [1,2,3,4,5]
    cat

Use ``copy=True`` to prevent such a behaviour or simply don't reuse `Categoricals`:

.. ipython:: python

    cat = Categorical([1,2,3,10], categories=[1,2,3,4,10])
    s = Series(cat, name="cat", copy=True)
    cat
    s.iloc[0:2] = 10
    cat

.. note::
    This also happens in some cases when you supply a `numpy` array instead of a `Categorical`:
    using an int array (e.g. ``np.array([1,2,3,4])``) will exhibit the same behaviour, while using
    a string array (e.g. ``np.array(["a","b","c","a"])``) will not.
