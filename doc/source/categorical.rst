.. _categorical:

.. currentmodule:: pandas

.. ipython:: python
   :suppress:

   import numpy as np
   import random
   import os
   np.random.seed(123456)
   from pandas import options
   import pandas as pd
   np.set_printoptions(precision=4, suppress=True)
   options.display.mpl_style='default'
   options.display.max_rows=15


***********
Categorical
***********

.. versionadded:: 0.15

.. note::
    While there was in `pandas.Categorical` in earlier versions, the ability to use
    `Categorical` data in `Series` and `DataFrame` is new.


This is a introduction to pandas :class:`pandas.Categorical` type, including a short comparison
with R's `factor`.

`Categoricals` are a pandas data type, which correspond to categorical variables in
statistics: a variable, which can take on only a limited, and usually fixed,
number of possible values (commonly called `levels`). Examples are gender, social class,
blood types, country affiliations, observation time or ratings via Likert scales.

In contrast to statistical categorical variables, a `Categorical` might have an order (e.g.
'strongly agree' vs 'agree' or 'first observation' vs. 'second observation'), but numerical
operations (additions, divisions, ...) are not possible.

All values of the `Categorical` are either in `levels` or `np.nan`. Order is defined by
the order of the `levels`, not lexical order of the values. Internally, the data structure
consists of a levels array and an integer array of `codes` which point to the real value in the
levels array.

`Categoricals` are useful in the following cases:

* A string variable consisting of only a few different values. Converting such a string
  variable to a categorical variable will save some memory.
* The lexical order of a variable is not the same as the logical order ("one", "two", "three").
  By converting to a categorical and specifying an order on the levels, sorting and
  min/max will use the logical order instead of the lexical order.
* As a signal to other python libraries that this column should be treated as a categorical
  variable (e.g. to use suitable statistical methods or plot types)

See also the :ref:`API docs on Categoricals<api.categorical>`.

Object Creation
---------------

Categorical `Series` or columns in a `DataFrame` can be crated in several ways:

By passing a `Categorical` object to a `Series` or assigning it to a `DataFrame`:

.. ipython:: python

    raw_cat = pd.Categorical(["a","b","c","a"])
    s = pd.Series(raw_cat)
    s
    df = pd.DataFrame({"A":["a","b","c","a"]})
    df["B"] = raw_cat
    df

By converting an existing `Series` or column to a ``category`` type:

.. ipython:: python

    df = pd.DataFrame({"A":["a","b","c","a"]})
    df["B"] = df["A"].astype('category')
    df

By using some special functions:

.. ipython:: python

    df = pd.DataFrame({'value': np.random.randint(0, 100, 20)})
    labels = [ "{0} - {1}".format(i, i + 9) for i in range(0, 100, 10) ]

    df['group'] = pd.cut(df.value, range(0, 105, 10), right=False, labels=labels)
    df.head(10)

See :ref:`documentation <reshaping.tile.cut>` for :func:`~pandas.cut`.

`Categoricals` have a specific ``category`` :ref:`dtype <basics.dtypes>`:

.. ipython:: python

    df.dtypes

.. note::

    In contrast to R's `factor` function, a `Categorical` is not converting input values to
    string and levels will end up the same data type as the original values.

.. note::

    In contrast to R's `factor` function, there is currently no way to assign/change labels at
    creation time. Use `levels` to change the levels after creation time.

To get back to the original Series or `numpy` array, use ``Series.astype(original_dtype)`` or
``np.asarray(categorical)``:

.. ipython:: python

    s = pd.Series(["a","b","c","a"])
    s
    s2 = s.astype('category')
    s2
    s3 = s2.astype('string')
    s3
    np.asarray(s2.cat)

If you have already `codes` and `levels`, you can use the :func:`~pandas.Categorical.from_codes`
constructor to save the factorize step during normal constructor mode:

.. ipython:: python

    splitter = np.random.choice([0,1], 5, p=[0.5,0.5])
    pd.Categorical.from_codes(splitter, levels=["train", "test"])

Description
-----------

Using ``.describe()`` on a ``Categorical(...)`` or a ``Series(Categorical(...))`` will show
different output.


As part of a `Dataframe` or as a `Series` a similar output as for a `Series` of type ``string`` is
shown. Calling ``Categorical.describe()`` will show the frequencies for each level, with NA for
unused levels.

.. ipython:: python

    cat = pd.Categorical(["a","c","c",np.nan], levels=["b","a","c",np.nan] )
    df = pd.DataFrame({"cat":cat, "s":["a","c","c",np.nan]})
    df.describe()
    cat.describe()

Working with levels
-------------------

`Categoricals` have a `levels` property, which list their possible values. If you don't
manually specify levels, they are inferred from the passed in values. `Series` of type
``category`` expose the same interface via their `cat` property.

.. ipython:: python

    raw_cat = pd.Categorical(["a","b","c","a"])
    raw_cat.levels
    raw_cat.ordered
    # Series of type "category" also expose these interface via the .cat property:
    s = pd.Series(raw_cat)
    s.cat.levels
    s.cat.ordered

.. note::
    New `Categorical` are automatically ordered if the passed in values are sortable or a
    `levels` argument is supplied. This is a difference to R's `factors`, which are unordered
    unless explicitly told to be ordered (``ordered=TRUE``).

It's also possible to pass in the levels in a specific order:

.. ipython:: python

    raw_cat = pd.Categorical(["a","b","c","a"], levels=["c","b","a"])
    s = pd.Series(raw_cat)
    s.cat.levels
    s.cat.ordered

.. note::

    Passing in a `levels` argument implies ``ordered=True``. You can of course overwrite that by
    passing in an explicit ``ordered=False``.

Any value omitted in the levels argument will be replaced by `np.nan`:

.. ipython:: python

    raw_cat = pd.Categorical(["a","b","c","a"], levels=["a","b"])
    s = pd.Series(raw_cat)
    s.cat.levels
    s

Renaming levels is done by assigning new values to the ``Category.levels`` or
``Series.cat.levels`` property:

.. ipython:: python

    s = pd.Series(pd.Categorical(["a","b","c","a"]))
    s
    s.cat.levels = ["Group %s" % g for g in s.cat.levels]
    s
    s.cat.levels = [1,2,3]
    s

.. note::

    I contrast to R's `factor`, a `Categorical` can have levels of other types than string.

Levels must be unique or a `ValueError` is raised:

.. ipython:: python

    try:
        s.cat.levels = [1,1,1]
    except ValueError as e:
        print("ValueError: " + str(e))

Appending levels can be done by assigning a levels list longer than the current levels:

.. ipython:: python

    s.cat.levels = [1,2,3,4]
    s.cat.levels
    s

.. note::
    Adding levels in other positions can be done with ``.reorder_levels(<levels_including_new>)``.

Removing a level is also possible, but only the last level(s) can be removed by assigning a
shorter list than current levels. Values which are omitted are replaced by ``np.nan``.

.. ipython:: python

    s.cat.levels = [1,2]
    s

.. note::

    It's only possible to remove or add a level at the last position. If that's not where you want
    to remove an old or add a new level, use ``Category.reorder_levels(new_order)`` or
    ``Series.cat.reorder_levels(new_order)`` methods before or after.

Removing unused levels can also be done:

.. ipython:: python

    raw = pd.Categorical(["a","b","a"], levels=["a","b","c","d"])
    c = pd.Series(raw)
    raw
    raw.remove_unused_levels()
    raw
    c.cat.remove_unused_levels()
    c

.. note::

    In contrast to R's `factor` function, passing a `Categorical` as the sole input to the
    `Categorical` constructor will *not* remove unused levels but create a new `Categorical`
    which is equal to the passed in one!


Ordered or not...
-----------------

If a `Categoricals` is ordered (``cat.ordered == True``), then the order of the levels has a
meaning and certain operations are possible. If the categorical is unordered, a `TypeError` is
raised.

.. ipython:: python

    s = pd.Series(pd.Categorical(["a","b","c","a"], ordered=False))
    try:
        s.sort()
    except TypeError as e:
        print("TypeError: " + str(e))
    s = pd.Series(pd.Categorical(["a","b","c","a"], ordered=True))
    s.sort()
    s
    print(s.min(), s.max())

.. note::
    ``ordered=True`` is not necessary needed in the second case, as lists of strings are sortable
    and so the resulting `Categorical` is ordered.

Sorting will use the order defined by levels, not any lexical order present on the data type.
This is even true for strings and numeric data:

.. ipython:: python

    s = pd.Series(pd.Categorical([1,2,3,1]))
    s.cat.levels = [2,3,1]
    s
    s.sort()
    s
    print(s.min(), s.max())

Reordering the levels is possible via the ``Categorical.reorder_levels(new_levels)``  or
``Series.cat.reorder_levels(new_levels)`` methods. All old levels must be included in the new
levels.

.. ipython:: python

    s2 = pd.Series(pd.Categorical([1,2,3,1]))
    s2.cat.reorder_levels([2,3,1])
    s2
    s2.sort()
    s2
    print(s2.min(), s2.max())


.. note::
    Note the difference between assigning new level names and reordering the levels: the first
    renames levels and therefore the individual values in the `Series`, but if the first
    position was sorted last, the renamed value will still be sorted last. Reordering means that the
    way values are sorted is different afterwards, but not that individual values in the
    `Series` are changed.

You can also add new levels with :func:`Categorical.reorder_levels`, as long as you include all
old levels:

.. ipython:: python

    s3 = pd.Series(pd.Categorical(["a","b","d"]))
    s3.cat.reorder_levels(["a","b","c","d"])
    s3


Operations
----------

The following operations are possible with categorical data:

Comparing `Categoricals` with other objects is possible in two cases:

 * comparing a `Categorical` to another `Categorical`, when `level` and `ordered` is the same or
 * comparing a `Categorical` to a scalar.

All other comparisons will raise a TypeError.

.. ipython:: python

    cat = pd.Series(pd.Categorical([1,2,3], levels=[3,2,1]))
    cat_base = pd.Series(pd.Categorical([2,2,2], levels=[3,2,1]))
    cat_base2 = pd.Series(pd.Categorical([2,2,2]))

    cat
    cat_base
    cat_base2

Comparing to a categorical with the same levels and ordering or to a scalar works:

.. ipython:: python

    cat > cat_base
    cat > 2

This doesn't work because the levels are not the same:

.. ipython:: python

    try:
        cat > cat_base2
    except TypeError as e:
         print("TypeError: " + str(e))

.. note::

    Comparisons with `Series`, `np.array` or a `Categorical` with different levels or ordering
    will raise an `TypeError` because custom level ordering would result in two valid results:
    one with taking in account the ordering and one without. If you want to compare a `Categorical`
    with such a type, you need to be explicit and convert the `Categorical` to values:

.. ipython:: python

    base = np.array([1,2,3])

    try:
        cat > base
    except TypeError as e:
         print("TypeError: " + str(e))

    np.asarray(cat) > base

Getting the minimum and maximum, if the categorical is ordered:

.. ipython:: python

    s = pd.Series(pd.Categorical(["a","b","c","a"], levels=["c","a","b","d"]))
    print(s.min(), s.max())

.. note::

    If the `Categorical` is not ordered, ``Categorical.min()`` and ``Categorical.max()`` and the
    corresponding operations on `Series` will raise `TypeError`.

The mode:

.. ipython:: python

    raw_cat = pd.Categorical(["a","b","c","c"], levels=["c","a","b","d"])
    s = pd.Series(raw_cat)
    raw_cat.mode()
    s.mode()

.. note::

    Numeric operations like ``+``, ``-``, ``*``, ``/`` and operations based on them (e.g.
    ``.median()``, which would need to compute the mean between two values if the length of an
    array is even) do not work and raise a `TypeError`.

`Series` methods like `Series.value_counts()` will use all levels, even if some levels are not
present in the data:

.. ipython:: python

    s = pd.Series(pd.Categorical(["a","b","c","c"], levels=["c","a","b","d"]))
    s.value_counts()

Groupby will also show "unused" levels:

.. ipython:: python

    cats = pd.Categorical(["a","b","b","b","c","c","c"], levels=["a","b","c","d"])
    df = pd.DataFrame({"cats":cats,"values":[1,2,2,2,3,4,5]})
    df.groupby("cats").mean()

    cats2 = pd.Categorical(["a","a","b","b"], levels=["a","b","c"])
    df2 = pd.DataFrame({"cats":cats2,"B":["c","d","c","d"], "values":[1,2,3,4]})
    df2.groupby(["cats","B"]).mean()


Pivot tables:

.. ipython:: python

    raw_cat = pd.Categorical(["a","a","b","b"], levels=["a","b","c"])
    df = pd.DataFrame({"A":raw_cat,"B":["c","d","c","d"], "values":[1,2,3,4]})
    pd.pivot_table(df, values='values', index=['A', 'B'])

Data munging
------------

The optimized pandas data access methods  ``.loc``, ``.iloc``, ``.ix`` ``.at``, and ``.iat``,
work as normal, the only difference is the return type (for getting) and
that only values already in the levels can be assigned.

Getting
~~~~~~~

If the slicing operation returns either a `DataFrame` or a column of type `Series`,
the ``category`` dtype is preserved.

.. ipython:: python

    cats = pd.Categorical(["a","b","b","b","c","c","c"], levels=["a","b","c"])
    idx = pd.Index(["h","i","j","k","l","m","n",])
    values= [1,2,2,2,3,4,5]
    df = pd.DataFrame({"cats":cats,"values":values}, index=idx)
    df.iloc[2:4,:]
    df.iloc[2:4,:].dtypes
    df.loc["h":"j","cats"]
    df.ix["h":"j",0:1]
    df[df["cats"] == "b"]

An example where the `Categorical` is not preserved is if you take one single row: the
resulting `Series` is of dtype ``object``:

.. ipython:: python

    # get the complete "h" row as a Series
    df.loc["h", :]

Returning a single item from a `Categorical` will also return the value, not a `Categorical`
of length "1".

.. ipython:: python

    df.iat[0,0]
    df["cats"].cat.levels = ["x","y","z"]
    df.at["h","cats"] # returns a string

.. note::
    This is a difference to R's `factor` function, where ``factor(c(1,2,3))[1]``
    returns a single value `factor`.

To get a single value `Series` of type ``category`` pass in a single value list:

.. ipython:: python

    df.loc[["h"],"cats"]

Setting
~~~~~~~

Setting values in a categorical column (or `Series`) works as long as the value is included in the
`levels`:

.. ipython:: python

    cats = pd.Categorical(["a","a","a","a","a","a","a"], levels=["a","b"])
    idx = pd.Index(["h","i","j","k","l","m","n"])
    values = [1,1,1,1,1,1,1]
    df = pd.DataFrame({"cats":cats,"values":values}, index=idx)

    df.iloc[2:4,:] = [["b",2],["b",2]]
    df
    try:
        df.iloc[2:4,:] = [["c",3],["c",3]]
    except ValueError as e:
        print("ValueError: " + str(e))

Setting values by assigning a `Categorical` will also check that the `levels` match:

.. ipython:: python

    df.loc["j":"k","cats"] = pd.Categorical(["a","a"], levels=["a","b"])
    df
    try:
        df.loc["j":"k","cats"] = pd.Categorical(["b","b"], levels=["a","b","c"])
    except ValueError as e:
        print("ValueError: " + str(e))

Assigning a `Categorical` to parts of a column of other types will use the values:

.. ipython:: python

    df = pd.DataFrame({"a":[1,1,1,1,1], "b":["a","a","a","a","a"]})
    df.loc[1:2,"a"] = pd.Categorical(["b","b"], levels=["a","b"])
    df.loc[2:3,"b"] = pd.Categorical(["b","b"], levels=["a","b"])
    df
    df.dtypes


Merging
~~~~~~~

You can concat two `DataFrames` containing categorical data together,
but the levels of these `Categoricals` need to be the same:

.. ipython:: python

    cat = pd.Categorical(["a","b"], levels=["a","b"])
    vals = [1,2]
    df = pd.DataFrame({"cats":cat, "vals":vals})
    res = pd.concat([df,df])
    res
    res.dtypes

In this case the levels are not the same and so an error is raised:

.. ipython:: python

    df_different = df.copy()
    df_different["cats"].cat.levels = ["a","b","c"]
    try:
        pd.concat([df,df_different])
    except ValueError as e:
        print("ValueError: " + str(e))

The same applies to ``df.append(df)``.

Getting Data In/Out
-------------------

Writing data (`Series`, `Frames`) to a HDF store that contains a ``category`` dtype will currently
raise ``NotImplementedError``.

Writing to a CSV file will convert the data, effectively removing any information about the
`Categorical` (levels and ordering). So if you read back the CSV file you have to convert the
relevant columns back to `category` and assign the right levels and level ordering.

.. ipython:: python
    :suppress:

    from pandas.compat import StringIO

.. ipython:: python

    s = pd.Series(pd.Categorical(['a', 'b', 'b', 'a', 'a', 'd']))
    # rename the levels
    s.cat.levels = ["very good", "good", "bad"]
    # reorder the levels and add missing levels
    s.cat.reorder_levels(["very bad", "bad", "medium", "good", "very good"])
    df = pd.DataFrame({"cats":s, "vals":[1,2,3,4,5,6]})
    csv = StringIO()
    df.to_csv(csv)
    df2 = pd.read_csv(StringIO(csv.getvalue()))
    df2.dtypes
    df2["cats"]
    # Redo the category
    df2["cats"] = df2["cats"].astype("category")
    df2["cats"].cat.reorder_levels(["very bad", "bad", "medium", "good", "very good"])
    df2.dtypes
    df2["cats"]


Missing Data
------------

pandas primarily uses the value `np.nan` to represent missing data. It is by
default not included in computations. See the :ref:`Missing Data section
<missing_data>`

There are two ways a `np.nan` can be represented in `Categorical`: either the value is not
available ("missing value") or `np.nan` is a valid level.

.. ipython:: python

    s = pd.Series(pd.Categorical(["a","b",np.nan,"a"]))
    s
    # only two levels
    s.cat.levels
    s2 = pd.Series(pd.Categorical(["a","b","c","a"]))
    s2.cat.levels = [1,2,np.nan]
    s2
    # three levels, np.nan included
    # Note: as int arrays can't hold NaN the levels were converted to object
    s2.cat.levels

.. note::
    Missing value methods like ``isnull`` and ``fillna`` will take both missing values as well as
    `np.nan` levels into account:

.. ipython:: python

    c = pd.Categorical(["a","b",np.nan])
    c.levels = ["a","b",np.nan]
    # will be inserted as a NA level:
    c[0] = np.nan
    s = pd.Series(c)
    s
    pd.isnull(s)
    s.fillna("a")


Gotchas
-------

`Categorical` is not a `numpy` array
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently, `Categorical` and the corresponding ``category`` `Series` is implemented as a python
object and not as a low level `numpy` array dtype. This leads to some problems.

`numpy` itself doesn't know about the new `dtype`:

.. ipython:: python

    try:
        np.dtype("category")
    except TypeError as e:
        print("TypeError: " + str(e))

    dtype = pd.Categorical(["a"]).dtype
    try:
        np.dtype(dtype)
    except TypeError as e:
         print("TypeError: " + str(e))

Dtype comparisons work:

.. ipython:: python

    dtype == np.str_
    np.str_ == dtype

Using `numpy` functions on a `Series` of type ``category`` should not work as `Categoricals`
are not numeric data (even in the case that ``.levels`` is numeric).

.. ipython:: python

    s = pd.Series(pd.Categorical([1,2,3,4]))
    try:
        np.sum(s)
        #same with np.log(s),..
    except TypeError as e:
         print("TypeError: " + str(e))

.. note::
    If such a function works, please file a bug at https://github.com/pydata/pandas!


Side effects
~~~~~~~~~~~~

Constructing a `Series` from a `Categorical` will not copy the input `Categorical`. This
means that changes to the `Series` will in most cases change the original `Categorical`:

.. ipython:: python

    cat = pd.Categorical([1,2,3,10], levels=[1,2,3,4,10])
    s = pd.Series(cat, name="cat")
    cat
    s.iloc[0:2] = 10
    cat
    df = pd.DataFrame(s)
    df["cat"].cat.levels = [1,2,3,4,5]
    cat

Use ``copy=True`` to prevent such a behaviour:

.. ipython:: python

    cat = pd.Categorical([1,2,3,10], levels=[1,2,3,4,10])
    s = pd.Series(cat, name="cat", copy=True)
    cat
    s.iloc[0:2] = 10
    cat

.. note::
    This also happens in some cases when you supply a `numpy` array instea dof a `Categorical`:
    using an int array (e.g. ``np.array([1,2,3,4])``) will exhibit the same behaviour, but using
    a string array (e.g. ``np.array(["a","b","c","a"])``) will not.


Danger of confusion
~~~~~~~~~~~~~~~~~~~

Both `Series` and `Categorical` have a method ``.reorder_levels()`` but for different things. For
Series of type ``category`` this means that there is some danger to confuse both methods.

.. ipython:: python

    s = pd.Series(pd.Categorical([1,2,3,4]))
    print(s.cat.levels)
    # wrong and raises an error:
    try:
        s.reorder_levels([4,3,2,1])
    except Exception as e:
        print("Exception: " + str(e))
    # right
    s.cat.reorder_levels([4,3,2,1])
    print(s.cat.levels)

See also the API documentation for :func:`pandas.Series.reorder_levels` and
:func:`pandas.Categorical.reorder_levels`

Old style constructor usage
~~~~~~~~~~~~~~~~~~~~~~~~~~~

I earlier versions, a `Categorical` could be constructed by passing in precomputed `codes`
(called then `labels`) instead of values with levels. The `codes` are interpreted as pointers
to the levels with `-1` as `NaN`. This usage is now deprecated and not available unless
``compat=True`` is passed to the constructor of `Categorical`.

.. ipython:: python
    :okwarning:

    # This raises a FutureWarning:
    cat = pd.Categorical([1,2], levels=[1,2,3], compat=True)
    cat.get_values()

In the default case (``compat=False``) the first argument is interpreted as values.

.. ipython:: python

    cat = pd.Categorical([1,2], levels=[1,2,3], compat=False)
    cat.get_values()

.. warning::
    Using Categorical with precomputed codes and levels is deprecated and a `FutureWarning`
    is raised. Please change your code to use the :func:`~pandas.Categorical.from_codes`
    constructor instead of adding ``compat=False``.

No categorical index
~~~~~~~~~~~~~~~~~~~~

There is currently no index of type ``category``, so setting the index to a `Categorical` will
convert the `Categorical` to a normal `numpy` array first and therefore remove any custom
ordering of the levels:

.. ipython:: python

    cats = pd.Categorical([1,2,3,4], levels=[4,2,3,1])
    strings = ["a","b","c","d"]
    values = [4,2,3,1]
    df = pd.DataFrame({"strings":strings, "values":values}, index=cats)
    df.index
    # This should sort by levels but does not as there is no CategoricalIndex!
    df.sort_index()

.. note::
    This could change if a `CategoricalIndex` is implemented (see
    https://github.com/pydata/pandas/issues/7629)

dtype in apply
~~~~~~~~~~~~~~

Pandas currently does not preserve the dtype in apply functions: If you apply along rows you get
a `Series` of ``object`` `dtype` (same as getting a row -> getting one element will return a
basic type) and applying along columns will also convert to object.

.. ipython:: python

    df = pd.DataFrame({"a":[1,2,3,4], "b":["a","b","c","d"], "cats":pd.Categorical([1,2,3,2])})
    df.apply(lambda row: type(row["cats"]), axis=1)
    df.apply(lambda col: col.dtype, axis=0)


Future compatibility
~~~~~~~~~~~~~~~~~~~~

As `Categorical` is not a native `numpy` dtype, the implementation details of
`Series.cat` can change if such a `numpy` dtype is implemented.
