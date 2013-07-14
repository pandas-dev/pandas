.. currentmodule:: pandas
.. _missing_data:

*************************
Working with missing data
*************************

In this section, we will discuss missing (also referred to as NA) values in
pandas.

.. ipython:: python
   :suppress:

   import numpy as np; randn = np.random.randn; randint =np.random.randint
   from pandas import *
   import matplotlib.pyplot as plt

.. note::

    The choice of using ``NaN`` internally to denote missing data was largely
    for simplicity and performance reasons. It differs from the MaskedArray
    approach of, for example, :mod:`scikits.timeseries`. We are hopeful that
    NumPy will soon be able to provide a native NA type solution (similar to R)
    performant enough to be used in pandas.

See the :ref:`cookbook<cookbook.missing_data>` for some advanced strategies

Missing data basics
-------------------

When / why does data become missing?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Some might quibble over our usage of *missing*. By "missing" we simply mean
**null** or "not present for whatever reason". Many data sets simply arrive with
missing data, either because it exists and was not collected or it never
existed. For example, in a collection of financial time series, some of the time
series might start on different dates. Thus, values prior to the start date
would generally be marked as missing.

In pandas, one of the most common ways that missing data is **introduced** into
a data set is by reindexing. For example

.. ipython:: python

   df = DataFrame(randn(5, 3), index=['a', 'c', 'e', 'f', 'h'],
                  columns=['one', 'two', 'three'])
   df['four'] = 'bar'
   df['five'] = df['one'] > 0
   df
   df2 = df.reindex(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
   df2

Values considered "missing"
~~~~~~~~~~~~~~~~~~~~~~~~~~~

As data comes in many shapes and forms, pandas aims to be flexible with regard
to handling missing data. While ``NaN`` is the default missing value marker for
reasons of computational speed and convenience, we need to be able to easily
detect this value with data of different types: floating point, integer,
boolean, and general object. In many cases, however, the Python ``None`` will
arise and we wish to also consider that "missing" or "null".

Until recently, for legacy reasons ``inf`` and ``-inf`` were also
considered to be "null" in computations. This is no longer the case by
default; use the ``mode.use_inf_as_null`` option to recover it.

.. _missing.isnull:

To make detecting missing values easier (and across different array dtypes),
pandas provides the :func:`~pandas.core.common.isnull` and
:func:`~pandas.core.common.notnull` functions, which are also methods on
``Series`` objects:

.. ipython:: python

   df2['one']
   isnull(df2['one'])
   df2['four'].notnull()

**Summary:** ``NaN`` and ``None`` (in object arrays) are considered
missing by the ``isnull`` and ``notnull`` functions. ``inf`` and
``-inf`` are no longer considered missing by default.

Datetimes
---------

For datetime64[ns] types, ``NaT`` represents missing values. This is a pseudo-native
sentinal value that can be represented by numpy in a singular dtype (datetime64[ns]).
Pandas objects provide intercompatibility between ``NaT`` and ``NaN``.

.. ipython:: python

   df2 = df.copy()
   df2['timestamp'] = Timestamp('20120101')
   df2
   df2.ix[['a','c','h'],['one','timestamp']] = np.nan
   df2
   df2.get_dtype_counts()


Calculations with missing data
------------------------------

Missing values propagate naturally through arithmetic operations between pandas
objects.

.. ipython:: python
   :suppress:

   df = df2.ix[:, ['one', 'two', 'three']]
   a = df2.ix[:5, ['one', 'two']].fillna(method='pad')
   b = df2.ix[:5, ['one', 'two', 'three']]

.. ipython:: python

   a
   b
   a + b

The descriptive statistics and computational methods discussed in the
:ref:`data structure overview <basics.stats>` (and listed :ref:`here
<api.series.stats>` and :ref:`here <api.dataframe.stats>`) are all written to
account for missing data. For example:

  * When summing data, NA (missing) values will be treated as zero
  * If the data are all NA, the result will be NA
  * Methods like **cumsum** and **cumprod** ignore NA values, but preserve them
    in the resulting arrays

.. ipython:: python

   df
   df['one'].sum()
   df.mean(1)
   df.cumsum()

NA values in GroupBy
~~~~~~~~~~~~~~~~~~~~

NA groups in GroupBy are automatically excluded. This behavior is consistent
with R, for example.



Cleaning / filling missing data
--------------------------------

pandas objects are equipped with various data manipulation methods for dealing
with missing data.

.. _missing_data.fillna:

Filling missing values: fillna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The **fillna** function can "fill in" NA values with non-null data in a couple
of ways, which we illustrate:

**Replace NA with a scalar value**

.. ipython:: python

   df2
   df2.fillna(0)
   df2['four'].fillna('missing')

**Fill gaps forward or backward**

Using the same filling arguments as :ref:`reindexing <basics.reindexing>`, we
can propagate non-null values forward or backward:

.. ipython:: python

   df
   df.fillna(method='pad')

.. _missing_data.fillna.limit:

**Limit the amount of filling**

If we only want consecutive gaps filled up to a certain number of data points,
we can use the `limit` keyword:

.. ipython:: python
   :suppress:

   df.ix[2:4, :] = np.nan

.. ipython:: python

   df
   df.fillna(method='pad', limit=1)

To remind you, these are the available filling methods:

.. csv-table::
    :header: "Method", "Action"
    :widths: 30, 50

    pad / ffill, Fill values forward
    bfill / backfill, Fill values backward

With time series data, using pad/ffill is extremely common so that the "last
known value" is available at every time point.

.. _missing_data.dropna:

Dropping axis labels with missing data: dropna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may wish to simply exclude labels from a data set which refer to missing
data. To do this, use the **dropna** method:

.. ipython:: python
   :suppress:

   df['two'] = df['two'].fillna(0)
   df['three'] = df['three'].fillna(0)

.. ipython:: python

   df
   df.dropna(axis=0)
   df.dropna(axis=1)
   df['one'].dropna()

**dropna** is presently only implemented for Series and DataFrame, but will be
eventually added to Panel. Series.dropna is a simpler method as it only has one
axis to consider. DataFrame.dropna has considerably more options, which can be
examined :ref:`in the API <api.dataframe.missing>`.

.. _missing_data.interpolate:

Interpolation
~~~~~~~~~~~~~

A linear **interpolate** method has been implemented on Series. The default
interpolation assumes equally spaced points.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   idx = date_range('1/1/2000', periods=100, freq='BM')
   ts = Series(randn(100), index=idx)
   ts[1:20] = np.nan
   ts[60:80] = np.nan
   ts = ts.cumsum()

.. ipython:: python

   ts.count()

   ts.head()

   ts.interpolate().count()

   ts.interpolate().head()

   @savefig series_interpolate.png
   ts.interpolate().plot()

Index aware interpolation is available via the ``method`` keyword:

.. ipython:: python
   :suppress:

   ts = ts[[0, 1, 30, 60, 99]]

.. ipython:: python

   ts

   ts.interpolate()

   ts.interpolate(method='time')

For a floating-point index, use ``method='values'``:

.. ipython:: python
   :suppress:

   idx = [0., 1., 10.]
   ser = Series([0., np.nan, 10.], idx)

.. ipython:: python

   ser

   ser.interpolate()

   ser.interpolate(method='values')

.. _missing_data.replace:

Replacing Generic Values
~~~~~~~~~~~~~~~~~~~~~~~~
Often times we want to replace arbitrary values with other values. New in v0.8
is the ``replace`` method in Series/DataFrame that provides an efficient yet
flexible way to perform such replacements.

For a Series, you can replace a single value or a list of values by another
value:

.. ipython:: python

   ser = Series([0., 1., 2., 3., 4.])

   ser.replace(0, 5)

You can replace a list of values by a list of other values:

.. ipython:: python

   ser.replace([0, 1, 2, 3, 4], [4, 3, 2, 1, 0])

You can also specify a mapping dict:

.. ipython:: python

   ser.replace({0: 10, 1: 100})

For a DataFrame, you can specify individual values by column:

.. ipython:: python

   df = DataFrame({'a': [0, 1, 2, 3, 4], 'b': [5, 6, 7, 8, 9]})

   df.replace({'a': 0, 'b': 5}, 100)

Instead of replacing with specified values, you can treat all given values as
missing and interpolate over them:

.. ipython:: python

   ser.replace([1, 2, 3], method='pad')

.. _missing_data.replace_expression:

String/Regular Expression Replacement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   Python strings prefixed with the ``r`` character such as ``r'hello world'``
   are so-called "raw" strings. They have different semantics regarding
   backslashes than strings without this prefix. Backslashes in raw strings
   will be interpreted as an escaped backslash, e.g., ``r'\' == '\\'``. You
   should `read about them
   <http://docs.python.org/2/reference/lexical_analysis.html#string-literals>`_
   if this is unclear.

Replace the '.' with ``nan`` (str -> str)

.. ipython:: python
   :suppress:

   from numpy.random import rand, randn
   from numpy import nan
   from pandas import DataFrame

.. ipython:: python

   d = {'a': range(4), 'b': list('ab..'), 'c': ['a', 'b', nan, 'd']}
   df = DataFrame(d)
   df.replace('.', nan)

Now do it with a regular expression that removes surrounding whitespace
(regex -> regex)

.. ipython:: python

   df.replace(r'\s*\.\s*', nan, regex=True)

Replace a few different values (list -> list)

.. ipython:: python

   df.replace(['a', '.'], ['b', nan])

list of regex -> list of regex

.. ipython:: python

   df.replace([r'\.', r'(a)'], ['dot', '\1stuff'], regex=True)

Only search in column ``'b'`` (dict -> dict)

.. ipython:: python

   df.replace({'b': '.'}, {'b': nan})

Same as the previous example, but use a regular expression for
searching instead (dict of regex -> dict)

.. ipython:: python

   df.replace({'b': r'\s*\.\s*'}, {'b': nan}, regex=True)

You can pass nested dictionaries of regular expressions that use ``regex=True``

.. ipython:: python

   df.replace({'b': {'b': r''}}, regex=True)

or you can pass the nested dictionary like so

.. ipython:: python

   df.replace(regex={'b': {r'\s*\.\s*': nan}})

You can also use the group of a regular expression match when replacing (dict
of regex -> dict of regex), this works for lists as well

.. ipython:: python

   df.replace({'b': r'\s*(\.)\s*'}, {'b': r'\1ty'}, regex=True)

You can pass a list of regular expressions, of which those that match
will be replaced with a scalar (list of regex -> regex)

.. ipython:: python

   df.replace([r'\s*\.\s*', r'a|b'], nan, regex=True)

All of the regular expression examples can also be passed with the
``to_replace`` argument as the ``regex`` argument. In this case the ``value``
argument must be passed explicity by name or ``regex`` must be a nested
dictionary. The previous example, in this case, would then be

.. ipython:: python

   df.replace(regex=[r'\s*\.\s*', r'a|b'], value=nan)

This can be convenient if you do not want to pass ``regex=True`` every time you
want to use a regular expression.

.. note::

   Anywhere in the above ``replace`` examples that you see a regular expression
   a compiled regular expression is valid as well.

Numeric Replacement
~~~~~~~~~~~~~~~~~~~

Similiar to ``DataFrame.fillna``

.. ipython:: python
   :suppress:

   from numpy.random import rand, randn
   from numpy import nan
   from pandas import DataFrame
   from pandas.util.testing import assert_frame_equal

.. ipython:: python

   df = DataFrame(randn(10, 2))
   df[rand(df.shape[0]) > 0.5] = 1.5
   df.replace(1.5, nan)

Replacing more than one value via lists works as well

.. ipython:: python

   df00 = df.values[0, 0]
   df.replace([1.5, df00], [nan, 'a'])
   df[1].dtype

You can also operate on the DataFrame in place

.. ipython:: python

   df.replace(1.5, nan, inplace=True)

Missing data casting rules and indexing
---------------------------------------

While pandas supports storing arrays of integer and boolean type, these types
are not capable of storing missing data. Until we can switch to using a native
NA type in NumPy, we've established some "casting rules" when reindexing will
cause missing data to be introduced into, say, a Series or DataFrame. Here they
are:

.. csv-table::
    :header: "data type", "Cast to"
    :widths: 40, 40

	integer, float
    boolean, object
    float, no cast
    object, no cast

For example:

.. ipython:: python

   s = Series(randn(5), index=[0, 2, 4, 6, 7])
   s > 0
   (s > 0).dtype
   crit = (s > 0).reindex(range(8))
   crit
   crit.dtype

Ordinarily NumPy will complain if you try to use an object array (even if it
contains boolean values) instead of a boolean array to get or set values from
an ndarray (e.g. selecting values based on some criteria). If a boolean vector
contains NAs, an exception will be generated:

.. ipython:: python
   :okexcept:

   reindexed = s.reindex(range(8)).fillna(0)
   reindexed[crit]

However, these can be filled in using **fillna** and it will work fine:

.. ipython:: python

   reindexed[crit.fillna(False)]
   reindexed[crit.fillna(True)]

