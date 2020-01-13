.. _missing_data:

{{ header }}

*************************
Working with missing data
*************************

In this section, we will discuss missing (also referred to as NA) values in
pandas.

.. note::

    The choice of using ``NaN`` internally to denote missing data was largely
    for simplicity and performance reasons.
    Starting from pandas 1.0, some optional data types start experimenting
    with a native ``NA`` scalar using a mask-based approach. See
    :ref:`here <missing_data.NA>` for more.

See the :ref:`cookbook<cookbook.missing_data>` for some advanced strategies.

Values considered "missing"
~~~~~~~~~~~~~~~~~~~~~~~~~~~

As data comes in many shapes and forms, pandas aims to be flexible with regard
to handling missing data. While ``NaN`` is the default missing value marker for
reasons of computational speed and convenience, we need to be able to easily
detect this value with data of different types: floating point, integer,
boolean, and general object. In many cases, however, the Python ``None`` will
arise and we wish to also consider that "missing" or "not available" or "NA".

.. note::

   If you want to consider ``inf`` and ``-inf`` to be "NA" in computations,
   you can set ``pandas.options.mode.use_inf_as_na = True``.

.. _missing.isna:

.. ipython:: python

   df = pd.DataFrame(np.random.randn(5, 3), index=['a', 'c', 'e', 'f', 'h'],
                     columns=['one', 'two', 'three'])
   df['four'] = 'bar'
   df['five'] = df['one'] > 0
   df
   df2 = df.reindex(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])
   df2

To make detecting missing values easier (and across different array dtypes),
pandas provides the :func:`isna` and
:func:`notna` functions, which are also methods on
Series and DataFrame objects:

.. ipython:: python

   df2['one']
   pd.isna(df2['one'])
   df2['four'].notna()
   df2.isna()

.. warning::

   One has to be mindful that in Python (and NumPy), the ``nan's`` don't compare equal, but ``None's`` **do**.
   Note that pandas/NumPy uses the fact that ``np.nan != np.nan``, and treats ``None`` like ``np.nan``.

   .. ipython:: python

      None == None                                                 # noqa: E711
      np.nan == np.nan

   So as compared to above, a scalar equality comparison versus a ``None/np.nan`` doesn't provide useful information.

   .. ipython:: python

      df2['one'] == np.nan

Integer dtypes and missing data
-------------------------------

Because ``NaN`` is a float, a column of integers with even one missing values
is cast to floating-point dtype (see :ref:`gotchas.intna` for more). Pandas
provides a nullable integer array, which can be used by explicitly requesting
the dtype:

.. ipython:: python

   pd.Series([1, 2, np.nan, 4], dtype=pd.Int64Dtype())

Alternatively, the string alias ``dtype='Int64'`` (note the capital ``"I"``) can be
used.

See :ref:`integer_na` for more.

Datetimes
---------

For datetime64[ns] types, ``NaT`` represents missing values. This is a pseudo-native
sentinel value that can be represented by NumPy in a singular dtype (datetime64[ns]).
pandas objects provide compatibility between ``NaT`` and ``NaN``.

.. ipython:: python

   df2 = df.copy()
   df2['timestamp'] = pd.Timestamp('20120101')
   df2
   df2.loc[['a', 'c', 'h'], ['one', 'timestamp']] = np.nan
   df2
   df2.dtypes.value_counts()

.. _missing.inserting:

Inserting missing data
~~~~~~~~~~~~~~~~~~~~~~

You can insert missing values by simply assigning to containers. The
actual missing value used will be chosen based on the dtype.

For example, numeric containers will always use ``NaN`` regardless of
the missing value type chosen:

.. ipython:: python

   s = pd.Series([1, 2, 3])
   s.loc[0] = None
   s

Likewise, datetime containers will always use ``NaT``.

For object containers, pandas will use the value given:

.. ipython:: python

   s = pd.Series(["a", "b", "c"])
   s.loc[0] = None
   s.loc[1] = np.nan
   s

.. _missing_data.calculations:

Calculations with missing data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Missing values propagate naturally through arithmetic operations between pandas
objects.

.. ipython:: python
   :suppress:

   df = df2.loc[:, ['one', 'two', 'three']]
   a = df2.loc[df2.index[:5], ['one', 'two']].fillna(method='pad')
   b = df2.loc[df2.index[:5], ['one', 'two', 'three']]

.. ipython:: python

   a
   b
   a + b

The descriptive statistics and computational methods discussed in the
:ref:`data structure overview <basics.stats>` (and listed :ref:`here
<api.series.stats>` and :ref:`here <api.dataframe.stats>`) are all written to
account for missing data. For example:

* When summing data, NA (missing) values will be treated as zero.
* If the data are all NA, the result will be 0.
* Cumulative methods like :meth:`~DataFrame.cumsum` and :meth:`~DataFrame.cumprod` ignore NA values by default, but preserve them in the resulting arrays. To override this behaviour and include NA values, use ``skipna=False``.

.. ipython:: python

   df
   df['one'].sum()
   df.mean(1)
   df.cumsum()
   df.cumsum(skipna=False)


.. _missing_data.numeric_sum:

Sum/prod of empties/nans
~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   This behavior is now standard as of v0.22.0 and is consistent with the default in ``numpy``; previously sum/prod of all-NA or empty Series/DataFrames would return NaN.
   See :ref:`v0.22.0 whatsnew <whatsnew_0220>` for more.

The sum of an empty or all-NA Series or column of a DataFrame is 0.

.. ipython:: python

   pd.Series([np.nan]).sum()

   pd.Series([], dtype="float64").sum()

The product of an empty or all-NA Series or column of a DataFrame is 1.

.. ipython:: python

   pd.Series([np.nan]).prod()

   pd.Series([], dtype="float64").prod()


NA values in GroupBy
~~~~~~~~~~~~~~~~~~~~

NA groups in GroupBy are automatically excluded. This behavior is consistent
with R, for example:

.. ipython:: python

    df
    df.groupby('one').mean()

See the groupby section :ref:`here <groupby.missing>` for more information.

Cleaning / filling missing data
--------------------------------

pandas objects are equipped with various data manipulation methods for dealing
with missing data.

.. _missing_data.fillna:

Filling missing values: fillna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`~DataFrame.fillna` can "fill in" NA values with non-NA data in a couple
of ways, which we illustrate:

**Replace NA with a scalar value**

.. ipython:: python

   df2
   df2.fillna(0)
   df2['one'].fillna('missing')

**Fill gaps forward or backward**

Using the same filling arguments as :ref:`reindexing <basics.reindexing>`, we
can propagate non-NA values forward or backward:

.. ipython:: python

   df
   df.fillna(method='pad')

.. _missing_data.fillna.limit:

**Limit the amount of filling**

If we only want consecutive gaps filled up to a certain number of data points,
we can use the `limit` keyword:

.. ipython:: python
   :suppress:

   df.iloc[2:4, :] = np.nan

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

:meth:`~DataFrame.ffill` is equivalent to ``fillna(method='ffill')``
and :meth:`~DataFrame.bfill` is equivalent to ``fillna(method='bfill')``

.. _missing_data.PandasObject:

Filling with a PandasObject
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also fillna using a dict or Series that is alignable. The labels of the dict or index of the Series
must match the columns of the frame you wish to fill. The
use case of this is to fill a DataFrame with the mean of that column.

.. ipython:: python

        dff = pd.DataFrame(np.random.randn(10, 3), columns=list('ABC'))
        dff.iloc[3:5, 0] = np.nan
        dff.iloc[4:6, 1] = np.nan
        dff.iloc[5:8, 2] = np.nan
        dff

        dff.fillna(dff.mean())
        dff.fillna(dff.mean()['B':'C'])

Same result as above, but is aligning the 'fill' value which is
a Series in this case.

.. ipython:: python

        dff.where(pd.notna(dff), dff.mean(), axis='columns')


.. _missing_data.dropna:

Dropping axis labels with missing data: dropna
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may wish to simply exclude labels from a data set which refer to missing
data. To do this, use :meth:`~DataFrame.dropna`:

.. ipython:: python
   :suppress:

   df['two'] = df['two'].fillna(0)
   df['three'] = df['three'].fillna(0)

.. ipython:: python

   df
   df.dropna(axis=0)
   df.dropna(axis=1)
   df['one'].dropna()

An equivalent :meth:`~Series.dropna` is available for Series.
DataFrame.dropna has considerably more options than Series.dropna, which can be
examined :ref:`in the API <api.dataframe.missing>`.

.. _missing_data.interpolate:

Interpolation
~~~~~~~~~~~~~

.. versionadded:: 0.23.0

  The ``limit_area`` keyword argument was added.

Both Series and DataFrame objects have :meth:`~DataFrame.interpolate`
that, by default, performs linear interpolation at missing data points.

.. ipython:: python
   :suppress:

   np.random.seed(123456)
   idx = pd.date_range('1/1/2000', periods=100, freq='BM')
   ts = pd.Series(np.random.randn(100), index=idx)
   ts[1:5] = np.nan
   ts[20:30] = np.nan
   ts[60:80] = np.nan
   ts = ts.cumsum()

.. ipython:: python

   ts
   ts.count()
   @savefig series_before_interpolate.png
   ts.plot()

.. ipython:: python

   ts.interpolate()
   ts.interpolate().count()

   @savefig series_interpolate.png
   ts.interpolate().plot()

Index aware interpolation is available via the ``method`` keyword:

.. ipython:: python
   :suppress:

   ts2 = ts[[0, 1, 30, 60, 99]]

.. ipython:: python

   ts2
   ts2.interpolate()
   ts2.interpolate(method='time')

For a floating-point index, use ``method='values'``:

.. ipython:: python
   :suppress:

   idx = [0., 1., 10.]
   ser = pd.Series([0., np.nan, 10.], idx)

.. ipython:: python

   ser
   ser.interpolate()
   ser.interpolate(method='values')

You can also interpolate with a DataFrame:

.. ipython:: python

   df = pd.DataFrame({'A': [1, 2.1, np.nan, 4.7, 5.6, 6.8],
                      'B': [.25, np.nan, np.nan, 4, 12.2, 14.4]})
   df
   df.interpolate()

The ``method`` argument gives access to fancier interpolation methods.
If you have scipy_ installed, you can pass the name of a 1-d interpolation routine to ``method``.
You'll want to consult the full scipy interpolation documentation_ and reference guide_ for details.
The appropriate interpolation method will depend on the type of data you are working with.

* If you are dealing with a time series that is growing at an increasing rate,
  ``method='quadratic'`` may be appropriate.
* If you have values approximating a cumulative distribution function,
  then ``method='pchip'`` should work well.
* To fill missing values with goal of smooth plotting, consider ``method='akima'``.

.. warning::

   These methods require ``scipy``.

.. ipython:: python

   df.interpolate(method='barycentric')

   df.interpolate(method='pchip')

   df.interpolate(method='akima')

When interpolating via a polynomial or spline approximation, you must also specify
the degree or order of the approximation:

.. ipython:: python

   df.interpolate(method='spline', order=2)

   df.interpolate(method='polynomial', order=2)

Compare several methods:

.. ipython:: python

   np.random.seed(2)

   ser = pd.Series(np.arange(1, 10.1, .25) ** 2 + np.random.randn(37))
   missing = np.array([4, 13, 14, 15, 16, 17, 18, 20, 29])
   ser[missing] = np.nan
   methods = ['linear', 'quadratic', 'cubic']

   df = pd.DataFrame({m: ser.interpolate(method=m) for m in methods})
   @savefig compare_interpolations.png
   df.plot()

Another use case is interpolation at *new* values.
Suppose you have 100 observations from some distribution. And let's suppose
that you're particularly interested in what's happening around the middle.
You can mix pandas' ``reindex`` and ``interpolate`` methods to interpolate
at the new values.

.. ipython:: python

   ser = pd.Series(np.sort(np.random.uniform(size=100)))

   # interpolate at new_index
   new_index = ser.index | pd.Index([49.25, 49.5, 49.75, 50.25, 50.5, 50.75])
   interp_s = ser.reindex(new_index).interpolate(method='pchip')
   interp_s[49:51]

.. _scipy: http://www.scipy.org
.. _documentation: http://docs.scipy.org/doc/scipy/reference/interpolate.html#univariate-interpolation
.. _guide: http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html

.. _missing_data.interp_limits:

Interpolation limits
--------------------

Like other pandas fill methods, :meth:`~DataFrame.interpolate` accepts a ``limit`` keyword
argument. Use this argument to limit the number of consecutive ``NaN`` values
filled since the last valid observation:

.. ipython:: python

   ser = pd.Series([np.nan, np.nan, 5, np.nan, np.nan,
                    np.nan, 13, np.nan, np.nan])
   ser

   # fill all consecutive values in a forward direction
   ser.interpolate()

   # fill one consecutive value in a forward direction
   ser.interpolate(limit=1)

By default, ``NaN`` values are filled in a ``forward`` direction. Use
``limit_direction`` parameter to fill ``backward`` or from ``both`` directions.

.. ipython:: python

   # fill one consecutive value backwards
   ser.interpolate(limit=1, limit_direction='backward')

   # fill one consecutive value in both directions
   ser.interpolate(limit=1, limit_direction='both')

   # fill all consecutive values in both directions
   ser.interpolate(limit_direction='both')

By default, ``NaN`` values are filled whether they are inside (surrounded by)
existing valid values, or outside existing valid values. Introduced in v0.23
the ``limit_area`` parameter restricts filling to either inside or outside values.

.. ipython:: python

   # fill one consecutive inside value in both directions
   ser.interpolate(limit_direction='both', limit_area='inside', limit=1)

   # fill all consecutive outside values backward
   ser.interpolate(limit_direction='backward', limit_area='outside')

   # fill all consecutive outside values in both directions
   ser.interpolate(limit_direction='both', limit_area='outside')

.. _missing_data.replace:

Replacing generic values
~~~~~~~~~~~~~~~~~~~~~~~~
Often times we want to replace arbitrary values with other values.

:meth:`~Series.replace` in Series and :meth:`~DataFrame.replace` in DataFrame provides an efficient yet
flexible way to perform such replacements.

For a Series, you can replace a single value or a list of values by another
value:

.. ipython:: python

   ser = pd.Series([0., 1., 2., 3., 4.])

   ser.replace(0, 5)

You can replace a list of values by a list of other values:

.. ipython:: python

   ser.replace([0, 1, 2, 3, 4], [4, 3, 2, 1, 0])

You can also specify a mapping dict:

.. ipython:: python

   ser.replace({0: 10, 1: 100})

For a DataFrame, you can specify individual values by column:

.. ipython:: python

   df = pd.DataFrame({'a': [0, 1, 2, 3, 4], 'b': [5, 6, 7, 8, 9]})

   df.replace({'a': 0, 'b': 5}, 100)

Instead of replacing with specified values, you can treat all given values as
missing and interpolate over them:

.. ipython:: python

   ser.replace([1, 2, 3], method='pad')

.. _missing_data.replace_expression:

String/regular expression replacement
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   Python strings prefixed with the ``r`` character such as ``r'hello world'``
   are so-called "raw" strings. They have different semantics regarding
   backslashes than strings without this prefix. Backslashes in raw strings
   will be interpreted as an escaped backslash, e.g., ``r'\' == '\\'``. You
   should `read about them
   <https://docs.python.org/3/reference/lexical_analysis.html#string-literals>`__
   if this is unclear.

Replace the '.' with ``NaN`` (str -> str):

.. ipython:: python

   d = {'a': list(range(4)), 'b': list('ab..'), 'c': ['a', 'b', np.nan, 'd']}
   df = pd.DataFrame(d)
   df.replace('.', np.nan)

Now do it with a regular expression that removes surrounding whitespace
(regex -> regex):

.. ipython:: python

   df.replace(r'\s*\.\s*', np.nan, regex=True)

Replace a few different values (list -> list):

.. ipython:: python

   df.replace(['a', '.'], ['b', np.nan])

list of regex -> list of regex:

.. ipython:: python

   df.replace([r'\.', r'(a)'], ['dot', r'\1stuff'], regex=True)

Only search in column ``'b'`` (dict -> dict):

.. ipython:: python

   df.replace({'b': '.'}, {'b': np.nan})

Same as the previous example, but use a regular expression for
searching instead (dict of regex -> dict):

.. ipython:: python

   df.replace({'b': r'\s*\.\s*'}, {'b': np.nan}, regex=True)

You can pass nested dictionaries of regular expressions that use ``regex=True``:

.. ipython:: python

   df.replace({'b': {'b': r''}}, regex=True)

Alternatively, you can pass the nested dictionary like so:

.. ipython:: python

   df.replace(regex={'b': {r'\s*\.\s*': np.nan}})

You can also use the group of a regular expression match when replacing (dict
of regex -> dict of regex), this works for lists as well.

.. ipython:: python

   df.replace({'b': r'\s*(\.)\s*'}, {'b': r'\1ty'}, regex=True)

You can pass a list of regular expressions, of which those that match
will be replaced with a scalar (list of regex -> regex).

.. ipython:: python

   df.replace([r'\s*\.\s*', r'a|b'], np.nan, regex=True)

All of the regular expression examples can also be passed with the
``to_replace`` argument as the ``regex`` argument. In this case the ``value``
argument must be passed explicitly by name or ``regex`` must be a nested
dictionary. The previous example, in this case, would then be:

.. ipython:: python

   df.replace(regex=[r'\s*\.\s*', r'a|b'], value=np.nan)

This can be convenient if you do not want to pass ``regex=True`` every time you
want to use a regular expression.

.. note::

   Anywhere in the above ``replace`` examples that you see a regular expression
   a compiled regular expression is valid as well.

Numeric replacement
~~~~~~~~~~~~~~~~~~~

:meth:`~DataFrame.replace` is similar to :meth:`~DataFrame.fillna`.

.. ipython:: python

   df = pd.DataFrame(np.random.randn(10, 2))
   df[np.random.rand(df.shape[0]) > 0.5] = 1.5
   df.replace(1.5, np.nan)

Replacing more than one value is possible by passing a list.

.. ipython:: python

   df00 = df.iloc[0, 0]
   df.replace([1.5, df00], [np.nan, 'a'])
   df[1].dtype

You can also operate on the DataFrame in place:

.. ipython:: python

   df.replace(1.5, np.nan, inplace=True)

.. warning::

   When replacing multiple ``bool`` or ``datetime64`` objects, the first
   argument to ``replace`` (``to_replace``) must match the type of the value
   being replaced. For example,

   .. code-block:: python

      >>> s = pd.Series([True, False, True])
      >>> s.replace({'a string': 'new value', True: False})  # raises
      TypeError: Cannot compare types 'ndarray(dtype=bool)' and 'str'

   will raise a ``TypeError`` because one of the ``dict`` keys is not of the
   correct type for replacement.

   However, when replacing a *single* object such as,

   .. ipython:: python

      s = pd.Series([True, False, True])
      s.replace('a string', 'another string')

   the original ``NDFrame`` object will be returned untouched. We're working on
   unifying this API, but for backwards compatibility reasons we cannot break
   the latter behavior. See :issue:`6354` for more details.

Missing data casting rules and indexing
---------------------------------------

While pandas supports storing arrays of integer and boolean type, these types
are not capable of storing missing data. Until we can switch to using a native
NA type in NumPy, we've established some "casting rules". When a reindexing
operation introduces missing data, the Series will be cast according to the
rules introduced in the table below.

.. csv-table::
    :header: "data type", "Cast to"
    :widths: 40, 40

    integer, float
    boolean, object
    float, no cast
    object, no cast

For example:

.. ipython:: python

   s = pd.Series(np.random.randn(5), index=[0, 2, 4, 6, 7])
   s > 0
   (s > 0).dtype
   crit = (s > 0).reindex(list(range(8)))
   crit
   crit.dtype

Ordinarily NumPy will complain if you try to use an object array (even if it
contains boolean values) instead of a boolean array to get or set values from
an ndarray (e.g. selecting values based on some criteria). If a boolean vector
contains NAs, an exception will be generated:

.. ipython:: python
   :okexcept:

   reindexed = s.reindex(list(range(8))).fillna(0)
   reindexed[crit]

However, these can be filled in using :meth:`~DataFrame.fillna` and it will work fine:

.. ipython:: python

   reindexed[crit.fillna(False)]
   reindexed[crit.fillna(True)]

Pandas provides a nullable integer dtype, but you must explicitly request it
when creating the series or column. Notice that we use a capital "I" in
the ``dtype="Int64"``.

.. ipython:: python

   s = pd.Series([0, 1, np.nan, 3, 4], dtype="Int64")
   s

See :ref:`integer_na` for more.


.. _missing_data.NA:

Experimental ``NA`` scalar to denote missing values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   Experimental: the behaviour of ``pd.NA`` can still change without warning.

.. versionadded:: 1.0.0

Starting from pandas 1.0, an experimental ``pd.NA`` value (singleton) is
available to represent scalar missing values. At this moment, it is used in
the nullable :doc:`integer <integer_na>`, boolean and
:ref:`dedicated string <text.types>` data types as the missing value indicator.

The goal of ``pd.NA`` is provide a "missing" indicator that can be used
consistently across data types (instead of ``np.nan``, ``None`` or ``pd.NaT``
depending on the data type).

For example, when having missing values in a Series with the nullable integer
dtype, it will use ``pd.NA``:

.. ipython:: python

    s = pd.Series([1, 2, None], dtype="Int64")
    s
    s[2]
    s[2] is pd.NA

Currently, pandas does not yet use those data types by default (when creating
a DataFrame or Series, or when reading in data), so you need to specify
the dtype explicitly.

Propagation in arithmetic and comparison operations
---------------------------------------------------

In general, missing values *propagate* in operations involving ``pd.NA``. When
one of the operands is unknown, the outcome of the operation is also unknown.

For example, ``pd.NA`` propagates in arithmetic operations, similarly to
``np.nan``:

.. ipython:: python

   pd.NA + 1
   "a" * pd.NA

There are a few special cases when the result is known, even when one of the
operands is ``NA``.


================ ======
Operation        Result
================ ======
``pd.NA ** 0``   0
``1 ** pd.NA``   1
================ ======

In equality and comparison operations, ``pd.NA`` also propagates. This deviates
from the behaviour of ``np.nan``, where comparisons with ``np.nan`` always
return ``False``.

.. ipython:: python

   pd.NA == 1
   pd.NA == pd.NA
   pd.NA < 2.5

To check if a value is equal to ``pd.NA``, the :func:`isna` function can be
used:

.. ipython:: python

   pd.isna(pd.NA)

An exception on this basic propagation rule are *reductions* (such as the
mean or the minimum), where pandas defaults to skipping missing values. See
:ref:`above <missing_data.calculations>` for more.

Logical operations
------------------

For logical operations, ``pd.NA`` follows the rules of the
`three-valued logic <https://en.wikipedia.org/wiki/Three-valued_logic>`__ (or
*Kleene logic*, similarly to R, SQL and Julia). This logic means to only
propagate missing values when it is logically required.

For example, for the logical "or" operation (``|``), if one of the operands
is ``True``, we already know the result will be ``True``, regardless of the
other value (so regardless the missing value would be ``True`` or ``False``).
In this case, ``pd.NA`` does not propagate:

.. ipython:: python

   True | False
   True | pd.NA
   pd.NA | True

On the other hand, if one of the operands is ``False``, the result depends
on the value of the other operand. Therefore, in this case ``pd.NA``
propagates:

.. ipython:: python

   False | True
   False | False
   False | pd.NA

The behaviour of the logical "and" operation (``&``) can be derived using
similar logic (where now ``pd.NA`` will not propagate if one of the operands
is already ``False``):

.. ipython:: python

   False & True
   False & False
   False & pd.NA

.. ipython:: python

   True & True
   True & False
   True & pd.NA


``NA`` in a boolean context
---------------------------

Since the actual value of an NA is unknown, it is ambiguous to convert NA
to a boolean value. The following raises an error:

.. ipython:: python
   :okexcept:

   bool(pd.NA)

This also means that ``pd.NA`` cannot be used in a context where it is
evaluated to a boolean, such as ``if condition: ...`` where ``condition`` can
potentially be ``pd.NA``. In such cases, :func:`isna` can be used to check
for ``pd.NA`` or ``condition`` being ``pd.NA`` can be avoided, for example by
filling missing values beforehand.

A similar situation occurs when using Series or DataFrame objects in ``if``
statements, see :ref:`gotchas.truth`.

NumPy ufuncs
------------

:attr:`pandas.NA` implements NumPy's ``__array_ufunc__`` protocol. Most ufuncs
work with ``NA``, and generally return ``NA``:

.. ipython:: python

   np.log(pd.NA)
   np.add(pd.NA, 1)

.. warning::

   Currently, ufuncs involving an ndarray and ``NA`` will return an
   object-dtype filled with NA values.

   .. ipython:: python

      a = np.array([1, 2, 3])
      np.greater(a, pd.NA)

   The return type here may change to return a different array type
   in the future.

See :ref:`dsintro.numpy_interop` for more on ufuncs.
