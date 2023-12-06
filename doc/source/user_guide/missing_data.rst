.. _missing_data:

{{ header }}

*************************
Working with missing data
*************************

Values considered "missing"
~~~~~~~~~~~~~~~~~~~~~~~~~~~

pandas uses different sentinel values to represent a missing (also referred to as NA)
depending on the data type.

``numpy.nan`` for NumPy data types. The disadvantage of using NumPy data types
is that the original data type will be coerced to ``np.float64`` or ``object``.

.. ipython:: python

   pd.Series([1, 2], dtype=np.int64).reindex([0, 1, 2])
   pd.Series([True, False], dtype=np.bool_).reindex([0, 1, 2])

:class:`NaT` for NumPy ``np.datetime64``, ``np.timedelta64``, and :class:`PeriodDtype`. For typing applications,
use :class:`api.types.NaTType`.

.. ipython:: python

   pd.Series([1, 2], dtype=np.dtype("timedelta64[ns]")).reindex([0, 1, 2])
   pd.Series([1, 2], dtype=np.dtype("datetime64[ns]")).reindex([0, 1, 2])
   pd.Series(["2020", "2020"], dtype=pd.PeriodDtype("D")).reindex([0, 1, 2])

:class:`NA` for :class:`StringDtype`, :class:`Int64Dtype` (and other bit widths),
:class:`Float64Dtype`(and other bit widths), :class:`BooleanDtype` and :class:`ArrowDtype`.
These types will maintain the original data type of the data.
For typing applications, use :class:`api.types.NAType`.

.. ipython:: python

   pd.Series([1, 2], dtype="Int64").reindex([0, 1, 2])
   pd.Series([True, False], dtype="boolean[pyarrow]").reindex([0, 1, 2])

To detect these missing value, use the :func:`isna` or :func:`notna` methods.

.. ipython:: python

   ser = pd.Series([pd.Timestamp("2020-01-01"), pd.NaT])
   ser
   pd.isna(ser)


.. note::

   :func:`isna` or :func:`notna` will also consider ``None`` a missing value.

   .. ipython:: python

      ser = pd.Series([1, None], dtype=object)
      ser
      pd.isna(ser)

.. warning::

   Equality compaisons between ``np.nan``, :class:`NaT`, and :class:`NA`
   do not act like ``None``

   .. ipython:: python

      None == None  # noqa: E711
      np.nan == np.nan
      pd.NaT == pd.NaT
      pd.NA == pd.NA

   Therefore, an equality comparison between a :class:`DataFrame` or :class:`Series`
   with one of these missing values does not provide the same information as
   :func:`isna` or :func:`notna`.

   .. ipython:: python

      ser = pd.Series([True, None], dtype="boolean[pyarrow]")
      ser == pd.NA
      pd.isna(ser)


.. _missing_data.NA:

:class:`NA` semantics
~~~~~~~~~~~~~~~~~~~~~

.. warning::

   Experimental: the behaviour of :class:`NA`` can still change without warning.

Starting from pandas 1.0, an experimental :class:`NA` value (singleton) is
available to represent scalar missing values. The goal of :class:`NA` is provide a
"missing" indicator that can be used consistently across data types
(instead of ``np.nan``, ``None`` or ``pd.NaT`` depending on the data type).

For example, when having missing values in a :class:`Series` with the nullable integer
dtype, it will use :class:`NA`:

.. ipython:: python

    s = pd.Series([1, 2, None], dtype="Int64")
    s
    s[2]
    s[2] is pd.NA

Currently, pandas does not yet use those data types using :class:`NA` by default
a :class:`DataFrame` or :class:`Series`, so you need to specify
the dtype explicitly. An easy way to convert to those dtypes is explained in the
:ref:`conversion section <missing_data.NA.conversion>`.

Propagation in arithmetic and comparison operations
---------------------------------------------------

In general, missing values *propagate* in operations involving :class:`NA`. When
one of the operands is unknown, the outcome of the operation is also unknown.

For example, :class:`NA` propagates in arithmetic operations, similarly to
``np.nan``:

.. ipython:: python

   pd.NA + 1
   "a" * pd.NA

There are a few special cases when the result is known, even when one of the
operands is ``NA``.

.. ipython:: python

   pd.NA ** 0
   1 ** pd.NA

In equality and comparison operations, :class:`NA` also propagates. This deviates
from the behaviour of ``np.nan``, where comparisons with ``np.nan`` always
return ``False``.

.. ipython:: python

   pd.NA == 1
   pd.NA == pd.NA
   pd.NA < 2.5

To check if a value is equal to :class:`NA`, use :func:`isna`

.. ipython:: python

   pd.isna(pd.NA)


.. note::

   An exception on this basic propagation rule are *reductions* (such as the
   mean or the minimum), where pandas defaults to skipping missing values. See the
   :ref:`calculation section <missing_data.calculations>` for more.

Logical operations
------------------

For logical operations, :class:`NA` follows the rules of the
`three-valued logic <https://en.wikipedia.org/wiki/Three-valued_logic>`__ (or
*Kleene logic*, similarly to R, SQL and Julia). This logic means to only
propagate missing values when it is logically required.

For example, for the logical "or" operation (``|``), if one of the operands
is ``True``, we already know the result will be ``True``, regardless of the
other value (so regardless the missing value would be ``True`` or ``False``).
In this case, :class:`NA` does not propagate:

.. ipython:: python

   True | False
   True | pd.NA
   pd.NA | True

On the other hand, if one of the operands is ``False``, the result depends
on the value of the other operand. Therefore, in this case :class:`NA`
propagates:

.. ipython:: python

   False | True
   False | False
   False | pd.NA

The behaviour of the logical "and" operation (``&``) can be derived using
similar logic (where now :class:`NA` will not propagate if one of the operands
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
to a boolean value.

.. ipython:: python
   :okexcept:

   bool(pd.NA)

This also means that :class:`NA` cannot be used in a context where it is
evaluated to a boolean, such as ``if condition: ...`` where ``condition`` can
potentially be :class:`NA`. In such cases, :func:`isna` can be used to check
for :class:`NA` or ``condition`` being :class:`NA` can be avoided, for example by
filling missing values beforehand.

A similar situation occurs when using :class:`Series` or :class:`DataFrame` objects in ``if``
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

.. _missing_data.NA.conversion:

Conversion
^^^^^^^^^^

If you have a :class:`DataFrame` or :class:`Series` using ``np.nan``,
:meth:`Series.convert_dtypes` and :meth:`DataFrame.convert_dtypes`
in :class:`DataFrame` that can convert data to use the data types that use :class:`NA`
such as :class:`Int64Dtype` or :class:`ArrowDtype`. This is especially helpful after reading
in data sets from IO methods where data types were inferred.

In this example, while the dtypes of all columns are changed, we show the results for
the first 10 columns.

.. ipython:: python

   import io
   data = io.StringIO("a,b\n,True\n2,")
   df = pd.read_csv(data)
   df.dtypes
   df_conv = df.convert_dtypes()
   df_conv
   df_conv.dtypes

.. _missing.inserting:

Inserting missing data
~~~~~~~~~~~~~~~~~~~~~~

You can insert missing values by simply assigning to a :class:`Series` or :class:`DataFrame`.
The missing value sentinel used will be chosen based on the dtype.

.. ipython:: python

   ser = pd.Series([1., 2., 3.])
   ser.loc[0] = None
   ser

   ser = pd.Series([pd.Timestamp("2021"), pd.Timestamp("2021")])
   ser.iloc[0] = np.nan
   ser

   ser = pd.Series([True, False], dtype="boolean[pyarrow]")
   ser.iloc[0] = None
   ser

For ``object`` types, pandas will use the value given:

.. ipython:: python

   s = pd.Series(["a", "b", "c"], dtype=object)
   s.loc[0] = None
   s.loc[1] = np.nan
   s

.. _missing_data.calculations:

Calculations with missing data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Missing values propagate through arithmetic operations between pandas objects.

.. ipython:: python

   ser1 = pd.Series([np.nan, np.nan, 2, 3])
   ser2 = pd.Series([np.nan, 1, np.nan, 4])
   ser1
   ser2
   ser1 + ser2

The descriptive statistics and computational methods discussed in the
:ref:`data structure overview <basics.stats>` (and listed :ref:`here
<api.series.stats>` and :ref:`here <api.dataframe.stats>`) are all
account for missing data.

When summing data, NA values or empty data will be treated as zero.

.. ipython:: python

   pd.Series([np.nan]).sum()
   pd.Series([], dtype="float64").sum()

When taking the product, NA values or empty data will be treated as 1.

.. ipython:: python

   pd.Series([np.nan]).prod()
   pd.Series([], dtype="float64").prod()

Cumulative methods like :meth:`~DataFrame.cumsum` and :meth:`~DataFrame.cumprod`
ignore NA values by default preserve them in the result. This behavior can be changed
with ``skipna``

* Cumulative methods like :meth:`~DataFrame.cumsum` and :meth:`~DataFrame.cumprod` ignore NA values by default, but preserve them in the resulting arrays. To override this behaviour and include NA values, use ``skipna=False``.


.. ipython:: python

   ser = pd.Series([1, np.nan, 3, np.nan])
   ser
   ser.cumsum()
   ser.cumsum(skipna=False)

.. _missing_data.dropna:

Dropping missing data
~~~~~~~~~~~~~~~~~~~~~

:meth:`~DataFrame.dropna` dropa rows or columns with missing data.

.. ipython:: python

   df = pd.DataFrame([[np.nan, 1, 2], [1, 2, np.nan], [1, 2, 3]])
   df
   df.dropna()
   df.dropna(axis=1)

   ser = pd.Series([1, pd.NA], dtype="int64[pyarrow]")
   ser.dropna()

Filling missing data
~~~~~~~~~~~~~~~~~~~~

.. _missing_data.fillna:

Filling by value
----------------

:meth:`~DataFrame.fillna` replaces NA values with non-NA data.

Replace NA with a scalar value

.. ipython:: python

   data = {"np": [1.0, np.nan, np.nan, 2], "arrow": pd.array([1.0, pd.NA, pd.NA, 2], dtype="float64[pyarrow]")}
   df = pd.DataFrame(data)
   df
   df.fillna(0)

Fill gaps forward or backward

.. ipython:: python

   df.ffill()
   df.bfill()

.. _missing_data.fillna.limit:

Limit the number of NA values filled

.. ipython:: python

   df.ffill(limit=1)

NA values can be replaced with corresponding value from a :class:`Series` or :class:`DataFrame`
where the index and column aligns between the original object and the filled object.

.. ipython:: python

   dff = pd.DataFrame(np.arange(30, dtype=np.float64).reshape(10, 3), columns=list("ABC"))
   dff.iloc[3:5, 0] = np.nan
   dff.iloc[4:6, 1] = np.nan
   dff.iloc[5:8, 2] = np.nan
   dff
   dff.fillna(dff.mean())

.. note::

   :meth:`DataFrame.where` can also be used to fill NA values.Same result as above.

   .. ipython:: python

      dff.where(pd.notna(dff), dff.mean(), axis="columns")


.. _missing_data.interpolate:

Interpolation
-------------

:meth:`DataFrame.interpolate` and :meth:`Series.interpolate` fills NA values
using various interpolation methods.

.. ipython:: python

   df = pd.DataFrame(
       {
           "A": [1, 2.1, np.nan, 4.7, 5.6, 6.8],
           "B": [0.25, np.nan, np.nan, 4, 12.2, 14.4],
       }
   )
   df
   df.interpolate()

   idx = pd.date_range("2020-01-01", periods=10, freq="D")
   data = np.random.default_rng(2).integers(0, 10, 10).astype(np.float64)
   ts = pd.Series(data, index=idx)
   ts.iloc[[1, 2, 5, 6, 9]] = np.nan

   ts
   @savefig series_before_interpolate.png
   ts.plot()

.. ipython:: python

   ts.interpolate()
   @savefig series_interpolate.png
   ts.interpolate().plot()

Interpolation relative to a :class:`Timestamp` in the :class:`DatetimeIndex`
is available by setting ``method="time"``

.. ipython:: python

   ts2 = ts.iloc[[0, 1, 3, 7, 9]]
   ts2
   ts2.interpolate()
   ts2.interpolate(method="time")

For a floating-point index, use ``method='values'``:

.. ipython:: python

   idx = [0.0, 1.0, 10.0]
   ser = pd.Series([0.0, np.nan, 10.0], idx)
   ser
   ser.interpolate()
   ser.interpolate(method="values")

If you have scipy_ installed, you can pass the name of a 1-d interpolation routine to ``method``.
as specified in the scipy interpolation documentation_ and reference guide_.
The appropriate interpolation method will depend on the data type.

.. tip::

   If you are dealing with a time series that is growing at an increasing rate,
   use ``method='barycentric'``.

   If you have values approximating a cumulative distribution function,
   use ``method='pchip'``.

   To fill missing values with goal of smooth plotting use ``method='akima'``.

   .. ipython:: python

      df = pd.DataFrame(
         {
            "A": [1, 2.1, np.nan, 4.7, 5.6, 6.8],
            "B": [0.25, np.nan, np.nan, 4, 12.2, 14.4],
         }
      )
      df
      df.interpolate(method="barycentric")
      df.interpolate(method="pchip")
      df.interpolate(method="akima")

When interpolating via a polynomial or spline approximation, you must also specify
the degree or order of the approximation:

.. ipython:: python

   df.interpolate(method="spline", order=2)
   df.interpolate(method="polynomial", order=2)

Comparing several methods.

.. ipython:: python

   np.random.seed(2)

   ser = pd.Series(np.arange(1, 10.1, 0.25) ** 2 + np.random.randn(37))
   missing = np.array([4, 13, 14, 15, 16, 17, 18, 20, 29])
   ser.iloc[missing] = np.nan
   methods = ["linear", "quadratic", "cubic"]

   df = pd.DataFrame({m: ser.interpolate(method=m) for m in methods})
   @savefig compare_interpolations.png
   df.plot()

Interpolating new observations from expanding data with :meth:`Series.reindex`.

.. ipython:: python

   ser = pd.Series(np.sort(np.random.uniform(size=100)))

   # interpolate at new_index
   new_index = ser.index.union(pd.Index([49.25, 49.5, 49.75, 50.25, 50.5, 50.75]))
   interp_s = ser.reindex(new_index).interpolate(method="pchip")
   interp_s.loc[49:51]

.. _scipy: https://scipy.org/
.. _documentation: https://docs.scipy.org/doc/scipy/reference/interpolate.html#univariate-interpolation
.. _guide: https://docs.scipy.org/doc/scipy/tutorial/interpolate.html

.. _missing_data.interp_limits:

Interpolation limits
^^^^^^^^^^^^^^^^^^^^

:meth:`~DataFrame.interpolate` accepts a ``limit`` keyword
argument to limit the number of consecutive ``NaN`` values
filled since the last valid observation

.. ipython:: python

   ser = pd.Series([np.nan, np.nan, 5, np.nan, np.nan, np.nan, 13, np.nan, np.nan])
   ser
   ser.interpolate()
   ser.interpolate(limit=1)

By default, ``NaN`` values are filled in a ``forward`` direction. Use
``limit_direction`` parameter to fill ``backward`` or from ``both`` directions.

.. ipython:: python

   ser.interpolate(limit=1, limit_direction="backward")
   ser.interpolate(limit=1, limit_direction="both")
   ser.interpolate(limit_direction="both")

By default, ``NaN`` values are filled whether they are surrounded by
existing valid values or outside existing valid values. The ``limit_area``
parameter restricts filling to either inside or outside values.

.. ipython:: python

   # fill one consecutive inside value in both directions
   ser.interpolate(limit_direction="both", limit_area="inside", limit=1)

   # fill all consecutive outside values backward
   ser.interpolate(limit_direction="backward", limit_area="outside")

   # fill all consecutive outside values in both directions
   ser.interpolate(limit_direction="both", limit_area="outside")

.. _missing_data.replace:

Replacing values
----------------

:meth:`Series.replace` and :meth:`DataFrame.replace` can be used similar to
:meth:`Series.fillna` and :meth:`DataFrame.fillna` to replace or insert missing values.

.. ipython:: python

   df = pd.DataFrame(np.eye(3))
   df
   df_missing = df.replace(0, np.nan)
   df_missing
   df_filled = df_missing.replace(np.nan, 2)
   df_filled

Replacing more than one value is possible by passing a list.

.. ipython:: python

   df_filled.replace([1, 44], [2, 28])

Replacing using a mapping dict.

.. ipython:: python

   df_filled.replace({1: 44, 2: 28})

.. _missing_data.replace_expression:

Regular expression replacement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   Python strings prefixed with the ``r`` character such as ``r'hello world'``
   are `"raw" strings <https://docs.python.org/3/reference/lexical_analysis.html#string-and-bytes-literals>`_.
   They have different semantics regarding backslashes than strings without this prefix.
   Backslashes in raw strings  will be interpreted as an escaped backslash, e.g., ``r'\' == '\\'``.

Replace the '.' with ``NaN``

.. ipython:: python

   d = {"a": list(range(4)), "b": list("ab.."), "c": ["a", "b", np.nan, "d"]}
   df = pd.DataFrame(d)
   df.replace(".", np.nan)

Replace the '.' with ``NaN`` with regular expression that removes surrounding whitespace

.. ipython:: python

   df.replace(r"\s*\.\s*", np.nan, regex=True)

Replace with a list of regexes.

.. ipython:: python

   df.replace([r"\.", r"(a)"], ["dot", r"\1stuff"], regex=True)

Replace with a regex in a mapping dict.

.. ipython:: python

   df.replace({"b": r"\s*\.\s*"}, {"b": np.nan}, regex=True)

Pass nested dictionaries of regular expressions that use the ``regex`` keyword.

.. ipython:: python

   df.replace({"b": {"b": r""}}, regex=True)
   df.replace(regex={"b": {r"\s*\.\s*": np.nan}})
   df.replace({"b": r"\s*(\.)\s*"}, {"b": r"\1ty"}, regex=True)

Pass a list of regular expressions that will replace matches with a scalar.

.. ipython:: python

   df.replace([r"\s*\.\s*", r"a|b"], "placeholder", regex=True)

All of the regular expression examples can also be passed with the
``to_replace`` argument as the ``regex`` argument. In this case the ``value``
argument must be passed explicitly by name or ``regex`` must be a nested
dictionary.

.. ipython:: python

   df.replace(regex=[r"\s*\.\s*", r"a|b"], value="placeholder")

.. note::

   A regular expression object from ``re.compile`` is a valid input as well.
