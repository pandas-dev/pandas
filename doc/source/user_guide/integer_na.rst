.. currentmodule:: pandas

{{ header }}

.. _integer_na:

**************************
Nullable integer data type
**************************

.. note::

   IntegerArray is currently experimental. Its API or implementation may
   change without warning. Uses :attr:`pandas.NA` as the missing value.

In :ref:`missing_data`, we saw that pandas primarily uses ``NaN`` to represent
missing data. Because ``NaN`` is a float, this forces an array of integers with
any missing values to become floating point. In some cases, this may not matter
much. But if your integer column is, say, an identifier, casting to float can
be problematic. Some integers cannot even be represented as floating point
numbers.

Construction
------------

pandas can represent integer data with possibly missing values using
:class:`arrays.IntegerArray`. This is an :ref:`extension type <extending.extension-types>`
implemented within pandas.

.. ipython:: python

   arr = pd.array([1, 2, None], dtype=pd.Int64Dtype())
   arr

Or the string alias ``"Int64"`` (note the capital ``"I"``) to differentiate from
NumPy's ``'int64'`` dtype:

.. ipython:: python

   pd.array([1, 2, np.nan], dtype="Int64")

All NA-like values are replaced with :attr:`pandas.NA`.

.. ipython:: python

   pd.array([1, 2, np.nan, None, pd.NA], dtype="Int64")

This array can be stored in a :class:`DataFrame` or :class:`Series` like any
NumPy array.

.. ipython:: python

   pd.Series(arr)

You can also pass the list-like object to the :class:`Series` constructor
with the dtype.

.. warning::

   Currently :meth:`pandas.array` and :meth:`pandas.Series` use different
   rules for dtype inference. :meth:`pandas.array` will infer a
   nullable-integer dtype

   .. ipython:: python

      pd.array([1, None])
      pd.array([1, 2])

   For backwards-compatibility, :class:`Series` infers these as either
   integer or float dtype.

   .. ipython:: python

      pd.Series([1, None])
      pd.Series([1, 2])

   We recommend explicitly providing the dtype to avoid confusion.

   .. ipython:: python

      pd.array([1, None], dtype="Int64")
      pd.Series([1, None], dtype="Int64")

   In the future, we may provide an option for :class:`Series` to infer a
   nullable-integer dtype.

If you create a column of ``NA`` values (for example to fill them later)
with ``df['new_col'] = pd.NA``, the ``dtype`` would be set to ``object`` in the
new column. The performance on this column will be worse than with
the appropriate type. It's better to use
``df['new_col'] = pd.Series(pd.NA, dtype="Int64")``
(or another ``dtype`` that supports ``NA``).

.. ipython:: python

   df = pd.DataFrame()
   df['objects'] = pd.NA
   df.dtypes

Operations
----------

Operations involving an integer array will behave similar to NumPy arrays.
Missing values will be propagated, and the data will be coerced to another
dtype if needed.

.. ipython:: python

   s = pd.Series([1, 2, None], dtype="Int64")

   # arithmetic
   s + 1

   # comparison
   s == 1

   # slicing operation
   s.iloc[1:3]

   # operate with other dtypes
   s + s.iloc[1:3].astype("Int8")

   # coerce when needed
   s + 0.01

These dtypes can operate as part of a ``DataFrame``.

.. ipython:: python

   df = pd.DataFrame({"A": s, "B": [1, 1, 3], "C": list("aab")})
   df
   df.dtypes


These dtypes can be merged, reshaped & casted.

.. ipython:: python

   pd.concat([df[["A"]], df[["B", "C"]]], axis=1).dtypes
   df["A"].astype(float)

Reduction and groupby operations such as :meth:`~DataFrame.sum` work as well.

.. ipython:: python

   df.sum(numeric_only=True)
   df.sum()
   df.groupby("B").A.sum()

Scalar NA value
---------------

:class:`arrays.IntegerArray` uses :attr:`pandas.NA` as its scalar
missing value. Slicing a single element that's missing will return
:attr:`pandas.NA`

.. ipython:: python

   a = pd.array([1, None], dtype="Int64")
   a[1]
