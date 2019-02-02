.. currentmodule:: pandas

{{ header }}

.. _integer_na:

**************************
Nullable Integer Data Type
**************************

.. versionadded:: 0.24.0

.. note::

   IntegerArray is currently experimental. Its API or implementation may
   change without warning.


In :ref:`missing_data`, we saw that pandas primarily uses ``NaN`` to represent
missing data. Because ``NaN`` is a float, this forces an array of integers with
any missing values to become floating point. In some cases, this may not matter
much. But if your integer column is, say, an identifier, casting to float can
be problematic. Some integers cannot even be represented as floating point
numbers.

Pandas can represent integer data with possibly missing values using
:class:`arrays.IntegerArray`. This is an :ref:`extension types <extending.extension-types>`
implemented within pandas. It is not the default dtype for integers, and will not be inferred;
you must explicitly pass the dtype into :meth:`array` or :class:`Series`:

.. ipython:: python

   arr = pd.array([1, 2, np.nan], dtype=pd.Int64Dtype())
   arr

Or the string alias ``"Int64"`` (note the capital ``"I"``, to differentiate from
NumPy's ``'int64'`` dtype:

.. ipython:: python

   pd.array([1, 2, np.nan], dtype="Int64")

This array can be stored in a :class:`DataFrame` or :class:`Series` like any
NumPy array.

.. ipython:: python

   pd.Series(arr)

You can also pass the list-like object to the :class:`Series` constructor
with the dtype.

.. ipython:: python

   s = pd.Series([1, 2, np.nan], dtype="Int64")
   s

By default (if you don't specify ``dtype``), NumPy is used, and you'll end
up with a ``float64`` dtype Series:

.. ipython:: python

   pd.Series([1, 2, np.nan])

Operations involving an integer array will behave similar to NumPy arrays.
Missing values will be propagated, and and the data will be coerced to another
dtype if needed.

.. ipython:: python

   # arithmetic
   s + 1

   # comparison
   s == 1

   # indexing
   s.iloc[1:3]

   # operate with other dtypes
   s + s.iloc[1:3].astype('Int8')

   # coerce when needed
   s + 0.01

These dtypes can operate as part of of ``DataFrame``.

.. ipython:: python

   df = pd.DataFrame({'A': s, 'B': [1, 1, 3], 'C': list('aab')})
   df
   df.dtypes


These dtypes can be merged & reshaped & casted.

.. ipython:: python

   pd.concat([df[['A']], df[['B', 'C']]], axis=1).dtypes
   df['A'].astype(float)

Reduction and groupby operations such as 'sum' work as well.

.. ipython:: python

   df.sum()
   df.groupby('B').A.sum()
